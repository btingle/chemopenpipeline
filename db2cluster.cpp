// standard library (cpdkit's devs should check this out!)
#include <vector>
#include <algorithm>
#include <cmath>
#include <map>
#include <stack>
#include <cassert>
#include <iostream>
#include <set>
#include <fstream>
#include <chrono>
#include <functional>
#include <sstream>
#include "db2cluster.hpp"

using namespace std;

struct vec3entry {
    vec3 v;
    int idx;
    vec3entry(vec3 _v, int _idx) : v(_v), idx(_idx) {}
    vec3entry() {}
};

float dist(vec3 a, vec3 b) {
    vec3 dx = vec3(a.x-b.x, a.y-b.y, a.z-b.z);
    return dx.x*dx.x + dx.y*dx.y + dx.z*dx.z;
}

struct groupentry {
    int gid;
    int idx;
    groupentry(int _gid, int _idx) : gid(_gid), idx(_idx) {}
    groupentry() {}
};

// input: list of 3d vectors, tolerance distance for clustering
// output: group affinity map for each vector- 0 = clustered to vec 0, 1 = clustered to vec1 etc...
// ex:
// tol = 0.5, l = {(0.7, 1, 1), (0.9, 2, 2), (1, 1, 1), (1.3, 2, 2)} => {0, 1, 0, 1}
vector<int> fast3dcluster(vector<vec3> inlist, int n, map<int, int>& group_map, float tol=0.01) {
    float tol2 = tol*tol;

    vector<vec3entry> l;
    for (int i = 0; i < n; i++) {
        l.push_back(vec3entry(inlist[i], i));
    }
    function<bool(groupentry, groupentry)> gcompidx = [](groupentry a, groupentry b) {
        return a.idx < b.idx;
    };
    // something to note:
    // to speed up clustering, we sort by x and only cluster subsets of positions that could plausibly be within tolerance distance of one another
    // for my expected data, this should be a sufficient speed up
    // but one could go further, and for each plausible x subset, sort by y and create plausible subsets, and so on with z
    // then you would be left with the smallest possible groups of positions to cluster
    function<bool(vec3entry, vec3entry)> vcompx = [](vec3entry a, vec3entry b) {
        return a.v.x < b.v.x;
    };
    function<bool(vec3entry, vec3entry)> vcompidx = [](vec3entry a, vec3entry b) {
        return a.idx < b.idx;
    };
    sort(l.begin(), l.end(), vcompx);

    int i = 1;
    int ci = 0;
    int clustercount = 0;
    vector<groupentry> s(n);
    for (int j = 0; j < n; j++) {
        s[j] = groupentry(0, l[j].idx);
    }

    auto _fclust = [&] () { // & means this lambda inherits references to the current local variables in scope
        stack<int> visited;
        // sort by original index to ensure that group hashes are always calculated the same way
        sort(l.begin()+ci, l.begin()+i, vcompidx);
        for (int j = ci; j < i; j++) {
            if (s[j].gid) { // already grouped
                continue;
            }
            int chash = l[j].idx;
            clustercount += 1;
            visited.push(j);
            for (int k = j+1; k < i; k++) { // simple greedy O(N^2) clustering
                if (dist(l[j].v, l[k].v) <= tol2) {
                    s[k].idx = l[k].idx;
                    chash = chash ^ (l[k].idx << 1);
                    visited.push(k);
                }
            }
            int gid;
            if (visited.size() == l.size()) { // this is only true if there is but one cluster in l
                gid = 1; // assign special gid 1
            } else if (!group_map[chash]) {
                gid = group_map.size()+1; // minimum id for non-rigid = 2
                group_map[chash] = gid;
            } else {
                gid = group_map[chash];
            }
            while (!visited.empty()) {
                int k = visited.top();
                visited.pop();
                s[k].gid = gid;
            }
        }
    };

    while (i < n) {
        vec3 v = l[i].v;
        vec3 pv = l[i-1].v;
        if (!(abs(v.x-pv.x) <= tol)) {
            _fclust();
            ci = i;
        }
        i++;
    }
    if (ci == 0) {
        ci = 0;
    }
    _fclust();

    vector<int> sout;
    sort(s.begin(), s.end(), gcompidx); // return to original sorting
    for (int i = 0; i < s.size(); i++) {
        sout.push_back(s[i].gid);
    }
    return sout;
}

// db should already have sets filled out (without group_id info) as well as atoms and bonds
db2 db2cluster(db2& db, vector<vector<vec3>> conformationsxyz) {

    int natoms = db.atoms.size();
    int nconfs = db.sets.size();
    vector<vector<int>> groups_by_atom;
    map<int, int> group_map; // keeps track of unique groups between clustering operations

    for (int i = 0; i < natoms; i++) {
        vector<int> atom_groups = fast3dcluster(conformationsxyz[i], nconfs, group_map);
        groups_by_atom.push_back(atom_groups);
    }

    int ngroups = group_map.size() + 1; // assuming rigid group was found- will throw error if none is detected later
    db.groups.resize(ngroups);

    map<int, bool> group_visited;
    for (int i = 0; i < nconfs; i++) {
        vector<db2coord> setdata;
        db2set* thisset = &db.sets[i];
        for (int j = 0; j < natoms; j++) { // transpose column (per-atom) data to row (per-conformation) data
            setdata.push_back(db2coord(conformationsxyz[j][i], j, groups_by_atom[j][i]));
        }
        /*struct {
            bool operator()(db2coord a, db2coord b) { return a.groupi < b.groupi; }
        } compcoordgid;*/
        function<bool(db2coord, db2coord)> compcoordgid = [](db2coord a, db2coord b) {
            return a.groupi < b.groupi;
        };
        // now sort by group id- besides having better memory access patterns, it puts our rigid group 1 first
        sort(setdata.begin(), setdata.end(), compcoordgid);
        assert(setdata[0].groupi == 1); // here's where we detect if we have a rigid component
        int coords_index_start = db.coords.size();
        int g, pg = 0;
        function<void()> pushgroup = [&] () { // just to condense the code a little
            thisset->group_ids.push_back(pg);
            // if (!group_visited[pg])
            if (!group_visited[pg]) {
                group_visited[pg] = true;
                db.groups[pg-1] = db2group(coords_index_start, db.coords.size());
                coords_index_start = db.coords.size();
            }
        };
        // assemble groups list for this set, assemble any newly encountered groups, push coordinates of any new groups
        for (int j = 0; j < setdata.size(); j++) {
            g = setdata[j].groupi;
            if (j != 0 && g != pg) {
                pushgroup();
            }
            if (!group_visited[g]) {
                db.coords.push_back(setdata[j]);
            }
            pg = g;
        }
        pushgroup();
    }

    cout << ngroups << " groups found for " << nconfs << " confs" << endl;

    return db;
}