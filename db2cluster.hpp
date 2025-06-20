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

using namespace std;

struct vec3 {
    float x, y, z;
    vec3(float _x, float _y, float _z) : x(_x), y(_y), z(_z) {}
    vec3() : x(0), y(0), z(0) {}
};

float dist(vec3 a, vec3 b);

struct db2group {
    int coords_start;
    int coords_end;
    db2group(int _coords_start, int _coords_end) : coords_start(_coords_start), coords_end(_coords_end) {}
    db2group() {}
};

struct db2set {
    vector<int> group_ids;  
    int hydro;            
    float energy_1;
    float energy_2;
    db2set() {}
};

struct db2coord {
    vec3 xyz;
    int atomi;
    int groupi;
    db2coord(vec3 _xyz, int _atomi, int _groupi) : xyz(_xyz), atomi(_atomi), groupi(_groupi) {}
};

enum dockcolortype {
    POSITIVE = 1,
    NEGATIVE = 2,
    ACCEPTOR = 3,
    DONOR    = 4,
    ESTER_O  = 5,
    AMIDE_O  = 6,
    NEUTRAL  = 7
};

struct db2atom { 
    string name;
    string type;
    int docktype;
    int color;
    float charge;
    float polarsolv;
    float apolarsolv;
    float solv;
    float surface;
    float sigma;
    db2atom() {}
};

struct db2bond {
    int atoma;
    int atomb;
    string bondtype;
    db2bond(int _atoma, int _atomb, string _bondtype) : atoma(_atoma), atomb(_atomb), bondtype(_bondtype) {}
};

struct db2 {
    string name;
    vector<db2atom>         atoms;  // various per-atom data, besides total count this information is not used by db2cluster
    vector<db2bond>         bonds;  // bond data, this information is not used at all by db2cluster
    vector<db2group>        groups; // all groups of atoms that are in an identical position between 1 or more input conformations
                                    // in the original db2 scripts groups are referred to as "confs", not to be mistaken as shorthand for "conformations"
                                    // i found this an odd choice- so i am calling "confs" "groups" now- still vague, but less confusing
    vector<db2set>          sets;   // a set is composed of groups- each set corresponds to an input conformation
    vector<db2coord>        coords; // actual coordinates for reconstructing conformations from sets/groups
    float charge;
    float polarsolv;
    float apolarsolv;
    float solv;
    float surface;
};

db2 db2cluster(db2& db, vector<vector<vec3>> conformationsxyz);