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
#include <iomanip>
#include "db2cluster.hpp"
// openbabel
#include "openbabel/obconversion.h"
// CDPKit
// "standard" functions
#include "CDPL/Chem/Atom.hpp"
#include "CDPL/Chem/AtomContainerFunctions.hpp"
#include "CDPL/Chem/AtomFunctions.hpp"
#include "CDPL/Chem/AtomFunctions.hpp"
#include "CDPL/Chem/BasicMolecule.hpp"
#include "CDPL/Chem/Bond.hpp"
#include "CDPL/Chem/BondFunctions.hpp"
#include "CDPL/Chem/ControlParameterFunctions.hpp"
#include "CDPL/Chem/Entity3DContainerFunctions.hpp"
#include "CDPL/Chem/MolecularGraph.hpp"
#include "CDPL/Chem/MolecularGraphFunctions.hpp"
#include "CDPL/Chem/MoleculeReader.hpp"
#include "CDPL/Chem/MOL2MolecularGraphWriter.hpp"
#include "CDPL/Chem/UtilityFunctions.hpp"
#include "CDPL/Chem/SybylAtomType.hpp"
#include "CDPL/Math/VectorArray.hpp"
// tautomer specific functions
#include "CDPL/Chem/TautomerGenerator.hpp"
#include "CDPL/Chem/TautomerScore.hpp"
#include "CDPL/Chem/KetoEnolTautomerization.hpp"  
#include "CDPL/Chem/ImineEnamineTautomerization.hpp"  
#include "CDPL/Chem/NitrosoOximeTautomerization.hpp"  
#include "CDPL/Chem/NitroAciTautomerization.hpp"  
#include "CDPL/Chem/AmideImidicAcidTautomerization.hpp"  
#include "CDPL/Chem/LactamLactimTautomerization.hpp"  
#include "CDPL/Chem/KeteneYnolTautomerization.hpp"  
#include "CDPL/Chem/PhosphinicAcidTautomerization.hpp"  
#include "CDPL/Chem/SulfenicAcidTautomerization.hpp"  
#include "CDPL/Chem/GenericHydrogen13ShiftTautomerization.hpp"
// conformer specific functions
#include "CDPL/ConfGen/ConformerGenerator.hpp"
#include "CDPL/ConfGen/MoleculeFunctions.hpp"
#include "CDPL/ConfGen/StructureGenerator.hpp"
// protonation specific functions
#include "CDPL/Chem/ProtonationStateStandardizer.hpp"

using namespace CDPL;
using namespace Chem;

bool enumtautomers(BasicMolecule mol, long unsigned int maxtauts, vector<BasicMolecule>& out) {

    // tautomerize, protonate
    function<void(BasicMolecule&)> initmolecule = [] (BasicMolecule& molgraph) {
        setRingFlags(molgraph, false);
        perceiveHybridizationStates(molgraph, false);
        perceiveSSSR(molgraph, false);
        setAromaticityFlags(molgraph, false);
    };
    calcImplicitHydrogenCounts(mol, false);
    initmolecule(mol);
    PatternBasedTautomerizationRule::SharedPointer h13_shift(new GenericHydrogen13ShiftTautomerization());
    PatternBasedTautomerizationRule::SharedPointer keto_enol(new KetoEnolTautomerization());
    PatternBasedTautomerizationRule::SharedPointer imine_enamine(new ImineEnamineTautomerization());  
    PatternBasedTautomerizationRule::SharedPointer nitroso_oxime(new NitrosoOximeTautomerization());  
    PatternBasedTautomerizationRule::SharedPointer nitro_aci(new NitroAciTautomerization());  
    PatternBasedTautomerizationRule::SharedPointer amide_imidic(new AmideImidicAcidTautomerization());  
    PatternBasedTautomerizationRule::SharedPointer lactam_lactim(new LactamLactimTautomerization());  
    PatternBasedTautomerizationRule::SharedPointer ketene_ynol(new KeteneYnolTautomerization());
    PatternBasedTautomerizationRule::SharedPointer posph_acid(new PhosphinicAcidTautomerization());
    PatternBasedTautomerizationRule::SharedPointer sulf_acid(new SulfenicAcidTautomerization());

    h13_shift->addExcludePatterns(*keto_enol);
    h13_shift->addExcludePatterns(*imine_enamine);
    h13_shift->addExcludePatterns(*nitroso_oxime);
    h13_shift->addExcludePatterns(*amide_imidic);
    h13_shift->addExcludePatterns(*lactam_lactim);
    h13_shift->addExcludePatterns(*nitro_aci);
    
    TautomerGenerator tautgen;

    tautgen.addTautomerizationRule(h13_shift);
    tautgen.addTautomerizationRule(keto_enol);
    tautgen.addTautomerizationRule(imine_enamine);
    tautgen.addTautomerizationRule(nitroso_oxime);
    tautgen.addTautomerizationRule(nitro_aci);
    tautgen.addTautomerizationRule(amide_imidic);
    tautgen.addTautomerizationRule(lactam_lactim);
    tautgen.addTautomerizationRule(ketene_ynol);
    tautgen.addTautomerizationRule(posph_acid);
    tautgen.addTautomerizationRule(sulf_acid);

    struct tautentry {
        BasicMolecule mol;
        double score;
        tautentry(BasicMolecule _mol, double _score) : mol(_mol), score(_score) {}
    };
    TautomerScore tautScore;
    vector<tautentry> tautbuffer;

    function<bool(MolecularGraph&)> tautcallback = [&] (MolecularGraph& taut) {
        initmolecule((BasicMolecule&)taut);
        double score = tautScore(taut);
        tautbuffer.push_back(tautentry((BasicMolecule)taut, score));
        if (tautbuffer.size() > 100) {
            return false;
        }
        return true;
    };

    tautgen.setCallbackFunction(tautcallback);
    tautgen.generate(mol);

    function<bool(tautentry, tautentry)> sortscore = [] (tautentry a, tautentry b) {
        return a.score < b.score;
    };

    sort(tautbuffer.begin(), tautbuffer.end(), sortscore);
    for (int i = 0; i < min(maxtauts, tautbuffer.size()); i++) {
        ProtonationStateStandardizer prot_state_gen;
        prot_state_gen.standardize(tautbuffer[i].mol, ProtonationStateStandardizer::PHYSIOLOGICAL_CONDITION_STATE);
        perceiveComponents(tautbuffer[i].mol);
        out.push_back(tautbuffer[i].mol);
    }
    return true;
}

bool generateconformers(Molecule& mol) {

    ConfGen::ConformerGenerator struct_gen;
    ConfGen::ConformerGeneratorSettings struct_gen_settings;
    struct_gen_settings.setTimeout(1000*5);
    struct_gen_settings.setSamplingMode(1);
    struct_gen_settings.setMaxPoolSize(2000);
    struct_gen_settings.setMacrocycleRotorBondCountThreshold(10);
    // this is a cdpkit function frustratingly locked behind an inacecssible header file
    // it is not too complicated though, so just define it here
    // usually conformer generation switches to DG stochastic mode beyond this threshold
    // we would rather avoid processing molecules like this, so just skip them if this is the case
    function<int(FragmentList&)> getMaxNonAromaticSingleBondCount = [](FragmentList& frags) {
        int maxcount = 0;
        for (auto& ring : frags) {
            int count = 0;
            for (auto& bond : ring.getBonds()) {
                if (getOrder(bond) == 1 && !getAromaticityFlag(bond)) {
                    count++;
                }
            }
            maxcount = max(count, maxcount);
        }
        return maxcount;
    };
    if (getMaxNonAromaticSingleBondCount(*getSSSR(mol)) > struct_gen_settings.getMacrocycleRotorBondCountThreshold()) {
        cout << "large macrocycle detected, skipping this molecule" << endl;
        return false;
    }
    ConfGen::prepareForConformerGeneration(mol);
    struct_gen.generate(mol);
    struct_gen.setConformers(mol);
    if (!has3DCoordinatesArray(mol.getAtom(0))) {
        cout << "failed to generate conformers for this molecule" << endl;
        return false;
    }
    set3DCoordinates(mol, struct_gen.getConformer(0));
    return true;
}

bool calcsolvation(MolecularGraph& mol, vector<db2atom>& atom_solv_data_out, db2& outtotal) {
    using namespace OpenBabel;
    stringstream mol2_string;
    stringstream mopin_wat_string;
    stringstream mopin_hex_string;
    MOL2MolecularGraphWriter mol_writer(mol2_string);
    setMultiConfExportParameter(mol_writer, false);
    mol_writer.write(mol);
    OBConversion conv;
    conv.SetInAndOutFormats("mol2", "mopin");
    stringstream watkeywords;
    watkeywords <<
    "CHARGE=-9.999999999983633e-05 AM1 1SCF TLIMIT=3 GEO-OK SM5.42R\n"
    "& SOLVNT=WATER\n" << "dummy" << " " << mol.getNumAtoms();
    conv.AddOption("k", OBConversion::OUTOPTIONS, watkeywords.str().c_str());
    conv.Convert(&mol2_string, &mopin_wat_string);
    stringstream hexkeywords;
    hexkeywords <<
    "CHARGE=-9.999999999983633e-05 AM1 1SCF TLIMIT=3 GEO-OK SM5.42R\n"
    "& SOLVNT=GENORG IOFR=1.4345 ALPHA=0.00 BETA=0.00 GAMMA=38.93\n"
    "& DIELEC=2.06 FACARB=0.00 FEHALO=0.00 DEV\n" << "dummy" << " " << mol.getNumAtoms();
    mol2_string = stringstream(mol2_string.str()); // seek() doesnt work, but this does!
    conv.RemoveOption("k", OBConversion::OUTOPTIONS);
    conv.AddOption("k", OBConversion::OUTOPTIONS, hexkeywords.str().c_str());
    conv.Convert(&mol2_string, &mopin_hex_string);

    // i had an idea to do the stdin and stdout via fancy stringstream
    function<stringstream(stringstream&)> run_amsol = [] (stringstream& indata) {
        filebuf fb;
        fb.open("temp.in", ios::out);
        ostream os(&fb);
        os << indata.str();
        fb.close();
        system("/home/btingle/chemopenpipeline/amsol7.1/amsol < temp.in > temp.ot 2>/dev/null");
        stringstream outdata;
        fb.open("temp.ot", ios::in);
        istream is(&fb);
        outdata << is.rdbuf();
        return outdata;
    };

    // very scrungly function to grok the useful data from amsol's overly verbose output
    // 
    function<void(stringstream&, vector<db2atom>&, db2atom&)> extract_data = 
    [] (stringstream& indata, vector<db2atom>& outdata, db2atom& outtotal) {
        string line;
        while (getline(indata, line)) {
            stringstream linestream(line);
            string token;
            linestream >> token;
            // we are interested in some of these stats
            if (token == "Total:") {
                float dummy;
                linestream >> outtotal.charge;
                linestream >> dummy;
                linestream >> outtotal.surface;
                linestream >> outtotal.apolarsolv;
                linestream >> dummy;
                break; // we are done once we find this line
            }
            if (token.empty() || !isdigit(token[0]) || (token.size() > 1 && !isdigit(token[1]))) {
                continue;
            }
            linestream >> token;
            if (token.empty() || !isalpha(token[0])) {
                continue;
            }
            db2atom atom;
            atom.name[0] = token[0];
            linestream >> atom.charge;
            linestream >> atom.polarsolv;
            if (linestream.eof()) { // there should be more fields, if there aren't we're at the wrong place
                continue;
            }
            linestream >> atom.surface;
            linestream >> atom.sigma;
            linestream >> atom.apolarsolv;
            linestream >> atom.solv;
            outdata.push_back(atom);
        }
        return outdata;
    };

    stringstream solv_wat = run_amsol(mopin_wat_string);
    stringstream solv_hex = run_amsol(mopin_hex_string);
    vector<db2atom> atoms_wat;
    db2atom total_wat;
    extract_data(solv_wat, atoms_wat, total_wat);
    vector<db2atom> atoms_hex;
    db2atom total_hex;
    extract_data(solv_hex, atoms_hex, total_hex);

    if (atoms_hex.size() != mol.getNumAtoms()) {
        cout << "DETECTED INCOMPLETE DATA! " << atoms_hex.size() << "<->" << mol.getNumAtoms() << endl;
        return false;
    }
    // this is all copied from the original UCSF implementation
    float sum_apolar = 0;
    for (db2atom& atom_solv_hex : atoms_hex) {
        sum_apolar += atom_solv_hex.apolarsolv;
    }
    float cs_coeff = (total_hex.apolarsolv - sum_apolar) / total_hex.surface;
    float tot_diff_polar = 0;
    float tot_diff_apolar = 0;
    float tot_diff_allsolv = 0;
    for (int i = 0; i < mol.getNumAtoms(); i++) {
        db2atom atom_solv_wat = atoms_wat[i];
        db2atom atom_solv_hex = atoms_hex[i];
        atom_solv_wat.polarsolv  -= atom_solv_hex.polarsolv;
        atom_solv_wat.apolarsolv -= atom_solv_hex.apolarsolv + cs_coeff * atom_solv_hex.surface;
        atom_solv_wat.solv        = atom_solv_wat.polarsolv + atom_solv_wat.apolarsolv;
        atom_solv_wat.surface     = atom_solv_hex.surface;
        tot_diff_polar += atom_solv_wat.polarsolv;
        tot_diff_apolar += atom_solv_wat.apolarsolv;
        tot_diff_allsolv += atom_solv_wat.solv;
    }
    outtotal.charge = total_hex.charge;
    outtotal.polarsolv = tot_diff_polar;
    outtotal.apolarsolv = tot_diff_apolar;
    outtotal.surface = total_hex.surface;
    outtotal.solv = tot_diff_allsolv;
    for (auto atom : atoms_wat) {
        atom_solv_data_out.push_back(atom);
    }
    return true;
}

bool generatedb2(MolecularGraph& mol, vector<db2atom> atoms_solv_data, db2& db) {
    using namespace SybylAtomType;
    int natoms = mol.getNumAtoms();
    int nbonds = mol.getNumBonds();
    int nconfs = getNumConformations(mol);
    vector<vector<vec3>> conformationsxyz;
    vector<int> dock_colors;
    vector<int> dock_numbers;
    // a rule determines what atoms are assigned which dock color types
    struct dockcolorrule {
        int atom_sybyl_type;
        // bond dist <0 = only apply to atoms not bonded to bonded_sybyl_type within abs(bond_dist)
        // bond dist 0  = dont consider bonded atoms
        // bond dist 1+ = apply to atoms bonded to bonded_sybyl_type within bond_dist
        int bond_dist;
        int bonded_sybyl_type;
        int dock_color_type; // if rule conditions met, apply this color type
        dockcolorrule(int _atom_sybyl_type, int _bond_dist, int _bonded_sybyl_type, int _dock_color_type) :
            atom_sybyl_type(_atom_sybyl_type), bond_dist(_bond_dist), bonded_sybyl_type(_bonded_sybyl_type), dock_color_type(_dock_color_type) {}
    };
    vector<dockcolorrule> rules = {
        dockcolorrule(N_4  ,  0, 0    , dockcolortype::POSITIVE),
        dockcolorrule(O_co2,  0, 0    , dockcolortype::NEGATIVE),
        dockcolorrule(O_2  ,  0, 0    , dockcolortype::ACCEPTOR),
        dockcolorrule(O_3  ,  0, 0    , dockcolortype::ACCEPTOR),
        dockcolorrule(S_2  ,  0, 0    , dockcolortype::ACCEPTOR),
        dockcolorrule(N_ar ,  0, 0    , dockcolortype::ACCEPTOR),
        dockcolorrule(P_3  ,  1, O_co2, dockcolortype::NEGATIVE),
        dockcolorrule(S_O2 ,  1, O_co2, dockcolortype::NEGATIVE),
        dockcolorrule(N_2  ,  1, H    , dockcolortype::DONOR),
        dockcolorrule(N_am ,  1, H    , dockcolortype::DONOR),
        dockcolorrule(N_pl3,  1, H    , dockcolortype::DONOR),
        dockcolorrule(O_3  ,  1, H    , dockcolortype::DONOR),
        dockcolorrule(N_ar , -1, H    , dockcolortype::ACCEPTOR),
        dockcolorrule(N_ar , -1, C_3  , dockcolortype::ACCEPTOR),
        dockcolorrule(N_ar ,  1, H    , dockcolortype::DONOR),
        dockcolorrule(O_3  ,  1, H    , dockcolortype::DONOR),
        dockcolorrule(O_2  ,  2, O_3  , dockcolortype::ESTER_O),
        dockcolorrule(O_2  ,  2, N_pl3, dockcolortype::AMIDE_O),
        dockcolorrule(O_2  ,  2, N_am , dockcolortype::AMIDE_O),
        dockcolorrule(O_2  ,  2, N_3  , dockcolortype::AMIDE_O)
    };
    // special types for DOCK
    map<int, int> sybyl_to_docktype = {
        {C_3, 5},    {C_2, 1},    {C_ar, 1},   {C_1, 1},    {N_3, 10},   {N_2, 8},    
        {N_1, 8},    {O_3, 12},   {O_2, 11},   {S_3, 14},   {N_ar, 8},   {P_3, 13},   
        {H, 6},      {Co_oh, 2},  {Br, 17},    {Cl, 16},    {F, 15},     {I, 18},     
        {S_2, 14},   {N_pl3, 8},  {LP, 25},    {Na, 19},    {K, 19},     {Ca, 21},    
        {Li, 20},    {Al, 20},    {Du, 25},    {Du_C, 25},  {Si, 24},    {N_am, 8},   
        {S_O, 14},   {S_O2, 14},  {N_4, 9},    {O_co2, 11}, {C_cat, 1},  {H_spc, 6},  
        {O_spc, 11}, {H_t3p, 6},  {O_t3p, 11}, {Any, 25},   {Hev, 25},   {Het, 25},   
        {Hal, 25},   {Mg, 20},    {Cr_oh, 25}, {Cr_th, 25}, {Se, 25},    {Fe, 25},    
        {Cu, 25},    {Zn, 26},    {Sn, 25},    {Mo, 25},    {Mn, 25}
    };
    // find color types & dock types
    for (auto& atom : mol.getAtoms()) {
        function<bool(Atom&, int, int, set<int>&)> eval_atoms_in_range;
        // graph search all atoms within dmax range
        eval_atoms_in_range = [&] (Atom& atom, int cdist, int dmax, set<int>& type_query) {
            int sybyl_type = perceiveSybylType(atom, mol);
            if (type_query.count(sybyl_type))
                return true;
            else if (cdist >= dmax)
                return false;
            for (auto& connected_atom : atom.getAtoms()) {
                if (eval_atoms_in_range(connected_atom, cdist+1, dmax, type_query))
                    return true;
            }
            return false;
        };
        int sybyl_type = perceiveSybylType(atom, mol);
        int lastcolor = 0; 
        for (auto rule : rules) {
            set<int> query {rule.dock_color_type};
            if (!(sybyl_type == rule.atom_sybyl_type))
                continue;
            if (rule.bond_dist == 0) {
                lastcolor = rule.dock_color_type;
                continue;
            }
            if (eval_atoms_in_range(atom, 0, abs(rule.bond_dist), query)) {
                if (rule.bond_dist > 0)
                    lastcolor = rule.dock_color_type; // positive rules
            } else if (rule.bond_dist < 0) {
                lastcolor = rule.dock_color_type; // negative rules
            }
                
        }
        dock_colors.push_back(lastcolor);
        int dock_type = sybyl_to_docktype[sybyl_type];
        // special case for hydrogen atoms connected to carbons 
        // they are not recognized by SYBYL, but are by DOCK- so find them here
        if (dock_type == 6) {
            set<int> c_query {C_1, C_2, C_3, C_ar, C_cat};
            if (eval_atoms_in_range(atom, 0, 1, c_query)) {
                dock_type = 7; // special dock type 7 for these cases
            }
        }
        dock_numbers.push_back(dock_type);
    } 

    for (auto& atom : mol.getAtoms()) {
        vector<vec3> atomxyz;
        db2atom atom_solv_data = atoms_solv_data[atom.getIndex()];
        Math::Vector3DArray::SharedPointer atom_coordinates;
        atom_coordinates = get3DCoordinatesArray(atom);
        for (int j = 0; j < nconfs; j++) {
            Math::Vector3D xyz = atom_coordinates->getElement(j);
            atomxyz.push_back(vec3(xyz[0], xyz[1], xyz[2]));

        }
        db2atom dbatom = db2atom(); 
        // todo: properly fill out all these fields
        dbatom.name = getSymbolForType(atom)+to_string(atom.getIndex()+1);
        dbatom.type = getSybylAtomTypeString(perceiveSybylType(atom, mol));
        dbatom.docktype     = dock_numbers[atom.getIndex()];
        dbatom.color        = dock_colors[atom.getIndex()];
        dbatom.charge       = atom_solv_data.charge;
        dbatom.polarsolv    = atom_solv_data.polarsolv;  
        dbatom.apolarsolv   = atom_solv_data.apolarsolv;
        dbatom.solv         = atom_solv_data.solv;
        dbatom.surface      = atom_solv_data.surface;
        db.atoms.push_back(dbatom);
        conformationsxyz.push_back(atomxyz);
    }

    for (auto& bond : mol.getBonds()) {
        string bondtype = getSybylBondTypeString(perceiveSybylType(bond, mol));
        db2bond db_bond(bond.getBegin().getIndex(), bond.getEnd().getIndex(), bondtype);
        db.bonds.push_back(db_bond);
    }

    auto energies = getConformerEnergies(mol);
    for (int i = 0; i < nconfs; i++) {
        db2set dbset = db2set();
        dbset.hydro = 0; // weird flag, ignoring it
        dbset.energy_1 = energies->getElement(i);
        dbset.energy_2 = 0;
        db.sets.push_back(dbset);
    }
    // finally, actually cluster the db2
    db2cluster(db, conformationsxyz);
    return true;
}

bool writedb2(db2& db, ostream& os) {
    char rigbuf[512]; // save this for a little later, we want to calculate something first
    int numrigidheavy = 0; // this value, for M header
    int buff_offset = 0;
    for (int i = db.groups[0].coords_start; i < db.groups[0].coords_end; i++) {
        auto coord = db.coords[i];
        auto atom = db.atoms[coord.atomi];
        if (atom.name[0] == 'H')
            continue;
        buff_offset += sprintf(rigbuf+buff_offset, "R %6d %2d %+9.4f %+9.4f %+9.4f\n", numrigidheavy+1, 
            atom.color, coord.xyz.x, coord.xyz.y, coord.xyz.z);
        numrigidheavy++;
    }
    char buf[512];
    sprintf(buf, "M %16s %9s %3d %3d %6d %6d %6d %6d %6d %6d\n",
        db.name.c_str(), "...", (int)db.atoms.size(), (int)db.bonds.size(), 
        (int)db.coords.size(), (int)db.groups.size(), (int)db.sets.size(), numrigidheavy, 5, 0);
    os << buf;
    sprintf(buf, "M %9.4f %+10.3f %+10.3f %+10.3f %9.3f\n",
        db.charge, db.polarsolv, db.apolarsolv, db.solv, db.surface);
    os << buf;
    for (int i = 0; i < db.atoms.size(); i++) {
        auto atom = db.atoms[i];
        sprintf(buf, "A %3d %-4s %-5s %2d %2d %+9.4f %+10.3f %+10.3f %+10.3f %9.3f\n", i+1,
            atom.name.c_str(), atom.type.c_str(), atom.docktype, atom.color,
            atom.charge, atom.polarsolv, atom.apolarsolv, atom.solv, atom.surface);
        os << buf;
    }
    for (int i = 0; i < db.bonds.size(); i++) {
        auto bond = db.bonds[i];
        sprintf(buf, "B %3d %3d %3d %-2s\n", i+1, bond.atoma, bond.atomb, bond.bondtype.c_str());
        os << buf;
    }
    for (int i = 0; i < db.coords.size(); i++) {
        auto coord = db.coords[i];
        sprintf(buf, "X %9d %3d %6d %+9.4f %+9.4f %+9.4f\n", i+1, 
            coord.atomi+1, coord.groupi, coord.xyz.x, coord.xyz.y, coord.xyz.z);
        os << buf;
    }
    os << rigbuf; // output that data we were saving here
    for (int i = 0; i < db.groups.size(); i++) {
        auto group = db.groups[i];
        sprintf(buf, "C %6d %9d %9d\n", i+1, group.coords_start, group.coords_end);
        os << buf;
    }
    // antiquated set format
    for (int i = 0; i < db.sets.size(); i++) {
        auto set = db.sets[i];
        int numlines = (int)ceil(set.group_ids.size()/8.0);
        int lastlinelength = set.group_ids.size() % 8;
        sprintf(buf, "S %6d %6d %3d %1d %1d %+11.3f\n", i+1,
            numlines, (int)set.group_ids.size(), 0, 0, set.energy_1);
        os << buf;
        for (int j = 0; j < numlines; j++) {
            sprintf(buf, "S %6d %6d %1d", i+1, numlines, lastlinelength);
            os << buf;
            for (int k = 0; k < 8; k++) {
                int gidx = j*8+k;
                if (gidx >= set.group_ids.size())
                    break;
                sprintf(buf, " %6d", set.group_ids[gidx]);
                os << buf;
            }
            os << endl;
        }
    }
    os << "E" << endl; // done!
    return true;
}

// evil macro function for more prettiness and compactness
// times FUNC_CALL with std::chrono, continues if function returns false
#define TIME_FUNCTION_OR_CONTINUE(FUNC_CALL, TOTAL_TIME)                    \
{                                                                           \
    auto tstart = high_resolution_clock::now();                             \
    if (!FUNC_CALL) {                                                       \
        continue;                                                           \
    }                                                                       \
    auto tend = high_resolution_clock::now();                               \
    long dur = duration_cast<milliseconds>(tend-tstart).count();            \
    TOTAL_TIME += dur;                                                      \
    cout << dur << "ms elapsed for "#FUNC_CALL << endl;                     \
}            

void mol2db2_CDPKit(string fname) {

    using namespace chrono;
    
    MoleculeReader mol_file_reader(fname);
    BasicMolecule basemol;
    int moltot = 0;
    int maxtauts = 4;
    long taut_exec_tot = 0;
    long conf_exec_tot = 0;
    long solv_exec_tot = 0;
    long clst_exec_tot = 0;
    long writ_exec_tot = 0;
    auto tstart0 = high_resolution_clock::now();

    while (mol_file_reader.read(basemol)) {
        int tautid = 0;
        vector<BasicMolecule> tautomers;
        cout << "===============================" << endl;

        TIME_FUNCTION_OR_CONTINUE(enumtautomers(basemol, maxtauts, tautomers), taut_exec_tot);

        for (BasicMolecule mol : tautomers) {
            cout << "......" << getName(basemol) << "." << tautid << "::" << moltot << endl;
            db2 db;
            vector<db2atom> atom_solv_data;
            filebuf fb;
            db.name = getName(basemol)+"."+to_string(tautid);
            fb.open(db.name+".db2", ios::out);
            ostream os(&fb);

            TIME_FUNCTION_OR_CONTINUE(generateconformers(mol), conf_exec_tot)

            TIME_FUNCTION_OR_CONTINUE(calcsolvation(mol, atom_solv_data, db), solv_exec_tot)

            TIME_FUNCTION_OR_CONTINUE(generatedb2(mol, atom_solv_data, db), clst_exec_tot)

            TIME_FUNCTION_OR_CONTINUE(writedb2(db, os), writ_exec_tot);

            moltot++;
            tautid++;

            cout << db.sets.size() << " " << db.groups.size() << " " << db.coords.size() << endl;
            cout << " charge=" << db.charge << " psolve=" << db.polarsolv << " surf=" << db.surface << " apsolv=" << db.apolarsolv << " solv=" << db.solv << " minE=" << db.sets[0].energy_1 << endl;
        }
    }
    auto tend0 = high_resolution_clock::now();
    long exec_tot = duration_cast<milliseconds>(tend0-tstart0).count();
    long tracked_exec_tot = taut_exec_tot + conf_exec_tot + solv_exec_tot + clst_exec_tot;

    cout << "+++++++++++++++++++++++++++++++" << endl;
    cout << "evaluated " << moltot << " mols" << endl;
    cout << taut_exec_tot << "ms for tautomerization (n=" << maxtauts << ")" << endl;
    cout << conf_exec_tot << "ms for 3d conformer (n=100)" << endl;
    cout << solv_exec_tot << "ms for solvation" << endl;
    cout << clst_exec_tot << "ms for db2 clustering" << endl;
    cout << exec_tot - tracked_exec_tot << "ms for other" << endl;
    cout << exec_tot << "ms total" << endl;
}

int main(int argc, char** argv) {
    string fname = argv[1];
    mol2db2_CDPKit(fname);
}
