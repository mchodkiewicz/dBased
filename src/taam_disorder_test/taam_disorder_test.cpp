#include "discamb/BasicChemistry/periodic_table.h"
#include "discamb/BasicUtilities/on_error.h"
#include "discamb/BasicUtilities/parse_cmd.h"
#include "discamb/BasicUtilities/string_utilities.h"
#include "discamb/CrystalStructure/crystal_structure_utilities.h"
#include "discamb/IO/hkl_io.h"
#include "discamb/IO/structure_io.h"
#include "discamb/IO/tsc_io.h"
#include "discamb/Scattering/SfCalculator.h"
#include "json.hpp"


using namespace discamb;
using namespace std;

#include <fstream>
#include <iomanip>

void mix_tsc(
    const string& config)
{
    Crystal crystal;

    // read config
    ifstream in(config);
    string line, structureFile;
    vector<string> words;
    in >> structureFile;
    structure_io::read_structure(structureFile, crystal);
    getline(in, line);

    map<string, int> atomLabelToIndexInCrystal;

    vector<string> crystalAtomLabels;
    int idx = 0;
    for (auto& atom : crystal.atoms)
    {
        crystalAtomLabels.push_back(atom.label);
        atomLabelToIndexInCrystal[atom.label] = idx;
        idx++;
    }

    vector<string> atomLabels, tsc_files;
    vector<vector<pair<string, double> > > tsc_weights;

    while (getline(in, line))
    {
        string_utilities::split(line, words);
        if (words.size() == 0)
            continue;
        if (words.size() == 1)
        {
            tsc_files.push_back(words[0]);
            tsc_weights.resize(tsc_weights.size() + 1);
        }
        if (words.size() == 2)
        {
            double wght;
            if (isdigit(words[1][0]))
                wght = stod(words[1]);
            else
                wght = crystal.atoms[atomLabelToIndexInCrystal[words[1]]].occupancy;
            tsc_weights.back().push_back({ words[0], wght });
        }
    }


    vector<Vector3i> hkls;
    vector< vector<complex<double> > > tsc_ff, total_ff;
    vector<double> weights_sum(crystal.atoms.size(), 0.0);
    for (int i = 0; i < tsc_files.size(); i++)
    {
        tsc_io::read_tsc(tsc_files[i], atomLabels, hkls, tsc_ff);

        if (i == 0)
        {
            total_ff.resize(hkls.size());
            for (int j = 0; j < hkls.size(); j++)
                total_ff[j].resize(crystal.atoms.size(), complex<double>(0.0, 0.0));
        }

        double other_atom_weights = 0.0;
        map<string, double> atomWeightInTsc;
        for (auto& p : tsc_weights[i])
            if (p.first == "other")
                other_atom_weights = p.second;
        for (auto& label : atomLabels)
            atomWeightInTsc[label] = other_atom_weights;
        for (auto& p : tsc_weights[i])
            if (p.first != "other")
                atomWeightInTsc[p.first] = p.second;

        for (int idxInTsc = 0; idxInTsc < atomLabels.size(); idxInTsc++)
        {
            int idxInCrystal = atomLabelToIndexInCrystal[atomLabels[idxInTsc]];
            double atom_weight = atomWeightInTsc[atomLabels[idxInTsc]];
            weights_sum[idxInCrystal] += atom_weight;
            for (int hklIdx = 0; hklIdx < hkls.size(); hklIdx++)
                total_ff[hklIdx][idxInCrystal] += tsc_ff[hklIdx][idxInTsc] * atom_weight;
        }
    }

    for (int idxInCrystal = 0; idxInCrystal < crystal.atoms.size(); idxInCrystal++)
        for (int hklIdx = 0; hklIdx < hkls.size(); hklIdx++)
            if (weights_sum[idxInCrystal] > 0.0)
                total_ff[hklIdx][idxInCrystal] /= weights_sum[idxInCrystal];
    tsc_io::write_tsc("total_ff.tsc", crystalAtomLabels, hkls, total_ff);
}

void test_next_gen_taam(
    const string& structureFile,
    const string& hklFile)
{
    Crystal crystal;
    structure_io::read_structure(structureFile, crystal);

    nlohmann::json json_data;
    std::ifstream json_file("aspher.json");
    json_file >> json_data;
    /*
    "macromolecular structural information": {
    "atom data": [
        {
            "altloc": "",
            "bonding": [
                1,
                8,
                9,
                10
            ],
            "chain_id": "B",
            "name": "N",
            "plane": null,
            "r_idx": 1,
            "r_name": "LEU"
        },
    */
    /*
    auto& atom_list = json_data["macromolecular structural information"]["atom data"];
    int atomIdx = 0;
    for (auto& atom : atom_list)
    {
        string name = atom["name"];
        string chain_id = atom["chain_id"];
        string altloc = atom["altloc"];
        int r_idx = atom["r_idx"];
        string r_name = atom["r_name"];
        vector<int> bonding;
        for (auto& b : atom["bonding"])
            bonding.push_back(b);
        cout << "atom " << atomIdx << " " << name << " " << chain_id << " " << altloc << " " << r_idx << " " << r_name << "\n";
        cout << "bonding with: ";
        for (auto& b : bonding)
            cout << b << " ";
        cout << "\n";
        atomIdx++;
    }
    */
    vector<Vector3i> hkl;
    hkl_io::readHklIndices(hklFile, hkl);

    SfCalculator* sf_calc = SfCalculator::create(crystal, string("aspher.json"));
    vector<complex<double> > sf;
    sf_calc->calculateStructureFactors(crystal.atoms, hkl, sf);
    int i, n = hkl.size();
    ofstream out("sf");
    for (i = 0; i < n; i++)
        out << hkl[i] << " " << sf[i] << "\n";
    out.close();
}


int main(int argc, char *argv[])
{
    
    try{
        vector<string> arguments, options;
        parse_cmd::get_args_and_options(argc, argv, arguments, options);


        if (parse_cmd::hasOption(options, "-calc_sf"))
        {
            if (arguments.size() < 2)
                on_error::throwException("expected structure and hkl file\n", __FILE__, __LINE__);

            test_next_gen_taam(arguments[0], arguments[1]);
            return 0;
        }

        if (parse_cmd::hasOption(options, "-mix_tsc"))
        {
            if (arguments.size() != 1)
                on_error::throwException("expected config file for mix_tsc\n", __FILE__, __LINE__);
            mix_tsc(arguments[0]);
            return 0;
        }

        cout << "Usage: \n";
        cout << "taam_disorder_test -calc_sf structure_file hkl_file\n";
        cout << "taam_disorder_test -mix_tsc config_file\n";
        return 0;
    }
    catch (exception &e)
    {
        cout << e.what() << endl;
    }
}


