#include "discamb/BasicChemistry/periodic_table.h"
#include "discamb/BasicUtilities/on_error.h"
#include "discamb/BasicUtilities/parse_cmd.h"
#include "discamb/CrystalStructure/crystal_structure_utilities.h"
#include "discamb/IO/structure_io.h"
#include "pdb_io.h"

using namespace discamb;
using namespace std;

#include <fstream>
#include <iomanip>


int main(int argc, char *argv[])
{
    
    try{
        vector<string> arguments, options;
        parse_cmd::get_args_and_options(argc, argv, arguments, options);

        if (arguments.size() != 3)
            on_error::throwException("expected pdb input_res output_res\n", __FILE__, __LINE__);
        pdb_io::PdbData data;
        pdb_io::read(arguments[0], data);
        Crystal crystal;
        structure_io::read_structure(arguments[1], crystal);
        crystal.atoms.clear();
        crystal.unitCell = data.unitCell;
        bool split_by_altloc = parse_cmd::hasOption(options, "-split");
        int nAtoms = 0;
        set<char> altLocs;
        map<char, vector<int> > altLocAtomIdx;
        for (auto& chain : data.chains)
            for (auto& residue : chain.residues)
                for (auto& atom : residue.atoms)
                {
                    altLocAtomIdx[atom.alternateLocation].push_back(nAtoms);
                    AtomInCrystal atomInCrystal;
                    nAtoms++;
                    atomInCrystal.adp = atom.temperatureFactor;
                    data.unitCell.cartesianToFractional(atom.position, atomInCrystal.coordinates);
                    atomInCrystal.label = periodic_table::symbol(atom.atomicNumber) + to_string(nAtoms);
                    if (atom.alternateLocation != ' ')
                    {
                        atomInCrystal.label += atom.alternateLocation;
                        altLocs.insert(atom.alternateLocation);
                    }
                    atomInCrystal.occupancy = atom.occupancy;
                    atomInCrystal.type = periodic_table::symbol(atom.atomicNumber);
                    crystal.atoms.push_back(atomInCrystal);
                    cout << atomInCrystal.label << " " << atom.serialNumber << endl;
                }
        crystal_structure_utilities::set_atoms_multiplicity(crystal);
        structure_io::write_structure(arguments[2], crystal);
        if (split_by_altloc)
        {
            Crystal splitCrystal;
            splitCrystal = crystal;
            splitCrystal.atoms.clear();
            vector<AtomInCrystal> orderedAtoms;
            for(int idx: altLocAtomIdx[' '])
                orderedAtoms.push_back(crystal.atoms[idx]);
            for (char altLoc : altLocs)
                if (altLoc != ' ')
                {
                    Crystal splitCrystal;
                    splitCrystal = crystal;
                    splitCrystal.atoms = orderedAtoms;

                    for (int idx : altLocAtomIdx[altLoc])
                        splitCrystal.atoms.push_back(crystal.atoms[idx]);
                    string name_core = arguments[2].substr(0, arguments[2].find_last_of('.'));
                    string file_extension = arguments[2].substr(arguments[2].find_last_of('.'));
                    structure_io::write_structure(name_core + "_" + altLoc + file_extension, splitCrystal);
                }
        }
    }
    catch (exception &e)
    {
        cout << e.what() << endl;
    }
}


