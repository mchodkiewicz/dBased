#include "pdb_io.h"
#include "discamb/BasicChemistry/periodic_table.h"
#include "discamb/MathUtilities/math_utilities.h"
#include "gemmi/pdb.hpp"

#include <fstream>
#include <iomanip>

using namespace std;

namespace {

    vector<string> standardAminoacidCodes({ "GLY", "ALA", "VAL", "LEU", "ILE", "PRO", "SER", "THR", "ASN", "GLN", "CYS", "MET",
    "PHE", "TYR", "TRP", "ASP", "GLU", "HIS", "LYS", "ARG", "SEC", "PYL"});

}

namespace discamb {

    namespace pdb_io {

        void write(
            const std::string& fileName, 
            const  std::vector<std::string>& symbols, 
            const std::vector<Vector3d>& positions)
        {
            int nAtoms = symbols.size();
            ofstream out(fileName);
            for (int i = 0; i < nAtoms; i++)
            {
                out << "HETATM" << setw(5) << i + 1 << setw(4) << symbols[i]
                    << "   LIG     1  ";
                for (int j = 0; j < 3; j++)                   
                    out << setw(8) << setprecision(3) << fixed << positions[i][j];
                
                out << "                       ";
                if (symbols[i].size() == 1) 
                    out << " ";
                out << symbols[i] << "\n";
            }
            out.close();
        }

        void write(
            const std::string& fileName, 
            const std::vector<int>& atomicNumber, 
            const std::vector<Vector3d>& position)
        {

            vector<string> symbols;
            for (auto const& z : atomicNumber)
                symbols.push_back(periodic_table::symbol(z));
            write(fileName, symbols, position);
        }


        void read(const std::string& fileName, PdbData& data)
        {
            data.chains.clear();
            gemmi::Structure structure = gemmi::read_pdb_file(fileName);
            data.unitCell.set(structure.cell.a, structure.cell.b, structure.cell.c, structure.cell.alpha, structure.cell.beta, structure.cell.gamma);
            int nChains = structure.models[0].chains.size();
            PdbChain chain;
            PdbResidue residue;
            PdbAtom atom;
            for (int chainIdx = 0; chainIdx < nChains; chainIdx++)
            {
                chain.residues.clear();
                chain.name = structure.models[0].chains[chainIdx].name;
                int nResidues = structure.models[0].chains[chainIdx].residues.size();
                for (int residueIdx = 0; residueIdx < nResidues; residueIdx++)
                {
                    residue.name = structure.models[0].chains[chainIdx].residues[residueIdx].name;
                    residue.atoms.clear();
                    int nAtoms = structure.models[0].chains[chainIdx].residues[residueIdx].atoms.size();
                    for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
                    {
                        auto const& pdbAtom = structure.models[0].chains[chainIdx].residues[residueIdx].atoms[atomIdx];
                        atom.alternateLocation = pdbAtom.altloc;
                        if(atom.alternateLocation == '\0')
                            atom.alternateLocation = ' ';
                        atom.atomicNumber = pdbAtom.element.atomic_number();
                        atom.serialNumber = pdbAtom.serial;
                        atom.name = pdbAtom.name;
                        atom.charge = pdbAtom.charge;
                        atom.occupancy = pdbAtom.occ;
                        atom.position.set(pdbAtom.pos.x, pdbAtom.pos.y, pdbAtom.pos.z);
                        if(pdbAtom.aniso.nonzero())
                            atom.temperatureFactor = { pdbAtom.aniso.u11, pdbAtom.aniso.u22, pdbAtom.aniso.u33, pdbAtom.aniso.u12, pdbAtom.aniso.u13, pdbAtom.aniso.u23 };
                        else
                            atom.temperatureFactor = { pdbAtom.b_iso/(8.0*M_PI*M_PI)};
                        residue.atoms.push_back(atom);
                    }
                    chain.residues.push_back(residue);
                }
                data.chains.push_back(chain);
            }
        }

        bool is_standard_aminoacid(
            const std::string& code)
        {
            return (find(standardAminoacidCodes.begin(), standardAminoacidCodes.end(), code) != standardAminoacidCodes.end());
        }

        //[chain][idx] residue in chain which is standard amino acid
        void standardAminoAcidResidues(
            const PdbData& data,
            std::vector< std::vector<int> >& ids)
        {
            int residueIdx, nResidues, nChains = data.chains.size();

            ids.resize(data.chains.size());

            for (int chainIdx = 0; chainIdx < nChains; chainIdx++)
            {
                nResidues = data.chains[chainIdx].residues.size();
                for (residueIdx = 0; residueIdx < nResidues; residueIdx++)
                    if (is_standard_aminoacid(data.chains[chainIdx].residues[residueIdx].name))
                        ids[chainIdx].push_back(residueIdx);

            }
        }


        void standardAminoAcidResidueAtoms(
            const PdbData& data,
            std::vector<int>& ids)
        {
            for(auto const & chain: data.chains)
                for (auto const& residue : chain.residues)
                {
                    if (is_standard_aminoacid(residue.name))
                        for (auto const& atom : residue.atoms)
                            ids.push_back(atom.serialNumber - 1);
                }
        }
    }

}

