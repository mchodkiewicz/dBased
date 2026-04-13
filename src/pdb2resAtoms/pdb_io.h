#pragma once

#include "discamb/MathUtilities/Vector3.h"
#include "discamb/CrystalStructure/UnitCell.h"

#include <vector>
#include <string>

namespace discamb {

    /**
    * \addtogroup IO IO
    * @{
    */


    namespace pdb_io {
        
        struct PdbAtom {
            std::string name;
            int atomicNumber;
            Vector3d position;
            int serialNumber = 0;
            char alternateLocation = ' ';
            std::vector<double> temperatureFactor;
            double occupancy = 1.0;
            int charge = 0;
        };

        struct PdbResidue {
            std::vector<PdbAtom> atoms;
            std::string name;
        };

        struct PdbChain {
            std::vector<PdbResidue> residues;
            std::string name;
        };


        struct PdbData {
            std::vector<PdbChain> chains;
            UnitCell unitCell;
        };

        void read(const std::string& fileName, PdbData& data);
        void write(const std::string& fileName, const std::vector<int>& atomicNumber, const std::vector<Vector3d>& position);
        void write(const std::string& fileName, const  std::vector<std::string>& symbol, const std::vector<Vector3d>& position);
        /**
        numeration starts with 0
        */
        void standardAminoAcidResidueAtoms(const PdbData& data, std::vector<int> &ids);
        //[chain][idx] residue in chain which is standard amino acid
        void standardAminoAcidResidues(const PdbData& data, std::vector< std::vector<int> > & ids);

        bool is_standard_aminoacid(const std::string& code);
    }
    /**@}*/
}

