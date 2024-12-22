#include "discamb/BasicChemistry/periodic_table.h"
#include "discamb/BasicUtilities/on_error.h"
#include "discamb/CrystalStructure/crystal_structure_utilities.h"
#include "discamb/IO/structure_io.h"
#include "discamb/StructuralProperties/structural_properties.h"

using namespace discamb;
using namespace std;

#include <fstream>


int main(int argc, char *argv[])
{
    
    try{

        if (argc != 2)
            on_error::throwException("expected structure file on input", __FILE__, __LINE__);

        Crystal crystal;
        UnitCellContent ucContent;

        structure_io::read_structure(argv[1], crystal);
        
        ucContent.set(crystal);

        vector<int> atomicNumbers;
        crystal_structure_utilities::atomicNumbers(crystal, atomicNumbers);

        map<string, int> atomLabel2Idx;
        for (int i = 0; i < crystal.atoms.size(); i++)
            atomLabel2Idx[crystal.atoms[i].label] = i;

        vector<vector<UnitCellContent::AtomID> > molecules;
        vector< vector< pair< UnitCellContent::AtomID, UnitCellContent::AtomID > > > networkBonds;

        structural_properties::splitAsymmetricUnitIntoMolecules(ucContent, molecules, networkBonds, 0.4);

        ofstream out("out.txt");
        
        
        for (auto& molecule : molecules)
        {
            for (auto& atomId : molecule)
            {
                string label, symmOpStr;
                ucContent.interpreteAtomID(atomId, label, symmOpStr);
                int atomIdx = atomLabel2Idx[label];
                out << periodic_table::symbol(atomicNumbers[atomIdx]);
                SpaceGroupOperation symmOp(symmOpStr);
                Vector3d positionCart, positionFrac;
                symmOp.apply(crystal.atoms[atomIdx].coordinates, positionFrac);
                crystal.unitCell.fractionalToCartesian(positionFrac, positionCart);

            }
        }
    }
    catch (exception &e)
    {
        cout << e.what() << endl;
    }
}


