#include "CrystalStructure.h"
#include "taam_utilities.h"

//using namespace discamb;
using namespace std;

#include <fstream>


int main(int argc, char *argv[])
{
    
    if (argc != 2)
    {
        cerr << "expecting structure file\n";
        exit(0);
    }

    CrystalStructure crystalStructure(argv[1]);
    vector<vector<double> > adps;
    adps = crystalStructure.getAdps();
    vector<string> symbols;
    symbols = crystalStructure.getAtomSymbols();
    vector<vector<double> > xyz;
    xyz = crystalStructure.getPositions();
    vector<int> bonds;
    bonds = crystalStructure.getBonds();
    vector<string> names;
    names = crystalStructure.getAtomNames();

    vector< vector<vector<int> > > fragments{ {{0, 0, 0, 0}, {1, 0, 0, 0}, {4, 0, 0, 0}, {6, 0, 0, 0}, {8, 0, 0, 0}, {9, 0, 0, 0}, {7, -1, -1, 0}, {3, 0, -1, -1}, {5, -1, 0, -1}, {2, 0, 0, 0}}, {{0, 0, 0, 0}, {1, 0, 0, 0}, {4, 0, 0, 0}} };

    auto types = taam_utilities::find_atom_types_for_fragments(crystalStructure, "MATTS2021databank.txt", fragments);
    cout << types.size() << endl;
            
}


