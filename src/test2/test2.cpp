#include "CrystalStructure.h"
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
            
}


