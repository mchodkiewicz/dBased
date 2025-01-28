#include "discamb/BasicChemistry/periodic_table.h"
#include "discamb/BasicUtilities/on_error.h"
#include "discamb/CrystalStructure/crystal_structure_utilities.h"
#include "discamb/CrystalStructure/StructuralParametersConverter.h"
#include "discamb/IO/structure_io.h"
#include "discamb/MathUtilities/algebra3d.h"
#include "discamb/StructuralProperties/structural_properties.h"

using namespace discamb;
using namespace std;

#include <fstream>
#include <iomanip>

int colors[] =
{
      0xFFFFFFFF, // H  1
      0xFFD9FFFF, // He 2
      0xFFCC80FF, // Li 3
      0xFFC2FF00, // Be 4
      0xFFFFB5B5, // B  5
      0xFF909090, // C  6 - changed from ghemical
      0xFF3050F8, // N  7 - changed from ghemical
      0xFFFF0D0D, // O  8
      0xFF90E050, // F  9 - changed from ghemical
      0xFFB3E3F5, // Ne 10
      0xFFAB5CF2, // Na 11
      0xFF8AFF00, // Mg 12
      0xFFBFA6A6, // Al 13
      0xFFF0C8A0, // Si 14 - changed from ghemical
      0xFFFF8000, // P  15
      0xFFFFFF30, // S  16
      0xFF1FF01F, // Cl 17
      0xFF80D1E3, // Ar 18
      0xFF8F40D4, // K  19
      0xFF3DFF00, // Ca 20
      0xFFE6E6E6, // Sc 21
      0xFFBFC2C7, // Ti 22
      0xFFA6A6AB, // V  23
      0xFF8A99C7, // Cr 24
      0xFF9C7AC7, // Mn 25
      0xFFE06633, // Fe 26 - changed from ghemical
      0xFFF090A0, // Co 27 - changed from ghemical
      0xFF50D050, // Ni 28 - changed from ghemical
      0xFFC88033, // Cu 29 - changed from ghemical
      0xFF7D80B0, // Zn 30
      0xFFC28F8F, // Ga 31
      0xFF668F8F, // Ge 32
      0xFFBD80E3, // As 33
      0xFFFFA100, // Se 34
      0xFFA62929, // Br 35
      0xFF5CB8D1, // Kr 36
      0xFF702EB0, // Rb 37
      0xFF00FF00, // Sr 38
      0xFF94FFFF, // Y  39
      0xFF94E0E0, // Zr 40
      0xFF73C2C9, // Nb 41
      0xFF54B5B5, // Mo 42
      0xFF3B9E9E, // Tc 43
      0xFF248F8F, // Ru 44
      0xFF0A7D8C, // Rh 45
      0xFF006985, // Pd 46
      0xFFC0C0C0, // Ag 47 - changed from ghemical
      0xFFFFD98F, // Cd 48
      0xFFA67573, // In 49
      0xFF668080, // Sn 50
      0xFF9E63B5, // Sb 51
      0xFFD47A00, // Te 52
      0xFF940094, // I  53
      0xFF429EB0, // Xe 54
      0xFF57178F, // Cs 55
      0xFF00C900, // Ba 56
      0xFF70D4FF, // La 57
      0xFFFFFFC7, // Ce 58
      0xFFD9FFC7, // Pr 59
      0xFFC7FFC7, // Nd 60
      0xFFA3FFC7, // Pm 61
      0xFF8FFFC7, // Sm 62
      0xFF61FFC7, // Eu 63
      0xFF45FFC7, // Gd 64
      0xFF30FFC7, // Tb 65
      0xFF1FFFC7, // Dy 66
      0xFF00FF9C, // Ho 67
      0xFF00E675, // Er 68
      0xFF00D452, // Tm 69
      0xFF00BF38, // Yb 70
      0xFF00AB24, // Lu 71
      0xFF4DC2FF, // Hf 72
      0xFF4DA6FF, // Ta 73
      0xFF2194D6, // W  74
      0xFF267DAB, // Re 75
      0xFF266696, // Os 76
      0xFF175487, // Ir 77
      0xFFD0D0E0, // Pt 78 - changed from ghemical
      0xFFFFD123, // Au 79 - changed from ghemical
      0xFFB8B8D0, // Hg 80 - changed from ghemical
      0xFFA6544D, // Tl 81
      0xFF575961, // Pb 82
      0xFF9E4FB5, // Bi 83
      0xFFAB5C00, // Po 84
      0xFF754F45, // At 85
      0xFF428296, // Rn 86
      0xFF420066, // Fr 87
      0xFF007D00, // Ra 88
      0xFF70ABFA, // Ac 89
      0xFF00BAFF, // Th 90
      0xFF00A1FF, // Pa 91
      0xFF008FFF, // U  92
      0xFF0080FF, // Np 93
      0xFF006BFF, // Pu 94
      0xFF545CF2, // Am 95
      0xFF785CE3, // Cm 96
      0xFF8A4FE3, // Bk 97
      0xFFA136D4, // Cf 98
      0xFFB31FD4, // Es 99
      0xFFB31FBA, // Fm 100
      0xFFB30DA6, // Md 101
      0xFFBD0D87, // No 102
      0xFFC70066, // Lr 103
      0xFFCC0059, // Rf 104
      0xFFD1004F, // Db 105
      0xFFD90045, // Sg 106
      0xFFE00038, // Bh 107
      0xFFE6002E, // Hs 108
      0xFFEB0026  // Mt 109

};

int main(int argc, char *argv[])
{
    
    try{
        int i = 0xF242;
        int k = (i & 0xFF00);
        int l = k >> 8;
        cout << std::hex << l;
        ofstream out("out.xyz");
        out << "118\nperiodic table\n";
        for (int z = 1; z < 119; z++)
            out << periodic_table::symbol(z) << "  " << 2 * z << "  0.0   0.0\n";
        out.close();
        out.open("colors.txt");
        for (int z = 1; z <= 109; z++)
        {
            int r = (colors[z - 1] & 0x00FF0000) >> 16;
            int g = (colors[z - 1] & 0x0000FF00) >> 8;
            int b = (colors[z - 1] & 0x000000FF);
            out << "\"" << periodic_table::symbol(z) << "\": [" << r << ", " << g << ", " << b << "],\n";
        }
    }
    catch (exception &e)
    {
        cout << e.what() << endl;
    }
}


