#include "CrystalStructure.h"
#include "discamb/IO/structure_io.h"
#include "discamb/BasicChemistry/periodic_table.h"
#include "discamb/CrystalStructure/crystal_structure_utilities.h"
#include "discamb/CrystalStructure/UnitCellContent.h"
#include "discamb/MathUtilities/algebra3d.h"
#include "discamb/StructuralProperties/structural_properties.h"

using namespace std;
using namespace discamb;

CrystalStructure::CrystalStructure()
{
    
}

CrystalStructure::CrystalStructure(
    const std::string& fileName)
{
    Crystal crystal;
    structure_io::read_structure(fileName, crystal);
    mUnitCellContent.set(crystal);

    structural_properties::calcUnitCellConnectivity(mUnitCellContent, mUnitCellConnectivity, 0.4);
    
    //crystal_structure_utilities::ge

    vector<vector<UnitCellContent::AtomID> > molecules;
    vector< vector< pair< UnitCellContent::AtomID, UnitCellContent::AtomID > > > networkBonds;

    structural_properties::splitAsymmetricUnitIntoMolecules(mUnitCellContent, molecules, networkBonds, 0.4);
    mAtoms.clear();
    
    for (auto& molecule : molecules)
        mAtoms.insert(mAtoms.end(), molecule.begin(), molecule.end());

    for (auto const& atom : crystal.atoms)
    {
        int adpSize = atom.adp.size();
        mThermalEllispoids.resize(mThermalEllispoids.size() + 1);
        if (adpSize > 0)
            mThermalEllispoids.back().defined = true;
        if (adpSize == 1)
        {
            mThermalEllispoids.back().isotropic = true;
            mThermalEllispoids.back().u_iso = atom.adp[0];
        }
        else if (adpSize == 6)
        {
            auto& ellipsoid = mThermalEllispoids.back();
            ellipsoid.adpEigenvalues.resize(3);
            ellipsoid.axesDirections.resize(3);
            StructuralParametersConverter converter(crystal.unitCell);
            vector<double> adpCart(6);
            converter.convertADP(
                atom.adp,
                adpCart,
                structural_parameters_convention::AdpConvention::U_cif,
                structural_parameters_convention::AdpConvention::U_cart);

            Matrix3d u_cart(adpCart[0], adpCart[3], adpCart[4],
                adpCart[3], adpCart[1], adpCart[5],
                adpCart[4], adpCart[5], adpCart[2]);
            Vector3d v1, v2, v3, v_frac, v_frac_conv;
            double lambda1, lambda2, lambda3;
            Matrix3d rotation;

            algebra3d::eigensystemRealSymm(
                u_cart,
                ellipsoid.axesDirections[0],
                ellipsoid.axesDirections[1],
                ellipsoid.axesDirections[2],
                ellipsoid.adpEigenvalues[0],
                ellipsoid.adpEigenvalues[1],
                ellipsoid.adpEigenvalues[2]);
        }
    }

    vector<vector<double> > positions = getPositions();
    vector<Vector3d> r;
    for (auto& position : positions)
        r.push_back(position);
    vector<vector<int> > connectivity;
    vector<int> z;
    crystal_structure_utilities::atomicNumbers(crystal, mAtomicNumbers);
    

    for (auto& atom : mAtoms)
    {
        int atomInAsymmUnit = mUnitCellContent.indexOfSymmetryEquivalentAtomInCrystal(atom.atomIndex);
        z.push_back(mAtomicNumbers[atomInAsymmUnit]);
    }

    structural_properties::calculateConnectivity(r, z, connectivity);
    int nBonds = 0;
    for (auto& bonds : connectivity)
        nBonds += bonds.size();
    
    for (int i = 0; i < connectivity.size(); i++)
        for (int j = 0; j < connectivity[i].size(); j++)
            if (connectivity[i][j] < i)
            {
                mBonds.push_back(i);
                mBonds.push_back(connectivity[i][j]);
            }

}

std::vector<int> CrystalStructure::getBonds() const
{
    return mBonds;
}

std::vector<std::string> CrystalStructure::getAtomSymbols()
const
{
    vector<int> z;
    const Crystal& c = mUnitCellContent.getCrystal();
    crystal_structure_utilities::atomicNumbers(c, z);

    vector<string> symbols;
    for (auto const& atom : mAtoms)
    {
        int idx = mUnitCellContent.indexOfSymmetryEquivalentAtomInCrystal(atom.atomIndex);
        symbols.push_back(periodic_table::symbol(z[idx]));
    }
    return symbols;
}

vector<string> CrystalStructure::getAtomNames()
const
{
    vector<string> atomNames;
    Crystal const & crystal = mUnitCellContent.getCrystal();
    for (auto const& atom : crystal.atoms)
        atomNames.push_back(atom.label);
    return atomNames;
}

vector<vector<double> > CrystalStructure::getPositions()
const
{
    vector<vector<double> > positions;
    vector<double> r(3);
    UnitCellContent uc;
    
    for (auto const& atom : mAtoms)
    {
        //crystal_structure_utilities::atomPosition()
        Vector3d v = mUnitCellContent.getAtomPositionCart(atom);
        r[0] = v[0];
        r[1] = v[1];
        r[2] = v[2];
        positions.push_back(r);
    }
    return positions;

}

vector<vector<double> > CrystalStructure::getAdps()
const
{
    vector<vector<double> > adps;
//    double c05 = 1.538;
//    for (auto const& atom : mAtoms)
//    {
//        adps.resize(adps.size() + 1);
//        auto& adp = adps.back();
//        auto ellipsoid = getThermalEllispoid(atom);
//        if (!ellipsoid.defined)
//            continue;
//        if (ellipsoid.isotropic)
//            adp.push_back(sqrt(ellipsoid.u_iso) * c05);
//        else
//        {
//            for (int i = 0; i < 3; i++)
//                for (int j = 0; j < 3; j++)
//                    adp.push_back(ellipsoid.axesDirections[i][j]);
//            for (int i = 0; i < 3; i++)
//                adp.push_back(sqrt(ellipsoid.adpEigenvalues[i]) * c05);
//        }
//    }
    for (auto const& atom : mAtoms)
        adps.push_back(getThermalEllipsoid(atom.atomIndex));
    return adps;

}

ThermalEllispoid CrystalStructure::getThermalEllispoid(
    int atomInUnitCellIdx)
    //discamb::UnitCellContent::AtomID const& atomId)
    const
{
    ThermalEllispoid ellipsoid;

    int idx = mUnitCellContent.indexOfSymmetryEquivalentAtomInCrystal(atomInUnitCellIdx);// atomId.atomIndex);
    ellipsoid = mThermalEllispoids[idx];
    if (mThermalEllispoids[idx].adpEigenvalues.size() == 3)
    {
        auto const& crystal = mUnitCellContent.getCrystal();
        ellipsoid.axesDirections.resize(3);
        Vector3d v[3], v_frac, v_frac_conv;
        Matrix3d rotation;
        mUnitCellContent.getGeneratingOperation(atomInUnitCellIdx, 0).getRotation(rotation);

        for (int i = 0; i < 3; i++)
        {
            crystal.unitCell.cartesianToFractional(mThermalEllispoids[idx].axesDirections[i], v_frac);
            v_frac_conv = rotation * v_frac;
            crystal.unitCell.fractionalToCartesian(v_frac_conv, ellipsoid.axesDirections[i]);
        }
    }
    return ellipsoid;

    //out << " " << c05 * sqrt(lambda1) << " " << c05 * sqrt(lambda2) << " " << c05 * sqrt(lambda3) << "\n";

}

int CrystalStructure::atomicNumber(int indexInUnitCell)
{
    return mAtomicNumbers[mUnitCellContent.indexOfSymmetryEquivalentAtomInCrystal(indexInUnitCell)];
}

std::string CrystalStructure::getAtomLabel(
    int indexInUnitCell)
    const
{
    UnitCellContent::AtomID id(indexInUnitCell);
    string symmOp, label;
    mUnitCellContent.interpreteAtomID(id, label, symmOp);
    return label;
}

std::string CrystalStructure::getSymmetryOperationStr(
    const std::vector<int>& atomId)
    const
{
    UnitCellContent::AtomID id;
    id.atomIndex = atomId[0];
    id.unitCellPosition.set(atomId[1], atomId[2], atomId[3]);
    string symmOp, label;
    mUnitCellContent.interpreteAtomID(id, label, symmOp);
    return symmOp;
}



std::vector< std::vector<std::vector<int> > > CrystalStructure::getGrouppedAsymetricUnitBondedAtoms()
const
{

    vector<vector<UnitCellContent::AtomID> > molecules;
    vector< vector< pair< UnitCellContent::AtomID, UnitCellContent::AtomID > > > networkBonds;

    structural_properties::splitAsymmetricUnitIntoMolecules(mUnitCellContent, molecules, networkBonds, 0.4);
    vector<vector<vector<int> > > _molecules;
    vector<int> atomId(4);
    for (auto& molecule : molecules)
    {
        _molecules.resize(_molecules.size() + 1);
        for (auto const& atom : molecule)
        {
            atomId[0] = atom.atomIndex;
            atomId[1] = atom.unitCellPosition[0];
            atomId[2] = atom.unitCellPosition[1];
            atomId[3] = atom.unitCellPosition[2];
            _molecules.back().push_back(atomId);
        }
    }

    return _molecules;
}

std::vector<double> CrystalStructure::getAtomPositionCart(
    int indexInUnitCell,
    int nA,
    int nB,
    int nC)
    const
{
    Vector3d r = mUnitCellContent.getAtomPositionCart(UnitCellContent::AtomID(indexInUnitCell, Vector3i(nA, nB, nC)));
    return vector<double>{r[0], r[1], r[2]};
}


std::vector<double> CrystalStructure::getThermalEllipsoid(int indexInUnitCell)
const
{
    vector<double> result;
    double c05 = 1.538;

    auto ellipsoid = getThermalEllispoid(indexInUnitCell);
    if (!ellipsoid.defined)
        return result;
    if (ellipsoid.isotropic)
        result.push_back(sqrt(ellipsoid.u_iso) * c05);
    else
    {
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                result.push_back(ellipsoid.axesDirections[i][j]);
        for (int i = 0; i < 3; i++)
            result.push_back(sqrt(ellipsoid.adpEigenvalues[i]) * c05);
    }
    
    return result;

}

std::vector<std::pair<int, int> > CrystalStructure::getBonds(
    const vector<vector<int> >& atomIds)
    const
{
    vector<pair<int, int> > bonds;
    map<int, vector<int> > unitCellIdxAtoms;
    vector<Vector3i> unitCellShifts(atomIds.size());

    for (int i = 0; i < atomIds.size(); i++)
    {
        unitCellIdxAtoms[atomIds[i][0]].push_back(i);
        unitCellShifts[i].set(atomIds[i][1], atomIds[i][2], atomIds[i][3]);
    }

    for (int i = 0; i < atomIds.size(); i++)
    {
        int atomIdxInUnitCell = atomIds[i][0];
        Vector3i atomUnitCell(atomIds[i][1], atomIds[i][2], atomIds[i][3]);

        for (int j = 0; j < mUnitCellConnectivity[atomIdxInUnitCell].size(); j++)
        {
            auto neighbourId = mUnitCellConnectivity[atomIdxInUnitCell][j];
            neighbourId.unitCellPosition += atomUnitCell;
            for (int k = 0; k < unitCellIdxAtoms[neighbourId.atomIndex].size(); k++)
            {
                int candidateIdx = unitCellIdxAtoms[neighbourId.atomIndex][k];
                if(candidateIdx < i) // to avoid including bond twice
                    if (neighbourId.unitCellPosition == unitCellShifts[candidateIdx])
                    {
                        bonds.push_back({ i, candidateIdx });
                    }
            }
        }
    }
    return bonds;
}
