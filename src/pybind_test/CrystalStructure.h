#pragma once

#include "discamb/CrystalStructure/UnitCellContent.h"
#include <string>
#include <vector>

struct ThermalEllispoid {
    std::vector<discamb::Vector3d> axesDirections;
    std::vector<double> adpEigenvalues;
    double u_iso = 0.0;
    bool isotropic = false;
    bool defined = false;
};

// commnet
class CrystalStructure
{
public:
    CrystalStructure();
    CrystalStructure(const std::string &fileName);
    
    std::vector<std::string> getAtomNames() const;
    std::vector<std::string> getAtomSymbols() const;
    std::vector<std::vector<double> > getPositions() const;
    std::vector<std::vector<double> > getAdps() const;
    std::vector<int> getBonds() const;
    std::vector<std::pair<int,int> > getBonds(const std::vector<std::vector<int> > &atomIds) const;
    std::vector<std::vector<int> > getAsymetricUnitBondedAtoms() const;
    std::vector< std::vector<std::vector<int> > > getGrouppedAsymetricUnitBondedAtoms() const;
    std::vector<std::vector<int> > getAsymetricUnitAtoms() const;
    std::vector< std::vector<std::vector<int> > > getGrouppedAsymetricUnitAtoms() const;
    std::string getAtomLabel(int indexInUnitCell) const;
    std::vector<double> getThermalEllipsoid(int indexInUnitCell) const;
    std::vector<double> getAtomPositionCart(int indexInUnitCell, int nA, int nB, int nC) const;
    //std::string getSymmetryOperationStr(int indexInUnitCell) const;
    std::string getSymmetryOperationStr(const std::vector<int> &atomId) const;
    int atomicNumber(int indexInUnitCell);

    
    
private:
    discamb::UnitCellContent mUnitCellContent;
    std::vector<discamb::UnitCellContent::AtomID> mAtoms;
    std::vector<int> mAtomicNumbers;
    ThermalEllispoid getThermalEllispoid(discamb::UnitCellContent::AtomID const& atomId) const;
    ThermalEllispoid getThermalEllispoid(int idx) const;
    std::vector<ThermalEllispoid> mThermalEllispoids;
    std::vector<int> mBonds;
    std::vector<std::vector<discamb::UnitCellContent::AtomID> > mUnitCellConnectivity;
    
};

