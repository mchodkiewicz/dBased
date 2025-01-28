#pragma once

#include "discamb/CrystalStructure/UnitCellContent.h"
#include <string>
#include <vector>

struct ThermalEllispoid {
    std::vector<discamb::Vector3d> axesDirections;
    std::vector<double> adpEigenvalues;
    double u_iso;
    bool isotropic;
};

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
    
private:
    discamb::UnitCellContent mUnitCellContent;
    std::vector<discamb::UnitCellContent::AtomID> mAtoms;
    ThermalEllispoid getThermalEllispoid(discamb::UnitCellContent::AtomID const& atomId) const;
    std::vector<ThermalEllispoid> mThermalEllispoids;
    std::vector<int> mBonds;
    
};

