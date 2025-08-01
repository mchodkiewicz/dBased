#include "har_utilities.h"
#include "discamb/Scattering/gar_utilities.h"

#include <iostream>

using namespace std;

namespace har_utilities {

    std::vector<Representative>  find_default_representatives(
        const CrystalStructure& crystal_structure,
        const std::vector<std::vector<std::vector<int> > >& fragments)
        {
            std::vector<Representative> rerpresenataives;
            //rerpresenataives.push_back({ 1, 2, 0.23 });
            //rerpresenataives.push_back({ 5, 4, 0.12 });
            //return rerpresenataives;

            int nFragments = fragments.size();
            vector<vector<pair<string, string> > > subsystemAtoms(nFragments);
            for (int fragIdx = 0; fragIdx < nFragments; fragIdx++)
            {
                for (auto const& atom : fragments[fragIdx])
                {
                    discamb::UnitCellContent::AtomID atomId(atom[0], discamb::Vector3i(atom[1], atom[2], atom[3]));
                    string label, symmOpStr;
                    crystal_structure.getUnitCellContent().interpreteAtomID(atomId, label, symmOpStr);
                    subsystemAtoms[fragIdx].push_back({ label, symmOpStr });
                    cout << "har_utilities.cpp line " << __LINE__ << " " << subsystemAtoms[fragIdx].back().first << " " << subsystemAtoms[fragIdx].back().second << endl;
                }
            }

            vector<vector<discamb::AtomRepresentativeInfo> > reps;
            discamb::gar_utilities::findDefaultRepresentatives(
                crystal_structure.getCrystal(), 
                subsystemAtoms,
                reps);
            for (int atomIdx=0; atomIdx< reps.size(); atomIdx++)
            {
                for (auto& rep : reps[atomIdx])
                {
                    Representative representative;
                    representative.atomIdx = atomIdx;
                    representative.substructureIdx = rep.fragmentIdx;
                    representative.weight = rep.fixedWeightValue;
                    cout << "representative.weight = " << representative.weight << endl;
                    rerpresenataives.push_back(representative);
                }
            }
            return rerpresenataives;
        }

}
