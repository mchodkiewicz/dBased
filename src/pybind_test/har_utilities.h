#include "discamb/Scattering/HirshfeldAtomModelSettings.h"

namespace har_utilities {

    struct Representatives {

        int substructureIdx;
        int atomIdx;
        double weight;
    };

    void find_default_representatives(const std::vector<std::vector<int> >& fragments);

}