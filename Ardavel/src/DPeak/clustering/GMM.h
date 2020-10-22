//
// Created by Wojciech 'Ardavel' Szałapski
//

#ifndef DPEAK_GMM_H
#define DPEAK_GMM_H

#include "math/NormallyDistributedClusterFit.h"

#include <array>
#include <vector>

class GMM
{
public:

    struct ModelFit
    {
        std::array<NormallyDistributedClusterFit, 2> m_clusters;
        bool m_valid = false;
    };

    ModelFit clusterize(
        const std::vector<real>& samples,
        const std::array<real, 2>& initialCenters,
        int maxNumberOfIterations,
        real minMoveOfMeans,
        std::vector<std::array<real, 2>>& sampleProbability) const;
};

#endif //DPEAK_GMM_H
