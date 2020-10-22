//
// Created by Wojciech 'Ardavel' Szałapski
//

#ifndef DPEAK_PIPELINE_H
#define DPEAK_PIPELINE_H

#include "math/NormallyDistributedClusterFit.h"

#include <vector>

class Pipeline
{
public:

    void run();

private:

    int findIndexOfFirstSampleAssignedToSecondCluster(
        const std::vector<real>& samples,
        const NormallyDistributedClusterFit& cluster1,
        const NormallyDistributedClusterFit& cluster2,
        const real& thresholdValue) const;
};

#endif //DPEAK_PIPELINE_H
