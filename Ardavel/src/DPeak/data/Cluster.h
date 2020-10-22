//
// Created by Wojciech 'Ardavel' Szałapski
//

#ifndef DPEAK_CLUSTER_H
#define DPEAK_CLUSTER_H

#include "data/types.h"

#include <array>

typedef std::array<real, 2> Cluster;

struct ClusterWithVariance
{
    Cluster m_cluster;
    real m_variance;
};

#endif //DPEAK_CLUSTER_H
