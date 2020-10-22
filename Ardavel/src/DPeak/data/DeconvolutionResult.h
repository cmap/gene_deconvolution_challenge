//
// Created by Wojciech 'Ardavel' Szałapski
//

#ifndef DPEAK_DECONVOLUTIONRESULT_H
#define DPEAK_DECONVOLUTIONRESULT_H

#include "data/Cluster.h"

#include <map>

using DeconvolutionResult = std::map<int, std::array<ClusterWithVariance, 2>>;

#endif //DPEAK_DECONVOLUTIONRESULT_H