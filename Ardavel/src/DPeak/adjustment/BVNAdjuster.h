//
// Created by Wojciech 'Ardavel' Szałapski
//

#ifndef DPEAK_BVNADJUSTER_H
#define DPEAK_BVNADJUSTER_H

#include "data/DeconvolutionResult.h"
#include "data/UnsureResults.h"
#include "math/BivariateNormalDistribution.h"

#include <vector>

class BVNAdjuster
{
public:

    std::vector<UnsureResults> adjust(
        std::vector<DeconvolutionResult>& deconvolutionResults,
        std::vector<std::array<BivariateNormalDistribution, 2>>& distributionsPerBarcode,
        int maximumNumberOfIterations) const;
};

#endif //DPEAK_BVNADJUSTER_H
