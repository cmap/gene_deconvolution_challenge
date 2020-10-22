//
// Created by Wojciech 'Ardavel' Szałapski
//

#ifndef DPEAK_BIVARIATENORMALDISTRIBUTION_H
#define DPEAK_BIVARIATENORMALDISTRIBUTION_H

#include "data/types.h"

#include <array>
#include <vector>

class BivariateNormalDistribution
{
public:

    void initialize(
        const std::vector<std::array<real, 2>> &samples);

    real density(
        const std::array<real, 2> &value) const;

    bool isInitialized() const;

private:

    std::array<real, 2> m_standardDeviation;
    std::array<real, 2> m_mean;
    real m_correlation;

    bool m_initialized = false;
};

#endif //DPEAK_BIVARIATENORMALDISTRIBUTION_H
