//
// Created by Wojciech 'Ardavel' Szałapski
//

#ifndef DPEAK_UNSURERESULTS_H
#define DPEAK_UNSURERESULTS_H

#include "data/types.h"

#include <array>
#include <vector>

struct UnsureResults
{
    std::vector<int> m_indices;

    std::array<real, 2> m_means;
    std::array<real, 2> m_stdDevs;
    std::array<real, 2> m_medianStdDevs;
};

#endif //DPEAK_UNSURERESULTS_H
