//
// Created by Wojciech 'Ardavel' Szałapski
//

#ifndef DPEAK_BASIC_STATISTICS_H
#define DPEAK_BASIC_STATISTICS_H

#include "data/types.h"

#include <cmath>

template<typename Iterator>
real mean(
    const Iterator begin,
    const Iterator end)
{
    int count = 0;
    real sum = 0;

    for (Iterator it = begin; it != end; ++it)
    {
        ++count;
        sum += *it;
    }

    return sum / count;
}

template<typename Iterator>
real sampleVariance(
    const Iterator begin,
    const Iterator end,
    const real mean)
{
    int count = 0;
    real sum = 0;

    for (Iterator it = begin; it != end; ++it)
    {
        ++count;
        sum += std::pow(*it - mean, 2);
    }

    if (count == 1)
    {
        return 1;
    }

    return sum / (count - 1);
}

#endif //DPEAK_BASIC_STATISTICS_H