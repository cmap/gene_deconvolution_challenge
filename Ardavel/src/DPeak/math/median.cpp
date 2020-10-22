//
// Created by Wojciech 'Ardavel' Szałapski
//

#include "median.h"

real calculateMedianOfSortedData(
    const std::vector<real>& sortedData)
{
    const auto n = static_cast<int>(sortedData.size());

    if (n % 2)
    {
        return sortedData[n / 2];
    }
    else
    {
        return (sortedData[n / 2 - 1] + sortedData[n / 2]) / 2;
    }
}