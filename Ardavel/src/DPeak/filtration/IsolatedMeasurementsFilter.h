//
// Created by Wojciech 'Ardavel' Szałapski
//

#ifndef DPEAK_ISOLATEDMEASUREMENTSFILTER_H
#define DPEAK_ISOLATEDMEASUREMENTSFILTER_H

#include "data/types.h"

#include <vector>

class IsolatedMeasurementsFilter
{
public:

    void filter(
        const std::vector<real>& samples,
        std::vector<real>& filteredSamples,
        real oneSideMargin,
        real requiredProportionOfSamples) const;
};

#endif //DPEAK_ISOLATEDMEASUREMENTSFILTER_H
