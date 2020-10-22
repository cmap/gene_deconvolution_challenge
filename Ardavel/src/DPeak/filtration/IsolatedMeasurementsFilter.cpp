//
// Created by Wojciech 'Ardavel' Szałapski
//

#include "IsolatedMeasurementsFilter.h"

#include <algorithm>
#include <cmath>

void IsolatedMeasurementsFilter::filter(
    const std::vector<real>& samples,
    std::vector<real>& filteredSamples,
    const real oneSideMargin,
    const real requiredProportionOfSamples) const
{
    const auto n = static_cast<int>(samples.size());
    const auto minimumNumberOfSamplesWithinMargin = static_cast<int>(std::round(n * requiredProportionOfSamples));

    filteredSamples.clear();

    int lowerBound = 0;
    int upperBound = 0;

    for (int i = 0; i < n; ++i)
    {
        const real value = samples[i];
        const real leftMargin = value - oneSideMargin;
        const real rightMargin = value + oneSideMargin;

        while (samples[lowerBound] < leftMargin)
        {
            ++lowerBound;
        }

        while (upperBound < n && samples[upperBound] <= rightMargin)
        {
            ++upperBound;
        }

        const int samplesWithinMargin = upperBound - lowerBound;

        if (samplesWithinMargin >= minimumNumberOfSamplesWithinMargin)
        {
            filteredSamples.push_back(value);
        }
    }
}
