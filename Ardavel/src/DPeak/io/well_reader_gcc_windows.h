//
// Created by Wojciech 'Ardavel' Szałapski
//

#ifndef DPEAK_WELL_READER_GCC_WINDOWS_H
#define DPEAK_WELL_READER_GCC_WINDOWS_H

#include "io/fast_input.h"
#include "params/Parameters.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <vector>

void readWells(
    const std::string& path,
    std::vector<std::vector<real>>& result)
{
    std::ifstream input(path);

    {
        std::string placeholder;
        input >> placeholder >> placeholder;
    }

    for (auto& entry : result)
    {
        entry.clear();
    }

    int barcode, expIntensity;
    while (input >> barcode >> expIntensity)
    {
        if (expIntensity >= Parameters::MIN_EXP_EXPRESSION &&
            expIntensity <= Parameters::MAX_EXP_EXPRESSION)
        {
            result[barcode].push_back(std::log2(expIntensity));
        }
    }

    for (auto& entry : result)
    {
        std::sort(entry.begin(), entry.end());
    }
}

#endif //DPEAK_WELL_READER_GCC_WINDOWS_H
