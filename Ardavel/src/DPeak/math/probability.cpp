//
// Created by Wojciech 'Ardavel' Szałapski
//

#include "math/probability.h"

#include <cmath>

real dnorm(
    const real value,
    const real mean,
    const real variance,
    const real denominator)
{
    return static_cast<real>(std::exp(-std::pow(value - mean, 2) / (2 * variance)) / denominator);
}