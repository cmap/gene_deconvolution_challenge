//
// Created by Wojciech 'Ardavel' Szałapski
//

#ifndef DPEAK_NORMALLYDISTRIBUTEDCLUSTERFIT_H
#define DPEAK_NORMALLYDISTRIBUTEDCLUSTERFIT_H

#include "data/types.h"

#include <limits>

struct NormallyDistributedClusterFit
{
    real m_mean;
    real m_variance;
    real m_apriori;
};

#endif //DPEAK_NORMALLYDISTRIBUTEDCLUSTERFIT_H
