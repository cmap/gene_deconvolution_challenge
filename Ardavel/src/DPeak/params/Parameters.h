//
// Created by Wojciech 'Ardavel' Szałapski
//

#ifndef DPEAK_PARAMETERS_H
#define DPEAK_PARAMETERS_H

#include "data/types.h"

#include <string>

class Parameters
{
public:

    static std::string m_inputDir;

    static std::string m_outputDir;
    static std::string m_plateName;

    static constexpr int MIN_EXP_EXPRESSION = 16;
    static constexpr int MAX_EXP_EXPRESSION = 32768;

    static constexpr int FIRST_RELEVANT_BARCODE = 12;
    static constexpr int MAXIMUM_BARCODE = 500;

    static constexpr int FALLBACK_EXPRESSION = 1024;

    static constexpr int MIN_BEADS = 20;

    static constexpr real MIN_CLUSTER_RELATIVE_SIZE = 0.1;

    static constexpr real FILTER_ONE_SIDE_MARGIN = 0.625;
    static constexpr real FILTER_MIN_PROPORTION = 0.15;

    static constexpr int GMM_EM_ITERATIONS = 10;
    static constexpr real GMM_EM_MIN_MOVE_OF_MEANS = 0.001;

    static constexpr real MIN_CLUSTER_VARIANCE = 0.0001;

    static constexpr real CLUSTER_ASSIGNMENT_MAX_Z_SCORE = 3;

    static constexpr int BVN_ADJUSTMENT_ITERATIONS = 20;
    static constexpr real HIGH_OVERLAP_MAX_Z_SCORE = 1.25;
    static constexpr real FAR_FROM_LOW_CLUSTER_MIN_Z_SCORE = 4;
    static constexpr real FAR_FROM_LOW_CLUSTER_MIN_COMBINED_Z_SCORE = 3;
    static constexpr real NEAR_TO_HIGH_CLUSTER_MAX_Z_SCORE = 1;
    static constexpr real NEAR_TO_HIGH_CLUSTER_MAX_Z_SCORE_SUM = 2;
    static constexpr real CONFIDENT_PAIR_MAX_Z_SCORE = 1.5;

    static constexpr real MIN_GAP_BETWEEN_DISTRIBUTIONS_TO_REFINE_ADJUSTMENT = 1;
    static constexpr real LOW_CLUSTER_SEARCH_RANGE = 1.5;
    static constexpr real MIN_RATIO_TO_CREATE_LOW_CLUSTER = 0.1;
    static constexpr real FORBIDDEN_Z_SCORE_FOR_LOW_CLUSTER_GENERATION = 2;
};

#endif //DPEAK_PARAMETERS_H