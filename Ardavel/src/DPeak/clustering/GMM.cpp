//
// Created by Wojciech 'Ardavel' Szałapski
//

#include "clustering/GMM.h"

#include "math/basic_statistics.h"
#include "math/probability.h"
#include "params/Parameters.h"

#include <algorithm>

GMM::ModelFit GMM::clusterize(
    const std::vector<real>& samples,
    const std::array<real, 2>& initialCenters,
    const int maxNumberOfIterations,
    const real minMoveOfMeans,
    std::vector<std::array<real, 2>>& sampleProbability) const
{
    const auto n = static_cast<int>(samples.size());

    ModelFit bestFit;
    std::array<NormallyDistributedClusterFit, 2> currentFit;

    {
        const int firstSampleIdxInSecondCluster = static_cast<int>(
            std::distance(samples.begin(),
                          std::upper_bound(samples.begin(),
                                           samples.end(),
                                           (initialCenters[0] + initialCenters[1]) / 2)));

        if (firstSampleIdxInSecondCluster < 2 || firstSampleIdxInSecondCluster > n - 2)
        {
            return bestFit;
        }

        std::array<int, 3> initialClusterBoundaryIndices = {0, firstSampleIdxInSecondCluster, n};

        for (int i = 0; i < 2; ++i)
        {
            const auto beginIt = samples.begin() + initialClusterBoundaryIndices[i];
            const auto endIt = samples.begin() + initialClusterBoundaryIndices[i + 1];

            currentFit[i].m_mean = mean(beginIt, endIt);
            currentFit[i].m_variance = sampleVariance(beginIt, endIt, currentFit[i].m_mean);
            currentFit[i].m_apriori =
                static_cast<real>(initialClusterBoundaryIndices[i + 1] - initialClusterBoundaryIndices[i]) / n;

            if (!std::isfinite(currentFit[i].m_mean) || !std::isfinite(currentFit[i].m_variance) ||
                currentFit[i].m_variance <= 0)
            {
                return bestFit;
            }
        }
    }

    bestFit.m_valid = true;
    bestFit.m_clusters = currentFit;

    sampleProbability.resize(samples.size());

    std::array<real, 2> unscaledSampleProbability;
    std::array<real, 2> sumOfProbabilities;
    std::array<real, 2> denominator;
    std::array<real, 2> weighedMean;
    std::array<real, 2> varianceAccumulator;

    for (int iteration = 0; iteration < maxNumberOfIterations; ++iteration)
    {
        for (int i = 0; i < 2; ++i)
        {
            denominator[i] = std::sqrt(2 * PI * currentFit[i].m_variance);
        }

        std::fill(sumOfProbabilities.begin(), sumOfProbabilities.end(), 0);

        for (int sampleIdx = 0; sampleIdx < n; ++sampleIdx)
        {
            const real sample = samples[sampleIdx];

            real sampleProbabilityNormalizingFactor = 0;

            for (int i = 0; i < 2; ++i)
            {
                const NormallyDistributedClusterFit& clusterFit = currentFit[i];

                unscaledSampleProbability[i] = clusterFit.m_apriori
                    * dnorm(sample, clusterFit.m_mean, clusterFit.m_variance, denominator[i]);

                sampleProbabilityNormalizingFactor += unscaledSampleProbability[i];
            }

            for (int i = 0; i < 2; ++i)
            {
                sampleProbability[sampleIdx][i] =
                    unscaledSampleProbability[i] / sampleProbabilityNormalizingFactor;

                sumOfProbabilities[i] += sampleProbability[sampleIdx][i];
            }
        }

        std::fill(weighedMean.begin(), weighedMean.end(), 0);

        for (int sampleIdx = 0; sampleIdx < n; ++sampleIdx)
        {
            const real sample = samples[sampleIdx];

            for (int i = 0; i < 2; ++i)
            {
                weighedMean[i] += sample * sampleProbability[sampleIdx][i];
            }
        }

        real moveOfMeans = 0;

        for (int i = 0; i < 2; ++i)
        {
            const real newMean = weighedMean[i] / sumOfProbabilities[i];
            moveOfMeans += std::fabs(newMean - currentFit[i].m_mean);
            currentFit[i].m_mean = newMean;
        }

        if (moveOfMeans < minMoveOfMeans)
        {
            break;
        }

        std::fill(varianceAccumulator.begin(), varianceAccumulator.end(), 0);

        for (int sampleIdx = 0; sampleIdx < n; ++sampleIdx)
        {
            for (int i = 0; i < 2; ++i)
            {
                varianceAccumulator[i] +=
                    sampleProbability[sampleIdx][i] * std::pow(samples[sampleIdx] - currentFit[i].m_mean, 2);
            }
        }

        {
            bool invalidVarianceAppeared = false;

            for (int i = 0; i < 2; ++i)
            {
                NormallyDistributedClusterFit& clusterFit = currentFit[i];

                clusterFit.m_variance = varianceAccumulator[i] / sumOfProbabilities[i];
                clusterFit.m_apriori = sumOfProbabilities[i] / n;

                if (std::isnan(clusterFit.m_variance) || std::isinf(clusterFit.m_variance)
                    || clusterFit.m_variance < Parameters::MIN_CLUSTER_VARIANCE)
                {
                    invalidVarianceAppeared = true;
                    break;
                }
            }

            if (invalidVarianceAppeared)
            {
                break;
            }
        }

        bestFit.m_clusters = currentFit;
    }

    return bestFit;
}
