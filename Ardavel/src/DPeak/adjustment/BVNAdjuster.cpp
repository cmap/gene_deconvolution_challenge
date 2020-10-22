//
// Created by Wojciech 'Ardavel' Szałapski
//

#include "BVNAdjuster.h"

#include "math/basic_statistics.h"
#include "math/median.h"
#include "params/Parameters.h"

#include "omp.h"

#include <algorithm>

std::vector<UnsureResults> BVNAdjuster::adjust(
    std::vector<DeconvolutionResult>& deconvolutionResults,
    std::vector<std::array<BivariateNormalDistribution, 2>>& distributionsPerBarcode,
    const int maximumNumberOfIterations) const
{
    std::vector<UnsureResults> unsureResults(Parameters::MAXIMUM_BARCODE + 1);

#ifdef NDEBUG
#pragma omp parallel for num_threads(omp_get_max_threads())
#endif
    for (int barcode = Parameters::FIRST_RELEVANT_BARCODE; barcode <= Parameters::MAXIMUM_BARCODE; ++barcode)
    {
        std::vector<std::map<int, std::array<ClusterWithVariance, 2>>::iterator> singleClusterResultsRef;
        std::vector<int> singleClusterResultsWells;

        std::vector<std::map<int, std::array<ClusterWithVariance, 2>>::iterator>
            unsureTwoClustersRef;
        std::vector<int> unsureTwoClustersWells;

        std::vector<std::map<int, std::array<ClusterWithVariance, 2>>::iterator>
            candidatesForMergeLowClusterGenerationRef;
        std::vector<int> candidatesForMergeLowClusterGenerationWells;

        std::vector<std::map<int, std::array<ClusterWithVariance, 2>>::iterator> candidatesRefs;
        std::array<std::vector<real>, 2> candidates;
        std::vector<int> candidatesWells;

        {
            int wellIdx = -1;
            for (DeconvolutionResult& wellResult : deconvolutionResults)
            {
                ++wellIdx;

                auto barcodeResultIt = wellResult.find(barcode);

                if (barcodeResultIt != wellResult.end() &&
                    std::isfinite(barcodeResultIt->second[0].m_variance) &&
                    std::isfinite(barcodeResultIt->second[1].m_variance) &&
                    barcodeResultIt->second[0].m_variance >= 0 && barcodeResultIt->second[1].m_variance >= 0 &&
                    barcodeResultIt->second[0].m_cluster[1] > 0 && barcodeResultIt->second[1].m_cluster[1] > 0)
                {
                    if (barcodeResultIt->second[0].m_cluster[0] != barcodeResultIt->second[1].m_cluster[0])
                    {
                        const real margin =
                            Parameters::HIGH_OVERLAP_MAX_Z_SCORE * (std::sqrt(barcodeResultIt->second[0].m_variance) +
                                std::sqrt(barcodeResultIt->second[1].m_variance));
                        const real meanDifference = std::fabs(barcodeResultIt->second[0].m_cluster[0] -
                            barcodeResultIt->second[1].m_cluster[0]);

                        if (meanDifference < margin)
                        {
                            unsureTwoClustersRef.push_back(barcodeResultIt);
                            unsureTwoClustersWells.push_back(wellIdx);
                        }
                        else
                        {
                            candidatesRefs.push_back(barcodeResultIt);

                            for (int i = 0; i < 2; ++i)
                            {
                                candidates[i].push_back(barcodeResultIt->second[i].m_cluster[0]);
                            }

                            candidatesWells.push_back(wellIdx);
                        }
                    }
                    else
                    {
                        singleClusterResultsRef.push_back(barcodeResultIt);
                        singleClusterResultsWells.push_back(wellIdx);
                    }
                }
            }
        }

        if (candidates[0].size() < 50)
        {
            continue;
        }

        std::array<real, 2> means;

        for (int i = 0; i < 2; ++i)
        {
            means[i] = mean(candidates[i].begin(), candidates[i].end());
        }

        std::array<real, 2> vars;

        for (int i = 0; i < 2; ++i)
        {
            vars[i] = sampleVariance(candidates[i].begin(), candidates[i].end(), means[i]);
        }

        if (!std::isfinite(vars[0]) || vars[0] <= 0 || !std::isfinite(vars[1]) || vars[1] <= 0)
        {
            continue;
        }

        std::vector<std::map<int, std::array<ClusterWithVariance, 2>>::iterator> resultsRef;
        std::vector<std::array<std::array<real, 2>, 2>> results;

        for (int candidateIdx = 0; candidateIdx < static_cast<int>(candidatesRefs.size()); ++candidateIdx)
        {
            real zScores[2][2];

            for (int i = 0; i < 2; ++i)
            {
                for (int j = 0; j < 2; ++j)
                {
                    zScores[i][j] = std::fabs(candidates[i][candidateIdx] - means[j]) / std::sqrt(vars[j]);
                }
            }

            if (std::max(zScores[0][0], zScores[1][1]) <= Parameters::CONFIDENT_PAIR_MAX_Z_SCORE ||
                std::max(zScores[0][1], zScores[1][0]) <= Parameters::CONFIDENT_PAIR_MAX_Z_SCORE)
            {
                if (std::max(zScores[0][0], zScores[1][1]) > Parameters::CONFIDENT_PAIR_MAX_Z_SCORE)
                {
                    std::swap(candidatesRefs[candidateIdx]->second[0], candidatesRefs[candidateIdx]->second[1]);
                }

                resultsRef.push_back(candidatesRefs[candidateIdx]);
                results
                    .push_back(std::array<std::array<real, 2>, 2>{candidatesRefs[candidateIdx]->second[0].m_cluster,
                                                                  candidatesRefs[candidateIdx]->second[1].m_cluster});
            }
            else
            {
                if (std::min(zScores[0][0], zScores[1][0]) > Parameters::FAR_FROM_LOW_CLUSTER_MIN_Z_SCORE &&
                    std::max(zScores[0][1], zScores[1][1]) <= Parameters::NEAR_TO_HIGH_CLUSTER_MAX_Z_SCORE &&
                    zScores[0][1] + zScores[1][1] <= Parameters::NEAR_TO_HIGH_CLUSTER_MAX_Z_SCORE_SUM)
                {
                    candidatesForMergeLowClusterGenerationRef.push_back(candidatesRefs[candidateIdx]);
                    candidatesForMergeLowClusterGenerationWells.push_back(candidatesWells[candidateIdx]);
                }
            }
        }

        const auto numberOfResults = static_cast<int>(resultsRef.size());

        if (numberOfResults < 40)
        {
            continue;
        }

        std::vector<std::array<real, 2>> samples[2];
        for (auto& entry : samples)
        {
            entry.reserve(static_cast<size_t>(numberOfResults));
        }

        std::vector<int> assignment[2];

        for (int clusterIdx = 0; clusterIdx < 2; ++clusterIdx)
        {
            assignment[clusterIdx].resize(static_cast<size_t>(numberOfResults), clusterIdx);
        }

        std::array<BivariateNormalDistribution, 2>& bvns = distributionsPerBarcode[barcode];

        for (int iteration = 0; iteration < maximumNumberOfIterations; ++iteration)
        {
            samples[0].clear();
            samples[1].clear();

            for (int originalClusterIdx = 0; originalClusterIdx < 2; ++originalClusterIdx)
            {
                for (int sampleIdx = 0; sampleIdx < numberOfResults; ++sampleIdx)
                {
                    samples[assignment[originalClusterIdx][sampleIdx]]
                        .push_back(results[sampleIdx][originalClusterIdx]);
                }
            }

            for (int i = 0; i < 2; ++i)
            {
                bvns[i].initialize(samples[i]);
            }

            if (!bvns[0].isInitialized() || !bvns[1].isInitialized())
            {
                break;
            }

            bool anyChange = false;

            for (int sampleIdx = 0; sampleIdx < numberOfResults; ++sampleIdx)
            {
                std::array<real, 2> stayVsSwapProbability = {1, 1};

                for (int originalClusterIdx = 0; originalClusterIdx < 2; ++originalClusterIdx)
                {
                    const std::array<real, 2>& sample = results[sampleIdx][originalClusterIdx];

                    for (int candidateClusterIdx = 0; candidateClusterIdx < 2; ++candidateClusterIdx)
                    {
                        stayVsSwapProbability[originalClusterIdx != candidateClusterIdx] *=
                            bvns[candidateClusterIdx].density(sample);
                    }
                }

                if (std::isfinite(stayVsSwapProbability[0]) && std::isfinite(stayVsSwapProbability[1]))
                {
                    if (assignment[0][sampleIdx] != (stayVsSwapProbability[1] > stayVsSwapProbability[0]))
                    {
                        std::swap(assignment[0][sampleIdx], assignment[1][sampleIdx]);
                        anyChange = true;
                    }
                }
            }

            if (!anyChange)
            {
                break;
            }
        }

        if (!bvns[0].isInitialized() || !bvns[1].isInitialized())
        {
            continue;
        }

        for (int sampleIdx = 0; sampleIdx < numberOfResults; ++sampleIdx)
        {
            const real stayProbability = bvns[0].density(resultsRef[sampleIdx]->second[0].m_cluster)
                * bvns[1].density(resultsRef[sampleIdx]->second[1].m_cluster);
            const real swapProbability = bvns[0].density(resultsRef[sampleIdx]->second[1].m_cluster)
                * bvns[1].density(resultsRef[sampleIdx]->second[0].m_cluster);

            if (std::isfinite(swapProbability) && std::isfinite(stayProbability) && swapProbability > stayProbability)
            {
                std::swap(resultsRef[sampleIdx]->second[0], resultsRef[sampleIdx]->second[1]);
            }
        }

        std::array<std::vector<real>, 2> buffer;

        for (int sampleIdx = 0; sampleIdx < numberOfResults; ++sampleIdx)
        {
            for (int i = 0; i < 2; ++i)
            {
                buffer[i].push_back(resultsRef[sampleIdx]->second[i].m_cluster[0]);
            }
        }

        for (int i = 0; i < 2; ++i)
        {
            means[i] = mean(buffer[i].begin(), buffer[i].end());
        }

        for (int i = 0; i < 2; ++i)
        {
            vars[i] = sampleVariance(buffer[i].begin(), buffer[i].end(), means[i]);
        }

        if (!std::isfinite(vars[0]) || vars[0] <= 0 || !std::isfinite(vars[1]) || vars[1] <= 0)
        {
            continue;
        }

        std::vector<std::map<int, std::array<ClusterWithVariance, 2>>::iterator> usedResultsRef;

        for (int sampleIdx = 0; sampleIdx < numberOfResults; ++sampleIdx)
        {
            usedResultsRef.push_back(resultsRef[sampleIdx]);

            int smallIdx = 0;
            int bigIdx = 1;
            if (resultsRef[sampleIdx]->second[0].m_cluster[1] > resultsRef[sampleIdx]->second[1].m_cluster[1])
            {
                std::swap(smallIdx, bigIdx);
            }

            constexpr std::array<std::array<real, 2>, 1> ratiosWithZScoreTresholds = {{{0.6, 4}}};

            const real smallClusterSize = resultsRef[sampleIdx]->second[smallIdx].m_cluster[1];
            const real bigClusterSize = resultsRef[sampleIdx]->second[bigIdx].m_cluster[1];
            const real bigSizeRatio = bigClusterSize / (smallClusterSize + bigClusterSize);

            int bigSizeRatioIndex = 0;

            while (bigSizeRatioIndex < static_cast<int>(ratiosWithZScoreTresholds.size()) &&
                bigSizeRatio < ratiosWithZScoreTresholds[bigSizeRatioIndex][0])
            {
                ++bigSizeRatioIndex;
            }

            if (bigSizeRatioIndex < static_cast<int>(ratiosWithZScoreTresholds.size()))
            {
                const real smallZScore = std::fabs(resultsRef[sampleIdx]->second[smallIdx].m_cluster[0] - means[0]) /
                    std::sqrt(vars[0]);
                const real bigZScore = std::fabs(resultsRef[sampleIdx]->second[bigIdx].m_cluster[0] - means[1]) /
                    std::sqrt(vars[1]);

                if (std::max(smallZScore, bigZScore) > ratiosWithZScoreTresholds[bigSizeRatioIndex][1])
                {
                    if (smallIdx == 1)
                    {
                        std::swap(resultsRef[sampleIdx]->second[0], resultsRef[sampleIdx]->second[1]);
                        usedResultsRef.pop_back();
                    }
                }
            }
        }

        const auto numberOfUsedResults = static_cast<int>(usedResultsRef.size());

        if (numberOfUsedResults < 30)
        {
            continue;
        }

        for (int i = 0; i < 2; ++i)
        {
            buffer[i].clear();
            samples[i].clear();
        }

        for (int sampleIdx = 0; sampleIdx < numberOfUsedResults; ++sampleIdx)
        {
            for (int i = 0; i < 2; ++i)
            {
                buffer[i].push_back(usedResultsRef[sampleIdx]->second[i].m_cluster[0]);
                samples[i].push_back(usedResultsRef[sampleIdx]->second[i].m_cluster);
            }
        }

        for (int i = 0; i < 2; ++i)
        {
            bvns[i].initialize(samples[i]);
        }

        for (int i = 0; i < 2; ++i)
        {
            means[i] = mean(buffer[i].begin(), buffer[i].end());
        }

        for (int i = 0; i < 2; ++i)
        {
            vars[i] = sampleVariance(buffer[i].begin(), buffer[i].end(), means[i]);
        }

        if (!std::isfinite(vars[0]) || vars[0] <= 0 || !std::isfinite(vars[1]) || vars[1] <= 0)
        {
            continue;
        }

        for (int idx = 0; idx < static_cast<int>(unsureTwoClustersRef.size()); ++idx)
        {
            const auto& entry = unsureTwoClustersRef[idx];

            real zScores[2][2];
            real combinedZScores[2][2];

            for (int j = 0; j < 2; ++j)
            {
                for (int i = 0; i < 2; ++i)
                {
                    zScores[i][j] = std::fabs(entry->second[i].m_cluster[0] - means[j]) / std::sqrt(vars[j]);
                    combinedZScores[i][j] = std::fabs(entry->second[i].m_cluster[0] - means[j])
                        / (std::sqrt(vars[j]) + std::sqrt(entry->second[i].m_variance));
                }
            }

            const real estimatedClustersOverlapZScore =
                std::fabs(entry->second[0].m_cluster[0] - entry->second[1].m_cluster[0])
                    / (std::sqrt(entry->second[0].m_variance) + std::sqrt(entry->second[1].m_variance));

            if (std::min(zScores[0][0], zScores[1][0]) > Parameters::FAR_FROM_LOW_CLUSTER_MIN_Z_SCORE &&
                ((estimatedClustersOverlapZScore <= Parameters::HIGH_OVERLAP_MAX_Z_SCORE &&
                    std::min(combinedZScores[0][0], combinedZScores[1][0])
                        > Parameters::FAR_FROM_LOW_CLUSTER_MIN_COMBINED_Z_SCORE) ||
                    (std::max(zScores[0][1], zScores[1][1]) <= Parameters::NEAR_TO_HIGH_CLUSTER_MAX_Z_SCORE &&
                        zScores[0][1] + zScores[1][1] <= Parameters::NEAR_TO_HIGH_CLUSTER_MAX_Z_SCORE_SUM)))
            {
                candidatesForMergeLowClusterGenerationRef.push_back(entry);
                candidatesForMergeLowClusterGenerationWells.push_back(unsureTwoClustersWells[idx]);
            }
            else
            {
                const real stayProbability =
                    bvns[0].density(entry->second[0].m_cluster) * bvns[1].density(entry->second[1].m_cluster);
                const real swapProbability =
                    bvns[0].density(entry->second[1].m_cluster) * bvns[1].density(entry->second[0].m_cluster);

                if (std::isfinite(swapProbability) && std::isfinite(stayProbability)
                    && swapProbability > stayProbability)
                {
                    std::swap(entry->second[0], entry->second[1]);
                }
            }
        }

        unsureResults[barcode].m_means = means;
        unsureResults[barcode].m_stdDevs = vars;

        for (auto& sd : unsureResults[barcode].m_stdDevs)
        {
            sd = std::sqrt(sd);
        }

        for (int idx = 0; idx < static_cast<int>(candidatesForMergeLowClusterGenerationRef.size()); ++idx)
        {
            unsureResults[barcode].m_indices.push_back(candidatesForMergeLowClusterGenerationWells[idx]);
        }

        for (int idx = 0; idx < static_cast<int>(singleClusterResultsRef.size()); ++idx)
        {
            const auto& entry = singleClusterResultsRef[idx];

            real zScores[2];

            for (int j = 0; j < 2; ++j)
            {
                zScores[j] = std::fabs(entry->second[0].m_cluster[0] - means[j]) / std::sqrt(vars[j]);
            }

            if (zScores[1] <= 3 && zScores[0] > 3)
            {
                unsureResults[barcode].m_indices.push_back(singleClusterResultsWells[idx]);
            }
        }

        if (!unsureResults[barcode].m_indices.empty())
        {
            for (int i = 0; i < 2; ++i)
            {
                buffer[i].clear();
            }

            for (const auto& entry : usedResultsRef)
            {
                for (int i = 0; i < 2; ++i)
                {
                    buffer[i].push_back(std::sqrt(entry->second[i].m_variance));
                }
            }

            for (int i = 0; i < 2; ++i)
            {
                std::sort(buffer[i].begin(), buffer[i].end());
                unsureResults[barcode].m_medianStdDevs[i] = calculateMedianOfSortedData(buffer[i]);
            }
        }
    }

    return unsureResults;
}
