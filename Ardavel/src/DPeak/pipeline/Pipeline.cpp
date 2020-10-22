//
// Created by Wojciech 'Ardavel' Szałapski
//

#include "Pipeline.h"

#include "adjustment/BVNAdjuster.h"
#include "clustering/GMM.h"
#include "data/BarcodeToGeneMap.h"
#include "data/DeconvolutionResult.h"
#include "filesystem/filesystem.h"
#include "filesystem/GCTFile.h"
#include "filesystem/utils.h"
#include "filtration/IsolatedMeasurementsFilter.h"
#include "io/well_reader.h"
#include "math/basic_statistics.h"
#include "math/median.h"

#include "omp.h"

#include <algorithm>
#include <cmath>
#include <fstream>

void Pipeline::run()
{
    std::vector<std::string> wellsPaths = listDirectoryContent(Parameters::m_inputDir);
    std::sort(wellsPaths.begin(), wellsPaths.end());

    std::vector<DeconvolutionResult> deconvolutionResultPerWell(wellsPaths.size());

    const int numberOfThreads = omp_get_max_threads();

    std::vector<std::vector<real>> filteredSamplesBuffer(static_cast<size_t>(numberOfThreads));

    std::vector<std::vector<std::vector<real>>> allWellData(static_cast<int>(wellsPaths.size()),
                                                            std::vector<std::vector<real>>(
                                                                Parameters::MAXIMUM_BARCODE + 1));

#ifdef NDEBUG
#pragma omp parallel for num_threads(numberOfThreads)
#endif
    for (int wellIdx = 0; wellIdx < static_cast<int>(wellsPaths.size()); ++wellIdx)
    {
        const std::string& wellPath = wellsPaths[wellIdx];

        std::vector<std::vector<real>>& wellData = allWellData[wellIdx];
        readWells(wellPath, wellData);
        DeconvolutionResult& deconvolutionResult = deconvolutionResultPerWell[wellIdx];

        std::vector<std::array<real, 2>> sampleProbabilityBuffer;

        for (int barcode = Parameters::FIRST_RELEVANT_BARCODE; barcode <= Parameters::MAXIMUM_BARCODE; ++barcode)
        {
            if (barcode == 499)
            {
                continue;
            }

            const std::vector<real>& rawSamples = wellData[barcode];

            if (rawSamples.empty())
            {
                deconvolutionResult[barcode] =
                    std::array<ClusterWithVariance, 2>{
                        ClusterWithVariance{Cluster{10, 0}, -1},
                        ClusterWithVariance{Cluster{10, 0}, -1}};

                continue;
            }

            std::vector<real>& samples = filteredSamplesBuffer[omp_get_thread_num()];

            IsolatedMeasurementsFilter isolatedMeasurementsFilter;
            isolatedMeasurementsFilter
                .filter(rawSamples, samples, Parameters::FILTER_ONE_SIDE_MARGIN, Parameters::FILTER_MIN_PROPORTION);
            const auto validBeadsCount = static_cast<int>(samples.size());

            bool fallbackNeeded = true;

            if (validBeadsCount >= Parameters::MIN_BEADS)
            {
                GMM::ModelFit gmmFit = GMM().clusterize(samples,
                                                        std::array<real, 2>{samples[0], samples.back()},
                                                        Parameters::GMM_EM_ITERATIONS,
                                                        Parameters::GMM_EM_MIN_MOVE_OF_MEANS,
                                                        sampleProbabilityBuffer);

                if (gmmFit.m_clusters[0].m_mean > gmmFit.m_clusters[1].m_mean)
                {
                    std::swap(gmmFit.m_clusters[0], gmmFit.m_clusters[1]);
                }

                bool gmmFitValid = gmmFit.m_valid;

                for (int i = 0; i < 2; ++i)
                {
                    const auto& cluster = gmmFit.m_clusters[i];

                    gmmFitValid &= std::isfinite(cluster.m_mean);
                    gmmFitValid &= std::isfinite(cluster.m_variance);
                    gmmFitValid &= (cluster.m_variance > 0);
                    gmmFitValid &= std::isfinite(cluster.m_apriori);
                    gmmFitValid &= cluster.m_apriori > 0;
                    gmmFitValid &= cluster.m_apriori <= 1;
                }

                if (gmmFitValid)
                {
                    std::array<real, 2> clusterSizes;

                    {
                        std::array<int, 3> indexOfFirstSampleInCluster;
                        indexOfFirstSampleInCluster[0] = 0;

                        const real assignmentThreshold =
                            2 * std::log(
                                std::sqrt(gmmFit.m_clusters[0].m_variance / gmmFit.m_clusters[1].m_variance) *
                                    gmmFit.m_clusters[1].m_apriori / gmmFit.m_clusters[0].m_apriori);

                        indexOfFirstSampleInCluster[1] = findIndexOfFirstSampleAssignedToSecondCluster(
                            samples,
                            gmmFit.m_clusters[0],
                            gmmFit.m_clusters[1],
                            assignmentThreshold);

                        indexOfFirstSampleInCluster[2] = validBeadsCount;

                        {
                            for (int i = 0; i < 2; ++i)
                            {
                                const real margin = Parameters::CLUSTER_ASSIGNMENT_MAX_Z_SCORE
                                    * std::sqrt(gmmFit.m_clusters[i].m_variance);
                                const int firstAllowedSample = static_cast<int>(
                                    std::distance(samples.begin(), std::lower_bound(samples.begin(), samples.end(),
                                                                                    gmmFit.m_clusters[i].m_mean
                                                                                        - margin)));
                                const int firstNotAllowedSample = static_cast<int>(
                                    std::distance(samples.begin(), std::upper_bound(samples.begin(), samples.end(),
                                                                                    gmmFit.m_clusters[i].m_mean
                                                                                        + margin)));

                                clusterSizes[i] = std::min(indexOfFirstSampleInCluster[i + 1], firstNotAllowedSample) -
                                    std::max(indexOfFirstSampleInCluster[i], firstAllowedSample);
                            }

                            if (clusterSizes[0] > clusterSizes[1])
                            {
                                std::swap(clusterSizes[0], clusterSizes[1]);
                                std::swap(gmmFit.m_clusters[0], gmmFit.m_clusters[1]);
                            }
                        }
                    }

                    const real smallerClusterSize = clusterSizes[0];
                    const real biggerClusterSize = clusterSizes[1];

                    const real secondBiggestClusterProportion = smallerClusterSize / validBeadsCount;

                    if (secondBiggestClusterProportion >= Parameters::MIN_CLUSTER_RELATIVE_SIZE)
                    {
                        deconvolutionResult[barcode] =
                            std::array<ClusterWithVariance, 2>{
                                ClusterWithVariance{Cluster{gmmFit.m_clusters[0].m_mean, smallerClusterSize},
                                                    gmmFit.m_clusters[0].m_variance},
                                ClusterWithVariance{Cluster{gmmFit.m_clusters[1].m_mean, biggerClusterSize},
                                                    gmmFit.m_clusters[1].m_variance}};

                        fallbackNeeded = false;
                    }

                    if (fallbackNeeded)
                    {
                        const real biggestClusterProportion =
                            static_cast<real>(biggerClusterSize) / validBeadsCount;

                        if (biggestClusterProportion >= Parameters::MIN_CLUSTER_RELATIVE_SIZE)
                        {
                            const real fallbackClusterSize = static_cast<real>(samples.size()) / 2;
                            deconvolutionResult[barcode] = std::array<ClusterWithVariance, 2>{
                                ClusterWithVariance{Cluster{gmmFit.m_clusters[1].m_mean, fallbackClusterSize},
                                                    gmmFit.m_clusters[1].m_variance},
                                ClusterWithVariance{Cluster{gmmFit.m_clusters[1].m_mean, fallbackClusterSize},
                                                    gmmFit.m_clusters[1].m_variance}};

                            fallbackNeeded = false;
                        }
                    }
                }
            }

            if (fallbackNeeded)
            {
                const real fallbackClusterSize = static_cast<real>(validBeadsCount) / 2;
                const real estimatedCenter = (validBeadsCount > 0) ? mean(samples.begin(), samples.end()) : 10;
                const real variance =
                    (validBeadsCount > 1) ? sampleVariance(samples.begin(), samples.end(), estimatedCenter) : -1;

                deconvolutionResult[barcode] =
                    std::array<ClusterWithVariance, 2>{
                        ClusterWithVariance{Cluster{estimatedCenter, fallbackClusterSize}, variance},
                        ClusterWithVariance{Cluster{estimatedCenter, fallbackClusterSize}, variance}};
            }
        }
    }

    std::vector<std::array<BivariateNormalDistribution, 2>> distributionsPerBarcode(Parameters::MAXIMUM_BARCODE + 1);
    std::vector<UnsureResults> unsureResultsPerBarcode = BVNAdjuster()
        .adjust(deconvolutionResultPerWell, distributionsPerBarcode, Parameters::BVN_ADJUSTMENT_ITERATIONS);

#ifdef NDEBUG
#pragma omp parallel for schedule(dynamic, 1) num_threads(numberOfThreads)
#endif
    for (int barcode = Parameters::FIRST_RELEVANT_BARCODE; barcode <= Parameters::MAXIMUM_BARCODE; ++barcode)
    {
        const UnsureResults& unsureResults = unsureResultsPerBarcode[barcode];

        bool invalidDistributions = false;

        for (int i = 0; i < 2; ++i)
        {
            if (!std::isfinite(unsureResults.m_means[i]) ||
                !std::isfinite(unsureResults.m_stdDevs[i]) ||
                unsureResults.m_stdDevs[i] <= 0 ||
                !std::isfinite(unsureResults.m_medianStdDevs[i]) ||
                unsureResults.m_medianStdDevs[i] <= 0)
            {
                invalidDistributions = true;
            }
        }

        if (invalidDistributions)
        {
            continue;
        }

        const std::vector<int>& unsureWellsIndices = unsureResults.m_indices;

        if (unsureWellsIndices.empty())
        {
            continue;
        }

        const real distanceBetweenMeans = std::fabs(unsureResults.m_means[0] - unsureResults.m_means[1]);
        const real sumOfStdDevs = unsureResults.m_stdDevs[0] + unsureResults.m_stdDevs[1];

        if (Parameters::MIN_GAP_BETWEEN_DISTRIBUTIONS_TO_REFINE_ADJUSTMENT * sumOfStdDevs > distanceBetweenMeans)
        {
            continue;
        }

        const real lowClusterSearchMargin =
            Parameters::CLUSTER_ASSIGNMENT_MAX_Z_SCORE * unsureResults.m_medianStdDevs[0];
        const real lowClusterSearchSpan = 2 * lowClusterSearchMargin;

        for (const int& unsureWellIdx : unsureWellsIndices)
        {
            bool lowClusterGenerated = false;
            bool sameMeans = (deconvolutionResultPerWell[unsureWellIdx][barcode][0].m_cluster[0] ==
                deconvolutionResultPerWell[unsureWellIdx][barcode][1].m_cluster[0]);

            const std::vector<real>& rawSamples = allWellData[unsureWellIdx][barcode];
            const auto numberOfRawSamples = static_cast<int>(rawSamples.size());

            if (numberOfRawSamples == 0)
            {
                continue;
            }

            const real start =
                unsureResults.m_means[0] - Parameters::LOW_CLUSTER_SEARCH_RANGE * unsureResults.m_stdDevs[0]
                    - lowClusterSearchMargin;
            const real stop =
                unsureResults.m_means[0] + Parameters::LOW_CLUSTER_SEARCH_RANGE * unsureResults.m_stdDevs[0]
                    + lowClusterSearchMargin;

            int highestNumberOfSamples = 0;
            int firstSampleIdxInNewLowCluster = 0;
            int pastLastSampleIdxInNewLowCluster = 0;

            for (int sampleIdx = 0; sampleIdx < numberOfRawSamples; ++sampleIdx)
            {
                const real value = rawSamples[sampleIdx];
                if (value <= start)
                {
                    continue;
                }

                const real maximum = std::min(value + lowClusterSearchSpan, stop);

                const int upperBoundIdx = static_cast<int>(std::distance(rawSamples.begin(),
                                                                         std::upper_bound(rawSamples.begin(),
                                                                                          rawSamples.end(),
                                                                                          maximum)));

                const int numberOfSamples = upperBoundIdx - sampleIdx;
                if (numberOfSamples > highestNumberOfSamples)
                {
                    highestNumberOfSamples = numberOfSamples;
                    firstSampleIdxInNewLowCluster = sampleIdx;
                    pastLastSampleIdxInNewLowCluster = upperBoundIdx;
                }
            }

            const real ratio = static_cast<real>(highestNumberOfSamples) / numberOfRawSamples;

            if (ratio > Parameters::MIN_RATIO_TO_CREATE_LOW_CLUSTER)
            {
                real newLowClusterMean;

                if ((pastLastSampleIdxInNewLowCluster - firstSampleIdxInNewLowCluster) % 2)
                {
                    newLowClusterMean = rawSamples[firstSampleIdxInNewLowCluster
                        + (pastLastSampleIdxInNewLowCluster - firstSampleIdxInNewLowCluster) / 2];
                }
                else
                {
                    newLowClusterMean = (rawSamples[firstSampleIdxInNewLowCluster
                        + (pastLastSampleIdxInNewLowCluster - firstSampleIdxInNewLowCluster) / 2 - 1]
                        + rawSamples[firstSampleIdxInNewLowCluster
                            + (pastLastSampleIdxInNewLowCluster - firstSampleIdxInNewLowCluster) / 2]) / 2;
                }

                const int firstForbiddenSampleIdx = static_cast<int>(std::min(
                    std::distance(rawSamples.begin(),
                                  std::lower_bound(rawSamples.begin(), rawSamples.end(),
                                                   deconvolutionResultPerWell[unsureWellIdx][barcode]
                                                   [0].m_cluster[0]
                                                       - Parameters::FORBIDDEN_Z_SCORE_FOR_LOW_CLUSTER_GENERATION
                                                           * std::sqrt(
                                                               deconvolutionResultPerWell[unsureWellIdx][barcode]
                                                               [0].m_variance))),
                    std::distance(rawSamples.begin(),
                                  std::lower_bound(rawSamples.begin(), rawSamples.end(),
                                                   deconvolutionResultPerWell[unsureWellIdx][barcode]
                                                   [1].m_cluster[0]
                                                       - Parameters::FORBIDDEN_Z_SCORE_FOR_LOW_CLUSTER_GENERATION
                                                           * std::sqrt(
                                                               deconvolutionResultPerWell[unsureWellIdx][barcode]
                                                               [1].m_variance)))));

                const int pastLastForbiddenSampleIdx = static_cast<int>(std::max(
                    std::distance(rawSamples.begin(),
                                  std::upper_bound(rawSamples.begin(), rawSamples.end(),
                                                   deconvolutionResultPerWell[unsureWellIdx][barcode]
                                                   [0].m_cluster[0]
                                                       + Parameters::FORBIDDEN_Z_SCORE_FOR_LOW_CLUSTER_GENERATION
                                                           * std::sqrt(
                                                               deconvolutionResultPerWell[unsureWellIdx][barcode]
                                                               [0].m_variance))),
                    std::distance(rawSamples.begin(),
                                  std::upper_bound(rawSamples.begin(), rawSamples.end(),
                                                   deconvolutionResultPerWell[unsureWellIdx][barcode]
                                                   [1].m_cluster[0]
                                                       + Parameters::FORBIDDEN_Z_SCORE_FOR_LOW_CLUSTER_GENERATION
                                                           * std::sqrt(
                                                               deconvolutionResultPerWell[unsureWellIdx][barcode]
                                                               [1].m_variance)))));

                if (pastLastForbiddenSampleIdx <= firstForbiddenSampleIdx)
                {
                    deconvolutionResultPerWell[unsureWellIdx][barcode][0].m_cluster[0] = newLowClusterMean;
                    lowClusterGenerated = true;
                }
                else
                {
                    const real firstForbiddenSampleValue = rawSamples[firstForbiddenSampleIdx];
                    const real lastForbiddenSampleValue = rawSamples[pastLastForbiddenSampleIdx - 1];

                    if (newLowClusterMean < firstForbiddenSampleValue || newLowClusterMean > lastForbiddenSampleValue)
                    {
                        if (deconvolutionResultPerWell[unsureWellIdx][barcode][0].m_cluster[0] !=
                            deconvolutionResultPerWell[unsureWellIdx][barcode][1].m_cluster[0])
                        {
                            const real margin =
                                Parameters::HIGH_OVERLAP_MAX_Z_SCORE
                                    * (std::sqrt(deconvolutionResultPerWell[unsureWellIdx][barcode][0].m_variance) +
                                        std::sqrt(deconvolutionResultPerWell[unsureWellIdx][barcode][1].m_variance));

                            const real meanDifference = std::fabs(
                                deconvolutionResultPerWell[unsureWellIdx][barcode][0].m_cluster[0]
                                    - deconvolutionResultPerWell[unsureWellIdx][barcode][1].m_cluster[0]);

                            if (meanDifference < margin)
                            {
                                const int firstSampleIdxInAnyOfClusters = static_cast<int>(std::min(
                                    std::distance(rawSamples.begin(),
                                                  std::lower_bound(rawSamples.begin(), rawSamples.end(),
                                                                   deconvolutionResultPerWell[unsureWellIdx][barcode]
                                                                   [0].m_cluster[0] - 3 * std::sqrt(
                                                                       deconvolutionResultPerWell[unsureWellIdx][barcode]
                                                                       [0].m_variance))),
                                    std::distance(rawSamples.begin(),
                                                  std::lower_bound(rawSamples.begin(), rawSamples.end(),
                                                                   deconvolutionResultPerWell[unsureWellIdx][barcode]
                                                                   [1].m_cluster[0] - 3 * std::sqrt(
                                                                       deconvolutionResultPerWell[unsureWellIdx][barcode]
                                                                       [1].m_variance)))));

                                const int pastLastSampleIdxInAnyOfClusters = static_cast<int>(std::max(
                                    std::distance(rawSamples.begin(),
                                                  std::upper_bound(rawSamples.begin(), rawSamples.end(),
                                                                   deconvolutionResultPerWell[unsureWellIdx][barcode]
                                                                   [0].m_cluster[0] + 3 * std::sqrt(
                                                                       deconvolutionResultPerWell[unsureWellIdx][barcode]
                                                                       [0].m_variance))),
                                    std::distance(rawSamples.begin(),
                                                  std::upper_bound(rawSamples.begin(), rawSamples.end(),
                                                                   deconvolutionResultPerWell[unsureWellIdx][barcode]
                                                                   [1].m_cluster[0] + 3 * std::sqrt(
                                                                       deconvolutionResultPerWell[unsureWellIdx][barcode]
                                                                       [1].m_variance)))));

                                if (pastLastSampleIdxInAnyOfClusters > firstSampleIdxInAnyOfClusters)
                                {
                                    const real newHighClusterMean =
                                        mean(rawSamples.begin() + firstSampleIdxInAnyOfClusters,
                                             rawSamples.begin() + pastLastSampleIdxInAnyOfClusters);
                                    deconvolutionResultPerWell[unsureWellIdx][barcode][1].m_cluster[0] =
                                        newHighClusterMean;
                                }
                            }
                        }

                        deconvolutionResultPerWell[unsureWellIdx][barcode][0].m_cluster[0] = newLowClusterMean;
                        lowClusterGenerated = true;
                    }
                }
            }

            if (!lowClusterGenerated && !sameMeans)
            {
                std::array<BivariateNormalDistribution, 2>& bvns = distributionsPerBarcode[barcode];
                if (bvns[0].isInitialized() && bvns[1].isInitialized())
                {
                    const real stayProbability =
                        bvns[0].density(deconvolutionResultPerWell[unsureWellIdx][barcode][0].m_cluster)
                            * bvns[1].density(deconvolutionResultPerWell[unsureWellIdx][barcode][1].m_cluster);
                    const real swapProbability =
                        bvns[0].density(deconvolutionResultPerWell[unsureWellIdx][barcode][1].m_cluster)
                            * bvns[1].density(deconvolutionResultPerWell[unsureWellIdx][barcode][0].m_cluster);

                    if (swapProbability > stayProbability)
                    {
                        std::swap(deconvolutionResultPerWell[unsureWellIdx][barcode][0].m_cluster,
                                  deconvolutionResultPerWell[unsureWellIdx][barcode][1].m_cluster);
                    }
                }
            }
        }
    }

    std::vector<std::string> wellsNames;
    wellsNames.reserve(wellsPaths.size());

    for (const std::string& wellPath : wellsPaths)
    {
        wellsNames.push_back(getWellName(wellPath, Parameters::m_plateName));
    }

    GCTFile gexFile(wellsNames);

    for (int wellIdx = 0; wellIdx < static_cast<int>(wellsPaths.size()); ++wellIdx)
    {
        const DeconvolutionResult& deconvolutionResult = deconvolutionResultPerWell[wellIdx];

        for (const auto& entry : deconvolutionResult)
        {
            const int& barcode = entry.first;
            const std::array<ClusterWithVariance, 2>& clusters = entry.second;
            const PairedGenes& pairedGenes = barcodeToGeneMap[barcode];

            for (int sizeCategory = 0; sizeCategory <= 1; ++sizeCategory)
            {
                const int gene = pairedGenes.m_geneIds[sizeCategory];
                if (gene != -1)
                {
                    real expressionValue = std::pow(static_cast<real>(2), clusters[sizeCategory].m_cluster[0]);
                    if (!std::isfinite(expressionValue))
                    {
                        expressionValue = Parameters::FALLBACK_EXPRESSION;
                    }
                    gexFile.insertValue(barcode, sizeCategory, wellIdx, expressionValue);
                }
            }
        }
    }

    {
        filesystem::create_directories(Parameters::m_outputDir);
        const std::string outputPath =
            (filesystem::path(Parameters::m_outputDir) / (Parameters::m_plateName + ".gct")).string();
        gexFile.output(outputPath);
    }
}

int Pipeline::findIndexOfFirstSampleAssignedToSecondCluster(
    const std::vector<real>& samples,
    const NormallyDistributedClusterFit& cluster1,
    const NormallyDistributedClusterFit& cluster2,
    const real& thresholdValue) const
{
    if (samples.empty())
    {
        return 0;
    }

    int low = static_cast<int>(std::distance(samples.begin(),
                                             std::lower_bound(samples.begin(), samples.end(), cluster1.m_mean)));
    int high = static_cast<int>(std::distance(samples.begin(),
                                              std::lower_bound(samples.begin(), samples.end(), cluster2.m_mean)));
    const int end = high;

    if (high != 0)
    {
        --high;
    }

    while (low < high)
    {
        const int mid = (low + high) / 2;
        if (std::pow(samples[mid] - cluster2.m_mean, 2) / cluster2.m_variance -
            std::pow(samples[mid] - cluster1.m_mean, 2) / cluster1.m_variance < thresholdValue)
        {
            high = mid;
        }
        else
        {
            low = mid + 1;
        }
    }

    if (std::pow(samples[low] - cluster2.m_mean, 2) / cluster2.m_variance -
        std::pow(samples[low] - cluster1.m_mean, 2) / cluster1.m_variance < thresholdValue)
    {
        return low;
    }

    return end;
}
