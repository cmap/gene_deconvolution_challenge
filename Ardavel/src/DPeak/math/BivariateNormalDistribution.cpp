//
// Created by Wojciech 'Ardavel' Szałapski
//

#include "BivariateNormalDistribution.h"

#include <cmath>

void BivariateNormalDistribution::initialize(
    const std::vector<std::array<real, 2>>& samples)
{
    m_initialized = false;

    m_mean = {0, 0};
    std::array<real, 2> meanOfSquares{0, 0};
    real meanOfProducts = 0;

    for (const auto& dataPoint : samples)
    {
        for (int i = 0; i < 2; ++i)
        {
            m_mean[i] += dataPoint[i];
            meanOfSquares[i] += dataPoint[i] * dataPoint[i];
        }

        meanOfProducts += dataPoint[0] * dataPoint[1];
    }

    for (int i = 0; i < 2; ++i)
    {
        m_mean[i] /= samples.size();
        meanOfSquares[i] /= samples.size();
    }
    meanOfProducts /= samples.size();

    std::array<real, 2> squaredDeviation{0, 0};

    for (const auto& sample : samples)
    {
        const std::array<real, 2>& dataPoint = sample;

        for (int i = 0; i < 2; ++i)
        {
            squaredDeviation[i] += std::pow(dataPoint[i] - m_mean[i], 2);
        }
    }

    for (int i = 0; i < 2; ++i)
    {
        m_standardDeviation[i] = std::sqrt(squaredDeviation[i] / samples.size());
    }

    m_correlation = meanOfProducts - m_mean[0] * m_mean[1];
    for (int i = 0; i < 2; ++i)
    {
        m_correlation /= std::sqrt(meanOfSquares[i] - m_mean[i] * m_mean[i]);
    }

    if (!std::isfinite(m_correlation) ||
        m_correlation <= -1 ||
        m_correlation >= 1)
    {
        return;
    }

    for (int i = 0; i < 2; ++i)
    {
        if (!std::isfinite(m_standardDeviation[i]) ||
            m_standardDeviation[i] <= 0)
        {
            return;
        }

        if (!std::isfinite(m_mean[i]))
        {
            return;
        }
    }

    m_initialized = true;
}

real BivariateNormalDistribution::density(
    const std::array<real, 2>& value) const
{
    real result = 0;

    for (int i = 0; i < 2; ++i)
    {
        result += std::pow((value[i] - m_mean[i]) / m_standardDeviation[i], 2);
    }

    real subtrahend = 2 * m_correlation;
    for (int i = 0; i < 2; ++i)
    {
        subtrahend *= (value[i] - m_mean[i]) / m_standardDeviation[i];
    }

    result -= subtrahend;

    result /= -2 * (1 - std::pow(m_correlation, 2));
    result = std::exp(result);

    result /= 2 * PI * std::sqrt(1 - pow(m_correlation, 2));
    for (int i = 0; i < 2; ++i)
    {
        result /= m_standardDeviation[i];
    }

    return result;
}

bool BivariateNormalDistribution::isInitialized() const
{
    return m_initialized;
}