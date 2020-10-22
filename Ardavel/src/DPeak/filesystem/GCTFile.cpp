//
// Created by Wojciech 'Ardavel' Szałapski
//

#include "GCTFile.h"

#include "data/BarcodeToGeneMap.h"
#include "params/Parameters.h"

#include <cstdio>
#include <string>

GCTFile::GCTFile(
    const std::vector<std::string>& wellNames)
    : m_wellsNames(wellNames),
      m_wellsPerBarcode(Parameters::MAXIMUM_BARCODE + 1,
                        std::vector<std::array<real, 2>>(
                            wellNames.size(),
                            std::array<real, 2>{Parameters::FALLBACK_EXPRESSION,
                                                Parameters::FALLBACK_EXPRESSION}))
{
}

void GCTFile::insertValue(
    const int barcode,
    const int sizeCategory,
    const int wellIdx,
    const real expression)
{
    m_wellsPerBarcode[barcode][wellIdx][sizeCategory] = expression;
}

void GCTFile::output(
    const std::string& path)
{
    std::FILE* const file = fopen(path.c_str(), "w");

    std::fprintf(file, "#1.3\n976\t%d\t0\t0\nid", static_cast<int>(m_wellsNames.size()));

    for (const std::string& wellName : m_wellsNames)
    {
        std::fprintf(file, "\t%s", wellName.c_str());
    }

    std::fprintf(file, "\n");

    for (int barcode = Parameters::FIRST_RELEVANT_BARCODE; barcode <= Parameters::MAXIMUM_BARCODE; ++barcode)
    {
        if (barcode == 499)
        {
            continue;
        }

        for (int i = 0; i < 2; ++i)
        {
            std::fprintf(file, "%d", barcodeToGeneMap[barcode].m_geneIds[i]);

            for (int wellIdx = 0; wellIdx < static_cast<int>(m_wellsNames.size()); ++wellIdx)
            {
                std::fprintf(file, "\t%.1f", m_wellsPerBarcode[barcode][wellIdx][i]);
            }

            std::fprintf(file, "\n");
        }
    }

    std::fclose(file);
}
