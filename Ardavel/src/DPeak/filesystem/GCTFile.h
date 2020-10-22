//
// Created by Wojciech 'Ardavel' Szałapski
//

#ifndef DPEAK_GCTFILE_H
#define DPEAK_GCTFILE_H

#include "data/types.h"

#include <array>
#include <ostream>
#include <vector>

class GCTFile
{
public:

    explicit GCTFile(
        const std::vector<std::string>& wellsNames);

    void insertValue(
        int barcode,
        int sizeCategory,
        int wellIdx,
        real expression);

    void output(
        const std::string& path);

private:

    std::vector<std::string> m_wellsNames;
    std::vector<std::vector<std::array<real, 2>>> m_wellsPerBarcode;
};

#endif //DPEAK_GCTFILE_H