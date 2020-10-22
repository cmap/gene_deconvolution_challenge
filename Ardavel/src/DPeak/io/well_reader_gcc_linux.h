//
// Created by Wojciech 'Ardavel' Szałapski
//

#ifndef DPEAK_WELL_READER_MSVC_H
#define DPEAK_WELL_READER_MSVC_H

#include "io/fast_input.h"
#include "params/Parameters.h"

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

#include <algorithm>
#include <cmath>
#include <vector>

void readWells(
	const std::string& path,
    std::vector<std::vector<real>>& result)
{
    const int fileDescriptor = open(path.c_str(), O_RDONLY);

    struct stat fileStatus;
    fstat(fileDescriptor, &fileStatus);

    const size_t dataSize = fileStatus.st_size;

    const char* dataPtr = static_cast<const char*>(mmap(NULL, dataSize, PROT_READ, MAP_PRIVATE, fileDescriptor, 0u));
    const char* const endPtr = dataPtr + dataSize;

    for (auto& entry : result)
    {
        entry.clear();
    }

    int barcode = 0;
    int expIntensity;

    while (true)
    {
        barcode = readUnsignedInt(dataPtr, endPtr);

        if (barcode != -1)
        {
            expIntensity = readUnsignedInt(dataPtr, endPtr);

            if (expIntensity >= Parameters::MIN_EXP_EXPRESSION &&
                expIntensity <= Parameters::MAX_EXP_EXPRESSION)
            {
                result[barcode].push_back(std::log2(expIntensity));
            }
        }
        else
        {
            break;
        }
    }

    close(fileDescriptor);

    for (auto& entry : result)
    {
        std::sort(entry.begin(), entry.end());
    }
}

#endif //DPEAK_WELL_READER_MSVC_H
