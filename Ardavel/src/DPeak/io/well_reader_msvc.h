//
// Created by Wojciech 'Ardavel' Szałapski
//

#ifndef DPEAK_WELL_READER_MSVC_H
#define DPEAK_WELL_READER_MSVC_H

#include "io/fast_input.h"
#include "params/Parameters.h"

#include <boost/iostreams/device/mapped_file.hpp>

#include <algorithm>
#include <cmath>
#include <vector>

void readWells(
	const std::string& path,
	std::vector<std::vector<real>>& result)
{
	boost::iostreams::mapped_file mappedFile(path, boost::iostreams::mapped_file::readonly);
	const char* dataPtr = mappedFile.const_data();
	const char* const endPtr = dataPtr + mappedFile.size();

	for (auto& entry : result)
	{
		entry.clear();
	}

	int barcode = 0;
	int expIntensity;

	while (barcode != -1)
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
	}

	for (auto& entry : result)
	{
	    std::sort(entry.begin(), entry.end());
	}
}

#endif //DPEAK_WELL_READER_MSVC_H
