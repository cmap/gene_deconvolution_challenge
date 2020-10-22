//
// Created by Wojciech 'Ardavel' Szałapski
//

#include "ParameterParser.h"

void ParameterParser::parse(
    const int argc,
    const char* const* const argv)
{
    for (int idx = 1; idx < argc - 1; idx += 2)
    {
        m_parameterMap[argv[idx]] = argv[idx + 1];
    }

    Parameters::m_inputDir = get("--dspath");
    Parameters::m_outputDir = get("--out");
    Parameters::m_plateName = get("--plate");
}

std::string ParameterParser::get(
    const std::string& key)
{
    return m_parameterMap[key];
}
