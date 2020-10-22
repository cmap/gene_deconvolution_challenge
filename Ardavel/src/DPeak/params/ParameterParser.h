//
// Created by Wojciech 'Ardavel' Szałapski
//

#ifndef DPEAK_PARAMETERPARSER_H
#define DPEAK_PARAMETERPARSER_H

#include "Parameters.h"

#include <map>

class ParameterParser
{
public:

    void parse(
        int argc,
        const char* const* argv);

    std::string get(
        const std::string& key);

private:

    std::map<std::string, std::string> m_parameterMap;
};

#endif //DPEAK_PARAMETERPARSER_H
