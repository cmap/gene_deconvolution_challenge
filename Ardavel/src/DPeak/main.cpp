//
// Created by Wojciech 'Ardavel' Szałapski
//

#include "params/ParameterParser.h"
#include "pipeline/Pipeline.h"

int main(
    const int argc,
    const char* const* const argv)
{
    ParameterParser().parse(argc, argv);

    Pipeline pipeline;
    pipeline.run();

    return 0;
}