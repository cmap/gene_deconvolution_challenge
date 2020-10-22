//
// Created by Wojciech 'Ardavel' Szałapski
//

#include "filesystem/filesystem.h"
#include "utils.h"

std::vector<std::string> listDirectoryContent(
    const std::string& directoryPath)
{
    std::vector<std::string> result;

    for (const auto& entry : filesystem::directory_iterator(directoryPath))
    {
        result.push_back(entry.path().string());
    }

    return result;
}

std::string getWellName(
    const std::string& wellPath,
    const std::string& plateName)
{
    const std::string filename = filesystem::path(wellPath).stem().string();
    return filename.substr(plateName.size() + 1, std::string::npos);
}

