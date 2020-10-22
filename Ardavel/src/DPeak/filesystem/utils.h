//
// Created by Wojciech 'Ardavel' Szałapski
//

#ifndef DPEAK_UTILS_H
#define DPEAK_UTILS_H

#include <string>
#include <vector>

std::vector<std::string> listDirectoryContent(
    const std::string& directoryPath);

std::string getWellName(
    const std::string& wellPath,
    const std::string& plateName);

#endif //DPEAK_UTILS_H