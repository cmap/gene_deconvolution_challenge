//
// Created by Wojciech 'Ardavel' Szałapski
//

#ifndef DPEAK_WELL_READER_H
#define DPEAK_WELL_READER_H

#ifdef _MSC_VER
#include "io/well_reader_msvc.h"
#elif __linux__
#include "io/well_reader_gcc_linux.h"
#else
#include "io/well_reader_gcc_windows.h"
#endif

#endif //DPEAK_WELL_READER_H
