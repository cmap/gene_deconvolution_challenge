cmake_minimum_required(VERSION 3.12.1)
project(DPeak)

if (MSVC)
    set(BOOST_ROOT "" CACHE PATH "A path to the Boost root directory.")
endif ()

set(CMAKE_CXX_STANDARD 17)

set(BUILD_FLAGS "-std=c++17 -Wall -fopenmp")

if (NOT DEFINED DEBUG_BUILD OR DEBUG_BUILD EQUAL "0")
    set(BUILD_FLAGS "${BUILD_FLAGS} -O3 -DNDEBUG -s -ffast-math")
endif ()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${BUILD_FLAGS}")

set(SOURCES
        main.cpp params/ParameterParser.cpp params/ParameterParser.h params/Parameters.h pipeline/Pipeline.cpp
        pipeline/Pipeline.h filesystem/utils.h filesystem/utils.cpp filesystem/filesystem.h
        data/types.h data/DeconvolutionResult.h
        math/median.h math/median.cpp data/PairedGenes.h data/BarcodeToGeneMap.h
        filesystem/GCTFile.h filesystem/GCTFile.cpp clustering/GMM.h math/basic_statistics.h
        math/probability.h adjustment/BVNAdjuster.cpp adjustment/BVNAdjuster.h data/Cluster.h
        math/BivariateNormalDistribution.cpp math/BivariateNormalDistribution.h
        filtration/IsolatedMeasurementsFilter.cpp filtration/IsolatedMeasurementsFilter.h math/probability.cpp
        math/NormallyDistributedClusterFit.h
        io/fast_input.h)

if (MSVC)
    set(WELL_READER io/well_reader_msvc.h)
elseif (UNIX)
    set(WELL_READER io/well_reader_gcc_linux.h)
else ()
    set(WELL_READER io/well_reader_gcc_windows.h)
endif ()

set(SOURCES "${SOURCES};${WELL_READER}" io/well_reader.h data/UnsureResults.h params/Parameters.cpp clustering/GMM.cpp)

add_executable(DPeak ${SOURCES})

target_include_directories(DPeak PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

if (NOT MSVC)
    target_link_libraries(DPeak PRIVATE stdc++fs)
else ()
    set(Boost_USE_STATIC_LIBS ON)

    find_package(
            Boost
            COMPONENTS
            iostreams
            REQUIRED)

    target_include_directories(
            DPeak
            PRIVATE
            ${Boost_INCLUDE_DIR})

    target_link_libraries(
            DPeak
            PRIVATE
            ${Boost_LIBRARIES})
endif ()