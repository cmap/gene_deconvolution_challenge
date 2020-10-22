
CMap DPeak Challenge
Deployment Guide

// Docker instructions for making predictions //////////////////////////////////

to build for Docker: 
  docker build -t <image name> <solution dir with Dockerfile>

to run using Docker:
  docker run --rm -it \
    -v <input path>:/input \
    -v <output path>:/output \
	<image name> \
      --dspath /input/<input set directory> \
      --out /output \
      --create_subdir 0 \
      --plate <output plate name without .gct extension>

// Command-line instructions ///////////////////////////////////////////////////

requirements: 
  Java version 1.8 or later is required

to compile: 
  sh compile.sh

to train:
  sh train.sh <DPK input dir> <DPK truth file> <LITMUS input dir> <LITMUS truth file>
  
  (this is a shortcut for running the java .jar file)
  java -Xmx2G -jar gardn999_CMap_DPeak.jar <DPK input dir> <DPK truth file> 
                                     <LITMUS input dir> <LITMUS truth file>

to make predictions: 
  sh predict.sh <input dir> <output dir> <plate name>

  (this is a shortcut for running the java .jar file)
  java -Xmx2G -jar gardn999_CMap_DPeak.jar <input dir> <output dir> <plate name>

////////////////////////////////////////////////////////////////////////////////
