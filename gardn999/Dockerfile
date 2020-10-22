
FROM ubuntu:16.04

## Install General Requirements
RUN apt-get update && \
    apt-get install -y openjdk-8-jdk && \
    apt-get install -y --no-install-recommends \
    apt-utils \
    build-essential \
    cmake \
    git \
    wget \
    nano \
    software-properties-common
		
WORKDIR /work

# copy entire directory where docker file is into docker container at /work
COPY . /work/

RUN chmod 777 compile.sh
RUN chmod 777 train.sh
RUN chmod 777 test.sh

ENTRYPOINT [ "/work/test.sh" ]
