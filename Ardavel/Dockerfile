FROM gcc:8.2.0

RUN mkdir /usr/src/cmake
RUN cd /usr/src/cmake && wget https://github.com/Kitware/CMake/releases/download/v3.13.2/cmake-3.13.2-Linux-x86_64.sh && chmod u+rwx ./cmake-3.13.2-Linux-x86_64.sh && ./cmake-3.13.2-Linux-x86_64.sh --skip-license

COPY ./src/DPeak /usr/src/solution
WORKDIR /usr/src/build

RUN /usr/src/cmake/bin/cmake /usr/src/solution
RUN make

ENTRYPOINT ["./DPeak"]