FROM gcc:7.2
COPY script.cpp ./
ADD barcode_to_gene_map.txt ./
RUN g++ -std=c++11 -pthread -o myapp script.cpp
ENTRYPOINT ["./myapp"]
