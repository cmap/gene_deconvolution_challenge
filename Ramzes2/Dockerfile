FROM romanchernenko/cpp_parser

COPY predict.py .
COPY libs libs
COPY model.hdf5 .
RUN mkdir -p data/input
COPY data/input/barcode_to_gene_map.txt data/input

RUN python predict.py --python-cache true

ENTRYPOINT ["python", "predict.py"]
