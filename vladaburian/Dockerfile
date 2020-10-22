FROM ubuntu:cosmic

RUN apt-get update
RUN apt-get install -y gcc python3 python3-setuptools python3-dev python3-venv
RUN apt-get install -y pybind11-dev libeigen3-dev libtbb-dev

COPY gene /workspace/gene
COPY setup.py /workspace/

RUN python3 -m venv /app && /app/bin/pip3 install /workspace


# Prepare "clean" dist image
FROM ubuntu:cosmic

RUN apt-get update && apt-get install -y \
    python3 libtbb2 \
 && rm -rf /var/lib/apt/lists/*

COPY --from=0 /app /app

ENTRYPOINT ["/app/bin/python3", "-m", "gene.cmd_main"]
