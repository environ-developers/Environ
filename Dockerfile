# Use containerized installation
# 0) Install Docker (https://docs.docker.com/engine/install/) or Docker Desktop (https://www.docker.com/products/docker-desktop/)
# 1) Clone repo:     git clone https://github.com/environ-developers/Environ.git && cd Environ
# 2) Build image:    docker build -t environ-sandbox . && docker image ls
# 3) Run container:  docker run -i -p 8000:80 environ-sandbox && docker ps
# 4) Access container:
#    a) Docker Desktop
#    b) CLI: docker exec -it <container ID> bash

# Configure sandbox -- compiler, dependencies, python

FROM ubuntu
RUN apt-get update
RUN apt-get install --assume-yes --no-install-recommends \
    autoconf \
    build-essential \
    ca-certificates \
    gfortran \
    git-all \
    libblas3 \
    libc6 \
    libfftw3-dev \
    libgcc-s1 \
    liblapack-dev \
    python3 \
    pipx \
    wget 
RUN pipx install Cython numpy ase matplotlib pandas notebook --include-deps
WORKDIR /app
COPY . /app

# Set environment variables
# N               number of processors
# QE_VERSION      Quantum ESPRESSO source version
# QE_CONFIG       Quantum ESPRESSO configuration options
# ENVIRON_CONFIG  Environ configuration options
# BIN_DIR         absolute path to binary directory (QE variable)
# PSEUDO_DIR      absolute path to pseudopotential directory (QE variable)
# TMP_DIR         absolute path of placeholder directory (QE variable)
ENV N=1 \
    QE_VERSION=qe-7.2 \
    QE_CONFIG="--with-environ=/app/q-e/Environ" \
    ENVIRON_CONFIG="--with-qe=/app/q-e" \
    BIN_DIR="/app/q-e/bin" \
    PSEUDO_DIR="/app/sssp_efficiency" \
    TMP_DIR="/app"

# Install Environ
    
RUN git clone https://github.com/QEF/q-e.git
WORKDIR /app/q-e
RUN git checkout $QE_VERSION && git clone https://github.com/environ-developers/Environ.git
WORKDIR /app/q-e/Environ
RUN ./configure $ENVIRON_CONFIG && make -j $N compile
WORKDIR /app/q-e
RUN ./configure $QE_CONFIG && make -j $N pw
RUN alias pw=/app/q-e/bin/pw.x
WORKDIR /app

# Download SSSP pseudopotentials

RUN mkdir sssp_efficiency sssp_precision
RUN wget "https://archive.materialscloud.org/record/file?record_id=1732&filename=SSSP_1.3.0_PBE_efficiency.tar.gz" --output-document=efficiency.tar.gz && \
    tar -xvzf efficiency.tar.gz -C sssp_efficiency
RUN wget "https://archive.materialscloud.org/record/file?record_id=1732&filename=SSSP_1.3.0_PBE_precision.tar.gz" --output-document=precision.tar.gz && \
    tar -xvzf precision.tar.gz -C sssp_precision
RUN rm efficiency.tar.gz precision.tar.gz

# Run shell

CMD bash