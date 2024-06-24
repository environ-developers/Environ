# Use containerized installation
# 0) Install Docker (https://docs.docker.com/engine/install/) or Docker Desktop (https://www.docker.com/products/docker-desktop/)
# 1) Clone repo:     git clone https://github.com/environ-developers/Environ.git && cd Environ
# 2) Build image:    docker build -t environ-sandbox . && docker image ls
# 3) Run container:  docker run -i -p 8000:80 environ-sandbox && docker ps
# 4) Access container:
#    a) Docker Desktop
#    b) CLI: docker exec -it <container ID> bash
# 5) Review test and/or example logs:
#    a) cat /app/tests.log
#    b) cat /app/examples.log

# Configure sandbox -- compiler, dependencies, openmpi, python

FROM ubuntu:22.04
RUN apt-get update
RUN DEBIAN_FRONTEND=noninteractive apt-get install --assume-yes --no-install-recommends \
    autoconf \
    bc \
    build-essential \
    ca-certificates \
    gfortran \
    git-all \
    libblas3 \
    libc6 \
    libelpa17 \
    libfftw3-dev \
    libgcc-s1 \
    liblapack-dev \
    libopenmpi-dev \
    libscalapack-openmpi-dev \
    python3 \
    pipx \
    wget 
#RUN pipx install Cython numpy ase matplotlib pandas notebook --include-deps
WORKDIR /app

# Set environment variables
# N                               number of processors
# QE_VERSION                      Quantum ESPRESSO source version
# QE_CONFIG                       Quantum ESPRESSO configuration options
# ENVIRON_CONFIG                  Environ configuration options
# DOWNLOAD_EFFICIENCY             download efficiency SSSP pseudopotential during build
# DOWNLOAD_PRECISION              download precision SSSP pseudopotential during build
# RUN_TESTS                       run tests during build
# RUN_EXAMPLES                    run examples during build
# BIN_DIR                         absolute path to binary directory (QE variable)
# PSEUDO_DIR                      absolute path to pseudopotential directory (QE variable)
# TMP_DIR                         absolute path of placeholder directory (QE variable)
# OMPI_ALLOW_RUN_AS_ROOT          root-level execution permission (MPI variable)
# OMPI_ALLOW_RUN_AS_ROOT_CONFIRM  root-level execution confirmation (MPI variable)
ENV N=4 \
    QE_VERSION=qe-7.2 \
    QE_CONFIG="--enable-parallel --with-environ=/app/q-e/Environ" \
    ENVIRON_CONFIG="--with-qe=/app/q-e" \
    DOWNLOAD_EFFICIENCY=0 \
    DOWNLOAD_PRECISION=0 \
    RUN_TESTS=1 \
    RUN_EXAMPLES=0 \
    BIN_DIR="/app/q-e/bin" \
    PSEUDO_DIR="/app/q-e/pseudo" \
    TMP_DIR="/app" \
    OMPI_ALLOW_RUN_AS_ROOT=1 \
    OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1

# Compile & install QE with Environ
    
RUN git clone https://github.com/QEF/q-e.git
WORKDIR /app/q-e
RUN git checkout $QE_VERSION && git clone https://github.com/environ-developers/Environ.git
RUN grep -rl .*NETWORK_PSEUDO=http:.* | xargs sed -i "s/NETWORK_PSEUDO=http:/NETWORK_PSEUDO=https:/"
RUN grep -rl .*NETWORK_PSEUDO=https:\/\/www.quantum-espresso.org\/wp-content\/uploads\/upf_files\/.* | xargs sed -i "s/NETWORK_PSEUDO=https:\/\/www.quantum-espresso.org\/wp-content\/uploads\/upf_files\//NETWORK_PSEUDO=https:\/\/pseudopotentials.quantum-espresso.org\/upf_files\//"
WORKDIR /app/q-e/Environ
RUN ./configure $ENVIRON_CONFIG && make -j $N compile
WORKDIR /app/q-e
RUN ./configure $QE_CONFIG && make -j $N pw
RUN export PATH="$BIN_DIR:$PATH"
WORKDIR /app

# Download SSSP pseudopotentials

ENV EFFICIENCY_TARBALL="https://archive.materialscloud.org/record/file?record_id=1732&filename=SSSP_1.3.0_PBE_efficiency.tar.gz"
ENV PRECISION_TARBALL="https://archive.materialscloud.org/record/file?record_id=1732&filename=SSSP_1.3.0_PBE_precision.tar.gz"
RUN if [ "$DOWNLOAD_EFFICIENCY" = "1" ]; then \
    wget $EFFICIENCY_TARBALL --output-document=efficiency.tar.gz && \
    tar -xvzf efficiency.tar.gz -C /app/q-e/pseudo && \
    rm efficiency.tar.gz; fi
RUN if [ "$DOWNLOAD_PRECISION" = "1" ]; then \
    wget $PRECISION_TARBALL --output-document=precision.tar.gz && \
    tar -xvzf precision.tar.gz -C /app/q-e/pseudo && \
    rm precision.tar.gz; fi

# Run Environ tests and examples

WORKDIR /app/q-e/Environ/tests
RUN if [ "$RUN_TESTS" = "1" ]; then \
    make --ignore-errors run-tests-parallel | tee /app/tests.log; fi
WORKDIR /app/q-e/Environ/examples/qe/pw
RUN if [ "$RUN_EXAMPLES" = "1" ]; then \
    ./run_all.sh | tee /app/examples.log; fi
WORKDIR /app

# Run shell

CMD bash