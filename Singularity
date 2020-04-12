Bootstrap: docker
From: debian:10.2-slim

%post
    apt-get update -qq
    apt-get install -qqy --no-install-recommends \
        build-essential \
        ca-certificates \
        vim \
        git \
        cmake \
        wget \
        curl \
        htop \
        pkg-config \
        doxygen \
        libyaml-cpp-dev \
        libgsl-dev \
        libfftw3-dev \
        libopenmpi-dev \
        libboost-math1.67-dev
    rm -rf ~/.cache
    mkdir /build
    cd /build
    git clone --recursive --single-branch https://github.com/jeffmm/simcore.git .
    ./install.sh -otI
    rm -rf /build

%runscript
    simcore.exe $*

%test
    simcore.exe --version

%help
    Usage: `singularity run <simcore container> [flags...] <parameter file>`

