FROM debian:10.2-slim
WORKDIR /build
RUN apt-get update -qq && \
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
        libboost-math1.67-dev &&\
    rm -rf ~/.cache
RUN git clone --recursive --single-branch --branch jmm/directory_tree_refactor https://github.com/jeffmm/simcore.git . &&\
    ./install.sh -otI
RUN rm -rf /build
WORKDIR /mnt
ENTRYPOINT []
