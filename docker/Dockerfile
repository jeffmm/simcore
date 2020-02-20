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
        libboost-math1.67-dev &&\
        rm -rf ~/.cache
RUN git clone --recursive --single-branch --branch jmm/directory_tree_refactor http://github.com/jeffmm/simcore.git . &&\
    mkdir build &&\
    cd build &&\
    cmake -DTESTS=TRUE .. &&\
    make -j4 &&\
    make install
RUN rm -rf /build
WORKDIR /mnt
ENTRYPOINT ["bash"]
