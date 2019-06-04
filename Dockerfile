FROM frolvlad/alpine-gxx:latest

COPY CMakeLists.txt /app/CMakeLists.txt
COPY src /app/src
COPY docs /app/docs
COPY extern /app/extern
COPY install.sh /app/install.sh
COPY .CMake_FFTW /app/.CMake_FFTW

RUN apk add --no-cache --virtual=.build_dependencies build-base &&\
    apk add --no-cache \
                    bash \
                    cmake \
                    yaml-cpp-dev \
                    boost-dev \
                    fftw-dev \
                    gsl-dev \
                    doxygen \
                    boost-dev &&\
    cd /app && \
    bash ./install.sh clean &&\
    bash ./install.sh build &&\
    apk del .build_dependencies &&\
    mkdir /app/run &&\
    mkdir /app/run/outputs &&\
    mv /app/simcore /app/run &&\
    mkdir /app/run/src/ &&\
    mv /app/src/config_params.yaml /app/run/src/config_params.yaml &&\
    touch /app/run/params.yaml


WORKDIR /app/run

CMD ["bash", "-c", "./simcore params.yaml && cp *.y* outputs && if ls *.[^y]* 1> /dev/null 2>&1; then mv *.[^y]* outputs; fi"]
