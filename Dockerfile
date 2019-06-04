FROM frolvlad/alpine-gxx:latest

COPY CMakeLists.txt /app/CMakeLists.txt
COPY src /app/src
COPY docs /app/docs
COPY extern /app/extern
COPY cmake.sh /app/cmake.sh
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
    bash ./cmake.sh clean &&\
    bash ./cmake.sh build &&\
    apk del .build_dependencies &&\
    mkdir /app/run &&\
    mkdir /app/run/outputs &&\
    mv /app/simcore /app/run &&\
    mkdir /app/run/src/ &&\
    mv /app/src/config_params.yaml /app/run/src/config_params.yaml &&\
    touch /app/run/params.yaml


WORKDIR /app/run

CMD ["sh", "-c", "./simcore params.yaml && mv *.[^y]* outputs && cp *.y* outputs"]
