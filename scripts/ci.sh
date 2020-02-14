#!/bin/bash

ci_python() {
    black -l 88 --check . &&\
    flake8 --max-line-length=88 \
           --ignore E203 D202 \
           --max-complexity 16 . &&\
    pytest -s
}

ci_cpp() {
    if [ ! -d "build" ]; then
        mkdir build
    fi
    cd build || exit 1
    cmake -DUNIT_TESTS=1 ..
    make -j8
    make test && cd ..
}

#ci_python && ci_cpp
ci_cpp
