#!/bin/bash

ci_python() {
    black -l 88 --check . &&\
    flake8 --max-line-length=88 \
           --ignore E203 D202 \
           --max-complexity 10 . &&\
    pytest -s
}

ci_cpp() {
    if [ ! -d "build" ]; then
        mkdir build
    fi
    cd build || exit 1
    cmake ..
    make -j8
    make test
    cd ..
}

ci_python
ci_cpp
