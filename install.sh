#!/bin/bash

do_build() {
    mkdir build
    cd build || exit 1
    cmake ..
    make -j8
    cd ..
}

do_omp_build() {
    mkdir build
    cd build || exit 1
    cmake -DOMP=1 ..
    make -j8
    cd ..
}

do_graph_omp_build() {
    mkdir build
    cd build || exit 1
    cmake -DGRAPH=1 -DOMP=1 ..
    make -j8
    cd ..
}

do_graph_build() {
    mkdir build
    cd build || exit 1
    cmake -DGRAPH=1 ..
    make -j8
    cd ..
}

do_windows_graph_build() {
    mkdir build
    cd build || exit 1
    cmake -DGRAPH=1 -DWINDOWS=1 ..
    make -j8
    cd ..
}

do_test_build() {
    mkdir build
    cd build || exit 1
    cmake -DTESTS=1 ..
    make -j8
    make test
    cd ..
}

do_debug_build() {
    mkdir build
    cd build || exit 1
    cmake -DDEBUG=1 ..
    make -j8
    cd ..
}

do_debug_graph_build() {
    mkdir build
    cd build || exit 1
    cmake -DDEBUG=1 -DGRAPH=1 ..
    make -j8
    cd ..
}

do_trace_build() {
    mkdir build
    cd build || exit 1
    cmake -DDEBUG=1 -DTRACE=1 ..
    make -j8
    cd ..
}

do_trace_graph_build() {
    mkdir build
    cd build || exit 1
    cmake -DDEBUG=1 -DTRACE=1 -DGRAPH=1 ..
    make -j8
    cd ..
}

do_docs_build() {
    mkdir build
    cd build || exit 1
    cmake ..
    make docs
    cd ..
}

do_clean() {
    rm -rf build/CMake*
    rm -rf build/Makefile
    rm -rf build/cmake_install.cmake
    rm -rf build/src
    rm -rf build/Doxyfile
    rm -rf build/tests
}

do_usage() {
    echo "Usage: $0 [arg]"
    echo "arg must be one of:"
    echo "  clean   - remove temporary installation files"
    echo "  build   - build simcore without graphics"
    echo "  gbuild  - build simcore with graphics"
    echo "  omp     - build simcore with openmp without graphics"
    echo "  gomp    - build simcore with openmp with graphics"
    echo "  debug   - build simcore in debug mode without graphics"
    echo "  gdebug  - build simcore in debug mode with graphics"
    echo "  gwin    - build simcore in windows with graphics"
    echo "  trace   - build simcore in trace mode (verbose logging)"
    echo "  gtrace  - build simcore in trace mode with graphics (verbose logging)"
    echo "  test    - build simcore and run unit tests"
    echo "  docs    - build Doxygen documentation"
}

case $1 in
clean)
    do_clean
    ;;
gwin)
    do_windows_graph_build
    ;;
build)
    do_build
    ;;
debug)
    do_debug_build
    ;;
omp)
    do_omp_build
    ;;
gomp)
    do_graph_omp_build
    ;;
gbuild)
    do_graph_build
    ;;
gdebug)
    do_debug_graph_build
    ;;
test)
    do_test_build
    ;;
docs)
    do_docs_build
    ;;
trace)
    do_trace_build
    ;;
gtrace)
    do_trace_graph_build
    ;;
*)
    do_usage
    ;;
esac
