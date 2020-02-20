#!/bin/bash

do_build() {
    mkdir build
    cd build || exit 1
    cmake ${CMAKE_FLAGS} ..
    make -j8
    if $build_docs; then
        make docs
    fi
    if $run_tests; then
        make test
    fi
    if $install_packages; then
        make install
    fi
    cd ..
}

clean_cmake_files() {
    rm -rf build/CMake*
    rm -rf build/Makefile
    rm -rf build/cmake_install.cmake
    rm -rf build/src
    rm -rf build/Doxyfile
    rm -rf build/tests
}

show_help() {
    echo "USAGE:"
    echo "  $0 [-hcIgwodtbx]"
    echo "OPTIONS:"
    echo "  -h      show this menu"
    echo "  -c      clean build directory"
    echo "  -I      install simcore locally after building"
    echo "  -g      build simcore with graphics"
    echo "  -w      build simcore in Windows"
    echo "  -o      build simcore with OpenMP parallelization"
    echo "  -d      build Doxygen documentation"
    echo "  -t      build and run simcore unit tests"
    echo "  -b      build simcore in Debug mode"
    echo "  -x      build simcore in Trace mode (verbose logging)"
}

# A POSIX variable
OPTIND=1         # Reset in case getopts has been used previously in the shell.

CMAKE_FLAGS="-DCMAKE_EXPORT_COMPILE_COMMANDS=ON"
build_docs=false
run_tests=false
install_packages=false
while getopts "h?cIgwodtDT:" opt; do
    case "$opt" in
    h|\?)
        show_help
        exit 0
        ;;
    c)  
        clean_cmake_files
        exit 0
        ;;
    I)  install_packages=true
        ;;
    g)  CMAKE_FLAGS="${CMAKE_FLAGS} -DGRAPH=TRUE"
        ;;
    w)  CMAKE_FLAGS="${CMAKE_FLAGS} -DWINDOWS_MODE=TRUE"
        ;;
    o)  CMAKE_FLAGS="${CMAKE_FLAGS} -DOMP=TRUE"
        ;;
    d)  build_docs=true
        ;;
    t)  
        run_tests=true
        CMAKE_FLAGS="${CMAKE_FLAGS} -DTESTS=TRUE"
        ;;
    D)  CMAKE_FLAGS="${CMAKE_FLAGS} -DDEBUG=TRUE"
        ;;
    T)  CMAKE_FLAGS="${CMAKE_FLAGS} -DTRACE=TRUE"
        ;;
    esac
done

shift $((OPTIND-1))

[ "${1:-}" = "--" ] && shift

do_build
