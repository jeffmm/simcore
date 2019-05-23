#!/bin/bash

do_build() {
  mkdir build
  cd build
  cmake ..
  make -j8
  cd ..
}

do_clean() {
  rm -rf build/CMake*
  rm -rf build/Makefile
  rm -rf build/cmake_install.cmake
  rm -rf build/src
  rm -rf build/Doxyfile
}

do_usage() {
  echo "Usage: $0 [arg]"
  echo "arg must be one of clean, build, or install"
}

case $1 in
  clean)
    do_clean
    ;;
  build)
    do_build
    ;;
  install)
    do_install
    ;;
  *)
    do_usage
    ;;
esac
