#!/bin/bash
g++-6 -std=c++11 -g *.cpp -I/opt/X11/include -I/usr/X11R6/include -I/usr/include -I/usr/local/include -I/usr/local/include/gsl -L/opt/X11/lib -L/usr/local/lib -lglfw3 -framework OpenGL -lglew -framework Cocoa -framework IOKit -framework CoreVideo -lgsl -lgslcblas -Wno-deprecated-declarations
