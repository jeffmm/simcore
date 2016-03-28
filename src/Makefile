CC = g++-5
CXX = g++-5
CFLAGS = -I/opt/X11/include -I/usr/X11R6/include -I/usr/include -I/usr/local/include -g -DDEBUG -std=c++11

LNKFLAGS = -lyaml-cpp

ifeq ($(CC),icpc)
    OMPFLAGS = -openmp -DBOB_OMP -Wno-deprecated
else
    CFLAGS += -Wno-deprecated-declarations -Wno-deprecated
    OMPFLAGS =	-fopenmp -DBOB_OMP
endif
# for Ubuntu distros
#GLXFLAGS = -L/opt/X11/lib -lGL -lGLU -lGLEW -lglfw3 -lX11 -lXinerama -lXcursor -lXxf86vm -lXrandr -lpthread -lXi
GLXFLAGS = -L/opt/X11/lib -lglfw3 -framework OpenGL -lglew -lpthread
GSLFLAGS = -L/usr/local/lib/ -I/usr/local/include/gsl/ -lgsl -lgslcblas
FFTFLAGS = -L/usr/lib64 -lfftw3

CXX = $(CC)
CXXFLAGS = $(CFLAGS)

../configure_parameters:	configure_parameters.o
	$(CC) -o ../configure_parameters configure_parameters.o

CYTOSCORE_OBJS = cytoscore_main.o simulation_manager.o allocate.o vector_algebra.o cpu.o error_exit.o generate_random_unit_vector.o

../cytoscore:	$(CYTOSCORE_OBJS) $(AUX_OBJS)
	$(CC)  -o ../cytoscore $(CYTOSCORE_OBJS) $(GLXFLAGS) $(GSLFLAGS) -lm $(LNKFLAGS)

