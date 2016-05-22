CC = g++-5
CXX = g++-5
CFLAGS = -I/opt/X11/include -I/usr/X11R6/include -I/usr/include -I/usr/local/include -g -DDEBUG -std=c++11
LNKFLAGS = -lyaml-cpp -gnu

ifeq ($(CC),icpc)
	OMPFLAGS = -openmp -DBOB_OMP -Wno-deprecated
else
	CFLAGS += -Wno-deprecated-declarations -Wno-deprecated
	OMPFLAGS =	-fopenmp -DBOB_OMP
endif

GLXFLAGS = -L/opt/X11/lib -lglfw3 -framework OpenGL -lglew -lpthread
GSLFLAGS = -L/usr/local/lib/ -I/usr/local/include/gsl/ -lgsl -lgslcblas
FFTFLAGS = -L/usr/lib64 -lfftw3
CXX = $(CC)
CXXFLAGS = $(CFLAGS)

CYTOSCORE_OBJS = cytoscore_main.o simulation_manager.o simulation.o space.o allocate.o vector_algebra.o cpu.o error_exit.o generate_random_unit_vector.o graphics.o grabber.o writebmp.o

CONFIGURE_OBJS = configure_cytoscore.o

CYTOSCORE_OBJECT_OBJS = brownian_dimer.o

../cytoscore:	$(CYTOSCORE_OBJS) $(CYTOSCORE_OBJECT_OBJS)
	$(CC)  -o ../cytoscore $(CYTOSCORE_OBJS) $(CYTOSCORE_OBJECT_OBJS) $(GLXFLAGS) $(GSLFLAGS) -lm $(LNKFLAGS) 

../configure_cytoscore: $(CONFIGURE_OBJS)
	$(CC) -o ../configure_cytoscore $(CONFIGURE_OBJS) $(GLXFLAGS) $(GSLFLAGS) -lm $(LNKFLAGS)
