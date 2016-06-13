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

SIMCORE_OBJS = simcore_main.o simulation_manager.o simulation.o object.o forces.o cell_list.o space.o allocate.o vector_algebra.o cpu.o error_exit.o generate_random_unit_vector.o graphics.o grabber.o writebmp.o minimum_distance.o

CONFIGURE_OBJS = configure_simcore.o

SIMCORE_OBJECT_OBJS = md_bead.o br_dimer.o br_bead.o argon.o neon.o

../simcore:	$(SIMCORE_OBJS) $(SIMCORE_OBJECT_OBJS)
	$(CC)  -o ../simcore $(SIMCORE_OBJS) $(SIMCORE_OBJECT_OBJS) $(GLXFLAGS) $(GSLFLAGS) -lm $(LNKFLAGS) 

../configure_simcore: $(CONFIGURE_OBJS)
	$(CC) -o ../configure_simcore $(CONFIGURE_OBJS) $(GLXFLAGS) $(GSLFLAGS) -lm $(LNKFLAGS)
