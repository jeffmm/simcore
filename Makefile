CC = gcc-5
CXX = g++-5

SRCDIR = src
OBJDIR = obj
BINDIR = bin
SRCEXT = cpp

COMPILE_FLAGS = -std=c++11
RCOMPILE_FLAGS = -D NDEBUG -O2 -march=native -g
DCOMPILE_FLAGS = -D DEBUG -O0 -g
LINK_FLAGS = -gnu 
RLINK_FLAGS = 
DLINK_FLAGS =

ifeq ($(THREADING),eomp)
	COMPILE_FLAGS += -fopenmp
	LINK_FLAGS += -fopenmp
endif

ifeq ($(NOGRAPH),true)
	COMPILE_FLAGS += -D NOGRAPH
	LINK_FLAGS += -D NOGRAPH
endif

UNAME_S:=$(shell uname -s)

ifeq ($(UNAME_S),Darwin)
	INCLUDES = -I/opt/X11/include -I/usr/X11R6/include -I/usr/include -I/usr/local/include -I/usr/local/include/gsl -I/usr/local/include/freetype2
	GLXLIBS = -L/opt/X11/lib -L/usr/local/lib -lglfw3 -framework OpenGL -lglew -lfreetype
	GSLLIBS = -lgsl -lgslcblas
	FFTLIBS = -L/usr/lib64 -lfftw3
	YAMLLIBS = -L/Users/$(USER)/Projects/Lib/yaml-cpp-gcc5/yaml-cpp/build -lyaml-cpp
	LIBS = $(GLXLIBS) $(GSLLIBS) $(FFTLIBS) $(YAMLLIBS) -L/usr/local/lib
else
	GSLINCS = -I/usr/local/include
	GSLLIBS = -L/usr/local/lib -lgsl -lgslcblas -lm
	GLFW3INCS = -I/home/cedelmaier/common/glfw/include
	GLFW3LIBS = -L/home/cedelmaier/common/glfw/build/src -lglfw3 -lGLEW -lGLU -lGL -lX11 -lXxf86vm -lpthread -ldl -lXrandr -lXi -lXcursor -lXinerama
	YAMLINCS = -I/home/cedelmaier/common/yaml-cpp/include
	#YAMLLIBS = -L/home/cedelmaier/common/yaml-cpp/build -lyaml-cpp
	YAMLLIBS = -Wl,-rpath,/home/cedelmaier/common/yaml-cpp/build -L/home/cedelmaier/common/yaml-cpp/build -lyaml-cpp
	INCLUDES = $(GLFW3INCS) $(YAMLINCS) $(GSLINCS)
	LIBS = $(GLFW3LIBS) $(YAMLLIBS) $(GSLLIBS)
	#export LD_LIBRARY_PATH=/home/cedelmaier/common/glfw/build/src:/home/cedelmaier/common/yaml-cpp/build:$LD_LIBRARY_PATH
endif

print-%: ; @echo $*=$($*)

SHELL = /bin/bash

.SUFFIXES:

ifneq ($(LIBS),)
	COMPILE_FLAGS += $(shell pkg-config --cflags $(LIBS))
	LINK_FLAGS += $(shell pkg-config --libs $(LIBS))
endif

# Special stuff for intel compiler
CC=$(CXX)
ifeq ($(CC),icpc)
	COMPILE_FLAGS += -Wno-deprecated
	#RCOMPILE_FLAGS += -openmp -DBOB_OMP
else
	COMPILE_FLAGS += -Wno-deprecated-declarations -Wno-deprecated
	#RCOMPILE_FLAGS += -fopenmp -DBOB_OMP
endif

# Combine compiler and linker flags
ifeq ($(CFG),release)
	export CXXFLAGS := $(CXXFLAGS) $(COMPILE_FLAGS) $(RCOMPILE_FLAGS)
	export LDFLAGS := $(LDFLAGS) $(LINK_FLAGS) $(RLINK_FLAGS)
else
	export CXXFLAGS := $(CXXFLAGS) $(COMPILE_FLAGS) $(DCOMPILE_FLAGS)
	export LDFLAGS := $(LDFLAGS) $(LINK_FLAGS) $(DLINK_FLAGS)
endif

# build information on all sources
ifeq ($(UNAME_S),Darwin)
	SOURCES = $(shell find $(SRCDIR) -name '*.$(SRCEXT)' | sort -k 1nr | cut -f2-)
else
	SOURCES = $(shell find $(SRCDIR) -name '*.$(SRCEXT)' -printf '%T@\t%p\n' \
			  										| sort -k 1nr | cut -f2-)
endif

# fallback case
rwildcard = $(foreach d, $(wildcard $1*), $(call rwildcard,$d/,$2) \
									$(filter $(subst *,%,$2), $d))
ifeq ($(SOURCES),)
	SOURCES := $(call rwildcard, $(SRCDIR), *.$(SRCEXT))
endif

# Now we have to figure out which we are building of the program list, since that matters
# for things like not compiling more than one main
SIMCORE_SOURCE = $(SRCDIR)/simcore_main.cpp
CONFIGURE_SIMCORE_SOURCE = $(SRCDIR)/configure_simcore.cpp
EXCLUDE_SOURCES = $(SRCDIR)/integrator_manager.cpp $(SRCDIR)/make_params.cpp 

MAIN_SOURCES = $(SIMCORE_SOURCE) $(CONFIGURE_SIMCORE_SOURCE)

# These are the common sources
SRCS = $(filter-out $(MAIN_SOURCES) $(EXCLUDE_SOURCES), $(SOURCES))

OBJECTS = $(SRCS:$(SRCDIR)/%.$(SRCEXT)=$(OBJDIR)/%.o)
SIMCORE_MAIN_OBJ = $(SIMCORE_SOURCE:$(SRCDIR)/%.$(SRCEXT)=$(OBJDIR)/%.o)
CONFIGURE_SIMCORE_OBJ = $(CONFIGURE_SIMCORE_SOURCE:$(SRCDIR)/%.$(SRCEXT)=$(OBJDIR)/%.o)
DEPS = $(OBJECTS:.o=.d)

.PHONY: dirs
dirs:
	mkdir -p $(OBJDIR)
	mkdir -p $(BINDIR)

.PHONY: clean
clean:
	$(RM) -r $(OBJDIR)
	$(RM) -r $(BINDIR)

.PHONY : clean-output
clean-output :
	rm -f *.posit *.config *.thermo *.initial_config.* *.final_config *.final_config.* \
	*.checkpoint.* *.checkpoint *.thermo_ext sphero.crosslinks.* test-log

simcore: dirs $(BINDIR)/simcore;cp $(BINDIR)/simcore simcore

configure_simcore: dirs $(BINDIR)/configure_simcore;cp $(BINDIR)/configure_simcore configure_simcore

$(BINDIR)/simcore: $(OBJECTS) $(SIMCORE_MAIN_OBJ)
	$(CXX) $^ -o $@ $(LDFLAGS) $(LIBS) 

$(BINDIR)/configure_simcore: $(OBJECTS) $(CONFIGURE_SIMCORE_OBJ)
	$(CXX) $(LDFLAGS) $(LIBS) $^ -o $@

# add dependencies
-include $(DEPS)


# source file rules
$(OBJDIR)/%.o: $(SRCDIR)/%.$(SRCEXT)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -MP -MMD -c $< -o $@ 
