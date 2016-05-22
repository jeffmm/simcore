#Generic Makefile for animation projects

TARGET	= ArNe_sim

#Compilers and flags
CXX		= g++-5

#Compiler flags
CXXFLAGS	= -g -Wall -Wextra
CXXFLAGS	+= -I/usr/local/include
ifeq ($(CFG),test)
	CXXFLAGS	+= -g -pg -O0
else
	CXXFLAGS	+= -O3 -ffast-math -fexpensive-optimizations -funroll-all-loops -msse3
endif

ifeq ($(USETHREADS),yes)
	OPENMP		= -fopenmp
else
	OPENMP		=
endif
THREADS		= -std=c++11

LINKER	= g++-5 -o
#Linking flags
LFLAGS	= -Wall -I. -lm -I/usr/local/include
GLFW3FLAGS = `pkg-config --static --libs glfw3`

#Set directories for src, objs, and bin
SRCDIR	= src
OBJDIR	= obj
BINDIR	= bin

SOURCES		:= $(wildcard $(SRCDIR)/*.cpp)
INCLUDES	:= $(wildcard $(SRCDIR)/*.h)
OBJECTS		:= $(SOURCES:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o)
rm	= rm -f

$(BINDIR)/$(TARGET): $(OBJECTS)
	$(LINKER) $@ $(LFLAGS) $(GLFW3FLAGS) $(OPENMP) $(OBJECTS) $(THREADS)
	@echo "Linking complete!"

$(OBJECTS): $(OBJDIR)/%.o : $(SRCDIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(THREADS) $(OPENMP) -c $< -o $@
	@echo "Compiled "$<" successfully!"

.PHONEY: clean
clean:
	$(rm) $(OBJECTS)
	@echo "Cleanup complete!"

.PHONEY: remove
remove: clean
	$(rm) $(BINDIR)/$(TARGET)
	@echo "Executable removed!"

