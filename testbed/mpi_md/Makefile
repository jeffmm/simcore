# Maybe a generic makefile?

CC = gcc-5
CXX = g++-5
MPICC = mpicc
MPICXX = mpic++

TARGET := argon

SRCDIR = src
OBJDIR = obj
BINDIR = bin
SRCEXT = cpp

COMPILE_FLAGS = -std=c++11 -Wall -Wextra -fopenmp
RCOMPILE_FLAGS = -D NDEBUG -O3 -ffast-math -fexpensive-optimizations -funroll-all-loops -msse3 
DCOMPILE_FLAGS = -D DEBUG -g
LINK_FLAGS = -fopenmp
RLINK_FLAGS = 
DLINK_FLGAS =

INCLUDES = -I $(SRCDIR) -I/usr/local/include
LIBS = -L/usr/local/lib -lm -lboost_iostreams-mt -lboost_serialization-mt

UNAME_S:=$(shell uname -s)

print-%: ; echo $*=$($*)

SHELL = /bin/bash

.SUFFIXES:

ifneq ($(LIBS),)
	COMPILE_FLAGS += $(shell pkg-config --cflags $(LIBS))
	LINK_FLAGS += $(shell pkg-config --libs $(LIBS))
endif

# Combine compiler and linker flags
release: export CXXFLAGS := $(CXXFLAGS) $(COMPILE_FLAGS) $(RCOMPILE_FLAGS)
release: export LDFLAGS := $(LDFLAGS) $(LINK_FLAGS) $(RLINK_FLAGS)
debug: export CXXFLAGS := $(CXXFLAGS) $(COMPILE_FLAGS) $(DCOMPILE_FLAGS)
debug: export LDFLAGS := $(LDFLAGS) $(LINK_FLAGS) $(DLINK_FLAGS)

# find all the stupid files
ifeq ($(UNAME_S),Darwin)
	SOURCES = $(shell find $(SRCDIR) -name '*.$(SRCEXT)' | sort -k 1nr | cut -f2-)
else
	SOURCES = $(shell find $(SRCDIR) -name '*.$(SRCEXT)' -printf '%T@\t%p\n' \
			  								| sort -k 1nr | cut -f2-)
endif

# fallback
rwildcard = $(foreach d, $(wildcard $1*), $(call rwildcard,$d/,$2) \
											$(filter $(subst *,%,$2),$d))

ifeq ($(SOURCES),)
	SOURCES := $(call rwildcard, $(SRCDIR), *.$(SRCEXT))
endif

# Set the appropriate object filenames
OBJECTS = $(SOURCES:$(SRCDIR)/%.$(SRCEXT)=$(OBJDIR)/%.o)
DEPS = $(OBJECTS:.o=.d)

.PHONY: release
release: dirs
	@echo "Beginning release build"
	@$(MAKE) all --no-print-directory

.PHONY: debug
debug: dirs
	@echo "Beginning debug build"
	@$(MAKE) all --no-print-directory

.PHONY: dirs
dirs:
	mkdir -p $(dir $(OBJECTS))
	mkdir -p $(BINDIR)

.PHONY: clean
clean:
	$(RM) -r $(OBJDIR)
	$(RM) -r $(BINDIR)

# the whole shebang
all: $(BINDIR)/$(TARGET)
	@echo "Making all"

$(BINDIR)/$(TARGET): $(OBJECTS)
	$(MPICXX) $(OBJECTS) $(LDFLAGS) $(LIBS) -o $@

# add depedencies
-include $(DEPS)

# source file ruels
$(OBJDIR)/%.o: $(SRCDIR)/%.$(SRCEXT)
	$(MPICXX) $(CXXFLAGS) $(INCLUDES) -MP -MMD -c $< -o $@

