OBJ := obj
BIN := bin
INCLUDE := -Iinclude -Isrc/include

# Set default MPI compiler
MPICXX = mpic++

LFLAGS = 
CXXFLAGS = -O3 -std=c++1y -Wno-literal-suffix
DEBUG_CXXFLAGS = -g -std=c++1y -O0
DEBUG_LFLAGS =

# Add our own headers common to multiple modules here.
# Ensures that the modules are rebuilt when we change these headers
HEADERS = src/include/pair_index_set.h 

# Include modules here
include src/mpi-naive/module.mk

all: $(TARGETS)

debug: all
debug: CXXFLAGS = $(DEBUG_CXXFLAGS)
debug: LFLAGS = $(DEBUG_LFLAGS)

.PHONY: clean
clean:
	rm -f $(TARGETS)
	rm -fR $(OBJ)
	rmdir --ignore-fail-on-non-empty $(BIN)
