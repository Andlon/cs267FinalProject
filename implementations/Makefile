OBJ := obj
BIN := bin
INCLUDE := include

# Set default MPI compiler
MPICXX = CC

LFLAGS = 
CXXFLAGS = -O3 -std=c++1y -Wno-literal-suffix
DEBUG_CXXFLAGS = -g -std=c++1y -O0
DEBUG_LFLAGS = 

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
