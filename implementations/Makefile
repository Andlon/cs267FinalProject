OBJ := obj
BIN := bin
INCLUDE := include

LFLAGS = 
CXXFLAGS = -O3 -std=c++11
DEBUG_CXXFLAGS = -g -std=c++11 -O0
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
