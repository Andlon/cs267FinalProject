OBJ := obj
BIN := bin
INCLUDE := include

LFLAGS := 
CXXFLAGS := -O3 -std=c++11

# Include modules here
include src/mpi-naive/module.mk

all: $(TARGETS)

.PHONY: clean
clean:
	rm -f $(TARGETS)
	rm -fR $(OBJ)
	rmdir --ignore-fail-on-non-empty $(BIN)
