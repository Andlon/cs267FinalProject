OBJ := obj
BIN := bin
INCLUDE := -Iinclude -Isrc/include
TARGETS := 

# Set default MPI compiler
MPICXX = CC

LFLAGS = 
CXXFLAGS = -O3 -std=c++1y -Wno-literal-suffix
DEBUG_CXXFLAGS = -g -std=c++1y -O0
DEBUG_LFLAGS =

# Add our own headers common to multiple modules here.
# Ensures that the modules are rebuilt when we change these headers
HEADERS = src/include/pair_index_set.h src/include/mpi_util.h src/include/rect.h src/include/uniform_distribution.h

# Need to put 'all' target before includes so that it becomes the default target.
# Moreover, it depends on main-build, which uses a later definition of $(TARGETS).
all: main-build

# Include modules here
include src/mpi-naive/module.mk
include src/mpi-ca-all-pairs/module.mk
include src/mpi-all-unique-pairs/module.mk
include src/index-set-tool/module.mk

main-build: $(TARGETS)

debug: all
debug: CXXFLAGS = $(DEBUG_CXXFLAGS)
debug: LFLAGS = $(DEBUG_LFLAGS)

.PHONY: clean
clean:
	rm -f $(TARGETS)
	rm -fR $(OBJ)
	rmdir --ignore-fail-on-non-empty $(BIN)

# The following rule lets us print any variable by calling `make print-VARIABLE`. Useful for debuggin
print-%  : ; @echo $* = $($*)
