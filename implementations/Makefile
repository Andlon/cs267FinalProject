OBJ := obj
BIN := bin
INCLUDE := include

# Include modules here
include src/mpi-naive/module.mk

all: $(TARGETS)

.PHONY: clean
clean:
	rm -f $(TARGETS)
	rm -fR $(OBJ)
	rmdir --ignore-fail-on-non-empty $(BIN)
