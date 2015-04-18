OBJ := obj
BIN := bin
INCLUDE := include

# Include modules here
include mpi-naive/module.mk

all: $(TARGETS)

.PHONY: clean
clean:
	rm -f $(BIN)/$(TARGETS)
	rm -R $(OBJ)
