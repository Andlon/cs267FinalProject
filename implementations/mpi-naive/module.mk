# Get relative directory of this module.mk
DIR := $(dir $(lastword $(MAKEFILE_LIST)))
MODULE := $(BIN)/mpi-naive

$(OBJ)/$(DIR)/%.o: $(DIR)/%.cpp
	mkdir -p $(OBJ)/$(DIR)
	$(CXX) -c -I$(INCLUDE) $< -o $@

$(MODULE): $(OBJ)/$(DIR)/main.o
	mkdir -p $(BIN)
	$(CXX) $^ -o $@	

TARGETS += $(MODULE)
