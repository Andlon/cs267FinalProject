# Get relative directory of this module.mk
DIR := $(dir $(lastword $(MAKEFILE_LIST)))
MODULE := $(BIN)/mpi-naive
OBJDIR := $(OBJ)/$(DIR)

$(OBJDIR)/%.o: $(DIR)/%.cpp
	mkdir -p $(OBJ)/$(DIR)
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDE) $< -o $@

$(MODULE): $(OBJDIR)/main.o $(OBJDIR)/data.o
	mkdir -p $(BIN)
	$(CXX) $^ -o $@	$(LFLAGS)

TARGETS += $(MODULE)
