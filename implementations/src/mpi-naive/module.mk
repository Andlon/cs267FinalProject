# Get relative directory of this module.mk
DIR := $(dir $(lastword $(MAKEFILE_LIST)))
MODULE := $(BIN)/mpi-naive
OBJDIR := $(OBJ)/$(DIR)

# For now just add all headers as dependency of every object.
# Not optimal, but not a problem for such a small project.
HEADERS := $(DIR)/pair_index_set.h $(DIR)/data.h $(DIR)/variogram.h

# Compile object files
$(OBJDIR)/%.o:: $(DIR)/%.cpp $(HEADERS)
	mkdir -p $(OBJDIR)
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDE) $< -o $@

# Compile targets
$(MODULE): $(OBJDIR)/main.o $(OBJDIR)/data.o $(OBJDIR)/variogram.o
	mkdir -p $(BIN)
	$(CXX) $^ -o $@	$(LFLAGS)

TARGETS += $(MODULE)
