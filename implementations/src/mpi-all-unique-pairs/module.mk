# Get relative directory of this module.mk
DIR := $(dir $(lastword $(MAKEFILE_LIST)))
MODULE := $(BIN)/mpi-all-unique-pairs
OBJDIR := $(OBJ)/$(DIR)

# For now just add all headers as dependency of every object.
# Not optimal, but not a problem for such a small project.
LOCAL_HEADERS := $(DIR)/data.h $(DIR)/variogram.h $(DIR)/uniform_distribution.h $(DIR)/column_team.h $(DIR)/node_grid.h

# Compile object files
$(OBJDIR)/%.o:: $(DIR)/%.cpp $(LOCAL_HEADERS) $(HEADERS)
	mkdir -p `dirname $@`
	$(MPICXX) $(CXXFLAGS) -c $(INCLUDE) $< -o $@

# Compile targets
$(MODULE): $(OBJDIR)/main.o $(OBJDIR)/data.o $(OBJDIR)/variogram.o $(OBJDIR)/column_team.o $(OBJDIR)/node_grid.o
	mkdir -p $(BIN)
	$(MPICXX) $^ -o $@ $(LFLAGS)

mpi-all-unique-pairs: $(MODULE)

TARGETS += $(MODULE)
