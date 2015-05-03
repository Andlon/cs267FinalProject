# Get relative directory of this module.mk
DIR := $(dir $(lastword $(MAKEFILE_LIST)))
MODULE := $(BIN)/mpi-ca-all-pairs
OBJDIR := $(OBJ)/$(DIR)

# For now just add all headers as dependency of every object.
# Not optimal, but not a problem for such a small project.
LOCAL_HEADERS := $(DIR)/data.h $(DIR)/variogram.h $(DIR)/uniform_distribution.h $(DIR)/team.h $(DIR)/grid.h

# Compile object files
$(OBJDIR)/%.o:: $(DIR)/%.cpp $(LOCAL_HEADERS) $(HEADERS)
	mkdir -p `dirname $@`
	$(MPICXX) $(CXXFLAGS) -c $(INCLUDE) $< -o $@

# Compile targets
$(MODULE): $(OBJDIR)/main.o $(OBJDIR)/data.o $(OBJDIR)/variogram.o $(OBJDIR)/team.o $(OBJDIR)/grid.o
	mkdir -p $(BIN)
	$(MPICXX) $^ -o $@ $(LFLAGS)

mpi-ca-all-pairs: $(MODULE)

TARGETS += $(MODULE)
