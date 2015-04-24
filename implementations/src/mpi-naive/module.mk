# Get relative directory of this module.mk
DIR := $(dir $(lastword $(MAKEFILE_LIST)))
MODULE := $(BIN)/mpi-naive
OBJDIR := $(OBJ)/$(DIR)

# For now just add all headers as dependency of every object.
# Not optimal, but not a problem for such a small project.
LOCAL_HEADERS := $(DIR)/data.h $(DIR)/variogram.h

# Compile object files
$(OBJDIR)/%.o:: $(DIR)/%.cpp $(LOCAL_HEADERS) $(HEADERS)
	mkdir -p $(OBJDIR)
	$(MPICXX) $(CXXFLAGS) -c $(INCLUDE) $< -o $@

# Compile targets
$(MODULE): $(OBJDIR)/main.o $(OBJDIR)/data.o $(OBJDIR)/variogram.o
	mkdir -p $(BIN)
	$(MPICXX) $^ -o $@	$(LFLAGS)

TARGETS += $(MODULE)
