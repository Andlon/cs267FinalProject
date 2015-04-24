# Get relative directory of this module.mk
DIR := $(dir $(lastword $(MAKEFILE_LIST)))
MODULE := $(BIN)/index-set-tool
OBJDIR := $(OBJ)/$(DIR)

# For now just add all headers as dependency of every object.
# Not optimal, but not a problem for such a small project.
LOCAL_HEADERS :=

# Compile object files
$(OBJDIR)/%.o:: $(DIR)/%.cpp $(LOCAL_HEADERS) $(HEADERS)
	mkdir -p `dirname $@`
	$(CXX) $(CXXFLAGS) -c $(INCLUDE) $< -o $@

# Compile targets
$(MODULE): $(OBJDIR)/main.o
	mkdir -p $(BIN)
	$(CXX) $^ -o $@	$(LFLAGS)

# Add custom target
index-set-tool: $(MODULE)

TARGETS += $(MODULE)
