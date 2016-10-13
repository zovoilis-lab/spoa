CP = g++
LD = g++
DX = doxygen

NAME = spoa

OBJ_DIR = obj
SRC_DIR = src
VND_DIR = vendor
INC_DIR = include/$(NAME)
LIB_DIR = lib
DOC_DIR = doc
EXC_DIR = bin

I_CMD = $(addprefix -I, $(SRC_DIR) $(VND_DIR))
L_CMD = $(addprefix -L, )

CP_FLAGS = $(I_CMD) -O3 -Wall -std=c++11 -march=native -DSPOA_MAIN_
LD_FLAGS = $(I_CMD) $(L_CMD)

API = $(addprefix $(SRC_DIR)/, alignment.hpp edge.hpp graph.hpp node.hpp \
		simd_alignment.hpp sisd_alignment.hpp spoa.hpp)

SRC = $(shell find $(SRC_DIR) -type f -regex ".*\.cpp")
OBJ = $(subst $(SRC_DIR), $(OBJ_DIR), $(addsuffix .o, $(basename $(SRC))))
DEP = $(OBJ:.o=.d)
INC = $(subst $(SRC_DIR), $(INC_DIR), $(API))
LIB = $(LIB_DIR)/lib$(NAME).a
EXC = $(NAME)
BIN = $(EXC_DIR)/$(EXC)
DOC = $(DOC_DIR)/DoxyFile

all: $(EXC)

install: lib include

include: $(INC)

lib: $(LIB)

$(EXC): $(OBJ)
	@echo [LD] $@
	@mkdir -p $(dir $@)
	@$(LD) -o $@ $^ $(LD_FLAGS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	@echo [CP] $<
	@mkdir -p $(dir $@)
	@$(CP) $< -c -o $@ -MMD $(CP_FLAGS)

$(INC_DIR)/%.hpp: $(SRC_DIR)/%.hpp
	@echo [CP] $@
	@mkdir -p $(dir $@)
	@cp $< $@

$(LIB): $(OBJ)
	@echo [AR] $@
	@mkdir -p $(dir $@)
	@ar rcs $(LIB) $(OBJ)

docs:
	@echo [DX] generating documentation
	@$(DX) $(DOC)

clean:
	@echo [RM] cleaning
	@rm $(OBJ_DIR) $(EXC) -rf

remove:
	@echo [RM] removing
	@rm $(OBJ_DIR) $(EXC) $(LIB) $(INC_DIR) -rf

-include $(DEP)
