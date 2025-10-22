# ===== Makefile =====
# Compiler and flags
CXX := g++
CXXFLAGS := -std=c++17 -O2 -I./includes -I./solverLib -I./common_mixes

# Directories
COMMON_AIR_DIR 	:= common_mixes
SOLVER_DIR     	:= equilibrium
BUILD_DIR      	:= build
BIN_DIR		   	:= bin
SRC_DIR			:= source

# Files
COMMON_AIR_SRC 	:= $(COMMON_AIR_DIR)/commonMixes.cpp
SOLVER_SRC     	:= $(SOLVER_DIR)/equilibrium.cpp
GIBBS_SRC      	:= $(SRC_DIR)/gibbs.cpp
PLOT_SRC		:= $(SRC_DIR)/plot.cpp

COMMON_AIR_OBJ := $(BUILD_DIR)/commonMixes.o
SOLVER_OBJ     := $(BUILD_DIR)/equilibrium.o
GIBBS_OBJ       := $(BUILD_DIR)/gibbs.o
PLOT_OBJ		:= $(BUILD_DIR)/plot.o

GIBBS_TARGET    := $(BIN_DIR)/gibbs
PLOT_TARGET		:= $(BIN_DIR)/plot

# ===== Targets =====
all: $(GIBBS_TARGET) $(PLOT_TARGET)

# Create build directory if missing
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

# Compile source files
$(COMMON_AIR_OBJ): $(COMMON_AIR_SRC) | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(SOLVER_OBJ): $(SOLVER_SRC) | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(GIBBS_OBJ): $(GIBBS_SRC) | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(PLOT_OBJ): $(PLOT_SRC) | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Link final executable
$(GIBBS_TARGET): $(COMMON_AIR_OBJ) $(SOLVER_OBJ) $(GIBBS_OBJ)
	$(CXX) $(CXXFLAGS) $^ -o $@

$(PLOT_TARGET): $(COMMON_AIR_OBJ) $(SOLVER_OBJ) $(PLOT_OBJ)
	$(CXX) $(CXXFLAGS) $^ -o $@

# Clean up
clean:
	rm -rf $(BUILD_DIR) $(TEST_TARGET)

.PHONY: all clean
