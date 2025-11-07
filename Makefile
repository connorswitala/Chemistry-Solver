# ===== Makefile =====
# Compiler and flags
CXX := g++
CXXFLAGS := -std=c++17 -O3 -I./includes -I./solverLib -I./common_mixes

# Directories
COMMON_AIR_DIR 	:= common_mixes
SOLVER_DIR     	:= CESolver
BUILD_DIR      	:= build
BIN_DIR		   	:= bin
SRC_DIR			:= source

# Files
COMMON_AIR_SRC 	:= $(COMMON_AIR_DIR)/commonMixes.cpp
SOLVER_SRC     	:= $(SOLVER_DIR)/CESolver.cpp
MIN_SRC			:= $(SRC_DIR)/minimize.cpp
PLOT_SRC		:= $(SRC_DIR)/plot_composition.cpp
TIMING_SRC		:= $(SRC_DIR)/timing.cpp
# Objects
COMMON_AIR_OBJ 	:= $(BUILD_DIR)/commonMixes.o
SOLVER_OBJ     	:= $(BUILD_DIR)/CESolver.o
MIN_OBJ   		:= $(BUILD_DIR)/minimize.o
PLOT_OBJ   		:= $(BUILD_DIR)/plot.o
TIMING_OBJ   	:= $(BUILD_DIR)/timing.o

# Targets
MIN_TARGET := $(BIN_DIR)/minimize
PLOT_TARGET := $(BIN_DIR)/plot
TIMING_TARGET := $(BIN_DIR)/timing

# ===== Targets =====
all: $(MIN_TARGET) $(PLOT_TARGET)

# Create build directory if missing
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

# Compile source files
$(COMMON_AIR_OBJ): $(COMMON_AIR_SRC) | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(SOLVER_OBJ): $(SOLVER_SRC) | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(MIN_OBJ): $(MIN_SRC) | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(PLOT_OBJ): $(PLOT_SRC) | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(TIMING_OBJ): $(TIMING_SRC) | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Link final executable
$(MIN_TARGET): $(COMMON_AIR_OBJ) $(SOLVER_OBJ) $(MIN_OBJ)
	$(CXX) $(CXXFLAGS) $^ -o $@

$(PLOT_TARGET): $(COMMON_AIR_OBJ) $(SOLVER_OBJ) $(PLOT_OBJ)
	$(CXX) $(CXXFLAGS) $^ -o $@

$(TIMING_TARGET): $(COMMON_AIR_OBJ) $(SOLVER_OBJ) $(TIMING_OBJ)
	$(CXX) $(CXXFLAGS) $^ -o $@

# Clean up
clean:
	rm -rf $(BUILD_DIR) $(TEST_TARGET)

.PHONY: all clean
