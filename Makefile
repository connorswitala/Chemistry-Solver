# ===== Makefile =====
# Compiler and flags
CXX := g++
CXXFLAGS := -std=c++17 -O3 -I./includes -I./solverLib -I./common_mixes -I/readNASA.h

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
PLOT_SRC		:= $(SRC_DIR)/plot_composition_tec.cpp
CSV_SRC			:= $(SRC_DIR)/plot_composition_csv.cpp
TIMING_SRC		:= $(SRC_DIR)/timing.cpp
TEST_SRC		:= $(SRC_DIR)/test.cpp

# Objects
COMMON_AIR_OBJ 	:= $(BUILD_DIR)/commonMixes.o
SOLVER_OBJ     	:= $(BUILD_DIR)/CESolver.o
MIN_OBJ   		:= $(BUILD_DIR)/minimize.o
PLOT_OBJ   		:= $(BUILD_DIR)/plot.o
CSV_OBJ   		:= $(BUILD_DIR)/csv.o
TIMING_OBJ   	:= $(BUILD_DIR)/timing.o
TEST_OBJ   		:= $(BUILD_DIR)/test.o

# Targets
MIN_TARGET := $(BIN_DIR)/minimize
PLOT_TARGET := $(BIN_DIR)/tec
CSV_TARGET := $(BIN_DIR)/csv
TIMING_TARGET := $(BIN_DIR)/timing
TEST_TARGET := $(BIN_DIR)/test

# ===== Targets =====
all: $(MIN_TARGET) $(PLOT_TARGET) $(CSV_TARGET) $(TIMING_TARGET) $(TEST_TARGET)

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

$(CSV_OBJ): $(CSV_SRC) | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(TIMING_OBJ): $(TIMING_SRC) | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(TEST_OBJ): $(TEST_SRC) | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Link final executable
$(MIN_TARGET): $(COMMON_AIR_OBJ) $(SOLVER_OBJ) $(MIN_OBJ)
	$(CXX) $(CXXFLAGS) $^ -o $@

$(PLOT_TARGET): $(COMMON_AIR_OBJ) $(SOLVER_OBJ) $(PLOT_OBJ)
	$(CXX) $(CXXFLAGS) $^ -o $@

$(CSV_TARGET): $(COMMON_AIR_OBJ) $(SOLVER_OBJ) $(CSV_OBJ)
	$(CXX) $(CXXFLAGS) $^ -o $@

$(TIMING_TARGET): $(COMMON_AIR_OBJ) $(SOLVER_OBJ) $(TIMING_OBJ)
	$(CXX) $(CXXFLAGS) $^ -o $@

$(TEST_TARGET): $(COMMON_AIR_OBJ) $(SOLVER_OBJ) $(TEST_OBJ)
	$(CXX) $(CXXFLAGS) $^ -o $@

# Clean up
clean:
	rm -rf $(BUILD_DIR) 

.PHONY: all clean
