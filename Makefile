# ===== Makefile =====
# Compiler and flags
CXX := g++
CXXFLAGS := -std=c++17 -O3 -I./libraries/includes -I./libraries/solverLib -I./libraries/mixes -I/libraries/readNASA.h

# Directories
MIXES_DIR 		:= libraries/mixes
SOLVER_DIR     	:= libraries/CESolver
BUILD_DIR      	:= misc/build
BIN_DIR		   	:= bin
SRC_DIR			:= source

# Files
MIXES_SRC 		:= $(MIXES_DIR)/mixes.cpp
SOLVER_SRC     	:= $(SOLVER_DIR)/CESolver.cpp
MIN_SRC			:= $(SRC_DIR)/minimize.cpp
TIMING_SRC		:= $(SRC_DIR)/timing.cpp
EXAMPLE_SRC		:= $(SRC_DIR)/example.cpp

# Objects
MIXES_OBJ	 	:= $(BUILD_DIR)/mixes.o
SOLVER_OBJ     	:= $(BUILD_DIR)/CESolver.o
MIN_OBJ   		:= $(BUILD_DIR)/minimize.o
TIMING_OBJ   	:= $(BUILD_DIR)/timing.o
EXAMPLE_OBJ		:= $(BUILD_DIR)/example.o

# Targets
MIN_TARGET := $(BIN_DIR)/minimize
TIMING_TARGET := $(BIN_DIR)/timing
EXAMPLE_TARGET := $(BIN_DIR)/example


# ===== Targets =====
all: $(MIN_TARGET) $(TIMING_TARGET) $(EXAMPLE_TARGET)

# Create build directory if missing
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

# Compile source files
$(MIXES_OBJ): $(MIXES_SRC) | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(SOLVER_OBJ): $(SOLVER_SRC) | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(MIN_OBJ): $(MIN_SRC) | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(EXAMPLE_OBJ): $(EXAMPLE_SRC) | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(TIMING_OBJ): $(TIMING_SRC) | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Link final executable
$(MIN_TARGET): $(MIXES_OBJ) $(SOLVER_OBJ) $(MIN_OBJ)
	$(CXX) $(CXXFLAGS) $^ -o $@

$(EXAMPLE_TARGET): $(MIXES_OBJ) $(SOLVER_OBJ) $(EXAMPLE_OBJ)
	$(CXX) $(CXXFLAGS) $^ -o $@

$(TIMING_TARGET): $(MIXES_OBJ) $(SOLVER_OBJ) $(TIMING_OBJ)
	$(CXX) $(CXXFLAGS) $^ -o $@


# Clean up
clean:
	rm -rf $(BUILD_DIR) 

.PHONY: all clean
