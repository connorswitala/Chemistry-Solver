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
PLOT_SRC		:= $(SRC_DIR)/plot_composition_tec.cpp
CSV_SRC			:= $(SRC_DIR)/plot_composition_csv.cpp
TIMING_SRC		:= $(SRC_DIR)/timing.cpp
TEST_SRC		:= $(SRC_DIR)/test.cpp
PRINTNASA_SRC   := $(SRC_DIR)/printnasa.cpp

# Objects
MIXES_OBJ	 	:= $(BUILD_DIR)/mixes.o
SOLVER_OBJ     	:= $(BUILD_DIR)/CESolver.o
MIN_OBJ   		:= $(BUILD_DIR)/minimize.o
PLOT_OBJ   		:= $(BUILD_DIR)/plot.o
CSV_OBJ   		:= $(BUILD_DIR)/csv.o
TIMING_OBJ   	:= $(BUILD_DIR)/timing.o
TEST_OBJ   		:= $(BUILD_DIR)/test.o
PRINTNASA_OBJ   := $(BUILD_DIR)/printnasa.o

# Targets
MIN_TARGET := $(BIN_DIR)/minimize
PLOT_TARGET := $(BIN_DIR)/tec
CSV_TARGET := $(BIN_DIR)/csv
TIMING_TARGET := $(BIN_DIR)/timing
TEST_TARGET := $(BIN_DIR)/test
PRINTNASA_TARGET := $(BIN_DIR)/printnasa


# ===== Targets =====
all: $(MIN_TARGET) $(PLOT_TARGET) $(CSV_TARGET) $(TIMING_TARGET) $(TEST_TARGET) $(PRINTNASA_TARGET)

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

$(PLOT_OBJ): $(PLOT_SRC) | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(CSV_OBJ): $(CSV_SRC) | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(TIMING_OBJ): $(TIMING_SRC) | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(TEST_OBJ): $(TEST_SRC) | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(PRINTNASA_OBJ): $(PRINTNASA_SRC) | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Link final executable
$(MIN_TARGET): $(MIXES_OBJ) $(SOLVER_OBJ) $(MIN_OBJ)
	$(CXX) $(CXXFLAGS) $^ -o $@

$(PLOT_TARGET): $(MIXES_OBJ) $(SOLVER_OBJ) $(PLOT_OBJ)
	$(CXX) $(CXXFLAGS) $^ -o $@

$(CSV_TARGET): $(MIXES_OBJ) $(SOLVER_OBJ) $(CSV_OBJ)
	$(CXX) $(CXXFLAGS) $^ -o $@

$(TIMING_TARGET): $(MIXES_OBJ) $(SOLVER_OBJ) $(TIMING_OBJ)
	$(CXX) $(CXXFLAGS) $^ -o $@

$(TEST_TARGET): $(MIXES_OBJ) $(SOLVER_OBJ) $(TEST_OBJ)
	$(CXX) $(CXXFLAGS) $^ -o $@

$(PRINTNASA_TARGET): $(MIXES_OBJ) $(SOLVER_OBJ) $(PRINTNASA_OBJ)
	$(CXX) $(CXXFLAGS) $^ -o $@	

# Clean up
clean:
	rm -rf $(BUILD_DIR) 

.PHONY: all clean
