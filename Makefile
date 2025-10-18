# ===== Makefile =====
# Compiler and flags
CXX := g++
CXXFLAGS := -std=c++17 -O2 -I./includes -I./solverLib -I./common_mixes

# Directories
COMMON_AIR_DIR := common_mixes
SOLVER_DIR     := solverLib
BUILD_DIR      := build

# Files
COMMON_AIR_SRC := $(COMMON_AIR_DIR)/commonMixes.cpp
SOLVER_SRC     := $(SOLVER_DIR)/equilibrium.cpp
TEST_SRC       := main.cpp

COMMON_AIR_OBJ := $(BUILD_DIR)/commonMixes.o
SOLVER_OBJ     := $(BUILD_DIR)/equilibrium.o
TEST_OBJ       := $(BUILD_DIR)/main.o

TEST_TARGET    := gibbs

# ===== Targets =====
all: $(TEST_TARGET)

# Create build directory if missing
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

# Compile source files
$(COMMON_AIR_OBJ): $(COMMON_AIR_SRC) | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(SOLVER_OBJ): $(SOLVER_SRC) | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(TEST_OBJ): $(TEST_SRC) | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Link final executable
$(TEST_TARGET): $(COMMON_AIR_OBJ) $(SOLVER_OBJ) $(TEST_OBJ)
	$(CXX) $(CXXFLAGS) $^ -o $@

# Clean up
clean:
	rm -rf $(BUILD_DIR) $(TEST_TARGET)

.PHONY: all clean
