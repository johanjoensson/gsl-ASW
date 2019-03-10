################################################################################
# Makefile for building an implementation of the ASW method (someday this might#
# actually be a "full" DFT code)                                               #
################################################################################
################################################################################
# Author: Johan Jönsson                                                        #
# Email: johanjoensson@gmail.com                                               #
# Created: 2018-01-31                                                          #
# Last updated: 2018-02-25                                                     #
################################################################################

# Compilers to use
CXX = clang++
#CXX ?= g++
CC  = clang
#CC  ?= gcc

CXXCHECK = clang-tidy

# Source files directory
SRC_DIR = src

# Build directory
BUILD_DIR = build

# Test directory
TEST_DIR = test
GSLLIBROOT=../GSL-lib
#  -Werror 
WFLAGS = -Wall -Wextra -pedantic -Wshadow -Wnon-virtual-dtor -Wold-style-cast -Wcast-align -Wunused -Woverloaded-virtual -Wpedantic -Wconversion -Wsign-conversion -Wnull-dereference -Wdouble-promotion -Wformat=2 -Weffc++ -Wmisleading-indentation -Wduplicated-cond -Wduplicated-branches -Wlogical-op  -Wuseless-cast
# Flags for the above defined compilers
CXXFLAGS = -std=c++11 $(WFLAGS) -I $(SRC_DIR) -I $(GSLLIBROOT)/include -O0  -DDEBUG -g

CXXCHECKS =clang-analyzer-*,-clang-analyzer-cplusplus*,cppcoreguidelines-*,bugprone-* 
CXXCHECKFLAGS = -checks=$(CXXCHECKS) -header-filter=.* -- -std=c++11

# Libraries to link against
GSLLIBDIR=$(GSLLIBROOT)/lib/GSLpp
LDFLAGS = -pg -L$(GSLLIBDIR) -L. -Wl,-rpath=$(GSLLIBDIR) -lGSLpp -lxc -lm -lgsl # -O3 -flto

# List of all executables in this project
EXE = gsl-asw

ASW_OBJ =     main.o\
	      bloch_sum.o\
	      gaunt.o\
	      structure_const.o\
	      augmented_fun.o\
	      augmented_spherical_wave.o\
	      atomic_quantity.o\
	      atom.o\
	      crystal.o\
	      lattice.o\
	      ewald_int.o\
	      utils.o\
	      xc_func.o\
	      envelope_fun.o\
	      simulation.o\
	      k-mesh.o\
	      log_mesh.o\
	      spherical_fun.o\
              numerov_solver.o\

TB_OBJS =     tb.o\
	      bloch_sum.o\
	      gaunt.o\
	      structure_const.o\
	      augmented_fun.o\
	      augmented_spherical_wave.o\
	      atomic_quantity.o\
	      atom.o\
	      crystal.o\
	      lattice.o\
	      ewald_int.o\
	      utils.o\
	      xc_func.o\
	      envelope_fun.o\
	      simulation.o\
	      k-mesh.o\
	      log_mesh.o\
	      spherical_fun.o\
              numerov_solver.o\

TEST_OBJ =	


OBJS = $(addprefix $(BUILD_DIR)/, $(ASW_OBJ))
TEST_OBJS = $(addprefix $(TEST_DIR)/, $(TEST_OBJ))

# Targets to always execute, even if new files with the same names exists
.PHONY: build all clean cleanall 

# Build all executables
all: build $(EXE)

# Create object files from c++ sources
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp $(SRC_DIR)/%.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Create object files from c sources
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c $(SRC_DIR)/%.h
	$(CC) $(CFLAGS) -c $? -o $@

# Link numerov_test
gsl-asw: $(OBJS)
	$(CXX) $^ -o $@ $(LDFLAGS)

TB: $(addprefix $(BUILD_DIR)/, $(TB_OBJS))
	$(CXX) $^ -o $@ $(LDFLAGS)

checkall: $(addprefix $(SRC_DIR)/, $(ASW_OBJ:o=cpp))
	$(CXXCHECK) $^ $(CXXCHECKFLAGS)

travis: GSLLIBROOT = GSL-lib-master
travis: GSLLIBDIR = $(GSLLIBROOT)/lib/GSLpp
travis: CXXFLAGS = -std=c++11 -I $(SRC_DIR) -I $(GSLLIBROOT)/include -O0
travis: all

tests: 	CXXFLAGS = -std=c++11 -I $(SRC_DIR) -I $(GSLLIBROOT)/include -I$(TEST_DIR) -O0 -fprofile-arcs -ftest-coverage
tests:  LDFLAGS = -lgcov -lgtest -L$(GSLLIBDIR) -L. -Wl,-rpath=$(GSLLIBDIR) -lGSLpp -lxc -lm 
tests: 	clean $(TEST_OBJS)
	$(CXX) $(TEST_OBJS) -o $@  $(LDFLAGS)

test_mesh: $(BUILD_DIR)/test_mesh.o $(BUILD_DIR)/log_mesh.o
	$(CXX) $^ -o $@ $(LDFLAGS)

test_cubic: $(BUILD_DIR)/test_cubic.o $(BUILD_DIR)/spherical_fun.o $(BUILD_DIR)/ewald_int.o
	$(CXX) $^ -o $@ $(LDFLAGS)

test_numerov: $(BUILD_DIR)/test_new_numerov.o
	$(CXX) $^ -o $@ $(LDFLAGS)

build : 
	mkdir -p build

# Remove object files
clean:
	rm -f $(BUILD_DIR)/*.o

# Remove executables and object files
cleanall:
	rm -f $(EXE) $(BUILD_DIR)/*.o
