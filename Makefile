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
#CXX = clang++
CXX ?= g++
#CC  = clang
CC  ?= gcc

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
CXXFLAGS = -std=c++14 $(WFLAGS) -I $(SRC_DIR) -I $(GSLLIBROOT)/include -O0  -DDEBUG -g

CXXCHECKS =clang-analyzer-*,-clang-analyzer-cplusplus*,cppcoreguidelines-*,bugprone-* 
CXXCHECKFLAGS = -checks=$(CXXCHECKS) -header-filter=.* -- -std=c++11 -I$(GSLLIBROOT)/include

# Libraries to link against
GSLLIBDIR=$(GSLLIBROOT)/lib/GSLpp
LDFLAGS = -pg -L$(GSLLIBDIR) -L. -Wl,-rpath=$(GSLLIBDIR) -lGSLpp -lxc -lm -lgsl -lpthread # -O3 -flto

# List of all executables in this project
EXE = gsl-asw

ASW_OBJ =     bloch_sum.o\
	      gaunt.o\
	      structure_const.o\
	      augmented_fun.o\
	      augmented_spherical_wave.o\
	      atomic_quantity.o\
	      atom.o\
	      lattice.o\
	      ewald_int.o\
	      utils.o\
	      xc_func.o\
	      envelope_fun.o\
	      simulation.o\
	      k-mesh.o\
	      log_mesh.o\
	      spherical_fun.o\

ASW_HEADERS = numerov_solver.h\
	      schroedinger.h\


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

TB_HEADERS =  numerov_solver.h\
	      schroedinger.h\

TEST_OBJ = $(ASW_OBJ)\
	   bloch_summed_structure_constants.o\
	   bloch_sum_test.o\
	   spherical_functions_test.o


OBJS = $(addprefix $(BUILD_DIR)/, $(ASW_OBJ)) $(BUILD_DIR)/main.o
TEST_OBJS = $(addprefix $(BUILD_DIR)/, $(TEST_OBJ)) $(BUILD_DIR)/main_test.o

# Targets to always execute, even if new files with the same names exists
.PHONY: build all clean cleanall 

# Build all executables
all: build $(EXE)

# Create object files from c++ sources
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp $(SRC_DIR)/%.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BUILD_DIR)/%.o: $(TEST_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c $? -o $@

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
tests:  LDFLAGS = -lgcov -lgtest -L$(GSLLIBDIR) -L. -Wl,-rpath=$(GSLLIBDIR) -lGSLpp -lgsl -lxc -lm 
tests: 	clean $(TEST_OBJS)
	$(CXX) $(TEST_OBJS) -o $@  $(LDFLAGS)

test_mesh: $(BUILD_DIR)/test_mesh.o $(BUILD_DIR)/log_mesh.o
	$(CXX) $^ -o $@ $(LDFLAGS)

test_cubic: $(BUILD_DIR)/test_cubic.o $(BUILD_DIR)/spherical_fun.o $(BUILD_DIR)/ewald_int.o
	$(CXX) $^ -o $@ $(LDFLAGS)

test_numerov: $(BUILD_DIR)/test_new_numerov.o
	$(CXX) $^ -o $@ $(LDFLAGS)

test_schroedinger: $(BUILD_DIR)/test_schroedinger.o $(BUILD_DIR)/log_mesh.o
	$(CXX) $^ -o $@ $(LDFLAGS)

build : 
	mkdir -p build

# Remove object files
clean:
	rm -f $(BUILD_DIR)/*.o $(BUILD_DIR)/*.gcno $(BUILD_DIR)/*.gcda

# Remove executables and object files
cleanall:
	rm -f $(EXE) $(BUILD_DIR)/*.o
