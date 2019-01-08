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

GSLLIBROOT="../GSL-lib"
WFLAGS = -Werror -Wall -Wextra -pedantic -Wshadow -Wnon-virtual-dtor -Wold-style-cast -Wcast-align -Wunused -Woverloaded-virtual -Wpedantic -Wconversion -Wsign-conversion -Wmisleading-indentation -Wduplicated-cond -Wduplicated-branches -Wlogical-op -Wnull-dereference -Wuseless-cast -Wdouble-promotion -Wformat=2 -Weffc++
# Flags for the above defined compilers
CXXFLAGS = -std=c++11 $(WFLAGS) -I $(SRC_DIR) -I $(GSLLIBROOT)/include -O0 -g -DDEBUG

CXXCHECKS =clang-analyzer-*,-clang-analyzer-cplusplus*,cppcoreguidelines-*,bugprone-* 
CXXCHECKFLAGS = -checks=$(CXXCHECKS) -header-filter=.* -- -std=c++11

# Libraries to link against
GSLLIBDIR=$(GSLLIBROOT)/lib/GSLpp
LDFLAGS = -pg -L$(GSLLIBDIR) -L. -Wl,-rpath=$(GSLLIBDIR) -lm -lGSLpp -lxc

# List of all executables in this project
EXE = gsl-asw

NUMEROV_OBJ = numerov_solver.o\
	      spherical_fun.o\
	      log_mesh.o\
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
	      main.o

OBJS = $(addprefix $(BUILD_DIR)/, $(NUMEROV_OBJ))

# Targets to always execute, even if new files with the same names exists
.PHONY: build all clean cleanall 

# Build all executables
all: build $(EXE)

# Create object files from c++ sources
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c $? -o $@

# Create object files from c sources
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c
	$(CC) $(CFLAGS) -c $? -o $@

# Link numerov_test
gsl-asw: $(OBJS)
	$(CXX) $(LDFLAGS) $? -o $@

checkall: $(addprefix $(SRC_DIR)/, $(NUMEROV_OBJ:o=cpp))
	$(CXXCHECK) $^ $(CXXCHECKFLAGS)

travis: GSLLIBROOT = GSL-lib-master
travis: CXXFLAGS = -std=c++11 -I$(SRC_DIR) -I $(GSLLIBROOT)/include -O0
travis: all

build : 
	mkdir -p build

# Remove object files
clean:
	rm -f $(BUILD_DIR)/*.o

# Remove executables and object files
cleanall:
	rm -f $(EXE) $(BUILD_DIR)/*.o
