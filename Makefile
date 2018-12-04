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
CXX = g++
#CC  = clang
CC  = gcc

CXXCHECK = clang-tidy

# Source files directory
SRC_DIR = src

# Build directory
BUILD_DIR = build

# Flags for the above defined compilers
CXXFLAGS = -g -std=c++11 -Wall -Wextra -Werror -W -pedantic -I $(SRC_DIR) -DDEBUG
CFLAGS = -g -std=c11 -Wall -Wextra -Werror -W -pedantic -I $(SRC_DIR) -DDEBUG

CXXCHECKS =clang-analyzer-*,-clang-analyzer-cplusplus*,cppcoreguidelines-*,bugprone-* 
CXXCHECKFLAGS = -checks=$(CXXCHECKS) -header-filter=.* -- -std=c++11

# Libraries to link against
GSLLIBDIR="../GSL-lib"
LDFLAGS = -L$(GSLLIBDIR) -Wl,-rpath=$(GSLLIBDIR) -lm -lgsl-lib -lxc

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
	      main.o

OBJS = $(addprefix $(BUILD_DIR)/, $(NUMEROV_OBJ))

# Targets to always execute, even if new files with the same names exists
.PHONY: all clean cleanall

# Build all executables
all: $(EXE)

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

# Remove object files
clean:
	rm -f $(BUILD_DIR)/*.o

# Remove executables and object files
cleanall:
	rm -f $(EXE) $(BUILD_DIR)/*.o
