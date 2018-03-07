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
CC  = clang

CXXCHECK = clang-tidy

# Flags for the above defined compilers
CXXFLAGS = -g -std=c++11 -Wall -Wextra -Werror -W -pedantic
CFLAGS = -g -std=c11 -Wall -Wextra -Werror -W -pedantic

CXXCHECKS =clang-analyzer-*,-clang-analyzer-cplusplus*,cppcoreguidelines-*,bugprone-* 
CXXCHECKFLAGS = -checks=$(CXXCHECKS) -header-filter=.* -- -std=c++11

# Libraries to link against
LDFLAGS = -lgsl -lgslcblas -lm

# Source files directory
SRC_DIR = src

# Build directory
BUILD_DIR = build

# List of all executables in this project
EXE = numerov_test

NUMEROV_OBJ = numerov_solver.o\
	      spherical_fun.o\
	      log_mesh.o \
	      gaunt.o \
	      structure_const.o\
	      atom.o\
	      crystal.o\
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
numerov_test: $(OBJS)
	$(CXX) $(LDFLAGS) $? -o $@

checkall: $(NUMEROV_OBJ:o=cpp)
	$(CXXCHECK) $^ $(CXXCHECKFLAGS)

# Remove object files
clean:
	rm -f $(BUILD_DIR)/*.o

# Remove executables and object files
cleanall:
	rm -f $(EXE) $(BUILD_DIR)/*.o
