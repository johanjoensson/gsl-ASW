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
CXXCHECKFLAGS = -checks=$(CXXCHECKS) -- -std=c++11

# Libraries to link against
LDFLAGS = -lgsl -lgslcblas -lm

# List of all executables in this project
EXE = numerov_test

NUMEROV_OBJ = numerov_solver.o \
	      spherical_fun.o\
	      log_mesh.o \
	      gaunt.o \
	      structure_const.o\
	      main.o

# Targets to always execute, even if new files with the same names exists
.PHONY: all clean cleanall

# Build all executables
all: $(EXE)

# Create object files from c++ sources
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $<

# Create object files from c sources
%.o: %.c
	$(CC) $(CFLAGS) -c $<

# Link numerov_test
numerov_test: $(NUMEROV_OBJ)
	$(CXX) $(LDFLAGS) $^ -o $@

checkall: $(NUMEROV_OBJ:o=cpp)
	$(CXXCHECK) $^ $(CXXCHECKFLAGS)

# Remove object files
clean:
	rm -f *.o

# Remove executables and object files
cleanall:
	rm -f $(EXE) *.o
