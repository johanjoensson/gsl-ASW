################################################################################
# Makefile for building an implementation of the ASW method (someday this might#
# actually be a "full" DFT code)                                               #
################################################################################
################################################################################
# Author: Johan Jönsson                                                        #
# Email: johanjoensson@gmail.com                                               #
# Created: 2018-01-31                                                          #
# Last updated: 2018-01-31                                                     #
################################################################################

# Compilers to use
CXX = g++
CC  = gcc

# Flags for the above defined compilers
CXXFLAGS = -g -std=c++11 -Wall -Wextra -Werror -W -pedantic
CFLAGS = -g -std=c11 -Wall -Wextra -Werror -W -pedantic

# Libraries to link against
LDFLAGS = -lgsl -lm

# List of all executables in this project
EXE = numerov_test

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
numerov_test: numerov_solver.o log_mesh.o
	$(CXX) $(LDFLAGS) $^ -o $@

# Remove object files
clean:
	rm -f *.o

# Remove executables and object files
cleanall:
	rm -f $(EXE) *.o
