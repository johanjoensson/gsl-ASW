cmake_minimum_required(VERSION 3.0)

project(gsl-asw)

set(WFLAGS " -Werror -Wall -Wextra -pedantic -Wshadow -Wnon-virtual-dtor -Wold-style-cast -Wcast-align -Wunused -Woverloaded-virtual -Wpedantic -Wconversion -Wsign-conversion -Wmisleading-indentation -Wduplicated-cond -Wduplicated-branches -Wlogical-op -Wnull-dereference -Wuseless-cast -Wdouble-promotion -Wformat=2 -Weffc++")

set(CXXFLAGS  "-g -std=c++11 ${WFLAGS} -DDEBUG")
set(CFLAGS  "-g -std=c11 ${WFLAGS} -DDEBUG")

set(SOURCE_FILES    src/numerov_solver.cpp
		    src/spherical_fun.cpp
		    src/log_mesh.cpp
		    src/bloch_sum.cpp
		    src/gaunt.cpp
		    src/structure_const.cpp
		    src/augmented_fun.cpp
		    src/augmented_spherical_wave.cpp
		    src/atomic_quantity.cpp
		    src/atom.cpp
		    src/crystal.cpp
		    src/lattice.cpp
		    src/ewald_int.cpp
		    src/utils.cpp
		    src/xc_func.cpp
		    src/envelope_fun.cpp
		    src/simulation.cpp
		    src/k-mesh.cpp
		    src/main.cpp)		

	    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CXXFLAGS}")
	    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${CFLAGS}")

add_executable(gsl-asw ${SOURCE_FILES})
find_library(GSL_LIB GSLpp ../GSL-lib)
target_link_libraries(gsl-asw m ${GSL_LIB} xc)
