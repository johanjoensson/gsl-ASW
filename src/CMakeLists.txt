set(SOURCE_FILES    spherical_fun.cpp
		    log_mesh.cpp
		    bloch_sum.cpp
		    gaunt.cpp
		    structure_const.cpp
		    augmented_fun.cpp
		    augmented_spherical_wave.cpp
		    atomic_quantity.cpp
		    atom.cpp
		    lattice.cpp
		    ewald_int.cpp
		    utils.cpp
		    xc_func.cpp
		    envelope_fun.cpp
		    simulation.cpp
		    k-mesh.cpp
)		
add_library(gsl-asw-lib ${SOURCE_FILES})

target_include_directories(gsl-asw-lib PUBLIC ../include/)
target_include_directories(gsl-asw-lib PUBLIC ${GSLpp_ROOT}/include/)
find_library(GSL_LIB GSLpp ${GSLpp_ROOT}/lib/GSLpp/ REQUIRED)

# This depends on (header only) boost 
target_link_libraries(gsl-asw-lib PUBLIC m ${GSL_LIB} xc ) 
# All users of this library will need at least C++11 
target_compile_features(gsl-asw-lib PUBLIC cxx_std_11) 

target_compile_options(gsl-asw-lib PRIVATE $<$<CONFIG:DEBUG>:${WFLAGS} -g -O1 -DDEBUG>)
target_compile_options(gsl-asw-lib PRIVATE $<$<CONFIG:RELEASE>:${WFLAGS} -O2>)
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
	target_compile_options(gsl-asw-lib PRIVATE $<$<CONFIG:DEBUG>:-Og -pg>)
	target_compile_options(gsl-asw-lib PRIVATE $<$<CONFIG:RELEASE>:-Ofast>)
endif()

# IDEs should put the headers in a nice place 
source_group(   
	TREE "${PROJECT_SOURCE_DIR}/include"   
	PREFIX "Header Files"   
	FILES ${HEADER_LIST}
)
