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
		    # simulation.cpp
		    k-mesh.cpp
)
add_library(gsl-asw-lib ${SOURCE_FILES})

check_ipo_supported(RESULT ipo_result)
if(ipo_result)
	set_property(TARGET gsl-asw-lib PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()

target_include_directories(gsl-asw-lib PUBLIC ${PROJECT_SOURCE_DIR}/include/gsl-asw)
#target_include_directories(gsl-asw-lib PUBLIC ${GSLpp_ROOT}/include/)

target_link_libraries(gsl-asw-lib PUBLIC m xc ${GSL_LIB} Numerical-mesh::numerical-mesh-lib GSLpp::GSLpp)

# All users of this library will need at least C++11
target_compile_features(gsl-asw-lib PUBLIC cxx_std_11)

add_executable(test_brillouin_zone_integration test_brillouin_zone_integration.cpp)
target_include_directories(test_brillouin_zone_integration PUBLIC ../include/)
#target_include_directories(test_brillouin_zone_integration PUBLIC ${GSLpp_ROOT}/include/)
target_link_libraries(test_brillouin_zone_integration PRIVATE m ${GSL_LIB} Numerical-mesh::numerical-mesh-lib GSLpp::GSLpp)



if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
	target_compile_options(gsl-asw-lib BEFORE PUBLIC $<$<CONFIG:DEBUG>:${WFLAGS} -Og -pg -DDEBUG -march=native>)
	target_link_options(gsl-asw-lib BEFORE PUBLIC $<$<CONFIG:DEBUG>:-Og -pg -march=native>)

	target_compile_options(gsl-asw-lib BEFORE PUBLIC $<$<CONFIG:RELEASE>:${WFLAGS} -Ofast -march=native>)
	target_link_options(gsl-asw-lib BEFORE PUBLIC $<$<CONFIG:RELEASE>: -Ofast -march=native>)
else()
	target_compile_options(gsl-asw-lib BEFORE PUBLIC $<$<CONFIG:DEBUG>:${WFLAGS} -Og -pg -DDEBUG>)
	target_link_options(gsl-asw-lib BEFORE PUBLIC $<$<CONFIG:DEBUG>:-Og -pg>)

	target_compile_options(gsl-asw-lib BEFORE PUBLIC $<$<CONFIG:RELEASE>:${WFLAGS} -O2>)
	target_link_options(gsl-asw-lib BEFORE PUBLIC $<$<CONFIG:RELEASE>:-O2>)
endif()

# IDEs should put the headers in a nice place
source_group(
	TREE "${PROJECT_SOURCE_DIR}/include"
	PREFIX "Header Files"
	FILES ${HEADER_LIST}
)
