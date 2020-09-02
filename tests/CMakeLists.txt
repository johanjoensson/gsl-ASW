set(SOURCE_FILES	spherical_functions_test.cpp
	   		bloch_sum_test.cpp
	   		bloch_summed_structure_constants.cpp
			main_test.cpp
)

add_executable(test ${SOURCE_FILES})
target_link_libraries(test PRIVATE gsl-asw-lib)
