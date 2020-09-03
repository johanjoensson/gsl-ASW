set(TEST_FILES	spherical_functions_test.cpp
	   		bloch_sum_test.cpp
	   		bloch_summed_structure_constants.cpp
			main_test.cpp
)

enable_testing()
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})
add_executable(gsl-asw-test ${TEST_FILES})
target_link_libraries(gsl-asw-test gsl-asw-lib ${GTEST_LIBRARIES} pthread)
GTEST_ADD_TESTS(gsl-asw-test "" ${TEST_FILES})
message(STATUS "CXX flags : ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG}")
