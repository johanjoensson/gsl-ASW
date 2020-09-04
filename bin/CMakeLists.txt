add_executable(gsl-asw main.cpp) 

set_property(TARGET gsl-asw PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)

target_link_libraries(gsl-asw PRIVATE gsl-asw-lib m xc ${GSL_LIB}) 
target_compile_features(gsl-asw PRIVATE cxx_std_11) 

target_compile_options(gsl-asw PRIVATE $<$<CONFIG:DEBUG>:${WFLAGS} -g -O1 -DDEBUG>)
target_link_options(gsl-asw PRIVATE $<$<CONFIG:DEBUG>:-O1>)
target_compile_options(gsl-asw PRIVATE $<$<CONFIG:RELEASE>:${WFLAGS} -O2>)
target_link_options(gsl-asw PRIVATE $<$<CONFIG:RELEASE>:-O2>)
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
	target_compile_options(gsl-asw PRIVATE $<$<CONFIG:DEBUG>:-Og -pg>)
	target_link_options(gsl-asw PRIVATE $<$<CONFIG:DEBUG>:-Og -pg>)
	target_compile_options(gsl-asw PRIVATE $<$<CONFIG:RELEASE>:-Ofast>)
	target_link_options(gsl-asw PRIVATE $<$<CONFIG:RELEASE>:-Ofast>)
endif()
