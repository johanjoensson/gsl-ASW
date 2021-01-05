# find all header files
get_target_property(GSL_ASW_PUBLIC_HEADER_DIR gsl-asw-lib INTERFACE_INCLUDE_DIRECTORIES)
file(GLOB_RECURSE GSL_ASW_PUBLIC_HEADERS ${GSL_ASW_PUBLIC_HEADER_DIR}/*.h)
set(DOXYGEN_INPUT_DIR ${PROJECT_SOURCE_DIR}/include/gsl-asw)
set(DOXYGEN_OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/doxygen)
set(DOXYGEN_INDEX_FILE ${CMAKE_DOXYGEN_OUTPUT_DIR}/html/index.html)
set(DOXYFILE_IN ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in)
set(DOXYFILE_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

configure_file(${DOXYFILE_IN} ${DOXYFILE_OUT} @ONLY)
file(MAKE_DIRECTORY ${DOXYGEN_OUTPUT_DIR}) #Doxygen won't create this for us

add_custom_command(OUTPUT ${DOXYGEN_INDEX_FILE}
               DEPENDS ${NUMERICAL_MESH_PUBLIC_HEADERS}
               COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYFILE_OUT}
               WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
               MAIN_DEPENDENCY ${DOXYFILE_OUT} ${DOXYFILE_IN}
               COMMENT "Generating docs")


find_package(Sphinx)
if(Sphinx_FOUND)
    set(SPHINX_SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/source)
    set(SPHINX_BUILD ${CMAKE_CURRENT_BINARY_DIR}/sphinx)
    set(SPHINX_INDEX_FILE ${SPHINX_BUILD}/index.html)

    add_custom_command(OUTPUT ${SPHINX_INDEX_FILE}
                   COMMAND
                     ${SPHINX_EXECUTABLE} -b html
                     -Dbreathe_projects.GSL-ASW=${DOXYGEN_OUTPUT_DIR}/xml
                   ${SPHINX_SOURCE} ${SPHINX_BUILD}
                   WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                   DEPENDS
                   # Other docs files you want to track should go here (or in some variable)
                   ${CMAKE_CURRENT_SOURCE_DIR}/source/index.rst
                   ${DOXYGEN_INDEX_FILE}
                   MAIN_DEPENDENCY ${SPHINX_SOURCE}/conf.py
                   COMMENT "Generating documentation with Sphinx")

    add_custom_target(docs ALL DEPENDS ${SPHINX_INDEX_FILE})
else()
   add_custom_target(docs ALL DEPENDS ${DOXYGEN_INDEX_FILE})
endif()
