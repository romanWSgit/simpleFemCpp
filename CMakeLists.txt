cmake_minimum_required(VERSION 3.24)
project(simpleFemCpp)

message("cmake dir :${CMAKE_CURRENT_SOURCE_DIR}")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# first we can indicate the documentation build as an option and set it to ON by default
option(BUILD_DOC "Build documentation" ON)

find_package(Eigen3 REQUIRED)
if(NOT EXISTS ${EIGEN3_INCLUDE_DIR})
    message(FATAL_ERROR "Please check that the set Eigen directory is valid!")
endif()
include_directories(${EIGEN3_INCLUDE_DIR})



include(FetchContent)
FetchContent_Declare(
        fast_matrix_market
        GIT_REPOSITORY https://github.com/alugowski/fast_matrix_market
        GIT_TAG main
        GIT_SHALLOW TRUE
)
FetchContent_MakeAvailable(fast_matrix_market)



#FIND_PACKAGE(PLplot REQUIRED)
#SET(PLplot_LIBRARIES ${PLplot_LIBRARIES} ${PLplot_cxx_LIBRARY})
#if(NOT EXISTS ${PLPLOT_INCLUDE_DIR})
#    message(FATAL_ERROR "Please check that the set PLplot directory is valid!")
#endif()
#include_directories(${PLPLOT_INCLUDE_DIR})




set(VTK_DIR "/opt/homebrew/Cellar/vtk/9.2.6")
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(simpleFemCpp src/fem.cpp src/fem.h)
target_link_libraries(simpleFemCpp
        PRIVATE
        ${VTK_LIBRARIES}
#        ${PLPLOT_LIBRARIES}
        fast_matrix_market::fast_matrix_market
        )


message("cmake dir :${CMAKE_CURRENT_SOURCE_DIR}")
# check if Doxygen is installed
find_package(Doxygen)
if (DOXYGEN_FOUND)
    # set input and output files
    set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/docs/Doxyfile.in)
    set(DOXYGEN_OUT ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile)

    # request to configure the file
    configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
    message("Doxygen build started")
    message("Doxygen executable = ${DOXYGEN_EXECUTABLE}")
    # note the option ALL which allows to build the docs together with the application
    add_custom_target( doc_doxygen ALL
            COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
            COMMENT "Generating API documentation with Doxygen"
            VERBATIM )
else (DOXYGEN_FOUND)
    message("Doxygen need to be installed to generate the doxygen documentation")
endif (DOXYGEN_FOUND)
