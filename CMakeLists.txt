cmake_minimum_required(VERSION 3.5 FATAL_ERROR)

# Project's name
project(champ LANGUAGES Fortran)

set( CMAKE_VERBOSE_MAKEFILE on )

# ## Build type
# if(NOT CMAKE_BUILD_TYPE)
#   set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
# endif()


## Compiler FLAGS
# Set the output folder where your program will be created
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})

## External packages
# search for lapack
find_package(LAPACK COMPONENTS Fortran REQUIRED)

# Search for Blas
find_package(BLAS COMPONENTS Fortran REQUIRED)

## User options
# documentation
option(BUILD_DOC "build_doc" off)

# MKL library
option(ENABLE_MKL "enable_mkl" off)

#  documentation requirements
if (BUILD_DOC)
  find_package(DOXYGEN)
endif(BUILD_DOC)

if (ENABLE_MKL)
  find_package(MKL)
endif(ENABLE_MKL)

## Subdirectories
# Include the following subdirectories
add_subdirectory(include)
add_subdirectory(input)
add_subdirectory(interface)
add_subdirectory(lib)
add_subdirectory(lib2)
# add_subdirectory(tools)
# add_subdirectory(dmc)
add_subdirectory(vmc)


