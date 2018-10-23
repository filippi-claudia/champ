cmake_minimum_required(VERSION 3.0.2)

enable_language(Fortran)

# Project's name
project(champ)

set( CMAKE_VERBOSE_MAKEFILE on )

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
add_subdirectory(lib)
add_subdirectory(interface)
# add_subdirectory(tools)
# add_subdirectory(dmc)
# add_subdirectory(vmc)

## Compiler FLAGS
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set(dialect "-O2 -ffixed-form -ffixed-line-length-132 -fno-automatic")

elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  set(dialect "-O2 -xHost -mcmodel=large -shared-intel -132 -save -g")
  set(memory "-dyncom=grid3d_data,orbital_num_spl,orbital_num_lag,orbital_num_spl2,grid3d_data")
  set(check " -check bounds -check format ")

elseif(CMAKE_Fortran_COMPILER_ID MATCHES "PGI")
  set(dialect "-Mfreeform -Mdclchk -Mstandard -Mallocatable=03")

endif()

## set compiler flags
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${dialect} ${memory} ${check}")
# set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${check}")



