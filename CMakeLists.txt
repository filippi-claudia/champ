cmake_minimum_required(VERSION 3.1 FATAL_ERROR)

# Project's name
project(champ LANGUAGES Fortran CXX)

set(CMAKE_VERBOSE_MAKEFILE off)

# ## Build type
# if(NOT CMAKE_BUILD_TYPE)
#   set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
# endif()

# Cmake modules/macros are in a subdirectory to keep this file cleaner
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules)


# Set the output folder where your program will be created
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})
# Input directory
set(INPUT_DIR ${CMAKE_SOURCE_DIR}/src/input)

## User options
# documentation
option(BUILD_DOC "build_doc" off)

# MPI 
option(ENABLE_MPI "enable_mpi" on)

# QMMM option
option(ENABLE_QMMM "enable_qmmm" off)

# Periodic option
option(ENABLE_PERIODIC "enable_periodic" off)

## External packages
# Include function to Load MKL
include(${CMAKE_MODULE_PATH}/FindMKL.cmake)
# Find MKL
find_package(MKL)
if (MKL_FOUND)
  include_directories(${MKL_INCLUDE_DIRS})
  set(MKL_INCLUDE ${MKL_INCLUDE_DIRS}/intel64/ilp64)
  set(LINEAR_ALGEBRA ${MKL_LIBRARIES} pthread m dl)
  message(STATUS "Using MKL for the linear algebra calculations!")
  message(STATUS "MKL_INCLUDE_DIRS: " ${MKL_INCLUDE_DIRS})
  message(STATUS "MKL_LIBRARIES: " ${MKL_LIBRARIES})
else()
  # search for lapack
  find_package(LAPACK COMPONENTS Fortran REQUIRED)
  # Search for Blas
  find_package(BLAS COMPONENTS Fortran REQUIRED)
  set(LINEAR_ALGEBRA ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
  message(STATUS "MKL NOT FOUND! Using default BLAS and LAPACK")
endif(MKL_FOUND)

  
#  documentation requirements
if (BUILD_DOC)
  find_package(DOXYGEN)
endif(BUILD_DOC)


# Find MPI
if (ENABLE_MPI)
  find_package(MPI COMPONENTS Fortran REQUIRED)
endif(ENABLE_MPI)


## Compiler FLAGS
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  list(APPEND Fortran_FLAGS "-O2" "-mcmodel=large")
  set(CMAKE_Fortran_FORMAT_FIXED_FLAG "-ffixed-form -ffixed-line-length-132 -fno-automatic")
  if (MKL_FOUND)
    #list(APPEND Fortran_FLAGS "-fdefault-integer-8")
    list(APPEND Fortran_FLAGS "-m64")
  endif()  
  #set(CMAKE_Fortran_FORMAT_FREE_FLAG "-ffree-form")
  
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  list(APPEND Fortran_FLAGS "-O2" "-mcmodel=large" "-shared-intel"
    "-dyncom=grid3d_data,orbital_num_spl,orbital_num_lag,orbital_num_spl2,grid3d_data")
  set(CMAKE_Fortran_FORMAT_FIXED_FLAG "-fixed -132")
  #if (MKL_FOUND)
  #    list(APPEND Fortran_FLAGS "-i8")
  #endif()  
  
elseif(Fortran_COMPILER_ID MATCHES "PGI")
  list(APPEND CMAKE_Fortran_FLAGS "-Mfreeform -Mdclchk -Mstandard -Mallocatable=03")
endif()
message(STATUS "compiler flags are: " ${Fortran_FLAGS})


## Subdirectories
# Include the following subdirectories
add_subdirectory(src/include)
add_subdirectory(src/input)
add_subdirectory(tools/interface)
add_subdirectory(lib/cyrus)
add_subdirectory(lib/pspline)
add_subdirectory(src/vmc)
add_subdirectory(src/dmc)

if(ENABLE_TEST)
  MESSAGE(STATUS "enable test")
  enable_testing()
  add_subdirectory(tests)
endif(ENABLE_TEST)
