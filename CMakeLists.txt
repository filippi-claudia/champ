cmake_minimum_required(VERSION 3.5 FATAL_ERROR)

# Project's name
project(champ LANGUAGES Fortran)

set( CMAKE_VERBOSE_MAKEFILE on )

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
set(INPUT_DIR ${CMAKE_SOURCE_DIR}/input)

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
  set(LINEAR_ALGEBRA ${MKL_LIBRARIES})
  message(STATUS "Using MKL for the linear algebra calculations!")
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
  find_package(MPI)
endif(ENABLE_MPI)


## Compiler FLAGS
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  list(APPEND Fortran_FLAGS "-O2" "-ffixed-line-length-132" "-fno-automatic")
  
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  list(APPEND Fortran_FLAGS "-O2" "-xHost" "-mcmodel=large" "-shared-intel" "-132"
    "-dyncom=grid3d_data,orbital_num_spl,orbital_num_lag,orbital_num_spl2,grid3d_data")
  list(APPEND Fortan_FLAGS_DEBUG "-g")
  
elseif(Fortran_COMPILER_ID MATCHES "PGI")
  list(APPEND CMAKE_Fortran_FLAGS "-Mfreeform -Mdclchk -Mstandard -Mallocatable=03")
endif()
message(STATUS "compiler flags are: " ${Fortran_FLAGS})

# # are crucial for MKL to avoid SIGSEGV at runtime!
# if (MKL_FOUND AND (CMAKE_Fortran_COMPILER_ID MATCHES "GNU"))
#  list(APPEND LINK_FLAGS "-m64" "-Wl" "--no-as-needed")
# elseif(MKL_FOUND AND (CMAKE_Fortran_COMPILER_ID MATCHES "Intel"))
#  list(APPEND LINK_FLAGS "-liomp5" "-lpthread" "-lm" "-ldl")
# endif()




## Subdirectories
# Include the following subdirectories
add_subdirectory(include)
add_subdirectory(input)
add_subdirectory(interface)
add_subdirectory(lib)
add_subdirectory(lib2)
# add_subdirectory(tools)
add_subdirectory(vmc)
add_subdirectory(dmc)


