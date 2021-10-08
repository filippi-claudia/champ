cmake_minimum_required(VERSION 3.1 FATAL_ERROR)

# Project's name
project(CHAMP LANGUAGES Fortran CXX C)

# Needed for something which we don't know
cmake_policy(SET CMP0079 NEW)

set(CMAKE_VERBOSE_MAKEFILE OFF)

# ## Build type
# if(NOT CMAKE_BUILD_TYPE)
#   set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
# endif()

# Cmake modules/macros are in a subdirectory to keep this file cleaner
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules)

# useful info
message(STATUS "Compiling                        :: " ${CMAKE_PROJECT_NAME})
message(STATUS "Fortran compiler                 :: " ${CMAKE_Fortran_COMPILER} )
message(STATUS "Fortran compiler version         :: " ${CMAKE_Fortran_COMPILER_VERSION} )


# Set folder variables:
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(CMAKE_VMC_DIR ${CMAKE_SOURCE_DIR}/src/vmc)
set(CMAKE_DMC_DIR ${CMAKE_SOURCE_DIR}/src/dmc)
set(CMAKE_PARSER_DIR ${CMAKE_SOURCE_DIR}/src/parser)

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})

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
  message("")
  message(STATUS "Using MKL for the linear algebra calculations!")
  message(STATUS "MKL INCLUDE DIRS                 :: " ${MKL_INCLUDE_DIRS})
# message(STATUS "MKL_LIBRARIES            :: " ${MKL_LIBRARIES})
  # Printing individual libraries in the list
  message(STATUS "MKL_LIBRARIES                    : ")
  foreach(lib ${MKL_LIBRARIES})
    message("                                    :: " ${lib})
  endforeach(lib)
else()
  # search for lapack
  find_package(LAPACK COMPONENTS Fortran REQUIRED)
  # Search for Blas
  find_package(BLAS COMPONENTS Fortran REQUIRED)
  set(LINEAR_ALGEBRA ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
  message(STATUS "MKL NOT FOUND! Using default BLAS and LAPACK")
endif(MKL_FOUND)
message(STATUS "                                    ")


# LibFDF is built-in in the src/parser folder

# Find HDF5 library
find_package(HDF5)
message(STATUS "Is HDF5 library found            :: " HDF5_FOUND)
message(STATUS "Is HDF5 library parallel         :: " HDF5_IS_PARALLEL)
message(STATUS "HDF5 Library C Compiler exec     :: " ${HDF5_C_COMPILER_EXECUTABLE})
message(STATUS "HDF5 include dirs                :: " ${HDF5_INCLUDE_DIRS})
message(STATUS "HDF5 library directories         :  " )
foreach(lib ${HDF5_LIBRARIES})
  message("                                    :: " ${lib})
endforeach(lib)
message(STATUS "                                    ")
add_definitions(-DHDF5_VERSION=\"${HDF5_VERSION}\")

# Find TREXIO
include(${CMAKE_MODULE_PATH}/FindTREXIO.cmake)
find_package(TREXIO)
message(STATUS "Is TREXIO library found          :: " TREXIO_FOUND)
message(STATUS "TREXIO Library include dirs      :: " ${TREXIO_INCLUDE_DIR})
message(STATUS "TREXIO Library lib dirs          :: " ${TREXIO_LIBRARY})
message(STATUS "                                    ")

#  documentation requirements
if (BUILD_DOC)
  find_package(DOXYGEN)
endif(BUILD_DOC)


# Find MPI
if (ENABLE_MPI)
  find_package(MPI COMPONENTS Fortran REQUIRED)
endif(ENABLE_MPI)

## Know the architecture of the CPU and set compile flags accordingly
set(TARGET_ARCHITECTURE "generic")
if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
   file(READ "/proc/cpuinfo" _cpuinfo)
   string(REGEX REPLACE ".*flags[ \t]*:[ \t]+([^\n]+).*" "\\1" _cpu_flags "${_cpuinfo}")
endif(CMAKE_SYSTEM_NAME STREQUAL "Linux")

if(_cpu_flags MATCHES "avx512")
  list(APPEND Fortran_FLAGS "-xCORE-AVX512")
  set(TARGET_ARCHITECTURE "avx512")
  add_definitions(-DVECTORIZATION="avx512")
  elseif(_cpu_flags MATCHES "avx2")
  list(APPEND Fortran_FLAGS "-xCORE-AVX2")
  set(TARGET_ARCHITECTURE "avx2")
  add_definitions(-DVECTORIZATION="avx2")
endif()

add_definitions(-DTARGET_ARCHITECTURE=\"${TARGET_ARCHITECTURE}\")
## Architecture specific flags end here

# Git info
execute_process(COMMAND git rev-parse --abbrev-ref HEAD  OUTPUT_VARIABLE GIT_BRANCH OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND git rev-parse HEAD OUTPUT_VARIABLE GIT_HASH OUTPUT_STRIP_TRAILING_WHITESPACE)
#

## Compiler FLAGS
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  list(APPEND Fortran_FLAGS "-O2" "-cpp" "-mcmodel=large" "-ffree-line-length-512")
  list(APPEND Fortran_FLAGS "-DCOMPILER=\"${CMAKE_Fortran_COMPILER}\"")
  list(APPEND Fortran_FLAGS "-DCOMPILER_VERSION=\"${CMAKE_Fortran_COMPILER_VERSION}\"")
  list(APPEND Fortran_FLAGS "-DGIT_BRANCH=\"${GIT_BRANCH}\"")
  list(APPEND Fortran_FLAGS "-DGIT_HASH=\"${GIT_HASH}\"")
  list(APPEND Fortran_FLAGS "-D_MPI_")
  list(APPEND Fortran_FLAGS "-DCLUSTER")
  if (TREXIO_FOUND)
    list(APPEND Fortran_FLAGS "-ltrexio")
  endif()

  set(CMAKE_Fortran_FORMAT_FIXED_FLAG "-ffixed-form -ffixed-line-length-512 -fno-automatic")

  if (MKL_FOUND)
    #list(APPEND Fortran_FLAGS "-fdefault-integer-8")
    list(APPEND Fortran_FLAGS "-m64")
  endif()
  #set(CMAKE_Fortran_FORMAT_FREE_FLAG "-ffree-form")

elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  list(APPEND Fortran_FLAGS "-O2")
	list(APPEND Fortran_FLAGS "-fpp" "-mcmodel=small" "-shared-intel" "-dyncom=grid3d_data,orbital_num_spl,orbital_num_lag,orbital_num_spl2,grid3d_data" )
  list(APPEND Fortran_FLAGS "-DCOMPILER=\"${CMAKE_Fortran_COMPILER}\"")
  list(APPEND Fortran_FLAGS "-DCOMPILER_VERSION=\"${CMAKE_Fortran_COMPILER_VERSION}\"")
  list(APPEND Fortran_FLAGS "-DGIT_BRANCH=\"${GIT_BRANCH}\"")
  list(APPEND Fortran_FLAGS "-DGIT_HASH=\"${GIT_HASH}\"")
  list(APPEND Fortran_FLAGS "-D_MPI_")
  list(APPEND Fortran_FLAGS "-DCLUSTER")
  if (TREXIO_FOUND)
    list(APPEND Fortran_FLAGS "-ltrexio")
  endif()

  set(CMAKE_Fortran_FORMAT_FIXED_FLAG "-fixed -132")
  #if (MKL_FOUND)
  #    list(APPEND Fortran_FLAGS "-i8")
  #endif()

  # debug
  #list(APPEND Fortran_FLAGS "-g" "-O0" "-debug" "-traceback")


elseif(Fortran_COMPILER_ID MATCHES "PGI")
  list(APPEND CMAKE_Fortran_FLAGS "-Mfreeform -Mdclchk -Mstandard -Mallocatable=03")
endif()

message("")
message(STATUS "Fortran compiler flags   : " )
foreach(flag ${Fortran_FLAGS})
  message("                                    :: " ${flag})
endforeach(flag)


## Subdirectories
# Include the following subdirectories
add_subdirectory(src/parser)
add_subdirectory(tools/interface)
add_subdirectory(lib/pspline)
add_subdirectory(src/vmc)
add_subdirectory(src/dmc)

target_link_libraries(vmc.mov1 ${TREXIO_LIBRARY})
target_link_libraries(dmc.mov1 ${TREXIO_LIBRARY})

if(ENABLE_TEST)
  MESSAGE(STATUS "enable test")
  enable_testing()
  add_subdirectory(tests/CI_test)
#  add_test(example ${FORTRAN_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test.py)
endif(ENABLE_TEST)
