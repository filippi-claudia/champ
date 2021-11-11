## Ravindra Shinde (2021) (c) TREX-CoE

cmake_minimum_required(VERSION 3.1 FATAL_ERROR)

# Project's name
project(CHAMP VERSION 2.0.5 LANGUAGES Fortran CXX C)


# Needed for something which we don't know
cmake_policy(SET CMP0079 NEW)
cmake_policy(SET CMP0076 NEW)

set(CMAKE_VERBOSE_MAKEFILE OFF)

# ## Build type
# if(NOT CMAKE_BUILD_TYPE)
#   set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
# endif()

# Cmake modules/macros are in a subdirectory to keep this file cleaner
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules)

# Print the header text art
message( "____________________________________________________________________" )
message( "                                                                    " )
message( "                                                                    " )
message( " .d8888b.   888    888         d8888  888b     d888  8888888b.      " )
message( "d88P  Y88b  888    888        d88888  8888b   d8888  888   Y88b     " )
message( "888    888  888    888       d88P888  88888b.d88888  888    888     " )
message( "888         8888888888      d88P 888  888Y88888P888  888   d88P     " )
message( "888         888    888     d88P  888  888 Y888P 888  8888888P''     " )
message( "888    888  888    888    d88P   888  888  Y8P  888  888            " )
message( "Y88b  d88P  888    888   d8888888888  888   '   888  888            " )
message( " 'Y8888P'   888    888  d88P     888  888       888  888            " )
message( "                                                                    " )
message( "____________________________________________________________________" )
message(" ")
message(" ")
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
message("")
include(${CMAKE_MODULE_PATH}/FindMKL.cmake)
# Find MKL
message(CHECK_START "Finding Math Kernel Library (MKL)")
list(APPEND CMAKE_MESSAGE_INDENT "  ")
find_package(MKL)
if (MKL_FOUND)
  include_directories(${MKL_INCLUDE_DIRS})
  set(MKL_INCLUDE ${MKL_INCLUDE_DIRS}/intel64/ilp64)
  set(LINEAR_ALGEBRA ${MKL_LIBRARIES} pthread m dl)
  message(CHECK_PASS "Found")
  message(STATUS "Using MKL for the linear algebra calculations!")
  message(STATUS "MKL INCLUDE DIRS                 :: " ${MKL_INCLUDE_DIRS})
  # Printing individual libraries in the list
  message(STATUS "MKL_LIBRARIES                    : ")
  foreach(lib ${MKL_LIBRARIES})
    message("                                    :: " ${lib})
  endforeach(lib)
  list(POP_BACK CMAKE_MESSAGE_INDENT)
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
if (HDF5_FOUND)
  message(STATUS "Is HDF5 library found            :: " ${HDF5_VERSION})
  message(STATUS "Is HDF5 library parallel         :: " ${HDF5_IS_PARALLEL})
  message(STATUS "HDF5 Library C Compiler exec     :: " ${HDF5_C_COMPILER_EXECUTABLE})
  message(STATUS "HDF5 include dirs                :: " ${HDF5_INCLUDE_DIRS})
  message(STATUS "HDF5 library directories         :  " )
  foreach(lib ${HDF5_LIBRARIES})
    message("                                    :: " ${lib})
  endforeach(lib)
  message(STATUS "                                    ")
  add_definitions(-DHDF5_VERSION=\"${HDF5_VERSION}\")
endif(HDF5_FOUND)

# Find TREXIO
if (HDF5_FOUND)
  include(${CMAKE_MODULE_PATH}/FindTREXIO.cmake)
  find_package(TREXIO)
endif(HDF5_FOUND)
if (TREXIO_FOUND)
  message(STATUS "Is TREXIO library found          :: " ${TREXIO_FOUND})
  message(STATUS "TREXIO Library include dirs      :: " ${TREXIO_INCLUDE_DIR})
  message(STATUS "TREXIO Library lib dirs          :: " ${TREXIO_LIBRARY})
  message(STATUS "                                    ")
  add_definitions(-DTREXIO_FOUND=\"${TREXIO_FOUND}\")
endif(TREXIO_FOUND)

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
   string(REGEX REPLACE ".*\nflags[ \t]*:[ \t]+([^\n]+).*" "\\1" _cpu_flags "${_cpuinfo}")
endif(CMAKE_SYSTEM_NAME STREQUAL "Linux")

string( TOLOWER "${VECTORIZED}" VECTORIZED )
message(STATUS "Vectorization status  " ${VECTORIZED} )

if(VECTORIZED STREQUAL "yes")
	message(STATUS "Vectorization switched ON by user  " )
	list(APPEND Fortran_FLAGS "-xCORE-AVX512")
	set(TARGET_ARCHITECTURE "avx512")
	add_definitions(-DVECTORIZATION="avx512")
elseif(VECTORIZED STREQUAL "no")
	message(STATUS "Vectorization switched off by user  " )
else()
	message(STATUS "Vectorization is default  " )
	if(_cpu_flags MATCHES "avx512")
  		list(APPEND Fortran_FLAGS "-xCORE-AVX512")
  		set(TARGET_ARCHITECTURE "avx512")
  		add_definitions(-DVECTORIZATION="avx512")
	elseif(_cpu_flags MATCHES "avx2")
      if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  		  list(APPEND Fortran_FLAGS "-xCORE-AVX2")
  			set(TARGET_ARCHITECTURE "avx2")
  			add_definitions(-DVECTORIZATION="avx2")
		  else()
  			list(APPEND Fortran_FLAGS "-axCORE-AVX2")
  			set(TARGET_ARCHITECTURE "avx2")
  			add_definitions(-DVECTORIZATION="avx2")
		endif()
	endif()
endif()
unset(VECTORIZED CACHE)

add_definitions(-DTARGET_ARCHITECTURE=\"${TARGET_ARCHITECTURE}\")
## Architecture specific flags end here

# system information
include(${CMAKE_MODULE_PATH}/systeminfo.cmake)

# Git information
include(${CMAKE_MODULE_PATH}/gitinfo.cmake)
gitinfo(${CMAKE_SOURCE_DIR})

## Compiler FLAGS
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  list(APPEND Fortran_FLAGS "-O2"  "-cpp" "-mcmodel=large" "-ffree-line-length-132")
  list(APPEND Fortran_FLAGS "-D_MPI_")
  list(APPEND Fortran_FLAGS "-DCLUSTER")
  if (TREXIO_FOUND)
    list(APPEND Fortran_FLAGS "-ltrexio")
  endif()

  list(APPEND Fortran_FLAGS "-g" "-O0" "-debug" "-traceback")

  set(CMAKE_Fortran_FORMAT_FIXED_FLAG "-ffixed-form -ffixed-line-length-132 -fno-automatic")

  if (MKL_FOUND)
    #list(APPEND Fortran_FLAGS "-fdefault-integer-8")
    list(APPEND Fortran_FLAGS "-m64")
  endif()
  #set(CMAKE_Fortran_FORMAT_FREE_FLAG "-ffree-form")

elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  list(APPEND Fortran_FLAGS "-O2")
	list(APPEND Fortran_FLAGS "-fpp" "-mcmodel=small" "-shared-intel" "-dyncom=grid3d_data,orbital_num_spl,orbital_num_lag,orbital_num_spl2,grid3d_data" )
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
add_subdirectory(src)
add_subdirectory(tools)
add_subdirectory(lib)

## Testing
MESSAGE(STATUS " Tests are enabled ::  cd build; ctest")
enable_testing()
add_subdirectory(tests)

#Packaging
set(CPACK_SOURCE_GENERATOR "TBZ2")
set(CPACK_PACKAGE_VERSION "2.0.5")
set(CPACK_PACKAGE_NAME "CHAMP")
set(CPACK_ARCHIVE_FILE_NAME "CHAMP-v2.0.5.tar.bz2")
set(CPACK_PACKAGE_VERSION "v2.0.5")
set(CPACK_PACKAGE_VERSION_MAJOR "2")
set(CPACK_PACKAGE_VERSION_MINOR "0")
set(CPACK_PACKAGE_VERSION_PATCH "5")
include(CPack)

