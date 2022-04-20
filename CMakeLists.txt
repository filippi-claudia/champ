## Ravindra Shinde (2021) (c) TREX-CoE

cmake_minimum_required(VERSION 3.20 FATAL_ERROR)

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
find_package(Banner)

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

if((CMAKE_Fortran_COMPILER_ID MATCHES "Intel") OR ENABLE_MKL)
  if (BLA_STATIC MATCHES "ON")
    set(BLA_VENDOR Intel10_64lp)
  else()
    set(BLA_VENDOR Intel10_64_dyn)
  endif()
endif()

find_package(BLAS COMPONENTS Fortran REQUIRED QUIET)
find_package(LAPACK COMPONENTS Fortran REQUIRED QUIET)

if (BLAS_FOUND AND LAPACK_FOUND)
  message(STATUS "Using BLAS and LAPACK for the linear algebra calculations!")
  # Printing individual libraries in the list
  message(STATUS "BLAS and LAPACK LIBRARIES          : ")
  foreach(lib ${BLAS_LIBRARIES})
    message("                                    :: " ${lib})
  endforeach(lib)
endif()

set(LINEAR_ALGEBRA ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
link_libraries(${LINEAR_ALGEBRA})

# LibFDF is built-in in the src/parser folder

# Find HDF5 library
find_package(HDF5 QUIET)
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
#  add_definitions(-DHDF5_VERSION=\"${HDF5_VERSION}\")
endif(HDF5_FOUND)

# Find TREXIO
if (HDF5_FOUND)
  find_package(TREXIO QUIET)
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

# SIMD instructions
find_package(SIMD)

# system information
find_package(SystemInfo)

# Git information
find_package(GitInfo)
gitinfo(${CMAKE_SOURCE_DIR})

## Compiler FLAGS
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  list(APPEND Fortran_FLAGS "-O3" "-funroll-loops" "-cpp" "-ffree-line-length-none" "-g"
       "-Werror=implicit-interface" "-fimplicit-none" #"-flto"
        "-Wall" "-Wextra"
        "-Wno-unused" "-Wno-conversion" "-Wno-unused-dummy-argument"
        "-Wno-character-truncation" "-Wno-compare-reals" "-Wno-tabs"
        "-Wno-c-binding-type" "-Wno-function-elimination")
  list(APPEND Fortran_FLAGS "-D_MPI_")
  list(APPEND Fortran_FLAGS "-DCLUSTER")
  if(CMAKE_BUILD_TYPE MATCHES "DEBUG")
    list(APPEND Fortran_FLAGS "-fcheck=all" "-fbacktrace"
        )
  endif()
  if (TREXIO_FOUND)
    list(APPEND Fortran_FLAGS "-ltrexio")
  endif()

  set(CMAKE_Fortran_FORMAT_FIXED_FLAG "-ffixed-form -ffixed-line-length-132 -fno-automatic")

  if (MKL_FOUND)
    #list(APPEND Fortran_FLAGS "-fdefault-integer-8")
    list(APPEND Fortran_FLAGS "-m64")
  endif()
  #set(CMAKE_Fortran_FORMAT_FREE_FLAG "-ffree-form")

elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  list(APPEND Fortran_FLAGS "-O3" "-g") 
	list(APPEND Fortran_FLAGS "-fpp" "-shared-intel" "-g" 
      "-warn" "all" "-diag-disable=8889,7712")
  list(APPEND Fortran_FLAGS "-D_MPI_")
  list(APPEND Fortran_FLAGS "-DCLUSTER")
  if(CMAKE_BUILD_TYPE MATCHES "DEBUG")
    list(APPEND Fortran_FLAGS "-traceback")
  else()
  endif()
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
set(CPACK_ARCHIVE_COMPONENT_INSTALL 1)
set(CPACK_ARCHIVE_GENERATOR "TBZ2")
set(CPACK_SOURCE_GENERATOR "ZIP")
set(CPACK_PACKAGE_VERSION "2.0.7")
set(CPACK_PACKAGE_NAME "CHAMP")
set(CPACK_ARCHIVE_FILE_NAME "CHAMP-v2.0.7.zip")
set(CPACK_PACKAGE_VERSION "v2.0.7")
set(CPACK_PACKAGE_VERSION_MAJOR "2")
set(CPACK_PACKAGE_VERSION_MINOR "0")
set(CPACK_PACKAGE_VERSION_PATCH "7")
include(CPack)
