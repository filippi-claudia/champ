## Ravindra Shinde (2021) (c) TREX-CoE

cmake_minimum_required(VERSION 3.17 FATAL_ERROR)

project(CHAMP VERSION 2.5.0 LANGUAGES Fortran C)


# Allow targets to link to targets defined in other directories.
cmake_policy(SET CMP0079 NEW)
cmake_policy(SET CMP0076 NEW)

set(CMAKE_VERBOSE_MAKEFILE OFF)
set(Fortran_FLAGS_DEBUG "")

# CMake modules/macros are in a subdirectory to keep this file cleaner
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules)

# Print the header text art
find_package(Banner)

# Build summary
message(STATUS "Compiling                        :: " ${CMAKE_PROJECT_NAME})
message(STATUS "Fortran compiler                 :: " ${CMAKE_Fortran_COMPILER} )
message(STATUS "Fortran compiler version         :: " ${CMAKE_Fortran_COMPILER_VERSION} )

# Set folder variables:
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(CMAKE_PARSER_DIR ${CMAKE_SOURCE_DIR}/src/parser)

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})

# User options
option(ENABLE_QMMM "enable_qmmm" off)
option(ENABLE_TREXIO "enable_trexio" off)
option(ENABLE_QMCKL  "enable_qmckl" off)
option(ENABLE_GPU  "enable_gpu" off)

message("")

if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  if(BLA_STATIC MATCHES "ON")
    set(BLA_VENDOR Intel10_64lp)
  else()
    set(BLA_VENDOR Intel10_64_dyn)
  endif()
endif()

find_package(BLAS COMPONENTS Fortran REQUIRED QUIET)
find_package(LAPACK COMPONENTS Fortran REQUIRED QUIET)

if(BLAS_FOUND AND LAPACK_FOUND)
  message(STATUS "Using BLAS and LAPACK for the linear algebra calculations!")
  # Printing individual libraries in the list
  message(STATUS "BLAS and LAPACK LIBRARIES          : ")
  foreach(lib ${BLAS_LIBRARIES})
    message("                                    :: " ${lib})
  endforeach()
endif()

set(LINEAR_ALGEBRA ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
set(LINEAR_ALGEBRA_LINK_LIBRARIES ${LINEAR_ALGEBRA})
set(LINEAR_ALGEBRA_LINK_OPTIONS "")

if(APPLE)
  set(LINEAR_ALGEBRA_LINK_LIBRARIES "")
  foreach(lib ${LINEAR_ALGEBRA})
    if(lib MATCHES "\\.framework$")
      get_filename_component(framework_name "${lib}" NAME_WE)
      list(APPEND LINEAR_ALGEBRA_LINK_OPTIONS "SHELL:-framework ${framework_name}")
    else()
      list(APPEND LINEAR_ALGEBRA_LINK_LIBRARIES "${lib}")
    endif()
  endforeach()
  list(REMOVE_DUPLICATES LINEAR_ALGEBRA_LINK_LIBRARIES)
  list(REMOVE_DUPLICATES LINEAR_ALGEBRA_LINK_OPTIONS)
endif()

function(target_link_linear_algebra target)
  target_link_libraries(${target}
    PRIVATE
    ${LINEAR_ALGEBRA_LINK_LIBRARIES})
  target_link_options(${target}
    PRIVATE
    ${LINEAR_ALGEBRA_LINK_OPTIONS})
endfunction()

# LibFDF is built-in in the src/parser folder

# Find HDF5 library
if(ENABLE_TREXIO)
  find_package(HDF5 QUIET COMPONENTS Fortran)
  if(HDF5_FOUND)
    message(STATUS "HDF5 library version             :: " ${HDF5_VERSION})
    message(STATUS "HDF5 library parallel            :: " ${HDF5_IS_PARALLEL})
    message(STATUS "HDF5 C compiler executable       :: " ${HDF5_C_COMPILER_EXECUTABLE})
    message(STATUS "HDF5 include dirs                :: " ${HDF5_INCLUDE_DIRS})
    message(STATUS "HDF5 library directories         :  " )
    foreach(lib ${HDF5_LIBRARIES})
      message("                                    :: " ${lib})
    endforeach()
    message(STATUS "                                    ")
  endif()

  # Find TREXIO
  message(" ")
  message("Looking for TREXIO library:")
  find_package(TREXIO QUIET)

  if(TREXIO_FOUND)
    message(STATUS "TREXIO library found             :: " ${TREXIO_FOUND})
    message(STATUS "TREXIO Library include dirs      :: " ${TREXIO_INCLUDE_DIR})
    message(STATUS "TREXIO library                   :: " ${TREXIO_LIBRARY})
    message(STATUS "                                    ")
  endif()
endif()


# Find QMCkl
if(TREXIO_FOUND AND ENABLE_QMCKL)
  message(" ")
  message("Looking for QMCkl library:")
  find_package(QMCKL QUIET)
endif()
if(QMCKL_FOUND)
  message(STATUS "QMCkl library found              :: " ${QMCKL_FOUND})
  message(STATUS "QMCkl include dirs               :: " ${QMCKL_INCLUDE_DIR})
  message(STATUS "QMCkl library                    :: " ${QMCKL_LIBRARY})
  message("                                    ")
endif()

function(target_link_optional_libraries target)
  if(TREXIO_FOUND)
    target_compile_definitions(${target}
      PRIVATE
      TREXIO_FOUND="${TREXIO_FOUND}")
    target_include_directories(${target}
      PRIVATE
      ${TREXIO_INCLUDE_DIRS})
    target_link_libraries(${target}
      PRIVATE
      ${TREXIO_LIBRARIES})
  endif()

  if(QMCKL_FOUND)
    target_compile_definitions(${target}
      PRIVATE
      QMCKL_FOUND="${QMCKL_FOUND}")
    target_include_directories(${target}
      PRIVATE
      ${QMCKL_INCLUDE_DIRS})
    target_link_libraries(${target}
      PRIVATE
      ${QMCKL_LIBRARIES})
  endif()

  if(HDF5_FOUND)
    target_include_directories(${target}
      PRIVATE
      ${HDF5_INCLUDE_DIRS})
    target_link_libraries(${target}
      PRIVATE
      ${HDF5_LIBRARIES})
    if(TARGET HDF5::HDF5_Fortran)
      target_link_libraries(${target}
        PRIVATE
        HDF5::HDF5_Fortran)
    endif()
  endif()
endfunction()

function(target_link_champ_runtime target)
  target_link_libraries(${target}
    PRIVATE
    pspline
    parser
    MPI::MPI_Fortran)
  if(NOT TREXIO_FOUND)
    target_link_libraries(${target}
      PRIVATE
      ${CUDA_LIBRARIES})
  endif()
  target_link_linear_algebra(${target})
  target_link_optional_libraries(${target})
endfunction()

set(NVFORTRAN_PATH "" CACHE PATH "Path to CUDA executable nvfortran and nvcc")
if(ENABLE_GPU)
  find_package(CUDA)
endif()

if(CUDA_FOUND)
  message(" ")
  message(STATUS "CUDA found                      :: " ${CUDA_FOUND})
  message(STATUS "CUDA version                    :: " ${CUDA_VERSION})
  message(STATUS "CUDA toolkit directory          :: " ${CUDA_TOOLKIT_ROOT_DIR})
endif()
if(NOT DEFINED CUDA_LIBRARIES)
  set(CUDA_LIBRARIES "")
endif()
if(NOT DEFINED HDF5_VERSION)
  set(HDF5_VERSION "")
endif()

if(ENABLE_GPU)
  list(APPEND CMAKE_PREFIX_PATH ${NVFORTRAN_PATH})
endif()

# Find MPI
find_package(MPI COMPONENTS Fortran REQUIRED)

# SIMD instructions
find_package(SIMD)

function(target_apply_common_compile_definitions target)
  if(CHAMP_COMPILE_DEFINITIONS)
    target_compile_definitions(${target}
      PRIVATE
      ${CHAMP_COMPILE_DEFINITIONS})
  endif()
endfunction()

function(target_apply_champ_fortran_compile_options target)
  target_compile_options(${target}
    PRIVATE
    ${Fortran_FLAGS}
    "$<$<CONFIG:Debug>:${Fortran_FLAGS_DEBUG}>")
  target_apply_common_compile_definitions(${target})
endfunction()

function(target_apply_champ_fortran_link_options target)
  if(CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC" OR CMAKE_Fortran_COMPILER_ID MATCHES "Fujitsu")
    target_link_options(${target}
      PRIVATE
      ${Fortran_FLAGS}
      "$<$<CONFIG:Debug>:${Fortran_FLAGS_DEBUG}>")
  endif()
endfunction()

# system information
find_package(SystemInfo)

# Git information
find_package(GitInfo)
gitinfo(${CMAKE_SOURCE_DIR})

list(APPEND Fortran_FLAGS "-D_MPI_" "-DCLUSTER")

# List of Fortran flags for Nvidia's NVHPC compiler begins here
if(CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC")
  list(APPEND Fortran_FLAGS "-O2" "-cpp" "-mcmodel=medium")
  list(APPEND Fortran_FLAGS "-DPGI")
  list(APPEND Fortran_FLAGS "-DOPENACC")
  list(APPEND Fortran_FLAGS "-Mfree")
  list(APPEND Fortran_FLAGS "-mp=multicore,gpu")
  list(APPEND Fortran_FLAGS "-acc")
  list(APPEND Fortran_FLAGS "-gpu=cc80")
  list(APPEND Fortran_FLAGS "-Mcudalib=cublas")
  list(APPEND Fortran_FLAGS "-Mcuda=lineinfo")
  list(APPEND Fortran_FLAGS "-Minfo=mp,acc")
  list(APPEND Fortran_FLAGS "-gopt")
  list(APPEND Fortran_FLAGS "-traceback")
  list(APPEND Fortran_FLAGS "-Mextend")
  list(APPEND Fortran_FLAGS "-Mnoframe")
  if(MKL_FOUND)
    list(APPEND Fortran_FLAGS "-m64")
  endif()
endif()
# List of Fortran flags for Nvidia's NVHPC compiler ends here

# List of Fortran flags for Fujitsu compiler begins here
if(CMAKE_Fortran_COMPILER_ID MATCHES "Fujitsu")
  message(STATUS "Fujitsu compiler detected :: ")
  list(APPEND Fortran_FLAGS "-O3" "-Cpp" "-Nalloc_assign" "-Fwide" "-fs" "-Kopenmp" "-Kparallel" "-Kfast" )
endif()
# List of Fortran flags for Fujitsu compiler ends here

# List of Fortran flags for GNU compiler begins here
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  if(APPLE AND CMAKE_HOST_SYSTEM_PROCESSOR MATCHES "arm64")
    list(APPEND Fortran_FLAGS "-O2" "-cpp" "-ffree-line-length-none")
  else()
    list(APPEND Fortran_FLAGS "-O2" "-cpp" "-mcmodel=large" "-ffree-line-length-none")
  endif()
  if(CMAKE_BUILD_TYPE MATCHES "^[Dd][Ee][Bb][Uu][Gg]$")
    list(APPEND Fortran_FLAGS "-fcheck=all" "-fbacktrace"
        "-g" "-Wall" "-Wextra" "-fPIC"
        "-Wno-unused" "-Wno-conversion" "-Wno-unused-dummy-argument"
        "-Wno-character-truncation" "-Wno-compare-reals"
        "-Wno-function-elimination"
        "-fno-check-array-temporaries"
        "-std=f2018"
        "-std=legacy"
        # "-ffpe-trap=invalid,zero,overflow"
        )
  endif()

  if(MKL_FOUND)
    list(APPEND Fortran_FLAGS "-m64")
  endif()
  # List of Fortran flags for GNU compiler ends here

# List of Fortran flags for Intel's ifort compiler begins here
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
  list(APPEND Fortran_FLAGS "-O2" "-diag-disable=all")
  list(APPEND Fortran_FLAGS "-fPIC")
  list(APPEND Fortran_FLAGS "-implicitnone")
  list(APPEND Fortran_FLAGS "-ipo")
  list(APPEND Fortran_FLAGS "-finline" "-ip" "-align" "array64byte" "-fma" "-ftz" "-fomit-frame-pointer")
  list(APPEND Fortran_FLAGS "-fpp" "-mcmodel=small" "-shared-intel" "-dyncom=grid3d_data,orbital_num_spl,orbital_num_lag,orbital_num_spl2,grid3d_data" )
  if(CMAKE_BUILD_TYPE MATCHES "^[Dd][Ee][Bb][Uu][Gg]$")
    list(APPEND Fortran_FLAGS "-O0" "-diag-enable=all" "-traceback" "-g" "-check" "all" )
  endif()
  # List of Fortran flags for Intel's ifort compiler ends here

# List of Fortran flags for Intel's LLVM based compiler begins here
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM")
  list(APPEND Fortran_FLAGS "-O2")
  list(APPEND Fortran_FLAGS "-fPIC")
  list(APPEND Fortran_FLAGS "-implicitnone")
  list(APPEND Fortran_FLAGS "-finline" "-align" "array64byte" "-fma" "-ftz" "-fomit-frame-pointer")
  list(APPEND Fortran_FLAGS "-fpp" "-mcmodel=small" "-shared-intel" "-dyncom=grid3d_data,orbital_num_spl,orbital_num_lag,orbital_num_spl2,grid3d_data" )
  if(CMAKE_BUILD_TYPE MATCHES "^[Dd][Ee][Bb][Uu][Gg]$")
    list(APPEND Fortran_FLAGS "-O0" "-diag-enable=all" "-traceback" "-g" "-check" "all" )
  endif()
  # List of Fortran flags for Intel's LLVM based compiler ends here

# List of Fortran flags for LLVM "classic" Flang based compiler begins here
# This includes flang and aocc-flang and shouldn't be confused with flang-new (LLVMFlang in CMake)
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Flang")
  list(APPEND Fortran_FLAGS "-O2" "-cpp")
# List of Fortran flags for LLVM Flang based compiler ends here

elseif(CMAKE_Fortran_COMPILER_ID MATCHES "PGI")
  list(APPEND CMAKE_Fortran_FLAGS "-Mfreeform -Mdclchk -Mstandard -Mallocatable=03")
endif()

message("")
message(STATUS "Fortran compiler flags   : " )
foreach(flag ${Fortran_FLAGS})
  message("                                    :: " ${flag})
endforeach()


# Subdirectories
add_subdirectory(src)
add_subdirectory(tools)
add_subdirectory(lib)

# Testing
message(STATUS "Tests are enabled :: cd build; ctest")
enable_testing()
add_subdirectory(tests)

# Packaging
set(CPACK_ARCHIVE_COMPONENT_INSTALL 1)
set(CPACK_ARCHIVE_GENERATOR "TBZ2")
set(CPACK_SOURCE_GENERATOR "ZIP")
set(CPACK_PACKAGE_NAME "${PROJECT_NAME}")
set(CPACK_PACKAGE_VERSION "${PROJECT_VERSION}")
set(CPACK_PACKAGE_VERSION_MAJOR "${PROJECT_VERSION_MAJOR}")
set(CPACK_PACKAGE_VERSION_MINOR "${PROJECT_VERSION_MINOR}")
set(CPACK_PACKAGE_VERSION_PATCH "${PROJECT_VERSION_PATCH}")
set(CPACK_ARCHIVE_FILE_NAME "${PROJECT_NAME}-v${PROJECT_VERSION}")
set(CPACK_SOURCE_PACKAGE_FILE_NAME "${PROJECT_NAME}-v${PROJECT_VERSION}")
include(CPack)
