# the netherlands escience center (2022)

project(vmc)

option(UNIT_TESTS "Build Project Unit tests" OFF)
message(STATUS "Build Unit Tests: ${UNIT_TESTS}")

set(N_Fortran_FLAGS "${Fortran_FLAGS}")
list(REMOVE_ITEM N_Fortran_FLAGS "-Werror")

if(UNIT_TESTS)
  string(REPLACE ";" " " CMAKE_Fortran_FLAGS "${N_Fortran_FLAGS}")
  if(NOT CMAKE_BUILD_TYPE MATCHES "^[Dd][Ee][Bb][Uu][Gg]$")
    message(FATAL_ERROR "if UNIT_TESTS is specified, CMAKE_BUILD_TYPE must be Debug")
  endif()

  file(GLOB_RECURSE FORTUTF_PROJECT_SRC_FILES CONFIGURE_DEPENDS "*.f90" "*.F90")
  set(FORTUTF_PROJECT_TEST_DIR ${CMAKE_CURRENT_LIST_DIR})
  set(FORTUTF_PROJECT_SRC_LIBRARY "vmc;pspline;parser")
  set(FORTUTF_PROJECT_MOD_DIR "${CMAKE_BINARY_DIR}/src/module;${CMAKE_BINARY_DIR}/src/vmc")

  set(FORTUTF_ROOT ../FortUTF)
  include(${FORTUTF_ROOT}/cmake/fortutf.cmake)
  FortUTF_Find_Tests()
endif()
