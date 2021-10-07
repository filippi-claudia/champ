# Compile the trexio interface source files
file(GLOB TREXIO_SOURCE ${CMAKE_TREXIO_DIR}/trexio_f.f90)

add_library(trexio ${TREXIO_SOURCE})

target_compile_options(trexio
  PUBLIC
  ${Fortran_FLAGS})
