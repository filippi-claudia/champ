file(GLOB PSPLINE_SOURCES LIST_DIRECTORIES false  *.f)

add_library(pspline ${PSPLINE_SOURCES})

target_compile_options(pspline
  PRIVATE
  ${Fortran_FLAGS} "-fno-implicit-none" "-Wno-implicit-interface")

