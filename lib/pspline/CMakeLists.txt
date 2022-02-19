file(GLOB PSPLINE_SOURCES LIST_DIRECTORIES false  *.f)

add_library(pspline ${PSPLINE_SOURCES})

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
target_compile_options(pspline
  PRIVATE
  ${Fortran_FLAGS} "-fno-implicit-none" "-Wno-implicit-interface")

elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  set(N_Fortran_FLAGS "${Fortran_Flags}")
  list(REMOVE_ITEM N_Fortran_FLAGS "-warn" "all")
target_compile_options(pspline
  PRIVATE
  ${N_Fortran_FLAGS} )
endif()

