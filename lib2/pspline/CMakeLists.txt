file(GLOB PSPLINE_SOURCES LIST_DIRECTORIES false  *.f)

add_library(pspline ${PSPLINE_SOURCES})

## Compiler FLAGS
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  list(APPEND Fortran_FLAGS "-O2" "-ffixed-line-length-132" "-fno-automatic")
  
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  list(APPEND Fortran_FLAGS "-O2" "-xHost" "-mcmodel=large" "-shared-intel" "-132" "-save "
    "-check bounds" "-check format ")
endif()
 
target_compile_options(pspline
  PRIVATE
  ${Fortran_FLAGS})

