
file(GLOB_RECURSE CHAMP_SOURCES *.f *.f90)

add_library(cyrus STATIC ${CHAMP_SOURCES})

## Compiler FLAGS
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  list(APPEND Fortran_FLAGS "-O3" "-ffixed-form" "-ffixed-line-length-132" "-fno-automatic")
  

elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  list(APPEND Fortran_FLAGS "-O3 -xHost -mcmodel=large -shared-intel -132 -save -g -dyncom=grid3d_data,orbital_num_spl,orbital_num_lag,orbital_num_spl2,grid3d_data -check bounds -check format ")

elseif(Fortran_COMPILER_ID MATCHES "PGI")
  list(APPEND CMAKE_Fortran_FLAGS "-Mfreeform -Mdclchk -Mstandard -Mallocatable=03")
endif()

# set compiler flags
target_compile_options(cyrus
  PRIVATE
  ${Fortran_FLAGS})
  # "$<$<CONFIG:Debug>:${Fortan_FLAGS_DEBUG}>"
  # "$<$<CONFIG:Release>:${Fortan_FLAGS_RELEASE}>"

install(
  TARGETS cyrus
  ARCHIVE DESTINATION ${LIBRARY_OUTPUT_PATH}
  LIBRARY DESTINATION ${LIBRARY_OUTPUT_PATH})
