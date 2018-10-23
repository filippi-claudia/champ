
file(GLOB_RECURSE CHAMP_SOURCES *.f *.f90)

# Use the ifort definition of time as default
file(COPY "my_second.f.ifort" "my_second.f")

elseif(CMAKE_Fortran_COMPILER_ID MATCHES "XL")
  file(COPY "my_second.f.xlf" "my_second.f")
endif()

add_library(cyrus STATIC ${CHAMP_SOURCES})
