
file(GLOB_RECURSE CHAMP_SOURCES *.f *.f90)

add_library(cyrus STATIC ${CHAMP_SOURCES})
