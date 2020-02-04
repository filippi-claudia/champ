
file(GLOB_RECURSE CHAMP_SOURCE *.f)

add_library(cyrus STATIC ${CHAMP_SOURCE})

### Set Fixed or Free format
## Include function to change fortran format
#include(${CMAKE_MODULE_PATH}/FortranFormat.cmake)
#
#formatFortran("${CHAMP_SOURCE}" "FIXED")
foreach(_source ${CHAMP_SOURCE})
  set_property(SOURCE ${_source} PROPERTY Fortran_FORMAT FIXED)
endforeach(_source)


# set compiler flags
target_compile_options(cyrus
  PRIVATE
  ${Fortran_FLAGS}
  "$<$<CONFIG:Debug>:${Fortan_FLAGS_DEBUG}>"
  )
  # "$<$<CONFIG:Release>:${Fortan_FLAGS_RELEASE}>"

install(
  TARGETS cyrus
  ARCHIVE DESTINATION ${LIBRARY_OUTPUT_PATH}
  LIBRARY DESTINATION ${LIBRARY_OUTPUT_PATH})
