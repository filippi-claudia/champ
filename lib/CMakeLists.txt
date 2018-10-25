
file(GLOB_RECURSE CHAMP_SOURCES *.f)

add_library(cyrus STATIC ${CHAMP_SOURCES})

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
