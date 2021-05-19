# Source files
file(GLOB DMC_SOURCE *.f *.f90)

include(${CMAKE_MODULE_PATH}/FortranFormat.cmake)
formatFortran("${DMC_SOURCE}") # Set Fixed format property

# --------------------------------------------------------------------

# executable dependencies
add_executable(dmc.mov1
  $<TARGET_OBJECTS:shared_modules>
  $<TARGET_OBJECTS:shared_objects>
  # DMC source code
  ${UTILS_SOURCE}
  ${DMC_SOURCE})

#  Add headers
target_include_directories(dmc.mov1
  PUBLIC
  ${CMAKE_SOURCE_DIR}/build/src/vmc
  )

# compilation
target_compile_options(dmc.mov1
  PRIVATE
  ${Fortran_FLAGS}
  "$<$<CONFIG:Debug>:${Fortan_FLAGS_DEBUG}>"
  )

target_link_libraries(dmc.mov1
  cyrus pspline
  ${LINEAR_ALGEBRA})

set(DMC_SOURCE ${DMC_SOURCE} CACHE INTERNAL "Path to the source of vmc f77")
