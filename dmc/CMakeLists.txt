# Source files
file(GLOB DMC_SOURCE LIST_DIRECTORIES false  *.f)

include(${CMAKE_MODULE_PATH}/FortranFormat.cmake)
formatFortran("${DMC_SOURCE}") # Set Fixed format property

# --------------------------------------------

# executable dependencies
add_executable(dmc.mov1
  $<TARGET_OBJECTS:shared_objects_1>
  $<TARGET_OBJECTS:shared_objects_2>
  $<TARGET_OBJECTS:vmc_dmc_mpi_objects>
  # DMC source code
  ${DMC_SOURCE})

#  Add headers
target_include_directories(dmc.mov1
  PUBLIC
  ${CMAKE_SOURCE_DIR}/include)

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
