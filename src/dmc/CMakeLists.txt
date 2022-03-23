# Source files
file(GLOB DMC_SOURCE "*.f")

find_package(FortranFormat)
formatFortran("${DMC_SOURCE}") # Set Fixed format property

#-------------------------------------------------------------
# Compile shared library for VMC
#-------------------------------------------------------------
formatFortran("${SHARED_OBJECTS_F}") # Set Fixed format property
add_library(shared_objects_dmc OBJECT ${SHARED_OBJECTS})

#  Add headers
target_include_directories(shared_objects_dmc
  PUBLIC
  ${CMAKE_SOURCE_DIR}/build/src/module
  ${CMAKE_BINARY_DIR}/src/parser
)

target_compile_options(shared_objects_dmc
  PRIVATE
  ${Fortran_FLAGS}
  "$<$<CONFIG:Debug>:${Fortan_FLAGS_DEBUG}>"
)

add_dependencies(shared_objects_dmc shared_modules)

# --------------------------------------------------------------------

# executable dependencies
add_executable(dmc.mov1
  $<TARGET_OBJECTS:shared_modules>
  $<TARGET_OBJECTS:shared_objects_dmc>
  # DMC source code
  ${DMC_SOURCE})

#  Add headers
target_include_directories(dmc.mov1
  PUBLIC
  ${CMAKE_BINARY_DIR}/src/parser
  ${CMAKE_BINARY_DIR}/src/shared
  ${CMAKE_BINARY_DIR}/src/module
  )

# compilation
target_compile_options(dmc.mov1
  PRIVATE
  ${Fortran_FLAGS}
  "$<$<CONFIG:Debug>:${Fortan_FLAGS_DEBUG}>"
  )

if (TREXIO_FOUND)
  target_link_libraries(dmc.mov1
    PRIVATE
    pspline
    parser
    ${LINEAR_ALGEBRA}
    ${TREXIO_LIBRARY}
    MPI::MPI_Fortran)
else()
  target_link_libraries(dmc.mov1
    PRIVATE
    pspline
    parser
    ${LINEAR_ALGEBRA}
    MPI::MPI_Fortran)
endif()

set(DMC_SOURCE ${DMC_SOURCE} CACHE INTERNAL "Path to the source of vmc f77")
