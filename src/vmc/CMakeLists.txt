file(GLOB VMC_SOURCE_F "*.f")
file(GLOB VMC_SOURCE_F90 "*.f90")
set(VMC_SOURCE ${VMC_SOURCE_F} ${VMC_SOURCE_F90})

find_package(FortranFormat)
formatFortran("${VMC_SOURCE_F}") # Set Fixed format property

#-------------------------------------------------------------
# Compile shared library for VMC
#-------------------------------------------------------------
formatFortran("${SHARED_OBJECTS_F}") # Set Fixed format property
add_library(shared_objects_vmc OBJECT ${SHARED_OBJECTS})

#  Add headers
target_include_directories(shared_objects_vmc
  PUBLIC
  ${CMAKE_BINARY_DIR}/src/module
  ${CMAKE_BINARY_DIR}/src/parser
)

target_compile_options(shared_objects_vmc
  PRIVATE
  ${Fortran_FLAGS}
  "$<$<CONFIG:Debug>:${Fortan_FLAGS_DEBUG}>"
)

add_dependencies(shared_objects_vmc shared_modules)

# --------------------------------------------------------------------
# Compile and link executable:
# --------------------------------------------------------------------
add_executable(vmc.mov1
  $<TARGET_OBJECTS:shared_modules>
  $<TARGET_OBJECTS:shared_objects_vmc>
  # source that is not already included in the objects
  ${VMC_SOURCE})

#  Add headers
target_include_directories(vmc.mov1
  PUBLIC
  ${CMAKE_BINARY_DIR}/src/module
  ${CMAKE_BINARY_DIR}/src/shared
)


# compilation
target_compile_options(vmc.mov1
  PRIVATE
  ${Fortran_FLAGS}
  "$<$<CONFIG:Debug>:${Fortan_FLAGS_DEBUG}>"
)

if (TREXIO_FOUND)
  target_link_libraries(vmc.mov1
    PRIVATE
    pspline
    parser
    ${LINEAR_ALGEBRA}
    ${TREXIO_LIBRARY}
    MPI::MPI_Fortran)
else()
  target_link_libraries(vmc.mov1
    PRIVATE
    pspline
    parser
    ${LINEAR_ALGEBRA}
    MPI::MPI_Fortran)
endif()

