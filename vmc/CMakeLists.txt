# Source files
file(GLOB VMC_SOURCES LIST_DIRECTORIES false  *.f)
file(GLOB VMC90_SOURCES LIST_DIRECTORIES false convcartzmat.f90 regterg.f90)

if (NOT ENABLE_QMMM)
  message(STATUS "QMMM is NOT enable!")
endif()

if (NOT ENABLE_PERIODIC)
  message(STATUS "PERIODIC is NOT enable!")
  list(REMOVE_ITEM VMC_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/pw_find_image.f ${CMAKE_CURRENT_SOURCE_DIR}/pw_ewald.f
    ${CMAKE_CURRENT_SOURCE_DIR}/pw_orbitals.f ${CMAKE_CURRENT_SOURCE_DIR}/pw_orbitals_e.f
    ${CMAKE_CURRENT_SOURCE_DIR}/pw_read.f)
endif()

# Input directory
set(INPUT_DIR ${CMAKE_SOURCE_DIR}/input)

# Compile and link executable
add_executable(vmc.mov1 ${CMAKE_SOURCE_DIR}/include/vmc.h ${VMC_SOURCES} ${VMC90_SOURCES}
  ${INPUT_DIR}/p2_defaults.f ${INPUT_DIR}/p2etc.f ${INPUT_DIR}/p2nmcheck.f)

#  Add headers
target_include_directories(vmc.mov1
  PRIVATE ${CMAKE_SOURCE_DIR}/include
  PRIVATE ${CMAKE_SOURCE_DIR}/input)


## Set Fixed or Free format
# Include function to change fortran format
include(${CMAKE_MODULE_PATH}/FortranFormat.cmake)

formatFortran("${VMC90_SOURCES}" "FREE")
formatFortran("${VMC_SOURCES}" "FIXED")

# compilation
target_compile_options(vmc.mov1
  PRIVATE
  ${Fortran_FLAGS}
  "$<$<CONFIG:Debug>:${Fortan_FLAGS_DEBUG}>"
  )

target_link_libraries(vmc.mov1 cyrus ${LAPACK_LIBRARIES} pspline ${BLAS_LIBRARIES})
