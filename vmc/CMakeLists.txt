# Source files
file(GLOB VMC_SOURCE LIST_DIRECTORIES false  *.f)
file(GLOB VMC90_SOURCE LIST_DIRECTORIES false coords_int.f90 optgeo_hessian.f90 convcartzmat.f90 regterg.f90)

if (NOT ENABLE_QMMM)
  message(STATUS "QMMM is NOT enable!")
endif()

if (NOT ENABLE_PERIODIC)
  message(STATUS "PERIODIC is NOT enable!")
  list(REMOVE_ITEM VMC_SOURCE
    ${CMAKE_CURRENT_SOURCE_DIR}/pw_find_image.f ${CMAKE_CURRENT_SOURCE_DIR}/pw_ewald.f
    ${CMAKE_CURRENT_SOURCE_DIR}/pw_orbitals.f ${CMAKE_CURRENT_SOURCE_DIR}/pw_orbitals_e.f
    ${CMAKE_CURRENT_SOURCE_DIR}/pw_read.f)
endif()

# Compile and link executable
add_executable(vmc.mov1 ${CMAKE_SOURCE_DIR}/include/vmc.h ${VMC_SOURCE} ${VMC90_SOURCE}
  ${INPUT_DIR}/p2_defaults.f ${INPUT_DIR}/p2etc.f ${INPUT_DIR}/p2nmcheck.f)

#  Add headers
target_include_directories(vmc.mov1
  PUBLIC
  ${CMAKE_SOURCE_DIR}/include
  ${CMAKE_SOURCE_DIR}/input
)


## Set Fixed or Free format
# Include function to change fortran format
#include(${CMAKE_MODULE_PATH}/FortranFormat.cmake)

#formatFortran("${VMC90_SOURCE}" "FREE")
#formatFortran("${VMC_SOURCE}" "FIXED")
foreach(_source ${VMC_SOURCE})
  set_property(SOURCE ${_source} PROPERTY Fortran_FORMAT FIXED)
endforeach(_source)

# compilation
target_compile_options(vmc.mov1
  PRIVATE
  ${Fortran_FLAGS}
  "$<$<CONFIG:Debug>:${Fortan_FLAGS_DEBUG}>"
  )

target_link_libraries(vmc.mov1
  cyrus pspline
  ${LINEAR_ALGEBRA})

# Include MPI subdir
set(VMC_SOURCE ${VMC_SOURCE} CACHE INTERNAL "Path to the source of vmc f77")
add_subdirectory(MPI)


# Everything related to miniparser p2
add_custom_command(
  OUTPUT
    ${INPUT_DIR}/p2nmcheck.f
  OUTPUT
    ${CMAKE_CURRENT_SOURCE_DIR}/p2prog.f
  WORKING_DIRECTORY
    ${CMAKE_CURRENT_SOURCE_DIR}
  COMMAND
    ${INPUT_DIR}/inpt.awk -v p2file=commands.p2 -v logf=${INPUT_DIR}/commands.log ${CMAKE_CURRENT_SOURCE_DIR}/read_input.f
  COMMAND
    cat commands.p2 ${INPUT_DIR}/end.p2 ${INPUT_DIR}/auto.p2 > all.p2
  COMMAND
    ${INPUT_DIR}/genp2_defaults.awk -v ofile=${CMAKE_CURRENT_SOURCE_DIR}/p2prog.f all.p2 > genp2.log
  COMMAND
    ${INPUT_DIR}/mkp2nmcheck.pl namelist.p2 > ${INPUT_DIR}/p2nmcheck.f
  VERBATIM
)

add_custom_target(
  generate-p2
  ${CMAKE_COMMAND} -E echo "Generating p2 files."
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/p2prog.f
  DEPENDS ${INPUT_DIR}/p2nmcheck.f
)

add_dependencies(vmc.mov1 generate-p2)

set_source_files_properties(
  ${INPUT_DIR}/p2nmcheck.f
  PROPERTIES GENERATED TRUE
)

set_source_files_properties(
  ${CMAKE_CURRENT_SOURCE_DIR}/p2prog.f
  PROPERTIES GENERATED TRUE
)
