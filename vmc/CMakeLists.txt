# Source files
file(GLOB VMC_SOURCE LIST_DIRECTORIES false  *.f)
file(GLOB VMC90_SOURCE LIST_DIRECTORIES false convcartzmat.f90 regterg.f90)

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
  ${CMAKE_SOURCE_DIR}/input)


## Set Fixed or Free format
# Include function to change fortran format
include(${CMAKE_MODULE_PATH}/FortranFormat.cmake)

formatFortran("${VMC90_SOURCE}" "FREE")
formatFortran("${VMC_SOURCE}" "FIXED")

# compilation
target_compile_options(vmc.mov1
  PRIVATE
  ${Fortran_FLAGS}
  "$<$<CONFIG:Debug>:${Fortan_FLAGS_DEBUG}>"
  )

set(VAR "${MKLROOT}/lib/intel64/libmkl_blas95_ilp64.a ${MKLROOT}/lib/intel64/libmkl_lapack95_ilp64.a -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl")

target_link_libraries(vmc.mov1 cyrus pspline ${VAR})

# Include MPI subdir
set(VMC_SOURCE ${VMC_SOURCE} CACHE INTERNAL "Path to the source of vmc f77")
add_subdirectory(MPI)
