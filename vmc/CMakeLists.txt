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
  message(STATUS "VMC_SOURCES: " ${VMC_SOURCES})
endif()


## Compiler FLAGS
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  list(APPEND Fortran_FLAGS "-O2" "-ffixed-line-length-132" "-fno-automatic")
  
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  list(APPEND Fortran_FLAGS "-O2" "-xHost" "-mcmodel=large" "-shared-intel" "-132" "-save "
    "-dyncom=grid3d_data,orbital_num_spl,orbital_num_lag,orbital_num_spl2,grid3d_data"
    "-check bounds" "-check format ")
  list(APPEND Fortan_FLAGS_DEBUG "-g")
  
elseif(Fortran_COMPILER_ID MATCHES "PGI")
  list(APPEND CMAKE_Fortran_FLAGS "-Mfreeform -Mdclchk -Mstandard -Mallocatable=03")
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


# Set Fixed or Free format
message(STATUS "Setting source format properties")
  # Intel does not care about the fixed format
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  foreach(_source ${VMC90_SOURCES})
    set_source_files_properties(${_source} PROPERTIES FORTRAN_FORMAT FREE)
    message (STATUS "${_source} is Free format")
  endforeach()
  foreach(_source ${VMC_SOURCES})
    set_source_files_properties(${_source} PROPERTIES FORTRAN_FORMAT FIXED)
    message (STATUS "${_source} is Fixed format")
  endforeach()
endif()

# compilation
target_compile_options(vmc.mov1
  PRIVATE
  ${Fortran_FLAGS}
  "$<$<CONFIG:Debug>:${Fortan_FLAGS_DEBUG}>"
  )

target_link_libraries(vmc.mov1 cyrus ${LAPACK_LIBRARIES} pspline ${BLAS_LIBRARIES})
