if (MPI_FOUND)
  file(GLOB VMC_MPI_SOURCE LIST_DIRECTORIES false  *.f)
  
  add_executable(vmc.mov1.mpi ${CMAKE_SOURCE_DIR}/include/vmc.h  ${VMC_MPI_SOURCE}
    ${INPUT_DIR}/p2_defaults.f ${INPUT_DIR}/p2etc.f ${INPUT_DIR}/p2nmcheck.f)

  # Add depency with vmc
  add_dependencies(vmc.mov1.mpi vmc.mov1)
  
  ## Set Fixed or Free format
  # Include function to change fortran format
  include(${CMAKE_MODULE_PATH}/FortranFormat.cmake)
  formatFortran("${VMC_MPI_SOURCES}" "FIXED")
  
  #  Add headers
  target_include_directories(vmc.mov1
    PRIVATE ${CMAKE_SOURCE_DIR}/include
    PRIVATE ${CMAKE_SOURCE_DIR}/input)
  
  # compilation
  target_compile_options(vmc.mov1.mpi
    PRIVATE
    ${Fortran_FLAGS}
    "$<$<CONFIG:Debug>:${Fortan_FLAGS_DEBUG}>"
    )
  
  target_link_libraries(vmc.mov1.mpi cyrus ${LAPACK_LIBRARIES} pspline ${BLAS_LIBRARIES})
endif()
