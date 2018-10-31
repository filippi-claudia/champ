if (MPI_FOUND)
  file(GLOB DMC_MPI_SOURCE LIST_DIRECTORIES false  *.f)

  # Include function to change fortran format
  include(${CMAKE_MODULE_PATH}/FortranFormat.cmake)
  formatFortran("${DMC_MPI_SOURCE}" "FIXED")
  
  add_executable(dmc.mov1.mpi
    # VMC dependencies
    ${VMC_SOURCE}
    # MPI VMC code dependencies
    ${VMC_MPI_SOURCE}
    # code in the current directory 
    ${DMC_MPI_SOURCE}
    #Input parser
    ${INPUT_DIR}/p2_defaults.f ${INPUT_DIR}/p2etc.f ${INPUT_DIR}/p2nmcheck.f)
  
  # Add dependency to vmc mpi
  add_dependencies(dmc.mov1.mpi vmc.mov1.mpi dmc.mov1)

    #  Add headers
  target_include_directories(dmc.mov1.mpi
    PRIVATE
    ${CMAKE_SOURCE_DIR}/include)


  
  # compilation
  target_compile_options(dmc.mov1
    PRIVATE
    ${Fortran_FLAGS}
    "$<$<CONFIG:Debug>:${Fortan_FLAGS_DEBUG}>")

  # Link dependencies
  target_link_libraries(dmc.mov1.mpi
    # PRIVATE
    cyrus pspline
    # PUBLIC
    ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES}
    MPI::MPI_Fortran)

endif()
