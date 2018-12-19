if (MPI_Fortran_FOUND)
  file(GLOB DMC_MPI_SOURCE LIST_DIRECTORIES false  *.f)

  # Include function to change fortran format
  include(${CMAKE_MODULE_PATH}/FortranFormat.cmake)
  formatFortran("${DMC_MPI_SOURCE}" "FIXED")

  add_executable(dmc.mov1.mpi
    $<TARGET_OBJECTS:shared_objects_1>
    $<TARGET_OBJECTS:shared_objects_2>
    $<TARGET_OBJECTS:vmc_dmc_mpi_objects>
    $<TARGET_OBJECTS:dmc_mpi_objects>
    # code in the current directory 
    ${DMC_MPI_SOURCE})
  
    #  Add headers
  target_include_directories(dmc.mov1.mpi
    PRIVATE
    ${CMAKE_SOURCE_DIR}/include)

  # compilation
  target_compile_options(dmc.mov1.mpi
    PRIVATE
    ${Fortran_FLAGS}
    "$<$<CONFIG:Debug>:${Fortan_FLAGS_DEBUG}>")

  # Link dependencies
  target_link_libraries(dmc.mov1.mpi
    # PRIVATE
    cyrus pspline
    # PUBLIC
    ${LINEAR_ALGEBRA}
    ${MPI_Fortran_LIBRARIES})

endif()
