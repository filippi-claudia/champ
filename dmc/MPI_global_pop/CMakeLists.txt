if (MPI_Fortran_FOUND)
  list(APPEND DMC_MPI_GPOP_SOURCE main.f walksav_jas.f)

  # Remove unused dependencies from Global population
  set(PATH_DMC_MPI_GPOPB ${CMAKE_SOURCE_DIR}/dmc/MPI_global_pop_big)
  list(REMOVE_ITEM DMC_MPI_GPOPB_DEPENDENCIES ${PATH_DMC_MPI_GPOPB}/main.f ${PATH_DMC_MPI_GPOPB}/walksav_jas.f)
  
  # REmove unused files
  # Include function to change fortran format
  include(${CMAKE_MODULE_PATH}/FortranFormat.cmake)
  formatFortran("${DMC_MPI_SOURCE}" "FIXED")

  # message(STATUS "SOURCE: " ${DMC_MPI_GPOPB_DEPENDENCIES})
  add_executable(dmc.mov1.gpop.mpi
    # VMC Fortran 90
    ${CMAKE_SOURCE_DIR}/vmc/convcartzmat.f90
    # All the dependencies from Global pop big
    ${DMC_MPI_GPOPB_DEPENDENCIES}
    # code in the current directory 
    ${DMC_MPI_GPOP_SOURCE}
    #Input parser
    ${INPUT_DIR}/p2_defaults.f ${INPUT_DIR}/p2etc.f ${INPUT_DIR}/p2nmcheck.f)
  
  # Add dependency to vmc mpi
  add_dependencies(dmc.mov1.gpop.mpi dmc.mov1.gpopb.mpi)

    #  Add headers
  target_include_directories(dmc.mov1.gpop.mpi
    PRIVATE
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/input)


  # compilation
  target_compile_options(dmc.mov1.gpop.mpi
    PRIVATE
    ${Fortran_FLAGS}
    "$<$<CONFIG:Debug>:${Fortan_FLAGS_DEBUG}>")

  # Link dependencies
  target_link_libraries(dmc.mov1.gpop.mpi
    # PRIVATE
    cyrus pspline
    # PUBLIC
    ${LINEAR_ALGEBRA}
    ${MPI_Fortran_LIBRARIES})

endif()
