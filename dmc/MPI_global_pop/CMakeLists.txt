if (MPI_Fortran_FOUND)
  list(APPEND DMC_MPI_GPOP_SOURCE main.f walksav_jas.f)
  include(${CMAKE_MODULE_PATH}/FortranFormat.cmake)
  formatFortran("${DMC_MPI_GPOP_SOURCE}" "FIXED")

  # Remove unused dependencies from Global population
  set(PATH_DMC_MPI_GPOPB ${CMAKE_SOURCE_DIR}/dmc/MPI_global_pop_big)
  list(REMOVE_ITEM DMC_MPI_GPOPB_DEPENDENCIES ${PATH_DMC_MPI_GPOPB}/main.f ${PATH_DMC_MPI_GPOPB}/walksav_jas.f)
  formatFortran("${DMC_MPI_GPOPB_DEPENDENCIES}" "FIXED")
  
  # REmove unused files

  # message(STATUS "SOURCE: " ${DMC_MPI_GPOPB_DEPENDENCIES})
  add_executable(dmc.mov1.gpop.mpi
    $<TARGET_OBJECTS:shared_objects_1>
    $<TARGET_OBJECTS:vmc_dmc_mpi_objects>
    # All the dependencies from Global pop big
    ${DMC_MPI_GPOPB_DEPENDENCIES}
    # code in the current directory 
    ${DMC_MPI_GPOP_SOURCE})
  
  # Add dependency to vmc mpi
  add_dependencies(dmc.mov1.gpop.mpi dmc.mov1.gpopb.mpi)

    #  Add headers
  target_include_directories(dmc.mov1.gpop.mpi
    PRIVATE
    ${CMAKE_SOURCE_DIR}/include)

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
