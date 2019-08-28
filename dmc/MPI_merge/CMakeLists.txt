if (MPI_Fortran_FOUND)

  file(GLOB DMC_MPI_SOURCE LIST_DIRECTORIES false  *.f)
  include(${CMAKE_MODULE_PATH}/FortranFormat.cmake)
  formatFortran("${DMC_MPI_SOURCE}" "FIXED")

    # Remove the files that are not used from dmc
  set(PATH_DMC ${CMAKE_SOURCE_DIR}/dmc)
  list(REMOVE_ITEM DMC_SOURCE ${PATH_DMC}/acuest.f ${PATH_DMC}/acues1.f ${PATH_DMC}/acues1_reduce.f 
    ${PATH_DMC}/dumper.f ${PATH_DMC}/init.f ${PATH_DMC}/main.f ${PATH_DMC}/finwrt.f 
    ${PATH_DMC}/walksav_det.f ${PATH_DMC}/walksav_jas.f ${PATH_DMC}/mc_configs.f)
  formatFortran("${DMC_SOURCE}" "FIXED")

  add_executable(dmc.mov1.mpi
    $<TARGET_OBJECTS:shared_objects_1>
    $<TARGET_OBJECTS:shared_objects_2>
    $<TARGET_OBJECTS:vmc_dmc_mpi_objects>
    #
    ${DMC_SOURCE}
    #
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

  set(DMC_MPI_DEPENDENCIES ${DMC_SOURCE} ${DMC_MPI_SOURCE} 
    CACHE INTERNAL "source files used by targets depending on dmc.mov1.gpop.mpi")  

endif()
