if (MPI_Fortran_FOUND)
  file(GLOB VMC_MPI_SOURCE LIST_DIRECTORIES false  *.f)
  file(GLOB VMC90_MPI_SOURCE LIST_DIRECTORIES false  *.f90)

  # Remove unused files
  list(REMOVE_ITEM VMC_MPI_SOURCE
    ${CMAKE_CURRENT_SOURCE_DIR}/lin_d_more_no_semiortho.f
    ${CMAKE_CURRENT_SOURCE_DIR}/lin_d_more_semiortho.f
    )

  include(${CMAKE_MODULE_PATH}/FortranFormat.cmake)
  formatFortran("${VMC_MPI_SOURCE}") # Set Fixed format property

# --------------------------------------------
# creates object library for dmc.mov1.mpi
# --------------------------------------------

  set(SHARED_DMC_MPI_SOURCE
    ${CMAKE_CURRENT_SOURCE_DIR}/optx_jas_ci_reduce.f
    ${CMAKE_CURRENT_SOURCE_DIR}/optx_jas_orb_reduce.f
    ${CMAKE_CURRENT_SOURCE_DIR}/optx_orb_ci_reduce.f
    ${CMAKE_CURRENT_SOURCE_DIR}/optci_reduce.f
    ${CMAKE_CURRENT_SOURCE_DIR}/optjas_reduce.f
    ${CMAKE_CURRENT_SOURCE_DIR}/optorb_reduce.f
    ${CMAKE_CURRENT_SOURCE_DIR}/finwrt_more.f
    ${CMAKE_CURRENT_SOURCE_DIR}/error_mpi.f
    ${CMAKE_CURRENT_SOURCE_DIR}/optwf_write.f)

  add_library(vmc_dmc_mpi_objects OBJECT ${SHARED_DMC_MPI_SOURCE})

  target_include_directories(vmc_dmc_mpi_objects
    PUBLIC
    ${CMAKE_SOURCE_DIR}/include)

  target_compile_options(vmc_dmc_mpi_objects
    PRIVATE
    ${Fortran_FLAGS}
    "$<$<CONFIG:Debug>:${Fortan_FLAGS_DEBUG}>")

  foreach(s ${SHARED_DMC_MPI_SOURCE})
    list(REMOVE_ITEM VMC_MPI_SOURCE ${s})
  endforeach(s)

# --------------------------------------------
  
  add_executable(vmc.mov1.mpi
    # MPI source code
    ${VMC_MPI_SOURCE} ${VMC90_MPI_SOURCE}
    # shared VMC objects
    $<TARGET_OBJECTS:shared_objects_1>
    $<TARGET_OBJECTS:shared_objects_2>
    $<TARGET_OBJECTS:vmc_mpi_objects>
    $<TARGET_OBJECTS:vmc_dmc_mpi_objects>)

  
  #  Add headers
  target_include_directories(vmc.mov1.mpi
    PUBLIC
    ${CMAKE_SOURCE_DIR}/include)
  
  # compilation
  target_compile_options(vmc.mov1.mpi
    PRIVATE
    ${Fortran_FLAGS}
    "$<$<CONFIG:Debug>:${Fortan_FLAGS_DEBUG}>"
    )

  target_link_libraries(vmc.mov1.mpi
    # PRIVATE
    cyrus pspline
    # PUBLIC
    ${LINEAR_ALGEBRA}
    ${MPI_Fortran_LIBRARIES}
    )

  # # # Cache the results for targets that depends on these files
  set(VMC_MPI_SOURCE ${VMC_MPI_SOURCE} CACHE INTERNAL "MPI vmc source files")

endif()
