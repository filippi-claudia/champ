# Source files
file(GLOB DMC_SOURCE
  acues1_reduce.f90
  assignment.f90
  averages.f90
  branching.f90
  dumper.f90
  main.f90
  mmpol_reduce.f90
  nonloc_grid.f90
  pcm_reduce.f90
  redistribute.f90
  splitj.f90
  zerest.f90
  acues1.f90
  acuest.f90
  dmc.f90
  dumper_gpop.f90
  hpsiedmc.f90
  mc_configs.f90
  move_walker.f90
  optwf_dmc.f90
  prop_dmc.f90
  restart.f90
  walksav_det.f90
  acues1_gpop.f90
  acuest_gpop.f90
  dmc_ps_mov1.f90
  finwrt.f90
  init.f90
  mmpol_dmc.f90
  multideterminant_tmove.f90
  pcm_dmc.f90
  prop_reduce.f90
  restart_gpop.f90
  walksav_jas.f90
  )

# --------------------------------------------------------------------

# executable dependencies
add_executable(dmc.mov1
  $<TARGET_OBJECTS:shared_modules>
  $<TARGET_OBJECTS:shared_objects>
  $<TARGET_OBJECTS:rng>
  # DMC source code
  ${DMC_SOURCE})

#  Add headers
target_include_directories(dmc.mov1
  PUBLIC
  ${CMAKE_BINARY_DIR}/src/vmc
  ${CMAKE_BINARY_DIR}/src/module
  )

# compilation
target_compile_options(dmc.mov1
  PRIVATE
  ${Fortran_FLAGS}
  "$<$<CONFIG:Debug>:${Fortan_FLAGS_DEBUG}>"
  )

  if(CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC" OR CMAKE_Fortran_COMPILER_ID MATCHES "Fujitsu")
  	target_link_options(dmc.mov1
	  PRIVATE
	  ${Fortran_FLAGS}
	  "$<$<CONFIG:Debug>:${Fortan_FLAGS_DEBUG}>"
	  )
  endif()

if (TREXIO_FOUND)
  target_link_libraries(dmc.mov1
    PRIVATE
    pspline
    parser
    ${LINEAR_ALGEBRA}
    ${TREXIO_LIBRARY}
    MPI::MPI_Fortran)
else()
  target_link_libraries(dmc.mov1
    PRIVATE
    pspline
    parser
    ${LINEAR_ALGEBRA}
    ${CUDA_LIBRARIES}
    MPI::MPI_Fortran)
endif()


set(DMC_SOURCE ${DMC_SOURCE} CACHE INTERNAL "Path to the source of vmc f77")
