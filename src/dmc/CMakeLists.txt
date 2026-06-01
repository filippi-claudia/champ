# Source files
file(GLOB DMC_SOURCE
  acues1_reduce.f90
  assignment.f90
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
  dmc_ae_mov1.f90
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

# DMC executable
add_executable(dmc.mov1
  $<TARGET_OBJECTS:shared_modules>
  $<TARGET_OBJECTS:shared_objects>
  $<TARGET_OBJECTS:rng>
  # DMC sources
  ${DMC_SOURCE})

# if HDF5 is found compile with HDF5 support
if(HDF5_FOUND)
  target_sources(dmc.mov1 PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/dmc_store_hdf5.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/dmc_restore_hdf5.f90)

  set_property(
    SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/dmc_store_hdf5.f90
    APPEND
    PROPERTY COMPILE_DEFINITIONS
    COMPILER=\"${CMAKE_Fortran_COMPILER}\"
    COMPILER_VERSION=\"${CMAKE_Fortran_COMPILER_VERSION}\"
    GIT_HEAD_BRANCH=\"${GIT_HEAD_BRANCH}\"
    GIT_REVISION_HASH=\"${GIT_REVISION_HASH}\"
    HDF5_VERSION=\"${HDF5_VERSION}\"
  )

  set_property(
    SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/dmc.f90 ${CMAKE_CURRENT_SOURCE_DIR}/mc_configs.f90
    APPEND
    PROPERTY COMPILE_DEFINITIONS
    HDF5_FOUND=\"${HDF5_FOUND}\"
  )

  set_source_files_properties(
    ${CMAKE_CURRENT_SOURCE_DIR}/dmc_store_hdf5.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/dmc_restore_hdf5.f90
    APPEND PROPERTY COMPILE_FLAGS "-lhdf5_fortran -lhdf5"
  )
endif()

# Add headers
target_include_directories(dmc.mov1
  PUBLIC
  ${CMAKE_BINARY_DIR}/src/vmc
  ${CMAKE_BINARY_DIR}/src/module
)

# Compile options
target_apply_champ_fortran_compile_options(dmc.mov1)
target_apply_champ_fortran_link_options(dmc.mov1)

target_link_champ_runtime(dmc.mov1)
