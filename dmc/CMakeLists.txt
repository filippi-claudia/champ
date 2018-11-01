# Source files
set(DMC_SOURCE main.f dmc.f mc_configs.f gauss.f prop_dmc.f optwf_dmc.f pcm_dmc.f
  mmpol_dmc.f averages.f dmc_good_ps_mov1.f acuest.f 
  finwrt.f dumper.f splitj.f walksav_det.f walksav_jas.f hpsiedmc.f 
  nonloc_grid.f multideterminant_tmove.f optwf_dummy.f acues1_reduce.f)

# Remove the files that are not used from vmc
set(PATH_VMC ${CMAKE_SOURCE_DIR}/vmc)
list(REMOVE_ITEM VMC_SOURCE
  ${PATH_VMC}/acuest.f ${PATH_VMC}/acuest_reduce.f ${PATH_VMC}/acuest_write.f
  ${PATH_VMC}/deriv_jastrow.f ${PATH_VMC}/determinant_psig.f ${PATH_VMC}/dl.f
  ${PATH_VMC}/dl_iter.f ${PATH_VMC}/dl_more.f ${PATH_VMC}/dumper.f ${PATH_VMC}/dumper_more.f
  ${PATH_VMC}/fetch_parameters.f ${PATH_VMC}/fin_reduce.f ${PATH_VMC}/finwrt.f
  ${PATH_VMC}/finwrt_more.f ${PATH_VMC}/jdqz_lib.f ${PATH_VMC}/jdqz_main.f ${PATH_VMC}/lin_d.f
  ${PATH_VMC}/lin_d_more.f ${PATH_VMC}/main.f ${PATH_VMC}/mc_configs.f
  ${PATH_VMC}/metrop_mov1_slat.f ${PATH_VMC}/mmpol_vmc.f ${PATH_VMC}/multiple_states.f
  ${PATH_VMC}/optwf_vmc.f ${PATH_VMC}/pcm_vmc.f
  ${PATH_VMC}/prop_vmc.f ${PATH_VMC}/pw_ewald.f ${PATH_VMC}/pw_find_image.f ${PATH_VMC}/pw_orbitals.f
  ${PATH_VMC}/pw_orbitals_e.f ${PATH_VMC}/pw_read.f ${PATH_VMC}/sr.f ${PATH_VMC}/sr_more.f
  ${PATH_VMC}/vmc.f)

# executable dependencies
add_executable(dmc.mov1
  # shared code with VMC
  ${VMC_SOURCE}
  # Fortran 90 source
  ${PATH_VMC}/convcartzmat.f90
  # DMC source code
  ${DMC_SOURCE}  
  #Input parser
  ${INPUT_DIR}/p2_defaults.f ${INPUT_DIR}/p2etc.f ${INPUT_DIR}/p2nmcheck.f)

#  Add headers
target_include_directories(dmc.mov1
  PUBLIC
  ${CMAKE_SOURCE_DIR}/include
  ${CMAKE_SOURCE_DIR}/input)

# Add depency with vmc
add_dependencies(dmc.mov1 vmc.mov1)

## Set Fixed or Free format
# Include function to change fortran format
include(${CMAKE_MODULE_PATH}/FortranFormat.cmake)
formatFortran("${VMC90_DEPENDENCIES}" "FREE")
formatFortran("${VMC_DEPENDENCIES}" "FIXED")
formatFortran("${DMC_SOURCE}" "FIXED")

# compilation
target_compile_options(dmc.mov1
  PRIVATE
  ${Fortran_FLAGS}
  "$<$<CONFIG:Debug>:${Fortan_FLAGS_DEBUG}>"
  )

target_link_libraries(dmc.mov1 cyrus ${LAPACK_LIBRARIES} pspline ${BLAS_LIBRARIES})
add_subdirectory(MPI)
add_subdirectory(MPI_global_pop_big)
