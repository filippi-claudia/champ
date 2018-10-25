# Source files
set(DMC_SOURCE main.f dmc.f mc_configs.f gauss.f prop_dmc.f optwf_dmc.f pcm_dmc.f
  mmpol_dmc.f averages.f dmc_good_ps_mov1.f acuest.f 
  finwrt.f dumper.f splitj.f walksav_det.f walksav_jas.f hpsiedmc.f 
  nonloc_grid.f multideterminant_tmove.f optwf_dummy.f acues1_reduce.f)

#VMC shared
set(VMC_PATH ${CMAKE_SOURCE_DIR}/vmc)

set(VMC90_DEPENDENCIES ${VMC_PATH}/convcartzmat.f90)

set(VMC_DEPENDENCIES
  # OBJCOM1
  ${VMC_PATH}/p2prog.f
  ${VMC_PATH}/read_input.f ${VMC_PATH}/orbitals.f 
  ${VMC_PATH}/hpsi.f ${VMC_PATH}/determinant.f ${VMC_PATH}/determinant_psit.f
  ${VMC_PATH}/bxmatrices.f ${VMC_PATH}/multideterminant.f
  ${VMC_PATH}/jastrow.f ${VMC_PATH}/jastrow_num.f ${VMC_PATH}/jastrow4.f
  ${VMC_PATH}/deriv_jastrow4.f
  ${VMC_PATH}/cuspinit4.f ${VMC_PATH}/cuspexact4.f
  ${VMC_PATH}/psi.f ${VMC_PATH}/gammai.f
  ${VMC_PATH}/matinv.f ${VMC_PATH}/distances.f ${VMC_PATH}/scale_dist.f
  ${VMC_PATH}/basis_norm.f ${VMC_PATH}/basis_fns_vgl.f
  ${VMC_PATH}/basis_fnse_vgl.f ${VMC_PATH}/basis_fnse_v.f ${VMC_PATH}/basis_fnse_vg.f
  ${VMC_PATH}/slm.f ${VMC_PATH}/read_bas_num.f ${VMC_PATH}/write_orb_loc.f
  ${VMC_PATH}/3dgrid.f ${VMC_PATH}/3dgrid_orbitals.f ${VMC_PATH}/verify_orbitals.f
  ${VMC_PATH}/readps.f ${VMC_PATH}/readps_tm.f ${VMC_PATH}/readps_gauss.f ${VMC_PATH}/readps_champ.f
  ${VMC_PATH}/pot_local.f ${VMC_PATH}/pot.f ${VMC_PATH}/nonloc_pot.f ${VMC_PATH}/nonloc.f
  ${VMC_PATH}/nonlpsi.f ${VMC_PATH}/rotqua.f
  ${VMC_PATH}/deriv_nonloc.f ${VMC_PATH}/deriv_nonlpsi.f
  ${VMC_PATH}/properties.f ${VMC_PATH}/pcm.f ${VMC_PATH}/mmpol.f ${VMC_PATH}/pcm_3dgrid.f ${VMC_PATH}/efield.f 
  ${VMC_PATH}/spline.f ${VMC_PATH}/spline2.f ${VMC_PATH}/splfit.f
  ${VMC_PATH}/rannyu.f ${VMC_PATH}/sites.f ${VMC_PATH}/nodes_distance.f
  ${VMC_PATH}/misc_grdnts.f
  ${VMC_PATH}/force_analytic.f ${VMC_PATH}/multiply_slmi_mderiv.f

  # OBJCOM1
  ${VMC_PATH}/hpsie.f ${VMC_PATH}/multideterminante.f ${VMC_PATH}/determinante.f ${VMC_PATH}/determinante_psit.f
  ${VMC_PATH}/jastrowe.f ${VMC_PATH}/jastrow4e.f ${VMC_PATH}/jassav.f ${VMC_PATH}/detsav.f 

  # OBJPW
  ${VMC_PATH}/pw_placeholder.f

  # OBJQMMM
  ${VMC_PATH}/qmmm_placeholder.f

  # OBJVMC
  ${VMC_PATH}/error.f ${VMC_PATH}/strech.f 
  ${VMC_PATH}/optwf.f ${VMC_PATH}/optwf_lib.f ${VMC_PATH}/optjas.f ${VMC_PATH}/optorb.f ${VMC_PATH}/optci.f 
  ${VMC_PATH}/optx_jas_orb.f ${VMC_PATH}/optx_jas_ci.f ${VMC_PATH}/optx_orb_ci.f
  )

# Compile and link executable
add_executable(dmc.mov1
  ${INPUT_DIR}/p2_defaults.f ${INPUT_DIR}/p2etc.f ${INPUT_DIR}/p2nmcheck.f
  ${DMC_SOURCE} ${VMC90_DEPENDENCIES} ${VMC_DEPENDENCIES})

#  Add headers
target_include_directories(dmc.mov1
  PRIVATE ${CMAKE_SOURCE_DIR}/include
  PRIVATE ${CMAKE_SOURCE_DIR}/input)

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

