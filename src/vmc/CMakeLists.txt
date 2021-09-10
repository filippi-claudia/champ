if (MPI_Fortran_FOUND)

  # Source files
  file(GLOB VMC_SOURCE *.f)
  file(GLOB VMC90_SOURCE *.f90)

  set_source_files_properties(VMC90_SOURCE PROPERTIES COMPILE_FLAGS "-check bounds -g -traceback")

  # --------------------------------------------------------------------
  # Set files according to CMAKE options:
  # --------------------------------------------------------------------

  if (NOT ENABLE_QMMM)
    message(STATUS "QMMM is NOT enabled!")
    list(REMOVE_ITEM VMC90_SOURCE
      ${CMAKE_CURRENT_SOURCE_DIR}/qmmm_bspline.f90
      ${CMAKE_CURRENT_SOURCE_DIR}/qmmm_modules.f90
      ${CMAKE_CURRENT_SOURCE_DIR}/qmmm_pot.f90
      ${CMAKE_CURRENT_SOURCE_DIR}/qmmm_density.f90
      ${CMAKE_CURRENT_SOURCE_DIR}/qmmm_writecube.f90)
  else()
    message(STATUS "QMMM is enabled!")
    list(REMOVE_ITEM VMC_SOURCE
      ${CMAKE_CURRENT_SOURCE_DIR}/qmmm_placeholder.f)
  endif()

  if (NOT ENABLE_PERIODIC)
    message(STATUS "PERIODIC is NOT enabled!")
    list(REMOVE_ITEM VMC_SOURCE
      ${CMAKE_CURRENT_SOURCE_DIR}/pw_find_image.f
      ${CMAKE_CURRENT_SOURCE_DIR}/pw_ewald.f
      ${CMAKE_CURRENT_SOURCE_DIR}/pw_orbitals.f
      ${CMAKE_CURRENT_SOURCE_DIR}/pw_orbitals_e.f
      ${CMAKE_CURRENT_SOURCE_DIR}/pw_read.f)
  else()
    message(STATUS "PERIODIC is enabled!")
    list(REMOVE_ITEM VMC_SOURCE
      ${CMAKE_CURRENT_SOURCE_DIR}/pw_placeholder.f)
  endif()

  include(${CMAKE_MODULE_PATH}/FortranFormat.cmake)
  formatFortran("${VMC_SOURCE}") # Set Fixed format property

  # --------------------------------------------------------------------
  # Create two libraries to share objects between the two binaries,
  # vmc.mov1 and dmc.mov1:
  # - shared_modules: all .f90 modules that must be compiled first.
  # - shared_objects: rest of shared objects (.f90/.f).
  # --------------------------------------------------------------------

  set(SHARED_MODULES
    ${CMAKE_CURRENT_SOURCE_DIR}/m_kinds.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/m_basis.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/m_common.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/m_periodic_table.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/m_atom.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/m_timer.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/m_det.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/m_coefs.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/m_constant.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/m_control.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/m_force.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/m_mstates.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/m_mpi.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/m_optci.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/m_pseudo.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/m_sampling.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/m_sr.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/m_dmc.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/m_vmc.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/m_grid.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/m_pcm.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/m_ewald.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/m_jastrow.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/m_mmpol.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/m_optorb.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/m_prop.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/m_efield.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/m_gradhess.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/m_optwf.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/m_estimators.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/m_deriv.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/m_mixderiv.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/m_input.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/m_grdnt.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/m_allocation.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/m_state_avrg.f90
    ${CMAKE_PARSER_DIR}/fdf.F90
    ${CMAKE_PARSER_DIR}/io_fdf.F90
    ${CMAKE_PARSER_DIR}/parse.F90
    ${CMAKE_PARSER_DIR}/prec.F90
    ${CMAKE_PARSER_DIR}/utils.F90)

  set(SHARED_OBJECTS
    ${CMAKE_CURRENT_SOURCE_DIR}/optwf_sr.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/optgeo_convcartzmat.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/optgeo_coords_int.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/optgeo_hessian.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/array_resize_utils.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/new_parser_read_data.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/set_input_data.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/pot_local.f
    ${CMAKE_CURRENT_SOURCE_DIR}/pcm_3dgrid.f
    ${CMAKE_CURRENT_SOURCE_DIR}/strech.f
    ${CMAKE_CURRENT_SOURCE_DIR}/slm.f
    ${CMAKE_CURRENT_SOURCE_DIR}/jastrow.f
    ${CMAKE_CURRENT_SOURCE_DIR}/jassav.f
    ${CMAKE_CURRENT_SOURCE_DIR}/gammai.f
    ${CMAKE_CURRENT_SOURCE_DIR}/sites.f
    ${CMAKE_CURRENT_SOURCE_DIR}/deriv_nonloc.f
    ${CMAKE_CURRENT_SOURCE_DIR}/splfit.f
    ${CMAKE_CURRENT_SOURCE_DIR}/readps.f
    ${CMAKE_CURRENT_SOURCE_DIR}/nonloc.f
    ${CMAKE_CURRENT_SOURCE_DIR}/deriv_jastrow4.f
    ${CMAKE_CURRENT_SOURCE_DIR}/qmmm_placeholder.f    # This file should be refactored and eliminated.
    ${CMAKE_CURRENT_SOURCE_DIR}/efield.f
    ${CMAKE_CURRENT_SOURCE_DIR}/nodes_distance.f
    ${CMAKE_CURRENT_SOURCE_DIR}/hpsie.f
    ${CMAKE_CURRENT_SOURCE_DIR}/verify_orbitals.f
    ${CMAKE_CURRENT_SOURCE_DIR}/properties.f
    ${CMAKE_CURRENT_SOURCE_DIR}/multideterminante.f
    ${CMAKE_CURRENT_SOURCE_DIR}/multideterminant.f
    ${CMAKE_CURRENT_SOURCE_DIR}/pcm.f
    ${CMAKE_CURRENT_SOURCE_DIR}/jastrowe.f
    ${CMAKE_CURRENT_SOURCE_DIR}/write_orb_loc.f
    ${CMAKE_CURRENT_SOURCE_DIR}/cuspexact4.f
    ${CMAKE_CURRENT_SOURCE_DIR}/optwf_matrix_corsamp.f
    ${CMAKE_CURRENT_SOURCE_DIR}/optwf_lin_matrix.f
    ${CMAKE_CURRENT_SOURCE_DIR}/optwf_lib.f
    ${CMAKE_CURRENT_SOURCE_DIR}/matinv.f
    ${CMAKE_CURRENT_SOURCE_DIR}/scale_dist.f
    ${CMAKE_CURRENT_SOURCE_DIR}/orbitals.f
    ${CMAKE_CURRENT_SOURCE_DIR}/readps_gauss.f
    ${CMAKE_CURRENT_SOURCE_DIR}/determinante_psit.f
    ${CMAKE_CURRENT_SOURCE_DIR}/readps_tm.f
    ${CMAKE_CURRENT_SOURCE_DIR}/pot.f
    ${CMAKE_CURRENT_SOURCE_DIR}/cuspinit4.f
    ${CMAKE_CURRENT_SOURCE_DIR}/mmpol.f
    ${CMAKE_CURRENT_SOURCE_DIR}/basis_fnse_vgl.f
    ${CMAKE_CURRENT_SOURCE_DIR}/determinant_psit.f
    ${CMAKE_CURRENT_SOURCE_DIR}/determinant.f
    ${CMAKE_CURRENT_SOURCE_DIR}/determinante.f
    ${CMAKE_CURRENT_SOURCE_DIR}/nonloc_pot.f
    ${CMAKE_CURRENT_SOURCE_DIR}/pw_placeholder.f    # This file should be refactored and eliminated.
    ${CMAKE_CURRENT_SOURCE_DIR}/basis_fns_vgl.f
    ${CMAKE_CURRENT_SOURCE_DIR}/jastrow_num.f
    ${CMAKE_CURRENT_SOURCE_DIR}/bxmatrices.f
    ${CMAKE_CURRENT_SOURCE_DIR}/deriv_nonlpsi.f
    ${CMAKE_CURRENT_SOURCE_DIR}/force_analytic.f
    ${CMAKE_CURRENT_SOURCE_DIR}/nonlpsi.f
    ${CMAKE_CURRENT_SOURCE_DIR}/basis_norm.f
    ${CMAKE_CURRENT_SOURCE_DIR}/psi.f
    ${CMAKE_CURRENT_SOURCE_DIR}/rotqua.f
    ${CMAKE_CURRENT_SOURCE_DIR}/misc_grdnts.f
    ${CMAKE_CURRENT_SOURCE_DIR}/spline.f
    ${CMAKE_CURRENT_SOURCE_DIR}/spline2.f
    ${CMAKE_CURRENT_SOURCE_DIR}/read_bas_num.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/store_diag_hs.f
    ${CMAKE_CURRENT_SOURCE_DIR}/rannyu.f
    ${CMAKE_CURRENT_SOURCE_DIR}/3dgrid_orbitals.f
    ${CMAKE_CURRENT_SOURCE_DIR}/readps_champ.f
    ${CMAKE_CURRENT_SOURCE_DIR}/distances.f
    ${CMAKE_CURRENT_SOURCE_DIR}/multiply_slmi_mderiv.f
    ${CMAKE_CURRENT_SOURCE_DIR}/basis_fnse_vg.f
    ${CMAKE_CURRENT_SOURCE_DIR}/parser.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/jastrow4.f
    ${CMAKE_CURRENT_SOURCE_DIR}/hpsi.f
    ${CMAKE_CURRENT_SOURCE_DIR}/jastrow4e.f
    ${CMAKE_CURRENT_SOURCE_DIR}/3dgrid.f
    ${CMAKE_CURRENT_SOURCE_DIR}/detsav.f
    ${CMAKE_CURRENT_SOURCE_DIR}/basis_fnse_v.f
    ${CMAKE_CURRENT_SOURCE_DIR}/optgeo_lib.f
    ${CMAKE_CURRENT_SOURCE_DIR}/optjas.f
    ${CMAKE_CURRENT_SOURCE_DIR}/optorb.f
    ${CMAKE_CURRENT_SOURCE_DIR}/optci.f
    ${CMAKE_CURRENT_SOURCE_DIR}/optx_orb_ci.f
    ${CMAKE_CURRENT_SOURCE_DIR}/optx_jas_ci.f
    ${CMAKE_CURRENT_SOURCE_DIR}/optx_jas_orb.f
    ${CMAKE_CURRENT_SOURCE_DIR}/optwf_handle_wf.f
    ${CMAKE_CURRENT_SOURCE_DIR}/optwf_sr_more.f
    ${CMAKE_CURRENT_SOURCE_DIR}/optx_jas_ci_reduce.f
    ${CMAKE_CURRENT_SOURCE_DIR}/optx_jas_orb_reduce.f
    ${CMAKE_CURRENT_SOURCE_DIR}/optx_orb_ci_reduce.f
    ${CMAKE_CURRENT_SOURCE_DIR}/optci_reduce.f
    ${CMAKE_CURRENT_SOURCE_DIR}/optjas_reduce.f
    ${CMAKE_CURRENT_SOURCE_DIR}/optorb_reduce.f
    ${CMAKE_CURRENT_SOURCE_DIR}/finwrt_more.f
    ${CMAKE_CURRENT_SOURCE_DIR}/error.f)


  # --------------------------------------------------------------------
  # Remove objects to share depending on CMAKE options:
  # --------------------------------------------------------------------
  if (ENABLE_QMMM)
    list(REMOVE_ITEM SHARED_OBJECTS
      ${CMAKE_CURRENT_SOURCE_DIR}/qmmm_placeholder.f)
    list(APPEND SHARED_OBJECTS
      ${CMAKE_CURRENT_SOURCE_DIR}/qmmm_bspline.f90
      ${CMAKE_CURRENT_SOURCE_DIR}/qmmm_modules.f90
      ${CMAKE_CURRENT_SOURCE_DIR}/qmmm_pot.f90
      ${CMAKE_CURRENT_SOURCE_DIR}/qmmm_density.f90
      ${CMAKE_CURRENT_SOURCE_DIR}/qmmm_writecube.f90)
  endif()

  if (ENABLE_PERIODIC)
    list(REMOVE_ITEM SHARED_SOURCE
      ${CMAKE_CURRENT_SOURCE_DIR}/pw_placeholder.f)
    list(APPEND SHARED_OBJECTS
      ${CMAKE_CURRENT_SOURCE_DIR}/pw_ewald.f
      ${CMAKE_CURRENT_SOURCE_DIR}/pw_orbitals_e.f
      ${CMAKE_CURRENT_SOURCE_DIR}/pw_find_image.f
      ${CMAKE_CURRENT_SOURCE_DIR}/pw_orbitals.f
      ${CMAKE_CURRENT_SOURCE_DIR}/pw_read.f)
  endif()

  add_library(shared_modules OBJECT ${SHARED_MODULES})
  target_compile_options(shared_modules
    PRIVATE
    ${Fortran_FLAGS}
    "$<$<CONFIG:Debug>:${Fortan_FLAGS_DEBUG}>")

  add_library(shared_objects OBJECT ${SHARED_OBJECTS})
  target_include_directories(shared_objects
    PUBLIC
    ${CMAKE_SOURCE_DIR}/src/input)

  target_compile_options(shared_objects
    PRIVATE
    ${Fortran_FLAGS}
    "$<$<CONFIG:Debug>:${Fortan_FLAGS_DEBUG}>")

  #add_dependencies(shared_objects generate-p2 shared_modules)
  add_dependencies(shared_objects shared_modules)

  # --------------------------------------------------------------------
  # Remove code that has been already compiled:
  # --------------------------------------------------------------------

  foreach(s ${SHARED_OBJECTS})
    list(REMOVE_ITEM VMC_SOURCE ${s})
  endforeach(s)
  foreach(s ${SHARED_OBJECTS})
    list(REMOVE_ITEM VMC90_SOURCE ${s})
  endforeach(s)
  foreach(s ${SHARED_MODULES})
    list(REMOVE_ITEM VMC_SOURCE ${s})
  endforeach(s)
  foreach(s ${SHARED_MODULES})
    list(REMOVE_ITEM VMC90_SOURCE ${s})
  endforeach(s)

  # --------------------------------------------------------------------
  # Compile and link executable:
  # --------------------------------------------------------------------

  add_executable(vmc.mov1
    $<TARGET_OBJECTS:shared_modules>
    $<TARGET_OBJECTS:shared_objects>
    # source that is not already included in the objects
    ${VMC_SOURCE}
    ${VMC90_SOURCE})

  # compilation
  target_compile_options(vmc.mov1
    PRIVATE
    ${Fortran_FLAGS}
    "$<$<CONFIG:Debug>:${Fortan_FLAGS_DEBUG}>"
    )

  target_link_libraries(vmc.mov1
    cyrus pspline
    ${LINEAR_ALGEBRA})

endif()
