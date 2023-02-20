if (MPI_Fortran_FOUND)

  # Source files
  file(GLOB VMC_SOURCE
    3dgrid.f
    3dgrid_orbitals.f
    acuest.f
    acuest_reduce.f
    acuest_write.f
    basis_fns.f
    bxmatrices.f
    cuspexact4.f
    deriv_jastrow.f
    deriv_jastrow4.f
    deriv_nonloc.f
    deriv_nonlpsi.f
    determinant.f
    determinant_psig.f
    determinant_psit.f
    determinante.f
    determinante_psit.f
    detsav.f
    distances.f
    dumper.f
    dumper_more.f
    efield.f
    fetch_parameters.f
    fin_reduce.f
    finwrt.f
    finwrt_more.f
    force_analy_reduce.f
    force_analytic.f
    gammai.f
    hpsi.f
    hpsie.f
    jassav.f
    jastrow.f
    jastrow4.f
    jastrow4e.f
    jastrow_num.f
    jastrowe.f
    mc_configs.f
    metrop_mov1_slat.f
    misc_grdnts.f
    mmpol.f
    mmpol_reduce.f
    mmpol_vmc.f
    multideterminant.f
    multideterminante.f
    multiple_states.f
    multiply_slmi_mderiv.f
    nodes_distance.f
    nonloc.f
    nonloc_pot.f
    nonlpsi.f
    optci.f
    optci_reduce.f
    optgeo_lib.f
    optjas.f
    optjas_reduce.f
    optorb.f
    optorb_reduce.f
    optwf_handle_wf.f
    optwf_lib.f
    optwf_lin_dav.f
    optwf_lin_dav_more.f
    optwf_lin_dav_extra.f
    optwf_lin_matrix.f
    optwf_matrix_corsamp.f
    optwf_mix_sa.f
    optwf_olbfgs.f
    optwf_olbfgs_more.f
    optwf_sr_more.f
    optwf_vmc.f
    optx_jas_ci.f
    optx_jas_ci_reduce.f
    optx_jas_orb.f
    optx_jas_orb_reduce.f
    optx_orb_ci.f
    optx_orb_ci_reduce.f
    orbitals.f
    pcm.f
    pcm_3dgrid.f
    pcm_reduce.f
    pcm_vmc.f
    pot.f
    pot_local.f
    prop_reduce.f
    prop_vmc.f
    properties.f
    psi.f
    pw_ewald.f
    pw_find_image.f
    pw_orbitals.f
    pw_orbitals_e.f
    pw_placeholder.f
    pw_read.f
    qmmm_placeholder.f
    random.f
    readps_gauss.f
    readps_tm.f
    rotqua.f
    scale_dist.f
    sites.f
    splfit.f
    spline2.f
    store_diag_hs.f
    strech.f
    verify_orbitals.f
    vmc.f
	  )

  file(GLOB VMC90_SOURCE
    array_resize_utils.f90
    array_utils.f90
    cuspinit4.f90
    davidson.f90
    davidson_wrap.f90
    get_norbterm.f90
    lapack_wrapper.f90
    main.f90
    matinv.f90
    new_parser_read_data.f90
    read_trexio_data.f90
    fitting.f90
    olbfgs.f90
    optgeo_convcartzmat.f90
    optgeo_coords_int.f90
    optgeo_hessian.f90
    optwf.f90 optwf_dl.f90
    optwf_sr.f90
    parser.f90
    qmmm_bspline.f90
    qmmm_density.f90
    qmmm_modules.f90
    qmmm_pot.f90
    qmmm_writecube.f90
    read_bas_num.f90
    regterg.f90
    set_input_data.f90
    slm.f90
    write_orb_loc.f90
    )

  set_property(
  SOURCE new_parser_read_data.f90
  APPEND
  PROPERTY COMPILE_DEFINITIONS
  CMAKE_Fortran_COMPILER=\"${CMAKE_Fortran_COMPILER}\"
  CMAKE_Fortran_COMPILER_VERSION=\"${CMAKE_Fortran_COMPILER_VERSION}\"
  GIT_HEAD_BRANCH=\"${GIT_HEAD_BRANCH}\"
  GIT_REVISION_HASH=\"${GIT_REVISION_HASH}\"
  HDF5_VERSION=\"${HDF5_VERSION}\"
  )


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

  find_package(FortranFormat)
  formatFortran("${VMC_SOURCE}") # Set Fixed format property

  # --------------------------------------------------------------------
  # Create two libraries to share objects between the two binaries,
  # vmc.mov1 and dmc.mov1:
  # - shared_modules: all .f90 modules that must be compiled first.
  # - shared_objects: rest of shared objects (.f90/.f).
  # --------------------------------------------------------------------

  set(SHARED_OBJECTS
    ${CMAKE_CURRENT_SOURCE_DIR}/array_resize_utils.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/lapack_wrapper.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/array_utils.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/optwf_sr.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/optgeo_convcartzmat.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/optgeo_coords_int.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/optgeo_hessian.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/new_parser_read_data.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/read_trexio_data.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/set_input_data.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/pot_local.f
    ${CMAKE_CURRENT_SOURCE_DIR}/pcm_3dgrid.f
    ${CMAKE_CURRENT_SOURCE_DIR}/strech.f
    ${CMAKE_CURRENT_SOURCE_DIR}/slm.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/jastrow.f
    ${CMAKE_CURRENT_SOURCE_DIR}/jassav.f
    ${CMAKE_CURRENT_SOURCE_DIR}/gammai.f
    ${CMAKE_CURRENT_SOURCE_DIR}/sites.f
    ${CMAKE_CURRENT_SOURCE_DIR}/deriv_nonloc.f
    ${CMAKE_CURRENT_SOURCE_DIR}/splfit.f
    ${CMAKE_CURRENT_SOURCE_DIR}/fitting.f90
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
    ${CMAKE_CURRENT_SOURCE_DIR}/write_orb_loc.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/cuspexact4.f
    ${CMAKE_CURRENT_SOURCE_DIR}/optwf_matrix_corsamp.f
    ${CMAKE_CURRENT_SOURCE_DIR}/optwf_lin_matrix.f
    ${CMAKE_CURRENT_SOURCE_DIR}/optwf_lib.f
    ${CMAKE_CURRENT_SOURCE_DIR}/matinv.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/scale_dist.f
    ${CMAKE_CURRENT_SOURCE_DIR}/orbitals.f
    ${CMAKE_CURRENT_SOURCE_DIR}/readps_gauss.f
    ${CMAKE_CURRENT_SOURCE_DIR}/readps_tm.f
    ${CMAKE_CURRENT_SOURCE_DIR}/determinante_psit.f
    ${CMAKE_CURRENT_SOURCE_DIR}/pot.f
    ${CMAKE_CURRENT_SOURCE_DIR}/cuspinit4.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/mmpol.f
    ${CMAKE_CURRENT_SOURCE_DIR}/basis_fns.f
    ${CMAKE_CURRENT_SOURCE_DIR}/determinant_psit.f
    ${CMAKE_CURRENT_SOURCE_DIR}/determinant.f
    ${CMAKE_CURRENT_SOURCE_DIR}/determinante.f
    ${CMAKE_CURRENT_SOURCE_DIR}/nonloc_pot.f
    ${CMAKE_CURRENT_SOURCE_DIR}/pw_placeholder.f    # This file should be refactored and eliminated.
    ${CMAKE_CURRENT_SOURCE_DIR}/jastrow_num.f
    ${CMAKE_CURRENT_SOURCE_DIR}/bxmatrices.f
    ${CMAKE_CURRENT_SOURCE_DIR}/deriv_nonlpsi.f
    ${CMAKE_CURRENT_SOURCE_DIR}/force_analytic.f
    ${CMAKE_CURRENT_SOURCE_DIR}/nonlpsi.f
    ${CMAKE_CURRENT_SOURCE_DIR}/psi.f
    ${CMAKE_CURRENT_SOURCE_DIR}/rotqua.f
    ${CMAKE_CURRENT_SOURCE_DIR}/misc_grdnts.f
    ${CMAKE_CURRENT_SOURCE_DIR}/spline2.f
    ${CMAKE_CURRENT_SOURCE_DIR}/read_bas_num.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/store_diag_hs.f
    ${CMAKE_CURRENT_SOURCE_DIR}/random.f
    ${CMAKE_CURRENT_SOURCE_DIR}/3dgrid_orbitals.f
    ${CMAKE_CURRENT_SOURCE_DIR}/distances.f
    ${CMAKE_CURRENT_SOURCE_DIR}/multiply_slmi_mderiv.f
    ${CMAKE_CURRENT_SOURCE_DIR}/parser.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/jastrow4.f
    ${CMAKE_CURRENT_SOURCE_DIR}/hpsi.f
    ${CMAKE_CURRENT_SOURCE_DIR}/jastrow4e.f
    ${CMAKE_CURRENT_SOURCE_DIR}/3dgrid.f
    ${CMAKE_CURRENT_SOURCE_DIR}/detsav.f
    ${CMAKE_CURRENT_SOURCE_DIR}/optgeo_lib.f
    ${CMAKE_CURRENT_SOURCE_DIR}/optjas.f
    ${CMAKE_CURRENT_SOURCE_DIR}/optorb.f
    ${CMAKE_CURRENT_SOURCE_DIR}/get_norbterm.f90
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
    ${CMAKE_CURRENT_SOURCE_DIR}/finwrt_more.f)

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



  add_library(shared_objects OBJECT ${SHARED_OBJECTS})

  #  Add headers
  target_include_directories(shared_objects
  PUBLIC
  ${CMAKE_BINARY_DIR}/src/module
  )

  target_compile_options(shared_objects
    PRIVATE
    ${Fortran_FLAGS}
    "$<$<CONFIG:Debug>:${Fortan_FLAGS_DEBUG}>")

  # add_dependencies(shared_objects parser)
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
  add_library(rng OBJECT
    xoroshiro256starstar.c)
  
  target_compile_options(rng
    PUBLIC
    "-O3")


  add_executable(vmc.mov1
    $<TARGET_OBJECTS:shared_modules>
    $<TARGET_OBJECTS:shared_objects>
    $<TARGET_OBJECTS:rng>
    # source that is not already included in the objects
    ${VMC_SOURCE}
    ${VMC90_SOURCE}
    )


IF(UNIT_TESTS)
  add_library(vmc SHARED
    $<TARGET_OBJECTS:shared_modules>
    $<TARGET_OBJECTS:shared_objects>
    $<TARGET_OBJECTS:rng>
    $<TARGET_OBJECTS:vmc.mov1>)

  target_compile_options(vmc
    PRIVATE
    ${Fortran_FLAGS}
    "$<$<CONFIG:Debug>:${Fortan_FLAGS_DEBUG}>"
    )

ENDIF()

  #  Add headers
  target_include_directories(vmc.mov1
  PUBLIC
  ${CMAKE_BINARY_DIR}/src/module
  )

  # compilation
  target_compile_options(vmc.mov1
    PRIVATE
    ${Fortran_FLAGS}
    "$<$<CONFIG:Debug>:${Fortan_FLAGS_DEBUG}>"
    )

  if (TREXIO_FOUND)
    target_link_libraries(vmc.mov1
      PRIVATE
      pspline
      parser
      ${LINEAR_ALGEBRA}
      ${TREXIO_LIBRARY}
      MPI::MPI_Fortran)
  else()
    target_link_libraries(vmc.mov1
      PRIVATE
      pspline
      parser
      ${LINEAR_ALGEBRA}
      MPI::MPI_Fortran)
  endif()
  

endif()
