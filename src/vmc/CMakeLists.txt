if (MPI_Fortran_FOUND)

  # Source files
  file(GLOB VMC90_SOURCE
    3dgrid.f90
    3dgrid_orbitals.f90
    acuest.f90
    acuest_reduce.f90
    acuest_write.f90
    basis_fns.f90
    bxmatrices.f90
    cuspexact4.f90
    deriv_nonloc.f90
    deriv_nonlpsi.f90
    determinant.f90
    determinant_psig.f90
    determinant_psit.f90
    determinante.f90
    determinante_psit.f90
    detsav.f90
    distances.f90
    dumper.f90
    dumper_more.f90
    efield.f90
    fetch_parameters.f90
    fin_reduce.f90
    finwrt.f90
    finwrt_more.f90
    force_analy_reduce.f90
    force_analytic.f90
    gammai.f90
    gauss.f90
    hpsi.f90
    hpsie.f90
    jassav.f90
    jastrow.f90
    jastrow1.f90
    jastrow4.f90
    jastrow4e.f90
    jastrow1e.f90
    jastrowe.f90
    jastrow_qmckl.f90
    mc_configs.f90
    metrop_mov1_driftdif.f90
    metrop_mov1_slat.f90
    metrop_driftdif.f90
    metrop_slat.f90
    misc_grdnts.f90
    mmpol.f90
    mmpol_reduce.f90
    mmpol_vmc.f90
    multideterminant.f90
    multideterminante.f90
    multiple_states.f90
    multiply_slmi_mderiv.f90
    nodes_distance.f90
    nonloc.f90
    nonloc_pot.f90
    nonlpsi.f90
    optci.f90
    optci_reduce.f90
    optgeo_lib.f90
    optjas.f90
    optjas_reduce.f90
    optorb.f90
    optorb_reduce.f90
    optwf_handle_wf.f90
    optwf_lib.f90
    optwf_lin_dav.f90
    optwf_lin_dav_more.f90
    optwf_lin_dav_extra.f90
    optwf_lin_matrix.f90
    optwf_matrix_corsamp.f90
    optwf_mix_sa.f90
    optwf_olbfgs.f90
    optwf_olbfgs_more.f90
    optwf_sr_more.f90
    optwf_vmc.f90
    optx_jas_ci.f90
    optx_jas_ci_reduce.f90
    optx_jas_orb.f90
    optx_jas_orb_reduce.f90
    optx_orb_ci.f90
    optx_orb_ci_reduce.f90
    orbitals_no_qmckl.f90
    orbitals_qmckl_periodic.f90
    orbitals_qmckl.f90
    orbitals.f90
    pcm.f90
    pcm_3dgrid.f90
    pcm_reduce.f90
    pcm_vmc.f90
    pot.f90
    pot_local.f90
    prop_reduce.f90
    prop_vmc.f90
    properties.f90
    ewald.f90
    find_image.f90
    qmmm_placeholder.f90
    random.f90
    readps_gauss.f90
    rotqua.f90
    scale_dist.f90
    sites.f90
    splfit.f90
    spline2.f90
    store_diag_hs.f90
    strech.f90
    verify_orbitals.f90
    vmc.f90
    array_resize_utils.f90
    array_utils.f90
    cuspinit4.f90
    davidson.f90
    davidson_wrap.f90
    deriv_jastrow_qmckl.f90
    deriv_jastrow1.f90
    deriv_jastrow4.f90
    get_norbterm.f90
    lapack_wrapper.f90
    main.f90
    matinv.f90
    read_data_parser.f90
    read_data_trexio.f90
    fitting.f90
    olbfgs.f90
    optgeo_convcartzmat.f90
    optgeo_coords_int.f90
    optgeo_hessian.f90
    optwf.f90 optwf_dl.f90
    optwf_sr.f90
    ortho_orbitals.f90
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
  SOURCE read_data_parser.f90
  APPEND
  PROPERTY COMPILE_DEFINITIONS
  CMAKE_Fortran_COMPILER=\"${CMAKE_Fortran_COMPILER}\"
  CMAKE_Fortran_COMPILER_VERSION=\"${CMAKE_Fortran_COMPILER_VERSION}\"
  GIT_HEAD_BRANCH=\"${GIT_HEAD_BRANCH}\"
  GIT_REVISION_HASH=\"${GIT_REVISION_HASH}\"
  HDF5_VERSION=\"${HDF5_VERSION}\"
  )


  set(ENABLE_GPU ${ENABLE_GPU} PARENT_SCOPE)
  message (STATUS "GPU IS ENABLED")

  set(ENABLE_TREXIO ${ENABLE_TREXIO} PARENT_SCOPE)
  set(ENABLE_QMCKL  ${ENABLE_QMCKL} PARENT_SCOPE)

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
    list(REMOVE_ITEM VMC90_SOURCE
      ${CMAKE_CURRENT_SOURCE_DIR}/qmmm_placeholder.f90)
  endif()

  if (NOT QMCKL_FOUND)
	  list(REMOVE_ITEM VMC90_SOURCE
         ${CMAKE_CURRENT_SOURCE_DIR}/orbitals_qmckl.f90
         ${CMAKE_CURRENT_SOURCE_DIR}/orbitals_qmckl_periodic.f90
         ${CMAKE_CURRENT_SOURCE_DIR}/jastrow_qmckl.f90
	 ${CMAKE_CURRENT_SOURCE_DIR}/deriv_jastrow_qmckl.f90)
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
    ${CMAKE_CURRENT_SOURCE_DIR}/read_data_parser.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/read_data_trexio.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/set_input_data.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/pot_local.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/pcm_3dgrid.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/strech.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/slm.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/jastrow.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/jassav.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/gammai.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/gauss.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/sites.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/deriv_nonloc.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/splfit.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/fitting.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/nonloc.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/deriv_jastrow1.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/deriv_jastrow4.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/qmmm_placeholder.f90   # This file should be refactored and eliminated.
    ${CMAKE_CURRENT_SOURCE_DIR}/efield.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/nodes_distance.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/hpsie.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/verify_orbitals.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/properties.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/multideterminante.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/multideterminant.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/pcm.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/jastrowe.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/write_orb_loc.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/cuspexact4.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/optwf_matrix_corsamp.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/optwf_lin_matrix.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/optwf_lib.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/ortho_orbitals.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/matinv.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/scale_dist.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/orbitals_no_qmckl.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/orbitals.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/readps_gauss.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/determinante_psit.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/pot.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/cuspinit4.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/mmpol.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/basis_fns.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/determinant_psit.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/determinant.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/determinante.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/nonloc_pot.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/bxmatrices.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/deriv_nonlpsi.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/force_analytic.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/force_analy_reduce.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/nonlpsi.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/rotqua.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/misc_grdnts.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/spline2.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/read_bas_num.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/store_diag_hs.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/random.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/3dgrid_orbitals.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/find_image.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/ewald.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/distances.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/multiply_slmi_mderiv.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/parser.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/jastrow1.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/jastrow1e.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/jastrow4.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/hpsi.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/jastrow4e.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/3dgrid.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/detsav.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/optgeo_lib.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/optjas.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/optorb.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/get_norbterm.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/optci.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/optx_orb_ci.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/optx_jas_ci.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/optx_jas_orb.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/optwf_handle_wf.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/optwf_sr_more.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/optx_jas_ci_reduce.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/optx_jas_orb_reduce.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/optx_orb_ci_reduce.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/optci_reduce.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/optjas_reduce.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/optorb_reduce.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/finwrt_more.f90)

  # --------------------------------------------------------------------
  # Remove objects to share depending on CMAKE options:
  # --------------------------------------------------------------------
  if (QMCKL_FOUND)
    list(APPEND SHARED_OBJECTS
      ${CMAKE_CURRENT_SOURCE_DIR}/jastrow_qmckl.f90
      ${CMAKE_CURRENT_SOURCE_DIR}/deriv_jastrow_qmckl.f90
      ${CMAKE_CURRENT_SOURCE_DIR}/orbitals_qmckl.f90
      ${CMAKE_CURRENT_SOURCE_DIR}/orbitals_qmckl_periodic.f90
	 )
  endif()

  if (ENABLE_QMMM)
    list(REMOVE_ITEM SHARED_OBJECTS
      ${CMAKE_CURRENT_SOURCE_DIR}/qmmm_placeholder.f90)
    list(APPEND SHARED_OBJECTS
      ${CMAKE_CURRENT_SOURCE_DIR}/qmmm_bspline.f90
      ${CMAKE_CURRENT_SOURCE_DIR}/qmmm_modules.f90
      ${CMAKE_CURRENT_SOURCE_DIR}/qmmm_pot.f90
      ${CMAKE_CURRENT_SOURCE_DIR}/qmmm_density.f90
      ${CMAKE_CURRENT_SOURCE_DIR}/qmmm_writecube.f90)
  endif()


  add_library(shared_objects OBJECT ${SHARED_OBJECTS})

  #  Add headers
  target_include_directories(shared_objects
  PUBLIC
  ${CMAKE_BINARY_DIR}/src/module
  ${CMAKE_BINARY_DIR}/src/parser
  )

  target_compile_options(shared_objects
    PRIVATE
    ${Fortran_FLAGS}
    "$<$<CONFIG:Debug>:${Fortan_FLAGS_DEBUG}>")

  # add_dependencies(shared_objects parser)
  add_dependencies(shared_objects shared_modules)
  add_dependencies(shared_objects parser)

  # --------------------------------------------------------------------
  # Remove code that has been already compiled:
  # --------------------------------------------------------------------

  foreach(s ${SHARED_OBJECTS})
    list(REMOVE_ITEM VMC90_SOURCE ${s})
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

  if(CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC" OR CMAKE_Fortran_COMPILER_ID MATCHES "Fujitsu")
     target_link_options(vmc.mov1
       PRIVATE
       ${Fortran_FLAGS}
       "$<$<CONFIG:Debug>:${Fortan_FLAGS_DEBUG}>"
       )
  endif()

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
      ${CUDA_LIBRARIES}
      MPI::MPI_Fortran)
  endif()


endif()
