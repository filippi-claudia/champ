if (MPI_Fortran_FOUND)

  # Source files
  file(GLOB VMC_SOURCE LIST_DIRECTORIES false  *.f)
  file(GLOB VMC90_SOURCE LIST_DIRECTORIES false optgeo_coords_int.f90 optgeo_hessian.f90 optgeo_convcartzmat.f90 regterg.f90)
  
  if (NOT ENABLE_QMMM)
    message(STATUS "QMMM is NOT enabled!")
  endif()
  
  if (NOT ENABLE_PERIODIC)
    message(STATUS "PERIODIC is NOT enabled!")
    list(REMOVE_ITEM VMC_SOURCE
      ${CMAKE_CURRENT_SOURCE_DIR}/pw_find_image.f ${CMAKE_CURRENT_SOURCE_DIR}/pw_ewald.f
      ${CMAKE_CURRENT_SOURCE_DIR}/pw_orbitals.f ${CMAKE_CURRENT_SOURCE_DIR}/pw_orbitals_e.f
      ${CMAKE_CURRENT_SOURCE_DIR}/pw_read.f)
  endif()
  
  include(${CMAKE_MODULE_PATH}/FortranFormat.cmake)
  formatFortran("${VMC_SOURCE}") # Set Fixed format property
  
  # --------------------------------------------
  # creates object library shared for all binaries
  # --------------------------------------------
  
  set(SHARED_SOURCE_1
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
    ${CMAKE_CURRENT_SOURCE_DIR}/qmmm_placeholder.f
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
    ${CMAKE_CURRENT_SOURCE_DIR}/p2prog.f
    ${CMAKE_CURRENT_SOURCE_DIR}/basis_fnse_vgl.f
    ${CMAKE_CURRENT_SOURCE_DIR}/determinant_psit.f
    ${CMAKE_CURRENT_SOURCE_DIR}/determinant.f
    ${CMAKE_CURRENT_SOURCE_DIR}/determinante.f
    ${CMAKE_CURRENT_SOURCE_DIR}/nonloc_pot.f
    ${CMAKE_CURRENT_SOURCE_DIR}/pw_placeholder.f
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
    ${CMAKE_CURRENT_SOURCE_DIR}/read_bas_num.f
    ${CMAKE_CURRENT_SOURCE_DIR}/rannyu.f
    ${CMAKE_CURRENT_SOURCE_DIR}/3dgrid_orbitals.f
    ${CMAKE_CURRENT_SOURCE_DIR}/readps_champ.f
    ${CMAKE_CURRENT_SOURCE_DIR}/distances.f
    ${CMAKE_CURRENT_SOURCE_DIR}/multiply_slmi_mderiv.f
    ${CMAKE_CURRENT_SOURCE_DIR}/basis_fnse_vg.f
    ${CMAKE_CURRENT_SOURCE_DIR}/read_input.f
    ${CMAKE_CURRENT_SOURCE_DIR}/jastrow4.f
    ${CMAKE_CURRENT_SOURCE_DIR}/hpsi.f
    ${CMAKE_CURRENT_SOURCE_DIR}/jastrow4e.f
    ${CMAKE_CURRENT_SOURCE_DIR}/3dgrid.f
    ${CMAKE_CURRENT_SOURCE_DIR}/detsav.f
    ${CMAKE_CURRENT_SOURCE_DIR}/basis_fnse_v.f
    ${CMAKE_CURRENT_SOURCE_DIR}/optgeo_convcartzmat.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/optgeo_coords_int.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/optgeo_hessian.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/optgeo_lib.f
    ${INPUT_DIR}/p2_defaults.f
    ${INPUT_DIR}/p2etc.f
    ${INPUT_DIR}/p2nmcheck.f)
  
  add_library(shared_objects_1 OBJECT ${SHARED_SOURCE_1})
  
  target_include_directories(shared_objects_1
    PUBLIC
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/input)
  
  target_compile_options(shared_objects_1
    PRIVATE
    ${Fortran_FLAGS}
    "$<$<CONFIG:Debug>:${Fortan_FLAGS_DEBUG}>")
  
  add_dependencies(shared_objects_1 generate-p2)
  
  foreach(s ${SHARED_SOURCE_1})
    list(REMOVE_ITEM VMC_SOURCE ${s})
  endforeach(s)
  foreach(s ${SHARED_SOURCE_1})
    list(REMOVE_ITEM VMC90_SOURCE ${s})
  endforeach(s)
  
  # --------------------------------------------
  
  # --------------------------------------------
  # creates object library shared for all binaries
  # --------------------------------------------
  
  set(SHARED_SOURCE_2
    ${CMAKE_CURRENT_SOURCE_DIR}/optjas.f
    ${CMAKE_CURRENT_SOURCE_DIR}/optorb.f
    ${CMAKE_CURRENT_SOURCE_DIR}/optci.f
    ${CMAKE_CURRENT_SOURCE_DIR}/optx_orb_ci.f
    ${CMAKE_CURRENT_SOURCE_DIR}/optx_jas_ci.f
    ${CMAKE_CURRENT_SOURCE_DIR}/optx_jas_orb.f)
  
  add_library(shared_objects_2 OBJECT ${SHARED_SOURCE_2})
  
  target_include_directories(shared_objects_2
    PUBLIC
    ${CMAKE_SOURCE_DIR}/include)
  
  target_compile_options(shared_objects_2
    PRIVATE
    ${Fortran_FLAGS}
    "$<$<CONFIG:Debug>:${Fortan_FLAGS_DEBUG}>")
  
  foreach(s ${SHARED_SOURCE_2})
    list(REMOVE_ITEM VMC_SOURCE ${s})
  endforeach(s)
  
  # --------------------------------------------
  # creates object library for vmc.mov1.mpi
  # --------------------------------------------
  
  set(SHARED_VMC_MPI_SOURCE
    ${CMAKE_CURRENT_SOURCE_DIR}/fetch_parameters.f
    ${CMAKE_CURRENT_SOURCE_DIR}/acuest.f
    ${CMAKE_CURRENT_SOURCE_DIR}/metrop_mov1_slat.f
    ${CMAKE_CURRENT_SOURCE_DIR}/acuest_write.f
    ${CMAKE_CURRENT_SOURCE_DIR}/jdqz_lib.f
    ${CMAKE_CURRENT_SOURCE_DIR}/prop_vmc.f
    ${CMAKE_CURRENT_SOURCE_DIR}/optwf_lin_dav.f
    ${CMAKE_CURRENT_SOURCE_DIR}/optwf_mix_sa.f
    ${CMAKE_CURRENT_SOURCE_DIR}/optwf_sr.f
    ${CMAKE_CURRENT_SOURCE_DIR}/optwf_vmc.f
    ${CMAKE_CURRENT_SOURCE_DIR}/optwf_dl.f
    ${CMAKE_CURRENT_SOURCE_DIR}/optwf_dl_iter.f
    ${CMAKE_CURRENT_SOURCE_DIR}/vmc.f
    ${CMAKE_CURRENT_SOURCE_DIR}/pcm_vmc.f
    ${CMAKE_CURRENT_SOURCE_DIR}/finwrt.f
    ${CMAKE_CURRENT_SOURCE_DIR}/dumper_more.f
    ${CMAKE_CURRENT_SOURCE_DIR}/mmpol_vmc.f
    ${CMAKE_CURRENT_SOURCE_DIR}/multiple_states.f
    ${CMAKE_CURRENT_SOURCE_DIR}/determinant_psig.f)
  
  add_library(vmc_mpi_objects OBJECT ${SHARED_VMC_MPI_SOURCE})
  
  target_include_directories(vmc_mpi_objects
    PUBLIC
    ${CMAKE_SOURCE_DIR}/include)
  
  target_compile_options(vmc_mpi_objects
    PRIVATE
    ${Fortran_FLAGS}
    "$<$<CONFIG:Debug>:${Fortan_FLAGS_DEBUG}>")
  
  foreach(s ${SHARED_VMC_MPI_SOURCE})
    list(REMOVE_ITEM VMC_SOURCE ${s})
  endforeach(s)
  foreach(s ${SHARED_VMC_MPI_SOURCE})
    list(REMOVE_ITEM VMC90_SOURCE ${s})
  endforeach(s)
  
  # --------------------------------------------
  # creates object library for dmc.mov1
  # --------------------------------------------
  
  set(SHARED_DMC_SOURCE
     ${CMAKE_CURRENT_SOURCE_DIR}/error.f
     ${CMAKE_CURRENT_SOURCE_DIR}/optwf_handle_wf.f
     ${CMAKE_CURRENT_SOURCE_DIR}/optwf_sr_more.f)
  
  add_library(dmc_objects OBJECT ${SHARED_DMC_SOURCE})
  
  target_include_directories(dmc_objects
    PUBLIC
    ${CMAKE_SOURCE_DIR}/include)
  
  target_compile_options(dmc_objects
    PRIVATE
    ${Fortran_FLAGS}
    "$<$<CONFIG:Debug>:${Fortan_FLAGS_DEBUG}>")
  
  foreach(s ${SHARED_DMC_SOURCE})
    list(REMOVE_ITEM VMC_SOURCE ${s})
  endforeach(s)
  
  # --------------------------------------------
  # Fortran Modules
  # --------------------------------------------
  
  add_library(vmc_modules optgeo_convcartzmat.f90 optgeo_coords_int.f90 optgeo_hessian.f90)
  
  # --------------------------------------------
  # Compile and link executable
  # --------------------------------------------
  
  add_executable(vmc.mov1
    $<TARGET_OBJECTS:shared_objects_1>
    $<TARGET_OBJECTS:shared_objects_2>
    $<TARGET_OBJECTS:vmc_mpi_objects>
    $<TARGET_OBJECTS:dmc_objects>
    # source that is not already included in the objects
    ${VMC_SOURCE} 
    ${VMC90_SOURCE})
  
  #  Add headers
  target_include_directories(vmc.mov1
    PUBLIC
    ${CMAKE_SOURCE_DIR}/include)
  
  # compilation
  target_compile_options(vmc.mov1
    PRIVATE
    ${Fortran_FLAGS}
    "$<$<CONFIG:Debug>:${Fortan_FLAGS_DEBUG}>"
    )
  
  target_link_libraries(vmc.mov1
    vmc_modules
    cyrus pspline
    ${LINEAR_ALGEBRA})
  
  # TODO remove after I took care of DMC
  set(VMC_SOURCE ${VMC_SOURCE} CACHE INTERNAL "VMC source files")
  
  
  # --------------------------------------------
  # Everything related to miniparser p2
  # --------------------------------------------
  
  add_custom_command(
    OUTPUT
      ${INPUT_DIR}/p2nmcheck.f
    OUTPUT
      ${CMAKE_CURRENT_SOURCE_DIR}/p2prog.f
    WORKING_DIRECTORY
      ${CMAKE_CURRENT_SOURCE_DIR}
    COMMAND
      ${INPUT_DIR}/inpt.awk -v p2file=commands.p2 -v logf=${INPUT_DIR}/commands.log ${CMAKE_CURRENT_SOURCE_DIR}/read_input.f
    COMMAND
      cat commands.p2 ${INPUT_DIR}/end.p2 ${INPUT_DIR}/auto.p2 > all.p2
    COMMAND
      ${INPUT_DIR}/genp2_defaults.awk -v ofile=${CMAKE_CURRENT_SOURCE_DIR}/p2prog.f all.p2 > genp2.log
    COMMAND
      ${INPUT_DIR}/mkp2nmcheck.pl namelist.p2 > ${INPUT_DIR}/p2nmcheck.f
    VERBATIM
  )
  
  add_custom_target(
    generate-p2
    ${CMAKE_COMMAND} -E echo "Generating p2 files."
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/p2prog.f
    DEPENDS ${INPUT_DIR}/p2nmcheck.f
  )
  
  set_source_files_properties(
    ${INPUT_DIR}/p2nmcheck.f
    PROPERTIES GENERATED TRUE
  )
  
  set_source_files_properties(
    ${CMAKE_CURRENT_SOURCE_DIR}/p2prog.f
    PROPERTIES GENERATED TRUE
  )
  
  # --------------------------------------------
  # creates object library for dmc.mov1.mpi
  # --------------------------------------------
  
    set(SHARED_DMC_MPI_SOURCE
      ${CMAKE_CURRENT_SOURCE_DIR}/error.f
      ${CMAKE_CURRENT_SOURCE_DIR}/optx_jas_ci_reduce.f
      ${CMAKE_CURRENT_SOURCE_DIR}/optx_jas_orb_reduce.f
      ${CMAKE_CURRENT_SOURCE_DIR}/optx_orb_ci_reduce.f
      ${CMAKE_CURRENT_SOURCE_DIR}/optci_reduce.f
      ${CMAKE_CURRENT_SOURCE_DIR}/optjas_reduce.f
      ${CMAKE_CURRENT_SOURCE_DIR}/optorb_reduce.f
      ${CMAKE_CURRENT_SOURCE_DIR}/finwrt_more.f
      ${CMAKE_CURRENT_SOURCE_DIR}/optwf_handle_wf.f
      ${CMAKE_CURRENT_SOURCE_DIR}/optwf_sr_more.f)
  
    add_library(vmc_dmc_mpi_objects OBJECT ${SHARED_DMC_MPI_SOURCE})
  
    target_include_directories(vmc_dmc_mpi_objects
      PUBLIC
      ${CMAKE_SOURCE_DIR}/include)
  
    target_compile_options(vmc_dmc_mpi_objects
      PRIVATE
      ${Fortran_FLAGS}
      "$<$<CONFIG:Debug>:${Fortan_FLAGS_DEBUG}>")
  
    foreach(s ${SHARED_DMC_MPI_SOURCE})
      list(REMOVE_ITEM VMC_SOURCE ${s})
    endforeach(s)
  
  # --------------------------------------------
    
    add_executable(vmc.mov1.mpi
      # MPI source code
      ${VMC_SOURCE} ${VMC90_MPI_SOURCE}
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
    set(VMC_SOURCE ${VMC_SOURCE} CACHE INTERNAL "MPI vmc source files")

endif()
