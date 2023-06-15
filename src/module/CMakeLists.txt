if (MPI_Fortran_FOUND)

  # --------------------------------------------------------------------
  # Create two libraries to share objects between the two binaries,
  # vmc.mov1 and dmc.mov1:
  # - shared_modules: all .f90 modules that must be compiled first.
  # - shared_objects: rest of shared objects (.f90/.f).
  # --------------------------------------------------------------------
  set(SHARED_MODULES
    m_kinds.f90
    m_metropolis.f90
    m_slater.f90
    m_basis.f90
    m_common.f90
    m_periodic_table.F90
    m_atom.f90
    m_timer.f90
    m_det.f90
    m_coefs.f90
    m_constant.f90
    m_control.f90
    m_slater.f90
    m_force.f90
    m_mstates.f90
    m_mpi.f90
    m_optci.f90
    m_pseudo.f90
    m_sampling.f90
    m_sr.f90
    m_slater.f90
    m_metropolis.f90
    m_dmc.f90
    m_vmc.f90
    m_grid.f90
    m_pcm.f90
    m_ewald.f90
    m_jastrow.f90
    m_mmpol.f90
    m_optorb.f90
    m_prop.f90
    m_efield.f90
    m_gradhess.f90
    m_optwf.f90
    m_optwf_func.f90
    m_estimators.f90
    m_deriv.f90
    m_mixderiv.f90
    m_input.f90
    m_grdnt.f90
    m_allocation.f90
    m_state_avrg.f90
    m_error.f90
    ${CMAKE_PARSER_DIR}/fdf.F90
    ${CMAKE_PARSER_DIR}/io_fdf.F90
    ${CMAKE_PARSER_DIR}/parse.F90
    ${CMAKE_PARSER_DIR}/prec.F90
    ${CMAKE_PARSER_DIR}/utils.F90)

    if (HDF5_FOUND)
        list(APPEND SHARED_MODULES ${CMAKE_CURRENT_SOURCE_DIR}/utils_hdf5.f90)
    endif()

    if (TREXIO_FOUND)
      list(APPEND SHARED_MODULES ${TREXIO_INCLUDE_DIR}/trexio_f.f90)
    endif()

    if (QMCKL_FOUND)
       list(APPEND SHARED_MODULES m_qmckl.F90)
    endif()

    add_library(shared_modules OBJECT ${SHARED_MODULES})

    if (HDF5_FOUND)
	    target_include_directories(shared_modules PRIVATE ${HDF5_INCLUDE_DIRS})
      target_link_libraries(shared_modules ${HDF5_LIBRARIES})
    endif(HDF5_FOUND)

    target_compile_options(shared_modules
      PRIVATE
      ${Fortran_FLAGS}
      "$<$<CONFIG:Debug>:${Fortan_FLAGS_DEBUG}>")

endif()
