if (MPI_Fortran_FOUND)

  file(GLOB DMC_MPI_GPOPB_SOURCE LIST_DIRECTORIES false  *.f)

    # Remove the files that are not used from dmc
  set(PATH_DMC ${CMAKE_SOURCE_DIR}/dmc)
  list(REMOVE_ITEM DMC_SOURCE ${PATH_DMC}/acuest.f ${PATH_DMC}/acues1_reduce.f ${PATH_DMC}/dumper.f
    ${PATH_DMC}/main.f ${PATH_DMC}/finwrt.f ${PATH_DMC}/optwf_dummy.f
    ${PATH_DMC}/walksav_det.f ${PATH_DMC}/walksav_jas.f)

  # Remove files from vmc
  set(PATH_VMC ${CMAKE_SOURCE_DIR}/vmc)
  list(REMOVE_ITEM VMC_SOURCE ${PATH_VMC}/error.f ${PATH_VMC}/optjas.f ${PATH_VMC}/optorb.f ${PATH_VMC}/optci.f
    ${PATH_VMC}/optx_jas_orb.f ${PATH_VMC}/optx_jas_ci.f ${PATH_VMC}/optx_orb_ci.f) 
  
  #Remove the files that are not used from vmc/mpi
  set(PATH_VMC_MPI ${CMAKE_SOURCE_DIR}/vmc/MPI)
  list(REMOVE_ITEM VMC_MPI_SOURCE
    ${PATH_VMC_MPI}/main.f ${PATH_VMC_MPI}/acuest_reduce.f ${PATH_VMC_MPI}/dumper.f
    ${PATH_VMC_MPI}/dl_more.f ${PATH_VMC_MPI}/fin_reduce.f ${PATH_VMC_MPI}/force_analy_reduce.f
    ${PATH_VMC_MPI}/jdqz_main.f ${PATH_VMC_MPI}/lin_d_more.f ${PATH_VMC_MPI}/lin_d_more_no_semiortho.f
    ${PATH_VMC_MPI}/mc_configs.f ${PATH_VMC_MPI}/mmpol_reduce.f ${PATH_VMC_MPI}/pcm_reduce.f
    ${PATH_VMC_MPI}/prop_reduce.f ${PATH_VMC_MPI}/regterg.f90 ${PATH_VMC_MPI}/sr_more.f)
  
  # DMC MPI dependencies
  set(PATH_DMC_MPI ${CMAKE_SOURCE_DIR}/dmc/MPI)
  set(DMC_MPI_SOURCE ${PATH_DMC_MPI}/mc_configs.f)
  
  # Include function to change fortran format
  include(${CMAKE_MODULE_PATH}/FortranFormat.cmake)
  formatFortran("${DMC_MPI_SOURCE}" "FIXED")

  add_executable(dmc.mov1.gpopb.mpi
    # VMC dependencies
    ${VMC_SOURCE}
    # VMC Fortran 90
    ${CMAKE_SOURCE_DIR}/vmc/convcartzmat.f90
    # # MPI VMC code dependencies
    ${VMC_MPI_SOURCE}
    # DMC sequential code
    ${DMC_SOURCE}
    # VMC Fortran 90
    ${CMAKE_SOURCE_DIR}/vmc/convcartzmat.f90
    # code in the current directory 
    ${DMC_MPI_GPOPB_SOURCE}
    #Input parser
    ${INPUT_DIR}/p2_defaults.f ${INPUT_DIR}/p2etc.f ${INPUT_DIR}/p2nmcheck.f)

  foreach(_source ${VMC_SOURCE})
    set_property(SOURCE ${_source} PROPERTY Fortran_FORMAT FIXED)
  endforeach(_source)
  foreach(_source ${VMC_MPI_SOURCE})
    set_property(SOURCE ${_source} PROPERTY Fortran_FORMAT FIXED)
  endforeach(_source)
  foreach(_source ${DMC_SOURCE})
    set_property(SOURCE ${_source} PROPERTY Fortran_FORMAT FIXED)
  endforeach(_source)
  foreach(_source ${DMC_MPI_GPOPB_SOURCE})
    set_property(SOURCE ${_source} PROPERTY Fortran_FORMAT FIXED)
  endforeach(_source)
  
  # Add dependency to vmc mpi
  add_dependencies(dmc.mov1.gpopb.mpi dmc.mov1.mpi)

    #  Add headers
  target_include_directories(dmc.mov1.gpopb.mpi
    PRIVATE
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/input)

  # compilation
  target_compile_options(dmc.mov1.gpopb.mpi
    PRIVATE
    ${Fortran_FLAGS}
    "$<$<CONFIG:Debug>:${Fortan_FLAGS_DEBUG}>")

  # Link dependencies
  target_link_libraries(dmc.mov1.gpopb.mpi
    # PRIVATE
    cyrus pspline
    # PUBLIC
    ${LINEAR_ALGEBRA}
    ${MPI_Fortran_LIBRARIES})

  set(DMC_MPI_GPOPB_DEPENDENCIES ${VMC_SOURCE} ${VMC_MPI_SOURCE} ${DMC_SOURCE} ${DMC_MPI_GPOPB_SOURCE}
    CACHE INTERNAL "source files used by targets depending on dmc.mov1.gpopb.mpi")  

  set_source_files_properties(
    ${INPUT_DIR}/p2nmcheck.f
    PROPERTIES GENERATED TRUE
  )
  
  set_source_files_properties(
    ${CMAKE_CURRENT_SOURCE_DIR}/p2prog.f
    PROPERTIES GENERATED TRUE
  )

endif()
