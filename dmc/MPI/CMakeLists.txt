if (MPI_Fortran_FOUND)
  list(APPEND DMC_MPI_SOURCE main.f acuest.f acues1_reduce.f dumper.f finwrt.f
    mc_configs.f mmpol_reduce.f pcm_reduce.f prop_reduce.f)

  # Remove the files that are not used from dmc
  set(PATH_DMC ${CMAKE_SOURCE_DIR}/dmc)
  list(REMOVE_ITEM DMC_SOURCE ${PATH_DMC}/acuest.f ${PATH_DMC}/acues1_reduce.f ${PATH_DMC}/dumper.f
    ${PATH_DMC}/main.f ${PATH_DMC}/mc_configs.f ${PATH_DMC}/finwrt.f)


  # Remove files from vmc
  list(REMOVE_ITEM VMC_SOURCE ${CMAKE_SOURCE_DIR}/vmc/error.f)
  
  #Remove the files that are not used from vmc/mpi
  set(PATH_VMC_MPI ${CMAKE_SOURCE_DIR}/vmc/MPI)
  list(REMOVE_ITEM VMC_MPI_SOURCE
    ${PATH_VMC_MPI}/main.f ${PATH_VMC_MPI}/acuest_reduce.f ${PATH_VMC_MPI}/dumper.f
    ${PATH_VMC_MPI}/dl_more.f ${PATH_VMC_MPI}/fin_reduce.f ${PATH_VMC_MPI}/force_analy_reduce.f
    ${PATH_VMC_MPI}/jdqz_main.f ${PATH_VMC_MPI}/lin_d_more.f ${PATH_VMC_MPI}/lin_d_more_no_semiortho.f
    ${PATH_VMC_MPI}/mc_configs.f ${PATH_VMC_MPI}/mmpol_reduce.f ${PATH_VMC_MPI}/pcm_reduce.f
    ${PATH_VMC_MPI}/prop_reduce.f ${PATH_VMC_MPI}/regterg.f90 ${PATH_VMC_MPI}/sr_more.f)
  
  # REmove unused files
  # Include function to change fortran format
  include(${CMAKE_MODULE_PATH}/FortranFormat.cmake)
  formatFortran("${DMC_MPI_SOURCE}" "FIXED")

  add_executable(dmc.mov1.mpi
    # VMC dependencies
    ${VMC_SOURCE}
    # VMC Fortran 90
    ${CMAKE_SOURCE_DIR}/vmc/convcartzmat.f90
    # # MPI VMC code dependencies
    ${VMC_MPI_SOURCE}
    # DMC sequential code
    ${DMC_SOURCE}
    # code in the current directory 
    ${DMC_MPI_SOURCE}
    #Input parser
    ${INPUT_DIR}/p2_defaults.f ${INPUT_DIR}/p2etc.f ${INPUT_DIR}/p2nmcheck.f)
  
  # Add dependency to vmc mpi
  add_dependencies(dmc.mov1.mpi vmc.mov1.mpi dmc.mov1)

    #  Add headers
  target_include_directories(dmc.mov1.mpi
    PRIVATE
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/input)

  # compilation
  target_compile_options(dmc.mov1.mpi
    PRIVATE
    ${Fortran_FLAGS}
    "$<$<CONFIG:Debug>:${Fortan_FLAGS_DEBUG}>")

  # Link dependencies
  target_link_libraries(dmc.mov1.mpi
    # PRIVATE
    cyrus pspline
    # PUBLIC
    ${LINEAR_ALGEBRA}
    ${MPI_Fortran_LIBRARIES})

endif()
