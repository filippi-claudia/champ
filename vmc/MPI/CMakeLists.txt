if (MPI_Fortran_FOUND)
  file(GLOB VMC_MPI_SOURCE LIST_DIRECTORIES false  *.f)
  file(GLOB VMC90_MPI_SOURCE LIST_DIRECTORIES false  *.f90)

  # Remove unused files
  list(REMOVE_ITEM VMC_MPI_SOURCE
    ${CMAKE_CURRENT_SOURCE_DIR}/lin_d_more_no_semiortho.f
    ${CMAKE_CURRENT_SOURCE_DIR}/lin_d_more_semiortho.f
    )
  
  # Remove sequential source from dependencies
  set(PATH_VMC ${CMAKE_SOURCE_DIR}/vmc)
  list(REMOVE_ITEM VMC_SOURCE
    ${PATH_VMC}/acuest_reduce.f ${PATH_VMC}/deriv_jastrow.f 
    ${PATH_VMC}/dl_more.f ${PATH_VMC}/dumper.f ${PATH_VMC}/error.f ${PATH_VMC}/fin_reduce.f
    ${PATH_VMC}/finwrt_more.f ${PATH_VMC}/jdqz_main.f
    ${PATH_VMC}/lin_d_more.f ${PATH_VMC}/main.f ${PATH_VMC}/mc_configs.f
    ${PATH_VMC}/pw_ewald.f ${PATH_VMC}/pw_find_image.f
    ${PATH_VMC}/pw_orbitals.f ${PATH_VMC}/pw_orbitals_e.f ${PATH_VMC}/pw_read.f
    ${PATH_VMC}/regterg.f90 ${PATH_VMC}/sr_more.f)
  
  add_executable(vmc.mov1.mpi
    # Header
    ${CMAKE_SOURCE_DIR}/include/vmc.h
    # MPI source code
    ${VMC_MPI_SOURCE} ${VMC90_MPI_SOURCE}
    # shared code with VMC
    ${VMC_SOURCE}
    # shared VMC Fortran 90 source
    ${PATH_VMC}/convcartzmat.f90
    # Input parser
    ${INPUT_DIR}/p2_defaults.f ${INPUT_DIR}/p2etc.f ${INPUT_DIR}/p2nmcheck.f)

  # Add depency with vmc
  add_dependencies(vmc.mov1.mpi vmc.mov1)
  
  ## Set Fixed or Free format
  # Include function to change fortran format
  include(${CMAKE_MODULE_PATH}/FortranFormat.cmake)
  formatFortran("${VMC_MPI_SOURCES}" "FIXED")
  
  #  Add headers
  target_include_directories(vmc.mov1.mpi
    PUBLIC
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/input)
  
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
  set(VMC_MPI_SOURCE ${VMC_MPI_SOURCE} CACHE INTERNAL "MPI vmc source files")

endif()
