# Source files
set(DMC_SOURCE
  ${CMAKE_CURRENT_SOURCE_DIR}/main.f ${CMAKE_CURRENT_SOURCE_DIR}/dmc.f
  ${CMAKE_CURRENT_SOURCE_DIR}/mc_configs.f ${CMAKE_CURRENT_SOURCE_DIR}/gauss.f
  ${CMAKE_CURRENT_SOURCE_DIR}/prop_dmc.f ${CMAKE_CURRENT_SOURCE_DIR}/optwf_dmc.f
  ${CMAKE_CURRENT_SOURCE_DIR}/pcm_dmc.f  ${CMAKE_CURRENT_SOURCE_DIR}/mmpol_dmc.f
  ${CMAKE_CURRENT_SOURCE_DIR}/averages.f ${CMAKE_CURRENT_SOURCE_DIR}/dmc_good_ps_mov1.f
  ${CMAKE_CURRENT_SOURCE_DIR}/acuest.f ${CMAKE_CURRENT_SOURCE_DIR}/finwrt.f
  ${CMAKE_CURRENT_SOURCE_DIR}/dumper.f ${CMAKE_CURRENT_SOURCE_DIR}/splitj.f
  ${CMAKE_CURRENT_SOURCE_DIR}/walksav_det.f ${CMAKE_CURRENT_SOURCE_DIR}/walksav_jas.f
  ${CMAKE_CURRENT_SOURCE_DIR}/hpsiedmc.f ${CMAKE_CURRENT_SOURCE_DIR}/nonloc_grid.f
  ${CMAKE_CURRENT_SOURCE_DIR}/multideterminant_tmove.f ${CMAKE_CURRENT_SOURCE_DIR}/optwf_dummy.f
  ${CMAKE_CURRENT_SOURCE_DIR}/acues1_reduce.f)

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
#include(${CMAKE_MODULE_PATH}/FortranFormat.cmake)
#formatFortran("${VMC90_DEPENDENCIES}" "FREE")
#formatFortran("${VMC_DEPENDENCIES}" "FIXED")
#formatFortran("${DMC_SOURCE}" "FIXED")
foreach(_source ${DMC_SOURCE})
  set_property(SOURCE ${_source} PROPERTY Fortran_FORMAT FIXED)
endforeach(_source)
foreach(_source ${VMC_SOURCE})
  set_property(SOURCE ${_source} PROPERTY Fortran_FORMAT FIXED)
endforeach(_source)

# compilation
target_compile_options(dmc.mov1
  PRIVATE
  ${Fortran_FLAGS}
  "$<$<CONFIG:Debug>:${Fortan_FLAGS_DEBUG}>"
  )

target_link_libraries(dmc.mov1
  cyrus pspline
  ${LINEAR_ALGEBRA})

set(DMC_SOURCE ${DMC_SOURCE} CACHE INTERNAL "Path to the source of vmc f77")
add_subdirectory(MPI)
add_subdirectory(MPI_global_pop_big)
add_subdirectory(MPI_global_pop)

set_source_files_properties(
  ${INPUT_DIR}/p2nmcheck.f
  PROPERTIES GENERATED TRUE
)

set_source_files_properties(
  ${CMAKE_CURRENT_SOURCE_DIR}/p2prog.f
  PROPERTIES GENERATED TRUE
)
