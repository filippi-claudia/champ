# Source files
file(GLOB DMC_SOURCE LIST_DIRECTORIES false  *.f)

include(${CMAKE_MODULE_PATH}/FortranFormat.cmake)
formatFortran("${DMC_SOURCE}") # Set Fixed format property

# --------------------------------------------
# creates object library for dmc.mov1.mpi
# --------------------------------------------

set(SHARED_DMC_MPI_SOURCE
  ${CMAKE_CURRENT_SOURCE_DIR}/averages.f
  ${CMAKE_CURRENT_SOURCE_DIR}/dmc_good_ps_mov1.f
  ${CMAKE_CURRENT_SOURCE_DIR}/hpsiedmc.f
  ${CMAKE_CURRENT_SOURCE_DIR}/multideterminant_tmove.f
  ${CMAKE_CURRENT_SOURCE_DIR}/optwf_dmc.f
  ${CMAKE_CURRENT_SOURCE_DIR}/pcm_dmc.f
  ${CMAKE_CURRENT_SOURCE_DIR}/splitj.f
  ${CMAKE_CURRENT_SOURCE_DIR}/walksav_jas.f
  ${CMAKE_CURRENT_SOURCE_DIR}/dmc.f
  ${CMAKE_CURRENT_SOURCE_DIR}/gauss.f
  ${CMAKE_CURRENT_SOURCE_DIR}/mmpol_dmc.f
  ${CMAKE_CURRENT_SOURCE_DIR}/nonloc_grid.f
  ${CMAKE_CURRENT_SOURCE_DIR}/optwf_dummy.f
  ${CMAKE_CURRENT_SOURCE_DIR}/prop_dmc.f
  ${CMAKE_CURRENT_SOURCE_DIR}/walksav_det.f)

add_library(dmc_mpi_objects OBJECT ${SHARED_DMC_MPI_SOURCE})

target_include_directories(dmc_mpi_objects
  PUBLIC
  ${CMAKE_SOURCE_DIR}/include)

target_compile_options(dmc_mpi_objects
  PRIVATE
  ${Fortran_FLAGS}
  "$<$<CONFIG:Debug>:${Fortan_FLAGS_DEBUG}>")

set(DMC_SOURCE_REDUCED "${DMC_SOURCE}")
foreach(s ${SHARED_DMC_MPI_SOURCE})
  list(REMOVE_ITEM DMC_SOURCE_REDUCED ${s})
endforeach(s)

# --------------------------------------------

# executable dependencies
add_executable(dmc.mov1
  $<TARGET_OBJECTS:shared_objects_1>
  $<TARGET_OBJECTS:shared_objects_2>
  $<TARGET_OBJECTS:dmc_objects>
  $<TARGET_OBJECTS:dmc_mpi_objects>
  # DMC source code
  ${DMC_SOURCE_REDUCED})

#  Add headers
target_include_directories(dmc.mov1
  PUBLIC
  ${CMAKE_SOURCE_DIR}/include)

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
