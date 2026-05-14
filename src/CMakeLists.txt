
add_subdirectory(parser)
add_subdirectory(vmc)
add_subdirectory(module)

if (MPI_Fortran_FOUND)
  add_subdirectory(dmc)
endif()
