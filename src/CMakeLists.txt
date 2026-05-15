
add_subdirectory(parser)
add_subdirectory(module)
add_subdirectory(vmc)

if(MPI_Fortran_FOUND)
  add_subdirectory(dmc)
endif()
