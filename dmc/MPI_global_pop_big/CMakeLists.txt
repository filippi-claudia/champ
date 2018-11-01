if ((MPI_FOUND) AND (NOT CMAKE_Fortran_COMPILER_ID MATCHES "GNU"))
file(GLOB DMC_MPI_GLOBAL_BIG LIST_DIRECTORIES false *.f)


# Include function to change fortran format
include(${CMAKE_MODULE_PATH}/FortranFormat.cmake)
formatFortran("${DMC_MPI_GLOBAL_BIG}" "FIXED")


add_executable(dmc.mov1.gpopb.mpi
  ${CMAKE_SOURCE_DIR}/include
  ${DMC_MPI_GLOBAL_BIG})

target_include_directories(dmc.mov1.gpopb.mpi
  PUBLIC
  ${CMAKE_SOURCE_DIR}/include
)

target_link_libraries(dmc.mov1.gpopb.mpi
  cyrus pspline
  ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES}
  MPI::MPI_Fortran)
endif()
