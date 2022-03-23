
add_subdirectory(parser)
add_subdirectory(shared)
add_subdirectory(vmc)
add_subdirectory(dmc)
add_subdirectory(module)

set_property(
SOURCE shared/new_parser_read_data.f90 DIRECTORY vmc dmc
APPEND
PROPERTY COMPILE_DEFINITIONS
CMAKE_Fortran_COMPILER=\"${CMAKE_Fortran_COMPILER}\"
CMAKE_Fortran_COMPILER_VERSION=\"${CMAKE_Fortran_COMPILER_VERSION}\"
GIT_HEAD_BRANCH=\"${GIT_HEAD_BRANCH}\"
GIT_REVISION_HASH=\"${GIT_REVISION_HASH}\"
HDF5_VERSION=\"${HDF5_VERSION}\"
)

