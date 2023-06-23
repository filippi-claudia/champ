
# Source files
file(GLOB PARSER_SOURCE
	fdf.F90
	io_fdf.F90
	parse.F90
	prec.F90
	utils.F90
	)

add_library(parser SHARED ${PARSER_SOURCE})

set(N_Fortran_FLAGS "${Fortran_Flags}")
list(REMOVE_ITEM N_Fortran_FLAGS "-Werror")
target_compile_options(parser
  PRIVATE
  ${N_Fortran_FLAGS})

link_libraries(parser)
