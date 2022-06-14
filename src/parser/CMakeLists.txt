
# Source files
file(GLOB PARSER_SOURCE
	fdf.F90
	io_fdf.F90
	parse.F90
	prec.F90
	utils.F90
	)

add_library(parser SHARED ${PARSER_SOURCE})

target_compile_options(parser
  PRIVATE
  ${Fortran_FLAGS})

link_libraries(parser)
