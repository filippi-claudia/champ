# Source files
file(GLOB PARSER_SOURCE
	fdf.F90 
	io_fdf.F90 
	parse.F90 
	prec.F90 
	utils.F90	
	)

# add_library(parser OBJECT ${PARSER_SOURCE})
# target_compile_options(parser
#   PUBLIC
#   ${Fortran_FLAGS}
#   "$<$<CONFIG:Debug>:${Fortan_FLAGS_DEBUG}>")


