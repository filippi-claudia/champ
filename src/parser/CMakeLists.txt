add_library(parser)

target_sources(parser 
	PUBLIC
	fdf.F90 
	io_fdf.F90 
	parse.F90 
	prec.F90 
	utils.F90	
)


# Source files
file(GLOB PARSER_SOURCE
	fdf.F90 
	io_fdf.F90 
	parse.F90 
	prec.F90 
	utils.F90	
	)

