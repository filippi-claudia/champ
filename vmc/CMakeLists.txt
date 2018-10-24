# Source files
file(GLOB VMC_SOURCES LIST_DIRECTORIES false  *.f)
file(GLOB VMC90_SOURCES LIST_DIRECTORIES false *f90)

# SET periodic and qmmm to false
set(PERIODIC F)
set(QMMM F)

## BUILD pspline by hand
add_custom_command(
  OUTPUT libpspline.a
  COMMAND  make
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/lib2/pspline)

# # Added the pspline library to the link_directories
# INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/bin)
# LINK_DIRECTORIES(${CMAKE_SOURCE_DIR}/bin)
# message(STATUS "LINK_DIRECTORIES: " ${LINK_DIRECTORIES})

## Compiler FLAGS
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  list(APPEND Fortran_FLAGS "-O2" "-ffixed-line-length-132" "-fno-automatic")
  
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  list(APPEND Fortran_FLAGS "-O2" "-xHost" "-mcmodel=large" "-shared-intel" "-132" "-save "
    "-dyncom=grid3d_data,orbital_num_spl,orbital_num_lag,orbital_num_spl2,grid3d_data"
    "-check bounds" "-check format ")
  list(APPEND Fortan_FLAGS_DEBUG "-g")
  
elseif(Fortran_COMPILER_ID MATCHES "PGI")
  list(APPEND CMAKE_Fortran_FLAGS "-Mfreeform -Mdclchk -Mstandard -Mallocatable=03")
endif()


# Compile and link executable
add_executable(vmc.mov1 ${VMC_SOURCES} ${VMC90_SOURCES} ${CMAKE_SOURCE_DIR}/include/vmc.h)

#  Add headers
target_include_directories(vmc.mov1
  PRIVATE ${CMAKE_SOURCE_DIR}/include
  PRIVATE ${CMAKE_SOURCE_DIR}/input)


# Set Fixed or Free format
message(STATUS "Setting source format properties")
  # Intel does not care about the fixed format
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  foreach(_source ${VMC90_SOURCES})
    set_source_files_properties(${_source} PROPERTIES FORTRAN_FORMAT FREE)
    message (STATUS "${_source} is Free format")
  endforeach()
  foreach(_source ${VMC_SOURCES})
    set_source_files_properties(${_source} PROPERTIES FORTRAN_FORMAT FIXED)
    message (STATUS "${_source} is Fixed format")
  endforeach()
endif()

# compilation
target_compile_options(vmc.mov1
  PRIVATE
  ${Fortran_FLAGS}
  "$<$<CONFIG:Debug>:${Fortan_FLAGS_DEBUG}>"
  )

target_link_libraries(vmc.mov1 cyrus libpspline.a ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})


# ifndef FFLAGS


# INSTALLDIR = $(QMCHOME)/bin
# INC        = -I$(INPUT) -I$(QMCHOME)/include


# include $(QMCHOME)/settings.make.local
# endif

# FFLAGS := $(FFLAGS) $(INC)

# include $(QMCHOME)/common_objects


# commands.p2: $V/read_input.f 
# 	$(INPUT)/inpt.awk  -v p2file=$@ -v logf=$(INPUT)/commands.log $^

# $V/p2prog.f: commands.p2 $(INPUT)/end.p2 $(INPUT)/auto.p2
# 	cat commands.p2 $(INPUT)/end.p2 $(INPUT)/auto.p2  > all.p2
# 	$(INPUT)/genp2_defaults.awk  -v ofile=$@ all.p2 > genp2.log

# $(INPUT)/p2nmcheck.f: namelist.p2
# 	$(INPUT)/mkp2nmcheck.pl $^ > $@

# vmc.mov1: $(OVMC1)
# 	$(FC) $(FFLAGS) $(OVMC1) $(CYRUSLIB) $(LINPACK) $(LAPACK) $(PSPLINE) $(BLAS) -o $@

