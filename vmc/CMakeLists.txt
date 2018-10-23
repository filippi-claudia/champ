file(GLOB_RECURSE VMC_SOURCES *.f *.f90)


## BUILD pspline by hand
add_custom_command(
  OUTPUT libpspline.a
  COMMAND  make
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/lib2/pspline)

add_executable(vmc.mov1 ${VMC_SOURCES})

target_link_libraries(vmc.mov1 ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} libpspline.a)


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

