# This Makefile is to be used with gmake.  On Linux systems, gmake is
# the default make.

include settings.make

system := $(shell uname)
hostname := $(shell uname -n)

ifeq ($(FC),xlf90_r)
 my_sec = my_second.f.xlf
else 
 my_sec = my_second.f.ifort
endif

OBJS = my_second.o ranf.o shellnn.o intpol.o intpol_monot.o grater.o \
hartre.o deriv.o deriv2.o divergence.o lapla3.o laplac.o fderiv.o simson.o bode.o \
hartre_grid2.o deriv_grid2.o deriv2_grid2.o divergence_grid2.o laplac_grid2.o fderiv_grid2.o \
bode_grid2.o

libcyrus.a: link_my_second $(OBJS) settings.make
	ar rs $@ $(OBJS)

clean:
	rm -f $(OBJS)
	rm -f libcyrus.a

link_my_second:
	ln -sf $(my_sec) my_second.f
