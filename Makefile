include settings.make

system    := $(shell uname)
hostname  := $(shell uname -n)
QMCHOME    = $(shell pwd)
LIBHOME    = $(QMCHOME)
INPUT      = $(QMCHOME)/input

CYRUSLIB   = -L$(LIBHOME)/lib -lcyrus

BLASHOME   = $(LIBHOME)/lib2/blas-3.6.0
LAPACKHOME = $(LIBHOME)/lib2/lapack-3.6.0
BLAS       = -L$(BLASHOME) -lblas
LAPACK     = -L$(LAPACKHOME) -llapack

#BLAS        =
#LAPACK      = -mkl=sequential

PSPLINE    = -L$(LIBHOME)/lib2/pspline -lpspline

INSTALLDIR = $(QMCHOME)/bin
INC        = -I$(INPUT) -I$(QMCHOME)/include

export FC FC_90 F90FLAGS FFLAGS FC_MPI FFLAGS_MPI FLINKER_MPI \
       QMCHOME LIBHOME BLASHOME LAPACKHOME INPUT INC \
       PSPLINE BLAS LAPACK CYRUSLIB \
       PERIODIC QMMM

###################################
# Targets, Dependencies and Rules #
###################################

# PHONY is a Make keyword that prevents it from getting confused if there happens to be a file named "clean" in the directory.
 
.PHONY: vmc.mov1 dmc.mov1 vmc.MPI dmc.MPI dmc.MPI.global.pop.big dmc.MPI.global.pop clean clean_all
.PHONY: seq mpi all lib lib2

seq: lib lib2 vmc.mov1 dmc.mov1 

mpi: vmc.MPI dmc.MPI dmc.MPI.global.pop.big dmc.MPI.global.pop

all: seq mpi

lib:
	cd lib ; $(MAKE)

lib2:
	cd lib2/pspline ; $(MAKE) libpspline.a
	cd $(BLASHOME) ; $(MAKE)
	cd $(LAPACKHOME) ; $(MAKE) lapacklib

vmc.mov1: lib lib2 vmc
	cd vmc ; $(MAKE) vmc.mov1

dmc.mov1: lib lib2 vmc.mov1
	cd dmc ; $(MAKE) dmc.mov1

vmc.MPI: lib lib2 vmc.mov1 
	cd vmc/MPI ; $(MAKE) vmc.mov1

dmc.MPI: lib lib2 vmc.mov1 dmc.mov1
	cd dmc/MPI ; $(MAKE) dmc.mov1

dmc.MPI.global.pop.big: lib lib2 vmc.mov1 dmc.mov1
	cd dmc/MPI_global_pop_big ; $(MAKE) dmc.mov1

dmc.MPI.global.pop: lib lib2 vmc.mov1 dmc.mov1 dmc.MPI.global.pop.big
	cd dmc/MPI_global_pop ; $(MAKE) dmc.mov1

clean:
	-rm -f input/*.o vmc/*.o vmc/*.mod dmc/*.o vmc/MPI/*.o dmc/MPI*/*.o 
	if test -d $(INSTALLDIR); then rm -f $(INSTALLDIR)/*; fi

clean_all: clean
	-rm -f vmc/vmc.mov1 dmc/dmc.mov1 vmc/MPI/vmc.mov1 dmc/MPI*/dmc.mov1
	cd lib ; $(MAKE) clean
	cd $(BLASHOME) ; $(MAKE) clean
	cd $(LAPACKHOME) ; $(MAKE) clean
	cd lib2/pspline ; $(MAKE) clean

install: 
	cp -f vmc/vmc.mov1               $(INSTALLDIR)
	cp -f vmc/MPI/vmc.mov1           $(INSTALLDIR)/vmc.mov1.mpi
	cp -f dmc/dmc.mov1               $(INSTALLDIR)
	cp -f dmc/MPI/dmc.mov1           $(INSTALLDIR)/dmc.mov1.mpi
	cp -f dmc/MPI_global_pop/dmc.mov1  $(INSTALLDIR)/dmc.mov1.gpop.mpi
	cp -f dmc/MPI_global_pop_big/dmc.mov1  $(INSTALLDIR)/dmc.mov1.gpopb.mpi
