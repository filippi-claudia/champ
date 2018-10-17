ifndef FFLAGS

QMCHOME  = ../..

LIBHOME    = $(QMCHOME)
INPUT      = $(QMCHOME)/input

CYRUSLIB   = -L$(LIBHOME)/lib -lcyrus

BLASHOME   = $(LIBHOME)/lib2/blas-2007
LAPACKHOME = $(LIBHOME)/lib2/lapack-3.4.2
BLAS       = -L$(BLASHOME) -lblas
LAPACK     = -L$(LAPACKHOME) -llapack

PSPLINE    = -L$(LIBHOME)/lib2/pspline -lpspline

INSTALLDIR = $(QMCHOME)/bin
INC        = -I$(INPUT) -I$(QMCHOME)/include

QMMM = F
PERIODIC = F

include $(QMCHOME)/settings.make.local

endif

V       = $(QMCHOME)/vmc
FFLAGS := $(FFLAGS) $(INC)

include $(QMCHOME)/common_objects

OBJVMC = $(OBJCOM) $(OBJCOM1) $V/error.o $V/strech.o \
         $V/optwf.o $V/optwf_lib.o $V/optjas.o $V/optorb.o $V/optci.o \
         $V/optx_jas_orb.o $V/optx_jas_ci.o $V/optx_orb_ci.o 

OBJDMC = main.o dmc.o mc_configs.o gauss.o prop_dmc.o optwf_dmc.o pcm_dmc.o mmpol_dmc.o averages.o \
         dmc_good_ps_mov1.o acuest.o \
         finwrt.o dumper.o splitj.o walksav_det.o walksav_jas.o hpsiedmc.o \
         nonloc_grid.o multideterminant_tmove.o optwf_dummy.o \
         acues1_reduce.o \
         $(OBJVMC)

dmc.mov1: $(OBJDMC)
	$(FC) $(FFLAGS) $(OBJDMC) $(CYRUSLIB) $(PSPLINE) $(LAPACK) $(BLAS) -o $@
