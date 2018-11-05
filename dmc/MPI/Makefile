ifndef FFLAGS

QMCHOME  = ../../..
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
FC      = $(FC_MPI)
FFLAGS := $(FFLAGS_MPI) $(INC)

include $(QMCHOME)/common_objects

OBJVMC = $V/MPI/error_mpi.o $V/strech.o $(OBJCOM) $(OBJCOM1) \
         $V/MPI/optorb_reduce.o $V/MPI/optci_reduce.o $V/MPI/optjas_reduce.o \
         $V/MPI/optx_jas_orb_reduce.o $V/MPI/optx_jas_ci_reduce.o $V/MPI/optx_orb_ci_reduce.o \
         $V/MPI/finwrt_more.o \
         $V/optwf.o $V/optwf_lib.o $V/optjas.o $V/optorb.o $V/optci.o \
         $V/optx_jas_orb.o $V/optx_jas_ci.o $V/optx_orb_ci.o 

OBJDMC = main.o acuest.o dumper.o finwrt.o acues1_reduce.o prop_reduce.o \
         pcm_reduce.o mmpol_reduce.o mc_configs.o \
         ../dmc_good_ps_mov1.o \
         ../multideterminant_tmove.o ../nonloc_grid.o ../optwf_dummy.o \
         ../dmc.o \
x         ../gauss.o ../splitj.o ../walksav_det.o ../walksav_jas.o ../hpsiedmc.o \
         ../prop_dmc.o ../pcm_dmc.o ../mmpol_dmc.o ../optwf_dmc.o \
         ../averages.o \
         $(OBJVMC)

dmc.mov1: $(OBJDMC)
	$(FC) $(FFLAGS) $(OBJDMC) $(CYRUSLIB) $(PSPLINE) $(LAPACK) $(BLAS) $(LDFLAGS) -o $@ 
