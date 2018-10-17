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
M1      = $(QMCHOME)/dmc/MPI
M2      = $(QMCHOME)/dmc/MPI_global_pop_big
FC      = $(FC_MPI)
FFLAGS := $(FFLAGS_MPI) $(INC)

include $(QMCHOME)/common_objects

OBJVMC = $V/MPI/error_mpi.o $V/strech.o $(OBJCOM) $(OBJCOM1) \
         $V/MPI/optorb_reduce.o $V/MPI/optci_reduce.o $V/MPI/optjas_reduce.o \
         $V/MPI/optx_jas_orb_reduce.o $V/MPI/optx_jas_ci_reduce.o $V/MPI/optx_orb_ci_reduce.o \
         $V/MPI/finwrt_more.o 

OBJDMC = main.o walksav_jas.o \
         $(M1)/mc_configs.o \
         $(M2)/acuest.o $(M2)/dumper.o $(M2)/finwrt.o $(M2)/redistribute.o \
         $(M2)/move_walker.o $(M2)/walksav_det.o $(M2)/prop_mpi.o $(M2)/acues1_reduce.o \
         $(M2)/pcm_mpi.o $(M2)/mmpol_mpi.o $(M2)/optwf_dummy.o \
         ../dmc_good_ps_mov1.o \
         ../nonloc_grid.o ../multideterminant_tmove.o \
         ../dmc.o ../gauss.o ../splitj.o ../hpsiedmc.o \
         ../prop_dmc.o ../pcm_dmc.o ../mmpol_dmc.o \
         ../averages.o \
         $(OBJVMC)

dmc.mov1: $(OBJDMC)
	$(FC) $(FFLAGS) $(OBJDMC) $(CYRUSLIB) $(PSPLINE) $(LAPACK) $(BLAS) $(MPIETC) $(LDFLAGS) -o $@ 
