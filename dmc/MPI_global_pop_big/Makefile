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
         $V/MPI/finwrt_more.o 

OBJDMCMPI = main.o acuest.o dumper.o finwrt.o redistribute.o \
            move_walker.o walksav_det.o walksav_jas.o prop_mpi.o \
            pcm_mpi.o mmpol_mpi.o optwf_dummy.o \
            acues1_reduce.o ../MPI/mc_configs.o \

OBJDMC = $(OBJDMCMPI) $(OBJVMC) \
         ../dmc_good_ps_mov1.o \
         ../nonloc_grid.o ../multideterminant_tmove.o \
         ../dmc.o ../gauss.o ../splitj.o ../hpsiedmc.o \
         ../prop_dmc.o ../pcm_dmc.o ../mmpol_dmc.o \
         ../averages.o

dmc.mov1: $(OBJDMC)
	$(FC) $(FFLAGS) $(OBJDMC) $(CYRUSLIB) $(PSPLINE) $(LAPACK) $(BLAS) $(MPIETC) $(LDFLAGS) -o $@ 
