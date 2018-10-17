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
FC      = $(FC_MPI)
FFLAGS := $(FFLAGS_MPI) $(INC)

include $(QMCHOME)/common_objects

OBJVMC = main.o dumper.o mc_configs.o \
         acuest_reduce.o fin_reduce.o finwrt_more.o \
         prop_reduce.o pcm_reduce.o mmpol_reduce.o\
         optorb_reduce.o optci_reduce.o optjas_reduce.o \
         optx_jas_orb_reduce.o optx_jas_ci_reduce.o optx_orb_ci_reduce.o \
         sr_more.o regterg.o jdqz_main.o lin_d_more.o dl_more.o \
         force_analy_reduce.o error_mpi.o \
         $V/vmc.o $V/acuest.o $V/acuest_write.o $V/dumper_more.o \
         $V/finwrt.o $V/multiple_states.o $V/determinant_psig.o \
         $V/prop_vmc.o $V/strech.o $V/pcm_vmc.o $V/mmpol_vmc.o \
         $V/optwf.o $V/optwf_vmc.o $V/optwf_lib.o \
         $V/optjas.o $V/optorb.o $V/optci.o \
         $V/optx_jas_orb.o $V/optx_orb_ci.o $V/optx_jas_ci.o \
         $V/sr.o $V/lin_d.o $V/jdqz_lib.o $V/dl.o $V/dl_iter.o $V/fetch_parameters.o

OVMC1  = $(OBJVMC) $(OBJCOM) $(OBJCOM1) $V/metrop_mov1_slat.o 

vmc.mov1: $(OVMC1)
	$(FC) $(FFLAGS) $(OVMC1) $(CYRUSLIB) $(LAPACK) $(PSPLINE) $(BLAS) -o $@ $(LDFLAGS)
