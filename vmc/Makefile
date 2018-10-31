ifndef FFLAGS

QMCHOME  = ..
LIBHOME    = $(QMCHOME)
INPUT      = $(QMCHOME)/input

CYRUSLIB   = -L$(LIBHOME)/lib -lcyrus

BLASHOME   = $(LIBHOME)/lib2/blas-2007
LAPACKHOME = $(LIBHOME)/lib2/lapack-3.6.0
BLAS       = -L$(BLASHOME) -lblas
LAPACK     = -L$(LAPACKHOME) -llapack

PSPLINE    = -L$(LIBHOME)/lib2/pspline -lpspline

INSTALLDIR = $(QMCHOME)/bin
INC        = -I$(INPUT) -I$(QMCHOME)/include

QMMM = F
PERIODIC = F

include $(QMCHOME)/settings.make.local
endif

FFLAGS := -ggdb $(FFLAGS) $(INC)
F90FLAGS := -ggdb -fcheck=all $(F90FLAGS) $(INC)

include $(QMCHOME)/common_objects

OBJVMC = main.o vmc.o \
         acuest.o acuest_write.o acuest_reduce.o fin_reduce.o \
         mc_configs.o finwrt.o finwrt_more.o dumper.o dumper_more.o \
         multiple_states.o determinant_psig.o \
         prop_vmc.o strech.o error.o pcm_vmc.o mmpol_vmc.o \
         optjas.o optorb.o optci.o optx_jas_orb.o optx_jas_ci.o optx_orb_ci.o \
         optwf.o optwf_vmc.o optwf_lib.o \
	 coords_int.o \
         sr.o sr_more.o \
         lin_d.o lin_d_more.o regterg.o jdqz_main.o jdqz_lib.o \
         dl.o dl_more.o dl_iter.o fetch_parameters.o

OVMC1  = $(OBJCOM) $(OBJVMC) $(OBJCOM1) metrop_mov1_slat.o 

commands.p2: $V/read_input.f 
	$(INPUT)/inpt.awk  -v p2file=$@ -v logf=$(INPUT)/commands.log $^

$V/p2prog.f: commands.p2 $(INPUT)/end.p2 $(INPUT)/auto.p2
	cat commands.p2 $(INPUT)/end.p2 $(INPUT)/auto.p2  > all.p2
	$(INPUT)/genp2_defaults.awk  -v ofile=$@ all.p2 > genp2.log

$(INPUT)/p2nmcheck.f: namelist.p2
	$(INPUT)/mkp2nmcheck.pl $^ > $@

vmc.mov1: $(OVMC1)
	$(FC) $(FFLAGS) $(OVMC1) $(CYRUSLIB) $(LINPACK) $(LAPACK) $(PSPLINE) $(BLAS) -o $@
