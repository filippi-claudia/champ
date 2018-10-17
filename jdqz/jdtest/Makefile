FC = ifort
FFLAGS = -u -O

LIBHOME    = ../../
BLASHOME   = $(LIBHOME)/lib2/blas-3.6.0
LAPACKHOME = $(LIBHOME)/lib2/lapack-3.6.0

LIBPATH = -L../jdlib
LIBS = -ljdqz -L$(LAPACKHOME) -llapack -L$(BLASHOME) -lblas
SUBDIRS = ../jdlib
OBJ = example1.o

all: 	$(OBJ)
	$(FC) $(LDFLAGS) -o example1 $(OBJ) $(LIBPATH) $(LIBS)

clean:
	-rm -f example1
	-rm -f `find . -name '*.[ao]' -print`
