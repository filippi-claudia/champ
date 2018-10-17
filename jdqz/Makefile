FC = f77
FFLAGS = -u -O

LIBPATH = -L../jdlib
LIBS = -ljdqr -llapack -lblas
SUBDIRS = jdlib
TEST    = jdtest

all: 	libs jdqr

jdqr: 	
	@for i in $(TEST); do (cd $$i; $(MAKE)); done

libs:
	@for i in $(SUBDIRS); do (cd $$i; $(MAKE)); done

clean:
	-rm -f jdqz
	-rm -f `find . -name '*.[ao]' -print`
