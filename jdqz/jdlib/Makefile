#
#     Coded by Diederik R. Fokkema
#
#     $Id$
# 
#

export SHELL=/bin/sh
export FC = ifort
# export FFLAGS = -u -O
export FFLAGS = -O0 -CB

LIB = libjdqz.a

OBJS = $(patsubst %.f,%.o,$(wildcard *.f))

$(LIB): $(LIB)($(OBJS)) 
	ranlib $@

$(LIB)(jdqz.o): jdqz.f
	$(FC) -O0 -CB -c jdqz.f
	$(AR) rv $(LIB) jdqz.o
	-rm jdqz.o

# $(FC) -u -O1 -c jdqz.f
