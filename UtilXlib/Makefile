# Makefile for UtilXlib

include ../make.inc

OBJS = $(patsubst %.f90,%.o,$(wildcard *.f90))

C_OBJS = cptimer.o

OBJS += $(C_OBJS)

all : libutil.a

libutil.a: $(OBJS)
	$(AR) $(ARFLAGS) $@ $?       
	$(RANLIB) $@    

clean :
	@ /bin/rm -f *.o *.a *.d *.i *~ *_tmp.f90 *.mod *.L *.x

include make.depend
