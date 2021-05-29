# Makefile for FFTs

include ../make.inc

DFTI = $(findstring DFTI,$(DFLAGS))

MODFLAGS = $(BASEMOD_FLAGS)

OBJS = $(patsubst %.f90,%.o,$(wildcard *.f90))

C_OBJS = fft_stick.o fftw.o fftw_sp.o fftw_dp.o

OBJS += $(C_OBJS)

all: driver

driver:
	@ if test $(DFTI); then \
		  make enable-mkl-include; \
		  $(MAKE) libfft.a; \
		  make disable-mkl-include; \
		  $(MAKE) libfft.a; \
	  else \
	      $(MAKE) libfft.a; \
	  fi

libfft.a: $(OBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

##################################################################
# The following two targets are used to get around the inclusion 
# of mkl_dfti.f90 in fft_scalar.DFTI.f90 by both Environ and QE.
# The targets are called from the root Makefile by compile-fft.
#
# To avoid multiple definitions of DFTI items, we enable the 
# include in Environ for the initial compilation, then disable 
# it for a subsequent recompilation of the relevant files. This 
# effectively removes the mkl items from libfft, letting the DFTI
# items reference QE's definitions. 
#
# A better solution is more than welcomed :)
##################################################################

enable-mkl-include:
	@ file=fft_scalar.DFTI.f90; \
	  pattern='^! #include "mkl'; \
	  sed -i.tmp -e "/$$pattern/s/! //" $$file && rm $$file.tmp

disable-mkl-include:
	@ file=fft_scalar.DFTI.f90; \
	  pattern='^#include "mkl'; \
	  sed -i.tmp -e "/$$pattern/s/#/! #/" $$file && rm $$file.tmp
	  
##################################################################

clean :
	@ /bin/rm -f *.o *.a *.d *.i *~ *_tmp.f90 *.mod *.L 

include make.depend
