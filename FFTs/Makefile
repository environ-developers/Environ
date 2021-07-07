#----------------------------------------------------------------------------------------
#
# Copyright (C) 2021 ENVIRON (www.quantum-environ.org)
# Copyright (C) Quantum ESPRESSO (www.quantum-espresso.org)
#
#----------------------------------------------------------------------------------------
#
#	  This file is part of Environ version 2.0
#	  
#	  Environ 2.0 is free software: you can redistribute it and/or modify
#	  it under the terms of the GNU General Public License as published by
#	  the Free Software Foundation, either version 2 of the License, or
#	  (at your option) any later version.
#	  
#	  Environ 2.0 is distributed in the hope that it will be useful,
#	  but WITHOUT ANY WARRANTY; without even the implied warranty of
#	  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#	  GNU General Public License for more detail, either the file
#	  `License' in the root directory of the present distribution, or
#	  online at <http://www.gnu.org/licenses/>.
#
#----------------------------------------------------------------------------------------
#
# Authors: Modified by Edan Bainglass
#
#----------------------------------------------------------------------------------------
#
# Makefile for FFTs
#
#----------------------------------------------------------------------------------------

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
