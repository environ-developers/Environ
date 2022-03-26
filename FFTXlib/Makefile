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
# Makefile for FFTXlib
#
#----------------------------------------------------------------------------------------

include ../make.inc

DFTI = $(findstring DFTI,$(DFLAGS))

MODFLAGS = $(BASEMOD_FLAGS)

OBJS = $(patsubst %.f90,%.o,$(filter-out mkl_wrapper.f90,$(wildcard *.f90)))

all: libsdir libfft.a

libfft.a: $(OBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@
	@ /bin/mv libfft.a ../libs

fft_scalar.DFTI.o: mkl_wrapper.o

libsdir:
	@ test -d ../libs || mkdir ../libs

clean:
	@ /bin/rm -f *.o *.a *.d *.i *~ *_tmp.f90 *.mod *.L

include make.depend
