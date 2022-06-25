#----------------------------------------------------------------------------------------
#
# Copyright (C) 2018-2021 ENVIRON (www.quantum-environ.org)
#
#----------------------------------------------------------------------------------------
#
#     This file is part of Environ version 2.0
#     
#     Environ 2.0 is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 2 of the License, or
#     (at your option) any later version.
#     
#     Environ 2.0 is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more detail, either the file
#     `License' in the root directory of the present distribution, or
#     online at <http://www.gnu.org/licenses/>.
#
#----------------------------------------------------------------------------------------
#
# Authors: Oliviero Andreussi (Department of Physics, UNT)
#          Francesco Nattino  (THEOS and NCCR-MARVEL, EPFL)
#          Edan Bainglass     (Department of Physics, UNT)
#
#----------------------------------------------------------------------------------------
#
# Makefile for src
#
#----------------------------------------------------------------------------------------

include ../make.inc

MODFLAGS = $(BASEMOD_FLAGS)

OBJS = $(patsubst %.f90,%.o,$(wildcard *.f90))

all: libsdir libenvsrc.a

libenvsrc.a: $(OBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@
	@ /bin/mv libenvsrc.a ../libs

libsdir:
	@ test -d ../libs || mkdir ../libs

clean:
	@ /bin/rm -f *.o *~ *.a *.d *.i *.L *.mod *.F90 *.x *.i90

include make.depend

# DO NOT DELETE
