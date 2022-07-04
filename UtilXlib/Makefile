#----------------------------------------------------------------------------------------
#
# Copyright (C) 2018-2022 ENVIRON (www.quantum-environ.org)
# Copyright (C) Quantum ESPRESSO group
#
#----------------------------------------------------------------------------------------
#
#     This file is part of Environ version 3.0
#
#     Environ 3.0 is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 2 of the License, or
#     (at your option) any later version.
#
#     Environ 3.0 is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more detail, either the file
#     `License' in the root directory of the present distribution, or
#     online at <http://www.gnu.org/licenses/>.
#
#----------------------------------------------------------------------------------------
#
# Authors: Modified by Edan Bainglass
#
#----------------------------------------------------------------------------------------
#
# Makefile for UtilXlib
#
#----------------------------------------------------------------------------------------

include ../make.inc

OBJS = $(patsubst %.f90,%.o,$(wildcard *.f90))

C_OBJS = cptimer.o

OBJS += $(C_OBJS)

all: libsdir libenvutil.a

libenvutil.a: $(OBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@
	@ /bin/mv libenvutil.a ../libs

libsdir:
	@ test -d ../libs || mkdir ../libs

clean:
	@ /bin/rm -f *.o *.a *.d *.i *~ *_tmp.f90 *.mod *.L *.x *.i90

include make.depend
