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
# Authors: Edan Bainglass (Department of Physics, UNT)
#
#----------------------------------------------------------------------------------------
#
# Makefile for tester
#
#----------------------------------------------------------------------------------------

ifndef VERBOSE
.SILENT:
endif

include ../make.inc

MODFLAGS = $(BASEMOD_FLAGS)

all: tester.x

OBJS = utilities.o tester.o

tester.x: tester.o
	@ $(LD) $(LDFLAGS) -o $@ $(OBJS) ../libs/* $(ENVLIBS)

clean:
	@ /bin/rm -f *.o *~ *.a *.d *.i *.L *.mod *.F90 *.x fort.* environ.debug CRASH

tester.o : utilities.o
