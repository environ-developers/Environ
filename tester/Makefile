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

ifndef VERBOSE
.SILENT:
endif

include ../make.inc

MODFLAGS = $(BASEMOD_FLAGS)

all: tester.x

tester.x: tester.o
	@ $(LD) $(LDFLAGS) -o $@ tester.o ../libs/* $(ENVLIBS)

tester.o : ../utils/parallel_include.o
tester.o : ../utils/io.o
tester.o : ../src/constants.o
tester.o : ../src/setup.o
tester.o : ../src/environment.o
tester.o : ../src/calculator.o
tester.o : ../src/destructor.o
tester.o : ../src/input.o
tester.o : ../src/base_input.o

clean:
	@ /bin/rm -f *.o *~ *.a *.d *.i *.L *.mod *.F90 *.x fort.* environ.debug CRASH
