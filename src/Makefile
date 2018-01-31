# Copyright (C) 2018 ENVIRON (www.quantum-environment.org)
#
#    This file is part of Environ version 1.0
#
#    Environ 1.0 is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 2 of the License, or
#    (at your option) any later version.
#
#    Environ 1.0 is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more detail, either the file
#    `License' in the root directory of the present distribution, or
#    online at <http://www.gnu.org/licenses/>.
#
# Makefile for ENVIRON module
#
# Authors: Oliviero Andreussi (Department of Physics, University of North Thexas)
#          Francesco Nattino  (THEOS and NCCR-MARVEL, Ecole Polytechnique Federale de Lausanne)
#

include ../../make.sys

# location of needed modules
MODFLAGS= $(MOD_FLAG)../../iotk/src $(MOD_FLAG)../../Modules $(MOD_FLAG).

ENVIRONOBJS = \
boundary.o                  \
cavity.o                    \
charges.o                   \
dielectric.o                \
electrostatic.o             \
electrostatic_init.o        \
electrostatic_mod.o         \
electrostatic_types.o       \
environ_fd_gradient.o       \
environ_generate_function.o \
environ_init.o              \
environ_input.o             \
environ_main.o              \
environ_mod.o               \
environ_output.o            \
environ_types.o             \
extcharge.o                 \
functions.o                 \
generalized.o               \
generate_boundary.o         \
ioncc.o                     \
ions.o                      \
linearized_pb.o             \
periodic.o                  \
poisson.o                   \
pressure.o                  \
solvent_tddfpt.o            \
util.o

all : libenviron.a

libenviron.a: $(ENVIRONOBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

clean :
	- /bin/rm -f *.o *~ *.a *.d *.i *.L *.mod *.F90

include make.depend

# DO NOT DELETE
