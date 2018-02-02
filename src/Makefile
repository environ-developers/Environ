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
correction_periodic.o      \
electrostatic_init.o       \
electrostatic_mod.o        \
electrostatic_types.o      \
embedding_electrostatic.o  \
embedding_surface.o        \
embedding_volume.o         \
environ_init.o             \
environ_input.o            \
environ_main.o             \
environ_mod.o              \
environ_output.o           \
environ_types.o            \
problem_generalized.o      \
problem_linearized_pb.o    \
problem_poisson.o          \
solvent_tddfpt.o           \
tools_fd_gradient.o        \
tools_generate_boundary.o  \
tools_generate_functions.o \
utils_boundary.o           \
utils_charges.o            \
utils_dielectric.o         \
utils_electrolyte.o        \
utils_externals.o          \
utils_functions.o          \
utils_ions.o

all : libenviron.a

libenviron.a: $(ENVIRONOBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

clean :
	- /bin/rm -f *.o *~ *.a *.d *.i *.L *.mod *.F90

include make.depend

# DO NOT DELETE
