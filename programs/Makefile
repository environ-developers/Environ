#----------------------------------------------------------------------------------------
#
# Copyright (C) 2018-2022 ENVIRON (www.quantum-environ.org)
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

DRIVER_OBJS = cmdline_args.o parsers.o prog_utils.o programs.o driver.o

ENVIRON_LIBS = -L$(TOPDIR)/libs -lenvsrc -lenvfft -lenvutil

driver: driver.o
	@ $(LD) $(LDFLAGS) -o $@ $(DRIVER_OBJS) $(ENVIRON_LIBS) $(ENVLIBS)
	@ test -d ../bin || mkdir ../bin
	@ ( cd ../bin ; ln -fs ../programs/driver driver)

help:
	@ echo
	@ echo \'make driver\' \(default\) to compile the driver \(requires Environ compilation\)
	@ echo
	@ echo \'make avail\' for available calculations
	@ echo
	@ echo \'make usage name=name-of-calculation\' for specific execution instructions
	@ echo
	@ echo \'make clean\' to clean directory
	@ echo

avail:
	@ echo
	@ echo "Available calculations:"
	@ echo
	@ echo "tester"
	@ echo
	@ echo "  For testing custom features"
	@ echo
	@ echo "from_cube"
	@ echo
	@ echo "  Computes environment contributions to properties of a \"frozen\" density given as a cube file"
	@ echo
	@ echo "with_aims"
	@ echo
	@ echo "  For use with FHI-aims"
	@ echo
	@ echo "descriptors"
	@ echo
	@ echo "  An Environ calculation of the Soft-Sphere interface function surface, volume, and local surface"
	@ echo "  and volume. Additionally, the forces due to the interface function and energy contributions can"
	@ echo "  be calculated"
	@ echo
	@ echo "  For fine calculation control, see the README file or check out https://environ.readthedocs.io"
	@ echo

usage:
	@ echo
	@ case "$(name)" in \
	  	from_cube) usage='driver -n from_cube -i input-file -c cube-file';; \
	        with_aims) usage='<mpiexe> -np <N> driver -n with_aims -i input-file & <mpiexe> -np <N> <aims_binary>';; \
	  	descriptors) usage='driver -n descriptors -i input-file -c cube-file';; \
		*) usage='make usage name=name-of-calculation';; \
	  esac; \
	  echo Usage: $$usage
	@ echo

clean:
	@ /bin/rm -f *.o *~ *.a *.d *.i *.L *.mod *.F90 tester driver fort.* environ.debug CRASH

driver.o: cmdline_args.o programs.o
programs.o: cmdline_args.o prog_utils.o
prog_utils.o: cmdline_args.o parsers.o
parsers.o: cmdline_args.o
