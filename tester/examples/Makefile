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
# Makefile for tester examples
#
#----------------------------------------------------------------------------------------

ifndef VERBOSE
.SILENT:
endif

include ../../make.inc

PROC = 4

EXAMPLES = $(subst /,,$(wildcard */))

TESTER = $(TOPDIR)/tester/tester.x

help:
	@ echo 
	@ echo "    Usage: make run[-parallel] example=[name]"
	@ echo
	@ echo "    available examples = $(EXAMPLES)"
	@ echo

run: check-dir
	@ if test $$? = 0; then \
	    (cd $(example); $(TESTER)); \
	  else \
		exit; \
	  fi

run-parallel: check-dir
	@ if test $$? = 0; then \
	    (cd $(example); mpirun -n $(PROC) $(TESTER)); \
	  else \
		exit; \
	  fi

check-dir:
	@ if ! (test $(example) && test -d $(example)) ; then \
		make help; \
		exit 1; \
	  fi

clean:
	@ for dir in */; do \
	  	(cd $$dir; /bin/rm -f fort.* environ.debug CRASH); \
	  done
