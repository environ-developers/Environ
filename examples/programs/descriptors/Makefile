#----------------------------------------------------------------------------------------
#
# Copyright (C) 2018-2021 ENVIRON (www.quantum-environ.org)
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
# Makefile for tester examples
#
#----------------------------------------------------------------------------------------

ifndef VERBOSE
.SILENT:
endif

include ../../../make.inc

PROC = 4

EXAMPLES = $(subst .cube,,$(wildcard *.cube))

CMD = $(TOPDIR)/programs/driver -n descriptors -min 1 -max 10 -step 1

help:
	@ echo 
	@ echo "    usage: make run[-parallel] example=[name]"
	@ echo
	@ echo "    available examples = $(EXAMPLES)"
	@ echo

run:
	@ if test $(example); then \
		if test -f $(example).cube; then \
		  $(CMD) -c $(example).cube; \
		else \
		  echo; \
		  echo "    $(example) is not a valid example"; \
		  echo; \
		  echo "    available examples: $(EXAMPLES)"; \
		  echo; \
		fi; \
	  else \
	  	echo; \
	    echo "    missing example argument"; \
	  	echo; \
	  fi

run-parallel:
	@ if test $$? = 0; then \
	    (mpirun -n $(PROC) $(CMD) -c $(example).cube); \
	  else \
		exit; \
	  fi

clean:
	@ /bin/rm -f fort.* environ.debug CRASH
