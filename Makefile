#!/bin/bash
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
#          Edan Bainglass     (Department of Physics, UNT)
#
#----------------------------------------------------------------------------------------
#
# Makefile for Environ
#
#----------------------------------------------------------------------------------------

ifndef VERBOSE
.SILENT:
endif

.PHONY: install

ENVIRON_VERSION=2.0

################################################################################
# HELP/LIST
################################################################################

default:
	@ echo
	@ echo "- make help -> installation manual"
	@ echo
	@ echo "- make devs -> developer options"
	@ echo

# quick manual
help:
	@ echo
	@ sed -n '/^# QE + Environ/,/^---/p' README.md
	@ echo

# rundown of developer make routines
devs:
	@ echo
	@ echo "Installation"
	@ echo "------------"
	@ echo
	@ echo "* see 'make help'"
	@ echo
	@ echo "Compilation"
	@ echo "-----------"
	@ echo
	@ echo "* compile (requires Environ/make.inc)"
	@ echo
	@ echo "  - compiles Environ's utils, FFTs, and src"
	@ echo "  - compilation generates Environ/install/Environ_comp.log"
	@ echo
	@ echo "Updating Dependencies"
	@ echo "---------------------"
	@ echo
	@ echo "* depend"
	@ echo
	@ echo "  - updates dependencies in Environ's utils, FFTs, and src"
	@ echo
	@ echo "Cleaning"
	@ echo "--------"
	@ echo
	@ echo "* clean (requires Environ/make.inc)"
	@ echo
	@ echo "  - remove Environ objects, libraries, and compilation logs"
	@ echo
	@ echo "* veryclean (requires Environ/make.inc)"
	@ echo
	@ echo "  - 'make clean' + removes configuration files"
	@ echo

################################################################################
# COMPILATION ROUTINES
################################################################################

compile: check-makeinc
	@ printf "\nCompiling Environ $(ENVIRON_VERSION)...\n"; \
	  $(MAKE) compile-util; \
	  $(MAKE) compile-fft; \
	  $(MAKE) compile-src; \
	  ( \
		cd install; \
		cat utils_comp.log FFTs_comp.log src_comp.log > Environ_comp.log; \
		/bin/rm utils_comp.log FFTs_comp.log src_comp.log \
	  )
	  printf "\nEnviron $(ENVIRON_VERSION) compilation successful! \n"

compile-util: check-makeinc
	@ printf "\nCompiling utils...\n\n" 2>&1 | \
	tee install/utils_comp.log
	@ ( cd utils && $(MAKE) all || exit 1 ) 2>&1 | tee -a install/utils_comp.log
	@ $(MAKE) check-for-errors prog=utils

compile-fft: check-makeinc
	@ printf "\nCompiling FFTs...\n\n" 2>&1 | \
	tee install/FFTs_comp.log
	@ ( cd FFTs && $(MAKE) all || exit 1 ) 2>&1 | tee -a install/FFTs_comp.log
	@ $(MAKE) check-for-errors prog=FFTs
	 
compile-src: check-makeinc
	@ printf "\nCompiling src...\n\n" 2>&1 | \
	tee install/src_comp.log
	@ ( cd src && $(MAKE) all || exit 1 ) 2>&1 | tee -a install/src_comp.log
	@ $(MAKE) check-for-errors prog=src

compile-doc:
	@ if test -d Doc ; then (cd Doc; $(MAKE) TLDEPS=all || exit 1 ); fi

################################################################################
# CHECKS
################################################################################

check-makeinc:
	@ if [ ! -e make.inc ]; then \
		  printf "\nMissing make.inc. Please configure installation.\n\n"; \
		  exit 1; \
	  fi

check-for-errors:
	@ if grep -qE "^make.*Error [0-9]+" install/$(prog)_comp.log; then \
		  printf "\nErrors found. See install/$(prog)_comp.log\n\n"; \
		  exit 1; \
	  else \
		  printf "\n$(prog) compilation successful! \n"; \
		  exit; \
	  fi

################################################################################
# DEPENDENCY UPDATING ROUTINE
################################################################################

depend:
	@ printf "\nUpdating Environ dependencies...\n\n"
	@ ./install/makedeps.sh

################################################################################
# CLEANING
################################################################################

# remove executables and objects
clean: check-makeinc
	@ $(MAKE) clean-libs
	@ $(MAKE) clean-logs
	@ $(MAKE) clean-util
	@ $(MAKE) clean-fft
	@ $(MAKE) clean-src

clean-libs:
	@ printf "libs.........."
	@ if test -d libs; then /bin/rm -fr libs; fi
	@ printf " done! \n"

clean-logs:
	@ printf "logs.........."
	@ ( \
		cd install && \
		find . -type f -name '*log' -not -name config.log | xargs /bin/rm -f \
	)
	@ printf " done! \n"

clean-util:
	@ printf "utils........."
	@ (cd utils && $(MAKE) clean)
	@ printf " done! \n"

clean-fft:
	@ printf "FFTs.........."
	@ (cd FFTs && $(MAKE) clean)
	@ printf " done! \n"

clean-src:
	@ printf "src..........."
	@ (cd src && $(MAKE) clean)
	@ printf " done! \n"

clean-doc:
	@ printf "Docs.........."
	@ (cd Doc && $(MAKE) clean)
	@ printf " done! \n"

# also remove configuration files
veryclean: clean
	@ printf "config........"
	@ (cd install && \
	   /bin/rm -rf *.log configure.msg config.status)
	@ /bin/rm make.inc
	@ printf " done! \n"
