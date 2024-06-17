#!/bin/bash
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

ENVIRON_VERSION=3.0

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
	@ echo "  - compiles Environ's UtilXlib, FFTXlib, and src"
	@ echo "  - compilation generates Environ/install/Environ_comp.log"
	@ echo
	@ echo "Updating Dependencies"
	@ echo "---------------------"
	@ echo
	@ echo "* depend"
	@ echo
	@ echo "  - updates dependencies in Environ's UtilXlib, FFTXlib, and src"
	@ echo
	@ echo "* cannibalize-QE"
	@ echo
	@ echo "  - updates UtilXlib and FFTXlib to currently checked-out QE branch"
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
	  $(MAKE) compile-programs; \
	  ( \
		cd install; \
		cat UtilXlib_comp.log FFTXlib_comp.log src_comp.log > Environ_comp.log; \
		/bin/rm UtilXlib_comp.log FFTXlib_comp.log src_comp.log \
	)
	@ printf "\nEnviron $(ENVIRON_VERSION) compilation successful! \n"

compile-util: check-makeinc
	@ printf "\nCompiling UtilXlib...\n\n" 2>&1 | \
	tee install/UtilXlib_comp.log
	@ ( cd UtilXlib && $(MAKE) all || exit 1 ) 2>&1 | tee -a install/UtilXlib_comp.log
	@ $(MAKE) check-for-errors prog=UtilXlib

compile-fft: check-makeinc
	@ printf "\nCompiling FFTXlib...\n\n" 2>&1 | \
	tee install/FFTXlib_comp.log
	@ ( cd FFTXlib && $(MAKE) all || exit 1 ) 2>&1 | tee -a install/FFTXlib_comp.log
	@ $(MAKE) check-for-errors prog=FFTXlib

compile-src: check-makeinc
	@ printf "\nCompiling src...\n\n" 2>&1 | \
	tee install/src_comp.log
	@ ( cd src && $(MAKE) all || exit 1 ) 2>&1 | tee -a install/src_comp.log
	@ $(MAKE) check-for-errors prog=src

compile-programs: check-makeinc
	@ printf "\nCompiling programs...\n\n" 2>&1 | \
        tee install/programs_comp.log
	@ ( cd programs && $(MAKE) || exit 1 ) 2>&1 | tee -a install/programs_comp.log
	@ $(MAKE) check-for-errors prog=programs

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
	@ ./install/makedeps.sh; echo

cannibalize-QE:
	@ if test -d $(qedir) && test $(version); then \
		./update_libs.sh -l fft -p $(qedir) -m -v $(version); \
		./update_libs.sh -l util -p $(qedir) -m -v $(version); \
	  fi

################################################################################
# CLEANING
################################################################################

# remove executables and objects
clean: check-makeinc
	@ $(MAKE) clean-libs
	@ $(MAKE) clean-bin
	@ $(MAKE) clean-logs
	@ $(MAKE) clean-util
	@ $(MAKE) clean-fft
	@ $(MAKE) clean-src
	@ $(MAKE) clean-programs

clean-bin:
	@ printf "bin..........."
	@ if test -d bin; then /bin/rm -fr bin; fi
	@ printf " done! \n"

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
	@ printf "UtilXlib......"
	@ (cd UtilXlib && $(MAKE) clean)
	@ printf " done! \n"

clean-fft:
	@ printf "FFTXlib......."
	@ (cd FFTXlib && $(MAKE) clean)
	@ printf " done! \n"

clean-src:
	@ printf "src..........."
	@ (cd src && $(MAKE) clean)
	@ printf " done! \n"

clean-programs:
	@@ printf "programs......"
	@ (cd programs && $(MAKE) clean)
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
