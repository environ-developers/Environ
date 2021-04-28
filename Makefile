#!/bin/bash
#
# Copyright (C) 2018 ENVIRON (www.quantum-environment.org)
#
#    This file is part of Environ version 2.0
#
#    Environ 2.0 is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 2 of the License, or
#    (at your option) any later version.
#
#    Environ 2.0 is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more detail, either the file
#    `License' in the root directory of the present distribution, or
#    online at <http://www.gnu.org/licenses/>.
#
# Author: Oliviero Andreussi (Department of Physics, University of North Texas)
#	      Edan Bainglass (Department of Physics, University of North Texas)
#

ifndef VERBOSE
.SILENT:
endif

ENVIRON_VERSION=2.0
DIRS="PW/src CPV/src XSpectra/src TDDFPT/src"

default: all

all: doc libenviron

doc:
	if test -d Doc ; then (cd Doc; $(MAKE) || exit 1 ); fi

################################################################################
# COMPILATION ROUTINES
################################################################################

# for development purposes
recompile: compile-Environ compile-QE-pw

compile-Environ: check-Environ-makeinc libsdir update-Environ-dependencies
	@ $(MAKE) libfft 2>&1 | tee install/Environ_comp.log
	@ $(MAKE) libutil 2>&1 | tee -a install/Environ_comp.log
	@ $(MAKE) libenv 2>&1 | tee -a install/Environ_comp.log
	@ $(MAKE) check-for-errors prog=Environ

decompile-Environ:
	@ printf "\nCleaning up Environ...\n\n"; $(MAKE) clean

compile-QE-pw: check-QE-makeinc
	@ if test "$(title)"; then \
		  title="$(title)"; \
	  else \
		  title="Compiling QE"; \
	  fi; \
	  printf "\n$$title...\n\n" | tee install/QE_comp.log
	@ (cd ../ && $(MAKE) pw 2>&1 | tee -a Environ/install/QE_comp.log)
	@ $(MAKE) check-for-errors prog=QE

decompile-QE-pw:
	@ printf "\nCleaning up QE...\n\n"
	@ (cd ../ && $(MAKE) clean)

libfft:
	@ printf "\nCompiling FFTXlib...\n\n"
	@ ( \
		cd FFTXlib && $(MAKE) TLDEPS=all || exit 1; \
		mv *.a ../libs \
	 )

libutil:
	@ printf "\nCompiling UtilXlib...\n\n"
	@ ( \
		cd UtilXlib && $(MAKE) TLDEPS=all || exit 1; \
		mv *.a ../libs \
	)

libenv:
	@ printf "\nCompiling Environ/src...\n\n"
	@ ( \
		cd src && $(MAKE) TLDEPS=all || exit 1; \
	   	mv *.a ../libs \
	)

libsdir:
	@ test -d libs || mkdir libs

################################################################################
# CHECKS
################################################################################

check-Environ-makeinc:
	@ if [ ! -e make.inc ]; then \
		  printf "\nMissing make.inc. Please configure installation.\n\n"; \
		  exit 1; \
	  fi
	
check-QE-makeinc:
	@ if [ ! -e ../make.inc ]; then \
		  printf "\nMissing QE/make.inc. Please configure the QE installation.\n\n"; \
		  exit 1; \
	  fi

check-for-errors:
	@ if grep -qE "error #[0-9]+" install/$(prog)_comp.log; then \
		printf "\nErrors found. See install/$(prog)_comp.log\n\n"; \
		exit 1; \
	else \
		printf "\n$(prog) compilation successful!\n\n"; \
		exit; \
	fi

################################################################################
# PATCHING ROUTINES FOR QE+ENVIRON
################################################################################

patch-QE: check-QE-makeinc
	@ printf "\nApplying QE patches using Environ version ${ENVIRON_VERSION}\n"
	@ ./patches/environpatch.sh -patch

revert-QE-patches: check-QE-makeinc
	@ printf "\nReverting QE patches using Environ version ${ENVIRON_VERSION}\n"
	@ ./patches/environpatch.sh -revert

update-Environ-dependencies:
	@ printf "\nUpdating Environ dependencies...\n\n"
	@ ./install/makedeps.sh

update-QE-dependencies:
	@ printf "\nUpdating QE dependencies...\n\n"
	@ (cd ../ && ./install/makedeps.sh "$(DIRS)")

################################################################################
# INSTALL ROUTINES FOR QE+ENVIRON
################################################################################

install-QE+Environ: check-Environ-makeinc check-QE-makeinc
	@ printf "\nPreparing to install QE + Environ $(ENVIRON_VERSION)...\n"
	@ printf "\nDo you wish to proceed (y|n)? "; read c; \
	if [ "$$c" = "y" ]; then \
		printf "\nUse # cores (default = 1) -> "; read cores; \
		printf "\nWould you like to pre-compile QE (y|n)? "; read p; \
		if [ "$$p" = "y" ]; then \
			$(MAKE) -j$${cores:=1} compile-QE-pw title="Pre-compiling QE"; \
			if [ $$? != 0 ]; then exit; fi; \
			(cd install && mv QE_comp.log QE_precomp.log); \
			title="Re-compiling QE with Environ $(ENVIRON_VERSION)"; \
		else \
			printf "\nQE pre-compilation skipped!\n\n"; \
			title="Compiling QE with Environ $(ENVIRON_VERSION)"; \
		fi; \
		$(MAKE) -j$${cores:=1} compile-Environ; \
		if [ $$? != 0 ]; then exit; fi; \
		$(MAKE) patch-QE; \
		$(MAKE) update-QE-dependencies; \
		$(MAKE) -j$${cores:=1} compile-QE-pw title="$$title"; \
		if [ $$? != 0 ]; then exit; fi; \
		if [ $$p = "y" ]; then \
			( \
				cd install && \
				mv QE_comp.log temp; \
				cat QE_precomp.log temp > QE_comp.log; \
				rm temp QE_precomp.log; \
			); \
		fi; \
	else \
		echo; \
	fi

uninstall-QE+Environ:
	@ printf "\nPreparing to uninstall QE + Environ $(ENVIRON_VERSION)...\n"
	@ printf "\nDo you wish to proceed (y|n)? "; read c; \
	if [ "$$c" = "y" ]; then \
		$(MAKE) decompile-Environ; \
		$(MAKE) revert-QE-patches; \
		$(MAKE) update-QE-dependencies; \
		printf "\nPreparing to decompile QE...\n"; \
		printf "\nDo you wish to proceed (y|n)? "; read c; \
		if [ "$$c" = "y" ]; then \
			$(MAKE) decompile-QE-pw; \
			printf "\nDone!\n\n"; \
		else \
			printf "\nQE decompilation skipped!\n\n"; \
		fi; \
	else \
		echo; \
	fi

################################################################################
# CLEANING
################################################################################

clean: check-Environ-makeinc
	@ $(MAKE) clean-src
	@ $(MAKE) clean-libs
	@ $(MAKE) clean-fft
	@ $(MAKE) clean-util

clean-src:
	@ printf "src..........."
	@ (cd src && $(MAKE) clean)
	@ printf " done!\n"

clean-fft:
	@ printf "FFTXlib......."
	@ (cd FFTXlib && $(MAKE) clean)
	@ printf " done!\n"

clean-util:
	@ printf "UtilXlib......"
	@ (cd UtilXlib && $(MAKE) clean)
	@ printf " done!\n"

clean-libs:
	@ printf "libs.........."
	@ if test -d libs; then rm -fr libs; fi
	@ printf " done!\n"

clean-doc:
	@ printf "Docs.........."
	@ (cd Doc && $(MAKE) clean)
	@ printf " done!\n"

# remove files produced by "configure" as well
veryclean: clean
	@ printf "Config........"
	@ (cd install && \
	   rm -rf *.log configure.msg config.status)
	@ rm make.inc
	@ printf " done!\n"

distclean: clean
	@ $(MAKE) clean-doc