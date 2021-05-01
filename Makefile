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
DIRS="PW/src XSpectra/src TDDFPT/src" # TODO add CPV/src after PW/src when CP is fixed

################################################################################
# HELP/LIST
################################################################################

# quick manual
help:
	@ echo
	@ echo 'QE + Environ $(ENVIRON_VERSION) installaion:'
	@ echo
	@ echo 'From QE root:'
	@ echo
	@ echo '  1) ./configure [--prefix=]'
	@ echo
	@ echo 'From Environ root:'
	@ echo
	@ echo '  1) ./configure [--prefix=]'
	@ echo
	@ echo '  2) make install-QE+Environ'
	@ echo
	@ echo '     a) select a QE package to install:'
	@ echo
	@ echo '        - pw           basic code for scf, structure optimization, MD'
	@ # echo '	       - cp           CP code: Car-Parrinello molecular dynamics' # TODO add CP option when fixed
	@ echo '        - tddfpt       time dependent dft code'
	@ echo '        - xspectra     X-ray core-hole spectroscopy calculations'
	@ echo '        - all          same as "make pw tddfpt xspectra"' # TODO add cp after pw when CP is fixed
	@ echo
	@ echo '     b) select number of cores (for parallel compilation)'
	@ echo '        - if issues arise in parallel, run in serial (core=1)'
	@ echo
	@ echo '     c) (optional) pre-compile QE without Environ (for the more cautious)'
	@ echo
	@ echo 'To Uninstall:'
	@ echo
	@ echo '  1) make uninstall-QE+Environ'
	@ echo '     - QE decompilation step is optional (see documentaion for more detail)'
	@ echo

# list make targets (documented)
list:
	@ grep -B 1 "^[^.]\S*:" Makefile | cut -d ':' -f 1 

################################################################################
# COMPILATION ROUTINES
################################################################################

# compile Environ (updates dependencies by default)
compile-Environ: check-Environ-makeinc libsdir update-Environ-dependencies
	@ $(MAKE) libfft 2>&1 | tee install/Environ_comp.log
	@ $(MAKE) libutil 2>&1 | tee -a install/Environ_comp.log
	@ $(MAKE) libenv 2>&1 | tee -a install/Environ_comp.log
	@ $(MAKE) check-for-errors in=Environ

# decompile Environ
decompile-Environ:
	@ printf "\nCleaning up Environ...\n\n"; $(MAKE) clean

# compile QE (make compile-QE [prog=<pw|tddfpt|xspectra>])
compile-QE: check-QE-makeinc
	@ if test "$(prog)"; then prog="$(prog)"; else prog=pw; fi
	@ if test "$(title)"; then title="$(title)"; else title="Compiling QE"; fi; \
	  printf "\n$$title...\n\n" | tee install/QE_comp.log
	@ (cd ../ && $(MAKE) $$prog 2>&1 | tee -a Environ/install/QE_comp.log)
	@ $(MAKE) check-for-errors in=QE

# decompile QE
decompile-QE:
	@ printf "\nCleaning up QE...\n\n"
	@ (cd ../ && $(MAKE) clean)

libfft :
	@ printf "\nCompiling FFTXlib...\n\n"
	@ ( \
		cd FFTXlib && $(MAKE) TLDEPS=all || exit 1; \
		mv *.a ../libs \
	 )

libutil :
	@ printf "\nCompiling UtilXlib...\n\n"
	@ ( \
		cd UtilXlib && $(MAKE) TLDEPS=all || exit 1; \
		mv *.a ../libs \
	)

libenv :
	@ printf "\nCompiling Environ/src...\n\n"
	@ ( \
		cd src && $(MAKE) TLDEPS=all || exit 1; \
	   	mv *.a ../libs \
	)

libsdir :
	@ test -d libs || mkdir libs

# compile Environ documentation
compile-doc:
	if test -d Doc ; then (cd Doc; $(MAKE) TLDEPS=all || exit 1 ); fi

################################################################################
# CHECKS
################################################################################

check-Environ-makeinc :
	@ if [ ! -e make.inc ]; then \
		  printf "\nMissing make.inc. Please configure installation.\n\n"; \
		  exit 1; \
	  fi
	
check-QE-makeinc :
	@ if [ ! -e ../make.inc ]; then \
		  printf "\nMissing QE/make.inc. Please configure the QE installation.\n\n"; \
		  exit 1; \
	  fi

check-for-errors :
	@ if grep -qE "error #[0-9]+" install/$(in)_comp.log; then \
		printf "\nErrors found. See install/$(in)_comp.log\n\n"; \
		exit 1; \
	else \
		printf "\n$(in) compilation successful! \n\n"; \
		exit; \
	fi

################################################################################
# PATCHING ROUTINES FOR QE+ENVIRON
################################################################################

# patch QE
patch-QE: check-QE-makeinc
	@ printf "\nApplying QE patches using Environ version ${ENVIRON_VERSION}\n"
	@ ./patches/environpatch.sh -patch

# revert QE patches
revert-QE-patches: check-QE-makeinc
	@ printf "\nReverting QE patches using Environ version ${ENVIRON_VERSION}\n"
	@ ./patches/environpatch.sh -revert

# update Environ dependencies
update-Environ-dependencies:
	@ printf "\nUpdating Environ dependencies...\n\n"
	@ ./install/makedeps.sh

# update QE dependencies
update-QE-dependencies:
	@ printf "\nUpdating QE dependencies...\n\n"
	@ (cd ../ && ./install/makedeps.sh "$(DIRS)")

################################################################################
# INSTALL ROUTINES FOR QE+ENVIRON
################################################################################

print_menu :
	@ printf "\nSelect a package:\n\n"
	@ printf "%s\n%s\n%s\n%s\n\n%s" \
			 "   1 - PW" \
			 "   2 - TD" \
			 "   3 - XS" \
			 "   4 - ALL" \
			 "-> "

# TODO add CP option when fixed
# compile Environ, patch QE, update QE dependencies, and compile selected QE package(s)
install-QE+Environ: check-Environ-makeinc check-QE-makeinc
	@ printf "\nPreparing to install QE + Environ $(ENVIRON_VERSION)...\n"
	@ make print_menu; read c; \
	\
	case $$c in \
	1) opt=pw;; \
	2) opt=tddfpt;; \
	3) opt=xspectra;; \
	4) opt="pw tddfpt xspectra";; \
	*) exit;; \
	esac; \
	\
	printf "\nUse # cores (default = 1) -> "; read cores; \
	printf "\nWould you like to pre-compile QE (y|n)? "; read p; \
	\
	if [ "$$p" = "y" ]; then \
		$(MAKE) -j$${cores:=1} compile-QE prog="$$opt" title="Pre-compiling QE"; \
		if [ $$? != 0 ]; then exit; fi; \
		(cd install && mv QE_comp.log QE_precomp.log); \
		title="Re-compiling QE with Environ $(ENVIRON_VERSION)"; \
	else \
		printf "\nQE pre-compilation skipped! \n\n"; \
		title="Compiling QE with Environ $(ENVIRON_VERSION)"; \
	fi; \
	\
	$(MAKE) -j$${cores:=1} compile-Environ; \
	if [ $$? != 0 ]; then exit; fi; \
	$(MAKE) patch-QE; \
	$(MAKE) update-QE-dependencies; \
	$(MAKE) -j$${cores:=1} compile-QE prog="$$opt" title="$$title"; \
	if [ $$? != 0 ]; then exit; fi; \
	\
	if [ $$p = "y" ]; then \
		( \
			cd install && \
			mv QE_comp.log temp; \
			cat QE_precomp.log temp > QE_comp.log; \
			rm temp QE_precomp.log; \
		); \
	fi

# decompile Environ, revert QE patches, update QE dependencies, and (optional) decompile QE
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
			$(MAKE) decompile-QE; \
			printf "\nDone! \n\n"; \
		else \
			printf "\nQE decompilation skipped! \n\n"; \
		fi; \
	else \
		echo; \
	fi

################################################################################
# CLEANING
################################################################################

# remove executables and objects
clean: check-Environ-makeinc
	@ $(MAKE) clean-src
	@ $(MAKE) clean-libs
	@ $(MAKE) clean-fft
	@ $(MAKE) clean-util

clean-src :
	@ printf "src..........."
	@ (cd src && $(MAKE) clean)
	@ printf " done! \n"

clean-fft :
	@ printf "FFTXlib......."
	@ (cd FFTXlib && $(MAKE) clean)
	@ printf " done! \n"

clean-util :
	@ printf "UtilXlib......"
	@ (cd UtilXlib && $(MAKE) clean)
	@ printf " done! \n"

# remove compiled libraries
clean-libs :
	@ printf "libs.........."
	@ if test -d libs; then rm -fr libs; fi
	@ printf " done! \n"

# remove compiled documentation files
clean-doc:
	@ printf "Docs.........."
	@ (cd Doc && $(MAKE) clean)
	@ printf " done! \n"

# also remove configuration files
veryclean: clean
	@ printf "Config........"
	@ (cd install && \
	   rm -rf *.log configure.msg config.status)
	@ rm make.inc
	@ printf " done! \n"
