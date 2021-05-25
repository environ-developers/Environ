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
	@ echo "* compile-Environ (requires Environ/make.inc)"
	@ echo
	@ echo "  - compiles Environ's utils, FFTs, and src"
	@ echo "  - compilation generates Environ/install/Environ_comp.log"
	@ echo
	@ echo "  * NOTE: Environ is decoupled from QE. Changes to QE files"
	@ echo "          (exluding Environ-patched sections) do not require"
	@ echo "          Environ recompilation"
	@ echo
	@ echo "* compile-QE [prog=] (requires QE/make.inc)"
	@ echo
	@ echo "  - (re)compiles QE package -> prog = pw (default) | tddfpt | xspectra"
	@ echo "  - compilation generates Environ/install/QE_comp.log"
	@ echo
	@ echo "  * NOTE: If QE dependencies reflect an Environ-patched state,"
	@ echo "          Environ must be pre-compiled before QE is (re)compiled."
	@ echo
	@ echo "* recompile-QE+Environ"
	@ echo
	@ echo "  - wrapper on compile-Environ and compile-QE"
	@ echo "  - used when changes are made to Environ"
	@ echo "  - updates Environ dependencies before compilation"
	@ echo
	@ echo "  * DO NOT use this if changing the patching scheme!"
	@ echo
	@ echo "Patching & Dependencies"
	@ echo "-----------------------"
	@ echo
	@ echo "* patch-QE [prog=] (requires QE/make.inc)"
	@ echo
	@ echo "  - applies Environ patch to QE/install/makedeps.sh"
	@ echo "  - patches plugin files and Makefile of QE/<prog>"
	@ echo "  - prog = pw | tddfpt | xspectra | all (default)"
	@ echo "  - patches Makefiles of all pw.x-dependent QE packages"
	@ echo
	@ echo "* revert-QE-patches (requires QE/make.inc)"
	@ echo
	@ echo "  - reverts patches applied during patch-QE"
	@ echo
	@ echo "* update-Environ-dependencies"
	@ echo
	@ echo "  - updates dependencies in Environ's utils, FFTs, and src"
	@ echo
	@ echo "* update-QE-dependencies [prog=]"
	@ echo
	@ echo "  - updates dependencies in QE package (prog)"
	@ echo "  - prog = pw | tddfpt | xspectra | all (default)"
	@ echo
	@ echo "Cleaning"
	@ echo "--------"
	@ echo
	@ echo "* decompile-Environ"
	@ echo
	@ echo "  - wrapper for 'make clean-Environ'"
	@ echo
	@ echo "* decompile-QE (requires QE/make.inc)"
	@ echo
	@ echo "  - wrapper for QE/Makefile -> 'make clean'"
	@ echo
	@ echo "* clean-Environ (requires Environ/make.inc)"
	@ echo
	@ echo "  - remove Environ objects, libraries, and compilation logs"
	@ echo
	@ echo "* veryclean (requires Environ/make.inc)"
	@ echo
	@ echo "  - 'make clean-Environ' + removes configuration files"
	@ echo

################################################################################
# COMPILATION ROUTINES
################################################################################

compile-Environ: check-Environ-makeinc libsdir
	@ printf "\nCompiling Environ $(ENVIRON_VERSION)...\n\n"
	@ $(MAKE) compile-util
	@ $(MAKE) compile-fft
	@ $(MAKE) compile-src
	@ ( \
		cd install; \
		cat utils_comp.log FFTs_comp.log src_comp.log > Environ_comp.log; \
		rm utils_comp.log FFTs_comp.log src_comp.log \
	)
	@ printf "\nEnviron $(ENVIRON_VERSION) compilation successful! \n\n"

compile-QE: check-QE-makeinc
	@ if test "$(prog)"; then prog="$(prog)"; else prog=pw; fi; \
	  if [ $$prog = all ]; then prog="pw tddfpt xspectra"; fi; \
	  if test "$(title)"; then title="$(title)"; else title="Compiling QE"; fi; \
	  printf "\n$$title...\n\n" | tee install/QE_comp.log; \
	  (cd ../ && $(MAKE) $$prog 2>&1 | tee -a Environ/install/QE_comp.log)
	@ $(MAKE) check-for-errors prog=QE

# used after changes made to Environ
recompile-QE+Environ:
	@ make print_menu; read c; \
	\
	case $$c in \
	1) opt=pw;; \
	2) opt=tddfpt;; \
	3) opt=xspectra;; \
	4) opt=all;; \
	*) exit;; \
	esac; \
	\
	printf "\nUse # cores (default = 1) -> "; read cores; \
	make update-Environ-dependencies; \
	$(MAKE) compile-Environ; \
	if [ $$? -ne 0 ]; then exit; fi; \
	$(MAKE) compile-QE prog=$$opt

compile-util: libsdir
	@ printf "\nCompiling utils...\n\n" 2>&1 | \
	tee install/utils_comp.log
	@ ( \
		cd utils && $(MAKE) all || exit 1; \
		mv *.a ../libs \
	) 2>&1 | tee -a install/utils_comp.log
	@ $(MAKE) check-for-errors prog=utils

compile-fft: libsdir
	@ printf "\nCompiling FFTs...\n\n" 2>&1 | \
	tee install/FFTs_comp.log
	@ ( \
		cd FFTs && \
		$(MAKE) enable-mkl-include; \
		$(MAKE) all || exit 1; \
		$(MAKE) disable-mkl-include; \
		$(MAKE) all || exit 1; \
		mv *.a ../libs \
	) 2>&1 | tee -a install/FFTs_comp.log
	@ $(MAKE) check-for-errors prog=FFTs
	 
compile-src: libsdir
	@ printf "\nCompiling src...\n\n" 2>&1 | \
	tee install/src_comp.log
	@ ( \
		cd src && $(MAKE) all || exit 1; \
	   	mv *.a ../libs \
	) 2>&1 | tee -a install/src_comp.log
	@ $(MAKE) check-for-errors prog=src

libsdir:
	@ test -d libs || mkdir libs

compile-doc:
	@ if test -d Doc ; then (cd Doc; $(MAKE) TLDEPS=all || exit 1 ); fi

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
	@ if grep -qE "[Ee]rror #?[0-9]+" install/$(prog)_comp.log; then \
		  printf "\nErrors found. See install/$(prog)_comp.log\n\n"; \
		  exit 1; \
	  else \
		  printf "\n$(prog) compilation successful! \n\n"; \
		  exit; \
	  fi

################################################################################
# PATCHING ROUTINES FOR QE+ENVIRON
################################################################################

patch-QE: check-QE-makeinc
	@ printf "\nApplying QE patches using Environ version ${ENVIRON_VERSION}\n"
	@ if test $(prog); then prog=$(prog); else prog=all; fi; \
	  ./patches/environpatch.sh -patch $$prog

revert-QE-patches: check-QE-makeinc
	@ printf "\nReverting QE patches using Environ version ${ENVIRON_VERSION}\n"
	  ./patches/environpatch.sh -revert

update-Environ-dependencies:
	@ printf "\nUpdating Environ dependencies...\n\n"
	@ ./install/makedeps.sh

# TODO add CP option when fixed
update-QE-dependencies:
	@ printf "\nUpdating QE dependencies...\n\n"
	@ if test $(prog); then prog=$(prog); else prog=all; fi; \
	  case $$prog in \
	  pw) DIRS="PW/src";; \
	  tddfpt) DIRS="PW/src TDDFPT/src";; \
	  xspectra) DIRS="PW/src XSpectra/src";; \
	  *) DIRS="PW/src TDDFPT/src XSpectra/src";; \
	  esac; \
	  (cd ../ && ./install/makedeps.sh $$DIRS)

################################################################################
# INSTALL ROUTINES FOR QE+ENVIRON
################################################################################

print_menu:
	@ printf "\nSelect a package:\n\n"
	@ printf "%s\n%s\n%s\n%s\n\n%s" \
			 "   1 - PW" \
			 "   2 - TDDFPT" \
			 "   3 - XSpectra" \
			 "   4 - ALL" \
			 "-> "

# TODO add CP option when fixed
install-QE+Environ: check-Environ-makeinc check-QE-makeinc
	@ printf "\nPreparing to install QE + Environ $(ENVIRON_VERSION)...\n"
	@ make print_menu; read c; \
	\
	case $$c in \
	1) opt=pw;; \
	2) opt=tddfpt;; \
	3) opt=xspectra;; \
	4) opt=all;; \
	*) exit;; \
	esac; \
	\
	printf "\nUse # cores (default = 1) -> "; read cores; \
	\
	make -j$${cores:=1} compile-QE prog=$$opt title="Pre-compiling QE"; \
	if [ $$? -ne 0 ]; then exit; fi; \
	(cd install && mv QE_comp.log QE_precomp.log); \
	\
	make -j$${cores:=1} compile-Environ; \
	if [ $$? -ne 0 ]; then exit; fi; \
	make patch-QE prog=$$opt; \
	make update-QE-dependencies prog=$$opt; \
	\
	title="Re-compiling QE with Environ $(ENVIRON_VERSION)"; \
	make -j$${cores:=1} compile-QE prog=$$opt title="$$title"; \
	if [ $$? -ne 0 ]; then exit; fi; \
	\
	printf "\nQE + Environ $(ENVIRON_VERSION) installaion complete! \n\n"
	\
	( \
		cd install && \
		mv QE_comp.log temp; \
		cat QE_precomp.log temp > QE_comp.log; \
		rm temp QE_precomp.log; \
	)

uninstall-QE+Environ:
	@ printf "\nPreparing to uninstall QE + Environ $(ENVIRON_VERSION)...\n"
	@ printf "\nDo you wish to proceed (y|n)? "; read c; \
	if [ "$$c" = y ]; then \
		make decompile-Environ; \
		make revert-QE-patches; \
		make update-QE-dependencies; \
		printf "\nPreparing to decompile QE...\n"; \
		printf "\nDo you wish to proceed (y|n)? "; read c; \
		if [ "$$c" = y ]; then \
			make decompile-QE; \
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

decompile-Environ:
	@ printf "\nCleaning up Environ...\n\n"; $(MAKE) clean-Environ

decompile-QE: check-QE-makeinc
	@ printf "\nCleaning up QE...\n\n"
	@ (cd ../ && $(MAKE) clean)

# dummy routine called by QE
clean:

# remove executables and objects
clean-Environ: check-Environ-makeinc
	@ $(MAKE) clean-libs
	@ $(MAKE) clean-logs
	@ $(MAKE) clean-util
	@ $(MAKE) clean-fft
	@ $(MAKE) clean-src

clean-libs:
	@ printf "libs.........."
	@ if test -d libs; then rm -fr libs; fi
	@ printf " done! \n"

clean-logs:
	@ printf "logs.........."
	@ ( \
		cd install && \
		find . -type f -name '*log' -not -name config.log | xargs rm -f \
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
veryclean: clean-Environ
	@ printf "config........"
	@ (cd install && \
	   rm -rf *.log configure.msg config.status)
	@ rm make.inc
	@ printf " done! \n"
