# Copyright (C) 2018 ENVIRON (www.quantum-environment.org)
#
#    This file is part of Environ version 1.1
#
#    Environ 1.1 is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 2 of the License, or
#    (at your option) any later version.
#
#    Environ 1.1 is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more detail, either the file
#    `License' in the root directory of the present distribution, or
#    online at <http://www.gnu.org/licenses/>.
#
# Author: Oliviero Andreussi (Department of Physics, University of North Thexas)
#

ifndef VERBOSE
.SILENT:
endif

default: all

all: doc libenviron

doc:
	if test -d Doc ; then (cd Doc; $(MAKE) || exit 1 ); fi

################################################################################
# COMPILATION ROUTINES
################################################################################

compile-environ: check-environ-makeinc libsdir
	@ $(MAKE) libfft
	@ $(MAKE) libutil
	@ $(MAKE) libenv

decompile-environ: 
	@ echo "\nCleaning up Environ...\n"; $(MAKE) clean

compile-qe-pw: check-qe-makeinc
	@ echo "\nCompiling QE...\n"
	@ (cd ../ && $(MAKE) pw)

decompile-qe-pw:
	@ echo "\nCleaning up QE...\n"
	@ (cd ../ && $(MAKE) clean)

libfft:
	@ echo "\nCompiling FFTXlib...\n"
	@ ( \
		cd FFTXlib ; $(MAKE) TLDEPS= all || exit 1; \
		mv *.a ../libs \
	 )

libutil: 
	@ echo "\nCompiling UtilXlib...\n"
	@ ( \
		cd UtilXlib ; $(MAKE) TLDEPS= all || exit 1; \
		mv *.a ../libs \
	)

libenv:
	@ echo "\nCompiling Environ/src...\n"
	@ if test -d src ; then \
          (cd src && \
		  if test "$(MAKE)" = ""; then \
		  	  $(MAKE) $(MFLAGS) $@; \
          else \
			  $(MAKE) $(MFLAGS); \
		  fi; \
		  mv *.a ../libs \
		  ); \
	  fi

libsdir:
	@ test -d libs || mkdir libs

check-environ-makeinc: # TODO can the error be silent?
	@ if [ ! -e make.inc ]; then \
		  echo "\nMissing make.inc. Please configure installation.\n"; \
		  exit 1; \
	  fi
	
check-qe-makeinc: # TODO can the error be silent?
	@ if [ ! -e ../make.inc ]; then \
		  echo "\nMissing QE/make.inc. Please configure the QE installation.\n"; \
		  exit 1; \
	  fi

################################################################################
# INSTALL ROUTINES FOR QE+ENVIRON
################################################################################

patch-qe: check-qe-makeinc
	@ echo "\nApplying QE patches...\n"
	@ ./patches/environpatch.sh -patch
	@ $(MAKE) update-QE-dependencies

revert-qe-patches: check-qe-makeinc
	@ echo "\nReverting QE patches...\n"
	@ ./patches/environpatch.sh -revert
	@ $(MAKE) update-QE-dependencies

update-QE-dependencies:
	@ echo "\nUpdating QE dependencies...\n"
	@ (cd ../ && ./install/makedeps.sh)

install-QE+Environ: check-environ-makeinc check-qe-makeinc
	@ echo "\nThis will compile Environ, patch QE, then compile QE.\n"
	@ echo -n "Do you wish to proceed (y|n)? "; read c; \
	if [ "$$c" = "y" ]; then \
		echo -n "\nUse # cores (default = 1) -> "; read cores; \
		$(MAKE) -j$${cores:=1} compile-environ; \
		$(MAKE) -j$${cores:=1} patch-qe; \
		$(MAKE) -j$${cores:=1} compile-qe-pw; \
		echo "\nDone!\n"; \
	else \
		echo; \
	fi

uninstall-QE+Environ: 
	@ echo "\nThis will decompile Environ, revert QE patches, and decompile QE.\n"
	@ echo -n "Do you wish to proceed (y|n)? "; read c; \
	if [ "$$c" = "y" ]; then \
		$(MAKE) decompile-environ; \
		$(MAKE) revert-qe-patches; \
		$(MAKE) decompile-qe-pw; \
		echo "\nDone!\n"; \
	else \
		echo; \
	fi

################################################################################
# CLEANING
################################################################################

clean:
	@ touch make.inc
	@ $(MAKE) clean-src
	@ $(MAKE) clean-libs
	@ $(MAKE) clean-fft
	@ $(MAKE) clean-util

clean-src:
	@ echo -n "src..........."
	@ (cd src && $(MAKE) clean)
	@ echo -n "done!\n"

clean-fft:
	@ echo -n "FFTXlib......."
	@ (cd FFTXlib && $(MAKE) clean)
	@ echo -n "done!\n"

clean-util:
	@ echo -n "UtilXlib......"
	@ (cd UtilXlib && $(MAKE) clean)
	@ echo -n "done!\n"

clean-libs:
	@ echo -n "libs.........."
	@ if test -d libs; then rm -fr libs; fi
	@ echo -n "done!\n"

clean-doc:
	@ echo -n "Docs.........."
	@ (cd Doc && $(MAKE) clean)
	@ echo -n "done!\n"

# remove files produced by "configure" as well
veryclean: clean
	@ echo -n "Config........"
	@ (cd install && \
	   rm -rf config.log configure.msg config.status)
	@ rm make.inc
	@ echo -n "done!\n"

distclean: clean 
	@ $(MAKE) clean-doc