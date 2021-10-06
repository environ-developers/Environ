#!/bin/bash
#----------------------------------------------------------------------------------------
#
# Copyright (C) 2018-2021 ENVIRON (www.quantum-environ.org)
# Copyright (C) 2001-2017 Quantum ESPRESSO (www.quantum-espresso.org)
#
#----------------------------------------------------------------------------------------
#
# This file is part of Environ version 2.0
#
# Environ 2.0 is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# Environ 2.0 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more detail, either the file
# `License' in the root directory of the present distribution, or
# online at <http://www.gnu.org/licenses/>.
#
#----------------------------------------------------------------------------------------
#
# Authors: Oliviero Andreussi (Department of Physics, UNT)
#          Edan Bainglass     (Department of Physics, UNT)
#
#----------------------------------------------------------------------------------------
#
# This file was adapted from the following Quantum ESPRESSO v6.1 file:
#
# 		QE/test-suite/Makefile
#
#----------------------------------------------------------------------------------------

ifndef VERBOSE
.SILENT:
endif

include ./ENVIRONMENT

default :
	@echo 'To run Environ tests, type at the shell prompt:'
	@echo ' '
	@echo '  make target'
	@echo ' '
	@echo 'where <target> identifies an action'
	@echo ' '
	@echo '  run-tests                            : run-tests-serial (default)'
	@echo '  run-tests-serial                     : run all tests in SERIAL' 
	@echo '  run-tests-parallel                   : run all tests in PARALLEL (4 MPI)'
	@echo '  run-custom-test testdir=DIR          : run test in DIR only (serial)'
	@echo '  run-custom-test-parallel testdir=DIR : run test in DIR only (parallel 4 MPI)'
	@echo '  compare                              : compare last output with reference'
	@echo '  clean                                : clean stdout/sderr of all tests'
	@echo ' '
	@echo 'Any target can take an optional dblcell=true to trigger the doublecell test suite'
	@echo '- make clean by default cleans all tests including doublecell'
	@echo '- see Makefile for finer cleaning options'
	@echo ' '
	@echo 'For additional advanced commands and settings, please manually inspect'
	@echo 'ENVIRONMENT and Makefile files'

prolog : check-version
	@sed "s|XXXXXX|$(ENVIRON_ROOT)|g" < userconfig.tmp > userconfig

pseudo :
	@./check_pseudo.sh pw_

check-version:
	@ v=$$(grep -E 'New in.*[0-9]\.[0-9].*' $(ESPRESSO_ROOT)/Doc/release-notes | grep -oE '[0-9]\.[0-9]' | head -1); \
	  if [ $$(echo "$$v < 6.6" | bc) -eq 1 ]; then \
	  	  echo; \
		  echo "##       ##    ###     #######   ###   ##  ####  ###   ##   ###### "; \
		  echo "##   #   ##   ## ##    ##   ###  ## #  ##   ##   ## #  ##  ##      "; \
		  echo " ## ### ##   ## # ##   #######   ##  # ##   ##   ##  # ##  ##  ### "; \
		  echo "  ### ###   ##     ##  ##    ##  ##   ###  ####  ##   ###   #### ##"; \
		  echo; \
		  echo "* For QE < 6.6, some tests may fail with negligible errors."; \
		  echo; \
		  echo "  This is due to the following change in QE 6.6:"; \
		  echo; \
		  echo '      Default value of parameter "diago_david_ndim" changed to 2'; \
		  echo; \
		  echo "  (https://github.com/QEF/q-e/blob/master/Doc/release-notes)"; \
		  echo; \
	  fi

run-tests : run-tests-serial

run-tests-serial : run-tests-pw-serial

run-tests-parallel : run-tests-pw-parallel

run-tests-pw-serial : prolog pseudo
	@ if [ ${dblcell} ]; then if [ ${dblcell} = true ]; then cat=doublecell_; else cat=; fi; fi; \
	  env QE_USE_MPI=0 ${TESTCODE_DIR}/bin/testcode.py --verbose --category=$${cat}pw_all

run-tests-pw-parallel : prolog pseudo
	@ if [ ${dblcell} ]; then if [ ${dblcell} = true ]; then cat=doublecell_; else cat=; fi; fi; \
	  env QE_USE_MPI=1 ${TESTCODE_DIR}/bin/testcode.py --verbose --category=$${cat}pw_all

run-custom-test : run-custom-test-serial

run-custom-test-serial : prolog pseudo
	@ if [ ${dblcell} ]; then if [ ${dblcell} = true ]; then cat=doublecell/; else cat=; fi; fi; \
	  if [ $(testdir) ]; then \
		  if test -d $(testdir); then \
	    	  env QE_USE_MPI=0 ${TESTCODE_DIR}/bin/testcode.py --verbose --category=$${cat}$(testdir); \
		  else \
		      echo "${testdir} not found"; \
		  fi; \
	  else \
		  echo "Missing testdir argument"; \
	  fi

run-custom-test-parallel : prolog pseudo
	@ if [ ${dblcell} ]; then if [ ${dblcell} = true ]; then cat=doublecell/; else cat=; fi; fi; \
	  if [ $(testdir) ]; then \
		  if test -d $(testdir); then \
	    	  env QE_USE_MPI=1 ${TESTCODE_DIR}/bin/testcode.py --verbose --category=$${cat}$(testdir); \
		  else \
		      echo "${testdir} not found"; \
		  fi; \
	  else \
		  echo "Missing testdir argument"; \
	  fi

create-reference-pw : pseudo prolog
	@ if [ ${dblcell} ]; then if [ ${dblcell} = true ]; then cat=doublecell_; else cat=; fi; fi; \
	  env QE_USE_MPI=0 ${TESTCODE_DIR}/bin/testcode.py --category=$${cat}pw_all make-benchmarks

create-custom-reference-pw : pseudo prolog
	@ if [ ${dblcell} ]; then if [ ${dblcell} = true ]; then cat=doublecell/; else cat=; fi; fi; \
	  if [ $(testdir) ]; then \
		  if test -d $(testdir); then \
	    	  env QE_USE_MPI=0 ${TESTCODE_DIR}/bin/testcode.py --category=$${cat}$(testdir) make-benchmarks; \
		  else \
		      echo "${testdir} not found"; \
		  fi; \
	  else \
		  echo "Missing testdir argument"; \
	  fi

compare : compare-pw

compare-pw : prolog
	@ if [ ${dblcell} ]; then if [ ${dblcell} = true ]; then cat=doublecell_; else cat=; fi; fi; \
	  ${TESTCODE_DIR}/bin/testcode.py --verbose --category=$${cat}pw_all compare 2>&1 | \
	  tee out.PW.`date +%Y%m%d_%H%M%S`

time-diff :
	@ if test -d $(testdir); then ./time_diff.sh ${testdir}; fi

clean : clean-pw-tests clean-dblcell-tests

clean-pw-tests :
	@ for x in `find pw_* -name "test*"`; do rm -rf $$x; done
	@ for x in `find pw_* -name "CRASH"`; do rm -rf $$x; done
	@ for x in `find pw_* -name "environ.*"`; do rm -rf $$x; done
	@ for x in `find pw_* -name "*.cube"`; do rm -rf $$x; done
	@ for x in `find pw_* -name "fort.*"`; do rm -rf $$x; done

clean-dblcell-tests:
	@ for x in `find doublecell/pw_* -name "test*"`; do rm -rf $$x; done
	@ for x in `find doublecell/pw_* -name "CRASH"`; do rm -rf $$x; done
	@ for x in `find doublecell/pw_* -name "environ.*"`; do rm -rf $$x; done
	@ for x in `find doublecell/pw_* -name "*.cube"`; do rm -rf $$x; done
	@ for x in `find doublecell/pw_* -name "fort.*"`; do rm -rf $$x; done

# Calling 'purge' means remove all reference outputs... be careful!
purge: clean
	@ for x in `find pw_* -name "benchmark*"`; do rm -rf $$x; done
	@ for x in `find doublecell/pw_* -name "benchmark*"`; do rm -rf $$x; done
