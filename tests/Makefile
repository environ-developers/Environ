# Copyright (C) 2018 ENVIRON (www.quantum-environment.org)
# Copyright (C) 2001-2017 Quantum ESPRESSO group
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
# This file was adapted from the following Quantum ESPRESSO v6.1 file
#
#      QE/test-suite/Makefile
#
# Author: Oliviero Andreussi (Department of Physics, University of North Thexas)
#
#!/bin/bash

include ./ENVIRONMENT

default :
	@echo 'To run Environ tests, type at the shell prompt:'
	@echo ' '
	@echo '  make target'
	@echo ' '
	@echo 'where <target> identifies an action'
	@echo ' run-tests                   : run-tests-serial (default)'
	@echo ' run-tests-serial            : run all tests in SERIAL' 
	@echo ' run-tests-parallel          : run all tests in PARALLEL (4 MPI)'
	@echo ' run-custom-test testdir=DIR : run test in DIR only (serial)'
	@echo ' run-custom-test-parallel testdir=DIR : run test in DIR only (parallel 4 MPI)'
	@echo ' compare                     : compare last output with reference'
	@echo ' clean                       : clean stdout/sderr of all tests'
	@echo ' '
	@echo 'For additional advanced commands and settings please manually inspect'
	@echo 'ENVIRONMENT and Makefile files'

prolog :
	@sed "s|XXXXXX|$(ENVIRON_ROOT)|g" < userconfig.tmp > userconfig

pseudo :
	@./check_pseudo.sh pw_


run-tests : run-tests-serial

run-tests-serial : run-tests-pw-serial

run-tests-parallel : run-tests-pw-parallel

run-tests-pw-serial : prolog pseudo
	env QE_USE_MPI=0 ${TESTCODE_DIR}/bin/testcode.py --verbose --category=pw_all

run-tests-pw-parallel : prolog pseudo
	env QE_USE_MPI=1 ${TESTCODE_DIR}/bin/testcode.py --verbose --category=pw_all

run-custom-test : run-custom-test-serial

run-custom-test-serial : prolog pseudo
	@if test -d $(testdir); then \
	env QE_USE_MPI=0 ${TESTCODE_DIR}/bin/testcode.py --verbose --category=$(testdir) ; fi

run-custom-test-parallel : prolog pseudo
	@if test -d $(testdir); then \
	env QE_USE_MPI=1 ${TESTCODE_DIR}/bin/testcode.py --verbose --category=$(testdir) ; fi

create-reference-pw : pseudo prolog
	env QE_USE_MPI=0 ${TESTCODE_DIR}/bin/testcode.py --category=pw_all make-benchmarks

create-custom-reference-pw : pseudo prolog
	@if test -d $(testdir); then \
	env QE_USE_MPI=0 ${TESTCODE_DIR}/bin/testcode.py --category=$(testdir) make-benchmarks ; fi

compare : compare-pw

compare-pw : prolog
	${TESTCODE_DIR}/bin/testcode.py --category=pw_all  --verbose compare 2>&1 | tee out.PW.`date +%Y%m%d_%H%M%S`

clean:
	@for x in `find pw_* -name "test*"`; do rm -rf $$x; done
	@for x in `find pw_* -name "CRASH"`; do rm -rf $$x; done
#This command is removed because pw_ms-gcs directory has environ.in
#as input file.
#@for x in `find pw_* -name "environ.*"`; do rm -rf $$x; done

# Calling 'purge' means remove all reference outputs... be careful!
purge: clean
	@for x in `find pw_* -name "benchmark*"`; do rm -rf $$x; done
