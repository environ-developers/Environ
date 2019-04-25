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
#      QE/test-suite/run-pw.sh
#
# Author: Oliviero Andreussi (Department of Physics, University of North Thexas)
#
#!/bin/bash

#source ${ENVIRON_ROOT}/tests/ENVIRONMENT

if [ $QE_USE_MPI = 1 ]; then
  export PARA_PREFIX="mpirun -np ${TESTCODE_NPROCS}"
else
  unset PARA_PREFIX
fi

# Additional stuff before run special test-cases

${PARA_PREFIX} ${ESPRESSO_ROOT}/bin/pw.x --environ "$@"

if [ -f environ.debug ]; then
  cat environ.debug
  rm environ.debug
fi

rm -f input_tmp.in environ.in
