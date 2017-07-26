#!/bin/bash
#
# Copyright (C) 2017 ENVIRON (www.quantum-environment.org)
# Copyright (C) 2001-2016 Quantum ESPRESSO
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License. See the file `License' in the root directory
# of the present distribution.

#source ${ENVIRON_ROOT}/tests/ENVIRONMENT

if [ $QE_USE_MPI == 1 ]; then
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
