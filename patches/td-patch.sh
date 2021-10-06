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
#          Iurii Timrov       (THEOS and NCCR-MARVEL, EPFL)
#          Edan Bainglass     (Department of Physics, UNT)
#
#----------------------------------------------------------------------------------------
#
# PATCH script for plugin files and Makefile in TDDFPT/src
#
#----------------------------------------------------------------------------------------

cd $TD_SRC

patch_makefile

check_src_patched
if test "$PATCHED" == 1; then 
   return
else
   message "Patching"
fi

echo "#Please do not remove or modify this file"                          >  Environ_PATCH
echo "#It keeps track of patched versions of the Environ addson package" >> Environ_PATCH
echo "$ENVIRON_VERSION"                                                  >> Environ_PATCH

# The following patches were generated using the command
# diff -u lr_readin.f90.original lr_readin.f90.new > tmp


# plugin_int_forces
if [ -e plugin_tddfpt_potential.f90 ]; then

  sed '/Environ MODULES BEGIN/ a\
      !Environ patch\
      USE scf,                ONLY : rho\
      USE env_global_objects, ONLY : env\
      USE class_calculator,   ONLY : calc\
      !Environ patch
      ' plugin_tddfpt_potential.f90 > tmp.1

  sed '/Environ CALLS BEGIN/ a\
      !Environ patch\
      IF ( use_environ ) THEN\
      !\
      IF (.not.davidson) WRITE( stdout, 8200 )\
8200  FORMAT(5x,"Calculate Environ contribution to response potential")\
      !\
      CALL env%update_response( dfftp%nnr, drho(:,1) )\
      !\
      CALL calc%dpotential( env, dfftp%nnr, dv(:,1) )\
      !\
      END IF\
      !Environ patch
      ' tmp.1 > tmp.2

  mv tmp.2 plugin_tddfpt_potential.f90
  rm tmp.1

fi

printf " done!\n"

cd $QE_DIR
