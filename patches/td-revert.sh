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
# PATCH REVERT for plugin files and Makefile in TDDFPT/src
#
#----------------------------------------------------------------------------------------

cd $TD_SRC

revert_makefile

check_src_reverted
if test "$REVERTED" == 1; then 
   return
else
   message "Reverting"
fi

rm "Environ_PATCH"

# plugin_tddfpt_potential

sed '/Environ patch/,/Environ patch/d' plugin_tddfpt_potential.f90 > tmp.1

mv tmp.1 plugin_tddfpt_potential.f90

printf " done!\n"

cd $QE_DIR
