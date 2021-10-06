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
#          Edan Bainglass     (Department of Physics, UNT)
#
#----------------------------------------------------------------------------------------
#
# PATCH REVERT script for plugin files and Makefile in PW/src
#
#----------------------------------------------------------------------------------------

cd $PW_SRC

revert_makefile

check_src_reverted
if test "$REVERTED" == 1; then 
   return
else
   message "Reverting"
fi

rm "Environ_PATCH"

# plugin_int_forces

sed '/Environ patch/,/Environ patch/d' plugin_int_forces.f90 > tmp.1

mv tmp.1 plugin_int_forces.f90

# plugin_read_input

sed '/Environ patch/,/Environ patch/d' plugin_read_input.f90 > tmp.1

mv tmp.1 plugin_read_input.f90

# plugin_clean

sed '/Environ patch/,/Environ patch/d' plugin_clean.f90 > tmp.1

mv tmp.1 plugin_clean.f90

# plugin_summary

sed '/Environ patch/,/Environ patch/d' plugin_summary.f90 > tmp.1

mv tmp.1 plugin_summary.f90

# plugin_initbase

sed '/Environ patch/,/Environ patch/d' plugin_initbase.f90 > tmp.1

mv tmp.1 plugin_initbase.f90

# plugin_clock

sed '/Environ patch/,/Environ patch/d' plugin_clock.f90 > tmp.1

mv tmp.1 plugin_clock.f90

# plugin_print_energies

sed '/Environ patch/,/Environ patch/d' plugin_print_energies.f90 > tmp.1

mv tmp.1 plugin_print_energies.f90

# plugin_init_ions

sed '/Environ patch/,/Environ patch/d' plugin_init_ions.f90 > tmp.1

mv tmp.1 plugin_init_ions.f90

# plugin_init_cell

sed '/Environ patch/,/Environ patch/d' plugin_init_cell.f90 > tmp.1

mv tmp.1 plugin_init_cell.f90

# plugin_scf_energy

sed '/Environ patch/,/Environ patch/d' plugin_scf_energy.f90 > tmp.1

mv tmp.1 plugin_scf_energy.f90

# plugin_init_potential

sed '/Environ patch/,/Environ patch/d' plugin_init_potential.f90 > tmp.1

mv tmp.1 plugin_init_potential.f90

# plugin_scf_potential

sed '/Environ patch/,/Environ patch/d' plugin_scf_potential.f90 > tmp.1

mv tmp.1 plugin_scf_potential.f90

# plugin_check

sed '/Environ patch/,/Environ patch/d' plugin_check.f90 > tmp.1

mv tmp.1 plugin_check.f90

# makov_payne

sed '/Environ patch/,/Environ patch/d' makov_payne.f90 > tmp.1

mv tmp.1 makov_payne.f90

# force_lc

sed '/Environ patch/,/Environ patch/d' force_lc.f90 > tmp.1

mv tmp.1 force_lc.f90

# plugin_initialization

sed '/Environ patch/,/Environ patch/d' plugin_initialization.f90 > tmp.1

mv tmp.1 plugin_initialization.f90

# plugin_ext_forces

sed '/Environ patch/,/Environ patch/d' plugin_ext_forces.f90 > tmp.1

mv tmp.1 plugin_ext_forces.f90

printf " done!\n"

cd $QE_DIR
