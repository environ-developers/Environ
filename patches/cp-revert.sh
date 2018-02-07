# Copyright (C) 2018 ENVIRON (www.quantum-environment.org)
#
#    This file is part of Environ version 1.0
#
#    Environ 1.0 is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 2 of the License, or
#    (at your option) any later version.
#
#    Environ 1.0 is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more detail, either the file
#    `License' in the root directory of the present distribution, or
#    online at <http://www.gnu.org/licenses/>.
#
# PATCH REVERT script for plugin files in CPV/src
#
#!/bin/bash

QEDIR="$PWD"
cd $QEDIR/CPV/src

# plugin_int_forces

sed '/Environ patch/,/Environ patch/d' plugin_int_forces.f90 > tmp.1

mv tmp.1 plugin_int_forces.f90

# plugin_read_input

sed '/Environ patch/,/Environ patch/d' plugin_read_input.f90 > tmp.1

mv tmp.1 plugin_read_input.f90

# plugin_clean

sed '/Environ patch/,/Environ patch/d' plugin_clean.f90 > tmp.1

mv tmp.1 plugin_clean.f90

# plugin_print_info

sed '/Environ patch/,/Environ patch/d' plugin_print_info.f90 > tmp.1

mv tmp.1 plugin_print_info.f90

# plugin_init_base

sed '/Environ patch/,/Environ patch/d' plugin_init_base.f90 > tmp.1

mv tmp.1 plugin_init_base.f90

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

# plugin_energy

sed '/Environ patch/,/Environ patch/d' plugin_energy.f90 > tmp.1

mv tmp.1 plugin_energy.f90

# plugin_add_potential

sed '/Environ patch/,/Environ patch/d' plugin_add_potential.f90 > tmp.1

mv tmp.1 plugin_add_potential.f90

# plugin_get_potential

sed '/Environ patch/,/Environ patch/d' plugin_get_potential.f90 > tmp.1

mv tmp.1 plugin_get_potential.f90

# plugin_utilities

sed '/Environ patch/,/Environ patch/d' plugin_utilities.f90 > tmp.1

mv tmp.1 plugin_utilities.f90

cd $QEDIR
