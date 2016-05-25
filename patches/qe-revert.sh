#!/bin/bash

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

# makov_payne

sed '/Environ patch/,/Environ patch/d' makov_payne.f90 > tmp.1

mv tmp.1 makov_payne.f90
