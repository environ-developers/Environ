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
#          Francesco Nattino  (THEOS and NCCR-MARVEL, EPFL)
#          Ismaila Dabo       (DMSE, Penn State)
#          Edan Bainglass     (Department of Physics, UNT)
#
#----------------------------------------------------------------------------------------
#
# PATCH script for plugin files and Makefile in PW/src
#
#----------------------------------------------------------------------------------------

cd $PW_SRC

patch_makefile

check_src_patched
if test "$PATCHED" == 1; then 
   return
else
   message "Patching"
fi

echo "#Please do not remove or modify this file" >Environ_PATCH
echo "#It keeps track of patched versions of the Environ addson package" >>Environ_PATCH
echo "$ENVIRON_VERSION" >>Environ_PATCH

# plugin_int_forces

sed '/Environ MODULES BEGIN/ a\
!Environ patch\
  USE environ_QE_interface,  ONLY : calc_environ_force\
!Environ patch
' plugin_int_forces.f90 >tmp.1

sed '/Environ VARIABLES BEGIN/ a\
!Environ patch \
  REAL(DP), ALLOCATABLE :: force_environ(:,:)\
!Environ patch
' tmp.1 >tmp.2

sed '/Environ CALLS BEGIN/ a\
!Environ patch\
  IF (use_environ) THEN\
    !\
    ALLOCATE(force_environ(3,nat))\
    !\
    force_environ=0.0_dp\
    !\
    ! ... Add environment contributions\
    !\
    CALL calc_environ_force( nat, force_environ )\
    !\
    IF ( iverbosity > 0 ) THEN\
      WRITE( stdout, 9001 )\
      DO na = 1, nat\
         WRITE( stdout, 9002 ) na, ityp(na), ( force_environ(ipol,na), ipol = 1, 3 )\
      END DO\
      WRITE( stdout, * )\
    ENDIF\
    !\
    force = force_environ\
    !\
    DEALLOCATE(force_environ)\
    !\
  END IF\
  !\
9001 FORMAT(5x,"The global environment contribution to forces")\
9002 FORMAT(5X,"atom ",I4," type ",I2,"   force = ",3F14.8)\
!Environ patch
' tmp.2 >tmp.1

mv tmp.1 plugin_int_forces.f90

# plugin_read_input

sed '/Environ MODULES BEGIN/ a\
!Environ patch\
  USE io_global,            ONLY : ionode, ionode_id, stdout\
  USE mp_images,            ONLY : intra_image_comm\
  USE ions_base,            ONLY : nat, ntyp => nsp, atm\
  USE martyna_tuckerman,    ONLY : do_comp_mt\
  USE environ_QE_interface, ONLY : init_environ_io, init_environ_base_first\
!Environ patch
' plugin_read_input.f90 >tmp.1

sed '/Environ CALLS BEGIN/ a\
!Environ patch\
   IF (use_environ) THEN\
      CALL init_environ_io(prog, ionode, ionode_id, intra_image_comm, stdout)\
      CALL init_environ_base_first(1, nat, ntyp, atm, do_comp_mt)\
   ENDIF\
!Environ patch
' tmp.1 >tmp.2

mv tmp.2 plugin_read_input.f90

# plugin_clean

sed '/Environ MODULES BEGIN/ a\
!Environ patch\
USE environ_QE_interface, ONLY : environ_clean, environ_clean_first, & \
                                 environ_clean_second, is_tddfpt\
!Environ patch
' plugin_clean.f90 >tmp.1

sed '/Environ CALLS BEGIN/ a\
!Environ patch\
   IF (use_environ) THEN\
      !\
      IF (prog(1:2) == "PW") THEN\
         !\
         ! When called by PW, but inside a TD calculation\
         ! do not clean environ variables, they have been\
         ! already cleaned by TD. The lflag input is used\
         ! to fully clean the variable or to only clean\
         ! variables initialized during the PW run and not the\
         ! ones initialized while processing the input:\
         ! this allows NEB simulations\
         !\
         IF (.NOT. is_tddfpt()) CALL environ_clean(lflag)\
         !\
      ELSE IF ( prog(1:2) == "TD" ) THEN\
         !\
         ! When called by TD, use the flag input variable to\
         ! specify whether to clean the PW variables or\
         ! the TD variables. In both cases, the variables are\
         ! fully cleaned (no NEB with TD).\
         !\
         IF (.NOT. lflag) THEN\
            CALL environ_clean_first(.TRUE.)\
         ELSE\
            CALL environ_clean_second(.TRUE.)\
         END IF\
         !\
      END IF\
      !\
   END IF\
!Environ patch
' tmp.1 >tmp.2

mv tmp.2 plugin_clean.f90

# plugin_summary

sed '/Environ MODULES BEGIN/ a\
!Environ patch \
USE    io_global,            ONLY : stdout \
USE    environ_QE_interface, ONLY : print_environ_summary, & \
                                    update_output_program_unit\
!Environ patch
' plugin_summary.f90 >tmp.1

sed '/Environ CALLS BEGIN/ a\
!Environ patch \
   IF (use_environ) CALL update_output_program_unit( stdout )\
   IF (use_environ) CALL print_environ_summary()\
!Environ patch
' tmp.1 >tmp.2

mv tmp.2 plugin_summary.f90

# plugin_initbase

sed '/Environ MODULES BEGIN/ a\
!Environ patch \
USE    mp_bands,             ONLY : intra_bgrp_comm, me_bgrp, root_bgrp\
USE    cell_base,            ONLY : at, alat\
USE    gvect,                ONLY : gcutm\
USE    environ_QE_interface, ONLY : init_environ_base_second\
!Environ patch
' plugin_initbase.f90 >tmp.1

sed '/Environ VARIABLES BEGIN/ a\
!Environ patch\
INTEGER :: ir_end, idx0, j0, k0\
!Environ patch
' tmp.1 >tmp.2

sed '/Environ CALLS BEGIN/ a\
!Environ patch \
#if defined(__MPI)\
    j0 = dfftp%my_i0r2p\
    k0 = dfftp%my_i0r3p\
    ir_end = MIN(dfftp%nnr, dfftp%nr1x * dfftp%my_nr2p * dfftp%my_nr3p)\
#else\
    j0 = 0\
    k0 = 0\
    ir_end = dfftp%nnr\
#endif\
  IF (use_environ) & \
      CALL init_environ_base_second(alat, at, intra_bgrp_comm, me_bgrp, root_bgrp, gcutm)\
!Environ patch
' tmp.2 >tmp.1

mv tmp.1 plugin_initbase.f90

# plugin_clock

sed '/Environ MODULES BEGIN/ a\
!Environ patch \
USE    environ_QE_interface, ONLY : print_environ_clocks \
!Environ patch
' plugin_clock.f90 >tmp.1

sed '/Environ CALLS BEGIN/ a\
!Environ patch \
   if(use_environ) CALL print_environ_clocks() \
!Environ patch
' tmp.1 >tmp.2

mv tmp.2 plugin_clock.f90

# plugin_print_energies

sed '/Environ MODULES BEGIN/ a\
!Environ patch \
USE    control_flags,        ONLY : conv_elec \
USE    environ_QE_interface, ONLY : print_environ_energies, & \
                                    print_environ_potential_warning \
!Environ patch
' plugin_print_energies.f90 >tmp.1

sed '/Environ CALLS BEGIN/ a\
!Environ patch \
   if (use_environ) then \
     CALL print_environ_energies() \
     if (conv_elec) then \
       CALL print_environ_potential_warning() \
     end if \
   end if \
!Environ patch
' tmp.1 >tmp.2

mv tmp.2 plugin_print_energies.f90

# plugin_init_ions

sed '/Environ MODULES BEGIN/ a\
!Environ patch\
USE cell_base,            ONLY : alat\
USE ions_base,            ONLY : zv, nat, nsp, ityp, tau\
USE environ_QE_interface, ONLY : init_environ_ions\
!Environ patch
' plugin_init_ions.f90 >tmp.1

sed '/Environ CALLS BEGIN/ a\
!Environ patch\
IF (use_environ) CALL init_environ_ions(nat, nsp, ityp, zv, tau, alat)\
!Environ patch
' tmp.1 >tmp.2

mv tmp.2 plugin_init_ions.f90

# plugin_init_cell

sed '/Environ MODULES BEGIN/ a\
!Environ patch\
USE cell_base,            ONLY : at, alat\
USE environ_QE_interface, ONLY : init_environ_cell\
!Environ patch
' plugin_init_cell.f90 >tmp.1

sed '/Environ CALLS BEGIN/ a\
!Environ patch\
  IF ( use_environ ) call init_environ_cell( at, alat )\
!Environ patch
' tmp.1 >tmp.2

mv tmp.2 plugin_init_cell.f90

# plugin_scf_energy

sed '/Environ MODULES BEGIN/ a\
!Environ patch \
USE environ_QE_interface, ONLY : calc_environ_energy, & \
                                 calc_environ_denergy\
!Environ patch
' plugin_scf_energy.f90 >tmp.1

sed '/Environ CALLS BEGIN/ a\
!Environ patch \
  IF(use_environ) THEN \
        ! \
        ! compute environ contributions to total energy\
        ! \
        CALL calc_environ_denergy(plugin_etot)\
        ! \
        CALL calc_environ_energy(plugin_etot)\
        ! \
  END IF \
!Environ patch
' tmp.1 >tmp.2

mv tmp.2 plugin_scf_energy.f90

# plugin_init_potential

sed '/Environ MODULES BEGIN/ a\
!Environ patch \
USE environ_QE_interface, ONLY : init_environ_potential\
!Environ patch
' plugin_init_potential.f90 >tmp.1

sed '/Environ CALLS BEGIN/ a\
!Environ patch \
  IF(use_environ) CALL init_environ_potential( dfftp%nnr, vltot )\
!Environ patch
' tmp.1 >tmp.2

mv tmp.2 plugin_init_potential.f90

# plugin_scf_potential

sed '/Environ MODULES BEGIN/ a\
!Environ patch\
USE global_version,       ONLY : version_number\
USE klist,                ONLY : nelec\
USE control_flags,        ONLY : lscf\
USE lsda_mod,             ONLY : nspin\
USE environ_QE_interface, ONLY : init_environ_electrons, &\
                                 calc_environ_potential, &\
                                 get_environ_threshold, &\
                                 set_environ_restart, &\
                                 is_tddfpt, is_environ_restart, &\
                                 print_environ_potential_shift\
!Environ patch
' plugin_scf_potential.f90 >tmp.1

sed '/Environ VARIABLES BEGIN/ a\
!Environ patch\
LOGICAL :: update_venviron\
INTEGER :: local_verbose\
REAL(DP), ALLOCATABLE :: rhoaux(:)\
!Environ patch
' tmp.1 >tmp.2

sed '/Environ CALLS BEGIN/ a\
!Environ patch\
     IF(use_environ) THEN\
        !\
        ! reduce output at each scf iteration\
        !\
        local_verbose = 0\
        IF ( .NOT. lscf .OR. conv_elec ) local_verbose = 1\
        !\
        ! update electrons-related quantities in environ\
        !\
        ALLOCATE ( rhoaux(dfftp%nnr) )\
        rhoaux(:) = rhoin%of_r(:, 1)\
        !\
        IF ( version_number == "6.3" ) THEN\
            IF ( nspin == 2 ) rhoaux(:) = rhoaux(:) + rhoin%of_r(:, 2)\
        END IF\
        !\
        CALL init_environ_electrons( dfftp%nnr, rhoaux, nelec )\
        !\
        ! environ contribution to the local potential\
        !\
        IF ( dr2 .GT. 0.0_dp ) THEN\
           update_venviron = .NOT. conv_elec .AND. dr2 .LT. get_environ_threshold()\
        !\
        ELSE\
           update_venviron = is_environ_restart() .OR. is_tddfpt()\
           ! for subsequent steps of optimization or dynamics, compute\
           ! environ contribution during initialization\
           IF ( .NOT. is_environ_restart() ) CALL set_environ_restart(.TRUE.)\
        ENDIF\
        !\
        IF ( update_venviron ) WRITE( stdout, 9200 )\
        CALL calc_environ_potential( update_venviron, dfftp%nnr, vltot, local_verbose )\
        !\
        IF ( .NOT. lscf .OR. conv_elec ) THEN\
          CALL print_environ_potential_shift()\
        END IF\
        !\
9200 FORMAT(/"     add environment contribution to local potential")\
     ENDIF\
!Environ patch
' tmp.2 >tmp.1

mv tmp.1 plugin_scf_potential.f90

# plugin_check

sed '/Environ CALLS BEGIN/ a\
!Environ patch \
IF (use_environ) CALL errore( calling_subroutine, &\
   & "Calculation not compatible with Environ embedding", 1)\
!Environ patch
' plugin_check.f90 >tmp.1

mv tmp.1 plugin_check.f90

rm tmp.2

printf " done!\n"

cd $QE_DIR
