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
  USE env_global_objects, ONLY : env\
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
    CALL env%force( nat, force_environ )\
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
  USE io_global,          ONLY : ionode, ionode_id, stdout\
  USE mp_images,          ONLY : intra_image_comm\
  USE martyna_tuckerman,  ONLY : do_comp_mt\
  USE env_base_io,        ONLY : io\
  USE env_global_objects, ONLY : setup\
  USE environ_input,      ONLY : read_environ_input\
!Environ patch
' plugin_read_input.f90 >tmp.1

sed '/Environ VARIABLES BEGIN/ a\
!Environ patch \
  LOGICAL :: lstdout = .FALSE.\
!Environ patch
' tmp.1 >tmp.2

sed '/Environ CALLS BEGIN/ a\
!Environ patch\
   IF (use_environ) THEN\
      !\
      IF (ionode) lstdout = .TRUE.\
      !\
      CALL io%init(ionode, ionode_id, intra_image_comm, stdout, lstdout)\
      !\
      CALL read_environ_input()\
      !\
      CALL setup%init()\
      !\
      IF (prog == "TD") setup%ltddfpt = .TRUE.\
      !\
   ENDIF\
!Environ patch
' tmp.2 >tmp.1

mv tmp.1 plugin_read_input.f90

# plugin_clean

sed '/Environ MODULES BEGIN/ a\
!Environ patch\
USE env_destructor,     ONLY : clean\
USE env_global_objects, ONLY : setup\
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
         IF (.NOT. setup%ltddfpt) CALL clean%all()\
         !\
      ELSE IF ( prog(1:2) == "TD" ) THEN\
         !\
         ! When called by TD, use the flag input variable to\
         ! specify whether to clean the PW variables or\
         ! the TD variables. In both cases, the variables are\
         ! fully cleaned (no NEB with TD).\
         !\
         IF (.NOT. lflag) THEN\
            CALL clean%first()\
         ELSE\
            CALL clean%second()\
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
USE io_global,          ONLY : ionode, stdout \
USE env_global_objects, ONLY : setup\
USE env_base_io,        ONLY : io\
!Environ patch
' plugin_summary.f90 >tmp.1

sed '/Environ CALLS BEGIN/ a\
!Environ patch \
   IF (use_environ) CALL io%update_unit( stdout )\
   IF (use_environ .AND. ionode) CALL setup%print_summary()\
!Environ patch
' tmp.1 >tmp.2

mv tmp.2 plugin_summary.f90

# plugin_initbase

sed '/Environ MODULES BEGIN/ a\
!Environ patch \
USE kinds,              ONLY : DP\
USE mp_bands,           ONLY : intra_bgrp_comm, me_bgrp, root_bgrp\
USE cell_base,          ONLY : at, alat\
USE ions_base,          ONLY : nat, nsp, ityp, atm, zv, tau\
USE martyna_tuckerman,  ONLY : do_comp_mt\
USE gvect,              ONLY : gcutm\
USE env_global_objects, ONLY : env, setup\
!Environ patch
' plugin_initbase.f90 >tmp.1

sed '/Environ VARIABLES BEGIN/ a\
!Environ patch\
INTEGER :: ir_end, idx0, j0, k0\
REAL(DP), ALLOCATABLE :: at_scaled(:, :)\
REAL(DP) :: gcutm_scaled\
CHARACTER(LEN=80) :: sub_name = "plugin_initbase"\
!Environ patch
' tmp.1 >tmp.2

sed '/Environ CALLS BEGIN/ a\
!Environ patch \
  IF (use_environ) THEN\
      !\
      IF (alat < 1.D-8) CALL env_errore(sub_name, "Wrong alat", 1)\
      !\
      IF (alat < 1.0_DP) CALL env_warning("strange lattice parameter")\
      !\
      CALL env_allocate_mp_buffers()\
      !\
      ALLOCATE (at_scaled(3, 3))\
      at_scaled = at * alat\
      !\
      gcutm_scaled = gcutm / alat**2\
      !\
      CALL setup%init_cell(gcutm_scaled, intra_bgrp_comm, at_scaled)\
      !\
      DEALLOCATE (at_scaled)\
      !\
      CALL setup%init_cores(gcutm_scaled, do_comp_mt)\
      !\
      CALL env%init(setup, 1, nat, nsp, atm, ityp, zv)\
      !\
  END IF\
!Environ patch
' tmp.2 >tmp.1

mv tmp.1 plugin_initbase.f90

# plugin_clock

sed '/Environ MODULES BEGIN/ a\
!Environ patch \
USE env_global_objects, ONLY : setup\
!Environ patch
' plugin_clock.f90 >tmp.1

sed '/Environ CALLS BEGIN/ a\
!Environ patch \
   if(use_environ) CALL setup%print_clocks() \
!Environ patch
' tmp.1 >tmp.2

mv tmp.2 plugin_clock.f90

# plugin_print_energies

sed '/Environ MODULES BEGIN/ a\
!Environ patch \
USE control_flags,      ONLY : conv_elec\
USE env_global_objects, ONLY : env, setup\
!Environ patch
' plugin_print_energies.f90 >tmp.1

sed '/Environ CALLS BEGIN/ a\
!Environ patch \
   if (use_environ) then \
     CALL env%print_energies("PW") \
     if (conv_elec) then \
       CALL setup%print_potential_warning() \
     end if \
   end if \
!Environ patch
' tmp.1 >tmp.2

mv tmp.2 plugin_print_energies.f90

# plugin_init_ions

sed '/Environ MODULES BEGIN/ a\
!Environ patch\
USE cell_base,          ONLY : alat\
USE ions_base,          ONLY : nat, tau\
USE env_global_objects, ONLY : env, setup\
!Environ patch
' plugin_init_ions.f90 >tmp.1

sed '/Environ VARIABLES BEGIN/ a\
!Environ patch\
REAL(DP), ALLOCATABLE :: tau_scaled(:, :)\
!Environ patch
' tmp.1 >tmp.2

sed '/Environ CALLS BEGIN/ a\
!Environ patch\
IF (use_environ) THEN\
   ALLOCATE (tau_scaled(3, nat))\
   tau_scaled = tau * alat\
   !\
   CALL env%update_ions(nat, tau_scaled)\
   !\
   DEALLOCATE (tau_scaled)\
END IF\
!Environ patch
' tmp.2 >tmp.1

mv tmp.1 plugin_init_ions.f90

# plugin_init_cell

sed '/Environ MODULES BEGIN/ a\
!Environ patch\
USE cell_base,          ONLY : at, alat\
USE env_global_objects, ONLY : env, setup\
!Environ patch
' plugin_init_cell.f90 >tmp.1

sed '/Environ VARIABLES BEGIN/ a\
!Environ patch\
REAL(DP), ALLOCATABLE :: at_scaled(:, :)\
!Environ patch
' tmp.1 >tmp.2

sed '/Environ CALLS BEGIN/ a\
!Environ patch\
IF ( use_environ ) THEN\
   ALLOCATE (at_scaled(3, 3))\
   at_scaled = at * alat\
   !\
   CALL setup%update_cell(at_scaled)\
   !\
   CALL env%update_cell_dependent_quantities()\
   !\
   CALL setup%end_cell_update()\
   !\
   DEALLOCATE (at_scaled)\
END IF\
!\
!Environ patch
' tmp.2 >tmp.1

mv tmp.1 plugin_init_cell.f90

# plugin_scf_energy

sed '/Environ MODULES BEGIN/ a\
!Environ patch \
USE env_global_objects, ONLY : env\
!Environ patch
' plugin_scf_energy.f90 >tmp.1

sed '/Environ CALLS BEGIN/ a\
!Environ patch \
IF(use_environ) THEN \
   ! \
   ! compute environ contributions to total energy\
   ! \
   CALL env%denergy(plugin_etot)\
   ! \
   CALL env%energy(plugin_etot)\
   ! \
END IF \
!Environ patch
' tmp.1 >tmp.2

mv tmp.2 plugin_scf_energy.f90

# plugin_init_potential

sed '/Environ MODULES BEGIN/ a\
!Environ patch \
USE env_global_objects, ONLY : env\
!Environ patch
' plugin_init_potential.f90 >tmp.1

sed '/Environ CALLS BEGIN/ a\
!Environ patch \
  IF(use_environ) CALL env%update_potential( dfftp%nnr, vltot )\
!Environ patch
' tmp.1 >tmp.2

mv tmp.2 plugin_init_potential.f90

# plugin_scf_potential

sed '/Environ MODULES BEGIN/ a\
!Environ patch\
USE kinds,              ONLY : DP\
USE global_version,     ONLY : version_number\
USE klist,              ONLY : nelec\
USE control_flags,      ONLY : lscf\
USE lsda_mod,           ONLY : nspin\
USE env_global_objects, ONLY : env, setup\
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
        CALL env%update_electrons( dfftp%nnr, rhoaux, nelec )\
        !\
        ! environ contribution to the local potential\
        !\
        IF ( dr2 .GT. 0.0_dp ) THEN\
           update_venviron = .NOT. conv_elec .AND. dr2 .LT. setup%threshold\
        !\
        ELSE\
           update_venviron = setup%restart .OR. setup%ltddfpt\
           ! for subsequent steps of optimization or dynamics, compute\
           ! environ contribution during initialization\
           setup%restart = .TRUE.\
        ENDIF\
        !\
        IF ( update_venviron ) WRITE( stdout, 9200 )\
        !\
        CALL env%potential(update_venviron, local_verbose)\
        !\
        vltot = env%vzero%of_r + env%dvtot%of_r\
        !\
        IF ( .NOT. lscf .OR. conv_elec ) CALL env%print_potential_shift()\
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
