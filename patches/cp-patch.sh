#!/bin/bash
#
# Copyright (C) 2018 ENVIRON (www.quantum-environment.org)
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
# PATCH script for plugin files in CPV/src
#
# Authors: Oliviero Andreussi (Department of Physics, University of North Thexas)
#          Francesco Nattino  (THEOS and NCCR-MARVEL, Ecole Polytechnique Federale de Lausanne)
#

cd $CP_SRC

patch_makefile cp

check_src_patched
if test "$PATCHED" == 1; then 
   return
else
   message "Patching"
fi

echo "#Please do not remove or modify this file"                          >  Environ_PATCH
echo "#It keeps track of patched versions of the Environ addson package" >> Environ_PATCH
echo "$ENVIRON_VERSION"                                                  >> Environ_PATCH

#plugin_add_potential.f90

sed '/Environ MODULES BEGIN/ a\
!Environ patch \
USE base_environ, ONLY: vzero \
!Environ patch
' plugin_add_potential.f90 > tmp.1

sed '/Environ VARIABLES BEGIN/ a\
!Environ patch \
INTEGER :: ir \
!Environ patch
' tmp.1 > tmp.2

sed '/Environ CALLS BEGIN/ a\
!Environ patch \
  IF (use_environ) THEN \
     IF (nspin.EQ.1) THEN \
       ! \
!$omp parallel do \
       DO ir=1,dfftp%nnr \
         ! \
         v(ir,1)=v(ir,1)+vzero%of_r(ir) \
         ! \
       END DO \
!$omp end parallel do \
       ! \
     ELSE IF (nspin.EQ.2) THEN \
       ! \
!$omp parallel do \
       DO ir=1,dfftp%nnr \
         ! \
         v(ir,1)=v(ir,1)+vzero%of_r(ir) \
         v(ir,2)=v(ir,2)+vzero%of_r(ir) \
         ! \
       END DO \
!$omp end parallel do \
       ! \
     END IF \
  END IF \
!Environ patch
' tmp.2 > tmp.1

mv tmp.1 plugin_add_potential.f90

#plugin_clean.f90

sed '/Environ MODULES BEGIN/ a\
!Environ patch \
USE    init_environ, ONLY : environ_clean \
!Environ patch
' plugin_clean.f90 > tmp.1

sed '/Environ CALLS BEGIN/ a\
!Environ patch \
   if(use_environ) CALL environ_clean(.true.) \
!Environ patch
' tmp.1 > tmp.2

mv tmp.2 plugin_clean.f90

#plugin_clock.f90

sed '/Environ MODULES BEGIN/ a\
!Environ patch \
USE    environ_output, ONLY : environ_clock \
!Environ patch
' plugin_clock.f90 > tmp.1

sed '/Environ CALLS BEGIN/ a\
!Environ patch \
   if(use_environ) CALL environ_clock() \
!Environ patch
' tmp.1 > tmp.2

mv tmp.2 plugin_clock.f90

#plugin_energy.f90

sed '/Environ MODULES BEGIN/ a\
!Environ patch \
USE base_environ,          ONLY : deenviron, eelectrostatic, & \
  esurface, evolume, econfine, eelectrolyte \
USE environ_main,         ONLY : calc_eenviron \
!Environ patch
' plugin_energy.f90 > tmp.1

sed '/Environ CALLS BEGIN/ a\
!Environ patch \
  IF(use_environ) THEN \
     call calc_eenviron( deenviron, eelectrostatic, esurface, & \
                       & evolume, econfine, eelectrolyte ) \
     ! \
     plugin_etot = plugin_etot + eelectrostatic + esurface & \
                 & + evolume + econfine + eelectrolyte \
     ! \
  END IF \
!Environ patch
' tmp.1 > tmp.2

mv tmp.2 plugin_energy.f90

#plugin_get_potential.f90

sed '/Environ MODULES BEGIN/ a\
!Environ patch \
USE global_version,   ONLY : version_number\
USE lsda_mod,         ONLY : nspin\
USE base_environ,     ONLY : update_venviron, vzero,    & \
                             environ_nskip, environ_restart \
USE environ_output,   ONLY : verbose \
USE environ_main,     ONLY : calc_venviron \
USE init_environ,     ONLY : environ_initelectrons \
!Environ patch
' plugin_get_potential.f90 > tmp.1

sed '/Environ VARIABLES BEGIN/ a\
!Environ patch\
REAL(DP), ALLOCATABLE :: rhoaux(:)\
!Environ patch
' tmp.1 > tmp.2

sed '/Environ CALLS BEGIN/ a\
!Environ patch\
     IF(use_environ) THEN\
        !\
        ! update electrons-related quantities in environ\
        !\
        ALLOCATE ( rhoaux(dfftp%nnr) )\
        rhoaux(:) = rhoin(:, 1)\
        !\
        IF ( version_number == "6.3" ) THEN\
            IF ( nspin == 2 ) rhoaux(:) = rhoaux(:) + rhoin(:, 2)\
        END IF\
        !\
        CALL environ_initelectrons( dfftp%nnr, rhoaux, nelec )\
        !\
        ! environ contribution to the local potential, saved in vzero\
        !\
        vzero%of_r = 0.D0\
        !\
        update_venviron = ( nfi .GT. environ_nskip ) .OR. environ_restart\
        !\
        IF ( update_venviron .AND. verbose .GT. 1 ) WRITE( stdout, 9200 )\
        CALL calc_venviron( update_venviron, dfftp%nnr, vzero%of_r )\
        !\
9200 FORMAT(/"     add environment contribution to local potential")\
     ENDIF\
!Environ patch
' tmp.2 > tmp.3

mv tmp.3 plugin_get_potential.f90

#plugin_init_base.f90

sed '/Environ MODULES BEGIN/ a\
!Environ patch \
USE    init_environ, ONLY : environ_initbase\
USE    cell_base,    ONLY : at, alat\
USE    gvect,        ONLY : gcutm\
USE    mp_bands,     ONLY : intra_bgrp_comm, me_bgrp, root_bgrp\
!Environ patch
' plugin_init_base.f90 > tmp.1

sed '/Environ CALLS BEGIN/ a\
!Environ patch\
  IF ( use_environ ) CALL environ_initbase(alat, at, &\
                                           & intra_bgrp_comm, me_bgrp, root_bgrp, &\
                                           & gcutm, 1.D0)\
!Environ patch
' tmp.1 > tmp.2

mv tmp.2 plugin_init_base.f90

#plugin_init_cell.f90

sed '/Environ MODULES BEGIN/ a\
!Environ patch \
USE cell_base,        ONLY : at \
USE init_environ,     ONLY : environ_initcell \
!Environ patch
' plugin_init_cell.f90 > tmp.1

sed '/Environ CALLS BEGIN/ a\
!Environ patch \
  IF ( use_environ ) call environ_initcell( at ) \
!Environ patch
' tmp.1 > tmp.2

mv tmp.2 plugin_init_cell.f90

#plugin_init_ions.f90

sed '/Environ MODULES BEGIN/ a\
!Environ patch \
USE cell_base,        ONLY : alat \
USE ions_base,        ONLY : zv \
USE vlocal,           ONLY : vloc \
USE init_environ,     ONLY : environ_initions \
!Environ patch
' plugin_init_ions.f90 > tmp.1

sed '/Environ VARIABLES BEGIN/ a\
!Environ patch \
INTEGER :: i, is, ia \
REAL(DP) :: charge, shift \
REAL(DP) :: rhops, r2, fact \
INTEGER, ALLOCATABLE :: ityp_tmp(:) \
REAL(DP), ALLOCATABLE :: tau_tmp(:,:) \
!Environ patch
' tmp.1 > tmp.2

sed '/Environ CALLS BEGIN/ a\
!Environ patch \
  IF ( use_environ ) THEN \
     ! \
     ! need to rebuild ityp, as the atoms are reshuffled in CP wrt the input \
     ! \
     ALLOCATE(ityp_tmp(nat)) \
     i = 0 \
     DO is = 1, nsp \
       DO ia = 1, na(is) \
         i = i + 1 \
         ityp_tmp(i) = is \
       ENDDO \
     ENDDO \
     ! \
     ! need to convert atomic positions because Environ assumes the same units of PW \
     ! \
     ALLOCATE(tau_tmp(3,nat)) \
     tau_tmp = tau / alat \
     ! \
     call environ_initions( dfftp%nnr, nat, nsp, ityp_tmp, zv, tau_tmp, vloc ) \
     ! \
     DEALLOCATE(ityp_tmp) \
     DEALLOCATE(tau_tmp) \
     ! \
  ENDIF \
!Environ patch
' tmp.2 > tmp.1

mv tmp.1 plugin_init_ions.f90

#plugin_int_forces.f90

sed '/Environ MODULES BEGIN/ a\
!Environ patch \
  USE environ_main,    ONLY : calc_fenviron \
!Environ patch
' plugin_int_forces.f90 > tmp.1

sed '/Environ VARIABLES BEGIN/ a\
!Environ patch \
  REAL(DP), ALLOCATABLE :: force_environ(:,:) \
!Environ patch
' tmp.1 > tmp.2

sed '/Environ CALLS BEGIN/ a\
!Environ patch \
  IF (use_environ) THEN \
    ! \
    ALLOCATE(force_environ(3,nat)) \
    force_environ=0.0_dp \
    ! \
    ! ... Add the other environment contributions \
    ! \
    CALL calc_fenviron( nat, force_environ ) \
    ! \
    force = force + force_environ \
    ! \
    DEALLOCATE(force_environ) \
  END IF \
  ! \
9002 FORMAT(5x,"The dielectric solvent contribution to forces") \
9035 FORMAT(5X,"atom ",I4," type ",I2,"   force = ",3F14.8) \
!Environ patch 
' tmp.2 > tmp.1

mv tmp.1 plugin_int_forces.f90

#plugin_print_energies.f90

sed '/Environ MODULES BEGIN/ a\
!Environ patch \
USE    environ_output, ONLY : environ_print_energies \
!Environ patch
' plugin_print_energies.f90 > tmp.1

sed '/Environ CALLS BEGIN/ a\
!Environ patch \
   if(use_environ) CALL environ_print_energies() \
!Environ patch
' tmp.1 > tmp.2

mv tmp.2 plugin_print_energies.f90

#plugin_print_info.f90

sed '/Environ MODULES BEGIN/ a\
!Environ patch \
USE    environ_output, ONLY : environ_summary \
!Environ patch
' plugin_print_info.f90 > tmp.1

sed '/Environ CALLS BEGIN/ a\
!Environ patch \
   if(use_environ) CALL environ_summary() \
!Environ patch
' tmp.1 > tmp.2

mv tmp.2 plugin_print_info.f90

#plugin_read_input.f90

sed '/Environ MODULES BEGIN/ a\
!Environ patch \
USE io_global,        ONLY : ionode, ionode_id, stdout \
USE mp_images,        ONLY : intra_image_comm \
USE environ_input,    ONLY : read_environ \
USE input_parameters, ONLY : ion_radius, atom_label \
USE environ_output,   ONLY : set_environ_output \
!Environ patch
' plugin_read_input.f90 > tmp.1

sed '/Environ VARIABLES BEGIN/ a\
!Environ patch \
INTEGER :: is \
!Environ patch
' tmp.1 > tmp.2

sed '/Environ CALLS BEGIN/ a\
!Environ patch\
   IF ( use_environ ) THEN\
      CALL set_environ_output("CP", ionode, ionode_id, intra_image_comm, stdout)\
   ENDIF\
!Environ patch
' tmp.2 > tmp.1

mv tmp.1 plugin_read_input.f90

rm tmp.1 tmp.2

printf " done!\n"

cd $QE_DIR
