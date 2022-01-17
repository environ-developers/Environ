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
#          Edan Bainglass     (Department of Physics, UNT)
#
#----------------------------------------------------------------------------------------
#
# PATCH script for plugin files and Makefile in CPV/src
#
#----------------------------------------------------------------------------------------

cd $CP_SRC

patch_makefile

check_src_patched
if test "$PATCHED" == 1; then 
   return
else
   message "Patching"
fi

echo "#Please do not remove or modify this file"                          > Environ_PATCH
echo "#It keeps track of patched versions of the Environ addson package" >> Environ_PATCH
echo "$ENVIRON_VERSION"                                                  >> Environ_PATCH

#plugin_add_potential.f90

sed '/Environ MODULES BEGIN/ a\
!Environ patch\
USE environ_api, ONLY : environ\
!Environ patch
' plugin_add_potential.f90 > tmp.1

sed '/Environ VARIABLES BEGIN/ a\
!Environ patch\
INTEGER :: ir\
REAL(DP), ALLOCATABLE :: dvtot(:)\
!Environ patch
' tmp.1 > tmp.2

sed '/Environ CALLS BEGIN/ a\
!Environ patch\
  IF (use_environ) THEN\
     ALLOCATE (dvtot(dfftp%nnr))\
     dvtot = 0.5D0 * environ%main%get_dvtot(dfftp%nnr)\
     ! Rydberg to Hartree\
     !\
     IF (nspin.EQ.1) THEN\
       !\
!$omp parallel do\
       DO ir=1,dfftp%nnr\
         v(ir,1)=v(ir,1)+dvtot(ir)\
       END DO\
!$omp end parallel do\
       !\
     ELSE IF (nspin.EQ.2) THEN\
       !\
!$omp parallel do\
       DO ir=1,dfftp%nnr\
         v(ir,1)=v(ir,1)+dvtot(ir)\
         v(ir,2)=v(ir,2)+dvtot(ir)\
       END DO\
!$omp end parallel do\
       !\
     END IF\
     !\
     DEALLOCATE (dvtot)\
  END IF\
!Environ patch
' tmp.2 > tmp.1

mv tmp.1 plugin_add_potential.f90

#plugin_clean.f90

sed '/Environ MODULES BEGIN/ a\
!Environ patch\
USE environ_api, ONLY : environ\
!Environ patch
' plugin_clean.f90 > tmp.1

sed '/Environ CALLS BEGIN/ a\
!Environ patch\
   if(use_environ) CALL environ%destroy()\
!Environ patch
' tmp.1 > tmp.2

mv tmp.2 plugin_clean.f90

#plugin_clock.f90

sed '/Environ MODULES BEGIN/ a\
!Environ patch\
USE environ_api, ONLY : environ\
!Environ patch
' plugin_clock.f90 > tmp.1

sed '/Environ CALLS BEGIN/ a\
!Environ patch\
   if(use_environ) CALL environ%setup%print_clocks()\
!Environ patch
' tmp.1 > tmp.2

mv tmp.2 plugin_clock.f90

#plugin_energy.f90

sed '/Environ MODULES BEGIN/ a\
!Environ patch\
USE environ_api, ONLY : environ\
!Environ patch
' plugin_energy.f90 > tmp.1

sed '/Environ CALLS BEGIN/ a\
!Environ patch\
  IF(use_environ) THEN\
        !\
        ! compute environ contributions to total energy\
        !\
        CALL environ%calc%energy(plugin_etot)\
        !\
        plugin_etot = 0.5D0 * plugin_etot\
        ! Rydberg to Hartree\
        !\
  END IF\
!Environ patch
' tmp.1 > tmp.2

mv tmp.2 plugin_energy.f90

#plugin_get_potential.f90

sed '/Environ MODULES BEGIN/ a\
!Environ patch\
USE global_version, ONLY : version_number\
USE electrons_base, ONLY : nelt\
USE environ_api,    ONLY : environ\
!Environ patch
' plugin_get_potential.f90 > tmp.1

sed '/Environ VARIABLES BEGIN/ a\
!Environ patch\
LOGICAL :: update_venviron\
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
        CALL environ%main%update_electrons( dfftp%nnr, rhoaux, Real( nelt, DP) )\
        !\
        ! compute environ contribution to the local potential (dvtot)\
        !\
        update_venviron = ( nfi .GT. environ%setup%get_nskip() ) .OR. environ%setup%is_restart()\
        !\
        IF ( update_venviron .AND. environ%get_verbosity() .GT. 1) WRITE( stdout, 9200 )\
        !\
        CALL environ%calc%potential( update_venviron )\
        !\
9200 FORMAT(/"     add environment contribution to local potential")\
     ENDIF\
!Environ patch
' tmp.2 > tmp.3

mv tmp.3 plugin_get_potential.f90

#plugin_init_base.f90

sed '/Environ MODULES BEGIN/ a\
!Environ patch\
USE kinds,            ONLY : DP\
USE mp_bands,         ONLY : intra_bgrp_comm\
USE cell_base,        ONLY : at, alat\
USE input_parameters, ONLY : nat, ntyp, atom_label\
USE ions_base,        ONLY : ityp, zv\
USE gvect,            ONLY : gcutm\
USE environ_api,      ONLY : environ\
!Environ patch
' plugin_init_base.f90 > tmp.1

sed '/Environ VARIABLES BEGIN/ a\
!Environ patch\
REAL(DP), ALLOCATABLE :: at_scaled(:, :)\
REAL(DP) :: gcutm_scaled\
CHARACTER(LEN=80) :: sub_name = "plugin_init_base"\
!Environ patch
' tmp.1 >tmp.2

sed '/Environ CALLS BEGIN/ a\
!Environ patch\
  IF (use_environ) THEN\
      !\
      IF (alat < 1.D-8) CALL errore(sub_name, "Wrong alat", 1)\
      !\
      ALLOCATE (at_scaled(3, 3))\
      at_scaled = at * alat\
      !\
      gcutm_scaled = gcutm / alat**2\
      !\
      CALL environ%setup%init_cell(intra_bgrp_comm, at_scaled, gcutm=gcutm_scaled)\
      !\
      DEALLOCATE (at_scaled)\
      !\
      CALL environ%setup%init_cores()\
      !\
      CALL environ%main%init(1, nat, ntyp, atom_label, ityp, zv)\
      !\
  END IF\
!Environ patch
' tmp.2 >tmp.1

mv tmp.1 plugin_init_base.f90

# plugin_init_cell

sed '/Environ MODULES BEGIN/ a\
!Environ patch\
USE cell_base,   ONLY : at, alat\
USE environ_api, ONLY : environ\
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
   CALL environ%update_cell(at_scaled)\
   !\
   DEALLOCATE (at_scaled)\
END IF\
!\
!Environ patch
' tmp.2 >tmp.1

mv tmp.1 plugin_init_cell.f90

# plugin_init_ions

sed '/Environ MODULES BEGIN/ a\
!Environ patch\
USE environ_api, ONLY : environ\
!Environ patch
' plugin_init_ions.f90 >tmp.1

sed '/Environ CALLS BEGIN/ a\
!Environ patch\
IF (use_environ) CALL environ%main%update_ions(nat, tau)\
!Environ patch
' tmp.1 >tmp.2

mv tmp.2 plugin_init_ions.f90

#plugin_int_forces.f90

sed '/Environ MODULES BEGIN/ a\
!Environ patch\
USE environ_api, ONLY : environ\
!Environ patch
' plugin_int_forces.f90 > tmp.1

sed '/Environ VARIABLES BEGIN/ a\
!Environ patch\
  REAL(DP), ALLOCATABLE :: force_environ(:,:)\
!Environ patch
' tmp.1 > tmp.2

sed '/Environ CALLS BEGIN/ a\
!Environ patch\
  IF (use_environ) THEN\
    !\
    ALLOCATE(force_environ(3,nat))\
    force_environ=0.0_dp\
    !\
    ! ... Add the other environment contributions\
    !\
    CALL environ%calc%force( nat, force_environ )\
    !\
    force = force + 0.5D0 * force_environ\
    ! Rydberg to Hartree\
    !\
    DEALLOCATE(force_environ)\
  END IF\
  !\
9002 FORMAT(5x,"The dielectric solvent contribution to forces")\
9035 FORMAT(5X,"atom ",I4," type ",I2,"   force = ",3F14.8)\
!Environ patch 
' tmp.2 > tmp.1

mv tmp.1 plugin_int_forces.f90

#plugin_print_energies.f90

sed '/Environ MODULES BEGIN/ a\
!Environ patch\
USE environ_api, ONLY : environ\
!Environ patch
' plugin_print_energies.f90 > tmp.1

sed '/Environ CALLS BEGIN/ a\
!Environ patch\
   if(use_environ) CALL environ%main%print_energies("CP")\
!Environ patch
' tmp.1 > tmp.2

mv tmp.2 plugin_print_energies.f90

#plugin_print_info.f90

sed '/Environ MODULES BEGIN/ a\
!Environ patch\
USE environ_api, ONLY : environ\
!Environ patch
' plugin_print_info.f90 > tmp.1

sed '/Environ CALLS BEGIN/ a\
!Environ patch\
   if(use_environ) CALL environ%setup%print_summary()\
!Environ patch
' tmp.1 > tmp.2

mv tmp.2 plugin_print_info.f90

#plugin_read_input.f90

sed '/Environ MODULES BEGIN/ a\
!Environ patch\
USE io_global,   ONLY : ionode, ionode_id, stdout\
USE mp_images,   ONLY : intra_image_comm\
USE environ_api, ONLY : environ\
!Environ patch
' plugin_read_input.f90 > tmp.1

sed '/Environ CALLS BEGIN/ a\
!Environ patch\
   IF ( use_environ ) THEN\
      !\
      CALL environ%init_interface()\
      !\
      CALL environ%init_io(ionode, ionode_id, intra_image_comm, stdout, ionode)\
      !\
      CALL environ%read_input()\
      !\
      CALL environ%setup%init()\
      !\
   ENDIF\
!Environ patch
' tmp.1 > tmp.2

mv tmp.2 plugin_read_input.f90

rm tmp.1

printf " done!\n"

cd $QE_DIR
