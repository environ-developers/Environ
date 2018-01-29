#!/bin/bash

#plugin_add_potential.f90

sed '/Environ MODULES BEGIN/ a\
!Environ patch \
USE environ_base, ONLY: vzero \
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
USE    environ_init, ONLY : environ_clean \
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
USE    environ_info, ONLY : environ_clock \
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
USE environ_base,          ONLY : deenviron, eelectrostatic, & \
  ecavity, epressure, eelectrolyte \
USE environ_main,         ONLY : calc_eenviron \
!Environ patch
' plugin_energy.f90 > tmp.1

sed '/Environ CALLS BEGIN/ a\
!Environ patch \
  IF(use_environ) THEN \
     call calc_eenviron( deenviron, eelectrostatic, ecavity, epressure, eelectrolyte ) \
     ! \
     plugin_etot = plugin_etot + eelectrostatic + ecavity + epressure + eelectrolyte \
     ! \
  END IF \
!Environ patch
' tmp.1 > tmp.2

mv tmp.2 plugin_energy.f90

#plugin_get_potential.f90

sed '/Environ MODULES BEGIN/ a\
!Environ patch \
USE environ_base,     ONLY : update_venviron, vzero,    & \
                             environ_nskip, environ_restart \
USE environ_output,   ONLY : verbose \
USE environ_main,     ONLY : calc_venviron \
USE environ_init,     ONLY : environ_initelectrons \
!Environ patch
' plugin_get_potential.f90 > tmp.1

sed '/Environ CALLS BEGIN/ a\
!Environ patch \
     IF(use_environ) THEN \
        ! \
        ! update electrons-related quantities in environ \
        ! \
        CALL environ_initelectrons( nspin, dfftp%nnr, rhoin ) ! \
        ! \
        ! environ contribution to the local potential, saved in vzero \
        ! \
        vzero%of_r = 0.D0 \
        ! \
        update_venviron = ( nfi .GT. environ_nskip ) .OR. environ_restart \
        ! \
        IF ( update_venviron .AND. verbose .GT. 1 ) WRITE( stdout, 9200 ) \
        CALL calc_venviron( update_venviron, dfftp%nnr, vzero%of_r ) \
        ! \
9200 FORMAT(/"     add environment contribution to local potential") \
     ENDIF \
!Environ patch
' tmp.1 > tmp.2

mv tmp.2 plugin_get_potential.f90

#plugin_init_base.f90

sed '/Environ MODULES BEGIN/ a\
!Environ patch \
USE    environ_base, ONLY : ir_end \
USE    environ_init, ONLY : environ_initbase \
USE    cell_base,    ONLY : at, alat, omega, ibrav \
USE    mp_bands,     ONLY : intra_bgrp_comm, me_bgrp, root_bgrp_id \
!Environ patch
' plugin_init_base.f90 > tmp.1

sed '/Environ VARIABLES BEGIN/ a\
!Environ patch \
INTEGER :: ir_end \
!Environ patch
' tmp.1 > tmp.2

sed '/Environ CALLS BEGIN/ a\
!Environ patch \
  ir_end = MIN(dfftp%nnr,dfftp%nr1x*dfftp%nr2x*dfftp%npp(me_bgrp+1)) \
  IF ( use_environ ) CALL environ_initbase( dfftp%nr1, dfftp%nr2, dfftp%nr3, ibrav, alat, omega, at, & \
       & dfftp%nnr, ir_end, intra_bgrp_comm, me_bgrp, root_bgrp_id, 1.D0 ) \
!Environ patch
' tmp.2 > tmp.1

mv tmp.1 plugin_init_base.f90

#plugin_init_cell.f90

sed '/Environ MODULES BEGIN/ a\
!Environ patch \
USE cell_base,        ONLY : at, omega \
USE environ_init,     ONLY : environ_initcell \
!Environ patch
' plugin_init_cell.f90 > tmp.1

sed '/Environ CALLS BEGIN/ a\
!Environ patch \
  IF ( use_environ ) call environ_initcell( omega, at ) \
!Environ patch
' tmp.1 > tmp.2

mv tmp.2 plugin_init_cell.f90

#plugin_init_ions.f90

sed '/Environ MODULES BEGIN/ a\
!Environ patch \
USE cell_base,        ONLY : alat, omega, tpiba2 \
USE ions_base,        ONLY : zv \
USE environ_init,     ONLY : environ_initions \
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
     call environ_initions( dfftp%nnr, nat, nsp, ityp_tmp, zv, tau_tmp ) \
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
    CALL calc_fenviron( dfftp%nnr, nspin, nat, force_environ ) \
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
USE    environ_info, ONLY : environ_print_energies \
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
USE    environ_info, ONLY : environ_summary \
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
USE input_parameters, ONLY : ion_radius, atom_label, nspin \
USE environ_output,   ONLY : set_environ_output \
!Environ patch
' plugin_read_input.f90 > tmp.1

sed '/Environ VARIABLES BEGIN/ a\
!Environ patch \
INTEGER :: is \
!Environ patch
' tmp.1 > tmp.2

sed '/Environ CALLS BEGIN/ a\
!Environ patch \
   IF ( use_environ ) THEN \
      CALL set_environ_output("CP", ionode, ionode_id, intra_image_comm, stdout) \
      CALL read_environ("CP",1, nspin, nat, ntyp, atom_label, assume_isolated, ion_radius) \
   ENDIF \
!Environ patch
' tmp.2 > tmp.1

mv tmp.1 plugin_read_input.f90

rm tmp.2
