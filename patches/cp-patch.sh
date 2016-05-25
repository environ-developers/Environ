#!/bin/bash

#plugin_add_potential.f90

sed '/Environ MODULES BEGIN/ a\
!Environ patch \
USE environ_base, ONLY: vltot_zero \
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
         v(ir,1)=v(ir,1)+vltot_zero(ir) \
         ! \
       END DO \
!$omp end parallel do \
       ! \
     ELSE IF (nspin.EQ.2) THEN \
       ! \
!$omp parallel do \
       DO ir=1,dfftp%nnr \
         ! \
         v(ir,1)=v(ir,1)+vltot_zero(ir) \
         v(ir,2)=v(ir,2)+vltot_zero(ir) \
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
   if(use_environ) CALL environ_clock(stdout) \
!Environ patch
' tmp.1 > tmp.2

mv tmp.2 plugin_clock.f90

#plugin_energy.f90

sed '/Environ MODULES BEGIN/ a\
!Environ patch \
USE environ_base,         ONLY : deenviron, esolvent, ecavity, epressure, & \
                                 eperiodic, eioncc, eextcharge \
USE environ_main,         ONLY : calc_eenviron \
!Environ patch
' plugin_energy.f90 > tmp.1

sed '/Environ CALLS BEGIN/ a\
!Environ patch \
  IF(use_environ) THEN \
       call calc_eenviron( dfftp%nnr, nspin, rhoin, deenviron, esolvent, & \
                           ecavity, epressure, eperiodic, eioncc, eextcharge ) \
       ! \
       plugin_etot = plugin_etot + esolvent + ecavity + epressure + eperiodic + eioncc + eextcharge \
       ! \
  END IF \
!Environ patch
' tmp.1 > tmp.2

mv tmp.2 plugin_energy.f90

#plugin_get_potential.f90

sed '/Environ MODULES BEGIN/ a\
!Environ patch \
USE environ_base,     ONLY : update_venviron, vltot_zero,    & \
                             environ_nskip, environ_restart, & \
                             verbose \
USE environ_main,     ONLY : calc_venviron \
!Environ patch
' plugin_get_potential.f90 > tmp.1

sed '/Environ CALLS BEGIN/ a\
!Environ patch \
     IF(use_environ) THEN \
        ! \
        ! environ contribution to the local potential, saved in vltot_zero \
        ! \
        vltot_zero = 0.D0 \
        ! \
        update_venviron = ( nfi .GT. environ_nskip ) .OR. environ_restart \
        ! \
        IF ( update_venviron .AND. verbose .GT. 1 ) WRITE( stdout, 9200 ) \
        CALL calc_venviron( update_venviron, dfftp%nnr, nspin, 1.D0, rhoin, vltot_zero ) \
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
!Environ patch
' plugin_init_base.f90 > tmp.1

sed '/Environ CALLS BEGIN/ a\
!Environ patch \
  ir_end = MIN(dfftp%nnr,dfftp%nr1x*dfftp%nr2x*dfftp%npp(me_bgrp+1)) \
  IF ( use_environ ) CALL environ_initbase( dfftp%nnr, 1.D0 ) \
!Environ patch
' tmp.1 > tmp.2

mv tmp.2 plugin_init_base.f90

#plugin_init_cell.f90

sed '/Environ MODULES BEGIN/ a\
!Environ patch \
USE cell_base,        ONLY : at, alat, omega, ibrav \
USE environ_init,     ONLY : environ_initcell \
!Environ patch
' plugin_init_cell.f90 > tmp.1

sed '/Environ CALLS BEGIN/ a\
!Environ patch \
  IF ( use_environ ) call environ_initcell( dfftp%nnr, dfftp%nr1, dfftp%nr2, dfftp%nr3, & \
                           ibrav, omega, alat, at ) \
!Environ patch
' tmp.1 > tmp.2

mv tmp.2 plugin_init_cell.f90

#plugin_init_ions.f90

sed '/Environ MODULES BEGIN/ a\
!Environ patch \
USE cell_base,        ONLY : alat, omega, tpiba2 \
USE ions_base,        ONLY : zv \
USE environ_init,     ONLY : environ_initions \
USE environ_ions,     ONLY : rhoions \
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
     call environ_initions( dfftp%nnr, nat, nsp, ityp_tmp, zv, tau_tmp, alat ) \
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
  USE environ_base,    ONLY : env_static_permittivity, env_dielectric_regions, & \
                              env_external_charges, rhopol, rhoexternal \
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
    ! ... Computes here the contribution coming from the environ charges \
    ! \
    IF ( env_static_permittivity .GT. 1.D0 .OR. & \
         env_dielectric_regions  .GT. 0    .OR. & \
         env_external_charges    .GT. 0 )  THEN \
      ALLOCATE( auxr( dfftp%nnr ) ) \
      auxr = CMPLX(0.0,0.0, kind=DP) \
      if(env_static_permittivity .GT. 1.D0 .OR. & \
         env_dielectric_regions  .GT. 0 ) & \
        auxr = auxr + CMPLX(rhopol(:),0.0, kind=DP) \
      if(env_external_charges .GT. 0) & \
        auxr = auxr + CMPLX(rhoexternal(:),0.0, kind=DP) \
      CALL fwfft( "Dense", auxr, dfftp ) \
      ALLOCATE( auxg( ngm ) ) \
      auxg(:) = auxr( nl (:) ) \
      CALL force_h_of_rho_g( auxg, eigts1, eigts2, eigts3, omega, force_environ ) \
      DEALLOCATE( auxr, auxg ) \
    ENDIF \
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
USE environ_input,    ONLY : read_environ \
USE environ_base,     ONLY : atomicspread \
USE input_parameters, ONLY : ion_radius \
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
      CALL read_environ(nat, ntyp, assume_isolated, ibrav) \
      ! \
      ! ... Overwrites atomicspread with ion_radius from the CP input \
      !     to avoid inconsistency \
      ! \
      DO is = 1, ntyp \
        atomicspread(is) = ion_radius(is) \
      ENDDO \
   ENDIF \
!Environ patch
' tmp.2 > tmp.1

mv tmp.1 plugin_read_input.f90

rm tmp.2
