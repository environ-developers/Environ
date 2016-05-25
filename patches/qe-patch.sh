#!/bin/bash

# plugin_int_forces

sed '/Environ MODULES BEGIN/ a\
!Environ patch \
  USE environ_base,  ONLY : env_static_permittivity, env_dielectric_regions, rhopol \
  USE environ_base,  ONLY : env_external_charges, rhoexternal \
  USE environ_main,  ONLY : calc_fenviron \
!Environ patch
' plugin_int_forces.f90 > tmp.1

sed '/Environ VARIABLES BEGIN/ a\
!Environ patch \
  REAL(DP), ALLOCATABLE :: force_environ(:,:), force_tmp(:,:) \
!Environ patch
' tmp.1 > tmp.2

sed '/Environ CALLS BEGIN/ a\
!Environ patch \
  IF (use_environ) THEN \
    ! \
    ALLOCATE(force_tmp(3,nat)) \
    ALLOCATE(force_environ(3,nat)) \
    ! \
    force_environ=0.0_dp \
    ! \
    IF(do_comp_mt) THEN \
      force_tmp=0.0_dp \
      ALLOCATE(auxr(dfftp%nnr)) \
      ALLOCATE(auxg(ngm)) \
      auxg = CMPLX(0.0_dp,0.0_dp) \
      auxr = CMPLX(0.0_dp,0.0_dp) \
      if(env_static_permittivity .GT. 1.D0 .OR. env_dielectric_regions .GT. 0) & \
        auxr = CMPLX(rhopol(:),0.0, kind=DP) \
      if(env_external_charges .GT. 0) & \
        auxr = auxr + CMPLX(rhoexternal(:),0.0, kind=DP) \
      CALL fwfft ("Dense", auxr, dfftp) \
      auxg(:)=auxr(nl(:)) \
      CALL wg_corr_force(.false.,omega, nat, ntyp, ityp, ngm, g, tau, zv, strf, & \
                        1, auxg, force_tmp) \
      force_environ = force_environ + force_tmp \
      DEALLOCATE(auxr,auxg) \
      IF ( iverbosity > 0 ) THEN \
        WRITE( stdout, 9001 ) \
        DO na = 1, nat \
           WRITE( stdout, 9035 ) na, ityp(na), ( force_tmp(ipol,na), ipol = 1, 3 ) \
        END DO \
      ENDIF \
    ENDIF \
    ! ... Computes here the solvent contribution \
    ! \
    IF ( env_static_permittivity .GT. 1.D0 .OR. env_dielectric_regions .GT. 0 ) THEN \
      force_tmp=0.0_dp \
      CALL force_lc( nat, tau, ityp, alat, omega, ngm, ngl, igtongl, & \
                     g, rhopol, nl, 1, gstart, gamma_only, vloc, & \
                     force_tmp ) \
      force_environ = force_environ + force_tmp \
      IF ( iverbosity > 0 ) THEN \
        WRITE( stdout, 9002 ) \
        DO na = 1, nat \
           WRITE( stdout, 9035 ) na, ityp(na), ( force_tmp(ipol,na), ipol = 1, 3 ) \
        END DO \
      ENDIF \
    ENDIF \
    ! \
    ! ... Computes here the external charges contribution \
    ! \
    IF ( env_external_charges .GT. 0 ) THEN \
      force_tmp = 0.0_DP \
      CALL force_lc( nat, tau, ityp, alat, omega, ngm, ngl, igtongl, & \
                   g, rhoexternal, nl, 1, gstart, gamma_only, vloc, & \
                   force_tmp ) \
      force_environ = force_environ + force_tmp \
      IF ( iverbosity > 0 ) THEN \
        WRITE( stdout, 9003 ) \
        DO na = 1, nat \
           WRITE( stdout, 9035 ) na, ityp(na), ( force_tmp(ipol,na), ipol = 1, 3 ) \
        END DO \
      ENDIF \
    ENDIF \
    ! \
    ! ... Add the other environment contributions \
    ! \
    CALL calc_fenviron( dfftp%nnr, nspin, nat, force_environ ) \
    ! \
    DEALLOCATE(force_tmp) \
    ! \
    IF ( iverbosity > 0 ) THEN \
      WRITE( stdout, 9004 ) \
      DO na = 1, nat \
         WRITE( stdout, 9035 ) na, ityp(na), ( force_environ(ipol,na), ipol = 1, 3 ) \
      END DO \
      WRITE( stdout, * ) \
    ENDIF \
    ! \
    force = force_environ \
    ! \
    DEALLOCATE(force_environ) \
  END IF \
  ! \
9001 FORMAT(5x,"The MT-Environ correction contrib. to forces") \
9002 FORMAT(5x,"The dielectric solvent contribution to forces") \
9003 FORMAT(5x,"The external charges contribution to forces") \
9004 FORMAT(5x,"The global environment contribution to forces") \
9035 FORMAT(5X,"atom ",I4," type ",I2,"   force = ",3F14.8) \
!Environ patch
' tmp.2 > tmp.1

mv tmp.1 plugin_int_forces.f90

# plugin_read_input

sed '/Environ MODULES BEGIN/ a\
!Environ patch \
USE    environ_input, ONLY : read_environ \
!Environ patch
' plugin_read_input.f90 > tmp.1

sed '/Environ CALLS BEGIN/ a\
!Environ patch \
   if(use_environ) CALL read_environ(nat, ntyp, assume_isolated, ibrav) \
!Environ patch
' tmp.1 > tmp.2

mv tmp.2 plugin_read_input.f90

# plugin_clean

sed '/Environ MODULES BEGIN/ a\
!Environ patch \
USE    environ_init, ONLY : environ_clean \
!Environ patch
' plugin_clean.f90 > tmp.1

sed '/Environ CALLS BEGIN/ a\
!Environ patch \
   if(use_environ) CALL environ_clean(lflag) \
!Environ patch
' tmp.1 > tmp.2

mv tmp.2 plugin_clean.f90

# plugin_summary

sed '/Environ MODULES BEGIN/ a\
!Environ patch \
USE    environ_info, ONLY : environ_summary \
!Environ patch
' plugin_summary.f90 > tmp.1

sed '/Environ CALLS BEGIN/ a\
!Environ patch \
   if(use_environ) CALL environ_summary() \
!Environ patch
' tmp.1 > tmp.2

mv tmp.2 plugin_summary.f90

# plugin_initbase

sed '/Environ MODULES BEGIN/ a\
!Environ patch \
USE    environ_base, ONLY : ir_end \
USE    environ_init, ONLY : environ_initbase \
!Environ patch
' plugin_initbase.f90 > tmp.1

sed '/Environ CALLS BEGIN/ a\
!Environ patch \
  ir_end = MIN(dfftp%nnr,dfftp%nr1x*dfftp%nr2x*dfftp%npp(me_bgrp+1)) \
  IF ( use_environ ) CALL environ_initbase( dfftp%nnr ) \
!Environ patch
' tmp.1 > tmp.2

mv tmp.2 plugin_initbase.f90

# plugin_clock

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

# plugin_print_energies

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

# plugin_init_ions

sed '/Environ MODULES BEGIN/ a\
!Environ patch \
USE cell_base,            ONLY : alat \
USE ions_base,            ONLY : zv, nat, nsp, ityp, tau \
USE environ_init,         ONLY : environ_initions \
!Environ patch
' plugin_init_ions.f90 > tmp.1

sed '/Environ CALLS BEGIN/ a\
!Environ patch \
  IF ( use_environ ) call environ_initions( dfftp%nnr, nat, nsp, ityp, zv, tau, alat ) \
!Environ patch
' tmp.1 > tmp.2

mv tmp.2 plugin_init_ions.f90

# plugin_init_cell

sed '/Environ MODULES BEGIN/ a\
!Environ patch \
USE cell_base,            ONLY : at, alat, omega, ibrav \
USE environ_init,         ONLY : environ_initcell \
!Environ patch
' plugin_init_cell.f90 > tmp.1

sed '/Environ CALLS BEGIN/ a\
!Environ patch \
  IF ( use_environ ) call environ_initcell( dfftp%nnr, dfftp%nr1, dfftp%nr2, dfftp%nr3, & \
                           ibrav, omega, alat, at ) \
!Environ patch
' tmp.1 > tmp.2

mv tmp.2 plugin_init_cell.f90

# plugin_scf_energy

sed '/Environ MODULES BEGIN/ a\
!Environ patch \
USE environ_base,         ONLY : deenviron, esolvent,    & \
                                 ecavity, epressure, eperiodic, eioncc,   & \
                                 eextcharge \
USE environ_main,          ONLY : calc_eenviron \
!Environ patch
' plugin_scf_energy.f90 > tmp.1

sed '/Environ CALLS BEGIN/ a\
!Environ patch \
  IF(use_environ) THEN \
       call calc_eenviron( dfftp%nnr, nspin, rhoin%of_r, deenviron, esolvent, & \
                            ecavity, epressure, eperiodic, eioncc, eextcharge ) \
        ! \
        plugin_etot = plugin_etot + deenviron + esolvent + ecavity + epressure + eperiodic + eioncc + eextcharge \
        ! \
  END IF \
!Environ patch
' tmp.1 > tmp.2

mv tmp.2 plugin_scf_energy.f90

# plugin_init_potential

sed '/Environ MODULES BEGIN/ a\
!Environ patch \
USE environ_base,         ONLY : vltot_zero \
!Environ patch
' plugin_init_potential.f90 > tmp.1

sed '/Environ CALLS BEGIN/ a\
!Environ patch \
  IF(use_environ) vltot_zero = vltot \
!Environ patch
' tmp.1 > tmp.2

mv tmp.2 plugin_init_potential.f90

# plugin_scf_potential

sed '/Environ MODULES BEGIN/ a\
!Environ patch \
USE environ_base,         ONLY : update_venviron, vltot_zero,    & \
                                 environ_thr, environ_restart \
USE environ_main,          ONLY : calc_venviron \
!Environ patch
' plugin_scf_potential.f90 > tmp.1

sed '/Environ CALLS BEGIN/ a\
!Environ patch \
     IF(use_environ) THEN \
        ! \
        ! environ contribution to the local potential \
        ! \
        vltot = vltot_zero \
        !vltot_zero = 0.0_dp \
        ! \
        IF ( dr2 .GT. 0.0_dp ) THEN \
           update_venviron = .NOT. conv_elec .AND. dr2 .LT. environ_thr \
        ! \
        ELSE \
           update_venviron = environ_restart \
           ! for subsequent steps of optimization or dynamics, compute \
           ! environ contribution during initialization \
           IF ( .NOT. environ_restart ) environ_restart = .TRUE. \
        ENDIF \
        ! \
        IF ( update_venviron ) WRITE( stdout, 9200 ) \
        CALL calc_venviron( update_venviron, dfftp%nnr, nspin, dr2, rhoin%of_r, vltot ) \
        ! \
        !CALL sum_vrs( dfftp%nnr, nspin, vltot, vrs, vrs) \
        ! \
9200 FORMAT(/"     add environment contribution to local potential") \
     ENDIF \
!Environ patch
' tmp.1 > tmp.2

mv tmp.2 plugin_scf_potential.f90

# makov_payne

sed '/Environ MODULES BEGIN/ a\
!Environ patch \
USE environ_mp,  ONLY : environ_makov_payne \
!Environ patch
' makov_payne.f90 > tmp.1

sed '/Environ CALLS BEGIN/ a\
!Environ patch \
     IF(use_environ) THEN \
       CALL environ_makov_payne( dfftp%nnr, nspin, rho%of_r, x0, etot ) \
     ENDIF \
!Environ patch
' tmp.1 > tmp.2

mv tmp.2 makov_payne.f90

rm tmp.1
