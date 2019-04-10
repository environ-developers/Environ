#!/bin/bash
#
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
# PATCH script for plugin files in CPV/src
#
# Authors: Oliviero Andreussi (Department of Physics, University of North Thexas)
#          Francesco Nattino  (THEOS and NCCR-MARVEL, Ecole Polytechnique Federale de Lausanne)
#

cd $CP_SRC

if test -e "Environ_PATCH" ; then
    echo "-- File Environ_PATCH exists in CPV/src directory"
    echo "-- I guess you have already patched CPV/src with Environ $(tail -1 Environ_PATCH)"
    echo "-- Please unpatch it first, or start from a clean source tree"
    echo "-- See you later..."
    echo "* ABORT"
    exit
fi

echo "* I will try to patch CPV/src with Environ version $ENVIRON_VERSION ..."
echo "#Please do not remove or modify this file"                          >  Environ_PATCH
echo "#It keeps track of patched versions of the Environ addson package" >> Environ_PATCH
echo "$ENVIRON_VERSION"                                                  >> Environ_PATCH

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
  esurface, evolume, eelectrolyte \
USE environ_main,         ONLY : calc_eenviron \
!Environ patch
' plugin_energy.f90 > tmp.1

sed '/Environ CALLS BEGIN/ a\
!Environ patch \
  IF(use_environ) THEN \
     call calc_eenviron( deenviron, eelectrostatic, esurface, evolume, eelectrolyte ) \
     ! \
     plugin_etot = plugin_etot + eelectrostatic + esurface + evolume + eelectrolyte \
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

sed '/Environ VARIABLES BEGIN/ a\
!Environ patch\
! BACKWARD COMPATIBILITY\
! Compatible with QE-6.0 QE-6.1.X QE-6.2.X QE-6.3.X\
!\
! Compatible with QE-6.4.X QE-GIT\
REAL(DP), DIMENSION(:), ALLOCATABLE :: rhoaux\
! END BACKWARD COMPATIBILITY\
!Environ patch
' tmp.1 > tmp.2

sed '/Environ CALLS BEGIN/ a\
!Environ patch\
     IF(use_environ) THEN\
        !\
        ! update electrons-related quantities in environ\
        !\
! BACKWARD COMPATIBILITY\
! Compatible with QE-6.0 QE-6.1.X QE-6.2.X QE-6.3.X\
!        CALL environ_initelectrons( nspin, dfftp%nnr, rhoin )\
! Compatible with QE-6.4.X QE-GIT\
        ALLOCATE(rhoaux(dfftp%nnr))\
        rhoaux = rhoin(:,1)\
        IF ( nspin .EQ. 2) rhoaux = rhoaux + rhoin(:,2)\
        CALL environ_initelectrons( dfftp%nnr, rhoaux )\
        DEALLOCATE(rhoaux)\
! END BACKWARD COMPATIBILITY\
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
! BACKWARD COMPATIBILITY \
! Compatible with QE-5.X QE-6.1.X \
!  ir_end = dfftp%nr1x*dfftp%nr2x*dfftp%npl \
! Compatible with QE-6.2, QE-6.2.1 and QE-GIT \
#if defined (__MPI) \
    ir_end = MIN(dfftp%nnr,dfftp%nr1x*dfftp%my_nr2p*dfftp%my_nr3p) \
#else \
    ir_end = dfftp%nnr \
#endif \
! END BACKWARD COMPATIBILITY \
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
USE input_parameters, ONLY : ion_radius, atom_label \
! BACKWARD COMPATIBILITY\
! Compatible with QE-6.0 QE-6.1.X QE-6.2.X QE-6.3.X\
!USE input_parameters, ONLY : nspin\
! Compatible with QE-6.4.X QE-GIT\
!\
! END BACKWARD COMPATIBILITY\
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
! BACKWARD COMPATIBILITY\
! Compatible with QE-6.0 QE-6.1.X QE-6.2.X QE-6.3.X\
!      CALL read_environ("CP", 1, nspin, nat, ntyp, atom_label, .FALSE., ion_radius)\
! Compatible with QE-6.4.X QE-GIT\
       CALL read_environ("CP", 1, nat, ntyp, atom_label, .FALSE., ion_radius)\
! END BACKWARD COMPATIBILITY\
   ENDIF\
!Environ patch
' tmp.2 > tmp.1

mv tmp.1 plugin_read_input.f90

rm tmp.2

#plugin_utilities.f90

cat > tmp.1 <<EOF
!-----------------------------------------------------------------------
  SUBROUTINE external_laplacian( a, lapla )
!-----------------------------------------------------------------------
      !
      ! Interface for computing hessians in real space, to be called by
      ! an external module
      !
      USE kinds,            ONLY : DP
      USE fft_base,         ONLY : dfftp
      USE gvect,            ONLY : ngm, g
      USE fft_interfaces,   ONLY : fwfft, invfft
      !
      IMPLICIT NONE
      !
      REAL( DP ), INTENT(IN)  :: a( dfftp%nnr )
      REAL( DP ), INTENT(OUT) :: lapla( dfftp%nnr )
      !
      ! ... Locals
      !
      INTEGER :: is
      COMPLEX(DP), ALLOCATABLE :: auxr(:)
      COMPLEX(DP), ALLOCATABLE :: auxg(:)
      REAL(DP), ALLOCATABLE :: d2rho(:,:)
      REAL(DP), ALLOCATABLE :: dxdyrho(:), dxdzrho(:), dydzrho(:), grada(:)
      !
      ALLOCATE( auxg( ngm ) )
      ALLOCATE( auxr( dfftp%nnr ) )
      auxr(:) = CMPLX(a( : ),0.D0,kind=dp)
      CALL fwfft ('Dense', auxr, dfftp)
      auxg(:) = auxr(dfftp%nl(:))
      DEALLOCATE( auxr )
      !
      ALLOCATE( grada(dfftp%nnr) )
      ALLOCATE( d2rho(3,dfftp%nnr) )
      ALLOCATE( dxdyrho(dfftp%nnr) )
      ALLOCATE( dxdzrho(dfftp%nnr) )
      ALLOCATE( dydzrho(dfftp%nnr) )
      ! from G-space A compute R-space grad(A) and second derivatives
      CALL gradrho(1,auxg,grada,d2rho,dxdyrho,dxdzrho,dydzrho)
      DEALLOCATE( auxg )
      ! reorder second derivatives
      lapla(:) = d2rho(1,:)+d2rho(2,:)+d2rho(3,:)
      DEALLOCATE( grada, d2rho, dxdyrho, dxdzrho, dydzrho )

  RETURN

!-----------------------------------------------------------------------
  END SUBROUTINE external_laplacian
!-----------------------------------------------------------------------
EOF

cat > tmp.2 <<EOF
!-----------------------------------------------------------------------
  SUBROUTINE external_force_lc( rho, force )
!-----------------------------------------------------------------------

    USE kinds,             ONLY : DP
    USE cell_base,         ONLY : omega
    USE fft_base,          ONLY : dfftp
    USE fft_interfaces,    ONLY : fwfft
! BACKWARD COMPATIBILITY
! Compatible with QE-6.0 QE-6.1.X QE-6.2.X QE-6.3.X
!    USE electrons_base,    ONLY : nspin
! Compatible with QE-6.4.X QE-GIT
!
! END BACKWARD COMPATIBILITY
    USE gvect,             ONLY : ngm, eigts1, eigts2, eigts3
    USE ions_base,         ONLY : nat

    IMPLICIT NONE
    !
! BACKWARD COMPATIBILITY
! Compatible with QE-6.0 QE-6.1.X QE-6.2.X QE-6.3.X
!    REAL(DP), DIMENSION(dfftp%nnr,nspin), INTENT(IN) :: rho
! Compatible with QE-6.4.X QE-GIT
    REAL(DP), DIMENSION(dfftp%nnr), INTENT(IN) :: rho
! END BACKWARD COMPATIBILITY
    REAL(DP), DIMENSION(3,nat), INTENT(OUT) :: force
    !
    ! aux is used to store a possible additional density
    ! now defined in real space
    !
    COMPLEX(DP), ALLOCATABLE :: auxg(:), auxr(:)
    !
    force = 0.D0
    !
    ALLOCATE( auxr( dfftp%nnr ) )
! BACKWARD COMPATIBILITY
! Compatible with QE-6.0 QE-6.1.X QE-6.2.X QE-6.3.X
!    auxr = CMPLX(rho(:,1),0.0, kind=DP)
!    IF ( nspin .GE. 2 ) auxr = auxr + CMPLX(rho(:,2),0.0, kind=DP)
! Compatible with QE-6.4.X QE-GIT
    auxr = CMPLX(rho,0.0, kind=DP)
! END BACKWARD COMPATIBILITY
    CALL fwfft( "Dense", auxr, dfftp )
    ALLOCATE( auxg( ngm ) )
    auxg(:) = auxr( dfftp%nl (:) )
    !
    CALL force_h_of_rho_g( auxg, eigts1, eigts2, eigts3, omega, force )
    !
    DEALLOCATE( auxr, auxg )
    !
    RETURN
    !
!-----------------------------------------------------------------------
  END SUBROUTINE external_force_lc
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
  SUBROUTINE external_wg_corr_force( rhor, force )
!-----------------------------------------------------------------------

    USE kinds,             ONLY : DP
    USE ions_base,         ONLY : nat
    USE fft_base,          ONLY : dfftp
! BACKWARD COMPATIBILITY
! Compatible with QE-6.0 QE-6.1.X QE-6.2.X QE-6.3.X
!    USE electrons_base,    ONLY : nspin
! Compatible with QE-6.4.X QE-GIT
!
! END BACKWARD COMPATIBILITY
    !
    IMPLICIT NONE
    !
! BACKWARD COMPATIBILITY
! Compatible with QE-6.0 QE-6.1.X QE-6.2.X QE-6.3.X
!    REAL( DP ), INTENT(IN) ::  rhor ( dfftp%nnr, nspin )
! Compatible with QE-6.4.X QE-GIT
    REAL( DP ), INTENT(IN) ::  rhor ( dfftp%nnr )
! END BACKWARD COMPATIBILITY
    REAL( DP ), INTENT(IN) :: force (3, nat)
    !
    CHARACTER(LEN=80) :: sub_name = "external_wg_corr_force"
    !
    ! dummy routine, avoiding cross dependencies with pw objects
    !
    CALL errore(sub_name,&
       "the martyna-tuckerman correction is not available in cp",1)
    !
    RETURN
    !
!-----------------------------------------------------------------------
  END SUBROUTINE external_wg_corr_force
!-----------------------------------------------------------------------
EOF

echo "!Environ patch" >> plugin_utilities.f90
# BACKWARD COMPATIBILITY
# Compatible with QE-6.0 QE-6.1.X QE-6.2.X
#cat tmp.1             >> plugin_utilities.f90
# END BACKWARD COMPATIBILITY
cat tmp.2             >> plugin_utilities.f90
echo "!Environ patch" >> plugin_utilities.f90
rm tmp.1 tmp.2

echo "- DONE!"

cd $QE_DIR
