!
! Copyright (C) 2007-2008 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Module to include an environment pressure potential, proportional to 
! the quantum volume of the system. The method was developed by
! Cococcioni et al. PRL 94, 145501 (2005).
!
! original version by O. Andreussi, N. Marzari (MIT)
!
!----------------------------------------------------------------------------
MODULE pressure
!----------------------------------------------------------------------------
  !
  ! ... The variables needed to compute PV energy and potential corrections
  !
  USE kinds,         ONLY : DP
  USE mp,            ONLY : mp_sum
  USE mp_bands,      ONLY : intra_bgrp_comm
  USE environ_base,  ONLY : env_pressure, e2
  USE environ_cell,  ONLY : domega
  USE environ_debug, ONLY : write_cube
  USE generate_f_of_rho, ONLY : generate_volume, generate_dvoldrho
  !
  IMPLICIT NONE
  !
  SAVE

  PRIVATE

  PUBLIC :: calc_vpressure, calc_epressure

CONTAINS
!--------------------------------------------------------------------
  SUBROUTINE calc_vpressure( nnr, rho, vpressure )
!--------------------------------------------------------------------
    !
    ! ... Calculates the PV contribution to the potential
    !
    IMPLICIT NONE
    !
    ! ... Declares variables
    !
    INTEGER, INTENT(IN)     :: nnr
    REAL( DP ), INTENT(IN)  :: rho( nnr )
    REAL( DP ), INTENT(OUT) :: vpressure( nnr )
    !
    REAL( DP ), ALLOCATABLE :: dvoldrho(:)
    !
    CALL start_clock ('calc_vpre')
    !
    ! ... Computes derivative of step function
    !
    ALLOCATE( dvoldrho( nnr ) )
    !
    CALL generate_dvoldrho( nnr, rho, dvoldrho )
    !
    ! ... Multiply the derivative of the volume by the external pressure 
    !
    vpressure = env_pressure * dvoldrho
    !
    DEALLOCATE( dvoldrho )
    !
    CALL stop_clock ('calc_vpre')
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE calc_vpressure
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE calc_epressure( nnr, rho, epressure )
!--------------------------------------------------------------------
    !
    ! ... Calculates the PV contribution to the energy
    ! 
    IMPLICIT NONE
    !
    ! ... Declares variables
    !
    INTEGER, INTENT(IN)     :: nnr
    REAL( DP ), INTENT(IN)  :: rho( nnr )
    REAL( DP ), INTENT(OUT) :: epressure
    !
    REAL( DP )              :: volume
    REAL( DP ), ALLOCATABLE :: volofrho(:)
    !
    CALL start_clock ('calc_epre')
    !
    ! ... Computes step function
    !
    ALLOCATE( volofrho( nnr ) )
    !
    CALL generate_volume( nnr, rho, volofrho )
    !
    volume = SUM( volofrho ) * domega 
    !
    DEALLOCATE( volofrho )
    !
    CALL mp_sum( volume, intra_bgrp_comm )
    !
    ! ... Computes the PV energy 
    !
    epressure = env_pressure * volume * e2 / 2.D0
    !
    CALL stop_clock ('calc_epre')
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE calc_epressure
!--------------------------------------------------------------------
!--------------------------------------------------------------------
END MODULE pressure
!--------------------------------------------------------------------
