!
! Copyright (C) 2007-2008 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Module to include a step-like confining potential
!
! original version by O. Andreussi, N. Marzari (MIT)
!
!----------------------------------------------------------------------------
MODULE confine
!----------------------------------------------------------------------------
  !
  ! ... The variables needed to compute PV energy and potential corrections
  !
  USE kinds,         ONLY : DP
  USE mp,            ONLY : mp_sum
  USE mp_bands,      ONLY : intra_bgrp_comm
  USE environ_base,  ONLY : env_confine, e2
  USE environ_cell,  ONLY : domega
  USE environ_debug, ONLY : write_cube
  USE generate_f_of_rho, ONLY : generate_volume, generate_dvoldrho
  !
  IMPLICIT NONE
  !
  SAVE

  PRIVATE

  PUBLIC :: calc_vconfine, calc_econfine

CONTAINS
!--------------------------------------------------------------------
  SUBROUTINE calc_vconfine( nnr, rho, vconfine, vtheta )
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
    REAL( DP ), INTENT(OUT) :: vconfine( nnr ), vtheta( nnr )
    !
    CALL start_clock ('calc_vcon')
    !
    ! ... Compute step function and its derivative 
    !
    CALL generate_volume( nnr, rho, vconfine )
    CALL generate_dvoldrho( nnr, rho, vtheta ) 
    !
    ! ... Multiply by the input value of confinement
    !
    vconfine = env_confine * vconfine
    vtheta = env_confine * vtheta * rho
    !
    CALL stop_clock ('calc_vcon')
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE calc_vconfine
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE calc_econfine( nnr, rho, econfine )
!--------------------------------------------------------------------
    USE environ_base, ONLY : vconfine
    !
    ! ... Calculates the confining energy
    ! 
    IMPLICIT NONE
    !
    ! ... Declares variables
    !
    INTEGER, INTENT(IN)     :: nnr
    REAL( DP ), INTENT(IN)  :: rho( nnr )
    REAL( DP ), INTENT(OUT) :: econfine
    !
    CALL start_clock ('calc_econ')
    !
    ! ... Computes the confining energy 
    !
    econfine = SUM( vconfine( : ) * rho( : ) ) * domega
    CALL mp_sum( econfine, intra_bgrp_comm )
    !
    CALL stop_clock ('calc_econ')
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE calc_econfine
!--------------------------------------------------------------------
!--------------------------------------------------------------------
END MODULE confine
!--------------------------------------------------------------------
