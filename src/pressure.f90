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
  USE environ_types
  USE environ_output
  USE environ_base,  ONLY : env_pressure, e2
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC :: calc_depressure_dboundary, calc_epressure
  !
CONTAINS
!--------------------------------------------------------------------
  SUBROUTINE calc_depressure_dboundary( boundary, de_dboundary )
!--------------------------------------------------------------------
    !
    ! ... Calculates the PV contribution to the potential
    !
    IMPLICIT NONE
    !
    ! ... Declares variables
    !
    TYPE( environ_boundary ), INTENT(IN) :: boundary
    TYPE( environ_density ), INTENT(INOUT) :: de_dboundary
    !
    CHARACTER( LEN=80 ) :: sub_name = 'calc_dpressure_dboundary'
    !
    ! ... The functional derivative of the volume term is just unity
    !
    de_dboundary%of_r = de_dboundary%of_r + env_pressure
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE calc_depressure_dboundary
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE calc_epressure( boundary, epressure )
!--------------------------------------------------------------------
    !
    ! ... Calculates the PV contribution to the energy
    !
    IMPLICIT NONE
    !
    ! ... Declares variables
    !
    TYPE( environ_boundary ), TARGET, INTENT(IN) :: boundary
    REAL( DP ), INTENT(OUT) :: epressure
    !
    ! ... Computes the PV energy
    !
    epressure = env_pressure * boundary%volume * e2 / 2.D0
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE calc_epressure
!--------------------------------------------------------------------
!--------------------------------------------------------------------
END MODULE pressure
!--------------------------------------------------------------------
