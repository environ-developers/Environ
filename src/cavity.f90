!
! Copyright (C) 2007-2008 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Module to compute a cavitation potential, defined as the quantum surface
! of the system times the surface tension of the environment.
! Original method developed in Scherlis et al, J. Chem. Phys. (2006),
! but formulas for the quantum surface were not correct and have been
! rederived by Andreussi et al., J. Chem. Phys. 136, 064102 (2012)
!
! original version by O. Andreussi, M. Cococcioni, N. Marzari (MIT)
!
!----------------------------------------------------------------------------
MODULE cavity
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
  !
  ! ... The variables needed to introduce an environment with a surface tension
  !     around the system.
  !
  USE environ_types
  USE environ_output
  USE environ_base,  ONLY : env_surface_tension, e2
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC :: calc_decavity_dboundary, calc_ecavity
  !
CONTAINS
!--------------------------------------------------------------------
  SUBROUTINE calc_decavity_dboundary( boundary, de_dboundary )
!--------------------------------------------------------------------
    !
    ! ... Calculates the cavitation contribution to the potential
    !
    IMPLICIT NONE
    !
    ! ... Declares variables
    !
    TYPE( environ_boundary ), TARGET, INTENT(IN) :: boundary
    TYPE( environ_density ), TARGET, INTENT(INOUT) :: de_dboundary
    !
    CHARACTER( LEN=80 )     :: sub_name = 'calc_dcavity_dboundary'
    !
    de_dboundary % of_r = de_dboundary % of_r + &
         & env_surface_tension * boundary % dsurface % of_r
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE calc_decavity_dboundary
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE calc_ecavity( boundary, ecavity )
!--------------------------------------------------------------------
    !
    ! ... Calculates the cavitation contribution to the energy
    !
    IMPLICIT NONE
    !
    ! ... Declares variables
    !
    TYPE( environ_boundary ), TARGET, INTENT(IN) :: boundary
    REAL( DP ), INTENT(OUT) :: ecavity
    !
    ! ... Local variables
    !
    CHARACTER( LEN=80 ) :: sub_name = 'calc_ecavity'
    !
    ! ... Computes the cavitation energy
    !
    ecavity = env_surface_tension * boundary%surface * e2 / 2.D0
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE calc_ecavity
!--------------------------------------------------------------------
!--------------------------------------------------------------------
END MODULE cavity
!--------------------------------------------------------------------
