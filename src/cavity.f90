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
  PUBLIC :: calc_vcavity, calc_ecavity, calc_fcavity
  !
CONTAINS
!--------------------------------------------------------------------
  SUBROUTINE calc_vcavity( boundary, potential )
!--------------------------------------------------------------------
    !
    ! ... Calculates the cavitation contribution to the potential
    !
    IMPLICIT NONE
    !
    ! ... Declares variables
    !
    TYPE( environ_boundary ), TARGET, INTENT(IN) :: boundary
    TYPE( environ_density ), TARGET, INTENT(INOUT) :: potential
    !
    CHARACTER( LEN=80 )     :: sub_name = 'calc_vcavity'
    !
    REAL( DP ), DIMENSION(:), POINTER :: dsurface, dscaled, vcavity
    !
    CALL start_clock ('calc_vcav')
    !
    dsurface => boundary % dsurface % of_r
    vcavity => potential % of_r
    !
    IF ( boundary%need_electrons ) THEN
       !
       dscaled => boundary % dscaled % of_r
       vcavity = env_surface_tension * dsurface * dscaled
       !
    END IF
    !
    CALL stop_clock ('calc_vcav')
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE calc_vcavity
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
    REAL( DP ) :: surface
    !
    CHARACTER( LEN=80 ) :: sub_name = 'calc_ecavity'
    !
    TYPE( environ_cell ), POINTER :: cell
    TYPE( environ_density ) :: sqrt_modulus
    !
    CALL start_clock ('calc_ecav')
    !
    cell => boundary % gradient % cell
    CALL init_environ_density( cell, sqrt_modulus )
    sqrt_modulus % of_r = SQRT( boundary % gradient % modulus % of_r )
    !
    ! ... Computes the molecular surface
    !
    surface = integrate_environ_density( sqrt_modulus )
    !
    CALL destroy_environ_density( sqrt_modulus )
    !
    ! ... Computes the cavitation energy
    !
    ecavity = env_surface_tension * surface * e2 / 2.D0
    !
    CALL stop_clock ('calc_ecav')
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE calc_ecavity
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE calc_fcavity( nat, boundary, forces )
!--------------------------------------------------------------------
    !
    ! ... Calculates the cavitation contribution to the energy
    !
    IMPLICIT NONE
    !
    ! ... Declares variables
    !
    INTEGER, INTENT(IN) :: nat
    TYPE( environ_boundary ), INTENT(IN) :: boundary
    REAL( DP ), DIMENSION(3,nat), INTENT(INOUT) :: forces
    !
    REAL( DP ), DIMENSION(3,nat) :: fcavity
    !
    CHARACTER( LEN=80 ) :: sub_name = 'calc_fcavity'
    !
    CALL start_clock ('calc_fcav')
    !
    fcavity = 0.D0
    !
    ! ... Computes the molecular surface
    !
    IF ( boundary%need_ions ) THEN
       !
       ! ... Rigid cases to be implemented
       !
       CALL errore(sub_name,'Option not yet implemented',1)
       !
    ENDIF
    !
    CALL stop_clock ('calc_fcav')
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE calc_fcavity
!--------------------------------------------------------------------
!--------------------------------------------------------------------
END MODULE cavity
!--------------------------------------------------------------------
