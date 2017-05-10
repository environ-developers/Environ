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
    ! ... Local variables
    !
    REAL( DP )              :: grho2, factor
    REAL( DP ), ALLOCATABLE :: grho( :, : )
    REAL( DP ), ALLOCATABLE :: hrho( :, :, : )
    CHARACTER( LEN=80 )     :: sub_name = 'calc_vcavity'
    !
    ! ... Dummy variables
    !
    INTEGER                 :: ir, ipol, jpol
    !
    INTEGER, POINTER        :: ir_end, nnr
    REAL( DP ), POINTER     :: delta
    REAL( DP ), DIMENSION(:), POINTER :: theta, vcavity, rho
    !
    CALL start_clock ('calc_vcav')
    !
    ir_end => boundary % theta % cell % ir_end
    nnr => boundary % theta % cell % nnr
    delta => boundary % deltatheta
    theta => boundary % theta % of_r
    vcavity => potential % of_r
    !
    IF ( boundary%mode .EQ. 'electrons' ) THEN
       !
       ! ... Computes gradient and hessian of the density
       !
       rho => boundary%electrons%density%of_r
       ALLOCATE( grho ( 3, nnr ) )
       ALLOCATE( hrho ( 3, 3, nnr ) )
       CALL external_hessian ( rho, grho, hrho )
       !
       ! ... Computes the cavitation potential
       !
       DO ir = 1, ir_end
          grho2 = SUM( grho(:, ir) * grho (:, ir) )
          factor = 0.D0
          DO ipol = 1, 3
             DO jpol = 1, 3
                IF ( jpol .EQ. ipol ) CYCLE
                factor = factor + &
                     grho(ipol, ir) * grho(jpol, ir) * hrho(ipol, jpol, ir) - &
                     grho(ipol, ir) * grho(ipol, ir) * hrho(jpol, jpol, ir)
             END DO
          END DO
          IF ( grho2 .GT. 1.D-12 ) &
               vcavity( ir ) = ( theta( ir ) / grho2 / SQRT(grho2) ) * factor
       END DO
       !
       DEALLOCATE( grho )
       DEALLOCATE( hrho )
       !
    ELSE
       !
       ! ... Rigid cases to be implemented
       !
       CALL errore(sub_name,'Option not yet implemented',1)
       !
    ENDIF
    !
    ! ... Multiply the cavitation potential by the constant factor
    !
    vcavity = env_surface_tension / delta * vcavity
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
    REAL( DP )              :: surface
    TYPE( environ_density ) :: density
    TYPE( environ_gradient ) :: gradient
    !
    CHARACTER( LEN=80 ) :: sub_name = 'calc_ecavity'
    !
    REAL( DP ), POINTER :: delta
    TYPE( environ_cell ), POINTER :: cell
    REAL( DP ), DIMENSION(:), POINTER :: theta, rho
    !
    CALL start_clock ('calc_ecav')
    !
    cell => boundary % theta % cell
    theta => boundary % theta % of_r
    delta => boundary % deltatheta
    !
    ! ... Initializes the variables
    !
    surface  = 0.D0
    !
    ! ... Computes the molecular surface
    !
    IF ( boundary%mode .EQ. 'electrons' ) THEN
       !
       rho => boundary % electrons % density % of_r
       !
       CALL init_environ_gradient( cell, gradient )
       CALL external_gradient( rho, gradient%of_r )
       CALL update_gradient_modulus( gradient )
       CALL init_environ_density( cell, density )
       density%of_r = SQRT(gradient%modulus) * theta
       surface = integrate_environ_density( density ) / delta
       CALL destroy_environ_density( density )
       CALL destroy_environ_gradient( gradient )
       !
    ELSE
       !
       ! ... Rigid cases to be implemented
       !
       CALL errore(sub_name,'Option not yet implemented',1)
       !
    ENDIF
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
