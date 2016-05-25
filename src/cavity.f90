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
  USE kinds,         ONLY : DP
  USE mp,            ONLY : mp_sum
  USE mp_bands,      ONLY : intra_bgrp_comm
  USE environ_cell,  ONLY : domega
  USE environ_base,  ONLY : verbose, env_surface_tension, delta, e2
  USE environ_debug, ONLY : write_cube
  USE generate_f_of_rho, ONLY : generate_theta
  !
  IMPLICIT NONE
  !
  REAL( DP ), ALLOCATABLE :: theta( : )
  !
  SAVE

  PRIVATE

  PUBLIC :: cavity_initbase, cavity_clean, calc_vcavity, calc_ecavity

CONTAINS
!----------------------------------------------------------------------------
  SUBROUTINE cavity_initbase( nnr )
!----------------------------------------------------------------------------

    ! ... Local initialization

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nnr

    IF(ALLOCATED(theta)) DEALLOCATE(theta)
    ALLOCATE(theta(nnr))
    theta=0.D0
    
    RETURN

!----------------------------------------------------------------------------
  END SUBROUTINE cavity_initbase
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
  SUBROUTINE cavity_clean()
!----------------------------------------------------------------------------

    ! ... Local clean up

    IMPLICIT NONE

    IF(ALLOCATED(theta)) DEALLOCATE(theta)

    RETURN

!----------------------------------------------------------------------------
  END SUBROUTINE cavity_clean
!----------------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE calc_vcavity( nnr, rho, vcavity )
!--------------------------------------------------------------------
    !
    ! ... Calculates the cavitation contribution to the potential
    ! 
    IMPLICIT NONE
    !
    ! ... Declares variables
    !
    INTEGER, INTENT(IN)     :: nnr
    REAL( DP ), INTENT(IN)  :: rho( nnr )
    REAL( DP ), INTENT(OUT) :: vcavity( nnr )
    !
    ! ... Local variables
    !
    REAL( DP )              :: grho2, factor
    REAL( DP ), ALLOCATABLE :: grho( :, : )
    REAL( DP ), ALLOCATABLE :: hrho( :, :, : )
    !
    ! ... Dummy variables
    !
    INTEGER                 :: ir, ipol, jpol
    !
    CALL start_clock ('calc_vcav')
    !
    ! ... Computes step function
    !
    CALL generate_theta( nnr, rho, theta )
    IF ( verbose .GE. 4 ) CALL write_cube( nnr, theta, 'theta.cube' )
    !
    ! ... Computes gradient and hessian of the density
    ! 
    ALLOCATE( grho ( 3, nnr ) )
    ALLOCATE( hrho ( 3, 3, nnr ) )
    CALL external_hessian ( rho, grho, hrho )
    !
    ! ... Computes the cavitation potential 
    !
    DO ir = 1, nnr
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
  SUBROUTINE calc_ecavity( nnr, rho, ecavity )
!--------------------------------------------------------------------
    !
    ! ... Calculates the cavitation contribution to the energy
    !
    IMPLICIT NONE
    !
    ! ... Declares variables
    !
    INTEGER, INTENT(IN)     :: nnr
    REAL( DP ), INTENT(IN)  :: rho( nnr )
    REAL( DP ), INTENT(OUT) :: ecavity
    !
    ! ... Local variables
    !
    INTEGER                 :: ir
    REAL( DP )              :: surface, mod_grho
    REAL( DP ), ALLOCATABLE :: grho( :, : )
    !
    CALL start_clock ('calc_ecav')
    !
    ! ... Initializes the variables
    !
    surface  = 0.D0
    !
    ! ... Computes step function
    !
    CALL generate_theta( nnr, rho, theta )
    !
    ! ... Computes the molecular surface
    !
    ALLOCATE( grho( 3, nnr ) )
    grho = 0.D0
    CALL external_gradient ( rho, grho )
    DO ir = 1, nnr
      mod_grho = SQRT(SUM(grho(:,ir)*grho(:,ir)))
      surface = surface + theta(ir)*mod_grho/delta
    ENDDO
    DEALLOCATE( grho )
    surface = surface * domega 
    !
    CALL mp_sum( surface, intra_bgrp_comm )
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
END MODULE cavity
!--------------------------------------------------------------------
