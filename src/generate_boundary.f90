!
! Copyright (C) 2009-2014 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! original version by O. Andreussi and N. Marzari
!
MODULE generate_boundary
  !
  USE environ_types
  USE constants, ONLY : tpi
  !
  PRIVATE
  !
  PUBLIC :: epsilonfunct, depsilonfunct, d2epsilonfunct, boundary_of_density
  !
CONTAINS
!--------------------------------------------------------------------
      FUNCTION epsilonfunct( rho, rhomax, rhomin, tbeta, epszero,  &
                             ifunct )
!--------------------------------------------------------------------
      !
      ! ... Calculates the density-dependent dielectric constant
      ! ... ifunct = 0 => original Fattebert and Gygi function
      !
      IMPLICIT NONE
      !
      REAL( DP )             :: epsilonfunct
      REAL( DP )             :: rho
      REAL( DP )             :: rhomax
      REAL( DP )             :: rhomin
      REAL( DP )             :: tbeta
      REAL( DP )             :: epszero
      !
      INTEGER                :: ifunct
      !
      REAL( DP )             :: arg
      !
      SELECT CASE( ifunct )
      !
      CASE( 0 )
        !
        epsilonfunct = 1.D0 + 0.5D0 * ( epszero - 1.D0 ) *       &
          ( 1.D0 + ( 1.D0 - ( ABS( rho ) / rhomax ) ** tbeta )   &
          / ( 1.D0 + ( ABS( rho ) / rhomax ) ** tbeta ) )
        !
      CASE( 1 )
        !
        IF ( rho .LE. rhomin ) THEN
          epsilonfunct = epszero
        ELSE IF ( rho .LT. rhomax ) THEN
          arg = tpi * LOG(rhomax/ABS(rho)) / tbeta
          epsilonfunct = EXP( LOG( epszero ) *                     &
            ( arg - SIN( arg ) ) / tpi )
        ELSE
          epsilonfunct = 1.D0
        ENDIF
        !
      CASE DEFAULT
        !
        WRITE(*,*)'ERROR: solvent type unknown, stype=',ifunct
        STOP
        !
      END SELECT
      !
      RETURN
      !
!--------------------------------------------------------------------
      END FUNCTION epsilonfunct
!--------------------------------------------------------------------
!--------------------------------------------------------------------
      FUNCTION depsilonfunct( rho, rhomax, rhomin, tbeta, epszero, &
                              ifunct )
!--------------------------------------------------------------------
      !
      ! ... Calculates the derivative of the
      ! ... density-dependent dielectric constant
      ! ... ifunct = 0 => original Fattebert and Gygi function
      !
      IMPLICIT NONE
      !
      REAL( DP )             :: depsilonfunct
      REAL( DP )             :: rho
      REAL( DP )             :: rhomax
      REAL( DP )             :: rhomin
      REAL( DP )             :: tbeta
      REAL( DP )             :: epszero
      !
      INTEGER                :: ifunct
      !
      REAL( DP )             :: arg
      !
      SELECT CASE( ifunct )
      !
      CASE( 0 )
        !
        depsilonfunct = - tbeta * ( epszero - 1.D0 )               &
          * ABS( rho ) ** ( tbeta - 1.D0 ) / rhomax ** tbeta       &
          / ( 1.D0 + ( ABS( rho ) / rhomax ) ** tbeta ) ** 2
        !
      CASE( 1 )
        !
        IF ( rho .LE. rhomin ) THEN
          depsilonfunct = 0.D0
        ELSE IF ( rho .LT. rhomax ) THEN
          arg = tpi * log(rhomax/ABS(rho)) / tbeta
          depsilonfunct = - EXP( LOG( epszero ) *                  &
            ( arg - SIN( arg ) ) / tpi ) * LOG( epszero ) *        &
            ( 1.D0 - COS( arg ) ) / ABS(rho) / tbeta
        ELSE
          depsilonfunct = 0.D0
        ENDIF
        !
      CASE DEFAULT
        !
        WRITE(*,*)'ERROR: solvent type unknown, stype=',ifunct
        STOP
        !
      END SELECT
      !
      RETURN
      !
!--------------------------------------------------------------------
      END FUNCTION depsilonfunct
!--------------------------------------------------------------------
!--------------------------------------------------------------------
      FUNCTION d2epsilonfunct( rho, rhomax, rhomin, tbeta, epszero, &
                              ifunct )
!--------------------------------------------------------------------
      !
      ! ... Calculates the derivative of the
      ! ... density-dependent dielectric constant
      ! ... ifunct = 0 => original Fattebert and Gygi function
      !
      IMPLICIT NONE
      !
      REAL( DP )             :: d2epsilonfunct
      REAL( DP )             :: rho
      REAL( DP )             :: rhomax
      REAL( DP )             :: rhomin
      REAL( DP )             :: tbeta
      REAL( DP )             :: epszero
      !
      INTEGER                :: ifunct
      !
      REAL( DP )             :: arg, arg2
      !
      SELECT CASE( ifunct )
      !
      CASE( 0 )
        !
        WRITE(*,*)'ERROR: solvent type unknown, stype=',ifunct
        STOP
        !
      CASE( 1 )
        !
        IF ( rho .LE. rhomin ) THEN
          d2epsilonfunct = 0.D0
        ELSE IF ( rho .LT. rhomax ) THEN
          arg = tpi * LOG(rhomax/ABS(rho)) / tbeta
          arg2 =  1.D0 - COS( arg )
          d2epsilonfunct = EXP( LOG( epszero ) *                    &
            ( arg - SIN( arg ) ) / tpi ) * LOG( epszero ) / tbeta * &
            ( LOG( epszero ) / tbeta * arg2 ** 2 + arg2 +           &
              tpi / tbeta * SIN( arg ) ) / rho**2
        ELSE
          d2epsilonfunct = 0.D0
        ENDIF
        !
      CASE DEFAULT
        !
        WRITE(*,*)'ERROR: solvent type unknown, stype=',ifunct
        STOP
        !
      END SELECT
      !
      RETURN
      !
!--------------------------------------------------------------------
      END FUNCTION d2epsilonfunct
!--------------------------------------------------------------------
!--------------------------------------------------------------------
      SUBROUTINE boundary_of_density( density, boundary )
!--------------------------------------------------------------------
      !
      ! ... Calculates the dielectric constant as a function
      ! ... of the charge density, and the derivative of
      ! ... the dielectric constant wrt the charge density.
      !
      IMPLICIT NONE
      !
      TYPE( environ_density ), TARGET, INTENT(IN) :: density
      TYPE( environ_boundary ), TARGET, INTENT(INOUT) :: boundary
      !
      INTEGER, POINTER :: ir_end, stype
      REAL( DP ), POINTER :: constant, rhomax, rhomin, tbeta
      REAL( DP ), DIMENSION(:), POINTER :: rho, eps, deps, d2eps
      !
      REAL( DP ), POINTER :: delta, factor
      REAL( DP ), DIMENSION(:), POINTER :: theta
      !
      INTEGER :: ir
      REAL( DP ) :: rhotmp, theta_plus, theta_minus
      !
      CHARACTER( LEN=80 ) :: sub_name = 'boundary_of_density'
      !
      IF ( .NOT. ASSOCIATED(density%cell,boundary%scaled%cell) ) &
           & CALL errore(sub_name,'Inconsistent domains',1)
      !
      ir_end => density % cell % ir_end
      rho => density % of_r
      !
      constant => boundary % constant ! SHOULD BE REMOVED
      !
      stype => boundary % type
      eps => boundary % scaled % of_r
      deps => boundary % dscaled % of_r
      d2eps => boundary % d2scaled % of_r
      !
      IF ( stype .EQ. 1 ) THEN
         rhomax => boundary % rhomax
         rhomin => boundary % rhomin
         tbeta => boundary % fact
      ELSE IF ( stype .EQ. 2 ) THEN
         rhomax => boundary % rhozero
         rhomin => boundary % deltarho
         tbeta => boundary % tbeta
      ENDIF
      !
      IF ( boundary%need_theta ) THEN
         factor => boundary % scaling_factor ! SHOULD NOT BE HERE
         delta => boundary % deltatheta
         theta => boundary % theta % of_r
      ENDIF
      !
      DO ir = 1, ir_end
        !
        eps( ir )   = epsilonfunct( rho( ir ), rhomax, rhomin, tbeta, &
                                    constant, stype )
        deps( ir )  = depsilonfunct( rho( ir ), rhomax, rhomin, tbeta, &
                                     constant, stype )
        d2eps( ir ) = d2epsilonfunct( rho( ir ), rhomax, rhomin, tbeta, &
                                     constant, stype )
        !
        IF ( boundary % need_theta ) THEN
           !
           rhotmp = rho( ir ) - delta/2.D0
           theta_plus = epsilonfunct( rhotmp, rhomax, rhomin, tbeta, &
                                      constant, stype )
           !
           rhotmp = rho( ir ) + delta/2.D0
           theta_minus = epsilonfunct( rhotmp, rhomax, rhomin, tbeta, &
                                       constant, stype )
           !
           theta( ir ) = ( theta_minus - theta_plus ) * factor
           !
        ENDIF
        !
      END DO
      !
      RETURN
      !
!--------------------------------------------------------------------
    END SUBROUTINE boundary_of_density
!--------------------------------------------------------------------
END MODULE generate_boundary
