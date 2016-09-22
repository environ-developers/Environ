!
! Copyright (C) 2009-2014 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! original version by O. Andreussi and N. Marzari
!
MODULE generate_f_of_rho
!
PRIVATE
!
PUBLIC :: epsilonfunct, depsilonfunct, generate_dielectric, &
          generate_theta, generate_volume, generate_dvoldrho
!
CONTAINS
!--------------------------------------------------------------------
      FUNCTION sfunct0( x, xthr, fact )
!--------------------------------------------------------------------
      !
      ! ... Switching function 0: 1+(1-(x/xthr)^fact)/(1+(x/xthr)^fact)
      !     goes from 1 to 0 when passing through the treshold
      !
      USE kinds,              ONLY : DP
      USE constants,          ONLY : tpi
      !
      IMPLICIT NONE
      !
      REAL( DP )             :: sfunct0
      REAL( DP )             :: x, xthr, fact
      !
      ! ... Local variables
      !
      REAL( DP )             :: arg
      !
      arg = ( ABS( x ) / xthr ) ** fact
      sfunct0 = 0.5D0 * ( 1.D0 + ( 1.D0 - arg ) / ( 1.D0 + arg ) )
      !
      RETURN
      !
!--------------------------------------------------------------------
      END FUNCTION sfunct0
!--------------------------------------------------------------------
!--------------------------------------------------------------------
      FUNCTION dsfunct0( x, xthr, fact )
!--------------------------------------------------------------------
      !
      ! ... Derivative of switching function 0
      !
      USE kinds,              ONLY : DP
      USE constants,          ONLY : tpi
      !
      IMPLICIT NONE
      !
      REAL( DP )             :: dsfunct0
      REAL( DP )             :: x, xthr, fact
      !
      ! ... Local variables
      !
      REAL( DP )             :: arg
      !
      arg = ( ABS( x ) / xthr ) ** fact
      dsfunct0 = - fact * ABS( x ) ** ( fact - 1.D0 ) / xthr ** fact &
               &  / ( 1.D0 + arg ) ** 2 
      !
      RETURN
      !
!--------------------------------------------------------------------
      END FUNCTION dsfunct0
!--------------------------------------------------------------------
!--------------------------------------------------------------------
      FUNCTION sfunct1( x, xmax, xmin, fact )
!--------------------------------------------------------------------
      !
      ! ... Switching function 1: x - sin(x)
      !     goes from 1 to 0 when passing from xmin to xmax
      !
      !     NOTE: fact should be equal to LOG(xmax/xmin)
      !     but is passed in input to save time
      !
      USE kinds,              ONLY : DP
      USE constants,          ONLY : tpi
      !
      IMPLICIT NONE
      !
      REAL( DP )             :: sfunct1
      REAL( DP )             :: x, xmax, xmin, fact
      !
      ! ... Local variables
      !
      REAL( DP )             :: arg
      !
      IF ( x .LE. xmin ) THEN
        sfunct1 = 1.D0
      ELSE IF ( x .LT. xmax ) THEN
        arg = tpi * LOG( xmax / ABS( x ) ) / fact
        sfunct1 = ( arg - SIN( arg ) ) / tpi 
      ELSE
        sfunct1 = 0.D0
      ENDIF
      !
      RETURN
      !
!--------------------------------------------------------------------
      END FUNCTION sfunct1
!--------------------------------------------------------------------
!--------------------------------------------------------------------
      FUNCTION dsfunct1( x, xmax, xmin, fact )
!--------------------------------------------------------------------
      !
      ! ... Derivative of switching function 1
      !
      !     NOTE: fact should be equal to LOG(xmax/xmin)
      !     but is passed in input to save time
      !
      USE kinds,              ONLY : DP
      USE constants,          ONLY : tpi
      !
      IMPLICIT NONE
      !
      REAL( DP )             :: dsfunct1
      REAL( DP )             :: x, xmax, xmin, fact
      !
      ! ... Local variables
      !
      REAL( DP )             :: arg
      !
      IF ( x .LE. xmin ) THEN
        dsfunct1 = 0.D0
      ELSE IF ( x .LT. xmax ) THEN
        arg = tpi * LOG(xmax/ABS(x)) / fact
        dsfunct1 = ( COS( arg ) - 1.D0 ) / ABS(x) / fact
      ELSE
        dsfunct1 = 0.D0
      ENDIF
      !
      RETURN
      !
!--------------------------------------------------------------------
      END FUNCTION dsfunct1
!--------------------------------------------------------------------
!--------------------------------------------------------------------
      FUNCTION sfunct2( x, xthr, spread )
!--------------------------------------------------------------------
      !
      ! ... Switching function 2: erfc()
      !     goes from 1 to 0 when passing through xthr
      !
      USE kinds,              ONLY : DP
      !
      IMPLICIT NONE
      !
      REAL( DP )             :: sfunct2
      REAL( DP )             :: x, xthr, spread
      !
      ! ... Local variables
      !
      REAL( DP )             :: arg
      REAL( DP ), EXTERNAL   :: qe_erfc
      !
      arg = ( x - xthr ) / spread
      sfunct2 = 0.5D0 * qe_erfc(arg)
      !
      RETURN
      !
!--------------------------------------------------------------------
      END FUNCTION sfunct2
!--------------------------------------------------------------------
!--------------------------------------------------------------------
      FUNCTION dsfunct2( x, xthr, spread )
!--------------------------------------------------------------------
      !
      ! ... Derivative of switching function 2
      !
      USE kinds,              ONLY : DP
      USE constants,          ONLY : sqrtpi
      !
      IMPLICIT NONE
      !
      REAL( DP )             :: dsfunct2
      REAL( DP )             :: x, xthr, spread
      !
      ! ... Local variables
      !
      REAL( DP )             :: arg
      REAL( DP ), EXTERNAL   :: qe_erfc
      !
      arg = ( x - xthr ) / spread
      dsfunct2 = - EXP( -arg**2 ) / sqrtpi / spread
      !
      RETURN
      !
!--------------------------------------------------------------------
      END FUNCTION dsfunct2
!--------------------------------------------------------------------
!--------------------------------------------------------------------
      FUNCTION epsilonfunct( rho, rhomax, rhomin, tbeta, epszero,  &
                             ifunct )
!--------------------------------------------------------------------
      !
      ! ... Calculates the density-dependent dielectric constant
      ! ... ifunct = 0 => original Fattebert and Gygi function
      !
      USE kinds,              ONLY : DP
      USE constants,          ONLY : tpi
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
      USE kinds,              ONLY : DP
      USE constants,          ONLY : tpi
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
      USE kinds,              ONLY : DP
      USE constants,          ONLY : tpi
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
      SUBROUTINE generate_dielectric( nnr, rho, eps, deps, d2eps, optical_constant )
!--------------------------------------------------------------------
      !
      ! ... Calculates the dielectric constant as a function 
      ! ... of the charge density, and the derivative of 
      ! ... the dielectric constant wrt the charge density.
      !
      USE kinds,          ONLY : DP
      USE environ_base,   ONLY : env_static_permittivity,     &
                                 env_optical_permittivity,    &
                                 env_dielectric_regions,      &
                                 epsstatic, epsoptical,       &
                                 tbeta, rhomax, rhomin, stype
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN)         :: nnr
      !
      REAL( DP ), INTENT(IN)      :: rho( nnr )
      REAL( DP ), INTENT(OUT)     :: eps( nnr )
      REAL( DP ), INTENT(OUT)     :: deps( nnr )
      REAL( DP ), INTENT(OUT)     :: d2eps( nnr )
      LOGICAL, INTENT(IN)         :: optical_constant
      !
      INTEGER                     :: ir
      REAL( DP ), DIMENSION(nnr)  :: permittivity
      !
      IF (optical_constant) THEN
         !
         ! TDDFPT calculation
         !
         IF ( env_dielectric_regions .GT. 0 ) THEN
            permittivity = epsoptical
         ELSE ! omogeneous dielectric
            permittivity = env_optical_permittivity
         ENDIF
         !
      ELSE
         !
         ! Ground-state calculation
         !
         IF ( env_dielectric_regions .GT. 0 ) THEN
            permittivity = epsstatic
         ELSE ! omogeneous dielectric
            permittivity = env_static_permittivity
         ENDIF
         !
      ENDIF
      ! 
      DO ir = 1, nnr
        ! 
        eps( ir )  =  epsilonfunct( rho( ir ), rhomax, rhomin, tbeta, &
                                    permittivity( ir ), stype )
        deps( ir ) = depsilonfunct( rho( ir ), rhomax, rhomin, tbeta, &
                                    permittivity( ir ), stype )
        d2eps( ir ) = d2epsilonfunct( rho( ir ), rhomax, rhomin, tbeta, &
                                    permittivity( ir ), stype )
      END DO
      !
      RETURN
      !
!--------------------------------------------------------------------
      END SUBROUTINE generate_dielectric
!--------------------------------------------------------------------
!--------------------------------------------------------------------
      SUBROUTINE generate_theta( nnr, rho, theta )
!--------------------------------------------------------------------
      !
      ! ... Calculates the dielectric constant as a function
      ! ... of the charge density
      !
      USE kinds,          ONLY : DP
      USE environ_base,   ONLY : env_static_permittivity, &
                                 tbeta, rhomax, rhomin, stype, delta
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN)         :: nnr
      !
      REAL( DP ), INTENT(IN)      :: rho( nnr )
      REAL( DP ), INTENT(OUT)     :: theta( nnr )
      !
      INTEGER                     :: ir
      !
      REAL( DP )                  :: theta_plus, theta_minus
      REAL( DP )                  :: fact, epstmp, rhotmp
      !
      theta = 0.D0
      !
      epstmp = env_static_permittivity 
      IF ( env_static_permittivity .LE. 1.D0 ) epstmp = 2.D0
      fact = - 1.D0 / ( epstmp - 1.D0 )
      !
      DO ir = 1, nnr
        !
        rhotmp = rho( ir ) - delta/2.D0
        theta_plus = epsilonfunct( rhotmp, rhomax, rhomin,     &
                                   tbeta, epstmp, stype )
        !
        rhotmp = rho( ir ) + delta/2.D0
        theta_minus = epsilonfunct( rhotmp, rhomax, rhomin,    &
                                    tbeta, epstmp, stype )
        !
        theta( ir ) = ( theta_minus - theta_plus ) * fact
        !
      END DO
      !
      RETURN
      !
!--------------------------------------------------------------------
      END SUBROUTINE generate_theta
!--------------------------------------------------------------------
!--------------------------------------------------------------------
      SUBROUTINE generate_volume( nnr, rho, volofrho )
!--------------------------------------------------------------------
      !
      ! ... Calculates the dielectric constant as a function
      ! ... of the charge density
      !
      USE kinds,          ONLY : DP
      USE environ_base,   ONLY : env_static_permittivity, &
                                 tbeta, rhomax, rhomin, stype      
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN)         :: nnr
      !
      REAL( DP ), INTENT(IN)      :: rho( nnr )
      REAL( DP ), INTENT(OUT)     :: volofrho( nnr )
      !
      INTEGER                     :: ir
      !
      REAL( DP )                  :: fact, epstmp
      !
      volofrho = 0.D0
      !
      epstmp = env_static_permittivity 
      IF ( env_static_permittivity .LE. 1.D0 ) epstmp = 2.D0
      fact = - 1.D0 / ( epstmp - 1.D0 )
      !
      DO ir = 1, nnr
        !
        volofrho( ir ) = epsilonfunct( rho( ir ), rhomax, rhomin,     &
                                   tbeta, epstmp, stype )
        !
      END DO
      !
      volofrho = 1.D0 + ( volofrho - 1.D0 ) * fact
      !
      RETURN
      !
!--------------------------------------------------------------------
      END SUBROUTINE generate_volume
!--------------------------------------------------------------------
!--------------------------------------------------------------------
      SUBROUTINE generate_dvoldrho( nnr, rho, dvoldrho )
!--------------------------------------------------------------------
      !
      ! ... Calculates the dielectric constant as a function
      ! ... of the charge density
      !
      USE kinds,          ONLY : DP
      USE environ_base,   ONLY : env_static_permittivity, & 
                                 tbeta, rhomax, rhomin, stype 
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN)         :: nnr
      !
      REAL( DP ), INTENT(IN)      :: rho( nnr )
      REAL( DP ), INTENT(OUT)     :: dvoldrho( nnr )
      !
      INTEGER                     :: ir
      !
      REAL( DP )                  :: fact, epstmp
      !
      dvoldrho = 0.D0
      !
      epstmp = env_static_permittivity 
      IF ( env_static_permittivity .LE. 1.D0 ) epstmp = 2.D0
      fact = - 1.D0 / ( epstmp - 1.D0 )
      !
      DO ir = 1, nnr
        !
        dvoldrho( ir ) = depsilonfunct( rho( ir ), rhomax, rhomin, &
                                    tbeta, epstmp, stype )
        !
      END DO
      !
      dvoldrho = dvoldrho * fact
      !
      RETURN
      !
!--------------------------------------------------------------------
      END SUBROUTINE generate_dvoldrho
!--------------------------------------------------------------------
END MODULE generate_f_of_rho
