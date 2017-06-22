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
  !
  PRIVATE
  !
  PUBLIC :: boundary_of_density
  !
CONTAINS
!--------------------------------------------------------------------
      FUNCTION sfunct0( x, xthr, fact )
!--------------------------------------------------------------------
      !
      ! ... Switching function 0: 1+(1-(x/xthr)^fact)/(1+(x/xthr)^fact)
      !     goes from 1 to 0 when passing through the threshold
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
        arg = tpi * LOG( xmax / ABS( x ) ) / fact
        dsfunct1 = ( COS( arg ) - 1.D0 ) / ABS( x ) / fact ! in fact should not use ABS( x )
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
      FUNCTION d2sfunct1( x, xmax, xmin, fact )
!--------------------------------------------------------------------
      !
      ! ... Second derivative of switching function 1
      !
      !     NOTE: fact should be equal to LOG(xmax/xmin)
      !     but is passed in input to save time
      !
      USE kinds,              ONLY : DP
      USE constants,          ONLY : tpi
      !
      IMPLICIT NONE
      !
      REAL( DP )             :: d2sfunct1
      REAL( DP )             :: x, xmax, xmin, fact
      !
      ! ... Local variables
      !
      REAL( DP )             :: arg
      !
      IF ( x .LE. xmin ) THEN
        d2sfunct1 = 0.D0
      ELSE IF ( x .LT. xmax ) THEN
        arg = tpi * LOG( xmax / ABS( x ) ) / fact
        d2sfunct1 = ( tpi * SIN( arg ) + fact * ( COS( arg ) - 1.D0 ) ) &
                    / ( x * fact ) ** 2
      ELSE
        d2sfunct1 = 0.D0
      ENDIF
      !
      RETURN
      !
!--------------------------------------------------------------------
      END FUNCTION d2sfunct1
!--------------------------------------------------------------------
!--------------------------------------------------------------------
      FUNCTION sfunct2( x, xthr, spread )
!--------------------------------------------------------------------
      !
      ! ... Switching function 2: erfc()
      !     goes from 1 to 0 when passing through xthr
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
      FUNCTION boundfunct( rho, rhomax, rhomin, tbeta, const, ifunct )
!--------------------------------------------------------------------
      !
      ! ... Calculates the density-dependent dielectric constant
      ! ... ifunct = 0 => original Fattebert and Gygi function
      !
      IMPLICIT NONE
      !
      REAL( DP ) :: boundfunct
      REAL( DP ) :: rho
      REAL( DP ) :: rhomax
      REAL( DP ) :: rhomin
      REAL( DP ) :: tbeta
      REAL( DP ) :: const
      !
      INTEGER :: ifunct
      !
      REAL( DP ) :: arg
      !
      CHARACTER( LEN=80 ) :: fun_name = 'boundfunct'
      !
      SELECT CASE( ifunct )
      !
      CASE( 0 )
         !
         boundfunct = 1.D0 - sfunct0( rho, rhomax, tbeta )
         !
      CASE( 1 )
         !
         boundfunct = 1.D0 - sfunct1( rho, rhomax, rhomin, tbeta )
         !
      CASE( 2 )
         !
         boundfunct = ( const - EXP( LOG( const ) * sfunct1( rho, rhomax, rhomin, tbeta ) ) ) &
              &  / ( const - 1.D0 )
         !
      CASE DEFAULT
         !
         CALL errore(fun_name,'Unknown boundary type',1)
         !
      END SELECT
      !
      RETURN
      !
!--------------------------------------------------------------------
      END FUNCTION boundfunct
!--------------------------------------------------------------------
!--------------------------------------------------------------------
      FUNCTION dboundfunct( rho, rhomax, rhomin, tbeta, const, ifunct )
!--------------------------------------------------------------------
      !
      ! ... Calculates the derivative of the
      ! ... density-dependent dielectric constant
      ! ... ifunct = 0 => original Fattebert and Gygi function
      !
      IMPLICIT NONE
      !
      REAL( DP ) :: dboundfunct
      REAL( DP ) :: rho
      REAL( DP ) :: rhomax
      REAL( DP ) :: rhomin
      REAL( DP ) :: tbeta
      REAL( DP ) :: const
      !
      INTEGER :: ifunct
      !
      REAL( DP ) :: arg
      !
      CHARACTER( LEN=80 ) :: fun_name
      !
      SELECT CASE( ifunct )
      !
      CASE( 0 )
         !
         dboundfunct = - dsfunct0( rho, rhomax, tbeta )
         !
      CASE( 1 )
         !
         dboundfunct = - dsfunct1( rho, rhomax, rhomin, tbeta )
         !
      CASE( 2 )
         !
         dboundfunct = - EXP( LOG( const ) * sfunct1( rho, rhomax, rhomin, tbeta ) ) / &
              & ( const - 1.D0 ) * LOG( const ) * dsfunct1( rho, rhomax, rhomin, tbeta )
         !
      CASE DEFAULT
         !
         CALL errore(fun_name,'Unknown boundary type',1)
         !
      END SELECT
      !
      RETURN
      !
!--------------------------------------------------------------------
      END FUNCTION dboundfunct
!--------------------------------------------------------------------
!--------------------------------------------------------------------
      FUNCTION d2boundfunct( rho, rhomax, rhomin, tbeta, const, ifunct )
!--------------------------------------------------------------------
      !
      ! ... Calculates the derivative of the
      ! ... density-dependent dielectric constant
      ! ... ifunct = 0 => original Fattebert and Gygi function
      !
      IMPLICIT NONE
      !
      REAL( DP ) :: d2boundfunct
      REAL( DP ) :: rho
      REAL( DP ) :: rhomax
      REAL( DP ) :: rhomin
      REAL( DP ) :: tbeta
      REAL( DP ) :: const
      !
      INTEGER :: ifunct
      !
      REAL( DP ) :: arg, arg2
      !
      CHARACTER( LEN=80 ) :: fun_name
      !
      SELECT CASE( ifunct )
      !
      CASE( 0 )
         !
         CALL errore(fun_name,'Second derivative not implemented',1)
         !
      CASE( 1 )
         !
         d2boundfunct = - d2sfunct1( rho, rhomax, rhomin, tbeta )
         !
      CASE( 2 )
         !
         d2boundfunct = - EXP( LOG( const ) * sfunct1( rho, rhomax, rhomin, tbeta ) ) / &
              & ( const - 1.D0 ) * LOG( const ) * ( LOG( const ) * dsfunct1( rho, rhomax, rhomin, tbeta )**2 + &
              & d2sfunct1( rho, rhomax, rhomin, tbeta ) )
         !
      CASE DEFAULT
         !
         CALL errore(fun_name,'Unknown boundary type',1)
         !
      END SELECT
      !
      RETURN
      !
!--------------------------------------------------------------------
      END FUNCTION d2boundfunct
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
      REAL( DP ), POINTER :: const, rhomax, rhomin, tbeta
      REAL( DP ), DIMENSION(:), POINTER :: rho, eps, deps, d2eps
      !
      INTEGER :: ir
      !
      CHARACTER( LEN=80 ) :: sub_name = 'boundary_of_density'
      !
      IF ( .NOT. ASSOCIATED(density%cell,boundary%scaled%cell) ) &
           & CALL errore(sub_name,'Inconsistent domains',1)
      !
      ir_end => density % cell % ir_end
      rho => density % of_r
      !
      stype => boundary % type
      eps => boundary % scaled % of_r
      deps => boundary % dscaled % of_r
      d2eps => boundary % d2scaled % of_r
      !
      IF ( stype .EQ. 1 .OR. stype .EQ. 2 ) THEN
         rhomax => boundary % rhomax
         rhomin => boundary % rhomin
         tbeta => boundary % fact
         const => boundary % const
      ELSE IF ( stype .EQ. 0 ) THEN
         rhomax => boundary % rhozero
         rhomin => boundary % deltarho
         tbeta => boundary % tbeta
         const => boundary % const
      ENDIF
      !
      IF ( boundary%need_theta ) THEN
         delta => boundary % deltatheta
         theta => boundary % theta % of_r
      ENDIF
      !
      DO ir = 1, ir_end
        !
        eps( ir )   = boundfunct( rho( ir ), rhomax, rhomin, tbeta, const, stype )
        deps( ir )  = dboundfunct( rho( ir ), rhomax, rhomin, tbeta, const, stype )
        d2eps( ir ) = d2boundfunct( rho( ir ), rhomax, rhomin, tbeta, const, stype )
        !
      END DO
      !
      RETURN
      !
!--------------------------------------------------------------------
    END SUBROUTINE boundary_of_density
!--------------------------------------------------------------------
END MODULE generate_boundary
