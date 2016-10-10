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
USE environ_base,   ONLY : verbose, environ_unit
USE environ_debug,  ONLY : write_cube
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
      !     goes from 1 to 0 when passing through the threshold
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
        arg = tpi * LOG( xmax / ABS( x ) ) / fact
        dsfunct1 = ( COS( arg ) - 1.D0 ) / ABS( x ) / fact
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
        epsilonfunct = 1.D0 + ( epszero - 1.D0 ) * sfunct0( rho, rhomax, tbeta )
        !
      CASE( 1 )
        !
        epsilonfunct = EXP( LOG( epszero ) * sfunct1( rho, rhomax, rhomin, tbeta ) )
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
        depsilonfunct = ( epszero - 1.D0 ) * dsfunct0( rho, rhomax, tbeta )
        !
      CASE( 1 )
        !
        depsilonfunct = LOG( epszero ) * &
              dsfunct1( rho, rhomax, rhomin, tbeta ) * &
              epsilonfunct( rho, rhomax, rhomin, tbeta, epszero, ifunct )
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
                                 tbeta, rhomax, rhomin, stype,&
                                 solvent_radius
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
      REAL( DP ), DIMENSION(nnr)  :: rhoaug
      !
      IF ( optical_constant ) THEN
         !
         ! TDDFPT calculation
         !
!DEBUG         IF ( env_dielectric_regions .GT. 0 ) THEN
!DEBUG            permittivity = epsoptical
!DEBUG         ELSE ! omogeneous dielectric
            permittivity = env_optical_permittivity
!DEBUG         ENDIF
         !
      ELSE
         !
         ! Ground-state calculation
         !
!DEBUG         IF ( env_dielectric_regions .GT. 0 ) THEN
!DEBUG            permittivity = epsstatic
!DEBUG         ELSE ! omogeneous dielectric
            permittivity = env_static_permittivity
!DEBUG         ENDIF
         !
      ENDIF
      !
      rhoaug = rho
      !
      IF ( env_dielectric_regions .GT. 0 ) CALL fakerhodiel( nnr, rhoaug )
      !
      IF ( solvent_radius .GT. 0.D0 ) THEN
         !
         ! Augment dielectric density to empty
         ! environment pockets smaller than solvent radius
         !
         CALL empty_pockets( nnr, rhoaug )
         !
      ENDIF
      !
      IF ( verbose .GE. 3 ) CALL write_cube( nnr, rhoaug, 'rhoaug.cube' )
      !
      DO ir = 1, nnr
        ! 
        eps( ir ) = epsilonfunct( rhoaug( ir ), rhomax, rhomin, tbeta, &
                                  permittivity( ir ), stype )
        deps( ir ) = depsilonfunct( rhoaug( ir ), rhomax, rhomin, tbeta, &
                                    permittivity( ir ), stype )
        d2eps( ir ) = d2epsilonfunct( rhoaug( ir ), rhomax, rhomin, tbeta, &
                                      permittivity( ir ), stype )
        !
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
!--------------------------------------------------------------------
      SUBROUTINE empty_pockets( nnr, rho )
!--------------------------------------------------------------------
      !
      ! ... Calculates the dielectric constant as a function
      ! ... of the charge density
      !
      USE kinds,          ONLY : DP
      USE environ_base,   ONLY : tbeta, rhomax, rhomin, stype,     &
                                 solvent_radius, radial_spread,    &
                                 radial_scale, emptying_threshold, &
                                 emptying_spread
      USE environ_cell,   ONLY : ntot, domega, alat
      USE io_global,      ONLY : ionode
      USE mp,             ONLY : mp_sum
      USE mp_bands,       ONLY : intra_bgrp_comm
      USE fft_base,       ONLY : dfftp
! BACKWARD COMPATIBILITY !!!! NEED TO BUILD IT PROPERLY
      USE scatter_mod,    ONLY : gather_grid
! END BACKWARD COMPATIBILITY
      USE generate_function, ONLY : generate_axis, generate_erfc
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN)         :: nnr
      REAL( DP ), INTENT(INOUT)   :: rho( nnr )
      !
      ! Local variables
      !
      INTEGER :: ir
      INTEGER :: dim, axis
      REAL( DP ) :: charge, width, spread, ztmp
      REAL( DP ), DIMENSION( 3 ) :: origin = 0.D0
      REAL( DP ), DIMENSION( 3 ) :: pos
      !
      dim = 0
      axis = 1
      charge = 1.D0
      width = solvent_radius * radial_scale
      spread = radial_spread
      !
      RETURN
      !
!--------------------------------------------------------------------
      END SUBROUTINE empty_pockets
!--------------------------------------------------------------------
!--------------------------------------------------------------------
      SUBROUTINE fakerhodiel( nnr, rho )
!--------------------------------------------------------------------
      !
      ! ... Calculates the dielectric constant as a function 
      ! ... of the charge density, and the derivative of 
      ! ... the dielectric constant wrt the charge density.
      !
      USE kinds,          ONLY : DP
      USE environ_base,   ONLY : env_static_permittivity,     &
                                 epsstatic, rhomax, rhomin,   &
                                 verbose
      USE environ_debug,  ONLY : write_cube
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN)         :: nnr
      !
      REAL( DP ), INTENT(OUT)     :: rho( nnr )
      !
      INTEGER :: ir
      !
      rho = 0.D0
      !
      DO ir = 1, nnr
        !
        IF ( epsstatic( ir ) .GE. env_static_permittivity ) THEN
           rho( ir ) = rhomin
        ELSE IF ( epsstatic( ir ) .LE. 1.D0 ) THEN
           rho( ir ) = rhomax
        ELSE
           rho( ir ) = rhomin + ( rhomax - rhomin ) * &
                & ( env_static_permittivity - epsstatic( ir ) ) / &
                & ( env_static_permittivity - 1.D0 )
        END IF
        !
      END DO
      !
      IF ( verbose .GE. 3 ) CALL write_cube( nnr, rho, 'rhofake.cube' )
      !
      RETURN
      !
!--------------------------------------------------------------------
      END SUBROUTINE fakerhodiel
!--------------------------------------------------------------------
END MODULE generate_f_of_rho
