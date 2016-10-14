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
                                 solvent_radius, radial_scale,&
                                 radial_spread
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
      REAL( DP ), DIMENSION(nnr)  :: rhofake
      REAL( DP ), DIMENSION(nnr)  :: df1, df2, ftmp
      REAL( DP ) :: plus, minus, deltarho, depstmp
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
      IF ( env_dielectric_regions .GT. 0 ) THEN
         !
         CALL fakerhodiel( nnr, rhofake )
         rhoaug = rhofake
         !
      ENDIF
      !
      IF ( solvent_radius .GT. 0.D0 ) THEN
         !
         ! Augment dielectric density to empty
         ! environment pockets smaller than solvent radius
         !
         CALL empty_pockets( nnr, rhoaug, df1, df2 )
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
      IF ( verbose .GE. 3 ) CALL write_cube( nnr, deps, 'deps0.cube' )
      !
      IF ( solvent_radius .GT. 0.D0 ) THEN
         !
         ! Update the functional derivative of epsilon
         !
         ftmp = deps * df2
         CALL compute_convolution( nnr, solvent_radius*radial_scale, radial_spread, ftmp, df2 )
         IF ( verbose .GE. 3 ) CALL write_cube( nnr, deps, 'df2conv.cube' )
         deps(:) = deps(:) + df1(:) * df2(:)
         !
      ENDIF
      !
      IF ( verbose .GE. 3 ) CALL write_cube( nnr, deps, 'deps.cube' )
      !
      deltarho = 0.000000001D0
      DO ir = 62077, nnr, 14400 !nnr, 14400
         rhoaug = rhofake
         rhoaug( ir ) = rhofake( ir ) + deltarho
         CALL empty_pockets( nnr, rhoaug, df1, df2 )
!         DO ir = nnr
         plus = epsilonfunct( rhoaug( ir ), rhomax, rhomin, tbeta, permittivity(ir),stype)
         WRITE(environ_unit,*)'plus ',rhofake(ir),rhoaug(ir),plus
!         ENDDO
         rhoaug = rhofake
         rhoaug( ir ) = rhofake( ir ) - deltarho
         CALL empty_pockets( nnr, rhoaug, df1, df2 )
!         DO ir = nnr
            minus = epsilonfunct( rhoaug( ir ), rhomax, rhomin, tbeta, permittivity(ir),stype)
!         ENDDO
         WRITE(environ_unit,*)'minus ',rhofake(ir),rhoaug(ir),minus
         depstmp = ( plus - minus ) / 2.D0 / deltarho
         WRITE(environ_unit,'(1X,a,i8,3f20.10)')' ir = ',ir,deps(ir),depstmp,deps(ir)-depstmp 
         FLUSH(environ_unit)
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
      SUBROUTINE empty_pockets( nnr, rho, df1, df2 )
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
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN)         :: nnr
      REAL( DP ), INTENT(INOUT)   :: rho( nnr )
      REAL( DP ), INTENT(OUT)     :: df1( nnr ), df2( nnr )
      !
      ! Local variables
      !
      INTEGER :: ir, ir_end, idx0, i, j, k, ip, itmp, jtmp, jr
      REAL( DP ) :: inv_nr1, inv_nr2, inv_nr3
      REAL( DP ) :: width, spread
      REAL( DP ), DIMENSION( nnr ) :: f, g
      !
      ! Step 1: compute scaled dielectric function and its derivative
      !
      f = 0.D0
      df1 = 0.D0
      !
      DO ir = 1, nnr ! WARNING: CHECK THAT nnr INSTEAD OF ir_end IS OK
         !
         f( ir ) = sfunct1( rho( ir ), rhomax, rhomin, tbeta )
         df1( ir ) = dsfunct1( rho( ir ), rhomax, rhomin, tbeta )
         !
      ENDDO
      !
      IF ( verbose .GE. 3 ) CALL write_cube( nnr, f, 'scaledeps.cube' )
      IF ( verbose .GE. 3 ) CALL write_cube( nnr, df1, 'dscaledeps.cube' )
      !
      ! Step 2: compute filled fraction, i.e. convolution of scaledeps with erfc
      !
      width = solvent_radius * radial_scale
      spread = radial_spread
      CALL compute_convolution( nnr, width, spread, f, g )
      !
      IF ( verbose .GE. 3 ) CALL write_cube( nnr, g, 'filledfrac.cube' )
      !
      ! Step 3: compute the filling condition and its derivative
      !
      f = 0.D0
      df2 = 0.D0
      DO ir = 1, nnr ! WARNING: CHECK THAT nnr INSTEAD OF ir_end IS OK
         !
         f( ir ) = sfunct2( g( ir ), emptying_threshold, emptying_spread )
         df2( ir ) = dsfunct2( g( ir ), emptying_threshold, emptying_spread ) 
         !
      ENDDO
      !
      IF ( verbose .GE. 3 ) CALL write_cube( nnr, f, 'rhoholes.cube' )
      IF ( verbose .GE. 3 ) CALL write_cube( nnr, df2, 'drhoholes.cube' )
      !
      rho = rho + 2.D0 * rhomax * f
      df2 = 2.D0 * rhomax * df2
      !
      RETURN
      !
!--------------------------------------------------------------------
      END SUBROUTINE empty_pockets
!--------------------------------------------------------------------
!--------------------------------------------------------------------
      SUBROUTINE compute_convolution( nnr, width, spread, fin, fout )
!--------------------------------------------------------------------
      !
      ! ... Calculates on a subset of the real-space grid the 
      ! ... convolution of function fin with a fast-decaying function
      ! ... at this time the only fconv implemented is
      ! ... erfc( ( |r'-r| - width ) / spread )
      !
      USE kinds,          ONLY : DP
      USE cell_base,      ONLY : alat, at, omega
      USE mp,             ONLY : mp_sum
      USE mp_bands,       ONLY : intra_bgrp_comm
      USE fft_base,       ONLY : dfftp
! BACKWARD COMPATIBILITY !!!! NEED TO BUILD IT PROPERLY
      USE scatter_mod,    ONLY : gather_grid
! END BACKWARD COMPATIBILITY
      USE generate_function, ONLY : erfcvolume
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN)         :: nnr
      REAL( DP ), INTENT(IN)      :: width, spread
      REAL( DP ), INTENT(IN)      :: fin( nnr )
      REAL( DP ), INTENT(OUT)     :: fout( nnr )
      !
      ! Local variables
      !
      REAL( DP ) :: inv_nr1, inv_nr2, inv_nr3
      INTEGER :: ir, ir_end, idx0, i, j, k, ip, itmp, jtmp, jr
      INTEGER :: nr1c, nr2c, nr3c, ntot
      INTEGER :: dim, axis
      REAL( DP ) :: maxwidth, scale
      REAL( DP ) :: r, r2, arg
      REAL( DP ), EXTERNAL :: qe_erfc
      REAL( DP ), DIMENSION( : ), ALLOCATABLE :: ftot
      INTEGER, DIMENSION( : ), ALLOCATABLE :: i1, i2, i3
      REAL( DP ), DIMENSION( :, :, : ), ALLOCATABLE :: fconv
      !
      ! ... Initialise output function
      !
      fout = 0.D0
      !
      ! ... Each processor needs to have the whole function
      !
      ntot = dfftp%nr1x * dfftp%nr2x * dfftp%nr3x ! NOT SURE IT NEEDS THE X VALUES
      !
      ALLOCATE( ftot(ntot) )
      !
#if defined (__MPI)
      ftot = 0.D0
      CALL gather_grid ( dfftp, fin, ftot )
      CALL mp_sum( ftot, intra_bgrp_comm )
#else
      ftot = fin
#endif
      !
      ! ... Prepare to run on real-space grid
      !
      inv_nr1 = 1.D0 / DBLE( dfftp%nr1 )
      inv_nr2 = 1.D0 / DBLE( dfftp%nr2 )
      inv_nr3 = 1.D0 / DBLE( dfftp%nr3 )
      !
      idx0 = 0
      ir_end = dfftp%nnr
      !
#if defined (__MPI)
      DO i = 1, me_bgrp
        idx0 = idx0 + dfftp%nr1x*dfftp%nr2x*dfftp%npp(i)
      END DO
      ir_end = MIN(nnr,dfftp%nr1x*dfftp%nr2x*dfftp%npp(me_bgrp+1))
#endif
      !
      ! ... Compute the size of the convolution grid
      !
      maxwidth = ( width + 3.D0 * spread ) / alat ! factor 3.D0 corresponds to neglect fconv<=ERFC(3.D0)
      nr1c = INT(maxwidth/SQRT(SUM(at(:,1)**2))*dfftp%nr1x) + 1
      nr2c = INT(maxwidth/SQRT(SUM(at(:,2)**2))*dfftp%nr2x) + 1
      nr3c = INT(maxwidth/SQRT(SUM(at(:,3)**2))*dfftp%nr3x) + 1
      !
      ! ... Generate the convolution function on the conv. grid
      !
      ALLOCATE( fconv( 0 : nr1c, 0 : nr2c, 0 : nr3c ) )
      fconv = 0.D0
      DO i = 0, nr1c
         DO j = 0, nr2c
            DO k = 0, nr3c
               r2 = 0.D0
               DO ip = 1, 3
                  r = DBLE( i )*inv_nr1*at(ip,1) + &
                      DBLE( j )*inv_nr2*at(ip,2) + &
                      DBLE( k )*inv_nr3*at(ip,3)
                  r2 = r2 + r**2 
               ENDDO
               r = SQRT(r2) * alat
               ! ERFC convolution function BEGIN
               arg = ( r - width ) / spread
               fconv(i,j,k) = qe_erfc(arg)
               ! ERFC convolution function END
            ENDDO
         ENDDO         
      ENDDO
      !
      ! ... Normalize convolution function and multiply by finite difference volume
      !
      ! ERFC convolution function BEGIN
      dim = 0
      axis = 1
      scale = 0.5D0 / erfcvolume(dim,axis,width,spread,alat,omega,at) * omega / DBLE(ntot)
      ! ERFC convolution function END
      fconv = fconv * scale 
      !
      ! Step 4: Perform the convolution
      !
      ALLOCATE( i1( -nr1c : nr1c ) )
      ALLOCATE( i2( -nr2c : nr2c ) )
      ALLOCATE( i3( -nr3c : nr3c ) )
      DO ir = 1, ir_end
         !
         ! ... if the point is already empty we can skip it
         !
         IF ( fin( ir ) .LT. 1.D-8 ) CYCLE
         !
         ! ... three dimensional indexes
         !
         i = idx0 + ir - 1
         k = i / (dfftp%nr1x*dfftp%nr2x)
         i = i - (dfftp%nr1x*dfftp%nr2x)*k
         j = i / dfftp%nr1x
         i = i - dfftp%nr1x*j
         !
         i1(0) = i
         i2(0) = j
         i3(0) = k
         !
         ! ... map surrounding grid points
         !
         DO i = -nr1c, nr1c
            i1(i) = i1(0) + i - FLOOR( DBLE( i1(0) + i ) / dfftp%nr1x ) * dfftp%nr1x 
         ENDDO
         DO j = -nr2c, nr2c
            i2(j) = i2(0) + j - FLOOR( DBLE( i2(0) + j ) / dfftp%nr2x ) * dfftp%nr2x 
         ENDDO
         DO k = -nr3c, nr3c
            i3(k) = i3(0) + k - FLOOR( DBLE( i3(0) + k ) / dfftp%nr3x ) * dfftp%nr3x 
         ENDDO
         !
         ! ... integrate
         !
         DO i = -nr1c, nr1c
            itmp = 1 + i1(i)
            DO j = -nr2c, nr2c
               jtmp = itmp + i2(j) * dfftp%nr1x
               DO k = -nr3c, nr3c
                  jr = jtmp + i3(k) * dfftp%nr1x * dfftp%nr2x 
                  fout(ir) = fout(ir) + ftot(jr) * fconv(ABS(i),ABS(j),ABS(k))
               ENDDO
            ENDDO
         ENDDO
         !
      ENDDO
      !
      DEALLOCATE( fconv )
      DEALLOCATE( i1 )
      DEALLOCATE( i2 )
      DEALLOCATE( i3 )
      DEALLOCATE( ftot )
      !
      RETURN
      !
!--------------------------------------------------------------------
    END SUBROUTINE compute_convolution
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
