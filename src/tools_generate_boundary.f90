! Copyright (C) 2018 ENVIRON (www.quantum-environment.org)
!
!    This file is part of Environ version 1.0
!
!    Environ 1.0 is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 2 of the License, or
!    (at your option) any later version.
!
!    Environ 1.0 is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more detail, either the file
!    `License' in the root directory of the present distribution, or
!    online at <http://www.gnu.org/licenses/>.
!
!> This module contains all the procedures to generate the boundary
!! function and its derivative, either as a functional of the electronic
!! density (self-consistent boundary), or as a function of the ionic
!! positions (soft-sphere boundary), or a simple geometric surface
!! centered on the system position (system boundary)
!
! Authors: Oliviero Andreussi (Department of Physics, UNT)
!          Ismaila Dabo       (DMSE, Penn State)
!          Nicola Marzari     (THEOS and NCCR-MARVEL, EPFL)
!
!----------------------------------------------------------------------------
MODULE tools_generate_boundary
!----------------------------------------------------------------------------
  !
  USE environ_types
  USE environ_output
  USE electrostatic_base, ONLY : boundary_core, fd
  !
  PRIVATE
  !
  PUBLIC :: boundary_of_density, boundary_of_functions, boundary_of_system, &
       & invert_boundary, calc_dboundary_dions, solvent_aware_boundary, &
       & solvent_aware_de_dboundary, test_de_dboundary, &
       & compute_ion_field, compute_ion_field_partial, &
       & compute_normal_field, compute_dion_field_drho
  !
CONTAINS
!  Function: sfunct0
!
!> Switching function 0: goes from 1 to 0 when passing through the
!! threshold
!!
!! \f[
!!    1 + \frac{1 - (x/x_t)^k}{1 + (x/x_t)^k}
!! \f]
!! where \f$x_t\f$ is the threshold
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
!  Function: dsfunct0
!
!> Derivative of switching function 0
!--------------------------------------------------------------------
  FUNCTION dsfunct0( x, xthr, fact )
!--------------------------------------------------------------------
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
!  Function: sfunct1
!
!> Switching function 1 that goes from 1 to 0 when passing from
!! xmin to xmax. 
!!
!! NOTE: fact should be equal to LOG(xmax/xmin) but is
!! passed in input to save time
!!
!! \f[
!!    x - \sin(x)
!! \f]
!--------------------------------------------------------------------
  FUNCTION sfunct1( x, xmax, xmin, fact )
!--------------------------------------------------------------------
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
!  Function: dsfunct
!
!> Derivative of switching function 1
!! 
!! NOTE: fact should be equal to LOG(xmax/xmin) but is passed in
!! input to save time.
!--------------------------------------------------------------------
  FUNCTION dsfunct1( x, xmax, xmin, fact )
!--------------------------------------------------------------------
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
!  Function: d2sfunct1
!
!> Second derivative of switching function 1
!!
!! Note: fact should be equal to LOG(xmax/xmin) but is passed in
!! input to save time
!--------------------------------------------------------------------
  FUNCTION d2sfunct1( x, xmax, xmin, fact )
!--------------------------------------------------------------------
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
       d2sfunct1 = ( tpi * SIN( arg ) + fact * ( 1.D0 - COS( arg ) ) ) &
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
!  Function: sfunct2
!
!> Switching function 2, erfc() that goes from 1 to 0 when passing 
!! through xthr. 
!--------------------------------------------------------------------
  FUNCTION sfunct2( x, xthr, spread )
!--------------------------------------------------------------------
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
!  Function: dsfunct2
!
!> Derivative of switching function 2
!--------------------------------------------------------------------
  FUNCTION dsfunct2( x, xthr, spread )
!--------------------------------------------------------------------
    IMPLICIT NONE
    !
    REAL( DP )             :: dsfunct2
    REAL( DP )             :: x, xthr, spread
    !
    ! ... Local variables
    !
    REAL( DP )             :: arg
    !
    arg = ( x - xthr ) / spread
    IF ( abs(arg) .GT. 6.D0 ) THEN ! 6.D0 is the threshold of qe_erfc(x)
       dsfunct2 = 0.D0
    ELSE
       dsfunct2 = - EXP( -arg**2 ) / sqrtpi / spread
    END IF
    !
    RETURN
    !
!--------------------------------------------------------------------
  END FUNCTION dsfunct2
!--------------------------------------------------------------------
!  Function: d2sfunct2
!
!> Second derivative of switching function 2
!--------------------------------------------------------------------
  FUNCTION d2sfunct2( x, xthr, spread )
!--------------------------------------------------------------------
    IMPLICIT NONE
    !
    REAL( DP )             :: d2sfunct2
    REAL( DP )             :: x, xthr, spread
    !
    ! ... Local variables
    !
    REAL( DP )             :: arg
    !
    arg = ( x - xthr ) / spread
    IF ( abs(arg) .GT. 6.D0 ) THEN
       d2sfunct2 = 0.D0
    ELSE
       d2sfunct2 = EXP( -arg**2 ) / sqrtpi / spread**2 * 2.D0 * arg
    END IF
    !
    RETURN
    !
!--------------------------------------------------------------------
  END FUNCTION d2sfunct2
!--------------------------------------------------------------------
!  Function: boundfunct
!
!> Calculates the density-dependent dielectric constant
!! ifunct = 0 => original Fattebert and Gygi function
!--------------------------------------------------------------------
  FUNCTION boundfunct( rho, rhomax, rhomin, tbeta, const, ifunct )
!--------------------------------------------------------------------
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
!  Function: dboundfunct
!
!> Calculates the derivative of the density-dependent dielectric
!! constant
!! ifunct = 0 => original Fattebert and Gygi function
!--------------------------------------------------------------------
  FUNCTION dboundfunct( rho, rhomax, rhomin, tbeta, const, ifunct )
!--------------------------------------------------------------------
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
    CHARACTER( LEN=80 ) :: fun_name = 'dboundfunct'
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
!  Function: d2boundfunct
!
!> Calculates the second derivative of the density-dependent
!! dielectric constant
!!
!! ifunct = 0 => original Fattebery and Gygi function
!--------------------------------------------------------------------
  FUNCTION d2boundfunct( rho, rhomax, rhomin, tbeta, const, ifunct )
!--------------------------------------------------------------------
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
    CHARACTER( LEN=80 ) :: fun_name = 'd2boundfunct'
    !
    SELECT CASE( ifunct )
       !
    CASE( 0 )
       !
       CALL errore(fun_name,'Option not yet implemented',1)
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
!  Subroutine: boundary_of_density
!
!> Calculates the dielectric constant as a function of the charge
!! density, and the derivative of the the dielectric constant
!! with respect to the charge density. 
!--------------------------------------------------------------------
  SUBROUTINE boundary_of_density( density, boundary )
!--------------------------------------------------------------------
    !
    USE tools_fd_gradient, ONLY : calc_fd_gradient
    !
    IMPLICIT NONE
    !
    TYPE( environ_density ), TARGET, INTENT(IN) :: density
    TYPE( environ_boundary ), TARGET, INTENT(INOUT) :: boundary
    !
    INTEGER, POINTER :: ir_end, nnr, stype, deriv
    REAL( DP ), POINTER :: const, rhomax, rhomin, tbeta
    REAL( DP ), DIMENSION(:), POINTER :: rho, eps, deps, d2eps
    REAL( DP ), DIMENSION(:), POINTER :: lapleps, dsurface
    REAL( DP ), DIMENSION(:,:), POINTER :: gradeps
    REAL( DP ), DIMENSION(:,:,:), POINTER :: hesseps
    !
    INTEGER :: ir, ipol, jpol
    !
    CHARACTER( LEN=80 ) :: sub_name = 'boundary_of_density'
    !
    IF ( .NOT. ASSOCIATED(density%cell,boundary%scaled%cell) ) &
         & CALL errore(sub_name,'Inconsistent domains',1)
    !
    ir_end => density % cell % ir_end
    nnr => density % cell % nnr
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
    DO ir = 1, ir_end
       !
       eps( ir )   = boundfunct( rho( ir ), rhomax, rhomin, tbeta, const, stype )
       deps( ir )  = dboundfunct( rho( ir ), rhomax, rhomin, tbeta, const, stype )
       d2eps( ir ) = d2boundfunct( rho( ir ), rhomax, rhomin, tbeta, const, stype )
       !
    END DO
    !
    boundary % volume = integrate_environ_density( boundary % scaled )
    !
    ! ... Compute boundary derivatives, if needed
    !
    deriv => boundary % deriv
    !
    IF ( deriv .GE. 1 ) gradeps => boundary%gradient%of_r
    IF ( deriv .GE. 2 ) lapleps => boundary%laplacian%of_r
    IF ( deriv .GE. 3 ) THEN
       dsurface => boundary%dsurface%of_r
       IF ( boundary % solvent_aware ) THEN
          hesseps => boundary%hessian%of_r
       ELSE
          ALLOCATE( hesseps( 3, 3, nnr ) )
          hesseps = 0.D0
       ENDIF
    END IF
    !
    SELECT CASE ( boundary_core )
       !
    CASE ( 'fft' )
       !
       IF ( deriv .EQ. 1 .OR. deriv .EQ. 2 ) CALL external_gradient( eps, gradeps )
       IF ( deriv .EQ. 2 ) CALL external_laplacian( eps, lapleps )
       IF ( deriv .EQ. 3 ) CALL external_dsurface( nnr, ir_end, eps, gradeps, lapleps, hesseps, dsurface )
       !
    CASE ( 'analytic', 'fd' )
       !
       IF ( deriv .EQ. 1 .OR. deriv .EQ. 2 ) CALL external_gradient( rho , gradeps )
       IF ( deriv .EQ. 2 ) CALL external_laplacian( rho, lapleps )
       IF ( deriv .EQ. 3 ) THEN
          CALL external_dsurface( nnr, ir_end, rho, gradeps, lapleps, hesseps, dsurface )
          IF ( boundary % solvent_aware ) THEN
             DO ipol = 1, 3
                DO jpol = 1, 3
                   hesseps(ipol,jpol,:) = hesseps(ipol,jpol,:) * deps(:) + &
                        & gradeps(ipol,:) * gradeps(jpol,:) * d2eps(:)
                END DO
             END DO
          END IF
       END IF
       IF ( deriv .GT. 1 ) lapleps(:) = lapleps(:) * deps(:) + &
            & ( gradeps(1,:)**2 + gradeps(2,:)**2 + gradeps(3,:)**2 ) * d2eps(:)
       IF ( deriv .GE. 1 ) THEN
          IF ( boundary_core .EQ. 'analytic' ) THEN
             DO ipol = 1, 3
                gradeps(ipol,:) = gradeps(ipol,:) * deps(:)
             ENDDO
          ELSE IF ( boundary_core .EQ. 'fd' ) THEN
             CALL calc_fd_gradient( fd%nfdpoint, fd%icfd, fd%ncfd, nnr, eps, gradeps )
          ENDIF
       ENDIF
       !
    END SELECT
    !
    ! ... Final updates
    !
    IF ( deriv .GE. 1 ) THEN
       CALL update_gradient_modulus( boundary%gradient )
       boundary % surface = integrate_environ_density( boundary%gradient%modulus )
    END IF
    !
    IF ( deriv .GE. 3 .AND. .NOT. boundary % solvent_aware ) &
         & DEALLOCATE( hesseps )
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE boundary_of_density
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE external_dsurface( n, iend, x, grad, lapl, hess, dsurface )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: n, iend
    REAL( DP ), DIMENSION( n ), INTENT(IN) :: x
    REAL( DP ), DIMENSION( 3, n ), INTENT(OUT) :: grad
    REAL( DP ), DIMENSION( n ), INTENT(OUT) :: lapl
    REAL( DP ), DIMENSION( 3, 3, n ), INTENT(OUT) :: hess
    REAL( DP ), DIMENSION( n ), INTENT(OUT) :: dsurface
    !
    CALL external_hessian( x, grad, hess )
    lapl(:) = hess(1,1,:) + hess(2,2,:) + hess(3,3,:)
    CALL calc_dsurface( n, iend, grad, hess, dsurface )
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE external_dsurface
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE calc_dsurface( n, iend, grad, hess, dsurface )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    REAL( DP ), PARAMETER :: toldsurface = 1.D-50
    !
    INTEGER, INTENT(IN) :: n, iend
    REAL( DP ), DIMENSION( 3, n ), INTENT(IN) :: grad
    REAL( DP ), DIMENSION( 3, 3, n ), INTENT(IN) :: hess
    REAL( DP ), DIMENSION( n ), INTENT(OUT) :: dsurface
    !
    INTEGER :: ipol, jpol, i
    REAL( DP ) :: gmod
    !
    DO i = 1, iend
       dsurface(i) = 0.D0
       gmod = SUM( grad(:,i)**2 )
       IF ( gmod .LT. toldsurface ) CYCLE
       DO ipol = 1, 3
          DO jpol = 1,3
             IF ( ipol .EQ. jpol ) CYCLE
             dsurface(i) = dsurface(i) + &
                  grad(ipol, i) * grad(jpol, i) * hess(ipol, jpol, i) - &
                  grad(ipol, i) * grad(ipol, i) * hess(jpol, jpol, i)
          ENDDO
       ENDDO
       dsurface( i ) = dsurface( i ) / gmod / SQRT(gmod)
    END DO
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE calc_dsurface
!--------------------------------------------------------------------
!  Subroutine: boundary_of_functions
!
!> Calculates the dielectric constant as a function of the charge
!! density, and derivative of the dielectric constant with respect
!! to the charge density
!--------------------------------------------------------------------
  SUBROUTINE boundary_of_functions( nsoft_spheres, soft_spheres, boundary )
!--------------------------------------------------------------------
    USE utils_functions, ONLY: density_of_functions, gradient_of_functions, &
         & laplacian_of_functions, hessian_of_functions
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: nsoft_spheres
    TYPE( environ_functions ), DIMENSION(nsoft_spheres), INTENT(IN) :: soft_spheres
    TYPE( environ_boundary ), TARGET, INTENT(INOUT) :: boundary
    !
    INTEGER, POINTER :: nnr, ir_end, deriv
    TYPE( environ_cell ), POINTER :: cell
    !
    INTEGER :: i
    CHARACTER( LEN=80 ) :: sub_name = 'boundary_of_functions'
    !
    TYPE( environ_density ), DIMENSION(:), ALLOCATABLE :: local
    TYPE( environ_gradient ), DIMENSION(:), ALLOCATABLE :: gradlocal
    TYPE( environ_density ), DIMENSION(:), ALLOCATABLE :: lapllocal
    TYPE( environ_hessian ), DIMENSION(:), ALLOCATABLE :: hesslocal
    TYPE( environ_hessian ), POINTER :: hessian
    !
    ! Aliases and allocations
    !
    cell => boundary % scaled % cell
    nnr => cell % nnr
    ir_end => cell % ir_end
    !
    ALLOCATE( local( nsoft_spheres ) )
    !
    ! Compute soft spheres and generate boundary
    !
    boundary % scaled % of_r = 1.D0
    !
    DO i = 1, nsoft_spheres
       !
       CALL init_environ_density( cell, local(i) )
       !
       CALL density_of_functions( soft_spheres(i), local(i), .FALSE. )
       !
       boundary % scaled % of_r = boundary % scaled % of_r * local(i) % of_r
       !
    ENDDO
    !
    ! Generate boundary derivatives, if needed
    !
    deriv => boundary % deriv
    !
    IF ( deriv .EQ. 3 ) THEN
       !
       IF ( boundary%solvent_aware ) THEN
          !
          hessian => boundary%hessian
          !
       ELSE
          !
          ALLOCATE( hessian )
          CALL init_environ_hessian( cell, hessian )
          !
       END IF
       !
    END IF
    !
    SELECT CASE ( boundary_core )
       !
    CASE ( 'fft' )
       !
       IF ( deriv .EQ. 1 .OR. deriv .EQ. 2 ) CALL external_gradient( boundary%scaled%of_r, &
            & boundary%gradient%of_r )
       IF ( deriv .EQ. 2 ) CALL external_laplacian( boundary%scaled%of_r, boundary%laplacian%of_r )
       IF ( deriv .EQ. 3 ) CALL external_dsurface( nnr, ir_end, boundary%scaled%of_r, &
            & boundary%gradient%of_r, boundary%laplacian%of_r, hessian%of_r, boundary%dsurface%of_r )
       !
    CASE ( 'analytic', 'highmem' )
       !
       IF ( deriv .GE. 1 ) ALLOCATE( gradlocal( nsoft_spheres ) )
       IF ( deriv .EQ. 2 ) ALLOCATE( lapllocal( nsoft_spheres ) )
       IF ( deriv .EQ. 3 ) ALLOCATE( hesslocal( nsoft_spheres ) )
       !
       ! Compute and temporarily store soft spheres derivatives
       !
       DO i = 1, nsoft_spheres
          !
          IF ( deriv .GE. 1 ) CALL init_environ_gradient( cell, gradlocal(i) )
          IF ( deriv .EQ. 2 ) CALL init_environ_density( cell, lapllocal(i) )
          IF ( deriv .EQ. 3 ) CALL init_environ_hessian( cell, hesslocal(i) )
          !
          IF ( deriv .GE. 1 ) CALL gradient_of_functions( soft_spheres(i), gradlocal(i), .FALSE. )
          IF ( deriv .EQ. 2 ) CALL laplacian_of_functions( soft_spheres(i), lapllocal(i), .FALSE. )
          IF ( deriv .EQ. 3 ) CALL hessian_of_functions( soft_spheres(i), hesslocal(i), .FALSE. )
          !
       END DO
       !
       IF ( deriv .EQ. 1 .OR. deriv .EQ. 2 ) &
            & CALL calc_gradient_of_boundary_highmem(nsoft_spheres,local, &
            & gradlocal,boundary%gradient )
       IF ( deriv .EQ. 2 ) CALL calc_laplacian_of_boundary_highmem(nsoft_spheres, &
            & local,gradlocal,lapllocal,boundary%laplacian )
       IF ( deriv .EQ. 3 ) CALL calc_dsurface_of_boundary_highmem(nsoft_spheres, &
            & local,gradlocal,hesslocal,boundary%gradient,boundary%laplacian,hessian, &
            & boundary%dsurface )
       !
       ! Destroy and deallocate
       !
       DO i = 1, nsoft_spheres
          IF ( deriv .GE. 1 ) CALL destroy_environ_gradient( gradlocal(i) )
          IF ( deriv .EQ. 2 ) CALL destroy_environ_density( lapllocal(i) )
          IF ( deriv .EQ. 3 ) CALL destroy_environ_hessian( hesslocal(i) )
       ENDDO
       IF ( deriv .GE. 1 ) DEALLOCATE( gradlocal )
       IF ( deriv .EQ. 2 ) DEALLOCATE( lapllocal )
       IF ( deriv .EQ. 3 ) DEALLOCATE( hesslocal )
       !
    END SELECT
    !
    ! Final updates
    !
    boundary % scaled % of_r = 1.D0 - boundary % scaled % of_r
    !
    boundary % volume = integrate_environ_density( boundary % scaled )
    !
    IF ( deriv .GE. 1 ) THEN
       !
       boundary % gradient % of_r = - boundary % gradient % of_r
       !
       CALL update_gradient_modulus( boundary%gradient )
       !
       boundary % surface = integrate_environ_density( boundary%gradient%modulus )
       !
       IF ( deriv .GE. 2 ) boundary % laplacian % of_r = - boundary % laplacian % of_r
       !
       IF ( deriv .EQ. 3 ) THEN
          !
          boundary % dsurface % of_r = - boundary % dsurface % of_r
          !
          IF ( boundary % solvent_aware ) THEN
             !
             boundary % hessian % of_r = - boundary % hessian % of_r
             !
          ELSE
             !
             CALL destroy_environ_hessian( hessian )
             DEALLOCATE( hessian )
             !
          ENDIF
          !
       END IF
       !
    END IF
    !
    DO i = 1, nsoft_spheres
       CALL destroy_environ_density( local(i) )
    ENDDO
    DEALLOCATE( local )
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE boundary_of_functions
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE calc_partial_of_boundary_highmem(n,i,local,gradlocal,partial)
!--------------------------------------------------------------------
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: n, i
    TYPE( environ_density ), DIMENSION(n), INTENT(IN) :: local
    TYPE( environ_gradient ), DIMENSION(n), INTENT(IN) :: gradlocal
    TYPE( environ_gradient ), INTENT(INOUT) :: partial
    !
    INTEGER :: j, ipol
    CHARACTER( LEN=80 ) :: sub_name = 'calc_partial_of_boundary_highmem'
    !
    IF ( i .GT. n ) CALL errore(sub_name,'Index out of bound',1)
    !
    DO ipol = 1, 3
       partial % of_r( ipol, : ) = gradlocal(i) % of_r( ipol, : )
       DO j = 1, n
          IF ( j .EQ. i ) CYCLE
          partial % of_r( ipol, : ) = partial % of_r( ipol, : ) * local(j) % of_r( : )
       ENDDO
    ENDDO
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE calc_partial_of_boundary_highmem
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE calc_gradient_of_boundary_highmem(n,local,gradlocal,gradient)
!--------------------------------------------------------------------
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: n
    TYPE( environ_density ), DIMENSION(n), INTENT(IN) :: local
    TYPE( environ_gradient ), DIMENSION(n), INTENT(IN) :: gradlocal
    TYPE( environ_gradient ), INTENT(INOUT) :: gradient
    !
    INTEGER :: i, j, ipol
    TYPE( environ_cell ), POINTER :: cell
    TYPE( environ_gradient ) :: partial
    !
    cell => gradient % cell
    !
    CALL init_environ_gradient( cell, partial )
    !
    gradient % of_r = 0.D0
    DO i = 1, n
       CALL calc_partial_of_boundary_highmem(n,i,local,gradlocal,partial)
       gradient % of_r  = gradient % of_r + partial % of_r
    ENDDO
    !
    CALL destroy_environ_gradient( partial )
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE calc_gradient_of_boundary_highmem
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE calc_laplacian_of_boundary_highmem(n,local,gradlocal,lapllocal,laplacian)
!--------------------------------------------------------------------
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: n
    TYPE( environ_density ), DIMENSION(n), INTENT(IN) :: local
    TYPE( environ_gradient ), DIMENSION(n), INTENT(IN) :: gradlocal
    TYPE( environ_density ), DIMENSION(n), INTENT(IN) :: lapllocal
    TYPE( environ_density ), INTENT(INOUT) :: laplacian
    !
    INTEGER :: i, j, k, ipol
    TYPE( environ_cell), POINTER :: cell
    TYPE( environ_density ) :: tmp
    !
    cell => laplacian % cell
    CALL init_environ_density( cell, tmp )
    !
    laplacian % of_r = 0.D0
    DO i = 1, n
       DO j = 1, n
          IF ( j .EQ. i ) THEN
             tmp % of_r = lapllocal(i) % of_r
          ELSE
             CALL scalar_product_environ_gradient( gradlocal(i), gradlocal(j), tmp )
          ENDIF
          DO k = 1, n
             IF ( k .EQ. j .OR. k .EQ. i ) CYCLE
             tmp % of_r = tmp % of_r * local(k) % of_r
          ENDDO
          laplacian % of_r = laplacian % of_r + tmp % of_r
       ENDDO
    ENDDO
    !
    CALL destroy_environ_density( tmp )
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE calc_laplacian_of_boundary_highmem
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE calc_dsurface_of_boundary_highmem(n,local,gradlocal,hesslocal,gradient,laplacian,hessian,dsurface)
!--------------------------------------------------------------------
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: n
    TYPE( environ_density ), DIMENSION(n), INTENT(IN) :: local
    TYPE( environ_gradient ), DIMENSION(n), INTENT(IN) :: gradlocal
    TYPE( environ_hessian ), DIMENSION(n), INTENT(IN) :: hesslocal
    TYPE( environ_gradient ), INTENT(INOUT) :: gradient
    TYPE( environ_density ), INTENT(INOUT) :: laplacian, dsurface
    TYPE( environ_hessian ), INTENT(INOUT) :: hessian
    !
    INTEGER :: i, j, k, ipol, jpol
    TYPE( environ_cell), POINTER :: cell
    TYPE( environ_density ) :: dens
    TYPE( environ_gradient ) :: partial
    !
    cell => laplacian % cell
    CALL init_environ_density( cell, dens )
    CALL init_environ_gradient( cell, partial )
    !
    gradient % of_r = 0.D0
    DO i = 1, n
       CALL calc_partial_of_boundary_highmem( n, i, local, gradlocal, partial )
       gradient % of_r = gradient % of_r + partial % of_r
       DO j = 1, n
          DO ipol = 1, 3
             DO jpol = 1, 3
                IF ( j .EQ. i ) THEN
                   dens % of_r( : ) = hesslocal(i) % of_r( ipol, jpol, : )
                ELSE
                   dens % of_r( : ) = gradlocal(i) % of_r( ipol, : ) * gradlocal(j) % of_r( jpol, : )
                ENDIF
                DO k = 1, n
                   IF ( k .EQ. j .OR. k .EQ. i ) CYCLE
                   dens % of_r = dens % of_r * local(k) % of_r
                ENDDO
                hessian % of_r( ipol, jpol, : ) = hessian % of_r( ipol, jpol, : ) + dens % of_r( : )
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    !
    ! Final operations
    !
    laplacian % of_r = hessian % of_r( 1, 1, : ) + hessian % of_r( 2, 2, : ) + hessian % of_r( 3, 3, : )
    CALL calc_dsurface( cell%nnr, cell%ir_end, gradient%of_r, hessian%of_r, dsurface%of_r )
    !
    CALL destroy_environ_density( dens )
    CALL destroy_environ_gradient( partial )
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE calc_dsurface_of_boundary_highmem
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE calc_dboundary_dions( index, boundary, partial )
!--------------------------------------------------------------------
    !
    USE utils_functions, ONLY: density_of_functions, gradient_of_functions
    !
    IMPLICIT NONE
    !
    REAL( DP ), PARAMETER :: tolspuriousforce = 1.D-5
    !
    INTEGER, INTENT(IN) :: index
    TYPE( environ_boundary ), INTENT(IN), TARGET :: boundary
    TYPE( environ_gradient ), INTENT(INOUT) :: partial
    !
    INTEGER, POINTER :: number
    TYPE( environ_cell ), POINTER :: cell
    !
    INTEGER :: i, ipol
    REAL( DP ) :: spurious_force
    TYPE( environ_density ) :: local
    !
    CHARACTER( LEN=80 ) :: sub_name='calc_dboundary_dions'
    !
    ! ... Exit if boundary is only defined on electronic density
    !
    IF ( boundary % mode .EQ. 'electronic' ) RETURN
    !
    ! ... Aliases and sanity checks
    !
    cell => partial % cell
    IF ( boundary % need_ions ) THEN
       number => boundary % ions % number
    ELSE IF ( boundary % need_system ) THEN
       number => boundary % system % ions % number
    ELSE
       CALL errore(sub_name,'Missing details of ions',1)
    END IF
    !
    IF ( index .GT. number ) &
         & CALL errore(sub_name,'Index greater than number of ions',1)
    IF ( index .LE. 0 ) &
         & CALL errore(sub_name,'Index of ion is zero or lower',1)
    IF ( boundary % mode .EQ. 'ionic' .AND. .NOT. ALLOCATED( boundary % soft_spheres ) ) &
         & CALL errore(sub_name,'Missing details of ionic boundary',1)
    IF ( boundary % mode .EQ. 'full' .AND. .NOT. ALLOCATED( boundary % ions % core_electrons ) ) &
         & CALL errore(sub_name,'Missing details of core electrons',1)
    IF ( boundary % mode .EQ. 'full' .AND. .NOT. ASSOCIATED( boundary % dscaled % cell, cell ) ) &
         & CALL errore(sub_name,'Mismatch or unassociated boundary derivative',1)
    !
    IF ( boundary % mode .EQ. 'ionic' ) THEN
       !
       SELECT CASE ( boundary_core )
          !
       CASE( 'fft', 'fd', 'analytic', 'highmem' )
          !
          CALL gradient_of_functions( boundary%soft_spheres(index), partial, .TRUE. )
          !
          CALL init_environ_density( cell, local )
          !
          DO i = 1, number
             IF ( i .EQ. index ) CYCLE
             CALL density_of_functions( boundary%soft_spheres(i), local, .TRUE. )
             DO ipol = 1, 3
                partial % of_r( ipol, : ) = partial % of_r( ipol, : ) * local % of_r( : )
             ENDDO
          ENDDO
          !
          CALL destroy_environ_density( local )
          !
       END SELECT
       !
    ELSE IF ( boundary % mode .EQ. 'full' ) THEN
       !
       CALL gradient_of_functions( boundary%ions%core_electrons(index), partial, .TRUE. )
       DO ipol = 1, 3
          partial % of_r( ipol, : ) = - partial % of_r( ipol, : ) * boundary % dscaled % of_r( : )
       END DO
       CALL update_gradient_modulus( partial )
       spurious_force = integrate_environ_density( partial%modulus )
       IF ( spurious_force .GT. tolspuriousforce .AND. ionode ) WRITE( program_unit, 4001 )index, spurious_force
4001   FORMAT(1x,'WARNING: Unphysical forces due to core electrons are non-negligible '&
            /,1x,'atom type ',I3,' is subject to a spurious force of ',F12.6,' ')
       !
    ELSE IF ( boundary % mode .EQ. 'system' ) THEN
       !
       ! PROBABLY THERE IS A UNIFORM CONTRIBUTION TO THE FORCES
       ! WHICH SHOULD ONLY AFFECT THE COM OF THE SYSTEM, POSSIBLY NEED TO ADD
       ! A CHECK ON ATOMS THAT BELONG TO THE SYSTEM
       partial % of_r = 0.D0
       !
    END IF
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE calc_dboundary_dions
!--------------------------------------------------------------------
!  Subroutine: boundary_of_system
!
!> Calculates the dielectric constant as a function of the charge
!! density, and the derivative of the dielectric constant with
!! respect to the charge density.
!--------------------------------------------------------------------
  SUBROUTINE boundary_of_system( simple, boundary )
!--------------------------------------------------------------------
    USE utils_functions, ONLY: density_of_functions, gradient_of_functions, &
         & laplacian_of_functions, hessian_of_functions
    !
    IMPLICIT NONE
    !
    TYPE( environ_functions ), INTENT(IN) :: simple
    TYPE( environ_boundary ), TARGET, INTENT(INOUT) :: boundary
    !
    INTEGER, POINTER :: nnr, ir_end, deriv
    TYPE( environ_cell ), POINTER :: cell
    !
    INTEGER :: i
    CHARACTER( LEN=80 ) :: sub_name = 'boundary_of_system'
    !
    TYPE( environ_hessian ), POINTER :: hesslocal
    !
    ! Aliases and allocations
    !
    cell => boundary % scaled % cell
    nnr => cell % nnr
    ir_end => cell % ir_end
    !
    ! Compute soft spheres and generate boundary
    !
    CALL density_of_functions( simple, boundary % scaled, .TRUE. )
    !
    ! Generate boundary derivatives, if needed
    !
    deriv => boundary % deriv
    !
    IF ( deriv .GE. 3 ) THEN
       IF ( boundary % solvent_aware ) THEN
          hesslocal => boundary % hessian
       ELSE
          ALLOCATE( hesslocal )
          CALL init_environ_hessian( cell, hesslocal )
       ENDIF
    ENDIF
    !
    SELECT CASE ( boundary_core )
       !
    CASE ( 'fft' )
       !
       IF ( deriv .EQ. 1 .OR. deriv .EQ. 2 ) CALL external_gradient( boundary%scaled%of_r, boundary%gradient%of_r )
       IF ( deriv .EQ. 2 ) CALL external_laplacian( boundary%scaled%of_r, boundary%laplacian%of_r )
       IF ( deriv .EQ. 3 ) CALL external_dsurface( nnr, ir_end, boundary%scaled%of_r, boundary%gradient%of_r, &
            & boundary%laplacian%of_r, hesslocal%of_r, boundary%dsurface%of_r )
       !
    CASE ( 'analytic' )
       !
       IF ( deriv .GE. 1 ) CALL gradient_of_functions( simple, boundary%gradient, .TRUE. )
       IF ( deriv .GE. 2 ) CALL laplacian_of_functions( simple, boundary%laplacian, .TRUE. )
       IF ( deriv .GE. 3 ) THEN
          CALL hessian_of_functions( simple, hesslocal, .TRUE. )
          CALL calc_dsurface( nnr, ir_end, boundary%gradient%of_r, hesslocal%of_r, boundary%dsurface%of_r )
       ENDIF
       !
    END SELECT
    !
    IF ( deriv .GE. 3 ) THEN
       IF ( .NOT. boundary % solvent_aware ) THEN
          CALL destroy_environ_hessian( hesslocal )
          DEALLOCATE(hesslocal )
       END IF
    END IF
    !
    boundary % volume = integrate_environ_density( boundary % scaled )
    !
    IF ( deriv .GE. 1 ) THEN
       !
       CALL update_gradient_modulus( boundary%gradient )
       !
       boundary % surface = integrate_environ_density( boundary%gradient%modulus )
       !
    END IF
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE boundary_of_system
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE invert_boundary( boundary )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    TYPE( environ_boundary ), INTENT(INOUT) :: boundary
    !
    boundary % scaled % of_r = 1.D0 - boundary % scaled % of_r
    !
    boundary % volume = integrate_environ_density( boundary % scaled )
    !
    IF ( boundary % deriv .GE. 1 ) &
         & boundary % gradient % of_r = - boundary % gradient % of_r
    !
    IF ( boundary % deriv .GE. 2 ) &
         & boundary % laplacian % of_r = - boundary % laplacian % of_r
    !
    IF ( boundary % deriv .GE. 3 ) THEN
       boundary % dsurface % of_r = - boundary % dsurface % of_r
       IF ( boundary % solvent_aware ) &
            & boundary % hessian % of_r = - boundary % hessian % of_r
    END IF
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE invert_boundary
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE solvent_aware_boundary( boundary )
!--------------------------------------------------------------------
    !
    ! ... Fill voids of the continuum interface that are too small
    ! ... to fit a solvent molecule
    !
    USE utils_functions, ONLY: density_of_functions
    USE tools_generate_functions, ONLY : compute_convolution_fft
    !
    IMPLICIT NONE
    !
    TYPE( environ_boundary ), INTENT(INOUT), TARGET :: boundary
    !
    ! Local variables
    !
    INTEGER, POINTER :: nnr, ir_end, deriv
    REAL( DP ), POINTER :: thr, spr
    TYPE( environ_cell ), POINTER :: cell
    !
    INTEGER :: ir, ipol, jpol
!    REAL( DP ) :: probe_volume
    TYPE( environ_density ) :: filled_fraction
    TYPE( environ_density ) :: d2filling
    !
    TYPE( environ_density ) :: local
    TYPE( environ_gradient ) :: gradlocal
    TYPE( environ_density ) :: lapllocal
    TYPE( environ_hessian ) :: hesslocal
    !
    CHARACTER( LEN = 80 ) :: label
    !
    ! Aliases and sanity checks
    !
    cell => boundary % scaled % cell
    nnr => boundary % scaled % cell % nnr
    ir_end => boundary % scaled % cell % ir_end
    deriv => boundary % deriv
    !
    thr => boundary % filling_threshold
    spr => boundary % filling_spread
    !
    CALL init_environ_density( cell, filled_fraction )
    IF ( deriv .GE. 2 .AND. boundary_core .NE. 'fft' ) CALL init_environ_density( cell, d2filling )
    !
    ! Step 0: save local interface function for later use
    !
    boundary % local % of_r = boundary % scaled % of_r
    !
    ! Step 1: compute the convolution function, this may be made moved out of here
    !
    CALL density_of_functions( boundary%solvent_probe, boundary%probe, .TRUE. )
    !
!    probe_volume = integrate_environ_density( boundary%probe )
!    boundary%probe%of_r = boundary%probe%of_r / probe_volume
    boundary%probe%of_r = boundary%probe%of_r / integrate_environ_density( boundary%probe )
    !
    ! Step 2: compute filled fraction, i.e. convolution of local boundary with probe
    !
    CALL compute_convolution_fft( nnr, boundary%local%of_r, boundary%probe%of_r, filled_fraction%of_r )
    !
    ! Step 3: compute the filling function and its derivative
    !
    boundary % filling % of_r = 0.D0
    boundary % dfilling % of_r = 0.D0
    DO ir = 1, ir_end
       boundary % filling % of_r( ir ) = 1.D0 - sfunct2( filled_fraction % of_r( ir ), thr, spr )
       boundary % dfilling % of_r( ir ) = - dsfunct2( filled_fraction % of_r( ir ), thr, spr )
       IF ( deriv .GE. 2 .AND. boundary_core .NE. 'fft' ) &
            & d2filling % of_r( ir ) = - d2sfunct2( filled_fraction % of_r( ir ), thr, spr )
    ENDDO
    !
    ! Step 4: compute solvent-aware interface
    !
    boundary % scaled % of_r = boundary % local % of_r + &
         & ( 1.D0 - boundary % local % of_r ) * boundary % filling % of_r
    !
    ! Step 5: compute boundary derivatives, if needed
    !
    SELECT CASE ( boundary_core )
       !
    CASE ( 'fft' )
       !
       IF ( deriv .EQ. 1 .OR. deriv .EQ. 2 ) CALL external_gradient( boundary%scaled%of_r, &
            & boundary%gradient%of_r )
       IF ( deriv .EQ. 2 ) CALL external_laplacian( boundary%scaled%of_r, boundary%laplacian%of_r )
       IF ( deriv .EQ. 3 ) CALL external_dsurface( nnr, ir_end, boundary%scaled%of_r, &
            & boundary%gradient%of_r, boundary%laplacian%of_r, boundary%hessian%of_r, boundary%dsurface%of_r )
       !
    CASE( 'analytic', 'highmem' )
       !
       ! Allocate local fields for derivatives of convolution
       !
       IF ( deriv .GE. 1 ) CALL init_environ_gradient( cell, gradlocal )
       IF ( deriv .GE. 2 ) CALL init_environ_density( cell, lapllocal )
       IF ( deriv .GE. 3 ) CALL init_environ_hessian( cell, hesslocal )
       !
       ! Compute derivative of convolution with probe
       !
!       IF ( deriv .GE. 1 ) CALL compute_convolution_deriv( deriv, boundary%solvent_probe, &
!            & boundary%local, gradlocal, lapllocal, hesslocal, probe_volume )
       IF ( deriv .GE. 1 ) CALL compute_convolution_deriv( deriv, boundary, gradlocal, lapllocal, hesslocal )
       !
       ! Update derivatives of interface function in reverse order
       !
       IF ( deriv .GE. 3 ) THEN
          !
          DO ipol = 1, 3
             DO jpol = 1, 3
                !
                boundary % hessian % of_r( ipol, jpol, : ) = boundary % hessian % of_r( ipol, jpol, : ) * &
                     & ( 1.D0 - boundary % filling % of_r ) - boundary % dfilling % of_r * &
                     & ( boundary % gradient % of_r( ipol, : ) * gradlocal % of_r ( jpol, : ) + &
                     &   boundary % gradient % of_r( jpol, : ) * gradlocal % of_r ( ipol, : )  ) + &
                     & ( 1.D0 - boundary % local % of_r ) * &
                     & ( d2filling % of_r * gradlocal % of_r ( ipol, : ) * gradlocal % of_r ( jpol, : ) + &
                     &   boundary % dfilling % of_r * hesslocal % of_r( ipol, jpol, : ) )
                !
             ENDDO
          ENDDO
          !
          CALL destroy_environ_hessian( hesslocal )
          !
       ENDIF
       !
       IF ( deriv .GE. 2 ) THEN
          !
          CALL init_environ_density( cell, local )
          CALL scalar_product_environ_gradient( boundary%gradient, gradlocal, local )
          !
          boundary % laplacian % of_r = boundary % laplacian % of_r * ( 1.D0 - boundary % filling % of_r ) &
               & - 2.D0 * local % of_r * boundary % dfilling % of_r + ( 1.D0 - boundary % local % of_r ) * &
               & ( d2filling % of_r * gradlocal % modulus % of_r **2 + boundary % dfilling % of_r * lapllocal % of_r )
          !
          CALL destroy_environ_density( local )
          CALL destroy_environ_density( lapllocal )
          CALL destroy_environ_density( d2filling )
          !
       ENDIF
       !
       IF ( deriv .GE. 1 ) THEN
          !
          DO ipol = 1, 3
             !
             boundary % gradient % of_r(ipol,:) = boundary % gradient % of_r(ipol,:) * &
                  & ( 1.D0 - boundary % filling % of_r(:) ) + gradlocal%of_r(ipol,:) * &
                  & ( 1.D0 - boundary % local % of_r(:) ) * boundary % dfilling % of_r(:)
             !
          ENDDO
          !
          CALL destroy_environ_gradient( gradlocal )
          !
       ENDIF
       !
       ! if needed, now can recompute dsurface
       !
       IF ( deriv .GE. 3 ) THEN
          !
          CALL calc_dsurface( nnr, ir_end, boundary % gradient % of_r, &
               & boundary % hessian % of_r, boundary % dsurface % of_r )
          !
       END IF
       !
    END SELECT
    !
    ! Final updates
    !
    boundary % volume = integrate_environ_density( boundary % scaled )
    !
    IF ( deriv .GE. 1 ) THEN
       !
       CALL update_gradient_modulus( boundary % gradient )
       !
       boundary % surface = integrate_environ_density( boundary%gradient%modulus )
       !
    END IF
    !
    CALL destroy_environ_density( filled_fraction )
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE solvent_aware_boundary
!--------------------------------------------------------------------
!  Subroutine: solvent_aware_de_dboundary
!
!> Fill voids of the continuum interface that are too small to
!! fit a solvent molecule
!--------------------------------------------------------------------
  SUBROUTINE solvent_aware_de_dboundary( boundary, de_dboundary )
!--------------------------------------------------------------------
    USE tools_generate_functions, ONLY : compute_convolution_fft
    !
    IMPLICIT NONE
    !
    TYPE( environ_boundary ), INTENT(IN), TARGET :: boundary
    TYPE( environ_density ), INTENT(INOUT) :: de_dboundary
    !
    ! Local variables
    !
    INTEGER, POINTER :: nnr, ir_end
    REAL( DP ), POINTER :: thr, spr
    TYPE( environ_cell ), POINTER :: cell
    !
    TYPE( environ_density ) :: local
    !
    ! Aliases and sanity checks
    !
    cell => boundary % scaled % cell
    nnr => boundary % scaled % cell % nnr
    !
    CALL init_environ_density( cell, local )
    !
    ! Step 1: compute (1-s)*de_dboudary*dfilling
    !
    local % of_r = ( 1.D0 - boundary % local % of_r ) * de_dboundary % of_r * boundary % dfilling % of_r
    !
    ! Step 2: compute convolution with the probe function
    !
    CALL compute_convolution_fft( nnr, boundary%probe%of_r, local%of_r, local%of_r )
    !
    ! Step 3: update the functional derivative of the energy wrt boundary
    !
    de_dboundary % of_r = de_dboundary % of_r * ( 1.D0 - boundary % filling % of_r ) + &
         & local % of_r
    !
    CALL destroy_environ_density( local )
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE solvent_aware_de_dboundary
!--------------------------------------------------------------------
!!--------------------------------------------------------------------
!  SUBROUTINE compute_convolution_deriv( deriv, probe, f, grad, lapl, hess, probe_vol )
!!--------------------------------------------------------------------
!    !
!    USE utils_functions, ONLY : gradient_of_functions, laplacian_of_functions, hessian_of_functions
!    USE tools_generate_functions, ONLY : compute_convolution_fft
!    !
!    IMPLICIT NONE
!    !
!    INTEGER, INTENT(IN) :: deriv
!    REAL( DP ), INTENT( IN ) :: probe_vol
!    TYPE( environ_functions ), INTENT(IN) :: probe
!    TYPE( environ_density ), INTENT(IN) :: f
!    TYPE( environ_gradient ), INTENT(INOUT) :: grad
!    TYPE( environ_density ), INTENT(INOUT) :: lapl
!    TYPE( environ_hessian ), INTENT(INOUT) :: hess
!    !
!    INTEGER, POINTER :: nnr
!    !
!    INTEGER :: ipol, jpol
!    !
!    nnr => f % cell % nnr
!    !
!    IF ( deriv .LE. 0 ) RETURN
!    !
!    IF ( deriv .GE. 1 ) THEN
!       !
!       CALL gradient_of_functions( probe, grad, .FALSE. )
!       grad%of_r(:,:) = grad%of_r(:,:) / probe_vol
!       !
!       DO ipol = 1, 3
!          CALL compute_convolution_fft( nnr, f%of_r, grad%of_r(ipol,:), grad%of_r(ipol,:))
!       ENDDO
!       !
!       CALL update_gradient_modulus( grad )
!       !
!    ENDIF
!    !
!    IF ( deriv .GE. 2 ) THEN
!       !
!       CALL laplacian_of_functions( probe, lapl, .FALSE. )
!       lapl%of_r = lapl%of_r / probe_vol
!       !
!       CALL compute_convolution_fft( nnr, f%of_r, lapl%of_r, lapl%of_r )
!       !
!    END IF
!    !
!    IF ( deriv .GE. 3 ) THEN
!       !
!       CALL hessian_of_functions( probe, hess, .FALSE. )
!       hess%of_r(:,:,:) = hess%of_r(:,:,:) / probe_vol
!       !
!       DO ipol = 1, 3
!          DO jpol = 1, 3
!             CALL compute_convolution_fft( nnr, f%of_r, hess%of_r(ipol,jpol,:), &
!                  & hess%of_r(ipol,jpol,:) )
!          ENDDO
!       ENDDO
!       !
!    ENDIF
!    !
!    RETURN
!    !
!!--------------------------------------------------------------------
!  END SUBROUTINE compute_convolution_deriv
!!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE compute_convolution_deriv( deriv, bound, grad, lapl, hess )
!--------------------------------------------------------------------
    !
    USE tools_generate_functions, ONLY : compute_convolution_fft
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: deriv
    TYPE( environ_boundary ), INTENT(IN) :: bound
    TYPE( environ_gradient ), INTENT(INOUT) :: grad
    TYPE( environ_density ), INTENT(INOUT) :: lapl
    TYPE( environ_hessian ), INTENT(INOUT) :: hess
    !
    INTEGER, POINTER :: nnr
    !
    INTEGER :: ipol, jpol
    !
    nnr => bound % probe % cell % nnr
    !
    IF ( deriv .LE. 0 ) RETURN
    !
    IF ( deriv .GE. 1 ) THEN
       !
       DO ipol = 1, 3
          CALL compute_convolution_fft( nnr, bound%probe%of_r, bound%gradient%of_r(ipol,:), &
               & grad%of_r(ipol,:))
       ENDDO
       !
       CALL update_gradient_modulus( grad )
       !
    ENDIF
    !
    IF ( deriv .GE. 2 ) THEN
       !
       CALL compute_convolution_fft( nnr, bound%probe%of_r, bound%laplacian%of_r, lapl%of_r )
       !
    END IF
    !
    IF ( deriv .GE. 3 ) THEN
       !
       DO ipol = 1, 3
          DO jpol = 1, 3
             CALL compute_convolution_fft( nnr, bound%probe%of_r, bound%hessian%of_r(ipol,jpol,:), &
                  & hess%of_r(ipol,jpol,:) )
          ENDDO
       ENDDO
       !
    ENDIF
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE compute_convolution_deriv
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE compute_ion_field( nsoft_spheres, soft_spheres, ions, electrons, ion_field )
!--------------------------------------------------------------------
    !
    USE utils_functions, ONLY : density_of_functions, gradient_of_functions
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: nsoft_spheres
    TYPE( environ_functions), DIMENSION(nsoft_spheres), INTENT(IN) :: soft_spheres
    TYPE( environ_ions ), INTENT(IN) :: ions
    TYPE( environ_electrons ), INTENT(IN) :: electrons
    REAL( DP ), DIMENSION(nsoft_spheres), INTENT(OUT) :: ion_field
    !
    TYPE( environ_cell ), POINTER :: cell
    !
    INTEGER :: i, j
    CHARACTER( LEN=80 ) :: sub_name = 'compute_ion_field'
    !
    TYPE( environ_density ), DIMENSION(:), ALLOCATABLE :: local
    !
    TYPE( environ_density ) :: aux, prod
    TYPE( environ_gradient ) :: gradaux, field
    !
    ! Aliases and allocations
    !
    cell => ions % density % cell
    !
    ALLOCATE( local( nsoft_spheres ) )
    !
    ! Compute field-independent soft spheres and gradients
    !
    DO i = 1, nsoft_spheres
       !
       CALL init_environ_density( cell, local(i) )
       !
       CALL density_of_functions( soft_spheres(i), local(i), .FALSE. )
       !
    ENDDO  
    !
    ! Compute field
    !
    CALL init_environ_density( cell, aux )
    aux % of_r = electrons % density % of_r + ions % density % of_r
    !
    CALL init_environ_gradient( cell, field )
    CALL gradv_h_of_rho_r( aux%of_r, field%of_r )
    !
    ! Compute field flux
    !
    ion_field = 0.d0
    !
    CALL init_environ_density( cell, prod )
    CALL init_environ_gradient( cell, gradaux )
    !
    DO i = 1, nsoft_spheres
       !
       ! Compute product of other soft-spheres
       !
       prod % of_r = 1.D0
       DO j = 1, nsoft_spheres
          IF ( j .EQ. i ) CYCLE
          prod % of_r = prod % of_r * local(j) % of_r
       ENDDO
       !
       ! Compute field flux through soft-sphere interface
       !
       CALL gradient_of_functions( soft_spheres(i), gradaux, .TRUE. )
       !
       CALL scalar_product_environ_gradient( field, gradaux, aux ) ! here aux is the normal field
       !
       aux % of_r = - aux % of_r * prod % of_r
       !
       ion_field(i) = integrate_environ_density( aux )
       !
    END DO
    !
    CALL destroy_environ_gradient( gradaux )
    CALL destroy_environ_density( prod )
    !
    CALL destroy_environ_gradient( field )
    CALL destroy_environ_density( aux )
    !
    DO i = 1, nsoft_spheres
       CALL destroy_environ_density( local(i) )
    ENDDO
    DEALLOCATE( local )
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE compute_ion_field
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE compute_ion_field_partial( nsoft_spheres, soft_spheres, &
       & ions, electrons, ion_field, partial_of_ion_field )
!--------------------------------------------------------------------
    !
    USE utils_functions, ONLY : density_of_functions, gradient_of_functions, laplacian_of_functions, hessian_of_functions
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: nsoft_spheres
    TYPE( environ_functions), DIMENSION(nsoft_spheres), INTENT(IN) :: soft_spheres
    TYPE( environ_ions ), INTENT(IN) :: ions
    TYPE( environ_electrons ), INTENT(IN) :: electrons
    REAL( DP ), DIMENSION(nsoft_spheres), INTENT(OUT) :: ion_field
    REAL( DP ), DIMENSION(3,nsoft_spheres,nsoft_spheres), INTENT(OUT) :: partial_of_ion_field
    !
    INTEGER, POINTER :: nnr, ir_end, deriv
    TYPE( environ_cell ), POINTER :: cell
    !
    INTEGER :: i, j, k, ipol
    CHARACTER( LEN=80 ) :: sub_name = 'compute_ion_field'
    !
    TYPE( environ_density ), DIMENSION(:), ALLOCATABLE :: local
    TYPE( environ_gradient ), DIMENSION(:), ALLOCATABLE :: gradlocal
    TYPE( environ_hessian ) :: hesslocal
    !
    TYPE( environ_density ) :: aux, prod
    TYPE( environ_gradient ) :: gradaux, field
    TYPE( environ_hessian ) :: hessaux
    !
    ! Aliases and allocations
    !
    cell => ions % density % cell
    nnr => cell % nnr
    ir_end => cell % ir_end
    !
    ALLOCATE( local( nsoft_spheres ) )
    ALLOCATE( gradlocal( nsoft_spheres ) )
    !
    ! Compute field-independent soft spheres and gradients
    !
    DO i = 1, nsoft_spheres
       !
       CALL init_environ_density( cell, local(i) )
       !
       CALL init_environ_gradient( cell, gradlocal(i) )
       !
       CALL density_of_functions( soft_spheres(i), local(i), .FALSE. )
       !
       CALL gradient_of_functions( soft_spheres(i), gradlocal(i), .FALSE. )
       !
    ENDDO  
    !
    CALL init_environ_hessian( cell, hesslocal )
    !
    ! Compute field
    !
    CALL init_environ_density( cell, aux )
    aux % of_r = electrons % density % of_r + ions % density % of_r
    !
    CALL init_environ_gradient( cell, field )
    CALL gradv_h_of_rho_r( aux%of_r, field%of_r )
    !
    ! Compute field flux
    !
    ion_field = 0.d0
    partial_of_ion_field = 0.D0
    !
    CALL init_environ_density( cell, prod )
    CALL init_environ_gradient( cell, gradaux )
    CALL init_environ_hessian( cell, hessaux )
    !
    DO i = 1, nsoft_spheres
       !
       ! Compute product of other soft-spheres
       !
       prod % of_r = 1.D0
       DO j = 1, nsoft_spheres
          IF ( j .EQ. i ) CYCLE
          prod % of_r = prod % of_r * local(j) % of_r
       ENDDO
       !
       ! Compute field flux through soft-sphere interface
       !
       CALL scalar_product_environ_gradient( field, gradlocal(i), aux ) ! here aux is the normal field
       !
       aux % of_r = - aux % of_r * prod % of_r
       !
       ion_field(i) = integrate_environ_density( aux )
       !
       ! Compute partial derivatives of field flux wrt ionic positions
       !
       DO j = 1, nsoft_spheres
          !
          ! Compute hessian of poisson potential of individual nuclei
          !
          CALL density_of_functions( ions%smeared_ions(j), aux, .TRUE. ) ! THIS STEP SHOULD BE MOVED OUT OF THIS LOOP
          !
          CALL hessv_h_of_rho_r( aux%of_r, hesslocal%of_r ) ! THIS STEP SHOULD BE MOVED OUT OF THIS LOOP
          !
          CALL scalar_product_environ_hessian( hesslocal, gradlocal(i), gradaux )
          !
          partial_of_ion_field( :, i, j ) = partial_of_ion_field( :, i, j ) - &
               & scalar_product_environ_gradient_density( gradaux, prod )
          !
          IF ( j .EQ. i ) THEN
             !
             ! Hessian of soft-sphere times the field
             !
             CALL hessian_of_functions( soft_spheres(i), hessaux, .TRUE. )
             !
             CALL scalar_product_environ_hessian( hessaux, field, gradaux )
             !
             partial_of_ion_field( :, i, j ) = partial_of_ion_field( :, i, j ) + &
                  & scalar_product_environ_gradient_density( gradaux, prod )
             !
          ELSE
             !
             ! ion field times gradient of different soft-sphere
             !
             CALL scalar_product_environ_gradient( gradlocal(i), field, aux ) ! here aux is the normal field
             !
             DO k = 1, nsoft_spheres
                IF ( k .EQ. j .OR. k .EQ. i ) CYCLE
                aux % of_r = aux % of_r * local(k) % of_r
             ENDDO
             !
             partial_of_ion_field( :, i, j ) = partial_of_ion_field( :, i, j ) + &
                  & scalar_product_environ_gradient_density( gradlocal(j), aux )
             !
          ENDIF
          !
       ENDDO
       !
    END DO
    !
    CALL destroy_environ_gradient( field )
    CALL destroy_environ_density( prod )
    CALL destroy_environ_density( aux )
    CALL destroy_environ_gradient( gradaux )
    CALL destroy_environ_hessian( hessaux )
    !
    CALL destroy_environ_hessian( hesslocal )
    DO i = 1, nsoft_spheres
       CALL destroy_environ_gradient( gradlocal(i) )
       CALL destroy_environ_density( local(i) )
    ENDDO
    DEALLOCATE( gradlocal )
    DEALLOCATE( local )
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE compute_ion_field_partial
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE compute_dion_field_drho( nsoft_spheres, soft_spheres, dion_field_drho )
!--------------------------------------------------------------------
    !
    USE utils_functions, ONLY : density_of_functions, gradient_of_functions
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: nsoft_spheres
    TYPE( environ_functions), DIMENSION(nsoft_spheres), INTENT(IN) :: soft_spheres
    TYPE( environ_density ), DIMENSION(nsoft_spheres), INTENT(INOUT) :: dion_field_drho
    !
    TYPE( environ_cell ), POINTER :: cell
    !
    INTEGER :: i, j, ipol
    CHARACTER( LEN=80 ) :: sub_name = 'compute_dion_field_drho'
    !
    TYPE( environ_density ), DIMENSION(:), ALLOCATABLE :: local
    !
    TYPE( environ_density ) :: prod
    TYPE( environ_gradient ) :: gradaux
    !
    ! Aliases and allocations
    !
    IF ( nsoft_spheres .LT. 1 ) CALL errore(sub_name,'Missing soft-spheres',1)
    cell => dion_field_drho(1) % cell
    !
    ALLOCATE( local( nsoft_spheres ) )
    !
    ! Compute field-independent soft spheres and gradients
    !
    DO i = 1, nsoft_spheres
       !
       CALL init_environ_density( cell, local(i) )
       !
       CALL density_of_functions( soft_spheres(i), local(i), .FALSE. )
       !
    ENDDO  
    !
    ! Compute field flux
    !
    CALL init_environ_density( cell, prod )
    CALL init_environ_gradient( cell, gradaux )
    !
    DO i = 1, nsoft_spheres
       !
       ! Compute product of other soft-spheres
       !
       prod % of_r = 1.D0
       DO j = 1, nsoft_spheres
          IF ( j .EQ. i ) CYCLE
          prod % of_r = prod % of_r * local(j) % of_r
       ENDDO
       !
       ! Compute functional derivative of field flux wrt electronic density
       !
       CALL gradient_of_functions( soft_spheres(i), gradaux, .TRUE. )
       !
       DO ipol = 1, 3
          gradaux % of_r(ipol,:) = gradaux % of_r(ipol,:) * prod % of_r(:)
       ENDDO
       !
       CALL field_of_gradrho( gradaux%of_r, dion_field_drho(i)%of_r )
       !
    END DO
    !
    CALL destroy_environ_gradient( gradaux )
    !
    DO i = 1, nsoft_spheres
       CALL destroy_environ_density( local(i) )
    ENDDO
    DEALLOCATE( local )
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE compute_dion_field_drho
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE compute_normal_field( ions, electrons, normal_field )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    TYPE( environ_ions ), INTENT(IN) :: ions
    TYPE( environ_electrons ), INTENT(IN) :: electrons
    TYPE( environ_density ), INTENT(INOUT) :: normal_field
    !
    TYPE( environ_cell ), POINTER :: cell
    TYPE( environ_density ) :: rho
    TYPE( environ_gradient ) :: field, gradrho
    !
    cell => normal_field % cell
    !
    ! Compute field
    !
    CALL init_environ_density( cell, rho )
    rho % of_r = electrons % density % of_r + ions % density % of_r
    !
    CALL init_environ_gradient( cell, field )
    CALL gradv_h_of_rho_r( rho%of_r, field%of_r )
    !
    CALL init_environ_gradient( cell, gradrho )
    CALL external_gradient( electrons%density%of_r, gradrho%of_r )
    !
    CALL scalar_product_environ_gradient( field, gradrho, normal_field )
    !
    CALL destroy_environ_gradient( gradrho )
    CALL destroy_environ_gradient( field )
    CALL destroy_environ_density( rho )
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE compute_normal_field
!--------------------------------------------------------------------
!----------------------------------------------------------------------------
END MODULE tools_generate_boundary
!----------------------------------------------------------------------------
