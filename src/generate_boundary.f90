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
  USE environ_output
  USE electrostatic_base, ONLY : boundary_core, nfdpoint, icfd, ncfd
  !
  PRIVATE
  !
  PUBLIC :: boundary_of_density, boundary_of_functions, calc_dboundary_dions
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
      USE fd_gradient, ONLY : calc_fd_gradient
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
      INTEGER, POINTER :: ir_end, nnr, stype, deriv
      REAL( DP ), POINTER :: const, rhomax, rhomin, tbeta
      REAL( DP ), DIMENSION(:), POINTER :: rho, eps, deps, d2eps
      REAL( DP ), DIMENSION(:), POINTER :: lapleps, dsurface
      REAL( DP ), DIMENSION(:,:), POINTER :: gradeps
      !
      INTEGER :: ir, ipol
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
      IF ( deriv .GE. 3 ) dsurface => boundary%dsurface%of_r
      !
      SELECT CASE ( boundary_core )
         !
      CASE ( 'fft' )
         !
         IF ( deriv .EQ. 1 .OR. deriv .EQ. 2 ) CALL external_gradient( eps, gradeps )
         IF ( deriv .EQ. 2 ) CALL external_laplacian( eps, lapleps )
         IF ( deriv .EQ. 3 ) CALL external_dsurface( nnr, ir_end, eps, gradeps, lapleps, dsurface )
         !
      CASE ( 'analytic', 'fd' )
         !
         IF ( deriv .EQ. 1 .OR. deriv .EQ. 2 ) CALL external_gradient( rho , gradeps )
         IF ( deriv .EQ. 2 ) CALL external_laplacian( rho, lapleps )
         IF ( deriv .EQ. 3 ) CALL external_dsurface( nnr, ir_end, rho, gradeps, lapleps, dsurface )
         IF ( deriv .GT. 1 ) lapleps(:) = lapleps(:) * deps(:) + &
              & ( gradeps(1,:)**2 + gradeps(2,:)**2 + gradeps(3,:)**2 ) * d2eps(:)
         IF ( deriv .GE. 1 ) THEN
            IF ( boundary_core .EQ. 'analytic' ) THEN
               DO ipol = 1, 3
                  gradeps(ipol,:) = gradeps(ipol,:) * deps(:)
               ENDDO
            ELSE IF ( boundary_core .EQ. 'fd' ) THEN
               CALL calc_fd_gradient( nfdpoint, icfd, ncfd, nnr, eps, gradeps )
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
      RETURN
      !
!--------------------------------------------------------------------
    END SUBROUTINE boundary_of_density
!--------------------------------------------------------------------
!--------------------------------------------------------------------
    SUBROUTINE external_dsurface( n, iend, x, grad, lapl, dsurface )
!--------------------------------------------------------------------
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN) :: n, iend
      REAL( DP ), DIMENSION( n ), INTENT(IN) :: x
      REAL( DP ), DIMENSION( 3, n ), INTENT(OUT) :: grad
      REAL( DP ), DIMENSION( n ), INTENT(OUT) :: lapl
      REAL( DP ), DIMENSION( n ), INTENT(OUT) :: dsurface
      !
      REAL( DP ), DIMENSION( 3, 3, n ) :: hess
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
      REAL( DP ), PARAMETER :: toldsurface = 1.D-20
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
!--------------------------------------------------------------------
    SUBROUTINE boundary_of_functions( nsoft_spheres, soft_spheres, boundary )
!--------------------------------------------------------------------
      !
      ! ... Calculates the dielectric constant as a function
      ! ... of the charge density, and the derivative of
      ! ... the dielectric constant wrt the charge density.
      !
      USE functions, ONLY: density_of_functions, gradient_of_functions, &
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
      TYPE( environ_density), DIMENSION(:), ALLOCATABLE :: lapllocal
      TYPE( environ_hessian), DIMENSION(:), ALLOCATABLE :: hesslocal
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
         CALL density_of_functions( 1, soft_spheres(i), local(i) )
         !
         boundary % scaled % of_r = boundary % scaled % of_r * local(i) % of_r
         !
      ENDDO
      !
      boundary % scaled % of_r = 1.D0 - boundary % scaled % of_r
      !
      boundary % volume = integrate_environ_density( boundary % scaled )
      !
      ! Generate boundary derivatives, if needed
      !
      deriv => boundary % deriv
      !
      SELECT CASE ( boundary_core )
         !
      CASE ( 'fft' )
         !
         IF ( deriv .EQ. 1 .OR. deriv .EQ. 2 ) CALL external_gradient( boundary%scaled%of_r, &
              & boundary%gradient%of_r )
         IF ( deriv .EQ. 2 ) CALL external_laplacian( boundary%scaled%of_r, boundary%laplacian%of_r )
         IF ( deriv .EQ. 3 ) CALL external_dsurface( nnr, ir_end, boundary%scaled%of_r, &
              & boundary%gradient%of_r, boundary%laplacian%of_r, boundary%dsurface%of_r )
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
            IF ( deriv .GE. 1 ) CALL gradient_of_functions( 1, soft_spheres(i), gradlocal(i) )
            IF ( deriv .EQ. 2 ) CALL laplacian_of_functions( 1, soft_spheres(i), lapllocal(i) )
            IF ( deriv .EQ. 3 ) CALL hessian_of_functions( 1, soft_spheres(i), hesslocal(i) )
            !
         END DO
         !
         IF ( deriv .EQ. 1 .OR. deriv .EQ. 2 ) &
              & CALL calc_gradient_of_boundary_highmem(nsoft_spheres,local, &
              & gradlocal,boundary%gradient )
         IF ( deriv .EQ. 2 ) CALL calc_laplacian_of_boundary_highmem(nsoft_spheres, &
              & local,gradlocal,lapllocal,boundary%laplacian )
         IF ( deriv .EQ. 3 ) CALL calc_dsurface_of_boundary_highmem(nsoft_spheres, &
              & local,gradlocal,hesslocal,boundary%gradient,boundary%laplacian, &
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
      IF ( deriv .GE. 1 ) THEN
         CALL update_gradient_modulus( boundary%gradient )
         boundary % surface = integrate_environ_density( boundary%gradient%modulus )
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
    SUBROUTINE calc_dsurface_of_boundary_highmem(n,local,gradlocal,hesslocal,gradient,laplacian,dsurface)
!--------------------------------------------------------------------
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN) :: n
      TYPE( environ_density ), DIMENSION(n), INTENT(IN) :: local
      TYPE( environ_gradient ), DIMENSION(n), INTENT(IN) :: gradlocal
      TYPE( environ_hessian ), DIMENSION(n), INTENT(IN) :: hesslocal
      TYPE( environ_gradient ), INTENT(INOUT) :: gradient
      TYPE( environ_density ), INTENT(INOUT) :: laplacian, dsurface
      !
      INTEGER :: i, j, k, ipol, jpol
      TYPE( environ_cell), POINTER :: cell
      TYPE( environ_density ) :: dens
      TYPE( environ_gradient ) :: partial
      TYPE( environ_hessian ) :: hess
      !
      cell => laplacian % cell
      CALL init_environ_density( cell, dens )
      CALL init_environ_gradient( cell, partial )
      CALL init_environ_hessian( cell, hess )
      !
      gradient % of_r = 0.D0
      DO i = 1, n
         CALL calc_partial_of_boundary_highmem(n,i,local,gradlocal,partial)
         gradient % of_r  = gradient % of_r + partial % of_r
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
                  hess % of_r( ipol, jpol, : ) = hess % of_r( ipol, jpol, : ) + dens % of_r( : )
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      !
      ! Final operations
      !
      laplacian % of_r = hess % of_r( 1, 1, : ) + hess % of_r( 2, 2, : ) + hess % of_r( 3, 3, : )
      CALL calc_dsurface( cell%nnr, cell%ir_end, gradient%of_r, hess%of_r, dsurface%of_r )
      !
      CALL destroy_environ_density( dens )
      CALL destroy_environ_gradient( partial )
      CALL destroy_environ_hessian( hess )
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
      USE functions, ONLY: density_of_functions, gradient_of_functions
      !
      IMPLICIT NONE
      !
      REAL( DP ), PARAMETER :: tolspuriousforce = 1.D-5
      !
      INTEGER, INTENT(IN) :: index
      TYPE( environ_boundary ), INTENT(IN) :: boundary
      TYPE( environ_gradient ), INTENT(INOUT) :: partial
      !
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
      !
      IF ( index .GT. boundary % ions % number ) &
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
         CASE( 'fft', 'analytic', 'highmem' )
            !
            CALL gradient_of_functions( 1, boundary%soft_spheres(index), partial )
            !
            CALL init_environ_density( cell, local )
            !
            DO i = 1, boundary % ions % number
               IF ( i .EQ. index ) CYCLE
               CALL density_of_functions( 1, boundary%soft_spheres(i), local )
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
         CALL gradient_of_functions( 1, boundary%ions%core_electrons(index), partial )
         DO ipol = 1, 3
            partial % of_r( ipol, : ) = - partial % of_r( ipol, : ) * boundary % dscaled % of_r( : )
         END DO
         CALL update_gradient_modulus( partial )
         spurious_force = integrate_environ_density( partial%modulus )
         IF ( spurious_force .GT. tolspuriousforce ) &
              & CALL errore(sub_name,'Unphysical forces due to core electrons are too high',1)
         !
      END IF
      !
      RETURN
      !
!--------------------------------------------------------------------
    END SUBROUTINE calc_dboundary_dions
!--------------------------------------------------------------------

END MODULE generate_boundary
