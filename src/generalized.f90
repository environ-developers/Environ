!
! Copyright (C) Oliviero Andreussi
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Module for implicit solvation according to the Self-Consistent
! Continuum Solvation (SCCS) model of Andreussi et al.
!  J. Chem. Phys. 136, 064102 (2012).
!
! original version by O. Andreussi and N. Marzari
! includes improved algorithms from G. Fisicaro and S. Goedecker
!
!--------------------------------------------------------------------
MODULE generalized
!--------------------------------------------------------------------

  USE environ_types
  USE environ_output
  USE poisson, ONLY : poisson_direct, poisson_gradient_direct
  USE electrostatic_base, ONLY : auxiliary, preconditioner, solver, &
       maxiter, tolvelect
  USE environ_base, ONLY : e2

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: generalized_gradient

CONTAINS

!--------------------------------------------------------------------
SUBROUTINE generalized_gradient( charges, dielectric, potential )
!--------------------------------------------------------------------

  IMPLICIT NONE

  TYPE( environ_charges ), INTENT(INOUT) :: charges
  TYPE( environ_dielectric ), INTENT(IN) :: dielectric
  TYPE( environ_density ), INTENT(INOUT) :: potential

  CHARACTER*20 :: sub_name = 'generalized_gradient'

  CALL start_clock( 'calc_vsolv' )

  SELECT CASE ( auxiliary )

  CASE ( 'none' )

     SELECT CASE ( preconditioner )

     CASE ( 'none' )

        CALL generalized_gradient_none( charges, dielectric, potential )

     CASE ( 'sqrt' )

        CALL generalized_gradient_sqrt( charges, dielectric, potential )

     CASE ( 'left' )

        CALL errore(sub_name,'Option not yet implemented',1)
!        CALL generalized_gradient_left( charges, dielectric, potential )

     CASE DEFAULT

        CALL errore( sub_name, 'unexpected preconditioner keyword', 1 )

     END SELECT

  CASE ( 'full' )

     IF ( solver .EQ. 'iterative' ) THEN

        CALL generalized_iterative( charges, dielectric, potential )

     ELSE

        CALL errore(sub_name,'Option not yet implemented',1)
!        CALL generalized_gradient_rhoaux( charges, dielectric, potential )

     ENDIF

  CASE DEFAULT

     CALL errore( sub_name, 'unexpected auxiliary keyword', 1 )

  END SELECT

  CALL stop_clock( 'calc_vsolv' )

  RETURN

!--------------------------------------------------------------------
END SUBROUTINE generalized_gradient
!--------------------------------------------------------------------
!--------------------------------------------------------------------
SUBROUTINE generalized_iterative( charges, dielectric, potential )
!--------------------------------------------------------------------

  USE electrostatic_base, ONLY : mix, tolrhoaux

  IMPLICIT NONE

  TYPE( environ_charges ), TARGET, INTENT(INOUT) :: charges
  TYPE( environ_dielectric ), TARGET, INTENT(IN) :: dielectric
  TYPE( environ_density ), TARGET, INTENT(INOUT) :: potential

  TYPE( environ_cell ), POINTER :: cell
  TYPE( environ_density ), POINTER :: eps, rhoiter, rhozero
  TYPE( environ_gradient ), POINTER :: gradlogeps

  INTEGER :: iter
  REAL( DP ) :: total, totpol, totzero, totiter, deltar
  TYPE( environ_density ) :: residual
  TYPE( environ_gradient ) :: gradpoisson

  CHARACTER( LEN=80 ) :: sub_name = 'generalized_iterative'

  IF ( verbose .GE. 1 ) WRITE(environ_unit,9000)
9000 FORMAT(/,4('%'),' COMPUTE ELECTROSTATIC POTENTIAL ',43('%'))

  ! ... Check that fields have the same defintion domain

  IF ( .NOT. ASSOCIATED(charges%density%cell,dielectric%epsilon%cell) ) &
       & CALL errore(sub_name,'Inconsistent cells of input fields',1)
  IF ( .NOT. ASSOCIATED(charges%density%cell,potential%cell) ) &
       & CALL errore(sub_name,'Inconsistent cells for charges and potential',1)
  cell => charges%density%cell

  ! ... Aliases

  eps => dielectric % epsilon
  gradlogeps => dielectric % gradlog
  rhoiter => charges % auxiliary % iterative
  rhozero => charges % auxiliary % fixed

  ! ... Solute only

  charges % include_auxiliary = .FALSE.
  CALL update_environ_charges( charges )
  total = charges % charge

  ! ... Set up auxiliary charge

  totpol = total * ( 1.D0 - dielectric % constant ) / dielectric % constant
  rhozero % of_r = charges % density % of_r * ( 1.D0 - eps % of_r ) / eps % of_r
  totzero = integrate_environ_density( rhozero )
  totiter = integrate_environ_density( rhoiter )
  IF ( verbose .GE. 1 ) WRITE(environ_unit,9001) totiter
9001 FORMAT(' Starting from polarization: rhoiter = ',F13.6)
  charges % include_auxiliary = .TRUE.

  ! ... Create local variables

  CALL init_environ_density( cell, residual )
  CALL init_environ_gradient( cell, gradpoisson )

  ! ... Start iterative algorithm

  DO iter = 1, maxiter

       IF ( verbose .GE. 1 ) WRITE(environ_unit,9002) iter
9002 FORMAT(' Iteration # ',i10)

       CALL update_environ_auxiliary( charges % auxiliary )

       CALL update_environ_charges( charges )

       CALL poisson_gradient_direct( charges, gradpoisson )

       CALL scalar_product_environ_gradient( gradlogeps, gradpoisson, residual )

       residual % of_r = residual % of_r / fpi / e2 - rhoiter % of_r

       rhoiter % of_r = rhoiter % of_r + mix * residual % of_r

       ! ... If residual is small enough exit

       deltar = quadratic_mean_environ_density( residual )
       totiter = integrate_environ_density( rhoiter )
       IF ( verbose .GE. 1 ) WRITE(environ_unit,9004)deltar,tolrhoaux
9004   FORMAT(' deltarho = ',E14.6,' tol = ',E14.6)
       IF ( verbose .GE. 2 ) WRITE(environ_unit,9003)totiter,totzero,totpol,total
9003   FORMAT(' Total iterative polarization charge = ',4F13.6)
       IF ( deltar .LT. tolrhoaux .AND. iter .GT. 0 ) THEN
          IF ( verbose .GE. 1 ) WRITE(environ_unit,9005)
9005      FORMAT(' Charges are converged, exit!')
          EXIT
       ELSE IF ( iter .EQ. maxiter ) THEN
         WRITE(program_unit,9006)
9006     FORMAT(' Warning: Polarization charge not converged')
       ENDIF

    ENDDO

    IF (.not.tddfpt.AND.verbose.GE.1) WRITE(program_unit, 9007) deltar, iter
9007 FORMAT('     polarization accuracy =',1PE8.1,', # of iterations = ',i3)

    ! ... Compute total electrostatic potential

    CALL update_environ_auxiliary( charges % auxiliary )

    CALL update_environ_charges( charges )

    CALL poisson_direct( charges, potential )

    ! ... Destroy local variables

    CALL destroy_environ_density( residual )
    CALL destroy_environ_gradient( gradpoisson )

  RETURN

!--------------------------------------------------------------------
END SUBROUTINE generalized_iterative
!--------------------------------------------------------------------
!--------------------------------------------------------------------
SUBROUTINE generalized_gradient_none( charges, dielectric, potential )
!--------------------------------------------------------------------

  IMPLICIT NONE

  TYPE( environ_charges ), TARGET, INTENT(IN) :: charges
  TYPE( environ_dielectric ), TARGET, INTENT(IN) :: dielectric
  TYPE( environ_density ), TARGET, INTENT(INOUT) :: potential

  TYPE( environ_cell ), POINTER :: cell
  TYPE( environ_density ), POINTER :: x, b, eps
  TYPE( environ_gradient ), POINTER :: gradeps

  INTEGER :: iter
  REAL( DP ) :: rznew, rzold, alpha, beta, pAp, deltar
  TYPE( environ_density ) :: r, z, p, Ap, l
  TYPE( environ_gradient ) :: g

  CHARACTER( LEN=80 ) :: sub_name = 'generalized_gradient_none'

  IF ( verbose .GE. 1 ) WRITE(environ_unit,9000)
9000 FORMAT(/,4('%'),' COMPUTE ELECTROSTATIC POTENTIAL ',43('%'))

  ! ... Check that fields have the same defintion domain

  IF ( .NOT. ASSOCIATED(charges%density%cell,dielectric%epsilon%cell) ) &
       & CALL errore(sub_name,'Inconsistent cells of input fields',1)
  IF ( .NOT. ASSOCIATED(charges%density%cell,potential%cell) ) &
       & CALL errore(sub_name,'Inconsistent cells for charges and potential',1)
  cell => charges%density%cell

  ! ... Aliases

  b => charges % density
  eps => dielectric % epsilon
  gradeps => dielectric % gradient
  x => potential

  ! ... Create and initialize local variables

  CALL init_environ_density( cell, r )
  CALL init_environ_density( cell, z )
  CALL init_environ_density( cell, p )
  CALL init_environ_density( cell, Ap )

  CALL init_environ_gradient( cell, g )
  CALL init_environ_density( cell, l )

  ! ... Starting guess from new input and previous solution(s)

  x%of_r = 0.D0
  r = b
  z = r ! no preconditioner
  p = z
  rzold = scalar_product_environ_density( r, z )
  IF ( ABS(rzold) .LT. 1.D-30 ) &
       & CALL errore(sub_name,'Null step in gradient descent iteration',1)

  ! ... Start gradient descent

  DO iter = 1, maxiter

       IF ( verbose .GE. 1 ) WRITE(environ_unit,9002) iter
9002 FORMAT(' Iteration # ',i10)

       ! ... Apply operator to conjugate direction

       CALL external_gradient(p%of_r,g%of_r)
       CALL external_laplacian(p%of_r,l%of_r)
       Ap%of_r(:) = eps%of_r(:)*l%of_r(:) + &
                  & gradeps%of_r(1,:)*g%of_r(1,:) + &
                  & gradeps%of_r(2,:)*g%of_r(2,:) + &
                  & gradeps%of_r(3,:)*g%of_r(3,:)

       ! ... Step downhill

       pAp = scalar_product_environ_density( p, Ap )
       alpha = rzold / pAp

       x%of_r = x%of_r + alpha * p%of_r
       r%of_r = r%of_r - alpha * Ap%of_r

       ! ... If residual is small enough exit

       deltar = quadratic_mean_environ_density( r )
       IF ( verbose .GE. 1 ) WRITE(environ_unit,9004)deltar,tolvelect
9004   FORMAT(' deltarho = ',E14.6,' tol = ',E14.6)
       IF ( deltar .LT. tolvelect .AND. iter .GT. 0 ) THEN
          IF ( verbose .GE. 1 ) WRITE(environ_unit,9005)
9005      FORMAT(' Charges are converged, exit!')
          EXIT
       ELSE IF ( iter .EQ. maxiter ) THEN
         WRITE(program_unit,9006)
9006     FORMAT(' Warning: Polarization charge not converged')
       ENDIF

       ! ... Apply preconditioner to new state

       z = r ! no preconditioner

       rznew = scalar_product_environ_density( r, z )
       IF ( ABS(rznew) .LT. 1.D-30 ) &
            & CALL errore(sub_name,'Null step in gradient descent iteration',1)

       ! ... Conjugate gradient or steepest descent input

       IF ( solver .EQ. 'cg' ) THEN
          p%of_r = z%of_r + rznew / rzold * p%of_r
       ELSE
          p = z
       END IF
       rzold = rznew

    ENDDO

    IF (.not.tddfpt.AND.verbose.GE.1) WRITE(program_unit, 9007) deltar, iter
9007 FORMAT('     polarization accuracy =',1PE8.1,', # of iterations = ',i3)

    CALL destroy_environ_density( l )
    CALL destroy_environ_gradient( g )

    CALL destroy_environ_density( r )
    CALL destroy_environ_density( z )
    CALL destroy_environ_density( p )
    CALL destroy_environ_density( Ap )

  RETURN

!--------------------------------------------------------------------
END SUBROUTINE generalized_gradient_none
!--------------------------------------------------------------------
!--------------------------------------------------------------------
SUBROUTINE generalized_gradient_sqrt( charges, dielectric, potential )
!--------------------------------------------------------------------

  IMPLICIT NONE

  TYPE( environ_charges ), TARGET, INTENT(IN) :: charges
  TYPE( environ_dielectric ), TARGET, INTENT(IN) :: dielectric
  TYPE( environ_density ), TARGET, INTENT(INOUT) :: potential

  TYPE( environ_cell ), POINTER :: cell
  TYPE( environ_density ), POINTER :: x, b, eps, factsqrt
  TYPE( environ_gradient ), POINTER :: gradeps

  INTEGER :: iter
  REAL( DP ) :: rznew, rzold, alpha, beta, pAp, deltar
  TYPE( environ_density ) :: r, z, p, Ap, invsqrt

  CHARACTER( LEN=80 ) :: sub_name = 'generalized_gradient_sqrt'

  IF ( verbose .GE. 1 ) WRITE(environ_unit,9000)
9000 FORMAT(/,4('%'),' COMPUTE ELECTROSTATIC POTENTIAL ',43('%'))

  ! ... Check that fields have the same defintion domain

  IF ( .NOT. ASSOCIATED(charges%density%cell,dielectric%epsilon%cell) ) &
       & CALL errore(sub_name,'Inconsistent cells of input fields',1)
  IF ( .NOT. ASSOCIATED(charges%density%cell,potential%cell) ) &
       & CALL errore(sub_name,'Inconsistent cells for charges and potential',1)
  cell => charges%density%cell

  ! ... Aliases

  b => charges % density
  eps => dielectric % epsilon
  factsqrt => dielectric % factsqrt
  x => potential
  CALL init_environ_density( cell, invsqrt )
  invsqrt%of_r = 1.D0 / SQRT(eps%of_r)

  ! ... Create and initialize local variables

  CALL init_environ_density( cell, r )
  CALL init_environ_density( cell, z )
  CALL init_environ_density( cell, p )
  CALL init_environ_density( cell, Ap )

  ! ... Starting guess from new input and previous solution(s)

  IF ( x%update ) THEN

     r%of_r = b%of_r - factsqrt%of_r * x%of_r

     ! ... Preconditioning step

     z%of_r = r%of_r * invsqrt%of_r
     CALL poisson_direct( z, z )
     z%of_r = z%of_r * invsqrt%of_r

     rzold = scalar_product_environ_density( r, z )
     IF ( ABS(rzold) .LT. 1.D-30 ) &
          & CALL errore(sub_name,'Null step in gradient descent iteration',1)

     r%of_r = factsqrt%of_r * ( x%of_r - z%of_r )
     deltar = quadratic_mean_environ_density( r )
     IF ( deltar .LT. 1.D-02 ) THEN
        IF ( verbose .GE. 1 ) WRITE(environ_unit,9008)deltar
9008    FORMAT(' Sqrt-preconditioned input guess with residual norm = ',E14.6)
        x%of_r = z%of_r
     ELSE
        IF ( verbose .GE. 1 ) WRITE(environ_unit,9001)deltar
9001    FORMAT(' Warning: bad guess with residual norm = ',E14.6,', reset to no guess')
        x%update = .FALSE.
     ENDIF

  ENDIF

  IF ( .NOT. x%update ) THEN

     x%update = .TRUE.
     x%of_r = 0.D0
     r = b
     rzold = 0.D0

  ENDIF

  ! ... Start gradient descent

  DO iter = 1, maxiter

       IF ( verbose .GE. 1 ) WRITE(environ_unit,9002) iter
9002   FORMAT(' Iteration # ',i10)

       ! ... Apply preconditioner to new state

       z%of_r = r%of_r * invsqrt%of_r
       CALL poisson_direct( z, z )
       z%of_r = z%of_r * invsqrt%of_r

       rznew = scalar_product_environ_density( r, z )
       IF ( ABS(rznew) .LT. 1.D-30 ) &
            & CALL errore(sub_name,'Null step in gradient descent iteration',1)

       ! ... Conjugate gradient or steepest descent input

       IF ( solver .EQ. 'cg' .AND. ABS(rzold) .GT. 1.D-30 ) THEN
          beta = rznew / rzold
          p%of_r = z%of_r + beta * p%of_r
       ELSE
          beta = 0.D0
          p%of_r = z%of_r
       END IF
       rzold = rznew

       ! ... Apply operator to conjugate direction

       Ap%of_r = factsqrt%of_r * z%of_r + r%of_r + beta * Ap%of_r

       ! ... Step downhill

       pAp = scalar_product_environ_density( p, Ap )
       alpha = rzold / pAp

       x%of_r = x%of_r + alpha * p%of_r
       r%of_r = r%of_r - alpha * Ap%of_r

       IF ( verbose .GE. 1 ) WRITE(environ_unit,*)'alpha = ',alpha,' beta = ',beta
       IF ( verbose .GE. 2 ) WRITE(environ_unit,*)'rznew = ',rznew,' rzold = ',rzold,' pAp = ',pAp

       ! ... If residual is small enough exit

       deltar = quadratic_mean_environ_density( r )
       IF ( verbose .GE. 1 ) WRITE(environ_unit,9004)deltar,tolvelect
       9004 FORMAT(' deltarho = ',E14.6,' tol = ',E14.6)
       IF ( deltar .LT. tolvelect .AND. iter .GT. 0 ) THEN
          IF ( verbose .GE. 1 ) WRITE(environ_unit,9005)
          9005 FORMAT(' Charges are converged, exit!')
          EXIT
       ELSE IF ( iter .EQ. maxiter ) THEN
         WRITE(program_unit,9006)
9006     FORMAT(' Warning: Polarization charge not converged')
      ENDIF

    ENDDO

    IF (.not.tddfpt.AND.verbose.GE.1) WRITE(program_unit, 9007) deltar, iter
9007 FORMAT('     polarization accuracy =',1PE8.1,', # of iterations = ',i3)

    CALL destroy_environ_density( r )
    CALL destroy_environ_density( z )
    CALL destroy_environ_density( p )
    CALL destroy_environ_density( Ap )

    CALL destroy_environ_density(invsqrt)

  RETURN

!--------------------------------------------------------------------
END SUBROUTINE generalized_gradient_sqrt
!--------------------------------------------------------------------
!!!TEMPLATE!!!!--------------------------------------------------------------------
!!!TEMPLATE!!!SUBROUTINE generalized_gradient_template( charges, dielectric, potential )
!!!TEMPLATE!!!!--------------------------------------------------------------------
!!!TEMPLATE!!!
!!!TEMPLATE!!!  IMPLICIT NONE
!!!TEMPLATE!!!
!!!TEMPLATE!!!  TYPE( environ_charges ), INTENT(IN) :: charges
!!!TEMPLATE!!!  TYPE( environ_dielectric ), INTENT(IN) :: dielectric
!!!TEMPLATE!!!  TYPE( environ_density ), INTENT(INOUT) :: potential
!!!TEMPLATE!!!
!!!TEMPLATE!!!  TYPE( environ_cell ), POINTER :: cell
!!!TEMPLATE!!!  TYPE( environ_density ), POINTER :: x, b, eps
!!!TEMPLATE!!!  TYPE( environ_gradient ), POINTER :: gradeps
!!!TEMPLATE!!!
!!!TEMPLATE!!!  INTEGER :: iter
!!!TEMPLATE!!!  REAL( DP ) :: rznew, rzold, alpha, beta, pAp, deltar
!!!TEMPLATE!!!  TYPE( environ_density ) :: r, z, p, Ap
!!!TEMPLATE!!!
!!!TEMPLATE!!!  ! ... Check that fields have the same defintion domain
!!!TEMPLATE!!!
!!!TEMPLATE!!!  IF ( .NOT. ASSOCIATED(charges%density%cell,dielectric%epsilon%cell) ) &
!!!TEMPLATE!!!       & CALL errore(sub_name,'Inconsistent cells of input fields',1)
!!!TEMPLATE!!!  IF ( .NOT. ASSOCIATED(charges%density%cell,potential%cell) ) &
!!!TEMPLATE!!!       & CALL errore(sub_name,'Inconsistent cells for charges and potential',1)
!!!TEMPLATE!!!  cell => charges%density%cell
!!!TEMPLATE!!!
!!!TEMPLATE!!!  ! ... Aliases
!!!TEMPLATE!!!
!!!TEMPLATE!!!  b => charges % density
!!!TEMPLATE!!!  eps => dielectric % epsilon
!!!TEMPLATE!!!  gradeps => dielectric % gradient
!!!TEMPLATE!!!  x => potential
!!!TEMPLATE!!!
!!!TEMPLATE!!!  ! ... Create and initialize local variables
!!!TEMPLATE!!!
!!!TEMPLATE!!!  CALL init_environ_density( cell, r )
!!!TEMPLATE!!!  CALL init_environ_density( cell, z )
!!!TEMPLATE!!!  CALL init_environ_density( cell, p )
!!!TEMPLATE!!!  CALL init_environ_density( cell, Ap )
!!!TEMPLATE!!!
!!!TEMPLATE!!!  ! ... Starting guess from new input and previous solution(s)
!!!TEMPLATE!!!
!!!TEMPLATE!!!  x%of_r = 0.D0
!!!TEMPLATE!!!  r = b
!!!TEMPLATE!!!  z = r ! no preconditioner
!!!TEMPLATE!!!  p = z
!!!TEMPLATE!!!  rzold = scalar_product_environ_density( r, z )
!!!TEMPLATE!!!  IF ( rzold .LT. 1.D-30 ) THEN
!!!TEMPLATE!!!     WRITE(program_unit,*)'ERROR: null step in gradient descent iteration'
!!!TEMPLATE!!!     STOP
!!!TEMPLATE!!!  ENDIF
!!!TEMPLATE!!!
!!!TEMPLATE!!!  ! ... Start gradient descent
!!!TEMPLATE!!!
!!!TEMPLATE!!!  DO iter = 1, maxiter
!!!TEMPLATE!!!
!!!TEMPLATE!!!       IF ( verbose .GE. 1 ) WRITE(environ_unit,9002) iter
!!!TEMPLATE!!!
!!!TEMPLATE!!!       ! ... Apply operator to conjugate direction
!!!TEMPLATE!!!
!!!TEMPLATE!!!       CALL external_hessian(p,g,h)
!!!TEMPLATE!!!       Ap%of_r(:) = eps%of_r(:)*(h%of_r(1,1,:)**2+h%of_r(2,2,:)**2+h%of_r(3,3,:)**2) + &
!!!TEMPLATE!!!                  & gradeps%of_r(1,:)*g%of_r(1,:) + &
!!!TEMPLATE!!!                  & gradeps%of_r(2,:)*g%of_r(2,:) + &
!!!TEMPLATE!!!                  & gradeps%of_r(3,:)*g%of_r(3,:)
!!!TEMPLATE!!!
!!!TEMPLATE!!!       ! ... Step downhill
!!!TEMPLATE!!!
!!!TEMPLATE!!!       pAp = scalar_product_environ_density( p, Ap )
!!!TEMPLATE!!!       alpha = rzold / pAp
!!!TEMPLATE!!!
!!!TEMPLATE!!!       x%of_r = x%of_r + alpha * p%of_r
!!!TEMPLATE!!!       r%of_r = r%of_r - alpha * Ap%of_r
!!!TEMPLATE!!!
!!!TEMPLATE!!!       ! ... If residual is small enough exit
!!!TEMPLATE!!!
!!!TEMPLATE!!!       deltar = quadratic_mean_environ_density(r)
!!!TEMPLATE!!!       IF ( verbose .GE. 1 ) WRITE(environ_unit,9004)deltar,tolvelect
!!!TEMPLATE!!!       IF ( deltar .LT. tolvelect .AND. iter .GT. 0 ) THEN
!!!TEMPLATE!!!         IF ( verbose .GE. 1 ) WRITE(environ_unit,9005)
!!!TEMPLATE!!!         EXIT
!!!TEMPLATE!!!       ELSE IF ( iter .EQ. maxiter ) THEN
!!!TEMPLATE!!!         WRITE(program_unit,9006)
!!!TEMPLATE!!!       ENDIF
!!!TEMPLATE!!!
!!!TEMPLATE!!!       ! ... Apply preconditioner to new state
!!!TEMPLATE!!!
!!!TEMPLATE!!!       z = r ! no preconditioner
!!!TEMPLATE!!!
!!!TEMPLATE!!!       rznew = scalar_product_environ_density( r, z )
!!!TEMPLATE!!!       IF ( ABS(rznew) .LT. 1.D-30 ) THEN
!!!TEMPLATE!!!          WRITE(program_unit,*)'ERROR: null step in gradient descent iteration'
!!!TEMPLATE!!!          STOP
!!!TEMPLATE!!!       ENDIF
!!!TEMPLATE!!!
!!!TEMPLATE!!!       ! ... Conjugate gradient or steepest descent input
!!!TEMPLATE!!!
!!!TEMPLATE!!!       IF ( solver .EQ. 'cg' ) THEN
!!!TEMPLATE!!!          p%of_r = z%of_r + rznew / rzold * p%of_r
!!!TEMPLATE!!!       ELSE
!!!TEMPLATE!!!          p = z
!!!TEMPLATE!!!       END IF
!!!TEMPLATE!!!       rzold = rznew
!!!TEMPLATE!!!
!!!TEMPLATE!!!    ENDDO
!!!TEMPLATE!!!
!!!TEMPLATE!!!    IF (.not.tddfpt.AND.verbose.GE.1) WRITE(program_unit, 9000) deltar, iter
!!!TEMPLATE!!!
!!!TEMPLATE!!!  RETURN
!!!TEMPLATE!!!
!!!TEMPLATE!!!!--------------------------------------------------------------------
!!!TEMPLATE!!!END SUBROUTINE generalized_gradient_template
!!!TEMPLATE!!!!--------------------------------------------------------------------
!--------------------------------------------------------------------
END MODULE generalized
!--------------------------------------------------------------------
