! Copyright (C) 2018 ENVIRON (www.quantum-environment.org)
!
!    This file is part of Environ version 1.1
!
!    Environ 1.1 is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 2 of the License, or
!    (at your option) any later version.
!
!    Environ 1.1 is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more detail, either the file
!    `License' in the root directory of the present distribution, or
!    online at <http://www.gnu.org/licenses/>.
!
!> This module contains the main drivers and routines to compute the
!! electrostatic potential that is the solution of a
!! generalized Poisson equation:
!! \f[
!!      \nabla \cdot \epsilon (r) \nabla \phi = -4 \pi \rho
!! \f]
!!
!! Different algorithms (gradient descent on potential and iterative
!! on the polarization charge) are available and implemented.
!
! Authors: Oliviero Andreussi (Department of Physics, UNT)
!          Nicola Marzari     (THEOS and NCCR-MARVEL, EPFL)
!
! includes improved algorithms from G. Fisicaro and S. Goedecker
!
!----------------------------------------------------------------------------
MODULE problem_generalized
!----------------------------------------------------------------------------
  !
  USE modules_constants, ONLY : e2, fpi
  USE environ_types
  USE electrostatic_types
  USE environ_output
  USE problem_poisson, ONLY : poisson_direct, poisson_gradient_direct!, poisson_energy
  USE environ_base, ONLY : oldenviron, add_jellium, ltddfpt
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC :: generalized_gradient!, generalized_energy
  !
  INTERFACE generalized_gradient
     MODULE PROCEDURE generalized_gradient_charges, generalized_gradient_density
  END INTERFACE generalized_gradient
  !
CONTAINS
!--------------------------------------------------------------------
  SUBROUTINE generalized_gradient_charges( solver, core, charges, potential )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    TYPE( electrostatic_solver ), INTENT(IN) :: solver
    TYPE( electrostatic_core ), INTENT(IN) :: core
    TYPE( environ_charges ), INTENT(INOUT) :: charges
    TYPE( environ_density ), INTENT(INOUT) :: potential
    !
    CHARACTER*20 :: sub_name = 'generalized_gradient'
    !
    CALL start_clock( 'calc_vsolv' )
    !
    potential % of_r = 0.D0
    !
    IF ( solver % use_gradient ) THEN
       !
       IF ( solver % auxiliary .EQ. 'none' ) THEN
          !
          SELECT CASE ( solver % gradient % preconditioner )
             !
          CASE ( 'none' )
             !
             CALL generalized_gradient_none( solver % gradient, core, charges%density, &
                        & charges%dielectric, potential, charges%electrolyte, &
                        & charges%semiconductor )
             !
          CASE ( 'sqrt' )
             !
             CALL generalized_gradient_sqrt( solver % gradient, core, charges%density, &
                       & charges%dielectric, potential, charges%electrolyte, &
                       & charges%semiconductor )
             !
          CASE ( 'left' )
             !
             CALL generalized_gradient_left( solver % gradient, core, charges%density, &
                       & charges%dielectric, potential, charges%electrolyte, charges%semiconductor )
             !
          CASE DEFAULT
             !
             CALL errore( sub_name, 'unexpected preconditioner keyword', 1 )
             !
          END SELECT
          !
       ELSE
          !
          CALL errore( sub_name, 'Option not yet implemented', 1 )
          !        CALL generalized_gradient_rhoaux( charges, dielectric, potential )
          !
       END IF
       !
    ELSE IF ( solver % use_iterative ) THEN
       !
       IF ( solver % auxiliary .EQ. 'full' ) THEN
          !

          CALL generalized_iterative( solver % iterative, core, charges%density, charges%dielectric, &
                           potential, charges%electrolyte, charges%semiconductor )

          !
       ELSE
          !
          CALL errore(sub_name,'Option not yet implemented',1)
          !        CALL generalized_iterative_velect( charges, dielectric, potential )
          !
       ENDIF
       !
    ELSE
       !
       CALL errore( sub_name, 'unexpected auxiliary keyword', 1 )
       !
    END IF
    !
    CALL stop_clock( 'calc_vsolv' )
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE generalized_gradient_charges
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE generalized_gradient_density( solver, core, charges, dielectric, potential,&
             & electrolyte, semiconductor )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    TYPE( electrostatic_solver ), INTENT(IN) :: solver
    TYPE( electrostatic_core ), INTENT(IN) :: core
    TYPE( environ_density ), INTENT(IN) :: charges
    TYPE( environ_dielectric ), INTENT(IN) :: dielectric
    TYPE( environ_density ), INTENT(INOUT) :: potential
    TYPE( environ_electrolyte ), INTENT(IN), OPTIONAL :: electrolyte
    TYPE( environ_semiconductor ), INTENT(INOUT), OPTIONAL :: semiconductor
    !
    CHARACTER*20 :: sub_name = 'generalized_gradient'
    !
    CALL start_clock( 'calc_vsolv' )
    !
    potential % of_r = 0.D0
    !
    IF ( solver % use_gradient ) THEN
       !
       IF ( solver % auxiliary .EQ. 'none' ) THEN
          !
          SELECT CASE ( solver % gradient % preconditioner )
             !
          CASE ( 'none' )
             !
             CALL generalized_gradient_none( solver % gradient, core, charges, dielectric, potential, electrolyte )
             !
          CASE ( 'sqrt' )
             !
             CALL generalized_gradient_sqrt( solver % gradient, core, charges, dielectric, potential, electrolyte )
             !
          CASE ( 'left' )
             !
             CALL generalized_gradient_left( solver % gradient, core, charges, dielectric, potential, electrolyte )
             !
          CASE DEFAULT
             !
             CALL errore( sub_name, 'unexpected preconditioner keyword', 1 )
             !
          END SELECT
          !
       ELSE
          !
          CALL errore( sub_name, 'Option not yet implemented', 1 )
          !        CALL generalized_gradient_rhoaux( charges, dielectric, potential )
          !
       END IF
       !
    ELSE IF ( solver % use_iterative ) THEN
       !
       IF ( solver % auxiliary .EQ. 'full' ) THEN
          !
          CALL generalized_iterative( solver % iterative, core, charges, dielectric, potential, electrolyte, semiconductor )
          !
       ELSE
          !
          CALL errore( sub_name, 'Option not yet implemented', 1 )
          !        CALL generalized_iterative_velect( charges, dielectric, potential )
          !
       ENDIF
       !
    ELSE
       !
       CALL errore( sub_name, 'unexpected auxiliary keyword', 1 )
       !
    END IF
    !
    CALL stop_clock( 'calc_vsolv' )
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE generalized_gradient_density
!--------------------------------------------------------------------
!!--------------------------------------------------------------------
!  SUBROUTINE generalized_energy( core, charges, potential, energy )
!!--------------------------------------------------------------------
!    !
!    IMPLICIT NONE
!    !
!    TYPE( electrostatic_core ), INTENT(IN) :: core
!    TYPE( environ_charges ), INTENT(IN) :: charges
!    TYPE( environ_density ), INTENT(IN) :: potential
!    REAL( DP ), INTENT(OUT) :: energy
!    !
!    REAL( DP ) :: degauss, eself
!    CHARACTER*20 :: sub_name = 'generalized_energy'
!    !
!    CALL start_clock( 'calc_esolv' )
!    !
!    ! Aliases and sanity checks
!    !
!    IF ( .NOT. ASSOCIATED( charges % density % cell, potential % cell ) ) &
!         & CALL errore(sub_name,'Missmatch in charges and potential domains',1)
!    !
!    energy = 0.D0
!    !
!    CALL poisson_energy( core, charges, potential, energy )
!    !
!    ! Adding correction for point-like nuclei: only affects simulations of charged
!    ! systems, it does not affect forces, but shift the energy depending on the
!    ! fictitious Gaussian spread of the nuclei
!    !
!    eself = 0.D0
!    degauss = 0.D0
!    !
!    IF ( charges % include_ions .AND. charges % ions % use_smeared_ions ) THEN
!       !
!       ! Compute spurious self-polarization energy
!       !
!       eself = 0.D0
!       !
!       IF ( core % use_fft ) THEN
!          !
!          IF ( core % fft % use_internal_pbc_corr .OR. core % need_correction ) THEN
!             !
!             degauss = 0.D0
!             !
!          ELSE
!             !
!             degauss = - charges % ions % quadrupole_correction * charges % dielectric % charge * e2 * pi &
!                  & / charges % density % cell % omega
!             !
!          ENDIF
!          !
!       ENDIF
!       !
!    ENDIF
!    !
!    energy = energy + eself + degauss
!    !
!    CALL stop_clock( 'calc_esolv' )
!    !
!    RETURN
!    !
!!--------------------------------------------------------------------
!  END SUBROUTINE generalized_energy
!!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE generalized_iterative( iterative, core, charges, dielectric, potential, electrolyte,&
                   & semiconductor )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    TYPE( iterative_solver ), TARGET, INTENT(IN) :: iterative
    TYPE( electrostatic_core ), INTENT(IN) :: core
    TYPE( environ_density ), TARGET, INTENT(IN) :: charges
    TYPE( environ_dielectric ), TARGET, INTENT(IN) :: dielectric
    TYPE( environ_density ), TARGET, INTENT(INOUT) :: potential
    TYPE( environ_electrolyte ), INTENT(IN), OPTIONAL :: electrolyte
    TYPE( environ_semiconductor ), INTENT(INOUT), OPTIONAL :: semiconductor
    !
    TYPE( environ_cell ), POINTER :: cell
    TYPE( environ_density ), POINTER :: eps, rhoiter, rhotot
    TYPE( environ_gradient ), POINTER :: gradlogeps
    !
    INTEGER :: iter
    REAL( DP ) :: total, totpol, totzero, totiter, delta_qm, delta_en, jellium
    TYPE( environ_density ) :: rhozero
    TYPE( environ_density ) :: residual
    TYPE( environ_gradient ) :: gradpoisson
    !
    CHARACTER( LEN=80 ) :: sub_name = 'generalized_iterative'
    !
    INTEGER, POINTER :: maxiter
    REAL( DP ), POINTER :: tolrhoaux, mix
    !
    maxiter => iterative % maxiter
    mix => iterative % mix
    tolrhoaux => iterative % tol
    !
    IF ( verbose .GE. 1 .AND. ionode ) WRITE(environ_unit,9000)
9000 FORMAT(/,4('%'),' COMPUTE ELECTROSTATIC POTENTIAL ',43('%'))

    !
    ! ... Check that fields have the same defintion domain
    !
    IF ( .NOT. ASSOCIATED(charges%cell,dielectric%epsilon%cell) ) &
         & CALL errore( sub_name, 'Inconsistent cells of input fields', 1 )
    IF ( .NOT. ASSOCIATED(charges%cell,potential%cell) ) &
         & CALL errore( sub_name, 'Inconsistent cells for charges and potential', 1 )
    cell => charges%cell
    !
    ! ... If auxiliary charge is not passed, initialize it
    !
    rhoiter => dielectric % iterative
    rhotot => dielectric % density
    CALL init_environ_density( cell, rhozero )
    !
    ! ... Aliases
    !
    eps => dielectric % epsilon
    gradlogeps => dielectric % gradlog
    !
    ! ... Set up auxiliary charge
    !
    total = integrate_environ_density( charges )
    totpol = total * ( 1.D0 - dielectric % constant ) / dielectric % constant
    jellium = 0.D0
    IF ( add_jellium ) jellium =  total / cell % omega
    rhozero % of_r = ( charges % of_r - jellium ) * ( 1.D0 - eps % of_r ) / eps % of_r
    totzero = integrate_environ_density( rhozero )
    totiter = integrate_environ_density( rhoiter )
    IF ( verbose .GE. 1 .AND. ionode ) WRITE(environ_unit,9001) totiter, jellium
9001 FORMAT(' Starting from polarization: rhoiter = ',F13.6, ' jellium = ',F13.6)

    !
    ! ... Create local variables
    !
    CALL init_environ_density( cell, residual )
    CALL init_environ_gradient( cell, gradpoisson )
    !
    ! ... Start iterative algorithm
    !
    DO iter = 1, maxiter
       !
       IF ( verbose .GE. 1 .AND. ionode ) WRITE(environ_unit,9002) iter
9002   FORMAT(' Iteration # ',i10)
       !

       rhotot % of_r = ( charges % of_r - jellium ) + rhozero % of_r + rhoiter % of_r
       !
       WRITE(environ_unit,*)"calling poisson_gradient_direct"

       CALL poisson_gradient_direct( core, rhotot, gradpoisson, electrolyte, semiconductor )
       WRITE(environ_unit,*)"finished poisson_gradient_direct"
       !
       CALL scalar_product_environ_gradient( gradlogeps, gradpoisson, residual )
       WRITE(environ_unit,*)"finished scalar product environ gradient"
       !
       residual % of_r = residual % of_r / fpi / e2 - rhoiter % of_r
       !
       rhoiter % of_r = rhoiter % of_r + mix * residual % of_r
       !
       ! ... If residual is small enough exit
       !
       delta_en = euclidean_norm_environ_density( residual )
       IF ( oldenviron ) delta_en = quadratic_mean_environ_density_old( residual )
       delta_qm = quadratic_mean_environ_density( residual )
       totiter = integrate_environ_density( rhoiter )
       IF ( verbose .GE. 1 .AND. ionode ) WRITE(environ_unit,9004)delta_qm,delta_en,tolrhoaux
9004   FORMAT(' delta_qm = ',E14.6,' delta_en = ',E14.6,' tol = ',E14.6)
       IF ( verbose .GE. 3 .AND. ionode ) WRITE(environ_unit,9003)totiter,totzero,totpol,total
9003   FORMAT(' Total iterative polarization charge = ',4F13.6)
       IF ( delta_en .LT. tolrhoaux .AND. iter .GT. 0 ) THEN
          IF ( verbose .GE. 1 .AND. ionode ) WRITE(environ_unit,9005)
9005      FORMAT(' Charges are converged, EXIT')
          EXIT
       ELSE IF ( iter .EQ. maxiter ) THEN
          IF ( ionode ) WRITE(program_unit,9006)
9006      FORMAT(' Warning: Polarization charge not converged')
       ENDIF
       !
    ENDDO
    !
    IF (.not.ltddfpt.AND.verbose.GE.1.AND.ionode) WRITE(program_unit, 9007) delta_en, iter
9007 FORMAT('     polarization accuracy =',1PE8.1,', # of iterations = ',i3)
    !
    ! ... Compute total electrostatic potential
    !
    rhotot % of_r = ( charges % of_r - jellium ) + rhozero % of_r + rhoiter % of_r
    !
    CALL poisson_direct( core, rhotot, potential, electrolyte, semiconductor )
    !
    ! ... In rhotot store total polarization charge
    !
    rhotot % of_r = rhozero % of_r + rhoiter % of_r
    !
    ! ... Destroy local variables
    !
    CALL destroy_environ_density( rhozero )
    CALL destroy_environ_density( residual )
    CALL destroy_environ_gradient( gradpoisson )
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE generalized_iterative
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE generalized_gradient_none( gradient, core, charges, dielectric,&
           & potential, electrolyte, semiconductor )
!--------------------------------------------------------------------
    !
    USE core_fft, ONLY : gradient_fft, laplacian_fft
    !
    IMPLICIT NONE
    !
    TYPE( gradient_solver ), TARGET, INTENT(IN) :: gradient
    TYPE( electrostatic_core ), INTENT(IN) :: core
    TYPE( environ_density ), TARGET, INTENT(IN) :: charges
    TYPE( environ_dielectric ), TARGET, INTENT(IN) :: dielectric
    TYPE( environ_density ), TARGET, INTENT(INOUT) :: potential
    TYPE( environ_electrolyte ), INTENT(IN), OPTIONAL :: electrolyte
    TYPE( environ_semiconductor ), INTENT(INOUT), OPTIONAL :: semiconductor
    !
    TYPE( environ_cell ), POINTER :: cell
    TYPE( environ_density ), POINTER :: x, b, eps
    TYPE( environ_gradient ), POINTER :: gradeps
    !
    INTEGER :: iter
    REAL( DP ) :: rznew, rzold, alpha, beta, pAp, delta_qm, delta_en
    TYPE( environ_density ) :: r, z, p, Ap, l
    TYPE( environ_gradient ) :: g
    !
    CHARACTER( LEN=80 ) :: sub_name = 'generalized_gradient_none'
    !
    LOGICAL, POINTER :: lconjugate
    INTEGER, POINTER :: maxstep
    REAL( DP ), POINTER :: tolvelect
    !
    lconjugate => gradient % lconjugate
    maxstep => gradient % maxstep
    tolvelect => gradient % tol
    !
    IF ( verbose .GE. 1 .AND. ionode ) WRITE(environ_unit,9000)
9000 FORMAT(/,4('%'),' COMPUTE ELECTROSTATIC POTENTIAL ',43('%'))
    !
    ! ... Check that fields have the same defintion domain
    !
    IF ( .NOT. ASSOCIATED(charges%cell,dielectric%epsilon%cell) ) &
         & CALL errore( sub_name, 'Inconsistent cells of input fields', 1 )
    IF ( .NOT. ASSOCIATED(charges%cell,potential%cell) ) &
         & CALL errore( sub_name, 'Inconsistent cells for charges and potential', 1 )
    cell => charges%cell
    !
    ! ... Aliases
    !
    b => charges
    eps => dielectric % epsilon
    gradeps => dielectric % gradient
    x => potential
    !
    ! ... Create and initialize local variables
    !
    CALL init_environ_density( cell, r )
    CALL init_environ_density( cell, z )
    CALL init_environ_density( cell, p )
    CALL init_environ_density( cell, Ap )
    !
    CALL init_environ_gradient( cell, g )
    CALL init_environ_density( cell, l )
    !
    ! ... Starting guess from new input and previous solution(s)
    !
    IF ( x%update ) THEN
       !
       x%update = .FALSE.
       !
    ENDIF
    !
    IF ( .NOT. x%update ) THEN
       !
       x%update = .TRUE.
       x%of_r = 0.D0
       r = b
       rzold = 0.D0
       !
    ENDIF
    !
    ! ... Start gradient descent
    !
    DO iter = 1, maxstep
       !
       IF ( verbose .GE. 1 .AND. ionode ) WRITE(environ_unit,9002) iter
9002   FORMAT(' Iteration # ',i10)
       !
       ! ... Apply preconditioner to new state
       !
       z%of_r = r%of_r ! no preconditioner
       !
       rznew = scalar_product_environ_density( r, z )
       IF ( ABS(rznew) .LT. 1.D-30 ) &
            & CALL errore( sub_name, 'Null step in gradient descent iteration', 1 )
       !
       ! ... Conjugate gradient or steepest descent input
       !
       IF ( lconjugate .AND. ABS(rzold) .GT. 1.D-30 ) THEN
          beta = rznew / rzold
       ELSE
          beta = 0.D0
       END IF
       p%of_r = z%of_r + beta * p%of_r
       rzold = rznew
       !
       ! ... Apply operator to conjugate direction
       ! NOTE: the following steps should be extended to account for different cores
       CALL gradient_fft( core%fft, p, g )
       CALL laplacian_fft( core%fft, p, l )
       Ap%of_r(:) = eps%of_r(:)*l%of_r(:) + &
            & gradeps%of_r(1,:)*g%of_r(1,:) + &
            & gradeps%of_r(2,:)*g%of_r(2,:) + &
            & gradeps%of_r(3,:)*g%of_r(3,:)
       Ap%of_r = - Ap%of_r / fpi / e2
       !
       ! ... Step downhill
       !
       pAp = scalar_product_environ_density( p, Ap )
       alpha = rzold / pAp
       !
       x%of_r = x%of_r + alpha * p%of_r
       r%of_r = r%of_r - alpha * Ap%of_r
       !
       IF ( verbose .GE. 1 .AND. ionode ) WRITE(environ_unit,*)'alpha = ',alpha,' beta = ',beta
       IF ( verbose .GE. 3 .AND. ionode ) WRITE(environ_unit,*)'rznew = ',rznew,' rzold = ',rzold,' pAp = ',pAp
       !
       ! ... If residual is small enough exit
       !
       delta_qm = quadratic_mean_environ_density( r )
       delta_en = euclidean_norm_environ_density( r )
       IF ( verbose .GE. 1 .AND. ionode ) WRITE(environ_unit,9004)delta_qm,delta_en,tolvelect
9004   FORMAT(' delta_qm = ',E14.6,' delta_en = ',E14.6,' tol = ',E14.6)
       IF ( delta_en .LT. tolvelect .AND. iter .GT. 0 ) THEN
          IF ( verbose .GE. 1 .AND. ionode ) WRITE(environ_unit,9005)
9005      FORMAT(' Charges are converged, EXIT')
          EXIT
       ELSE IF ( iter .EQ. maxstep ) THEN
          IF ( ionode ) WRITE(program_unit,9006)
9006      FORMAT(' Warning: Polarization charge not converged')
       ENDIF
       !
    ENDDO
    !
    IF (.not.ltddfpt.AND.verbose.GE.1.AND.ionode) WRITE(program_unit, 9007) delta_en, iter
9007 FORMAT('     polarization accuracy =',1PE8.1,', # of iterations = ',i3)
    !
    CALL destroy_environ_density( l )
    CALL destroy_environ_gradient( g )
    !
    CALL destroy_environ_density( r )
    CALL destroy_environ_density( z )
    CALL destroy_environ_density( p )
    CALL destroy_environ_density( Ap )
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE generalized_gradient_none
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE generalized_gradient_sqrt( gradient, core, charges, dielectric, potential,&
               & electrolyte, semiconductor )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    TYPE( gradient_solver ), TARGET, INTENT(IN) :: gradient
    TYPE( electrostatic_core ), INTENT(IN) :: core
    TYPE( environ_density ), TARGET, INTENT(IN) :: charges
    TYPE( environ_dielectric ), TARGET, INTENT(IN) :: dielectric
    TYPE( environ_density ), TARGET, INTENT(INOUT) :: potential
    TYPE( environ_electrolyte ), INTENT(IN), OPTIONAL :: electrolyte
    TYPE( environ_semiconductor ), INTENT(INOUT), OPTIONAL :: semiconductor
    !
    TYPE( environ_cell ), POINTER :: cell
    TYPE( environ_density ), POINTER :: x, b, eps, factsqrt
    TYPE( environ_gradient ), POINTER :: gradeps
    !
    INTEGER :: iter
    REAL( DP ) :: rznew, rzold, alpha, beta, pAp, delta_qm, delta_en, jellium, shift
    TYPE( environ_density ) :: r, z, p, Ap, invsqrt
    !
    CHARACTER( LEN=80 ) :: sub_name = 'generalized_gradient_sqrt'
    !
    LOGICAL, POINTER :: lconjugate
    INTEGER, POINTER :: maxstep
    REAL( DP ), POINTER :: tolvelect
    !
    lconjugate => gradient % lconjugate
    maxstep => gradient % maxstep
    tolvelect => gradient % tol
    !
    IF ( verbose .GE. 1 .AND. ionode ) WRITE(environ_unit,9000)
9000 FORMAT(/,4('%'),' COMPUTE ELECTROSTATIC POTENTIAL ',43('%'))
    !
    ! ... Check that fields have the same defintion domain
    !
    IF ( .NOT. ASSOCIATED(charges%cell,dielectric%epsilon%cell) ) &
         & CALL errore( sub_name, 'Inconsistent cells of input fields', 1 )
    IF ( .NOT. ASSOCIATED(charges%cell,potential%cell) ) &
         & CALL errore( sub_name, 'Inconsistent cells for charges and potential', 1 )
    cell => charges%cell
    !
    ! ... Aliases
    !
    b => charges
    eps => dielectric % epsilon
    factsqrt => dielectric % factsqrt
    x => potential
    CALL init_environ_density( cell, invsqrt )
    invsqrt%of_r = 1.D0 / SQRT(eps%of_r)
    jellium = 0.D0
    IF ( add_jellium ) jellium = integrate_environ_density( charges ) / cell % omega
    !
    ! ... Create and initialize local variables
    !
    CALL init_environ_density( cell, r )
    CALL init_environ_density( cell, z )
    CALL init_environ_density( cell, p )
    CALL init_environ_density( cell, Ap )
    !
    ! ... Starting guess from new input and previous solution(s)
    !
    IF ( x%update ) THEN
       !
       r%of_r = ( b%of_r - jellium ) - factsqrt%of_r * x%of_r
       !
       ! ... Preconditioning step
       !
       z%of_r = r%of_r * invsqrt%of_r
       CALL poisson_direct( core, z, z, electrolyte, semiconductor )
       z%of_r = z%of_r * invsqrt%of_r
       !
       rzold = scalar_product_environ_density( r, z )
       IF ( ABS(rzold) .LT. 1.D-30 ) &
            & CALL errore( sub_name, 'Null step in gradient descent iteration', 1 )
       !
       r%of_r = factsqrt%of_r * ( x%of_r - z%of_r )
       delta_en = euclidean_norm_environ_density( r )
       delta_qm = quadratic_mean_environ_density( r )
       IF ( delta_en .LT. 1.D-02 ) THEN
          IF ( verbose .GE. 1 .AND. ionode ) WRITE(environ_unit,9008)delta_en
9008      FORMAT(' Sqrt-preconditioned input guess with residual norm = ',E14.6)
          x%of_r = z%of_r
       ELSE
          IF ( verbose .GE. 1 .AND. ionode ) WRITE(environ_unit,9001)delta_en
9001      FORMAT(' Warning: bad guess with residual norm = ',E14.6,', reset to no guess')
          x%update = .FALSE.
       ENDIF
       !
    ENDIF
    !
    IF ( .NOT. x%update ) THEN
       !
       x%update = .TRUE.
       x%of_r = 0.D0
       r%of_r = b%of_r - jellium
       rzold = 0.D0
       !
    ENDIF
    !
    ! ... Start gradient descent
    !
    DO iter = 1, maxstep
       !
       IF ( verbose .GE. 1 .AND. ionode ) WRITE(environ_unit,9002) iter
9002   FORMAT(' Iteration # ',i10)
       !
       ! ... Apply preconditioner to new state
       !
       z%of_r = r%of_r * invsqrt%of_r
       CALL poisson_direct( core, z, z, electrolyte, semiconductor )
       z%of_r = z%of_r * invsqrt%of_r
       !
       rznew = scalar_product_environ_density( r, z )
       IF ( ABS(rznew) .LT. 1.D-30 ) &
            & CALL errore( sub_name, 'Null step in gradient descent iteration', 1 )
       !
       ! ... Conjugate gradient or steepest descent input
       !
       IF ( lconjugate .AND. ABS(rzold) .GT. 1.D-30 ) THEN
          beta = rznew / rzold
       ELSE
          beta = 0.D0
       END IF
       IF ( verbose .GE. 3 .AND. ionode ) WRITE(environ_unit,*)'rznew = ',rznew,' rzold = ',rzold,' beta = ',beta
       rzold = rznew
       !
       p%of_r = z%of_r + beta * p%of_r
       !
       ! ... Apply operator to conjugate direction
       !
       Ap%of_r = factsqrt%of_r * z%of_r + r%of_r + beta * Ap%of_r
       !
       ! ... Step downhill
       !
       pAp = scalar_product_environ_density( p, Ap )
       alpha = rzold / pAp
       IF ( verbose .GE. 1 .AND. ionode ) WRITE(environ_unit,*)' pAp = ',pAp,' rzold = ',rzold,' alpha = ',alpha
       !
       x%of_r = x%of_r + alpha * p%of_r
       r%of_r = r%of_r - alpha * Ap%of_r
       !
       ! ... If residual is small enough exit
       !
       delta_qm = quadratic_mean_environ_density( r )
       delta_en = euclidean_norm_environ_density( r )
       IF ( verbose .GE. 1 .AND. ionode ) WRITE(environ_unit,9004)delta_qm,delta_en,tolvelect
9004   FORMAT(' delta_qm = ',E14.6,' delta_en = ',E14.6,' tol = ',E14.6)
       IF ( delta_en .LT. tolvelect .AND. iter .GT. 0 ) THEN
          IF ( verbose .GE. 1 .AND. ionode ) WRITE(environ_unit,9005)
9005      FORMAT(' Charges are converged, EXIT')
          EXIT
       ELSE IF ( iter .EQ. maxstep ) THEN
          IF ( ionode ) WRITE(program_unit,9006)
9006      FORMAT(' Warning: Polarization charge not converged')
       ENDIF
       !
    ENDDO
    !
    ! In PBC the potential need to have zero average
    !
    shift = 0.D0
    !
    IF ( core % use_fft ) THEN
       !
       IF ( .NOT. ( core % fft % use_internal_pbc_corr .OR. core % need_correction ) ) THEN
          !
          shift = - integrate_environ_density( x ) / cell % omega
          !
       END IF
       !
    END IF
    !
    x % of_r = x % of_r + shift
    !
    IF (.not.ltddfpt.AND.verbose.GE.1.AND.ionode) WRITE(program_unit, 9007) delta_en, iter
9007 FORMAT('     polarization accuracy =',1PE8.1,', # of iterations = ',i3)
    !
    CALL destroy_environ_density( r )
    CALL destroy_environ_density( z )
    CALL destroy_environ_density( p )
    CALL destroy_environ_density( Ap )
    !
    CALL destroy_environ_density(invsqrt)
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE generalized_gradient_sqrt
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE generalized_gradient_left( gradient, core, charges, dielectric, potential ,&
              & electrolyte, semiconductor )
!--------------------------------------------------------------------
    !
    USE core_fft, ONLY : gradient_fft
    !
    IMPLICIT NONE
    !
    TYPE( gradient_solver ), TARGET, INTENT(IN) :: gradient
    TYPE( electrostatic_core ), INTENT(IN) :: core
    TYPE( environ_density ), TARGET, INTENT(IN) :: charges
    TYPE( environ_dielectric ), TARGET, INTENT(IN) :: dielectric
    TYPE( environ_density ), TARGET, INTENT(INOUT) :: potential
    TYPE( environ_electrolyte ), INTENT(IN), OPTIONAL :: electrolyte
    TYPE( environ_semiconductor ), INTENT(INOUT), OPTIONAL :: semiconductor
    !
    TYPE( environ_cell ), POINTER :: cell
    TYPE( environ_density ), POINTER :: x, b, eps
    TYPE( environ_gradient ), POINTER :: gradeps
    !
    INTEGER :: iter
    REAL( DP ) :: rznew, rzold, alpha, beta, pAp, delta_en, delta_qm
    TYPE( environ_density ) :: r, z, p, Ap
    TYPE( environ_gradient ) :: g
    !
    CHARACTER( LEN=80 ) :: sub_name = 'generalized_gradient_left'
    !
    LOGICAL, POINTER :: lconjugate
    INTEGER, POINTER :: maxstep
    REAL( DP ), POINTER :: tolvelect
    !
    lconjugate => gradient % lconjugate
    maxstep => gradient % maxstep
    tolvelect => gradient % tol
    !
    IF ( verbose .GE. 1 .AND. ionode ) WRITE(environ_unit,9000)
9000 FORMAT(/,4('%'),' COMPUTE ELECTROSTATIC POTENTIAL ',43('%'))
    !
    ! ... Check that fields have the same defintion domain
    !
    IF ( .NOT. ASSOCIATED(charges%cell,dielectric%epsilon%cell) ) &
         & CALL errore( sub_name, 'Inconsistent cells of input fields', 1 )
    IF ( .NOT. ASSOCIATED(charges%cell,potential%cell) ) &
         & CALL errore( sub_name, 'Inconsistent cells for charges and potential', 1 )
    cell => charges%cell
    !
    ! ... Aliases
    !
    b => charges
    eps => dielectric % epsilon
    gradeps => dielectric % gradient
    x => potential
    !
    ! ... Create and initialize local variables
    !
    CALL init_environ_density( cell, r )
    CALL init_environ_density( cell, z )
    CALL init_environ_density( cell, p )
    CALL init_environ_density( cell, Ap )
    !
    CALL init_environ_gradient( cell, g )
    !
    ! ... Starting guess from new input and previous solution(s)
    !
!!!  IF ( x%update ) THEN
!!!
!!!     CALL gradient_fft( core%fft, x, g )
!!!     g%of_r = - g%of_r / fpi / e2
!!!     r%of_r(:) = b%of_r -  &
!!!                & gradeps%of_r(1,:)*g%of_r(1,:) + &
!!!                & gradeps%of_r(2,:)*g%of_r(2,:) + &
!!!                & gradeps%of_r(3,:)*g%of_r(3,:)
!!!
!!!     ! ... Preconditioning step
!!!
!!!     z%of_r = r%of_r / eps%of_r
!!!     CALL poisson_direct( core, z, z )
!!!
!!!     rzold = scalar_product_environ_density( r, z )
!!!     IF ( ABS(rzold) .LT. 1.D-30 ) &
!!!          & CALL errore(sub_name,'Null step in gradient descent iteration',1)
!!!
!!!     r%of_r = x%of_r - z%of_r
!!!     CALL gradient_fft( core%fft, r, g )
!!!     g%of_r = - g%of_r / fpi / e2
!!!     r%of_r(:) = gradeps%of_r(1,:)*g%of_r(1,:) + &
!!!               & gradeps%of_r(2,:)*g%of_r(2,:) + &
!!!               & gradeps%of_r(3,:)*g%of_r(3,:)
!!!     delta_qm = quadratic_mean_environ_density( r )
!!!     IF ( delta_qm .LT. 1.D-02 ) THEN
!!!        IF ( verbose .GE. 1 .AND. ionode ) WRITE(environ_unit,9008)delta_qm
!!!9008    FORMAT(' Sqrt-preconditioned input guess with residual norm = ',E14.6)
!!!        x%of_r = z%of_r
!!!     ELSE
!!!        IF ( verbose .GE. 1 .AND. ionode ) WRITE(environ_unit,9001)delta_qm
!!!9001    FORMAT(' Warning: bad guess with residual norm = ',E14.6,', reset to no guess')
!!!        x%update = .FALSE.
!!!     ENDIF
!!!
!!!  ENDIF
!!!
!!!  IF ( .NOT. x%update ) THEN
!!!
!!!     x%update = .TRUE.
    x%of_r = 0.D0
    r%of_r = b%of_r
    rzold = 0.D0
    !
!!!  ENDIF
    !
    ! ... Start gradient descent
    !
    DO iter = 1, maxstep
       !
       IF ( verbose .GE. 1 .AND. ionode ) WRITE(environ_unit,9002) iter
9002   FORMAT(' Iteration # ',i10)
       !
       ! ... Apply preconditioner to new state
       !
       z%of_r = r%of_r / eps%of_r
       CALL poisson_direct( core, z, z, electrolyte, semiconductor )
       !
       rznew = scalar_product_environ_density( r, z )
       IF ( ABS(rznew) .LT. 1.D-30 ) &
            & CALL errore( sub_name, 'Null step in gradient descent iteration', 1 )
       !
       ! ... Conjugate gradient or steepest descent input
       !
       IF ( lconjugate .AND. ABS(rzold) .GT. 1.D-30 ) THEN
          beta = rznew / rzold
       ELSE
          beta = 0.D0
       END IF
       IF ( verbose .GE. 3 .AND. ionode ) WRITE(environ_unit,*)'rznew = ',rznew,' rzold = ',rzold,' beta = ',beta
       rzold = rznew
       !
       p%of_r = z%of_r + beta * p%of_r
       !
       ! ... Apply operator to conjugate direction
       ! NOTE: the following steps should be extended to account for different cores
       CALL gradient_fft( core%fft, z, g )
       g%of_r = g%of_r / fpi / e2
       Ap%of_r(:) = beta * Ap%of_r(:) - r%of_r(:) + &
                  & gradeps%of_r(1,:)*g%of_r(1,:) + &
                  & gradeps%of_r(2,:)*g%of_r(2,:) + &
                  & gradeps%of_r(3,:)*g%of_r(3,:)
       !
       ! ... Step downhill
       !
       pAp = scalar_product_environ_density( p, Ap )
       alpha = rzold / pAp
       IF ( verbose .GE. 1 .AND. ionode ) WRITE(environ_unit,*)' pAp = ',pAp,' rzold = ',rzold,' alpha = ',alpha
       !
       x%of_r = x%of_r + alpha * p%of_r
       r%of_r = r%of_r - alpha * Ap%of_r
       !
       ! ... If residual is small enough exit
       !
       delta_qm = quadratic_mean_environ_density( r )
       delta_en = euclidean_norm_environ_density( r )
       IF ( verbose .GE. 1 .AND. ionode ) WRITE(environ_unit,9004)delta_qm,delta_en,tolvelect
9004   FORMAT(' delta_qm = ',E14.6,' delta_en = ',E14.6,' tol = ',E14.6)
       IF ( delta_en .LT. tolvelect .AND. iter .GT. 0 ) THEN
          IF ( verbose .GE. 1 .AND. ionode ) WRITE(environ_unit,9005)
9005      FORMAT(' Charges are converged, EXIT')
          EXIT
       ELSE IF ( iter .EQ. maxstep ) THEN
         IF ( ionode ) WRITE(program_unit,9006)
9006     FORMAT(' Warning: Polarization charge not converged')
       ENDIF
       !
    ENDDO
    !
    IF (.not.ltddfpt.AND.verbose.GE.1.AND.ionode) WRITE(program_unit, 9007) delta_en, iter
9007 FORMAT('     polarization accuracy =',1PE8.1,', # of iterations = ',i3)
    !
    CALL destroy_environ_gradient( g )
    !
    CALL destroy_environ_density( r )
    CALL destroy_environ_density( z )
    CALL destroy_environ_density( p )
    CALL destroy_environ_density( Ap )
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE generalized_gradient_left
!--------------------------------------------------------------------
!----------------------------------------------------------------------------
END MODULE problem_generalized
!----------------------------------------------------------------------------
