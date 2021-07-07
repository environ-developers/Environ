!----------------------------------------------------------------------------------------
!
! Copyright (C) 2018-2021 ENVIRON (www.quantum-environ.org)
!
!----------------------------------------------------------------------------------------
!
!     This file is part of Environ version 2.0
!     
!     Environ 2.0 is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation, either version 2 of the License, or
!     (at your option) any later version.
!     
!     Environ 2.0 is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more detail, either the file
!     `License' in the root directory of the present distribution, or
!     online at <http://www.gnu.org/licenses/>.
!
!----------------------------------------------------------------------------------------
!
! Authors: Oliviero Andreussi (Department of Physics, UNT)
!          Nicola Marzari     (THEOS and NCCR-MARVEL, EPFL)
!          Edan Bainglass     (Department of Physics, UNT)
!
! includes improved algorithms from G. Fisicaro and S. Goedecker
!
!----------------------------------------------------------------------------------------
!>
!! This module contains the main drivers and routines to compute the
!! electrostatic potential that is the solution of a
!! generalized Poisson equation:
!! \f[
!!      \nabla \cdot \epsilon (r) \nabla \phi = -4 \pi \rho
!! \f]
!!
!! Different algorithms (gradient descent on potential and iterative
!! on the polarization charge) are available and implemented.
!!
!----------------------------------------------------------------------------------------
MODULE solver_generalized
    !------------------------------------------------------------------------------------
    !
    USE env_base_io, ONLY: ionode, environ_unit, verbose, program_unit, lstdout
    !
    USE environ_param, ONLY: DP, e2, fpi
    !
    USE types_electrostatic, ONLY: electrostatic_solver, iterative_solver, &
                                   gradient_solver
    !
    USE types_physical, ONLY: environ_charges, environ_dielectric, &
                              environ_electrolyte, environ_semiconductor
    !
    USE types_core, ONLY: core_container
    USE types_representation, ONLY: environ_density, environ_gradient
    USE types_cell, ONLY: environ_cell
    !
    USE utils_density, ONLY: init_environ_density, destroy_environ_density
    USE utils_gradient, ONLY: init_environ_gradient, destroy_environ_gradient
    !
    USE tools_math, ONLY: scalar_product_environ_density, &
                          scalar_product_environ_gradient, &
                          euclidean_norm_environ_density, &
                          quadratic_mean_environ_density, &
                          integrate_environ_density
    !
    USE tools_fft, ONLY: gradient_fft, laplacian_fft
    !
    USE solver_poisson
    !
    !------------------------------------------------------------------------------------
    !
    INTERFACE generalized_gradient
        MODULE PROCEDURE generalized_gradient_charges, generalized_gradient_density
    END INTERFACE generalized_gradient
    !
    !------------------------------------------------------------------------------------
    !
    PRIVATE
    !
    PUBLIC :: generalized_gradient
    !
    !------------------------------------------------------------------------------------
CONTAINS
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE generalized_gradient_charges(solver, core, charges, potential)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(electrostatic_solver), INTENT(IN) :: solver
        TYPE(core_container), INTENT(IN) :: core
        !
        TYPE(environ_charges), INTENT(INOUT) :: charges
        TYPE(environ_density), INTENT(INOUT) :: potential
        !
        CHARACTER*20 :: sub_name = 'generalized_gradient_charges'
        !
        !--------------------------------------------------------------------------------
        !
        CALL env_start_clock(sub_name)
        !
        potential%of_r = 0.D0
        !
        IF (solver%use_gradient) THEN
            !
            IF (solver%auxiliary == 'none') THEN
                !
                SELECT CASE (solver%gradient%preconditioner) ! #TODO none and left need work
                    !
                CASE ('none')
                    !
                    CALL generalized_gradient_none(solver%gradient, core, &
                                                   charges%density, charges%dielectric, &
                                                   potential, charges%electrolyte, &
                                                   charges%semiconductor)
                    !
                CASE ('sqrt')
                    !
                    CALL generalized_gradient_sqrt(solver%gradient, core, &
                                                   charges%density, charges%dielectric, &
                                                   potential, charges%electrolyte, &
                                                   charges%semiconductor)
                    !
                CASE ('left')
                    !
                    CALL generalized_gradient_left(solver%gradient, core, &
                                                   charges%density, charges%dielectric, &
                                                   potential, charges%electrolyte, &
                                                   charges%semiconductor)
                    !
                CASE DEFAULT
                    CALL env_errore(sub_name, 'Unexpected preconditioner keyword', 1)
                    !
                END SELECT
                !
            ELSE
                !
                CALL env_errore(sub_name, 'Option not yet implemented', 1)
                !
                ! CALL generalized_gradient_rhoaux(charges, dielectric, potential) #TODO future-work
                !
            END IF
            !
        ELSE IF (solver%use_iterative) THEN
            !
            IF (solver%auxiliary == 'full') THEN
                !
                CALL generalized_iterative(solver%iterative, core, charges%density, &
                                           charges%dielectric, potential, &
                                           charges%electrolyte, charges%semiconductor)
                !
            ELSE
                !
                CALL env_errore(sub_name, 'Option not yet implemented', 1)
                !
                ! CALL generalized_iterative_velect(charges, dielectric, potential) #TODO future-work
                !
            END IF
            !
        ELSE
            CALL env_errore(sub_name, 'Unexpected auxiliary keyword', 1)
        END IF
        !
        CALL env_stop_clock(sub_name)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE generalized_gradient_charges
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE generalized_gradient_density(solver, core, charges, dielectric, &
                                            potential, electrolyte, semiconductor)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(electrostatic_solver), INTENT(IN) :: solver
        TYPE(core_container), INTENT(IN) :: core
        TYPE(environ_density), INTENT(IN) :: charges
        TYPE(environ_dielectric), INTENT(IN) :: dielectric
        TYPE(environ_electrolyte), INTENT(IN), OPTIONAL :: electrolyte
        !
        TYPE(environ_density), INTENT(INOUT) :: potential
        TYPE(environ_semiconductor), INTENT(INOUT), OPTIONAL :: semiconductor
        !
        CHARACTER*20 :: sub_name = 'generalized_gradient_density'
        !
        !--------------------------------------------------------------------------------
        !
        CALL env_start_clock(sub_name)
        !
        potential%of_r = 0.D0
        !
        IF (solver%use_gradient) THEN
            !
            IF (solver%auxiliary == 'none') THEN
                !
                SELECT CASE (solver%gradient%preconditioner)
                    !
                CASE ('none')
                    !
                    CALL generalized_gradient_none(solver%gradient, core, charges, &
                                                   dielectric, potential, electrolyte)
                    !
                CASE ('sqrt')
                    !
                    CALL generalized_gradient_sqrt(solver%gradient, core, charges, &
                                                   dielectric, potential, electrolyte)
                    !
                CASE ('left')
                    !
                    CALL generalized_gradient_left(solver%gradient, core, charges, &
                                                   dielectric, potential, electrolyte)
                    !
                CASE DEFAULT
                    CALL env_errore(sub_name, 'Unexpected preconditioner keyword', 1)
                    !
                END SELECT
                !
            ELSE
                !
                CALL env_errore(sub_name, 'Option not yet implemented', 1)
                !
                ! CALL generalized_gradient_rhoaux(charges, dielectric, potential) #TODO future-work
                !
            END IF
            !
        ELSE IF (solver%use_iterative) THEN
            !
            IF (solver%auxiliary == 'full') THEN
                !
                CALL generalized_iterative(solver%iterative, core, charges, &
                                           dielectric, potential, electrolyte, &
                                           semiconductor)
                !
            ELSE
                !
                CALL env_errore(sub_name, 'Option not yet implemented', 1)
                !
                ! CALL generalized_iterative_velect(charges, dielectric, potential) #TODO future-work
                !
            END IF
            !
        ELSE
            !
            CALL env_errore(sub_name, 'Unexpected auxiliary keyword', 1)
            !
        END IF
        !
        CALL env_stop_clock(sub_name)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE generalized_gradient_density
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE generalized_iterative(iterative, core, charges, dielectric, &
                                     potential, electrolyte, semiconductor)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(iterative_solver), TARGET, INTENT(IN) :: iterative
        TYPE(core_container), INTENT(IN) :: core
        TYPE(environ_density), TARGET, INTENT(IN) :: charges
        TYPE(environ_dielectric), TARGET, INTENT(IN) :: dielectric
        TYPE(environ_electrolyte), INTENT(IN), OPTIONAL :: electrolyte
        !
        TYPE(environ_density), TARGET, INTENT(INOUT) :: potential
        TYPE(environ_semiconductor), INTENT(INOUT), OPTIONAL :: semiconductor
        !
        TYPE(environ_cell), POINTER :: cell
        TYPE(environ_density), POINTER :: eps, rhoiter, rhotot
        TYPE(environ_gradient), POINTER :: gradlogeps
        !
        INTEGER :: iter
        REAL(DP) :: total, totpol, totzero, totiter, delta_qm, delta_en
        TYPE(environ_density) :: rhozero
        TYPE(environ_density) :: residual
        TYPE(environ_gradient) :: gradpoisson
        !
        INTEGER, POINTER :: maxiter
        REAL(DP), POINTER :: tolrhoaux, mix
        !
        CHARACTER(LEN=80) :: sub_name = 'generalized_iterative'
        !
        !--------------------------------------------------------------------------------
        !
        maxiter => iterative%maxiter
        mix => iterative%mix
        tolrhoaux => iterative%tol
        !
        IF (verbose >= 1 .AND. ionode) WRITE (environ_unit, 9000)
        !
9000    FORMAT(/, 4('%'), ' COMPUTE ELECTROSTATIC POTENTIAL ', 43('%'))
        !
        !--------------------------------------------------------------------------------
        ! Check that fields have the same defintion domain
        !
        IF (.NOT. ASSOCIATED(charges%cell, dielectric%epsilon%cell)) &
            CALL env_errore(sub_name, 'Inconsistent cells of input fields', 1)
        !
        IF (.NOT. ASSOCIATED(charges%cell, potential%cell)) &
            CALL env_errore(sub_name, 'Inconsistent cells for charges and potential', 1)
        !
        cell => charges%cell
        !
        !--------------------------------------------------------------------------------
        ! If auxiliary charge is not passed, initialize it
        !
        rhoiter => dielectric%iterative
        rhotot => dielectric%density
        !
        CALL init_environ_density(cell, rhozero)
        !
        eps => dielectric%epsilon
        gradlogeps => dielectric%gradlog
        !
        !--------------------------------------------------------------------------------
        ! Set up auxiliary charge
        !
        total = integrate_environ_density(charges)
        totpol = total * (1.D0 - dielectric%constant) / dielectric%constant
        rhozero%of_r = charges%of_r * (1.D0 - eps%of_r) / eps%of_r
        totzero = integrate_environ_density(rhozero)
        totiter = integrate_environ_density(rhoiter)
        !
        IF (verbose >= 1 .AND. ionode) WRITE (environ_unit, 9001) totiter
        !
9001    FORMAT(' Starting from polarization: rhoiter = ', F13.6)
        !
        !--------------------------------------------------------------------------------
        ! Create local variables
        !
        CALL init_environ_density(cell, residual)
        !
        CALL init_environ_gradient(cell, gradpoisson)
        !
        !--------------------------------------------------------------------------------
        ! Start iterative algorithm
        !
        DO iter = 1, maxiter
            !
            IF (verbose >= 1 .AND. ionode) WRITE (environ_unit, 9002) iter
            !
9002        FORMAT(' Iteration # ', i10)
            !
            rhotot%of_r = charges%of_r + rhozero%of_r + rhoiter%of_r
            !
            CALL poisson_gradient_direct(core, rhotot, gradpoisson, electrolyte, &
                                         semiconductor)
            !
            CALL scalar_product_environ_gradient(gradlogeps, gradpoisson, residual)
            !
            residual%of_r = residual%of_r / fpi / e2 - rhoiter%of_r
            !
            rhoiter%of_r = rhoiter%of_r + mix * residual%of_r
            !
            !----------------------------------------------------------------------------
            ! If residual is small enough exit
            !
            delta_en = euclidean_norm_environ_density(residual)
            delta_qm = quadratic_mean_environ_density(residual)
            totiter = integrate_environ_density(rhoiter)
            !
            IF (verbose >= 1 .AND. ionode) &
                WRITE (environ_unit, 9004) delta_qm, delta_en, tolrhoaux
            !
9004        FORMAT(' delta_qm = ', E14.6, ' delta_en = ', E14.6, ' tol = ', E14.6)
            !
            IF (verbose >= 3 .AND. ionode) &
                WRITE (environ_unit, 9003) totiter, totzero, totpol, total
            !
9003        FORMAT(' Total iterative polarization charge = ', 4F13.6)
            !
            IF (delta_en < tolrhoaux .AND. iter > 0) THEN
                !
                IF (verbose >= 1 .AND. ionode) WRITE (environ_unit, 9005)
                !
9005            FORMAT(' Charges are converged, EXIT')
                !
                EXIT
                !
            ELSE IF (iter == maxiter) THEN
                !
                IF (ionode) WRITE (program_unit, 9006)
                !
9006            FORMAT(' Warning: Polarization charge not converged')
            END IF
            !
        END DO
        !
        IF (lstdout .AND. verbose >= 1) WRITE (program_unit, 9007) delta_en, iter
        !
9007    FORMAT('     Polarization accuracy =', 1PE8.1, ', # of iterations = ', i3)
        !
        !--------------------------------------------------------------------------------
        ! Compute total electrostatic potential
        !
        rhotot%of_r = charges%of_r + rhozero%of_r + rhoiter%of_r
        !
        CALL poisson_direct(core, rhotot, potential, electrolyte, semiconductor)
        !
        !--------------------------------------------------------------------------------
        !
        rhotot%of_r = rhozero%of_r + rhoiter%of_r
        ! in rhotot store total polarization charge
        !
        !--------------------------------------------------------------------------------
        ! Destroy local variables
        !
        CALL destroy_environ_density(rhozero)
        !
        CALL destroy_environ_density(residual)
        !
        CALL destroy_environ_gradient(gradpoisson)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE generalized_iterative
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE generalized_gradient_none(gradient, core, charges, dielectric, &
                                         potential, electrolyte, semiconductor)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(gradient_solver), TARGET, INTENT(IN) :: gradient
        TYPE(core_container), INTENT(IN) :: core
        TYPE(environ_density), TARGET, INTENT(IN) :: charges
        TYPE(environ_dielectric), TARGET, INTENT(IN) :: dielectric
        TYPE(environ_electrolyte), INTENT(IN), OPTIONAL :: electrolyte
        !
        TYPE(environ_density), TARGET, INTENT(INOUT) :: potential
        TYPE(environ_semiconductor), INTENT(INOUT), OPTIONAL :: semiconductor
        !
        TYPE(environ_cell), POINTER :: cell
        TYPE(environ_density), POINTER :: x, b, eps
        TYPE(environ_gradient), POINTER :: gradeps
        !
        INTEGER :: iter
        REAL(DP) :: rznew, rzold, alpha, beta, pAp, delta_qm, delta_en
        TYPE(environ_density) :: r, z, p, Ap, l
        TYPE(environ_gradient) :: g
        !
        LOGICAL, POINTER :: lconjugate
        INTEGER, POINTER :: maxstep
        REAL(DP), POINTER :: tolvelect
        !
        CHARACTER(LEN=80) :: sub_name = 'generalized_gradient_none'
        !
        !--------------------------------------------------------------------------------
        !
        lconjugate => gradient%lconjugate
        maxstep => gradient%maxstep
        tolvelect => gradient%tol
        !
        IF (verbose >= 1 .AND. ionode) WRITE (environ_unit, 9000)
        !
9000    FORMAT(/, 4('%'), ' COMPUTE ELECTROSTATIC POTENTIAL ', 43('%'))
        !
        !--------------------------------------------------------------------------------
        ! Check that fields have the same defintion domain
        !
        IF (.NOT. ASSOCIATED(charges%cell, dielectric%epsilon%cell)) &
            CALL env_errore(sub_name, 'Inconsistent cells of input fields', 1)
        !
        IF (.NOT. ASSOCIATED(charges%cell, potential%cell)) &
            CALL env_errore(sub_name, 'Inconsistent cells for charges and potential', 1)
        !
        cell => charges%cell
        !
        b => charges
        eps => dielectric%epsilon
        gradeps => dielectric%gradient
        x => potential
        !
        CALL init_environ_density(cell, r)
        !
        CALL init_environ_density(cell, z)
        !
        CALL init_environ_density(cell, p)
        !
        CALL init_environ_density(cell, Ap)
        !
        CALL init_environ_gradient(cell, g)
        !
        CALL init_environ_density(cell, l)
        !
        !--------------------------------------------------------------------------------
        ! Starting guess from new input and previous solution(s)
        !
        IF (x%update) THEN
            x%update = .FALSE.
        END IF
        !
        IF (.NOT. x%update) THEN
            x%update = .TRUE.
            x%of_r = 0.D0
            r = b
            rzold = 0.D0
        END IF
        !
        !--------------------------------------------------------------------------------
        ! Start gradient descent
        !
        DO iter = 1, maxstep
            !
            IF (verbose >= 1 .AND. ionode) WRITE (environ_unit, 9002) iter
            !
9002        FORMAT(' Iteration # ', i10)
            !
            !----------------------------------------------------------------------------
            ! Apply preconditioner to new state
            !
            z%of_r = r%of_r ! no preconditioner
            !
            rznew = scalar_product_environ_density(r, z)
            !
            IF (ABS(rznew) < 1.D-30) &
                CALL env_errore(sub_name, 'Null step in gradient descent iteration', 1)
            !
            !----------------------------------------------------------------------------
            ! Conjugate gradient or steepest descent input
            !
            IF (lconjugate .AND. ABS(rzold) > 1.D-30) THEN
                beta = rznew / rzold
            ELSE
                beta = 0.D0
            END IF
            !
            p%of_r = z%of_r + beta * p%of_r
            rzold = rznew
            !
            !----------------------------------------------------------------------------
            ! Apply operator to conjugate direction
            ! NOTE: the following steps should be extended to account for different cores
            !
            CALL gradient_fft(core%fft, p, g)
            !
            CALL laplacian_fft(core%fft, p, l)
            !
            Ap%of_r(:) = eps%of_r(:) * l%of_r(:) + &
                         gradeps%of_r(1, :) * g%of_r(1, :) + &
                         gradeps%of_r(2, :) * g%of_r(2, :) + &
                         gradeps%of_r(3, :) * g%of_r(3, :)
            !
            Ap%of_r = -Ap%of_r / fpi / e2
            !
            !----------------------------------------------------------------------------
            ! Step downhill
            !
            pAp = scalar_product_environ_density(p, Ap)
            alpha = rzold / pAp
            !
            x%of_r = x%of_r + alpha * p%of_r
            r%of_r = r%of_r - alpha * Ap%of_r
            !
            IF (verbose >= 1 .AND. ionode) &
                WRITE (environ_unit, *) 'alpha = ', alpha, ' beta = ', beta
            !
            IF (verbose >= 3 .AND. ionode) &
                WRITE (environ_unit, *) &
                'rznew = ', rznew, ' rzold = ', rzold, ' pAp = ', pAp
            !
            !----------------------------------------------------------------------------
            ! If residual is small enough exit
            !
            delta_qm = quadratic_mean_environ_density(r)
            delta_en = euclidean_norm_environ_density(r)
            !
            IF (verbose >= 1 .AND. ionode) &
                WRITE (environ_unit, 9004) delta_qm, delta_en, tolvelect
            !
9004        FORMAT(' delta_qm = ', E14.6, ' delta_en = ', E14.6, ' tol = ', E14.6)
            !
            IF (delta_en < tolvelect .AND. iter > 0) THEN
                !
                IF (verbose >= 1 .AND. ionode) WRITE (environ_unit, 9005)
                !
9005            FORMAT(' Charges are converged, EXIT')
                !
                EXIT
                !
            ELSE IF (iter == maxstep) THEN
                !
                IF (ionode) WRITE (program_unit, 9006)
                !
9006            FORMAT(' Warning: Polarization charge not converged')
            END IF
            !
        END DO
        !
        IF (lstdout .AND. verbose >= 1) WRITE (program_unit, 9007) delta_en, iter
        !
9007    FORMAT('     Polarization accuracy =', 1PE8.1, ', # of iterations = ', i3)
        !
        CALL destroy_environ_density(l)
        !
        CALL destroy_environ_gradient(g)
        !
        CALL destroy_environ_density(r)
        !
        CALL destroy_environ_density(z)
        !
        CALL destroy_environ_density(p)
        !
        CALL destroy_environ_density(Ap)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE generalized_gradient_none
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE generalized_gradient_sqrt(gradient, core, charges, dielectric, &
                                         potential, electrolyte, semiconductor)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(gradient_solver), TARGET, INTENT(IN) :: gradient
        TYPE(core_container), INTENT(IN) :: core
        TYPE(environ_density), TARGET, INTENT(IN) :: charges
        TYPE(environ_dielectric), TARGET, INTENT(IN) :: dielectric
        TYPE(environ_electrolyte), INTENT(IN), OPTIONAL :: electrolyte
        !
        TYPE(environ_density), TARGET, INTENT(INOUT) :: potential
        TYPE(environ_semiconductor), INTENT(INOUT), OPTIONAL :: semiconductor
        !
        TYPE(environ_cell), POINTER :: cell
        TYPE(environ_density), POINTER :: x, b, eps, factsqrt
        TYPE(environ_gradient), POINTER :: gradeps
        !
        INTEGER :: iter
        REAL(DP) :: rznew, rzold, alpha, beta, pAp, delta_qm, delta_en, shift
        TYPE(environ_density) :: r, z, p, Ap, invsqrt
        !
        LOGICAL, POINTER :: lconjugate
        INTEGER, POINTER :: maxstep
        REAL(DP), POINTER :: tolvelect
        !
        CHARACTER(LEN=80) :: sub_name = 'generalized_gradient_sqrt'
        !
        !--------------------------------------------------------------------------------
        !
        lconjugate => gradient%lconjugate
        maxstep => gradient%maxstep
        tolvelect => gradient%tol
        !
        IF (verbose >= 1 .AND. ionode) WRITE (environ_unit, 9000)
        !
9000    FORMAT(/, 4('%'), ' COMPUTE ELECTROSTATIC POTENTIAL ', 43('%'))
        !
        !--------------------------------------------------------------------------------
        ! Check that fields have the same defintion domain
        !
        IF (.NOT. ASSOCIATED(charges%cell, dielectric%epsilon%cell)) &
            CALL env_errore(sub_name, 'Inconsistent cells of input fields', 1)
        !
        IF (.NOT. ASSOCIATED(charges%cell, potential%cell)) &
            CALL env_errore(sub_name, 'Inconsistent cells for charges and potential', 1)
        !
        cell => charges%cell
        !
        b => charges
        eps => dielectric%epsilon
        factsqrt => dielectric%factsqrt
        x => potential
        !
        CALL init_environ_density(cell, invsqrt)
        !
        invsqrt%of_r = 1.D0 / SQRT(eps%of_r)
        !
        CALL init_environ_density(cell, r)
        !
        CALL init_environ_density(cell, z)
        !
        CALL init_environ_density(cell, p)
        !
        CALL init_environ_density(cell, Ap)
        !
        !--------------------------------------------------------------------------------
        ! Starting guess from new input and previous solution(s)
        !
        IF (x%update) THEN
            r%of_r = b%of_r - factsqrt%of_r * x%of_r
            !
            !----------------------------------------------------------------------------
            ! Preconditioning step
            !
            z%of_r = r%of_r * invsqrt%of_r
            !
            CALL poisson_direct(core, z, z, electrolyte, semiconductor)
            !
            z%of_r = z%of_r * invsqrt%of_r
            !
            rzold = scalar_product_environ_density(r, z)
            !
            IF (ABS(rzold) < 1.D-30) &
                CALL env_errore(sub_name, 'Null step in gradient descent iteration', 1)
            !
            r%of_r = factsqrt%of_r * (x%of_r - z%of_r)
            delta_en = euclidean_norm_environ_density(r)
            delta_qm = quadratic_mean_environ_density(r)
            !
            IF (delta_en < 1.D-02) THEN
                !
                IF (verbose >= 1 .AND. ionode) WRITE (environ_unit, 9008) delta_en
                !
9008            FORMAT(' Sqrt-preconditioned input guess with residual norm = ', E14.6)
                x%of_r = z%of_r
            ELSE
                !
                IF (verbose >= 1 .AND. ionode) WRITE (environ_unit, 9001) delta_en
                !
9001            FORMAT(' Warning: Bad guess with residual norm = ', E14.6, &
                       ', reset to no guess')
                !
                x%update = .FALSE.
            END IF
            !
        END IF
        !
        IF (.NOT. x%update) THEN
            !
            x%update = .TRUE.
            x%of_r = 0.D0
            r%of_r = b%of_r
            rzold = 0.D0
            !
        END IF
        !
        !--------------------------------------------------------------------------------
        ! Start gradient descent
        !
        DO iter = 1, maxstep
            !
            IF (verbose >= 1 .AND. ionode) WRITE (environ_unit, 9002) iter
            !
9002        FORMAT(' Iteration # ', i10)
            !
            !----------------------------------------------------------------------------
            ! Apply preconditioner to new state
            !
            z%of_r = r%of_r * invsqrt%of_r
            !
            CALL poisson_direct(core, z, z, electrolyte, semiconductor)
            !
            z%of_r = z%of_r * invsqrt%of_r
            !
            rznew = scalar_product_environ_density(r, z)
            !
            IF (ABS(rznew) < 1.D-30) &
                CALL env_errore(sub_name, 'Null step in gradient descent iteration', 1)
            !
            !----------------------------------------------------------------------------
            ! Conjugate gradient or steepest descent input
            !
            IF (lconjugate .AND. ABS(rzold) > 1.D-30) THEN
                beta = rznew / rzold
            ELSE
                beta = 0.D0
            END IF
            !
            IF (verbose >= 3 .AND. ionode) &
                WRITE (environ_unit, *) 'rznew = ', rznew, ' rzold = ', &
                rzold, ' beta = ', beta
            !
            rzold = rznew
            p%of_r = z%of_r + beta * p%of_r
            !
            !----------------------------------------------------------------------------
            !
            Ap%of_r = factsqrt%of_r * z%of_r + r%of_r + beta * Ap%of_r
            ! apply operator to conjugate direction
            !
            !----------------------------------------------------------------------------
            ! Step downhill
            !
            pAp = scalar_product_environ_density(p, Ap)
            alpha = rzold / pAp
            !
            IF (verbose >= 1 .AND. ionode) &
                WRITE (environ_unit, *) ' pAp = ', pAp, ' rzold = ', &
                rzold, ' alpha = ', alpha
            !
            x%of_r = x%of_r + alpha * p%of_r
            r%of_r = r%of_r - alpha * Ap%of_r
            !
            !----------------------------------------------------------------------------
            ! If residual is small enough exit
            !
            delta_qm = quadratic_mean_environ_density(r)
            delta_en = euclidean_norm_environ_density(r)
            !
            IF (verbose >= 1 .AND. ionode) &
                WRITE (environ_unit, 9004) delta_qm, delta_en, tolvelect
            !
9004        FORMAT(' delta_qm = ', E14.6, ' delta_en = ', E14.6, ' tol = ', E14.6)
            !
            IF (delta_en < tolvelect .AND. iter > 0) THEN
                IF (verbose >= 1 .AND. ionode) WRITE (environ_unit, 9005)
9005            FORMAT(' Charges are converged, EXIT')
                !
                EXIT
                !
            ELSE IF (iter == maxstep) THEN
                IF (ionode) WRITE (program_unit, 9006)
9006            FORMAT(' Warning: Polarization charge not converged')
            END IF
            !
        END DO
        !
        !--------------------------------------------------------------------------------
        ! In PBC the potential need to have zero average
        !
        shift = 0.D0
        !
        IF (core%use_fft) THEN
            !
            IF (.NOT. (core%fft%use_internal_pbc_corr .OR. core%need_correction)) THEN
                shift = -integrate_environ_density(x) / cell%omega
            END IF
            !
        END IF
        !
        x%of_r = x%of_r + shift
        !
        IF (lstdout .AND. verbose >= 1) WRITE (program_unit, 9007) delta_en, iter
        !
9007    FORMAT('     Polarization accuracy =', 1PE8.1, ', # of iterations = ', i3)
        !
        CALL destroy_environ_density(r)
        !
        CALL destroy_environ_density(z)
        !
        CALL destroy_environ_density(p)
        !
        CALL destroy_environ_density(Ap)
        !
        CALL destroy_environ_density(invsqrt)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE generalized_gradient_sqrt
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE generalized_gradient_left(gradient, core, charges, dielectric, &
                                         potential, electrolyte, semiconductor)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(gradient_solver), TARGET, INTENT(IN) :: gradient
        TYPE(core_container), INTENT(IN) :: core
        TYPE(environ_density), TARGET, INTENT(IN) :: charges
        TYPE(environ_dielectric), TARGET, INTENT(IN) :: dielectric
        TYPE(environ_electrolyte), INTENT(IN), OPTIONAL :: electrolyte
        !
        TYPE(environ_density), TARGET, INTENT(INOUT) :: potential
        TYPE(environ_semiconductor), INTENT(INOUT), OPTIONAL :: semiconductor
        !
        TYPE(environ_cell), POINTER :: cell
        TYPE(environ_density), POINTER :: x, b, eps
        TYPE(environ_gradient), POINTER :: gradeps
        !
        INTEGER :: iter
        REAL(DP) :: rznew, rzold, alpha, beta, pAp, delta_en, delta_qm
        TYPE(environ_density) :: r, z, p, Ap
        TYPE(environ_gradient) :: g
        !
        LOGICAL, POINTER :: lconjugate
        INTEGER, POINTER :: maxstep
        REAL(DP), POINTER :: tolvelect
        !
        CHARACTER(LEN=80) :: sub_name = 'generalized_gradient_left'
        !
        !--------------------------------------------------------------------------------
        !
        lconjugate => gradient%lconjugate
        maxstep => gradient%maxstep
        tolvelect => gradient%tol
        !
        IF (verbose >= 1 .AND. ionode) WRITE (environ_unit, 9000)
        !
9000    FORMAT(/, 4('%'), ' COMPUTE ELECTROSTATIC POTENTIAL ', 43('%'))
        !
        !--------------------------------------------------------------------------------
        ! Check that fields have the same defintion domain
        !
        IF (.NOT. ASSOCIATED(charges%cell, dielectric%epsilon%cell)) &
            CALL env_errore(sub_name, 'Inconsistent cells of input fields', 1)
        !
        IF (.NOT. ASSOCIATED(charges%cell, potential%cell)) &
            CALL env_errore(sub_name, 'Inconsistent cells for charges and potential', 1)
        !
        cell => charges%cell
        !
        b => charges
        eps => dielectric%epsilon
        gradeps => dielectric%gradient
        x => potential
        !
        CALL init_environ_density(cell, r)
        !
        CALL init_environ_density(cell, z)
        !
        CALL init_environ_density(cell, p)
        !
        CALL init_environ_density(cell, Ap)
        !
        CALL init_environ_gradient(cell, g)
        !
        !--------------------------------------------------------------------------------
        ! Starting guess from new input and previous solution(s)
        !
!         IF (x%update) THEN

!             CALL gradient_fft(core%fft, x, g)
!             g%of_r = -g%of_r / fpi / e2
!             r%of_r(:) = b%of_r -  &
!                        & gradeps%of_r(1, :) * g%of_r(1, :) + &
!                        & gradeps%of_r(2, :) * g%of_r(2, :) + &
!                        & gradeps%of_r(3, :) * g%of_r(3, :)

!             ! ... Preconditioning step

!             z%of_r = r%of_r / eps%of_r
!             CALL poisson_direct(core, z, z)

!             rzold = scalar_product_environ_density(r, z)
!             IF (ABS(rzold) < 1.D-30) &
!                  & CALL env_errore(sub_name, 'Null step in gradient descent iteration', 1)

!             r%of_r = x%of_r - z%of_r
!             CALL gradient_fft(core%fft, r, g)
!             g%of_r = -g%of_r / fpi / e2
!             r%of_r(:) = gradeps%of_r(1, :) * g%of_r(1, :) + &
!                       & gradeps%of_r(2, :) * g%of_r(2, :) + &
!                       & gradeps%of_r(3, :) * g%of_r(3, :)
!             delta_qm = quadratic_mean_environ_density(r)
!             IF (delta_qm < 1.D-02) THEN
!                 IF (verbose >= 1 .AND. ionode) WRITE (environ_unit, 9008) delta_qm
! 9008            FORMAT(' Sqrt-preconditioned input guess with residual norm = ', E14.6)
!                 x%of_r = z%of_r
!             ELSE
!                 IF (verbose >= 1 .AND. ionode) WRITE (environ_unit, 9001) delta_qm
! 9001            FORMAT(' Warning: Bad guess with residual norm = ', E14.6, ', reset to no guess')
!                 x%update = .FALSE.
!             END IF

!         END IF

!         IF (.NOT. x%update) THEN

!             x%update = .TRUE.
        x%of_r = 0.D0
        r%of_r = b%of_r
        rzold = 0.D0
        !
!         END IF
        !
        !--------------------------------------------------------------------------------
        ! Start gradient descent
        !
        DO iter = 1, maxstep
            !
            IF (verbose >= 1 .AND. ionode) WRITE (environ_unit, 9002) iter
            !
9002        FORMAT(' Iteration # ', i10)
            !
            !----------------------------------------------------------------------------
            ! Apply preconditioner to new state
            !
            z%of_r = r%of_r / eps%of_r
            !
            CALL poisson_direct(core, z, z, electrolyte, semiconductor)
            !
            rznew = scalar_product_environ_density(r, z)
            !
            IF (ABS(rznew) < 1.D-30) &
                CALL env_errore(sub_name, 'Null step in gradient descent iteration', 1)
            !
            !----------------------------------------------------------------------------
            ! Conjugate gradient or steepest descent input
            !
            IF (lconjugate .AND. ABS(rzold) > 1.D-30) THEN
                beta = rznew / rzold
            ELSE
                beta = 0.D0
            END IF
            !
            IF (verbose >= 3 .AND. ionode) &
                WRITE (environ_unit, *) 'rznew = ', rznew, ' rzold = ', &
                rzold, ' beta = ', beta
            !
            rzold = rznew
            !
            p%of_r = z%of_r + beta * p%of_r
            !
            !----------------------------------------------------------------------------
            ! Apply operator to conjugate direction
            ! NOTE: the following steps should be extended to account for different cores
            !
            CALL gradient_fft(core%fft, z, g)
            !
            g%of_r = g%of_r / fpi / e2
            !
            Ap%of_r(:) = beta * Ap%of_r(:) - r%of_r(:) + &
                         gradeps%of_r(1, :) * g%of_r(1, :) + &
                         gradeps%of_r(2, :) * g%of_r(2, :) + &
                         gradeps%of_r(3, :) * g%of_r(3, :)
            !
            !----------------------------------------------------------------------------
            ! Step downhill
            !
            pAp = scalar_product_environ_density(p, Ap)
            alpha = rzold / pAp
            !
            IF (verbose >= 1 .AND. ionode) &
                WRITE (environ_unit, *) ' pAp = ', pAp, ' rzold = ', &
                rzold, ' alpha = ', alpha
            !
            x%of_r = x%of_r + alpha * p%of_r
            r%of_r = r%of_r - alpha * Ap%of_r
            !
            !----------------------------------------------------------------------------
            ! If residual is small enough exit
            !
            delta_qm = quadratic_mean_environ_density(r)
            delta_en = euclidean_norm_environ_density(r)
            !
            IF (verbose >= 1 .AND. ionode) &
                WRITE (environ_unit, 9004) delta_qm, delta_en, tolvelect
            !
9004        FORMAT(' delta_qm = ', E14.6, ' delta_en = ', E14.6, ' tol = ', E14.6)
            !
            IF (delta_en < tolvelect .AND. iter > 0) THEN
                IF (verbose >= 1 .AND. ionode) WRITE (environ_unit, 9005)
                !
9005            FORMAT(' Charges are converged, EXIT')
                !
                EXIT
                !
            ELSE IF (iter == maxstep) THEN
                !
                IF (ionode) WRITE (program_unit, 9006)
                !
9006            FORMAT(' Warning: Polarization charge not converged')
            END IF
            !
        END DO
        !
        IF (lstdout .AND. verbose >= 1) WRITE (program_unit, 9007) delta_en, iter
        !
9007    FORMAT('     Polarization accuracy =', 1PE8.1, ', # of iterations = ', i3)
        !
        CALL destroy_environ_gradient(g)
        !
        CALL destroy_environ_density(r)
        !
        CALL destroy_environ_density(z)
        !
        CALL destroy_environ_density(p)
        !
        CALL destroy_environ_density(Ap)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE generalized_gradient_left
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE solver_generalized
!----------------------------------------------------------------------------------------
