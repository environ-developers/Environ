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
! Authors: Francesco Nattino  (THEOS and NCCR-MARVEL, EPFL)
!          Oliviero Andreussi (Department of Physics, UNT)
!          Nicola Marzari     (THEOS and NCCR-MARVEL, EPFL)
!
!----------------------------------------------------------------------------------------
!>
!! This module contains the main drivers and routines to compute the
!! electrostatic potential that is the solution of a
!! linearized Poisson-Boltzmann equation, possibly in a dielectric medium:
!!
!! \f[
!!    ( \nabla \cdot \epsilon (r) \nabla + \gamma (r) * k^2 ) \phi = -4 \pi \rho
!! \f]
!!
!----------------------------------------------------------------------------------------
MODULE problem_linearized_pb
    !------------------------------------------------------------------------------------
    !
    USE modules_constants, ONLY: e2, fpi
    USE environ_types
    USE electrostatic_types
    USE environ_output
    USE problem_poisson, ONLY: poisson_direct
    !
    IMPLICIT NONE
    !
    PRIVATE
    !
    PUBLIC :: linearized_pb_gradient
    !
    INTERFACE linearized_pb_gradient
        MODULE PROCEDURE linearized_pb_gradient_charges, linearized_pb_gradient_density
    END INTERFACE linearized_pb_gradient
    !
    !------------------------------------------------------------------------------------
CONTAINS
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE linearized_pb_gradient_charges(solver, core, charges, potential, &
                                              screening)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(electrostatic_solver), INTENT(IN) :: solver
        TYPE(electrostatic_core), INTENT(IN) :: core
        TYPE(environ_charges), INTENT(IN) :: charges
        TYPE(environ_density), OPTIONAL, INTENT(IN) :: screening
        !
        TYPE(environ_density), INTENT(INOUT) :: potential
        !
        TYPE(environ_density) :: local_screening
        !
        CHARACTER(LEN=25) :: sub_name = 'linearized_pb_gradient'
        !
        !--------------------------------------------------------------------------------
        !
        CALL start_clock('calc_vlinpb')
        !
        IF (solver%use_gradient) THEN
            !
            CALL init_environ_density(potential%cell, local_screening)
            !
            IF (PRESENT(screening)) THEN
                local_screening%of_r = screening%of_r ! External screening
            ELSE
                !
                !------------------------------------------------------------------------
                ! Screening as in linearized pb problem
                !
                IF (charges%electrolyte%electrolyte_entropy == 'ions' &
                    .AND. charges%electrolyte%cionmax > 0.D0) THEN
                    !
                    local_screening%of_r = &
                        charges%electrolyte%k2 / e2 / fpi * &
                        charges%electrolyte%gamma%of_r / &
                        (1.D0 - SUM(charges%electrolyte%ioncctype(:)%cbulk) / &
                         charges%electrolyte%cionmax * &
                         (1.D0 - charges%electrolyte%gamma%of_r))
                    !
                ELSE
                    !
                    local_screening%of_r = charges%electrolyte%k2 / e2 / fpi * &
                                           charges%electrolyte%gamma%of_r
                    !
                END IF
                !
            END IF
            !
            IF (solver%auxiliary == 'none') THEN
                !
                IF (ASSOCIATED(charges%dielectric)) THEN
                    !
                    SELECT CASE (solver%gradient%preconditioner)
                    CASE ('sqrt')
                        !
                        CALL linearized_pb_gradient_sqrt(solver%gradient, core, &
                                                         charges%density, &
                                                         local_screening, potential, &
                                                         charges%dielectric)
                        !
                    CASE DEFAULT
                        CALL errore(sub_name, 'unexpected preconditioner keyword', 1)
                    END SELECT
                    !
                ELSE
                    !
                    CALL linearized_pb_gradient_sqrt(solver%gradient, core, &
                                                     charges%density, local_screening, &
                                                     potential)
                    !
                END IF
                !
            ELSE
                CALL errore(sub_name, 'Option not yet implemented', 1)
            END IF
            !
            CALL destroy_environ_density(local_screening)
            !
        ELSE
            CALL errore(sub_name, 'unexpected solver keyword', 1)
        END IF
        !
        CALL stop_clock('calc_vlinpb')
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE linearized_pb_gradient_charges
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE linearized_pb_gradient_density(solver, core, charges, electrolyte, &
                                              potential, dielectric, screening)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(electrostatic_solver), INTENT(IN) :: solver
        TYPE(electrostatic_core), INTENT(IN) :: core
        TYPE(environ_density), INTENT(IN) :: charges
        TYPE(environ_electrolyte), INTENT(IN) :: electrolyte
        TYPE(environ_dielectric), INTENT(IN), OPTIONAL :: dielectric
        TYPE(environ_density), INTENT(IN), OPTIONAL :: screening
        !
        TYPE(environ_density), INTENT(INOUT) :: potential
        !
        TYPE(environ_density) :: local_screening
        !
        CHARACTER(LEN=25) :: sub_name = 'linearized_pb_gradient'
        !
        !--------------------------------------------------------------------------------
        !
        CALL start_clock('calc_vlinpb')
        !
        IF (solver%use_gradient) THEN
            !
            CALL init_environ_density(potential%cell, local_screening)
            !
            IF (PRESENT(screening)) THEN
                local_screening%of_r = screening%of_r ! External screening
            ELSE
                !
                !------------------------------------------------------------------------
                ! Screening as in linearized pb problem
                !
                IF (electrolyte%electrolyte_entropy == 'ions' .AND. &
                    electrolyte%cionmax > 0.D0) THEN
                    !
                    local_screening%of_r = &
                        electrolyte%k2 / e2 / fpi * electrolyte%gamma%of_r / &
                        (1.D0 - SUM(electrolyte%ioncctype(:)%cbulk) / &
                         electrolyte%cionmax * (1.D0 - electrolyte%gamma%of_r))
                    !
                ELSE
                    !
                    local_screening%of_r = &
                        electrolyte%k2 / e2 / fpi * electrolyte%gamma%of_r
                    !
                END IF
                !
            END IF
            !
            IF (solver%auxiliary == 'none') THEN
                !
                IF (PRESENT(dielectric)) THEN
                    !
                    SELECT CASE (solver%gradient%preconditioner)
                    CASE ('sqrt')
                        !
                        CALL linearized_pb_gradient_sqrt(solver%gradient, core, &
                                                         charges, local_screening, &
                                                         potential, dielectric)
                        !
                    CASE DEFAULT
                        CALL errore(sub_name, 'unexpected preconditioner keyword', 1)
                    END SELECT
                    !
                ELSE
                    !
                    CALL linearized_pb_gradient_sqrt(solver%gradient, core, charges, &
                                                     local_screening, potential)
                    !
                END IF
                !
            ELSE
                CALL errore(sub_name, 'Option not yet implemented', 1)
            END IF
            !
            CALL destroy_environ_density(local_screening)
            !
        ELSE
            CALL errore(sub_name, 'unexpected solver keyword', 1)
        END IF
        !
        CALL stop_clock('calc_vlinpb')
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE linearized_pb_gradient_density
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE linearized_pb_gradient_sqrt(gradient, core, charges, screening, &
                                           potential, dielectric)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(gradient_solver), TARGET, INTENT(IN) :: gradient
        TYPE(electrostatic_core), INTENT(IN) :: core
        TYPE(environ_density), TARGET, INTENT(IN) :: charges
        TYPE(environ_density), TARGET, INTENT(IN) :: screening
        TYPE(environ_dielectric), OPTIONAL, TARGET, INTENT(IN) :: dielectric
        !
        TYPE(environ_density), TARGET, INTENT(INOUT) :: potential
        !
        TYPE(environ_cell), POINTER :: cell
        TYPE(environ_density), POINTER :: x, b, eps, factsqrt, scr
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
        CHARACTER(LEN=80) :: sub_name = 'linearized_pb_gradient_sqrt'
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
        IF (.NOT. ASSOCIATED(charges%cell, screening%cell)) &
            CALL errore(sub_name, 'Inconsistent cells of input fields', 1)
        !
        IF (.NOT. ASSOCIATED(charges%cell, potential%cell)) &
            CALL errore(sub_name, 'Inconsistent cells for charges and potential', 1)
        !
        IF (PRESENT(dielectric)) THEN
            !
            IF (.NOT. ASSOCIATED(charges%cell, dielectric%epsilon%cell)) &
                CALL errore(sub_name, 'Inconsistent cells of input fields', 1)
            !
        END IF
        !
        cell => charges%cell
        !
        b => charges
        x => potential
        scr => screening
        !
        IF (PRESENT(dielectric)) THEN
            eps => dielectric%epsilon
            factsqrt => dielectric%factsqrt
            !
            CALL init_environ_density(cell, invsqrt)
            !
            invsqrt%of_r = 1.D0 / SQRT(eps%of_r)
        END IF
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
            !
            IF (PRESENT(dielectric)) THEN
                r%of_r = b%of_r - (factsqrt%of_r + scr%of_r) * x%of_r
            ELSE
                r%of_r = b%of_r - scr%of_r * x%of_r
            END IF
            !
            !----------------------------------------------------------------------------
            ! Preconditioning step
            !
            IF (PRESENT(dielectric)) THEN
                z%of_r = r%of_r * invsqrt%of_r
                !
                CALL poisson_direct(core, z, z)
                !
                z%of_r = z%of_r * invsqrt%of_r
            ELSE
                z%of_r = r%of_r
                !
                CALL poisson_direct(core, z, z)
                !
            END IF
            !
            rzold = scalar_product_environ_density(r, z)
            !
            IF (ABS(rzold) < 1.D-30) &
                CALL errore(sub_name, 'Null step in gradient descent iteration', 1)
            !
            IF (PRESENT(dielectric)) THEN
                r%of_r = (factsqrt%of_r + scr%of_r) * (x%of_r - z%of_r)
            ELSE
                r%of_r = scr%of_r * (x%of_r - z%of_r)
            END IF
            !
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
9001            FORMAT(' Warning: bad guess with residual norm = ', &
                       E14.6, ', reset to no guess')
                !
                x%update = .FALSE.
            END IF
            !
        END IF
        !
        IF (.NOT. x%update) THEN
            x%update = .TRUE.
            x%of_r = 0.D0
            r%of_r = b%of_r
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
            IF (PRESENT(dielectric)) THEN
                z%of_r = r%of_r * invsqrt%of_r
                !
                CALL poisson_direct(core, z, z)
                !
                z%of_r = z%of_r * invsqrt%of_r
            ELSE
                z%of_r = r%of_r
                !
                CALL poisson_direct(core, z, z)
                !
            END IF
            !
            rznew = scalar_product_environ_density(r, z)
            !
            IF (ABS(rznew) < 1.D-30) &
                CALL errore(sub_name, 'Null step in gradient descent iteration', 1)
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
            !
            IF (PRESENT(dielectric)) THEN
                Ap%of_r = (factsqrt%of_r + scr%of_r) * z%of_r + r%of_r + beta * Ap%of_r
            ELSE
                Ap%of_r = scr%of_r * z%of_r + r%of_r + beta * Ap%of_r
            END IF
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
9007    FORMAT('     polarization accuracy =', 1PE8.1, ', # of iterations = ', i3)
        !
        CALL destroy_environ_density(r)
        !
        CALL destroy_environ_density(z)
        !
        CALL destroy_environ_density(p)
        !
        CALL destroy_environ_density(Ap)
        !
        IF (PRESENT(dielectric)) CALL destroy_environ_density(invsqrt)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE linearized_pb_gradient_sqrt
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE problem_linearized_pb
!----------------------------------------------------------------------------------------
