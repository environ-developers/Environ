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
! Authors: Francesco Nattino  (THEOS and NCCR-MARVEL, EPFL)
!          Oliviero Andreussi (Department of Physics, UNT)
!          Nicola Marzari     (THEOS and NCCR-MARVEL, EPFL)
!          Edan Bainglass     (Department of Physics, UNT)
!
!----------------------------------------------------------------------------------------
!>
!!
!----------------------------------------------------------------------------------------
MODULE class_solver_gradient
    !------------------------------------------------------------------------------------
    !
    USE class_io, ONLY: io
    !
    USE environ_param, ONLY: DP, e2, fpi
    !
    USE class_cell
    USE class_density
    USE class_gradient
    !
    USE class_core_container_electrostatics
    USE class_core_fft_electrostatics
    !
    USE class_solver_iterative
    !
    USE class_charges
    USE class_dielectric
    USE class_electrolyte
    USE class_semiconductor
    !
    !------------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    PRIVATE
    !
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    TYPE, EXTENDS(solver_iterative), PUBLIC :: solver_gradient
        !--------------------------------------------------------------------------------
        !
        LOGICAL :: lconjugate
        CHARACTER(LEN=80) :: step_type
        REAL(DP) :: step
        CHARACTER(LEN=80) :: preconditioner
        CHARACTER(LEN=80) :: screening_type
        REAL(DP) :: screening
        !
        !--------------------------------------------------------------------------------
    CONTAINS
        !--------------------------------------------------------------------------------
        !
        PROCEDURE :: init => init_solver_gradient
        !
        PROCEDURE, PRIVATE :: generalized_gradient_charges, generalized_gradient_density
        !
        GENERIC :: generalized => &
            generalized_gradient_charges, generalized_gradient_density
        !
        PROCEDURE, PRIVATE :: &
            linearized_pb_gradient_charges, linearized_pb_gradient_density
        !
        GENERIC :: linearized_pb => &
            linearized_pb_gradient_charges, linearized_pb_gradient_density
        !
        PROCEDURE, PRIVATE :: generalized_none => generalized_gradient_none
        PROCEDURE, PRIVATE :: generalized_sqrt => generalized_gradient_sqrt
        PROCEDURE, PRIVATE :: generalized_left => generalized_gradient_left
        !
        PROCEDURE, PRIVATE :: linearized_pb_sqrt => linearized_pb_gradient_sqrt
        !
        !--------------------------------------------------------------------------------
    END TYPE solver_gradient
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
CONTAINS
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !                                   ADMIN METHODS
    !
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_solver_gradient(this, lconjugate, step_type, step, preconditioner, &
                                    screening_type, screening, cores, maxiter, tol, &
                                    auxiliary)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(container_electrostatics), INTENT(IN) :: cores
        LOGICAL, INTENT(IN) :: lconjugate
        INTEGER, INTENT(IN) :: maxiter
        REAL(DP), INTENT(IN) :: tol, step, screening
        CHARACTER(LEN=80), INTENT(IN) :: step_type, preconditioner, screening_type
        CHARACTER(LEN=80), INTENT(IN), OPTIONAL :: auxiliary
        !
        CLASS(solver_gradient), INTENT(INOUT) :: this
        !
        !--------------------------------------------------------------------------------
        !
        CALL this%init_iterative(cores, maxiter, tol, auxiliary)
        !
        IF (lconjugate) THEN
            this%solver_type = 'cg'
        ELSE
            this%solver_type = 'sd'
        END IF
        !
        this%lconjugate = lconjugate
        this%step_type = step_type
        this%step = step
        this%preconditioner = preconditioner
        this%screening_type = screening_type
        this%screening = screening
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_solver_gradient
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !                                   SOLVER METHODS
    !
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE generalized_gradient_charges(this, charges, potential)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(solver_gradient), INTENT(IN) :: this
        !
        TYPE(environ_charges), INTENT(INOUT) :: charges
        TYPE(environ_density), INTENT(INOUT) :: potential
        !
        CHARACTER(LEN=80) :: sub_name = 'generalized_gradient_charges'
        !
        !--------------------------------------------------------------------------------
        !
        CALL env_start_clock(sub_name)
        !
        potential%of_r = 0.D0
        !
        IF (this%auxiliary == 'none') THEN
            !
            SELECT CASE (this%preconditioner) ! #TODO none and left need work
                !
            CASE ('none')
                !
                CALL this%generalized_none(charges%density, charges%dielectric, &
                                           potential, charges%electrolyte, &
                                           charges%semiconductor)
                !
            CASE ('sqrt')
                !
                CALL this%generalized_sqrt(charges%density, charges%dielectric, &
                                           potential, charges%electrolyte, &
                                           charges%semiconductor)
                !
            CASE ('left')
                !
                CALL this%generalized_left(charges%density, charges%dielectric, &
                                           potential, charges%electrolyte, &
                                           charges%semiconductor)
                !
            CASE DEFAULT
                CALL io%error(sub_name, 'Unexpected preconditioner keyword', 1)
                !
            END SELECT
            !
        ELSE
            !
            CALL io%error(sub_name, 'Option not yet implemented', 1)
            !
            ! CALL generalized_gradient_rhoaux(charges, dielectric, potential) #TODO future-work
            !
        END IF
        !
        CALL env_stop_clock(sub_name)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE generalized_gradient_charges
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE generalized_gradient_density(this, charges, dielectric, potential, &
                                            electrolyte, semiconductor)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(solver_gradient), INTENT(IN) :: this
        TYPE(environ_density), INTENT(IN) :: charges
        TYPE(environ_dielectric), INTENT(IN) :: dielectric
        TYPE(environ_electrolyte), INTENT(IN), OPTIONAL :: electrolyte
        !
        TYPE(environ_density), INTENT(INOUT) :: potential
        TYPE(environ_semiconductor), INTENT(INOUT), OPTIONAL :: semiconductor
        !
        CHARACTER(LEN=80) :: sub_name = 'generalized_gradient_density'
        !
        !--------------------------------------------------------------------------------
        !
        CALL env_start_clock(sub_name)
        !
        potential%of_r = 0.D0
        !
        IF (this%auxiliary == 'none') THEN
            !
            SELECT CASE (this%preconditioner) ! #TODO none and left need work
                !
            CASE ('none')
                CALL this%generalized_none(charges, dielectric, potential, electrolyte)
                !
            CASE ('sqrt')
                CALL this%generalized_sqrt(charges, dielectric, potential, electrolyte)
                !
            CASE ('left')
                CALL this%generalized_left(charges, dielectric, potential, electrolyte)
                !
            CASE DEFAULT
                CALL io%error(sub_name, 'Unexpected preconditioner keyword', 1)
                !
            END SELECT
            !
        ELSE
            !
            CALL io%error(sub_name, 'Option not yet implemented', 1)
            !
            ! CALL generalized_gradient_rhoaux(charges, dielectric, potential) #TODO future-work
            !
        END IF
        !
        CALL env_stop_clock(sub_name)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE generalized_gradient_density
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE linearized_pb_gradient_charges(this, charges, potential, screening)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(solver_gradient), INTENT(IN) :: this
        TYPE(environ_charges), INTENT(IN) :: charges
        TYPE(environ_density), OPTIONAL, INTENT(IN) :: screening
        !
        TYPE(environ_density), INTENT(INOUT) :: potential
        !
        TYPE(environ_density) :: local_screening
        !
        CHARACTER(LEN=80) :: sub_name = 'linearized_pb_gradient_charges'
        !
        !--------------------------------------------------------------------------------
        !
        CALL env_start_clock(sub_name)
        !
        CALL local_screening%init(potential%cell)
        !
        IF (PRESENT(screening)) THEN
            local_screening%of_r = screening%of_r ! external screening
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
        IF (this%auxiliary == 'none') THEN
            !
            IF (ASSOCIATED(charges%dielectric)) THEN
                !
                SELECT CASE (this%preconditioner)
                    !
                CASE ('sqrt')
                    !
                    CALL this%linearized_pb_sqrt(charges%density, local_screening, &
                                                 potential, charges%dielectric)
                    !
                CASE DEFAULT
                    CALL io%error(sub_name, 'Unexpected preconditioner keyword', 1)
                    !
                END SELECT
                !
            ELSE
                !
                CALL this%linearized_pb_sqrt(charges%density, local_screening, &
                                             potential)
                !
            END IF
            !
        ELSE
            CALL io%error(sub_name, 'Option not yet implemented', 1)
        END IF
        !
        CALL local_screening%destroy()
        !
        CALL env_stop_clock(sub_name)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE linearized_pb_gradient_charges
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE linearized_pb_gradient_density(this, charges, electrolyte, potential, &
                                              dielectric, screening)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(solver_gradient), INTENT(IN) :: this
        TYPE(environ_density), INTENT(IN) :: charges
        TYPE(environ_electrolyte), INTENT(IN) :: electrolyte
        TYPE(environ_dielectric), INTENT(IN), OPTIONAL :: dielectric
        TYPE(environ_density), INTENT(IN), OPTIONAL :: screening
        !
        TYPE(environ_density), INTENT(INOUT) :: potential
        !
        TYPE(environ_density) :: local_screening
        !
        CHARACTER(LEN=80) :: sub_name = 'linearized_pb_gradient_density'
        !
        !--------------------------------------------------------------------------------
        !
        CALL env_start_clock(sub_name)
        !
        CALL local_screening%init(potential%cell)
        !
        IF (PRESENT(screening)) THEN
            local_screening%of_r = screening%of_r ! external screening
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
                local_screening%of_r = electrolyte%k2 / e2 / fpi * electrolyte%gamma%of_r
            END IF
            !
        END IF
        !
        IF (this%auxiliary == 'none') THEN
            !
            IF (PRESENT(dielectric)) THEN
                !
                SELECT CASE (this%preconditioner)
                    !
                CASE ('sqrt')
                    !
                    CALL this%linearized_pb_sqrt(charges, local_screening, potential, &
                                                 dielectric)
                    !
                CASE DEFAULT
                    CALL io%error(sub_name, 'Unexpected preconditioner keyword', 1)
                    !
                END SELECT
                !
            ELSE
                CALL this%linearized_pb_sqrt(charges, local_screening, potential)
            END IF
            !
        ELSE
            CALL io%error(sub_name, 'Option not yet implemented', 1)
        END IF
        !
        CALL local_screening%destroy()
        !
        CALL env_stop_clock(sub_name)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE linearized_pb_gradient_density
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !                              PRIVATE SOLVER METHODS
    !
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE generalized_gradient_none(this, charges, dielectric, potential, &
                                         electrolyte, semiconductor)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(solver_gradient), TARGET, INTENT(IN) :: this
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
        INTEGER, POINTER :: maxiter
        REAL(DP), POINTER :: tolvelect
        !
        CHARACTER(LEN=80) :: sub_name = 'generalized_gradient_none'
        !
        !--------------------------------------------------------------------------------
        !
        lconjugate => this%lconjugate
        maxiter => this%maxiter
        tolvelect => this%tol
        !
        IF (io%verbosity >= 1 .AND. io%lnode) WRITE (io%debug_unit, 1000)
        !
        !--------------------------------------------------------------------------------
        ! Check that fields have the same defintion domain
        !
        IF (.NOT. ASSOCIATED(charges%cell, dielectric%epsilon%cell)) &
            CALL io%error(sub_name, 'Inconsistent cells of input fields', 1)
        !
        IF (.NOT. ASSOCIATED(charges%cell, potential%cell)) &
            CALL io%error(sub_name, 'Inconsistent cells for charges and potential', 1)
        !
        cell => charges%cell
        !
        b => charges
        eps => dielectric%epsilon
        gradeps => dielectric%gradient
        x => potential
        !
        CALL r%init(cell)
        !
        CALL z%init(cell)
        !
        CALL p%init(cell)
        !
        CALL Ap%init(cell)
        !
        CALL g%init(cell)
        !
        CALL l%init(cell)
        !
        !--------------------------------------------------------------------------------
        ! Starting guess from new input and previous solution(s)
        !
        IF (x%lupdate) THEN
            x%lupdate = .FALSE.
        END IF
        !
        IF (.NOT. x%lupdate) THEN
            x%lupdate = .TRUE.
            x%of_r = 0.D0
            r = b
            rzold = 0.D0
        END IF
        !
        !--------------------------------------------------------------------------------
        ! Write output table column headers
        !
        IF (io%lnode) THEN
            !
            IF (io%verbosity >= 3) THEN
                WRITE (io%debug_unit, 1001)
            ELSE IF (io%verbosity >= 1) THEN
                WRITE (io%debug_unit, 1002)
            END IF
            !
        END IF
        !
        !--------------------------------------------------------------------------------
        ! Start gradient descent
        !
        DO iter = 1, maxiter
            !
            !----------------------------------------------------------------------------
            ! Apply preconditioner to new state
            !
            z%of_r = r%of_r ! no preconditioner
            !
            rznew = r%scalar_product(z)
            !
            IF (ABS(rznew) < 1.D-30) &
                CALL io%error(sub_name, 'Null step in gradient descent iteration', 1)
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
            rzold = rznew
            p%of_r = z%of_r + beta * p%of_r
            !
            !----------------------------------------------------------------------------
            ! Apply operator to conjugate direction
            ! NOTE: the following steps should be extended to account for different cores
            !
            CALL this%cores%gradient(p, g)
            !
            CALL this%cores%laplacian(p, l)
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
            pAp = p%scalar_product(Ap)
            alpha = rzold / pAp
            !
            x%of_r = x%of_r + alpha * p%of_r
            r%of_r = r%of_r - alpha * Ap%of_r
            !
            delta_en = r%euclidean_norm()
            delta_qm = r%quadratic_mean()
            !
            !----------------------------------------------------------------------------
            ! Print iteration results
            !
            IF (io%lnode) THEN
                !
                IF (io%verbosity >= 3) THEN
                    !
                    WRITE (io%debug_unit, 1003) &
                        iter, alpha, beta, rznew, rzold, pAp, delta_qm, delta_en, &
                        tolvelect
                    !
                ELSE IF (io%verbosity >= 1) THEN
                    !
                    WRITE (io%debug_unit, 1004) &
                        iter, alpha, beta, delta_qm, delta_en, tolvelect
                    !
                END IF
                !
            END IF
            !
            !----------------------------------------------------------------------------
            ! If residual is small enough exit
            !
            IF (delta_en < tolvelect .AND. iter > 0) THEN
                !
                IF (io%verbosity >= 1 .AND. io%lnode) WRITE (io%debug_unit, 1005)
                !
                EXIT
                !
            ELSE IF (iter == maxiter) THEN
                IF (io%lnode) WRITE (io%unit, 1006)
            END IF
            !
        END DO
        !
        IF (io%lstdout .AND. io%verbosity >= 1) WRITE (io%unit, 1007) delta_en, iter
        !
        CALL l%destroy()
        !
        CALL g%destroy()
        !
        CALL r%destroy()
        !
        CALL z%destroy()
        !
        CALL p%destroy()
        !
        CALL Ap%destroy()
        !
        !--------------------------------------------------------------------------------
        !
1000    FORMAT(/, 4('%'), ' COMPUTE ELECTROSTATIC POTENTIAL ', 43('%'),/)
        !
1001    FORMAT('   i | ' &
               '     alpha     |      beta      |      rznew     |      rzold     | ' &
               '      pAp      |    delta_qm    |    delta_en    |       tol', /, &
               1X, 139('-'))
        !
1002    FORMAT('   i | ' &
               '     alpha     |      beta      |       pAp      |    delta_qm    | ' &
               '   delta_en    |       tol', /, 1X, 105('-'))
        !
1003    FORMAT(1X, I3, 8(' | ', E14.6))
        !
1004    FORMAT(1X, I3, 6(' | ', E14.6))
        !
1005    FORMAT(/, ' Charges are converged, EXIT')
        !
1006    FORMAT(' Warning: Polarization charge not converged',/)
        !
1007    FORMAT('     Polarization accuracy =', 1PE8.1, ', # of iterations = ', i3)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE generalized_gradient_none
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE generalized_gradient_sqrt(this, charges, dielectric, potential, &
                                         electrolyte, semiconductor)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(solver_gradient), TARGET, INTENT(IN) :: this
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
        INTEGER, POINTER :: maxiter
        REAL(DP), POINTER :: tolvelect
        !
        CHARACTER(LEN=80) :: sub_name = 'generalized_gradient_sqrt'
        !
        !--------------------------------------------------------------------------------
        !
        lconjugate => this%lconjugate
        maxiter => this%maxiter
        tolvelect => this%tol
        !
        IF (io%verbosity >= 1 .AND. io%lnode) WRITE (io%debug_unit, 1100)
        !
        !--------------------------------------------------------------------------------
        ! Check that fields have the same defintion domain
        !
        IF (.NOT. ASSOCIATED(charges%cell, dielectric%epsilon%cell)) &
            CALL io%error(sub_name, 'Inconsistent cells of input fields', 1)
        !
        IF (.NOT. ASSOCIATED(charges%cell, potential%cell)) &
            CALL io%error(sub_name, 'Inconsistent cells for charges and potential', 1)
        !
        cell => charges%cell
        !
        b => charges
        eps => dielectric%epsilon
        factsqrt => dielectric%factsqrt
        x => potential
        !
        CALL invsqrt%init(cell)
        !
        invsqrt%of_r = 1.D0 / SQRT(eps%of_r)
        !
        CALL r%init(cell)
        !
        CALL z%init(cell)
        !
        CALL p%init(cell)
        !
        CALL Ap%init(cell)
        !
        !--------------------------------------------------------------------------------
        ! Starting guess from new input and previous solution(s)
        !
        IF (x%lupdate) THEN
            r%of_r = b%of_r - factsqrt%of_r * x%of_r
            !
            !----------------------------------------------------------------------------
            ! Preconditioning step
            !
            z%of_r = r%of_r * invsqrt%of_r
            !
            CALL this%poisson(z, z, electrolyte, semiconductor)
            !
            z%of_r = z%of_r * invsqrt%of_r
            !
            rzold = r%scalar_product(z)
            !
            IF (ABS(rzold) < 1.D-30) &
                CALL io%error(sub_name, 'Null step in gradient descent iteration', 1)
            !
            r%of_r = factsqrt%of_r * (x%of_r - z%of_r)
            delta_en = r%euclidean_norm()
            delta_qm = r%quadratic_mean()
            !
            IF (delta_en < 1.D-02) THEN
                !
                IF (io%verbosity >= 1 .AND. io%lnode) WRITE (io%debug_unit, 1101) delta_en
                !
                x%of_r = z%of_r
            ELSE
                !
                IF (io%verbosity >= 1 .AND. io%lnode) WRITE (io%debug_unit, 1102) delta_en
                !
                x%lupdate = .FALSE.
            END IF
            !
        END IF
        !
        IF (.NOT. x%lupdate) THEN
            !
            x%lupdate = .TRUE.
            x%of_r = 0.D0
            r%of_r = b%of_r
            rzold = 0.D0
            !
        END IF
        !
        !--------------------------------------------------------------------------------
        ! Write output table column headers
        !
        IF (io%lnode) THEN
            !
            IF (io%verbosity >= 3) THEN
                WRITE (io%debug_unit, 1103)
            ELSE IF (io%verbosity >= 1) THEN
                WRITE (io%debug_unit, 1104)
            END IF
            !
        END IF
        !
        !--------------------------------------------------------------------------------
        ! Start gradient descent
        !
        DO iter = 1, maxiter
            !
            !----------------------------------------------------------------------------
            ! Apply preconditioner to new state
            !
            z%of_r = r%of_r * invsqrt%of_r
            !
            CALL this%poisson(z, z, electrolyte, semiconductor)
            !
            z%of_r = z%of_r * invsqrt%of_r
            !
            rznew = r%scalar_product(z)
            !
            IF (ABS(rznew) < 1.D-30) &
                CALL io%error(sub_name, 'Null step in gradient descent iteration', 1)
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
            pAp = p%scalar_product(Ap)
            alpha = rzold / pAp
            !
            x%of_r = x%of_r + alpha * p%of_r
            r%of_r = r%of_r - alpha * Ap%of_r
            !
            delta_en = r%euclidean_norm()
            delta_qm = r%quadratic_mean()
            !
            !----------------------------------------------------------------------------
            ! Print iteration results
            !
            IF (io%lnode) THEN
                !
                IF (io%verbosity >= 3) THEN
                    !
                    WRITE (io%debug_unit, 1105) &
                        iter, alpha, beta, rznew, rzold, pAp, delta_qm, delta_en, &
                        tolvelect
                    !
                ELSE IF (io%verbosity >= 1) THEN
                    !
                    WRITE (io%debug_unit, 1106) &
                        iter, alpha, pAp, rzold, delta_qm, delta_en, tolvelect
                    !
                END IF
                !
            END IF
            !
            !----------------------------------------------------------------------------
            ! If residual is small enough exit
            !
            IF (delta_en < tolvelect .AND. iter > 0) THEN
                !
                IF (io%verbosity >= 1 .AND. io%lnode) WRITE (io%debug_unit, 1107)
                !
                EXIT
                !
            ELSE IF (iter == maxiter) THEN
                IF (io%lnode) WRITE (io%unit, 1108)
            END IF
            !
        END DO
        !
        !--------------------------------------------------------------------------------
        ! In PBC the potential need to have zero average
        !
        shift = 0.D0
        !
        SELECT TYPE (core => this%cores%core)
            !
        TYPE IS (core_fft_electrostatics)
            !
            IF (.NOT. (core%use_internal_pbc_corr .OR. &
                       ASSOCIATED(this%cores%correction))) &
                shift = -x%integrate() / cell%omega
            !
        END SELECT
        !
        x%of_r = x%of_r + shift
        !
        IF (io%lstdout .AND. io%verbosity >= 1) WRITE (io%unit, 1109) delta_en, iter
        !
        CALL r%destroy()
        !
        CALL z%destroy()
        !
        CALL p%destroy()
        !
        CALL Ap%destroy()
        !
        CALL invsqrt%destroy()
        !
        !--------------------------------------------------------------------------------
        !
1100    FORMAT(/, 4('%'), ' COMPUTE ELECTROSTATIC POTENTIAL ', 43('%'),/)
        !
1101    FORMAT(' Sqrt-preconditioned input guess with residual norm = ', E14.6,/)
        !
1102    FORMAT(' Warning: Bad guess with residual norm = ', E14.6, &
               ', reset to no guess',/)
        !
1103    FORMAT('   i | ' &
               '     alpha     |      beta      |      rznew     |      rzold     | ' &
               '      pAp      |    delta_qm    |    delta_en    |       tol', /, &
               1X, 139('-'))
        !
1104    FORMAT('   i | ' &
               '     alpha     |       pAp      |      rzold     |    delta_qm    | ' &
               '   delta_en    |       tol', /, 1X, 105('-'))
        !
1105    FORMAT(1X, I3, 8(' | ', E14.6))
        !
1106    FORMAT(1X, I3, 6(' | ', E14.6))
        !
1107    FORMAT(/, ' Charges are converged, EXIT')
        !
1108    FORMAT(' Warning: Polarization charge not converged',/)
        !
1109    FORMAT('     Polarization accuracy =', 1PE8.1, ', # of iterations = ', i3)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE generalized_gradient_sqrt
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE generalized_gradient_left(this, charges, dielectric, potential, &
                                         electrolyte, semiconductor)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(solver_gradient), TARGET, INTENT(IN) :: this
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
        INTEGER, POINTER :: maxiter
        REAL(DP), POINTER :: tolvelect
        !
        CHARACTER(LEN=80) :: sub_name = 'generalized_gradient_left'
        !
        !--------------------------------------------------------------------------------
        !
        lconjugate => this%lconjugate
        maxiter => this%maxiter
        tolvelect => this%tol
        !
        IF (io%verbosity >= 1 .AND. io%lnode) WRITE (io%debug_unit, 1200)
        !
        !--------------------------------------------------------------------------------
        ! Check that fields have the same defintion domain
        !
        IF (.NOT. ASSOCIATED(charges%cell, dielectric%epsilon%cell)) &
            CALL io%error(sub_name, 'Inconsistent cells of input fields', 1)
        !
        IF (.NOT. ASSOCIATED(charges%cell, potential%cell)) &
            CALL io%error(sub_name, 'Inconsistent cells for charges and potential', 1)
        !
        cell => charges%cell
        !
        b => charges
        eps => dielectric%epsilon
        gradeps => dielectric%gradient
        x => potential
        !
        CALL r%init(cell)
        !
        CALL z%init(cell)
        !
        CALL p%init(cell)
        !
        CALL Ap%init(cell)
        !
        CALL g%init(cell)
        !
        !--------------------------------------------------------------------------------
        ! Starting guess from new input and previous solution(s)
        !
        x%of_r = 0.D0
        r%of_r = b%of_r
        rzold = 0.D0
        !
        !--------------------------------------------------------------------------------
        ! Write output table column headers
        !
        IF (io%lnode) THEN
            !
            IF (io%verbosity >= 3) THEN
                WRITE (io%debug_unit, 1201)
            ELSE IF (io%verbosity >= 1) THEN
                WRITE (io%debug_unit, 1202)
            END IF
            !
        END IF
        !
        !--------------------------------------------------------------------------------
        ! Start gradient descent
        !
        DO iter = 1, maxiter
            !
            !----------------------------------------------------------------------------
            ! Apply preconditioner to new state
            !
            z%of_r = r%of_r / eps%of_r
            !
            CALL this%poisson(z, z, electrolyte, semiconductor)
            !
            rznew = r%scalar_product(z)
            !
            IF (ABS(rznew) < 1.D-30) &
                CALL io%error(sub_name, 'Null step in gradient descent iteration', 1)
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
            rzold = rznew
            p%of_r = z%of_r + beta * p%of_r
            !
            !----------------------------------------------------------------------------
            ! Apply operator to conjugate direction
            ! NOTE: the following steps should be extended to account for different cores
            !
            CALL this%cores%gradient(z, g)
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
            pAp = p%scalar_product(Ap)
            alpha = rzold / pAp
            !
            x%of_r = x%of_r + alpha * p%of_r
            r%of_r = r%of_r - alpha * Ap%of_r
            !
            delta_en = r%euclidean_norm()
            delta_qm = r%quadratic_mean()
            !
            !----------------------------------------------------------------------------
            ! Print iteration results
            !
            IF (io%lnode) THEN
                !
                IF (io%verbosity >= 3) THEN
                    !
                    WRITE (io%debug_unit, 1203) &
                        iter, alpha, beta, rznew, rzold, pAp, delta_qm, delta_en, &
                        tolvelect
                    !
                ELSE IF (io%verbosity >= 1) THEN
                    !
                    WRITE (io%debug_unit, 1204) &
                        iter, alpha, pAp, rzold, delta_qm, delta_en, tolvelect
                    !
                END IF
                !
            END IF
            !
            !----------------------------------------------------------------------------
            ! If residual is small enough exit
            !
            IF (delta_en < tolvelect .AND. iter > 0) THEN
                !
                IF (io%verbosity >= 1 .AND. io%lnode) WRITE (io%debug_unit, 1205)
                !
                EXIT
                !
            ELSE IF (iter == maxiter) THEN
                IF (io%lnode) WRITE (io%unit, 1206)
            END IF
            !
        END DO
        !
        IF (io%lstdout .AND. io%verbosity >= 1) WRITE (io%unit, 1207) delta_en, iter
        !
        CALL g%destroy()
        !
        CALL r%destroy()
        !
        CALL z%destroy()
        !
        CALL p%destroy()
        !
        CALL Ap%destroy()
        !
        !--------------------------------------------------------------------------------
        !
1200    FORMAT(/, 4('%'), ' COMPUTE ELECTROSTATIC POTENTIAL ', 43('%'),/)
        !
1201    FORMAT('   i | ' &
               '     alpha     |      beta      |      rznew     |      rzold     | ' &
               '      pAp      |    delta_qm    |    delta_en    |       tol', /, &
               1X, 139('-'))
        !
1202    FORMAT('   i | ' &
               '     alpha     |       pAp      |      rzold     |    delta_qm    | ' &
               '   delta_en    |       tol', /, 1X, 105('-'))
        !
1203    FORMAT(1X, I3, 8(' | ', E14.6))
        !
1204    FORMAT(1X, I3, 6(' | ', E14.6))
        !
1205    FORMAT(/, ' Charges are converged, EXIT')
        !
1206    FORMAT(' Warning: Polarization charge not converged',/)
        !
1207    FORMAT('     Polarization accuracy =', 1PE8.1, ', # of iterations = ', i3)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE generalized_gradient_left
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE linearized_pb_gradient_sqrt(this, charges, screening, potential, &
                                           dielectric)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(solver_gradient), TARGET, INTENT(IN) :: this
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
        INTEGER, POINTER :: maxiter
        REAL(DP), POINTER :: tolvelect
        !
        CHARACTER(LEN=80) :: sub_name = 'linearized_pb_gradient_sqrt'
        !
        !--------------------------------------------------------------------------------
        !
        lconjugate => this%lconjugate
        maxiter => this%maxiter
        tolvelect => this%tol
        !
        IF (io%verbosity >= 1 .AND. io%lnode) WRITE (io%debug_unit, 1300)
        !
        !--------------------------------------------------------------------------------
        ! Check that fields have the same defintion domain
        !
        IF (.NOT. ASSOCIATED(charges%cell, screening%cell)) &
            CALL io%error(sub_name, 'Inconsistent cells of input fields', 1)
        !
        IF (.NOT. ASSOCIATED(charges%cell, potential%cell)) &
            CALL io%error(sub_name, 'Inconsistent cells for charges and potential', 1)
        !
        IF (PRESENT(dielectric)) THEN
            !
            IF (.NOT. ASSOCIATED(charges%cell, dielectric%epsilon%cell)) &
                CALL io%error(sub_name, 'Inconsistent cells of input fields', 1)
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
            CALL invsqrt%init(cell)
            !
            invsqrt%of_r = 1.D0 / SQRT(eps%of_r)
        END IF
        !
        CALL r%init(cell)
        !
        CALL z%init(cell)
        !
        CALL p%init(cell)
        !
        CALL Ap%init(cell)
        !
        !--------------------------------------------------------------------------------
        ! Starting guess from new input and previous solution(s)
        !
        IF (x%lupdate) THEN
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
                CALL this%poisson(z, z)
                !
                z%of_r = z%of_r * invsqrt%of_r
            ELSE
                z%of_r = r%of_r
                !
                CALL this%poisson(z, z)
                !
            END IF
            !
            rzold = r%scalar_product(z)
            !
            IF (ABS(rzold) < 1.D-30) &
                CALL io%error(sub_name, 'Null step in gradient descent iteration', 1)
            !
            IF (PRESENT(dielectric)) THEN
                r%of_r = (factsqrt%of_r + scr%of_r) * (x%of_r - z%of_r)
            ELSE
                r%of_r = scr%of_r * (x%of_r - z%of_r)
            END IF
            !
            delta_en = r%euclidean_norm()
            delta_qm = r%quadratic_mean()
            !
            IF (delta_en < 1.D-02) THEN
                !
                IF (io%verbosity >= 1 .AND. io%lnode) WRITE (io%debug_unit, 1301) delta_en
                !
                x%of_r = z%of_r
            ELSE
                !
                IF (io%verbosity >= 1 .AND. io%lnode) WRITE (io%debug_unit, 1302) delta_en
                !
                x%lupdate = .FALSE.
            END IF
            !
        END IF
        !
        IF (.NOT. x%lupdate) THEN
            x%lupdate = .TRUE.
            x%of_r = 0.D0
            r%of_r = b%of_r
            rzold = 0.D0
        END IF
        !
        !--------------------------------------------------------------------------------
        ! Write output table column headers
        !
        IF (io%lnode) THEN
            !
            IF (io%verbosity >= 3) THEN
                WRITE (io%debug_unit, 1303)
            ELSE IF (io%verbosity >= 1) THEN
                WRITE (io%debug_unit, 1304)
            END IF
            !
        END IF
        !
        !--------------------------------------------------------------------------------
        ! Start gradient descent
        !
        DO iter = 1, maxiter
            !
            !----------------------------------------------------------------------------
            ! Apply preconditioner to new state
            !
            IF (PRESENT(dielectric)) THEN
                z%of_r = r%of_r * invsqrt%of_r
                !
                CALL this%poisson(z, z)
                !
                z%of_r = z%of_r * invsqrt%of_r
            ELSE
                z%of_r = r%of_r
                !
                CALL this%poisson(z, z)
                !
            END IF
            !
            rznew = r%scalar_product(z)
            !
            IF (ABS(rznew) < 1.D-30) &
                CALL io%error(sub_name, 'Null step in gradient descent iteration', 1)
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
            rzold = rznew
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
            pAp = p%scalar_product(Ap)
            alpha = rzold / pAp
            !
            x%of_r = x%of_r + alpha * p%of_r
            r%of_r = r%of_r - alpha * Ap%of_r
            !
            delta_en = r%euclidean_norm()
            delta_qm = r%quadratic_mean()
            !
            !----------------------------------------------------------------------------
            ! Print iteration results
            !
            IF (io%lnode) THEN
                !
                IF (io%verbosity >= 3) THEN
                    !
                    WRITE (io%debug_unit, 1305) &
                        iter, alpha, beta, rznew, rzold, pAp, delta_qm, delta_en, &
                        tolvelect
                    !
                ELSE IF (io%verbosity >= 1) THEN
                    !
                    WRITE (io%debug_unit, 1306) &
                        iter, alpha, pAp, rzold, delta_qm, delta_en, tolvelect
                    !
                END IF
                !
            END IF
            !
            !----------------------------------------------------------------------------
            ! If residual is small enough exit
            !
            IF (delta_en < tolvelect .AND. iter > 0) THEN
                !
                IF (io%verbosity >= 1 .AND. io%lnode) WRITE (io%debug_unit, 1307)
                !
                EXIT
                !
            ELSE IF (iter == maxiter) THEN
                IF (io%lnode) WRITE (io%unit, 1308)
            END IF
            !
        END DO
        !
        !--------------------------------------------------------------------------------
        ! In PBC the potential need to have zero average
        !
        shift = 0.D0
        !
        SELECT TYPE (core => this%cores%core)
            !
        TYPE IS (core_fft_electrostatics)
            !
            IF (.NOT. (core%use_internal_pbc_corr .OR. &
                       ASSOCIATED(this%cores%correction))) THEN
                shift = -x%integrate() / cell%omega
            END IF
            !
        END SELECT
        !
        x%of_r = x%of_r + shift
        !
        IF (io%lstdout .AND. io%verbosity >= 1) WRITE (io%unit, 1309) delta_en, iter
        !
        CALL r%destroy()
        !
        CALL z%destroy()
        !
        CALL p%destroy()
        !
        CALL Ap%destroy()
        !
        IF (PRESENT(dielectric)) CALL invsqrt%destroy()
        !
        !--------------------------------------------------------------------------------
        !
1300    FORMAT(/, 4('%'), ' COMPUTE ELECTROSTATIC POTENTIAL ', 43('%'),/)
        !
1301    FORMAT(' Sqrt-preconditioned input guess with residual norm = ', E14.6,/)
        !
1302    FORMAT(' Warning: Bad guess with residual norm = ', E14.6, &
               ', reset to no guess',/)
        !
1303    FORMAT('   i | ' &
               '     alpha     |      beta      |      rznew     |      rzold     | ' &
               '      pAp      |    delta_qm    |    delta_en    |       tol', /, &
               1X, 139('-'))
        !
1304    FORMAT('   i | ' &
               '     alpha     |       pAp      |      rzold     |    delta_qm    | ' &
               '   delta_en    |       tol', /, 1X, 105('-'))
        !
1305    FORMAT(1X, I3, 8(' | ', E14.6))
        !
1306    FORMAT(1X, I3, 6(' | ', E14.6))
        !
1307    FORMAT(/, ' Charges are converged, EXIT')
        !
1308    FORMAT(' Warning: Polarization charge not converged',/)
        !
1309    FORMAT('     Polarization accuracy =', 1PE8.1, ', # of iterations = ', i3)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE linearized_pb_gradient_sqrt
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE class_solver_gradient
!----------------------------------------------------------------------------------------
