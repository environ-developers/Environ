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
    USE class_density
    USE class_gradient
    !
    USE class_core_container
    !
    USE class_core_fft
    !
    USE class_solver_direct
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
        PROCEDURE :: generalized_charges
        PROCEDURE :: generalized_density
        !
        PROCEDURE :: linearized_pb_charges
        PROCEDURE :: linearized_pb_density
        !
        PROCEDURE, PRIVATE :: generalized_none
        PROCEDURE, PRIVATE :: generalized_sqrt
        PROCEDURE, PRIVATE :: generalized_left
        !
        PROCEDURE, PRIVATE :: linearized_pb_sqrt
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
                                    screening_type, screening, cores, direct, maxiter, &
                                    tol, auxiliary)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(core_container), INTENT(IN) :: cores
        TYPE(solver_direct), TARGET, INTENT(IN) :: direct
        LOGICAL, INTENT(IN) :: lconjugate
        INTEGER, INTENT(IN) :: maxiter
        REAL(DP), INTENT(IN) :: tol, step, screening
        CHARACTER(LEN=*), INTENT(IN) :: step_type, preconditioner, screening_type
        CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: auxiliary
        !
        CLASS(solver_gradient), INTENT(INOUT) :: this
        !
        !--------------------------------------------------------------------------------
        !
        CALL this%init_iterative(cores, direct, maxiter, tol, auxiliary)
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
    SUBROUTINE generalized_charges(this, charges, v)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(solver_gradient), INTENT(IN) :: this
        !
        TYPE(environ_density), INTENT(INOUT) :: v
        TYPE(environ_charges), INTENT(INOUT) :: charges
        !
        CHARACTER(LEN=80) :: sub_name = 'generalized_charges'
        !
        !--------------------------------------------------------------------------------
        !
        CALL env_start_clock(sub_name)
        !
        v%of_r = 0.D0
        !
        IF (this%auxiliary == 'none') THEN
            !
            SELECT CASE (this%preconditioner) ! #TODO none and left need work
                !
            CASE ('none')
                !
                CALL this%generalized_none(charges%density, charges%dielectric, v, &
                                           charges%electrolyte, charges%semiconductor)
                !
            CASE ('sqrt')
                !
                CALL this%generalized_sqrt(charges%density, charges%dielectric, v, &
                                           charges%electrolyte, charges%semiconductor)
                !
            CASE ('left')
                !
                CALL this%generalized_left(charges%density, charges%dielectric, v, &
                                           charges%electrolyte, charges%semiconductor)
                !
            CASE DEFAULT
                CALL io%error(sub_name, "Unexpected 'preconditioner'", 1)
                !
            END SELECT
            !
        ELSE
            !
            CALL io%error(sub_name, 'Option not yet implemented', 1)
            !
            ! CALL generalized_rhoaux(charges, dielectric, v) #TODO future-work
            !
        END IF
        !
        CALL env_stop_clock(sub_name)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE generalized_charges
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE generalized_density(this, charges, dielectric, v, electrolyte, &
                                   semiconductor)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(solver_gradient), INTENT(IN) :: this
        TYPE(environ_density), INTENT(IN) :: charges
        TYPE(environ_electrolyte), INTENT(IN), OPTIONAL :: electrolyte
        TYPE(environ_semiconductor), INTENT(IN), OPTIONAL :: semiconductor
        !
        TYPE(environ_density), INTENT(INOUT) :: v
        TYPE(environ_dielectric), INTENT(INOUT) :: dielectric
        !
        CHARACTER(LEN=80) :: sub_name = 'generalized_density'
        !
        !--------------------------------------------------------------------------------
        !
        CALL env_start_clock(sub_name)
        !
        v%of_r = 0.D0
        !
        IF (this%auxiliary == 'none') THEN
            !
            SELECT CASE (this%preconditioner) ! #TODO none and left need work
                !
            CASE ('none')
                CALL this%generalized_none(charges, dielectric, v, electrolyte)
                !
            CASE ('sqrt')
                CALL this%generalized_sqrt(charges, dielectric, v, electrolyte)
                !
            CASE ('left')
                CALL this%generalized_left(charges, dielectric, v, electrolyte)
                !
            CASE DEFAULT
                CALL io%error(sub_name, "Unexpected 'preconditioner'", 1)
                !
            END SELECT
            !
        ELSE
            !
            CALL io%error(sub_name, 'Option not yet implemented', 1)
            !
            ! CALL generalized_rhoaux(charges, dielectric, v) #TODO future-work
            !
        END IF
        !
        CALL env_stop_clock(sub_name)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE generalized_density
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE linearized_pb_charges(this, charges, v, screening)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(solver_gradient), INTENT(IN) :: this
        TYPE(environ_density), INTENT(IN), OPTIONAL :: screening
        !
        TYPE(environ_density), INTENT(INOUT) :: v
        TYPE(environ_charges), INTENT(INOUT) :: charges
        !
        TYPE(environ_density) :: local_screening
        !
        CHARACTER(LEN=80) :: sub_name = 'linearized_pb_charges'
        !
        !--------------------------------------------------------------------------------
        !
        CALL env_start_clock(sub_name)
        !
        CALL local_screening%init(v%cell)
        !
        IF (PRESENT(screening)) THEN
            local_screening%of_r = screening%of_r ! external screening
        ELSE
            !
            !----------------------------------------------------------------------------
            ! Screening as in linearized pb problem
            !
            ASSOCIATE (base => charges%electrolyte%base, &
                       electrolyte => charges%electrolyte)
                !
                IF (base%electrolyte_entropy == 'ions' .AND. base%cionmax > 0.D0) THEN
                    !
                    local_screening%of_r = &
                        base%k2 / e2 / fpi * electrolyte%gamma%of_r / &
                        (1.D0 - SUM(base%ioncctype%cbulk) / base%cionmax * &
                         (1.D0 - electrolyte%gamma%of_r))
                    !
                ELSE
                    local_screening%of_r = base%k2 / e2 / fpi * electrolyte%gamma%of_r
                END IF
                !
            END ASSOCIATE
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
                    CALL this%linearized_pb_sqrt(charges%density, local_screening, v, &
                                                 charges%dielectric)
                    !
                CASE DEFAULT
                    CALL io%error(sub_name, "Unexpected 'preconditioner'", 1)
                    !
                END SELECT
                !
            ELSE
                !
                CALL this%linearized_pb_sqrt(charges%density, local_screening, v)
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
    END SUBROUTINE linearized_pb_charges
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE linearized_pb_density(this, charges, electrolyte, v, dielectric, &
                                     screening)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(solver_gradient), INTENT(IN) :: this
        TYPE(environ_density), INTENT(IN) :: charges
        TYPE(environ_electrolyte), INTENT(IN) :: electrolyte
        TYPE(environ_density), INTENT(IN), OPTIONAL :: screening
        !
        TYPE(environ_density), INTENT(INOUT) :: v
        TYPE(environ_dielectric), INTENT(INOUT), OPTIONAL :: dielectric
        !
        TYPE(environ_density) :: local_screening
        !
        CHARACTER(LEN=80) :: sub_name = 'linearized_pb_density'
        !
        !--------------------------------------------------------------------------------
        !
        CALL env_start_clock(sub_name)
        !
        CALL local_screening%init(v%cell)
        !
        IF (PRESENT(screening)) THEN
            local_screening%of_r = screening%of_r ! external screening
        ELSE
            !
            !----------------------------------------------------------------------------
            ! Screening as in linearized pb problem
            !
            ASSOCIATE (base => electrolyte%base)
                !
                IF (base%electrolyte_entropy == 'ions' .AND. base%cionmax > 0.D0) THEN
                    !
                    local_screening%of_r = &
                        base%k2 / e2 / fpi * electrolyte%gamma%of_r / &
                        (1.D0 - SUM(base%ioncctype%cbulk) / &
                         base%cionmax * (1.D0 - electrolyte%gamma%of_r))
                    !
                ELSE
                    local_screening%of_r = base%k2 / e2 / fpi * electrolyte%gamma%of_r
                END IF
                !
            END ASSOCIATE
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
                    CALL this%linearized_pb_sqrt(charges, local_screening, v, &
                                                 dielectric)
                    !
                CASE DEFAULT
                    CALL io%error(sub_name, "Unexpected 'preconditioner'", 1)
                    !
                END SELECT
                !
            ELSE
                CALL this%linearized_pb_sqrt(charges, local_screening, v)
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
    END SUBROUTINE linearized_pb_density
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
    SUBROUTINE generalized_none(this, charges, dielectric, v, electrolyte, semiconductor)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(solver_gradient), TARGET, INTENT(IN) :: this
        TYPE(environ_density), TARGET, INTENT(IN) :: charges
        TYPE(environ_dielectric), TARGET, INTENT(IN) :: dielectric
        TYPE(environ_electrolyte), INTENT(IN), OPTIONAL :: electrolyte
        TYPE(environ_semiconductor), INTENT(INOUT), OPTIONAL :: semiconductor
        !
        TYPE(environ_density), TARGET, INTENT(INOUT) :: v
        !
        INTEGER :: iter
        REAL(DP) :: rznew, rzold, alpha, beta, pAp, delta_qm, delta_en
        TYPE(environ_density) :: r, z, p, Ap, l
        TYPE(environ_gradient) :: g
        !
        CHARACTER(LEN=80) :: sub_name = 'generalized_none'
        !
        !--------------------------------------------------------------------------------
        !
        IF (io%verbosity >= 1 .AND. io%lnode) WRITE (io%debug_unit, 1000)
        !
        !--------------------------------------------------------------------------------
        !
        IF (.NOT. ASSOCIATED(charges%cell, dielectric%epsilon%cell)) &
            CALL io%error(sub_name, 'Inconsistent cells of input fields', 1)
        !
        IF (.NOT. ASSOCIATED(charges%cell, v%cell)) &
            CALL io%error(sub_name, 'Inconsistent cells for charges and potential', 1)
        !
        !--------------------------------------------------------------------------------
        !
        ASSOCIATE (cell => charges%cell, &
                   gradeps => dielectric%gradient, &
                   maxiter => this%maxiter, &
                   tolvelect => this%tol)
            !
            !----------------------------------------------------------------------------
            ! Initialize local densities
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
            !----------------------------------------------------------------------------
            ! Starting guess from new input and previous solution(s)
            !
            IF (v%lupdate) THEN
                v%lupdate = .FALSE.
            END IF
            !
            IF (.NOT. v%lupdate) THEN
                v%lupdate = .TRUE.
                v%of_r = 0.D0
                r = charges
                rzold = 0.D0
            END IF
            !
            !----------------------------------------------------------------------------
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
            !----------------------------------------------------------------------------
            ! Start gradient descent
            !
            DO iter = 1, maxiter
                !
                !------------------------------------------------------------------------
                ! Apply preconditioner to new state
                !
                z%of_r = r%of_r ! no preconditioner
                !
                rznew = r%scalar_product(z)
                !
                IF (ABS(rznew) < 1.D-30) &
                    CALL io%error(sub_name, 'Null step in gradient descent iteration', 1)
                !
                !------------------------------------------------------------------------
                ! Conjugate gradient or steepest descent input
                !
                IF (this%lconjugate .AND. ABS(rzold) > 1.D-30) THEN
                    beta = rznew / rzold
                ELSE
                    beta = 0.D0
                END IF
                !
                rzold = rznew
                p%of_r = z%of_r + beta * p%of_r
                !
                !------------------------------------------------------------------------
                ! Apply operator to conjugate direction
                ! NOTE: the following steps should be extended to account for
                !       different cores
                !
                CALL this%cores%derivatives%gradient(p, g)
                !
                CALL this%cores%derivatives%laplacian(p, l)
                !
                Ap%of_r = dielectric%epsilon%of_r * l%of_r + &
                          gradeps%of_r(1, :) * g%of_r(1, :) + &
                          gradeps%of_r(2, :) * g%of_r(2, :) + &
                          gradeps%of_r(3, :) * g%of_r(3, :)
                !
                Ap%of_r = -Ap%of_r / fpi / e2
                !
                !------------------------------------------------------------------------
                ! Step downhill
                !
                pAp = p%scalar_product(Ap)
                alpha = rzold / pAp
                !
                v%of_r = v%of_r + alpha * p%of_r
                r%of_r = r%of_r - alpha * Ap%of_r
                !
                delta_en = r%euclidean_norm()
                delta_qm = r%quadratic_mean()
                !
                !------------------------------------------------------------------------
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
                !------------------------------------------------------------------------
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
            !----------------------------------------------------------------------------
            ! Clean up local densities
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
        END ASSOCIATE
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
    END SUBROUTINE generalized_none
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE generalized_sqrt(this, charges, dielectric, v, electrolyte, semiconductor)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(solver_gradient), TARGET, INTENT(IN) :: this
        TYPE(environ_density), TARGET, INTENT(IN) :: charges
        TYPE(environ_dielectric), TARGET, INTENT(IN) :: dielectric
        TYPE(environ_electrolyte), INTENT(IN), OPTIONAL :: electrolyte
        TYPE(environ_semiconductor), INTENT(IN), OPTIONAL :: semiconductor
        !
        TYPE(environ_density), TARGET, INTENT(INOUT) :: v
        !
        INTEGER :: iter
        REAL(DP) :: rznew, rzold, alpha, beta, pAp, delta_qm, delta_en, shift
        TYPE(environ_density) :: r, z, p, Ap, invsqrt
        !
        CHARACTER(LEN=80) :: sub_name = 'generalized_sqrt'
        !
        !--------------------------------------------------------------------------------
        !
        IF (io%verbosity >= 1 .AND. io%lnode) WRITE (io%debug_unit, 1100)
        !
        !--------------------------------------------------------------------------------
        !
        IF (.NOT. ASSOCIATED(charges%cell, dielectric%epsilon%cell)) &
            CALL io%error(sub_name, 'Inconsistent cells of input fields', 1)
        !
        IF (.NOT. ASSOCIATED(charges%cell, v%cell)) &
            CALL io%error(sub_name, 'Inconsistent cells for charges and potential', 1)
        !
        !--------------------------------------------------------------------------------
        !
        ASSOCIATE (cell => charges%cell, &
                   b => charges, &
                   factsqrt => dielectric%factsqrt, &
                   maxiter => this%maxiter, &
                   tolvelect => this%tol)
            !
            !----------------------------------------------------------------------------
            ! Initialize local densities
            !
            CALL invsqrt%init(cell)
            !
            invsqrt%of_r = 1.D0 / SQRT(dielectric%epsilon%of_r)
            !
            CALL r%init(cell)
            !
            CALL z%init(cell)
            !
            CALL p%init(cell)
            !
            CALL Ap%init(cell)
            !
            !----------------------------------------------------------------------------
            ! Starting guess from new input and previous solution(s)
            !
            IF (v%lupdate) THEN
                r%of_r = b%of_r - factsqrt%of_r * v%of_r
                !
                !------------------------------------------------------------------------
                ! Preconditioning step
                !
                z%of_r = r%of_r * invsqrt%of_r
                !
                CALL this%direct%poisson(z, z, electrolyte, semiconductor)
                !
                z%of_r = z%of_r * invsqrt%of_r
                !
                rzold = r%scalar_product(z)
                !
                IF (ABS(rzold) < 1.D-30) &
                    CALL io%error(sub_name, 'Null step in gradient descent iteration', 1)
                !
                r%of_r = factsqrt%of_r * (v%of_r - z%of_r)
                delta_en = r%euclidean_norm()
                delta_qm = r%quadratic_mean()
                !
                IF (delta_en < 1.D-02) THEN
                    !
                    IF (io%verbosity >= 1 .AND. io%lnode) &
                        WRITE (io%debug_unit, 1101) delta_en
                    !
                    v%of_r = z%of_r
                ELSE
                    !
                    IF (io%verbosity >= 1 .AND. io%lnode) &
                        WRITE (io%debug_unit, 1102) delta_en
                    !
                    v%lupdate = .FALSE.
                END IF
                !
            END IF
            !
            IF (.NOT. v%lupdate) THEN
                v%lupdate = .TRUE.
                v%of_r = 0.D0
                r%of_r = b%of_r
                rzold = 0.D0
            END IF
            !
            !----------------------------------------------------------------------------
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
            !----------------------------------------------------------------------------
            ! Start gradient descent
            !
            DO iter = 1, maxiter
                !
                !------------------------------------------------------------------------
                ! Apply preconditioner to new state
                !
                z%of_r = r%of_r * invsqrt%of_r
                !
                CALL this%direct%poisson(z, z, electrolyte, semiconductor)
                !
                z%of_r = z%of_r * invsqrt%of_r
                !
                rznew = r%scalar_product(z)
                !
                IF (ABS(rznew) < 1.D-30) &
                    CALL io%error(sub_name, 'Null step in gradient descent iteration', 1)
                !
                !------------------------------------------------------------------------
                ! Conjugate gradient or steepest descent input
                !
                IF (this%lconjugate .AND. ABS(rzold) > 1.D-30) THEN
                    beta = rznew / rzold
                ELSE
                    beta = 0.D0
                END IF
                !
                rzold = rznew
                p%of_r = z%of_r + beta * p%of_r
                !
                !------------------------------------------------------------------------
                !
                Ap%of_r = factsqrt%of_r * z%of_r + r%of_r + beta * Ap%of_r
                ! apply operator to conjugate direction
                !
                !------------------------------------------------------------------------
                ! Step downhill
                !
                pAp = p%scalar_product(Ap)
                alpha = rzold / pAp
                !
                v%of_r = v%of_r + alpha * p%of_r
                r%of_r = r%of_r - alpha * Ap%of_r
                !
                delta_en = r%euclidean_norm()
                delta_qm = r%quadratic_mean()
                !
                !------------------------------------------------------------------------
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
                !------------------------------------------------------------------------
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
            !----------------------------------------------------------------------------
            ! In PBC the potential need to have zero average
            !
            shift = 0.D0
            !
            IF (.NOT. (this%cores%internal_correction .OR. this%cores%has_corrections)) &
                shift = -v%integrate() / cell%omega
            !
            v%of_r = v%of_r + shift
            !
            IF (io%lstdout .AND. io%verbosity >= 1) WRITE (io%unit, 1109) delta_en, iter
            !
            !----------------------------------------------------------------------------
            ! Clean up local densities
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
        END ASSOCIATE
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
    END SUBROUTINE generalized_sqrt
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE generalized_left(this, charges, dielectric, v, electrolyte, semiconductor)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(solver_gradient), TARGET, INTENT(IN) :: this
        TYPE(environ_density), TARGET, INTENT(IN) :: charges
        TYPE(environ_dielectric), TARGET, INTENT(IN) :: dielectric
        TYPE(environ_electrolyte), INTENT(IN), OPTIONAL :: electrolyte
        TYPE(environ_semiconductor), INTENT(IN), OPTIONAL :: semiconductor
        !
        TYPE(environ_density), TARGET, INTENT(INOUT) :: v
        !
        INTEGER :: iter
        REAL(DP) :: rznew, rzold, alpha, beta, pAp, delta_en, delta_qm
        TYPE(environ_density) :: r, z, p, Ap
        TYPE(environ_gradient) :: g
        !
        CHARACTER(LEN=80) :: sub_name = 'generalized_left'
        !
        !--------------------------------------------------------------------------------
        !
        IF (io%verbosity >= 1 .AND. io%lnode) WRITE (io%debug_unit, 1200)
        !
        !--------------------------------------------------------------------------------
        ! Check that fields have the same defintion domain
        !
        IF (.NOT. ASSOCIATED(charges%cell, dielectric%epsilon%cell)) &
            CALL io%error(sub_name, 'Inconsistent cells of input fields', 1)
        !
        IF (.NOT. ASSOCIATED(charges%cell, v%cell)) &
            CALL io%error(sub_name, 'Inconsistent cells for charges and potential', 1)
        !
        !--------------------------------------------------------------------------------
        !
        ASSOCIATE (cell => charges%cell, &
                   gradeps => dielectric%gradient, &
                   maxiter => this%maxiter, &
                   tolvelect => this%tol)
            !
            !----------------------------------------------------------------------------
            ! Initialize local densities
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
            !----------------------------------------------------------------------------
            ! Starting guess from new input and previous solution(s)
            !
            v%of_r = 0.D0
            r%of_r = charges%of_r
            rzold = 0.D0
            !
            !----------------------------------------------------------------------------
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
            !----------------------------------------------------------------------------
            ! Start gradient descent
            !
            DO iter = 1, maxiter
                !
                !------------------------------------------------------------------------
                ! Apply preconditioner to new state
                !
                z%of_r = r%of_r / dielectric%epsilon%of_r
                !
                CALL this%direct%poisson(z, z, electrolyte, semiconductor)
                !
                rznew = r%scalar_product(z)
                !
                IF (ABS(rznew) < 1.D-30) &
                    CALL io%error(sub_name, 'Null step in gradient descent iteration', 1)
                !
                !------------------------------------------------------------------------
                ! Conjugate gradient or steepest descent input
                !
                IF (this%lconjugate .AND. ABS(rzold) > 1.D-30) THEN
                    beta = rznew / rzold
                ELSE
                    beta = 0.D0
                END IF
                !
                rzold = rznew
                p%of_r = z%of_r + beta * p%of_r
                !
                !------------------------------------------------------------------------
                ! Apply operator to conjugate direction
                ! NOTE: the following steps should be extended to account for
                !       different cores
                !
                CALL this%cores%derivatives%gradient(z, g)
                !
                g%of_r = g%of_r / fpi / e2
                !
                Ap%of_r = beta * Ap%of_r - r%of_r + &
                          gradeps%of_r(1, :) * g%of_r(1, :) + &
                          gradeps%of_r(2, :) * g%of_r(2, :) + &
                          gradeps%of_r(3, :) * g%of_r(3, :)
                !
                !------------------------------------------------------------------------
                ! Step downhill
                !
                pAp = p%scalar_product(Ap)
                alpha = rzold / pAp
                !
                v%of_r = v%of_r + alpha * p%of_r
                r%of_r = r%of_r - alpha * Ap%of_r
                !
                delta_en = r%euclidean_norm()
                delta_qm = r%quadratic_mean()
                !
                !------------------------------------------------------------------------
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
                !------------------------------------------------------------------------
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
            !----------------------------------------------------------------------------
            ! Clean up local densities
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
        END ASSOCIATE
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
    END SUBROUTINE generalized_left
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE linearized_pb_sqrt(this, charges, screening, v, dielectric)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(solver_gradient), TARGET, INTENT(IN) :: this
        TYPE(environ_density), TARGET, INTENT(IN) :: charges
        TYPE(environ_density), TARGET, INTENT(IN) :: screening
        TYPE(environ_dielectric), TARGET, INTENT(IN), OPTIONAL :: dielectric
        !
        TYPE(environ_density), TARGET, INTENT(INOUT) :: v
        !
        TYPE(environ_density), POINTER :: eps, factsqrt
        !
        INTEGER :: iter
        REAL(DP) :: rznew, rzold, alpha, beta, pAp, delta_qm, delta_en, shift
        TYPE(environ_density) :: r, z, p, Ap, invsqrt
        !
        CHARACTER(LEN=80) :: sub_name = 'linearized_pb_sqrt'
        !
        !--------------------------------------------------------------------------------
        !
        IF (io%verbosity >= 1 .AND. io%lnode) WRITE (io%debug_unit, 1300)
        !
        !--------------------------------------------------------------------------------
        !
        IF (.NOT. ASSOCIATED(charges%cell, screening%cell)) &
            CALL io%error(sub_name, 'Inconsistent cells of input fields', 1)
        !
        IF (.NOT. ASSOCIATED(charges%cell, v%cell)) &
            CALL io%error(sub_name, 'Inconsistent cells for charges and potential', 1)
        !
        IF (PRESENT(dielectric)) THEN
            !
            IF (.NOT. ASSOCIATED(charges%cell, dielectric%epsilon%cell)) &
                CALL io%error(sub_name, 'Inconsistent cells of input fields', 1)
            !
        END IF
        !
        !--------------------------------------------------------------------------------
        !
        ASSOCIATE (cell => charges%cell, &
                   b => charges, &
                   scr => screening, &
                   maxiter => this%maxiter, &
                   tolvelect => this%tol)
            !
            !----------------------------------------------------------------------------
            ! Initialize local densities
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
            !----------------------------------------------------------------------------
            ! Starting guess from new input and previous solution(s)
            !
            IF (v%lupdate) THEN
                !
                IF (PRESENT(dielectric)) THEN
                    r%of_r = b%of_r - (factsqrt%of_r + scr%of_r) * v%of_r
                ELSE
                    r%of_r = b%of_r - scr%of_r * v%of_r
                END IF
                !
                !------------------------------------------------------------------------
                ! Preconditioning step
                !
                IF (PRESENT(dielectric)) THEN
                    z%of_r = r%of_r * invsqrt%of_r
                    !
                    CALL this%direct%poisson(z, z)
                    !
                    z%of_r = z%of_r * invsqrt%of_r
                ELSE
                    z%of_r = r%of_r
                    !
                    CALL this%direct%poisson(z, z)
                    !
                END IF
                !
                rzold = r%scalar_product(z)
                !
                IF (ABS(rzold) < 1.D-30) &
                    CALL io%error(sub_name, 'Null step in gradient descent iteration', 1)
                !
                IF (PRESENT(dielectric)) THEN
                    r%of_r = (factsqrt%of_r + scr%of_r) * (v%of_r - z%of_r)
                ELSE
                    r%of_r = scr%of_r * (v%of_r - z%of_r)
                END IF
                !
                delta_en = r%euclidean_norm()
                delta_qm = r%quadratic_mean()
                !
                IF (delta_en < 1.D-02) THEN
                    !
                    IF (io%verbosity >= 1 .AND. io%lnode) &
                        WRITE (io%debug_unit, 1301) delta_en
                    !
                    v%of_r = z%of_r
                ELSE
                    !
                    IF (io%verbosity >= 1 .AND. io%lnode) &
                        WRITE (io%debug_unit, 1302) delta_en
                    !
                    v%lupdate = .FALSE.
                END IF
                !
            END IF
            !
            IF (.NOT. v%lupdate) THEN
                v%lupdate = .TRUE.
                v%of_r = 0.D0
                r%of_r = b%of_r
                rzold = 0.D0
            END IF
            !
            !----------------------------------------------------------------------------
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
            !----------------------------------------------------------------------------
            ! Start gradient descent
            !
            DO iter = 1, maxiter
                !
                !------------------------------------------------------------------------
                ! Apply preconditioner to new state
                !
                IF (PRESENT(dielectric)) THEN
                    z%of_r = r%of_r * invsqrt%of_r
                    !
                    CALL this%direct%poisson(z, z)
                    !
                    z%of_r = z%of_r * invsqrt%of_r
                ELSE
                    z%of_r = r%of_r
                    !
                    CALL this%direct%poisson(z, z)
                    !
                END IF
                !
                rznew = r%scalar_product(z)
                !
                IF (ABS(rznew) < 1.D-30) &
                    CALL io%error(sub_name, 'Null step in gradient descent iteration', 1)
                !
                !------------------------------------------------------------------------
                ! Conjugate gradient or steepest descent input
                !
                IF (this%lconjugate .AND. ABS(rzold) > 1.D-30) THEN
                    beta = rznew / rzold
                ELSE
                    beta = 0.D0
                END IF
                !
                rzold = rznew
                p%of_r = z%of_r + beta * p%of_r
                !
                !------------------------------------------------------------------------
                ! Apply operator to conjugate direction
                !
                IF (PRESENT(dielectric)) THEN
                    !
                    Ap%of_r = &
                        (factsqrt%of_r + scr%of_r) * z%of_r + r%of_r + beta * Ap%of_r
                    !
                ELSE
                    Ap%of_r = scr%of_r * z%of_r + r%of_r + beta * Ap%of_r
                END IF
                !
                !------------------------------------------------------------------------
                ! Step downhill
                !
                pAp = p%scalar_product(Ap)
                alpha = rzold / pAp
                !
                v%of_r = v%of_r + alpha * p%of_r
                r%of_r = r%of_r - alpha * Ap%of_r
                !
                delta_en = r%euclidean_norm()
                delta_qm = r%quadratic_mean()
                !
                !------------------------------------------------------------------------
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
                !------------------------------------------------------------------------
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
            !----------------------------------------------------------------------------
            ! In PBC the potential need to have zero average
            !
            shift = 0.D0
            !
            IF (.NOT. (this%cores%internal_correction .OR. this%cores%has_corrections)) &
                shift = -v%integrate() / cell%omega
            !
            v%of_r = v%of_r + shift
            !
            IF (io%lstdout .AND. io%verbosity >= 1) WRITE (io%unit, 1309) delta_en, iter
            !
            !----------------------------------------------------------------------------
            ! Clean up local densities
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
        END ASSOCIATE
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
    END SUBROUTINE linearized_pb_sqrt
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE class_solver_gradient
!----------------------------------------------------------------------------------------
