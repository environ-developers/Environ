!----------------------------------------------------------------------------------------
!
! Copyright (C) 2018-2021 ENVIRON (www.quantum-environ.org)
!
!----------------------------------------------------------------------------------------
!
!     This file is part of Environ version 3.0
!
!     Environ 3.0 is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation, either version 2 of the License, or
!     (at your option) any later version.
!
!     Environ 3.0 is distributed in the hope that it will be useful,
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
MODULE class_solver_newton
    !------------------------------------------------------------------------------------
    !
    USE class_io, ONLY: io
    !
    USE environ_param, ONLY: DP, e2, K_BOLTZMANN_RY, fpi
    !
    USE class_density
    !
    USE class_core_container
    !
    USE class_solver
    USE class_solver_direct
    USE class_solver_gradient
    USE class_solver_iterative
    !
    USE class_charges
    USE class_dielectric
    USE class_electrolyte
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
    TYPE, EXTENDS(solver_iterative), PUBLIC :: solver_newton
        !--------------------------------------------------------------------------------
        !
        !--------------------------------------------------------------------------------
    CONTAINS
        !--------------------------------------------------------------------------------
        !
        PROCEDURE, PRIVATE :: create_solver_newton
        PROCEDURE :: init => init_solver_newton
        !
        PROCEDURE :: pb_nested_charges
        PROCEDURE :: pb_nested_density
        !
        PROCEDURE, PRIVATE :: pb_newton
        !
        !--------------------------------------------------------------------------------
    END TYPE solver_newton
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
    SUBROUTINE create_solver_newton(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(solver_newton), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: routine = 'create_solver_newton'
        !
        !--------------------------------------------------------------------------------
        !
        this%solver_type = 'newton'
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE create_solver_newton
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_solver_newton(this, cores, direct, maxiter, tol, auxiliary)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(core_container), INTENT(IN) :: cores
        TYPE(solver_direct), INTENT(IN) :: direct
        INTEGER, INTENT(IN) :: maxiter
        REAL(DP), INTENT(IN) :: tol
        CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: auxiliary
        !
        CLASS(solver_newton), INTENT(INOUT) :: this
        !
        !--------------------------------------------------------------------------------
        !
        CALL this%create_solver_newton()
        !
        CALL this%init_iterative(cores, direct, maxiter, tol, auxiliary)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_solver_newton
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
    SUBROUTINE pb_nested_charges(this, charges, v, inner)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(solver_newton), INTENT(IN) :: this
        CLASS(electrostatic_solver), OPTIONAL, INTENT(IN) :: inner
        !
        TYPE(environ_density), INTENT(INOUT) :: v
        TYPE(environ_charges), INTENT(INOUT) :: charges
        !
        CHARACTER(LEN=80) :: routine = 'pb_nested_charges'
        !
        !--------------------------------------------------------------------------------
        !
        CALL env_start_clock(routine)
        !
        IF (ASSOCIATED(charges%dielectric)) THEN
            !
            CALL this%pb_newton(v, charges%density, charges%electrolyte, &
                                charges%dielectric, inner=inner)
            !
        ELSE
            CALL this%pb_newton(v, charges%density, charges%electrolyte, inner=inner)
        END IF
        !
        CALL env_stop_clock(routine)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE pb_nested_charges
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE pb_nested_density(this, charges, v, electrolyte, dielectric, inner)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(solver_newton), INTENT(IN) :: this
        TYPE(environ_density), INTENT(IN) :: charges
        TYPE(environ_electrolyte), INTENT(IN) :: electrolyte
        CLASS(electrostatic_solver), OPTIONAL, INTENT(IN) :: inner
        !
        TYPE(environ_density), INTENT(INOUT) :: v
        TYPE(environ_dielectric), OPTIONAL, INTENT(INOUT) :: dielectric
        !
        CHARACTER(LEN=80) :: routine = 'pb_nested_density'
        !
        !--------------------------------------------------------------------------------
        !
        CALL env_start_clock(routine)
        !
        CALL this%pb_newton(v, charges, electrolyte, dielectric, inner=inner)
        !
        CALL env_stop_clock(routine)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE pb_nested_density
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
    SUBROUTINE pb_newton(this, v, charges, electrolyte, dielectric, inner)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(solver_newton), INTENT(IN) :: this
        TYPE(environ_density), INTENT(IN) :: charges
        TYPE(environ_electrolyte), INTENT(IN) :: electrolyte
        CLASS(electrostatic_solver), OPTIONAL, INTENT(IN) :: inner
        !
        TYPE(environ_density), INTENT(INOUT) :: v
        TYPE(environ_dielectric), OPTIONAL, INTENT(INOUT) :: dielectric
        !
        INTEGER :: i, j, k, l
        REAL(DP) :: totaux, delta_qm, delta_en, kT, arg
        !
        TYPE(environ_density) :: residual, rhotot, numerator, denominator, &
                                 cfactor, rhoaux, screening
        !
        REAL(DP), PARAMETER :: exp_arg_limit = 40.D0
        !
        CHARACTER(LEN=80) :: routine = 'pb_newton'
        !
        !--------------------------------------------------------------------------------
        !
        IF (io%verbosity >= 1 .AND. io%lnode) WRITE (io%debug_unit, 1000)
        !
        !--------------------------------------------------------------------------------
        !
        IF (.NOT. ASSOCIATED(charges%cell, electrolyte%gamma%cell)) &
            CALL io%error(routine, "Inconsistent cells of input fields", 1)
        !
        IF (.NOT. ASSOCIATED(charges%cell, v%cell)) &
            CALL io%error(routine, "Inconsistent cells for charges and potential", 1)
        !
        IF (PRESENT(dielectric)) THEN
            !
            IF (.NOT. ASSOCIATED(charges%cell, dielectric%epsilon%cell)) &
                CALL io%error(routine, "Inconsistent cells of input fields", 1)
            !
        END IF
        !
        !--------------------------------------------------------------------------------
        !
        ASSOCIATE (cell => charges%cell, &
                   maxiter => this%maxiter, &
                   tol => this%tol, &
                   base => electrolyte%base, &
                   cionmax => electrolyte%base%cionmax, &
                   gam => electrolyte%gamma)
            !
            !----------------------------------------------------------------------------
            ! Initialize local densities
            !
            CALL cfactor%init(cell)
            !
            CALL numerator%init(cell)
            !
            CALL denominator%init(cell)
            !
            CALL rhoaux%init(cell)
            !
            CALL rhotot%init(cell)
            !
            CALL residual%init(cell)
            !
            CALL screening%init(cell)
            !
            kT = K_BOLTZMANN_RY * base%temperature
            !
            v%of_r = 0.D0
            rhoaux%of_r = 0.D0
            screening%of_r = base%k2 / e2 / fpi * gam%of_r
            residual%of_r = 0.D0
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
            ! Start Newton's algorithm
            !
            DO i = 1, maxiter
                rhotot%of_r = charges%of_r + rhoaux%of_r + screening%of_r * v%of_r
                residual%of_r = v%of_r
                !
                CALL inner%linearized_pb(rhotot, electrolyte, v, dielectric, screening)
                !
                residual%of_r = v%of_r - residual%of_r
                rhoaux%of_r = 0.D0
                screening%of_r = 0.D0
                denominator%of_r = 1.D0
                !
                !------------------------------------------------------------------------
                ! General solution for symmetric & asymmetric electrolyte
                !
                DO j = 1, base%ntyp
                    !
                    ASSOCIATE (zi => electrolyte%base%ioncctype(j)%z, &
                               cbulki => electrolyte%base%ioncctype(j)%cbulk)
                        !
                        cfactor%of_r = 1.D0
                        !
                        DO k = 1, cell%ir_end
                            arg = -zi * v%of_r(k) / kT
                            !
                            IF (arg > exp_arg_limit) THEN
                                cfactor%of_r(k) = EXP(exp_arg_limit)
                            ELSE IF (arg < -exp_arg_limit) THEN
                                cfactor%of_r(k) = EXP(-exp_arg_limit)
                            ELSE
                                cfactor%of_r(k) = EXP(arg)
                            END IF
                            !
                        END DO
                        !
                        rhoaux%of_r = rhoaux%of_r + zi * cbulki * cfactor%of_r
                        numerator%of_r = 1.D0
                        !
                        IF (cionmax > 0.D0) THEN
                            !
                            SELECT CASE (base%electrolyte_entropy)
                                !
                            CASE ('full')
                                !
                                denominator%of_r = &
                                    denominator%of_r - &
                                    cbulki / cionmax * (1.D0 - cfactor%of_r)
                                !
                                DO l = 1, base%ntyp
                                    !
                                    ASSOCIATE (zj => electrolyte%base%ioncctype(l)%z, &
                                               cbulkj => electrolyte%base%ioncctype(l)%cbulk)
                                        !
                                        IF (l == j) THEN
                                            !
                                            numerator%of_r = numerator%of_r - &
                                                             cbulkj / cionmax
                                            !
                                        ELSE
                                            !
                                            numerator%of_r = &
                                                numerator%of_r - cbulkj / cionmax * &
                                                (1.D0 - (1.D0 - zj / zi) * &
                                                 cfactor%of_r**(zj / zi))
                                            !
                                        END IF
                                        !
                                    END ASSOCIATE
                                    !
                                END DO
                                !
                            CASE ('ions')
                                !
                                denominator%of_r = &
                                    denominator%of_r - cbulki / cionmax * &
                                    (1.D0 - gam%of_r * cfactor%of_r)
                                !
                                DO l = 1, base%ntyp
                                    !
                                    ASSOCIATE (zj => electrolyte%base%ioncctype(l)%z, &
                                               cbulkj => electrolyte%base%ioncctype(l)%cbulk)
                                        !
                                        IF (l == j) THEN
                                            !
                                            numerator%of_r = numerator%of_r - &
                                                             cbulkj / cionmax
                                            !
                                        ELSE
                                            !
                                            numerator%of_r = &
                                                numerator%of_r - cbulkj / cionmax * &
                                                (1.D0 - (1.D0 - zj / zi) * gam%of_r * &
                                                 cfactor%of_r**(zj / zi))
                                            !
                                        END IF
                                        !
                                    END ASSOCIATE
                                    !
                                END DO
                                !
                            CASE DEFAULT
                                CALL io%error(routine, "Unexpected electrolyte entropy", 1)
                                !
                            END SELECT
                            !
                        END IF
                        !
                        screening%of_r = &
                            screening%of_r + &
                            cbulki * zi**2 / kT * cfactor%of_r * numerator%of_r
                        !
                    END ASSOCIATE
                    !
                END DO
                !
                rhoaux%of_r = gam%of_r * rhoaux%of_r / denominator%of_r
                screening%of_r = screening%of_r * gam%of_r / denominator%of_r**2
                !
                delta_en = residual%euclidean_norm()
                delta_qm = residual%quadratic_mean()
                totaux = rhoaux%integrate()
                !
                !------------------------------------------------------------------------
                ! Print iteration results
                !
                IF (io%lnode) THEN
                    !
                    IF (io%verbosity >= 3) THEN
                        WRITE (io%debug_unit, 1003) i, delta_qm, delta_en, tol, totaux
                    ELSE IF (io%verbosity >= 1) THEN
                        WRITE (io%debug_unit, 1004) i, delta_qm, delta_en, tol
                    END IF
                    !
                END IF
                !
                !------------------------------------------------------------------------
                ! If residual is small enough exit
                !
                IF (delta_en < tol .AND. i > 0) THEN
                    !
                    IF (io%verbosity >= 1 .AND. io%lnode) WRITE (io%debug_unit, 1005)
                    !
                    EXIT
                    !
                ELSE IF (i == maxiter) THEN
                    IF (io%lnode) WRITE (io%unit, 1006)
                END IF
                !
            END DO
            !
            IF (io%lstdout .AND. io%verbosity >= 1) WRITE (io%unit, 1007) delta_en, i
            !
            !----------------------------------------------------------------------------
            ! Clean up local densities
            !
            CALL cfactor%destroy()
            !
            CALL numerator%destroy()
            !
            CALL denominator%destroy()
            !
            CALL rhoaux%destroy()
            !
            CALL rhotot%destroy()
            !
            CALL residual%destroy()
            !
            CALL screening%destroy()
            !
        END ASSOCIATE
        !
        !--------------------------------------------------------------------------------
        !
1000    FORMAT(/, 4('%'), " COMPUTE ELECTROSTATIC POTENTIAL ", 43('%'),/)
        !
1001    FORMAT("   i | outer delta_qm |    delta_en    |       tol      | " &
               "total iterative electrolyte charge", /, 1X, 91('-'))
        !
1002    FORMAT("   i | outer delta_qm |    delta_en    |       tol", /, 1X, 54('-'))
        !
1003    FORMAT(1X, I3, 4(" | ", E14.6))
        !
1004    FORMAT(1X, I3, 3(" | ", E14.6))
        !
1005    FORMAT(/, " Newton steps are converged, EXIT")
        !
1006    FORMAT("  Warning: Electrolyte charge not converged",/)
        !
1007    FORMAT("     Electrolyte accuracy =", 1PE8.1, ", # of iterations = ", i3)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE pb_newton
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE class_solver_newton
!----------------------------------------------------------------------------------------
