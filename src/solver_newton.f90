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
        PROCEDURE :: init => init_solver_newton
        !
        PROCEDURE, PRIVATE :: pb_nested_charges
        PROCEDURE, PRIVATE :: pb_nested_density
        !
        GENERIC :: pb_nested => pb_nested_charges, pb_nested_density
        !
        PROCEDURE, PRIVATE :: pb => pb_newton
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
    SUBROUTINE init_solver_newton(this, cores, maxiter, tol, auxiliary)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(core_container), TARGET, INTENT(IN) :: cores
        INTEGER, INTENT(IN) :: maxiter
        REAL(DP), INTENT(IN) :: tol
        CHARACTER(LEN=80), INTENT(IN), OPTIONAL :: auxiliary
        !
        CLASS(solver_newton), INTENT(INOUT) :: this
        !
        !--------------------------------------------------------------------------------
        !
        CALL this%init_iterative(cores, maxiter, tol, auxiliary)
        !
        this%solver_type = 'newton'
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
    SUBROUTINE pb_nested_charges(this, inner, charges, potential)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(electrostatic_solver), INTENT(IN) :: inner
        TYPE(environ_charges), INTENT(IN) :: charges
        !
        CLASS(solver_newton), INTENT(INOUT) :: this
        TYPE(environ_density), INTENT(INOUT) :: potential
        !
        CHARACTER(LEN=80) :: sub_name = 'pb_nested_charges'
        !
        !--------------------------------------------------------------------------------
        !
        CALL env_start_clock(sub_name)
        !
        IF (ASSOCIATED(charges%dielectric)) THEN
            !
            CALL this%pb(inner, potential, charges%density, charges%electrolyte, &
                         charges%dielectric)
            !
        ELSE
            CALL this%pb(inner, potential, charges%density, charges%electrolyte)
        END IF
        !
        CALL env_stop_clock(sub_name)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE pb_nested_charges
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE pb_nested_density(this, inner, charges, potential, electrolyte, &
                                 dielectric)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(electrostatic_solver), INTENT(IN) :: inner
        TYPE(environ_density), INTENT(IN) :: charges
        TYPE(environ_electrolyte), INTENT(IN) :: electrolyte
        TYPE(environ_dielectric), INTENT(IN), OPTIONAL :: dielectric
        !
        CLASS(solver_newton), INTENT(INOUT) :: this
        TYPE(environ_density), INTENT(INOUT) :: potential
        !
        CHARACTER(LEN=80) :: sub_name = 'pb_nested_density'
        !
        !--------------------------------------------------------------------------------
        !
        CALL env_start_clock(sub_name)
        !
        CALL this%pb(inner, potential, charges, electrolyte, dielectric)
        !
        CALL env_stop_clock(sub_name)
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
    SUBROUTINE pb_newton(this, inner, potential, charges, electrolyte, dielectric)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(electrostatic_solver), INTENT(IN) :: inner
        TYPE(environ_density), TARGET, INTENT(IN) :: charges
        TYPE(environ_electrolyte), TARGET, INTENT(IN) :: electrolyte
        TYPE(environ_dielectric), INTENT(IN), OPTIONAL :: dielectric
        !
        CLASS(solver_newton), TARGET, INTENT(INOUT) :: this
        TYPE(environ_density), TARGET, INTENT(INOUT) :: potential
        !
        INTEGER :: iter, itypi, itypj, ir
        REAL(DP) :: totaux, delta_qm, delta_en, kT, z, arg
        !
        TYPE(environ_density) :: residual, rhotot, numerator, denominator, &
                                 cfactor, rhoaux, screening
        !
        REAL(DP), POINTER :: cbulki, cbulkj, zi, zj
        !
        REAL(DP), PARAMETER :: exp_arg_limit = 40.D0
        !
        CHARACTER(LEN=80) :: sub_name = 'pb_newton'
        !
        !--------------------------------------------------------------------------------
        !
        IF (io%verbosity >= 1 .AND. io%lnode) WRITE (io%debug_unit, 1000)
        !
        !--------------------------------------------------------------------------------
        !
        IF (.NOT. ASSOCIATED(charges%cell, electrolyte%gamma%cell)) &
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
        !--------------------------------------------------------------------------------
        !
        ASSOCIATE (cell => charges%cell, &
                   maxiter => this%maxiter, &
                   tol => this%tol, &
                   base => electrolyte%base, &
                   cionmax => electrolyte%base%cionmax, &
                   x => potential, &
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
            x%of_r = 0.D0
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
            DO iter = 1, maxiter
                rhotot%of_r = charges%of_r + rhoaux%of_r + screening%of_r * x%of_r
                residual%of_r = x%of_r
                !
                SELECT TYPE (inner)
                    !
                TYPE IS (solver_gradient)
                    CALL inner%linearized_pb(rhotot, electrolyte, x, dielectric, &
                                             screening)
                    !
                CLASS DEFAULT
                    CALL io%error(sub_name, 'Unexpected inner solver', 1)
                    !
                END SELECT
                !
                residual%of_r = x%of_r - residual%of_r
                rhoaux%of_r = 0.D0
                screening%of_r = 0.D0
                denominator%of_r = 1.D0
                !
                !------------------------------------------------------------------------
                ! General solution for symmetric & asymmetric electrolyte
                !
                DO itypi = 1, base%ntyp
                    zi => base%ioncctype(itypi)%z
                    cbulki => base%ioncctype(itypi)%cbulk
                    !
                    cfactor%of_r = 1.D0
                    !
                    DO ir = 1, cell%ir_end
                        arg = -zi * x%of_r(ir) / kT
                        !
                        IF (arg > exp_arg_limit) THEN
                            cfactor%of_r(ir) = EXP(exp_arg_limit)
                        ELSE IF (arg < -exp_arg_limit) THEN
                            cfactor%of_r(ir) = EXP(-exp_arg_limit)
                        ELSE
                            cfactor%of_r(ir) = EXP(arg)
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
                            denominator%of_r = denominator%of_r - &
                                               cbulki / cionmax * (1.D0 - cfactor%of_r)
                            !
                            DO itypj = 1, base%ntyp
                                zj => base%ioncctype(itypj)%z
                                cbulkj => base%ioncctype(itypj)%cbulk
                                !
                                IF (itypj == itypi) THEN
                                    numerator%of_r = numerator%of_r - cbulkj / cionmax
                                ELSE
                                    !
                                    numerator%of_r = &
                                        numerator%of_r - cbulkj / cionmax * &
                                        (1.D0 - (1.D0 - zj / zi) * &
                                         cfactor%of_r**(zj / zi))
                                    !
                                END IF
                                !
                                NULLIFY (zj)
                                NULLIFY (cbulkj)
                            END DO
                            !
                        CASE ('ions')
                            !
                            denominator%of_r = denominator%of_r - cbulki / cionmax * &
                                               (1.D0 - gam%of_r * cfactor%of_r)
                            !
                            DO itypj = 1, base%ntyp
                                zj => base%ioncctype(itypj)%z
                                cbulkj => base%ioncctype(itypj)%cbulk
                                !
                                IF (itypj == itypi) THEN
                                    numerator%of_r = numerator%of_r - cbulkj / cionmax
                                ELSE
                                    !
                                    numerator%of_r = &
                                        numerator%of_r - cbulkj / cionmax * &
                                        (1.D0 - (1.D0 - zj / zi) * gam%of_r * &
                                         cfactor%of_r**(zj / zi))
                                    !
                                END IF
                                !
                                NULLIFY (zj)
                                NULLIFY (cbulkj)
                            END DO
                            !
                        CASE DEFAULT
                            CALL io%error(sub_name, 'Unexpected electrolyte entropy', 1)
                            !
                        END SELECT
                        !
                    END IF
                    !
                    screening%of_r = screening%of_r + &
                                     cbulki * zi**2 / kT * cfactor%of_r * numerator%of_r
                    !
                    NULLIFY (zi)
                    NULLIFY (cbulki)
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
                        WRITE (io%debug_unit, 1003) iter, delta_qm, delta_en, tol, totaux
                    ELSE IF (io%verbosity >= 1) THEN
                        WRITE (io%debug_unit, 1004) iter, delta_qm, delta_en, tol
                    END IF
                    !
                END IF
                !
                !------------------------------------------------------------------------
                ! If residual is small enough exit
                !
                IF (delta_en < tol .AND. iter > 0) THEN
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
1000    FORMAT(/, 4('%'), ' COMPUTE ELECTROSTATIC POTENTIAL ', 43('%'),/)
        !
1001    FORMAT('   i | outer delta_qm |    delta_en    |       tol      | ' &
               'total iterative electrolyte charge', /, 1X, 91('-'))
        !
1002    FORMAT('   i | outer delta_qm |    delta_en    |       tol', /, 1X, 54('-'))
        !
1003    FORMAT(1X, I3, 4(' | ', E14.6))
        !
1004    FORMAT(1X, I3, 3(' | ', E14.6))
        !
1005    FORMAT(/, ' Newton steps are converged, EXIT')
        !
1006    FORMAT('  Warning: Electrolyte charge not converged',/)
        !
1007    FORMAT('     Electrolyte accuracy =', 1PE8.1, ', # of iterations = ', i3)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE pb_newton
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE class_solver_newton
!----------------------------------------------------------------------------------------
