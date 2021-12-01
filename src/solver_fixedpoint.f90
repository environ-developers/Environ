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
MODULE class_solver_fixedpoint
    !------------------------------------------------------------------------------------
    !
    USE class_io, ONLY: io
    !
    USE environ_param, ONLY: DP, e2, K_BOLTZMANN_RY, fpi
    !
    USE class_density
    USE class_gradient
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
    TYPE, EXTENDS(solver_iterative), PUBLIC :: solver_fixedpoint
        !--------------------------------------------------------------------------------
        !
        CHARACTER(LEN=80) :: mix_type
        REAL(DP) :: mix
        INTEGER :: ndiis
        !
        !--------------------------------------------------------------------------------
    CONTAINS
        !--------------------------------------------------------------------------------
        !
        PROCEDURE :: init => init_solver_fixedpoint
        !
        PROCEDURE :: generalized_charges
        PROCEDURE :: generalized_density
        !
        PROCEDURE :: pb_nested_charges
        PROCEDURE :: pb_nested_density
        !
        PROCEDURE, PRIVATE :: pb_fixedpoint
        PROCEDURE, PRIVATE :: generalized_fixedpoint
        !
        !--------------------------------------------------------------------------------
    END TYPE solver_fixedpoint
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
    SUBROUTINE init_solver_fixedpoint(this, mix_type, mix, ndiis, cores, direct, &
                                      maxiter, tol, auxiliary)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(core_container), TARGET, INTENT(IN) :: cores
        TYPE(solver_direct), TARGET, INTENT(IN) :: direct
        INTEGER, INTENT(IN) :: ndiis, maxiter
        REAL(DP), INTENT(IN) :: tol, mix
        CHARACTER(LEN=*), INTENT(IN) :: mix_type
        CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: auxiliary
        !
        CLASS(solver_fixedpoint), INTENT(INOUT) :: this
        !
        !--------------------------------------------------------------------------------
        !
        CALL this%init_iterative(cores, direct, maxiter, tol, auxiliary)
        !
        this%solver_type = 'fixed-point'
        !
        this%mix_type = mix_type
        this%mix = mix
        this%ndiis = ndiis
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_solver_fixedpoint
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
        CLASS(solver_fixedpoint), INTENT(IN) :: this
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
        IF (this%auxiliary == 'full') THEN
            !
            CALL this%generalized_fixedpoint(charges%density, charges%dielectric, v, &
                                             charges%electrolyte, charges%semiconductor)
            !
        ELSE
            !
            CALL io%error(sub_name, 'Option not yet implemented', 1)
            !
            ! CALL generalized_iterative_velect(charges, dielectric, v) #TODO future-work
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
    SUBROUTINE generalized_density(this, charges, dielectric, v, &
                                   electrolyte, semiconductor)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(solver_fixedpoint), INTENT(IN) :: this
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
        IF (this%auxiliary == 'full') THEN
            !
            CALL this%generalized_fixedpoint(charges, dielectric, v, electrolyte, &
                                             semiconductor)
            !
        ELSE
            !
            CALL io%error(sub_name, 'Option not yet implemented', 1)
            !
            ! CALL generalized_iterative_velect(charges, dielectric, v) #TODO future-work
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
    SUBROUTINE pb_nested_charges(this, charges, v, inner)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(solver_fixedpoint), INTENT(IN) :: this
        CLASS(electrostatic_solver), OPTIONAL, INTENT(IN) :: inner
        !
        TYPE(environ_density), INTENT(INOUT) :: v
        TYPE(environ_charges), INTENT(INOUT) :: charges
        !
        CHARACTER(LEN=80) :: sub_name = 'pb_nested_charges'
        !
        !--------------------------------------------------------------------------------
        !
        CALL env_start_clock(sub_name)
        !
        IF (this%auxiliary == 'ioncc') THEN
            !
            IF (ASSOCIATED(charges%dielectric)) THEN
                !
                IF (.NOT. PRESENT(inner)) &
                    CALL io%error(sub_name, 'Missing inner solver', 1)
                !
                CALL this%pb_fixedpoint(v, charges%density, charges%electrolyte, &
                                        charges%dielectric, inner=inner)
                !
            ELSE
                CALL this%pb_fixedpoint(v, charges%density, charges%electrolyte)
            END IF
            !
        ELSE
            CALL io%error(sub_name, 'Option not available', 1)
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
    SUBROUTINE pb_nested_density(this, charges, v, electrolyte, dielectric, inner)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(solver_fixedpoint), INTENT(IN) :: this
        TYPE(environ_density), INTENT(IN) :: charges
        TYPE(environ_electrolyte), INTENT(IN) :: electrolyte
        CLASS(electrostatic_solver), INTENT(IN), OPTIONAL :: inner
        !
        TYPE(environ_density), INTENT(INOUT) :: v
        TYPE(environ_dielectric), INTENT(INOUT), OPTIONAL :: dielectric
        !
        CHARACTER(LEN=80) :: sub_name = 'pb_nested_density'
        !
        !--------------------------------------------------------------------------------
        !
        CALL env_start_clock(sub_name)
        !
        IF (this%auxiliary == 'ioncc') THEN
            !
            IF (PRESENT(dielectric)) THEN
                !
                IF (.NOT. PRESENT(inner)) &
                    CALL io%error(sub_name, 'Missing inner setup', 1)
                !
                CALL this%pb_fixedpoint(v, charges, electrolyte, dielectric, inner=inner)
                !
            ELSE
                CALL this%pb_fixedpoint(v, charges, electrolyte)
            END IF
            !
        ELSE
            CALL io%error(sub_name, 'Option not available', 1)
        END IF
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
    SUBROUTINE generalized_fixedpoint(this, charges, dielectric, v, electrolyte, &
                                      semiconductor)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(solver_fixedpoint), TARGET, INTENT(IN) :: this
        TYPE(environ_density), TARGET, INTENT(IN) :: charges
        TYPE(environ_electrolyte), INTENT(IN), OPTIONAL :: electrolyte
        TYPE(environ_semiconductor), INTENT(IN), OPTIONAL :: semiconductor
        !
        TYPE(environ_density), TARGET, INTENT(INOUT) :: v
        TYPE(environ_dielectric), TARGET, INTENT(INOUT) :: dielectric
        !
        INTEGER :: iter
        REAL(DP) :: total, totpol, totzero, totiter, delta_qm, delta_en
        TYPE(environ_density) :: rhozero
        TYPE(environ_density) :: residual
        TYPE(environ_gradient) :: gradpoisson
        !
        CHARACTER(LEN=80) :: sub_name = 'generalized_fixedpoint'
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
                   rhoiter => dielectric%iterative, &
                   rhotot => dielectric%density, &
                   eps => dielectric%epsilon, &
                   maxiter => this%maxiter, &
                   tolrhoaux => this%tol)
            !
            !----------------------------------------------------------------------------
            ! Initialize local densities
            !
            CALL rhozero%init(cell)

            CALL residual%init(cell)
            !
            CALL gradpoisson%init(cell)
            !
            !----------------------------------------------------------------------------
            ! Set up auxiliary charge
            !
            total = charges%integrate()
            totpol = total * (1.D0 - dielectric%constant) / dielectric%constant
            rhozero%of_r = charges%of_r * (1.D0 - eps%of_r) / eps%of_r
            totzero = rhozero%integrate()
            totiter = rhoiter%integrate()
            !
            !----------------------------------------------------------------------------
            ! Write output table column headers
            !
            IF (io%verbosity >= 1 .AND. io%lnode) WRITE (io%debug_unit, 1001) totiter
            !
            IF (io%lnode) THEN
                !
                IF (io%verbosity >= 3) THEN
                    WRITE (io%debug_unit, 1002)
                ELSE IF (io%verbosity >= 1) THEN
                    WRITE (io%debug_unit, 1003)
                END IF
                !
            END IF
            !
            !----------------------------------------------------------------------------
            ! Start iterative algorithm
            !
            DO iter = 1, maxiter
                rhotot%of_r = charges%of_r + rhozero%of_r + rhoiter%of_r
                !
                CALL this%direct%grad_poisson(rhotot, gradpoisson, electrolyte, &
                                              semiconductor)
                !
                CALL dielectric%gradlog%scalar_product(gradpoisson, residual)
                !
                residual%of_r = residual%of_r / fpi / e2 - rhoiter%of_r
                rhoiter%of_r = rhoiter%of_r + this%mix * residual%of_r
                !
                delta_en = residual%euclidean_norm()
                delta_qm = residual%quadratic_mean()
                totiter = rhoiter%integrate()
                !
                !------------------------------------------------------------------------
                ! Print iteration results
                !
                IF (io%lnode) THEN
                    !
                    IF (io%verbosity >= 3) THEN
                        !
                        WRITE (io%debug_unit, 1004) &
                            iter, delta_qm, delta_en, tolrhoaux, totiter, totzero, &
                            totpol, total
                        !
                    ELSE IF (io%verbosity >= 1) THEN
                        WRITE (io%debug_unit, 1005) iter, delta_qm, delta_en, tolrhoaux
                    END IF
                    !
                END IF
                !
                !------------------------------------------------------------------------
                ! If residual is small enough exit
                !
                IF (delta_en < tolrhoaux .AND. iter > 0) THEN
                    !
                    IF (io%verbosity >= 1 .AND. io%lnode) WRITE (io%debug_unit, 1006)
                    !
                    EXIT
                    !
                ELSE IF (iter == maxiter) THEN
                    IF (io%lnode) WRITE (io%unit, 1007)
                END IF
                !
            END DO
            !
            IF (io%lstdout .AND. io%verbosity >= 1) WRITE (io%unit, 1008) delta_en, iter
            !
            !----------------------------------------------------------------------------
            ! Compute total electrostatic potential
            !
            rhotot%of_r = charges%of_r + rhozero%of_r + rhoiter%of_r
            !
            CALL this%direct%poisson(rhotot, v, electrolyte, semiconductor)
            !
            rhotot%of_r = rhozero%of_r + rhoiter%of_r
            ! in rhotot store total polarization charge
            !
            !----------------------------------------------------------------------------
            ! Clean up local densities
            !
            CALL rhozero%destroy()
            !
            CALL residual%destroy()
            !
            CALL gradpoisson%destroy()
            !
        END ASSOCIATE
        !
        !--------------------------------------------------------------------------------
        !
1000    FORMAT(/, 4('%'), ' COMPUTE ELECTROSTATIC POTENTIAL ', 43('%'),/)
        !
1001    FORMAT(' Starting from polarization: rhoiter = ', F13.6,/)
        !
1002    FORMAT('   i |    delta_qm    |    delta_en    |       tol     | ' &
               '          total iterative polarization charge', /, 1X, 113('-'))
        !
1003    FORMAT('   i |    delta_qm    |    delta_en    |       tol', /, 1X, 54('-'))
        !
1004    FORMAT(1X, I3, 3(' | ', E14.6), ' | ', 4E14.6)
        !
1005    FORMAT(1X, I3, 3(' | ', E14.6))
        !
1006    FORMAT(/, ' Charges are converged, EXIT')
        !
1007    FORMAT(' Warning: Polarization charge not converged',/)
        !
1008    FORMAT('     Polarization accuracy =', 1PE8.1, ', # of iterations = ', i3)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE generalized_fixedpoint
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE pb_fixedpoint(this, v, charges, electrolyte, dielectric, inner)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(solver_fixedpoint), TARGET, INTENT(IN) :: this
        TYPE(environ_density), TARGET, INTENT(IN) :: charges
        TYPE(environ_electrolyte), TARGET, INTENT(IN) :: electrolyte
        CLASS(electrostatic_solver), INTENT(IN), OPTIONAL :: inner
        !
        TYPE(environ_density), TARGET, INTENT(INOUT) :: v
        TYPE(environ_dielectric), INTENT(INOUT), OPTIONAL :: dielectric
        !
        INTEGER :: iter, ityp, ir
        REAL(DP) :: total, totaux, delta_qm, delta_en, kT, factor, arg
        TYPE(environ_density) :: residual, rhotot, denominator, rhoaux, cfactor
        !
        REAL(DP), POINTER :: cbulk, z
        !
        REAL(DP), PARAMETER :: exp_arg_limit = 40.D0
        !
        CHARACTER(LEN=80) :: sub_name = 'pb_fixedpoint'
        !
        !--------------------------------------------------------------------------------
        !
        IF (io%verbosity >= 1 .AND. io%lnode) WRITE (io%debug_unit, 1100)
        !
        !--------------------------------------------------------------------------------
        !
        IF (.NOT. ASSOCIATED(charges%cell, electrolyte%gamma%cell)) &
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
            IF (.NOT. PRESENT(inner)) CALL io%error(sub_name, 'Missing inner solver', 1)
            !
        END IF
        !
        !--------------------------------------------------------------------------------
        !
        ASSOCIATE (cell => charges%cell, &
                   maxiter => this%maxiter, &
                   tolrhoaux => this%tol, &
                   base => electrolyte%base, &
                   cionmax => electrolyte%base%cionmax, &
                   gam => electrolyte%gamma)
            !
            !----------------------------------------------------------------------------
            ! Initialize local densities
            !
            CALL denominator%init(cell)
            !
            CALL rhoaux%init(cell)
            !
            CALL rhotot%init(cell)
            !
            CALL residual%init(cell)
            !
            CALL cfactor%init(cell)
            !
            !----------------------------------------------------------------------------
            !
            kT = K_BOLTZMANN_RY * base%temperature
            !
            v%of_r = 0.D0
            rhoaux%of_r = 0.D0
            !
            !----------------------------------------------------------------------------
            ! Write output table column headers
            !
            IF (io%lnode) THEN
                !
                IF (io%verbosity >= 3) THEN
                    WRITE (io%debug_unit, 1101)
                ELSE IF (io%verbosity >= 1) THEN
                    WRITE (io%debug_unit, 1102)
                END IF
                !
            END IF
            !
            !----------------------------------------------------------------------------
            ! Start iterative algorithm
            !
            DO iter = 1, maxiter
                !
                rhotot%of_r = charges%of_r + rhoaux%of_r
                !
                IF (PRESENT(dielectric)) THEN
                    CALL inner%generalized(rhotot, dielectric, v)
                ELSE
                    CALL this%direct%poisson(rhotot, v)
                END IF
                !
                !------------------------------------------------------------------------
                ! Calculate electrolyte charge
                !
                residual%of_r = 0.D0
                denominator%of_r = 1.D0
                !
                DO ityp = 1, base%ntyp
                    cbulk => base%ioncctype(ityp)%cbulk
                    z => base%ioncctype(ityp)%z
                    !
                    cfactor%of_r = 1.D0
                    !
                    DO ir = 1, cell%ir_end
                        arg = -z * v%of_r(ir) / kT
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
                    residual%of_r = residual%of_r + z * cbulk * cfactor%of_r
                    !
                    IF (cionmax > 0.D0) THEN
                        factor = cbulk / cionmax
                        !
                        SELECT CASE (base%electrolyte_entropy)
                            !
                        CASE ('full')
                            !
                            denominator%of_r = denominator%of_r - &
                                               factor * (1.D0 - cfactor%of_r)
                            !
                        CASE ('ions')
                            !
                            denominator%of_r = denominator%of_r - &
                                               factor * (1.D0 - gam%of_r * cfactor%of_r)
                            !
                        CASE DEFAULT
                            CALL io%error(sub_name, 'Unexpected electrolyte entropy', 1)
                            !
                        END SELECT
                        !
                    END IF
                    !
                    NULLIFY (z)
                    NULLIFY (cbulk)
                END DO
                !
                residual%of_r = gam%of_r * residual%of_r / denominator%of_r
                !
                !------------------------------------------------------------------------
                ! Residual is now the new electrolyte charge
                !
                residual%of_r = residual%of_r - rhoaux%of_r
                rhoaux%of_r = rhoaux%of_r + this%mix * residual%of_r
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
                        !
                        WRITE (io%debug_unit, 1103) &
                            iter, delta_qm, delta_en, tolrhoaux, totaux
                        !
                    ELSE IF (io%verbosity >= 1) THEN
                        WRITE (io%debug_unit, 1104) iter, delta_qm, delta_en, tolrhoaux
                    END IF
                    !
                END IF
                !
                !------------------------------------------------------------------------
                ! If residual is small enough exit
                !
                IF (delta_en < tolrhoaux .AND. iter > 0) THEN
                    !
                    IF (io%verbosity >= 1 .AND. io%lnode) WRITE (io%debug_unit, 1105)
                    !
                    EXIT
                    !
                ELSE IF (iter == maxiter) THEN
                    IF (io%lnode) WRITE (io%unit, 1106)
                END IF
                !
            END DO
            !
            IF (io%lstdout .AND. io%verbosity >= 1) WRITE (io%unit, 1107) delta_en, iter
            !
            !----------------------------------------------------------------------------
            ! Clean up local densities
            !
            CALL denominator%destroy()
            !
            CALL rhoaux%destroy()
            !
            CALL rhotot%destroy()
            !
            CALL residual%destroy()
            !
            CALL cfactor%destroy()
            !
        END ASSOCIATE
        !
        !--------------------------------------------------------------------------------
        !
1100    FORMAT(/, 4('%'), ' COMPUTE ELECTROSTATIC POTENTIAL ', 43('%'),/)
        !
1101    FORMAT('   i | outer delta_qm |    delta_en    |       tol      | ' &
               'total iterative electrolyte charge', /, 1X, 91('-'))
        !
1102    FORMAT('   i | outer delta_qm |    delta_en    |       tol', /, 1X, 54('-'))
        !
1103    FORMAT(1X, I3, 4(' | ', E14.6))
        !
1104    FORMAT(1X, I3, 3(' | ', E14.6))
        !
1105    FORMAT(/, ' Outer loop is converged, EXIT')
        !
1106    FORMAT(' Warning: Electrolyte charge not converged',/)
        !
1107    FORMAT('     Electrolyte accuracy =', 1PE8.1, ', # of iterations = ', i3)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE pb_fixedpoint
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE class_solver_fixedpoint
!----------------------------------------------------------------------------------------
