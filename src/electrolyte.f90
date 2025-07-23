!----------------------------------------------------------------------------------------
!
! Copyright (C) 2018-2022 ENVIRON (www.quantum-environ.org)
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
! Authors: Oliviero Andreussi (Department of Physics, UNT)
!          Francesco Nattino  (THEOS and NCCR-MARVEL, EPFL)
!          Nicola Marzari     (THEOS and NCCR-MARVEL, EPFL)
!          Edan Bainglass     (Department of Physics, UNT)
!
!----------------------------------------------------------------------------------------
!>
!! Module containing the main routines to handle environ_electrolyte
!! derived data types.
!!
!! Environ_electrolyte contains all the specifications and the details of
!! the electrolyte medium and of the ionic distribution in it.
!!
!----------------------------------------------------------------------------------------
MODULE class_electrolyte
    !------------------------------------------------------------------------------------
    !
    USE class_io, ONLY: io
    !
    USE environ_param, ONLY: DP, e2, BOHR_RADIUS_SI, AMU_SI, K_BOLTZMANN_RY, fpi
    !
    USE class_cell
    USE class_density
    !
    USE class_core_container
    !
    USE class_boundary
    USE class_boundary_electronic
    USE class_boundary_ionic
    USE class_boundary_system
    USE class_electrons
    USE class_electrolyte_base
    USE class_ions
    USE class_system
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
    TYPE, PUBLIC :: environ_electrolyte
        !--------------------------------------------------------------------------------
        !
        LOGICAL :: initialized = .FALSE.
        !
        !--------------------------------------------------------------------------------
        !
        LOGICAL :: lupdate = .FALSE.
        !
        TYPE(environ_electrolyte_base) :: base
        !
        !--------------------------------------------------------------------------------
        !
        CLASS(environ_boundary), ALLOCATABLE :: boundary
        TYPE(environ_density) :: density
        !
        !--------------------------------------------------------------------------------
        ! The electrolyte switch function and related quantities
        !
        TYPE(environ_density) :: gamma
        TYPE(environ_density) :: dgamma
        !
        TYPE(environ_density) :: de_dboundary_second_order
        REAL(DP) :: energy_second_order = 0.D0
        !
        REAL(DP) :: charge = 0.D0
        !
        !--------------------------------------------------------------------------------
    CONTAINS
        !--------------------------------------------------------------------------------
        !
        PROCEDURE, PRIVATE :: create => create_environ_electrolyte
        PROCEDURE :: init => init_environ_electrolyte
        PROCEDURE :: update => update_environ_electrolyte
        PROCEDURE :: destroy => destroy_environ_electrolyte
        !
        PROCEDURE :: of_boundary => electrolyte_of_boundary
        PROCEDURE :: of_potential => electrolyte_of_potential
        PROCEDURE :: energy => calc_eelectrolyte
        PROCEDURE :: de_dboundary => calc_deelectrolyte_dboundary
        !
        PROCEDURE :: printout => print_environ_electrolyte
        !
        !--------------------------------------------------------------------------------
    END TYPE environ_electrolyte
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
    SUBROUTINE create_environ_electrolyte(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_electrolyte), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: routine = 'create_environ_electrolyte'
        !
        !--------------------------------------------------------------------------------
        !
        this%lupdate = .FALSE.
        this%energy_second_order = 0.D0
        this%charge = 0.D0
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE create_environ_electrolyte
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_environ_electrolyte(this, ntyp, mode, rhomax, rhomin, &
                                        const, alpha, softness, distance, &
                                        spread, solvent_radius, radial_scale, &
                                        radial_spread, filling_threshold, &
                                        filling_spread, field_aware, field_factor, &
                                        field_asymmetry, field_max, field_min, &
                                        electrons, ions, system, temperature, cbulk, &
                                        cionmax, radius, z, electrolyte_entropy, &
                                        linearized, cores, deriv_method, cell)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        LOGICAL, INTENT(IN) :: linearized, field_aware
        INTEGER, INTENT(IN) :: ntyp
        CHARACTER(LEN=*), INTENT(IN) :: mode, electrolyte_entropy, deriv_method
        !
        REAL(DP), INTENT(IN) :: rhomax, rhomin, const, distance, spread, &
                                alpha, softness, temperature, solvent_radius, &
                                radial_scale, radial_spread, filling_threshold, &
                                filling_spread, field_factor, field_asymmetry, &
                                field_max, field_min, cionmax, radius
        !
        REAL(DP), DIMENSION(ntyp), INTENT(IN) :: cbulk, z
        !
        TYPE(environ_electrons), INTENT(IN) :: electrons
        TYPE(environ_ions), INTENT(IN) :: ions
        TYPE(environ_system), INTENT(IN) :: system
        TYPE(core_container), INTENT(IN) :: cores
        TYPE(environ_cell), INTENT(IN) :: cell
        !
        CLASS(environ_electrolyte), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: routine = 'init_environ_electrolyte'
        !
        !--------------------------------------------------------------------------------
        !
        CALL this%create()
        !
        CALL this%base%init(ntyp, const, distance, spread, temperature, cbulk, cionmax, &
                            radius, z, electrolyte_entropy, linearized, cell)
        !
        !--------------------------------------------------------------------------------
        ! Casting and general setup
        !
        SELECT CASE (mode)
            !
        CASE ('electronic', 'full')
            ALLOCATE (environ_boundary_electronic :: this%boundary)
            !
        CASE ('ionic')
            ALLOCATE (environ_boundary_ionic :: this%boundary)
            !
        CASE ('system')
            ALLOCATE (environ_boundary_system :: this%boundary)
            !
        CASE DEFAULT
            CALL io%error(routine, "Unrecognized boundary mode", 1)
            !
        END SELECT
        !
        CALL this%boundary%pre_init(mode, .TRUE., .TRUE., .FALSE., cores, deriv_method, &
                                    cell, 'electrolyte')
        !
        !--------------------------------------------------------------------------------
        ! Boundary awareness
        !
        IF (solvent_radius > 0.D0) &
            CALL this%boundary%init_solvent_aware(solvent_radius, radial_scale, &
                                                  radial_spread, filling_threshold, &
                                                  filling_spread)
        !
        IF (field_aware) &
            CALL this%boundary%init_field_aware(field_factor, field_asymmetry, &
                                                field_max, field_min)
        !
        !--------------------------------------------------------------------------------
        ! Specific setup
        !
        SELECT TYPE (boundary => this%boundary)
            !
        TYPE IS (environ_boundary_electronic)
            CALL boundary%init(rhomax, rhomin, electrons, ions)
            !
        TYPE IS (environ_boundary_ionic)
            CALL boundary%init(alpha, softness, ions, electrons)
            !
        TYPE IS (environ_boundary_system)
            CALL boundary%init(distance, spread, system)
            !
        CLASS DEFAULT
            CALL io%error(routine, "Unrecognized boundary mode", 1)
            !
        END SELECT
        !
        !--------------------------------------------------------------------------------
        ! Densities
        !
        CALL this%density%init(cell, 'electrolyte')
        !
        CALL this%gamma%init(cell, 'gamma')
        !
        CALL this%dgamma%init(cell, 'dgamma')
        !
        IF (this%base%linearized) CALL this%de_dboundary_second_order%init(cell)
        !
        this%initialized = .TRUE.
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_environ_electrolyte
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE update_environ_electrolyte(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_electrolyte), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: routine = 'update_environ_electrolyte'
        !
        !--------------------------------------------------------------------------------
        !
        CALL env_start_clock(routine)
        !
        !--------------------------------------------------------------------------------
        ! Check if the boundary is under update (status = 1) or
        ! has been fully updated (status = 2)
        !
        IF (this%boundary%update_status > 0) this%lupdate = .TRUE.
        !
        IF (this%lupdate) THEN
            !
            !----------------------------------------------------------------------------
            ! Update the electrolyte in space if the boundary is ready
            !
            IF (this%boundary%update_status == 2) THEN
                !
                CALL this%of_boundary()
                !
                this%lupdate = .FALSE.
                !
            END IF
            !
        END IF
        !
        CALL env_stop_clock(routine)
        !
        !--------------------------------------------------------------------------------
        ! Output current state
        !
        IF (.NOT. this%lupdate) CALL this%printout()
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE update_environ_electrolyte
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE destroy_environ_electrolyte(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_electrolyte), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: routine = 'destroy_environ_electrolyte'
        !
        !--------------------------------------------------------------------------------
        !
        CALL this%boundary%destroy()
        !
        CALL this%density%destroy()
        !
        CALL this%gamma%destroy()
        !
        CALL this%dgamma%destroy()
        !
        IF (this%base%linearized) CALL this%de_dboundary_second_order%destroy()
        !
        CALL this%base%destroy()
        !
        this%initialized = .FALSE.
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE destroy_environ_electrolyte
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !                                  GENERAL METHODS
    !
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE electrolyte_of_boundary(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_electrolyte), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: routine = 'electrolyte_of_boundary'
        !
        !--------------------------------------------------------------------------------
        ! Compute exclusion function gamma(r) and dgamma/ds(r)
        !
        this%gamma%of_r = 1.D0 - this%boundary%scaled%of_r
        this%dgamma%of_r = -1.D0
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE electrolyte_of_boundary
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE electrolyte_of_potential(this, potential)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_density), INTENT(IN) :: potential
        !
        CLASS(environ_electrolyte), INTENT(INOUT) :: this
        !
        INTEGER :: i, j
        REAL(DP) :: fact, sumcbulk, arg
        TYPE(environ_density) :: denominator
        !
        REAL(DP), PARAMETER :: exp_arg_limit = 40.D0
        !
        CHARACTER(LEN=80) :: routine = 'calc_electrolyte_density'
        !
        !--------------------------------------------------------------------------------
        !
        ASSOCIATE (cell => potential%cell, &
                   gamma => this%gamma%of_r, &
                   dgamma => this%dgamma%of_r, &
                   pot => potential%of_r, &
                   rho => this%density, &
                   base => this%base, &
                   k2 => this%base%k2, &
                   cionmax => this%base%cionmax, &
                   de_dboundary_second_order => this%de_dboundary_second_order%of_r, &
                   kT => K_BOLTZMANN_RY * this%base%temperature)
            !
            !----------------------------------------------------------------------------
            !
            rho%of_r = 0.D0
            !
            CALL denominator%init(cell)
            !
            denominator%of_r = 1.D0
            !
            DO i = 1, base%ntyp
                !
                ASSOCIATE (cfactor => base%ioncctype(i)%cfactor%of_r, &
                           cbulk => base%ioncctype(i)%cbulk, &
                           z => base%ioncctype(i)%z)
                    !
                    cfactor = 1.D0
                    !
                    IF (base%linearized) THEN ! Linearized PB and modified PB
                        cfactor = 1.D0 - z * pot / kT
                        !
                        IF (cionmax > 0.D0) THEN
                            fact = cbulk / cionmax
                            !
                            IF (base%electrolyte_entropy == 'ions') &
                                denominator%of_r = denominator%of_r - &
                                                   fact * (1.D0 - gamma)
                            !
                        END IF
                        !
                    ELSE ! full PB
                        !
                        DO j = 1, cell%ir_end
                            !
                            !------------------------------------------------------------
                            ! Numerical problems arise when computing exp( -z*pot/kT )
                            ! in regions close to the nuclei (exponent is too large).
                            !
                            arg = -z * pot(j) / kT
                            !
                            IF (arg > exp_arg_limit) THEN
                                cfactor(j) = EXP(exp_arg_limit)
                            ELSE IF (arg < -exp_arg_limit) THEN
                                cfactor(j) = EXP(-exp_arg_limit)
                            ELSE
                                cfactor(j) = EXP(arg)
                            END IF
                            !
                        END DO
                        !
                        IF (cionmax /= 0.D0) THEN ! full modified PB
                            fact = cbulk / cionmax
                            !
                            SELECT CASE (base%electrolyte_entropy)
                                !
                            CASE ('full')
                                denominator%of_r = denominator%of_r - &
                                                   fact * (1.D0 - cfactor)
                                !
                            CASE ('ions')
                                !
                                denominator%of_r = denominator%of_r - &
                                                   fact * (1.D0 - gamma * cfactor)
                                !
                            CASE DEFAULT
                                CALL io%error(routine, "Unexpected electrolyte entropy", 1)
                                !
                            END SELECT
                            !
                        END IF
                        !
                    END IF
                    !
                END ASSOCIATE
                !
            END DO
            !
            DO i = 1, base%ntyp
                !
                ASSOCIATE (c => base%ioncctype(i)%c%of_r, &
                           cbulk => base%ioncctype(i)%cbulk, &
                           cfactor => base%ioncctype(i)%cfactor%of_r, &
                           z => base%ioncctype(i)%z)
                    !
                    c = gamma * cbulk * cfactor / denominator%of_r
                    rho%of_r = rho%of_r + c * z
                    !
                END ASSOCIATE
                !
            END DO
            !
            this%charge = rho%integrate()
            !
            IF (base%linearized) THEN
                !
                !------------------------------------------------------------------------
                ! Compute energy and de_dboundary terms that depend on the potential.
                ! These are the second order terms, first order terms are equal to zero.
                !
                ! NOTE: the energy is equal to minus the electrostatic interaction of the
                !       electrolyte, and is identical for the various implementations of
                !       the entropy.
                !
                this%energy_second_order = 0.5D0 * rho%scalar_product(potential)
                !
                IF (base%electrolyte_entropy == 'ions' .AND. cionmax > 0.D0) THEN
                    sumcbulk = SUM(base%ioncctype%cbulk)
                    !
                    de_dboundary_second_order = &
                        -0.5D0 * k2 / e2 / fpi * (1.D0 - sumcbulk / cionmax) * &
                        (pot / denominator%of_r)**2 * dgamma
                    !
                ELSE
                    de_dboundary_second_order = -0.5D0 * k2 / e2 / fpi * dgamma * pot**2
                END IF
                !
            END IF
            !
            CALL denominator%destroy()
            !
        END ASSOCIATE
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE electrolyte_of_potential
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_eelectrolyte(this, energy)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_electrolyte), INTENT(IN) :: this
        !
        REAL(DP), INTENT(OUT) :: energy
        !
        INTEGER :: i
        REAL(DP) :: logterm, integral
        !
        TYPE(environ_density) :: arg, f
        !
        CHARACTER(LEN=80) :: routine = 'calc_eelectrolyte'
        !
        !--------------------------------------------------------------------------------
        !
        ASSOCIATE (cell => this%gamma%cell, &
                   omega => this%gamma%cell%omega, &
                   gamma => this%gamma%of_r, &
                   base => this%base, &
                   cionmax => this%base%cionmax, &
                   kT => K_BOLTZMANN_RY * this%base%temperature, &
                   sumcbulk => SUM(this%base%ioncctype%cbulk))
            !
            !----------------------------------------------------------------------------
            !
            energy = 0.D0
            !
            CALL arg%init(cell)
            !
            CALL f%init(cell)
            !
            IF (base%linearized) THEN
                !
                IF (cionmax == 0.D0) THEN ! Linearized PB
                    integral = this%gamma%integrate()
                    energy = kT * sumcbulk * (omega - integral)
                ELSE ! Linearized modified PB
                    !
                    SELECT CASE (base%electrolyte_entropy)
                        !
                    CASE ('full')
                        integral = this%gamma%integrate()
                        logterm = LOG(1.D0 - sumcbulk / cionmax)
                        energy = kT * cionmax * logterm * (integral - omega)
                        !
                    CASE ('ions')
                        arg%of_r = 1.D0 - sumcbulk / cionmax * (1.D0 - gamma)
                        f%of_r = LOG(arg%of_r)
                        integral = f%integrate()
                        energy = -kT * cionmax * integral
                        !
                    CASE DEFAULT
                        CALL io%error(routine, "Unexpected electrolyte entropy", 1)
                        !
                    END SELECT
                    !
                END IF
                !
                energy = energy + this%energy_second_order
            ELSE
                arg%of_r = 0.D0
                !
                DO i = 1, base%ntyp
                    !
                    arg%of_r = arg%of_r + base%ioncctype(i)%cfactor%of_r * &
                               base%ioncctype(i)%cbulk
                    !
                END DO
                !
                IF (cionmax == 0.D0) THEN ! Full PB
                    f%of_r = gamma * arg%of_r - sumcbulk
                    integral = f%integrate()
                    energy = -kT * integral
                ELSE ! full modified PB
                    !
                    SELECT CASE (base%electrolyte_entropy)
                        !
                    CASE ('full')
                        arg%of_r = arg%of_r / (cionmax - sumcbulk)
                        arg%of_r = arg%of_r + 1.D0
                        f%of_r = gamma * LOG(arg%of_r)
                        integral = f%integrate()
                        logterm = LOG(1.D0 - sumcbulk / cionmax)
                        energy = -kT * cionmax * (integral + logterm * omega)
                        !
                    CASE ('ions')
                        arg%of_r = arg%of_r / cionmax * gamma
                        arg%of_r = arg%of_r + 1.D0 - sumcbulk / cionmax
                        f%of_r = LOG(arg%of_r)
                        integral = f%integrate()
                        energy = -kT * cionmax * integral
                        !
                    CASE DEFAULT
                        CALL io%error(routine, "Unexpected electrolyte entropy", 1)
                        !
                    END SELECT
                    !
                END IF
                !
            END IF
            !
            CALL arg%destroy()
            !
            CALL f%destroy()
            !
        END ASSOCIATE
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_eelectrolyte
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_deelectrolyte_dboundary(this, de_dboundary)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_electrolyte), INTENT(IN) :: this
        !
        TYPE(environ_density), INTENT(INOUT) :: de_dboundary
        !
        INTEGER :: i
        !
        TYPE(environ_density) :: arg
        !
        CHARACTER(LEN=80) :: routine = 'calc_deelectrolyte_dboundary'
        !
        !--------------------------------------------------------------------------------
        !
        ASSOCIATE (gamma => this%gamma%of_r, &
                   dgamma => this%dgamma%of_r, &
                   base => this%base, &
                   cionmax => this%base%cionmax, &
                   kT => K_BOLTZMANN_RY * this%base%temperature, &
                   sumcbulk => SUM(this%base%ioncctype%cbulk))
            !
            !----------------------------------------------------------------------------
            !
            IF (base%linearized) THEN
                !
                IF (cionmax == 0.D0) THEN ! linearized PB
                    de_dboundary%of_r = de_dboundary%of_r - dgamma * kT * sumcbulk
                ELSE ! linearized modified PB
                    !
                    SELECT CASE (base%electrolyte_entropy)
                        !
                    CASE ('full')
                        !
                        de_dboundary%of_r = de_dboundary%of_r + &
                                            dgamma * kT * cionmax * &
                                            LOG(1.D0 - sumcbulk / cionmax)
                        !
                    CASE ('ions')
                        !
                        de_dboundary%of_r = de_dboundary%of_r - &
                                            dgamma * kT * sumcbulk / &
                                            (1.D0 - sumcbulk / cionmax * (1.D0 - gamma))
                        !
                    CASE DEFAULT
                        CALL io%error(routine, "Unexpected electrolyte entropy", 1)
                        !
                    END SELECT
                    !
                END IF
                !
                de_dboundary%of_r = de_dboundary%of_r + &
                                    this%de_dboundary_second_order%of_r
                !
            ELSE
                !
                CALL arg%init(de_dboundary%cell)
                !
                arg%of_r = 0.D0
                !
                DO i = 1, base%ntyp
                    !
                    arg%of_r = arg%of_r + &
                               base%ioncctype(i)%cfactor%of_r * base%ioncctype(i)%cbulk
                    !
                END DO
                !
                IF (cionmax == 0.D0) THEN ! full PB
                    de_dboundary%of_r = de_dboundary%of_r - dgamma * kT * arg%of_r
                ELSE ! full modified PB
                    !
                    SELECT CASE (base%electrolyte_entropy)
                        !
                    CASE ('full')
                        arg%of_r = arg%of_r / (cionmax - sumcbulk)
                        arg%of_r = arg%of_r + 1.D0
                        !
                        de_dboundary%of_r = de_dboundary%of_r - &
                                            dgamma * kT * cionmax * LOG(arg%of_r)
                        !
                    CASE ('ions')
                        !
                        de_dboundary%of_r = &
                            de_dboundary%of_r - &
                            dgamma * kT * arg%of_r / &
                            (1.D0 - (sumcbulk - arg%of_r * gamma) / cionmax)
                        !
                    CASE DEFAULT
                        CALL io%error(routine, "Unexpected electrolyte entropy", 1)
                        !
                    END SELECT
                    !
                END IF
                !
                CALL arg%destroy()
                !
            END IF
            !
        END ASSOCIATE
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_deelectrolyte_dboundary
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !                                   OUTPUT METHODS
    !
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !>
    !! Prints the details of the electrolyte
    !!
    !! Nested objects receive a decremented passed verbose to trigger block printing
    !!
    !! @param verbose       : (INTEGER) adds verbosity to global verbose
    !! @param debug_verbose : (INTEGER) replaces global verbose for debugging
    !! @param unit          : (INTEGER) output target (default = io%debug_unit)
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE print_environ_electrolyte(this, verbose, debug_verbose, unit)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_electrolyte), TARGET, INTENT(IN) :: this
        INTEGER, OPTIONAL, INTENT(IN) :: verbose, debug_verbose, unit
        !
        INTEGER :: i
        INTEGER :: base_verbose, local_verbose, passed_verbose, local_unit
        !
        TYPE(environ_electrolyte_base), POINTER :: base
        !
        CHARACTER(LEN=80) :: routine = 'print_environ_electrolyte'
        !
        !--------------------------------------------------------------------------------
        !
        IF (PRESENT(debug_verbose)) THEN
            base_verbose = debug_verbose
            !
            IF (PRESENT(verbose)) THEN
                local_verbose = verbose
            ELSE
                local_verbose = debug_verbose
            END IF
            !
            passed_verbose = verbose - 1
            !
        ELSE IF (io%verbosity > 0) THEN
            base_verbose = io%verbosity
            !
            IF (PRESENT(verbose)) THEN
                local_verbose = base_verbose + verbose
            ELSE
                local_verbose = base_verbose
            END IF
            !
            passed_verbose = local_verbose - base_verbose - 1
            !
        ELSE
            RETURN
        END IF
        !
        IF (PRESENT(unit)) THEN
            local_unit = unit
        ELSE
            local_unit = io%debug_unit
        END IF
        !
        base => this%base
        !
        IF (local_verbose >= 1) THEN
            !
            IF (io%lnode) THEN
                WRITE (local_unit, 1000)
                !
                WRITE (local_unit, 1001) &
                    base%ntyp, base%temperature, 1.D0 / SQRT(base%k2)
                !
                IF (base%cionmax > 0.D0) WRITE (local_unit, 1002) base%cionmax
                !
                WRITE (local_unit, 1003) base%linearized
                WRITE (local_unit, 1004) this%charge
            END IF
            !
            IF (local_verbose >= 3) THEN
                !
                CALL this%density%printout(passed_verbose, debug_verbose, local_unit)
                !
                CALL this%gamma%printout(passed_verbose, debug_verbose, local_unit)
                !
            END IF
            !
            IF (local_verbose >= 5) &
                CALL this%dgamma%printout(passed_verbose, debug_verbose, local_unit)
            !
            IF (io%lnode) THEN
                WRITE (local_unit, 1005)
                WRITE (local_unit, 1006) ! header
            END IF
            !
            DO i = 1, base%ntyp
                !
                IF (io%lnode) &
                    WRITE (local_unit, 1007) &
                    base%ioncctype(i)%index, base%ioncctype(i)%cbulk, &
                    base%ioncctype(i)%cbulk * AMU_SI / BOHR_RADIUS_SI**3, &
                    base%ioncctype(i)%z
                !
            END DO
            !
            IF (local_verbose >= 5) THEN
                !
                DO i = 1, base%ntyp
                    !
                    CALL base%ioncctype(i)%c%printout(passed_verbose, debug_verbose, &
                                                      local_unit)
                    !
                    CALL base%ioncctype(i)%cfactor%printout(passed_verbose, &
                                                            debug_verbose, local_unit)
                    !
                END DO
                !
            END IF
            !
        END IF
        !
        FLUSH (local_unit)
        !
        !--------------------------------------------------------------------------------
        !
1000    FORMAT(/, 4('%'), " ELECTROLYTE ", 64('%'))
        !
1001    FORMAT(/, " number electrolyte species = ", I14, /, &
                " solvent temperature        = ", F14.1, /, &
                " Debye length / sqrt(eps)   = ", F14.7)
        !
1002    FORMAT(/, " modified Poisson-Boltzmann:", /, &
                " maximum concentration      = ", F14.7)
        !
1003    FORMAT(/, " electrolyte flags:", /, &
                " linearized                 = ", L14)
        !
1004    FORMAT(/, " total electrolyte charge   = ", F14.7)
        !
1005    FORMAT(/, " species", /, 1X, 7('='))
        !
1006    FORMAT(/, "   i |  c_bulk (a.u.) | c_bulk (mol/L) | ionic charge", /, 1X, 52('-'))
        !
1007    FORMAT(1X, I3, " | ", E14.4, " | ", F14.7, " | ", F12.2)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE print_environ_electrolyte
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE class_electrolyte
!----------------------------------------------------------------------------------------
