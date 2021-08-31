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
!          Francesco Nattino  (THEOS and NCCR-MARVEL, EPFL)
!          Nicola Marzari     (THEOS and NCCR-MARVEL, EPFL)
!          Edan Bainglass     (Department of Physics, UNT)
!
!----------------------------------------------------------------------------------------
!> Module containing the main routines to handle environ_electrolyte
!! derived data types.
!!
!! Environ_electrolyte contains all the specifications and the details of
!! the electrolyte medium and of the ionic distribution in it.
!!
!----------------------------------------------------------------------------------------
MODULE class_electrolyte
    !------------------------------------------------------------------------------------
    !
    USE env_base_io, ONLY: ionode, environ_unit, verbose, depth
    !
    USE environ_param, ONLY: DP, e2, BOHR_RADIUS_SI, AMU_SI, K_BOLTZMANN_RY, fpi
    !
    USE class_cell
    USE class_density
    USE class_functions
    !
    USE class_core_container_derivatives
    !
    USE class_boundary
    USE class_electrons
    USE class_ions
    USE class_ioncctype
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
        LOGICAL :: lupdate
        CHARACTER(LEN=80) :: electrolyte_entropy
        LOGICAL :: linearized
        INTEGER :: ntyp
        TYPE(environ_ioncctype), ALLOCATABLE :: ioncctype(:)
        !
        REAL(DP) :: temperature
        REAL(DP) :: k2
        REAL(DP) :: cionmax
        REAL(DP) :: permittivity
        !
        TYPE(environ_boundary) :: boundary
        TYPE(environ_density) :: density
        !
        !--------------------------------------------------------------------------------
        ! The electrolyte switch function and related quantities
        !
        TYPE(environ_density) :: gamma
        TYPE(environ_density) :: dgamma
        !
        TYPE(environ_density) :: de_dboundary_second_order
        REAL(DP) :: energy_second_order
        REAL(DP) :: charge
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
        CHARACTER(LEN=80) :: sub_name = 'create_environ_electrolyte'
        !
        !--------------------------------------------------------------------------------
        !
        IF (ALLOCATED(this%ioncctype)) &
            CALL env_errore(sub_name, 'Trying to create an existing object', 1)
        !
        IF (ALLOCATED(this%density%of_r)) &
            CALL env_errore(sub_name, 'Trying to create an existing object', 1)
        !
        IF (ALLOCATED(this%gamma%of_r)) &
            CALL env_errore(sub_name, 'Trying to create an existing object', 1)
        !
        IF (ALLOCATED(this%dgamma%of_r)) &
            CALL env_errore(sub_name, 'Trying to create an existing object', 1)
        !
        IF (ALLOCATED(this%de_dboundary_second_order%of_r)) &
            CALL env_errore(sub_name, 'Trying to create an existing object', 1)
        !
        !--------------------------------------------------------------------------------
        !
        this%lupdate = .FALSE.
        this%linearized = .FALSE.
        !
        this%charge = 0.D0
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE create_environ_electrolyte
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_environ_electrolyte(this, ntyp, mode, stype, rhomax, rhomin, &
                                        tbeta, const, alpha, softness, distance, &
                                        spread, solvent_radius, radial_scale, &
                                        radial_spread, filling_threshold, &
                                        filling_spread, field_awareness, &
                                        charge_asymmetry, field_max, field_min, &
                                        electrons, ions, system, derivatives, &
                                        temperature, cbulk, cionmax, radius, z, &
                                        electrolyte_entropy, linearized, cell)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        LOGICAL, INTENT(IN) :: linearized
        INTEGER, INTENT(IN) :: ntyp, stype
        CHARACTER(LEN=80), INTENT(IN) :: mode, electrolyte_entropy
        !
        REAL(DP), INTENT(IN) :: rhomax, rhomin, tbeta, const, distance, spread, &
                                alpha, softness, temperature, solvent_radius, &
                                radial_scale, radial_spread, filling_threshold, &
                                filling_spread, field_awareness, charge_asymmetry, &
                                field_max, field_min, cionmax, radius
        !
        REAL(DP), DIMENSION(ntyp), INTENT(IN) :: cbulk, z
        !
        TYPE(environ_electrons), INTENT(IN) :: electrons
        TYPE(environ_ions), INTENT(IN) :: ions
        TYPE(environ_system), TARGET, INTENT(IN) :: system
        TYPE(container_derivatives), TARGET, INTENT(IN) :: derivatives
        TYPE(environ_cell), INTENT(IN) :: cell
        !
        CLASS(environ_electrolyte), INTENT(INOUT) :: this
        !
        INTEGER :: ityp
        REAL(DP) :: neutral, sumcbulk, sum_cz2, arg, kT, e
        !
        CHARACTER(LEN=80) :: ityps, local_label
        !
        CHARACTER(LEN=80) :: sub_name = 'init_environ_electrolyte'
        !
        !--------------------------------------------------------------------------------
        !
        CALL this%create()
        !
        local_label = 'electrolyte'
        !
        CALL this%boundary%init(.TRUE., .TRUE., .FALSE., mode, stype, rhomax, &
                                rhomin, tbeta, const, alpha, softness, distance, &
                                spread, solvent_radius, radial_scale, radial_spread, &
                                filling_threshold, filling_spread, field_awareness, &
                                charge_asymmetry, field_max, field_min, electrons, &
                                ions, system, derivatives, cell, local_label)
        !
        !--------------------------------------------------------------------------------
        ! Setup all electrolyte parameters (with checks)
        !
        this%linearized = linearized
        this%ntyp = ntyp
        this%electrolyte_entropy = TRIM(electrolyte_entropy)
        this%temperature = temperature
        this%cionmax = 0.D0
        this%permittivity = const
        !
        ALLOCATE (this%ioncctype(ntyp))
        !
        neutral = 0.D0
        sum_cz2 = 0.D0
        !
        DO ityp = 1, ntyp
            !
            CALL this%ioncctype(ityp)%init(ityp, cbulk(ityp), z(ityp), cell)
            !
            neutral = neutral + cbulk(ityp) * z(ityp)
            sum_cz2 = sum_cz2 + this%ioncctype(ityp)%cbulk * this%ioncctype(ityp)%z**2
        END DO
        !
        IF (neutral > 1.D-8) &
            CALL env_errore(sub_name, 'Bulk electrolyte is not neutral', 1)
        !
        kT = K_BOLTZMANN_RY * temperature
        !
        this%k2 = sum_cz2 / kT * e2 * fpi ! k^2 = eps / lambda_D^2
        !
        !--------------------------------------------------------------------------------
        ! If cionmax is not provided in input but radius is, calculate cionmax
        !
        this%cionmax = cionmax * BOHR_RADIUS_SI**3 / AMU_SI
        !
        IF (cionmax == 0.D0 .AND. radius > 0.D0) &
            this%cionmax = 0.64D0 * 3.D0 / fpi / radius**3
        !
        !--------------------------------------------------------------------------------
        ! Check suitability of cionmax value
        !
        sumcbulk = SUM(this%ioncctype(:)%cbulk)
        !
        IF (this%cionmax > 0.D0 .AND. this%cionmax <= sumcbulk) &
            CALL env_errore(sub_name, &
                            'cionmax should be larger than the sum of cbulks', 1)
        !
        this%energy_second_order = 0.D0
        !
        !--------------------------------------------------------------------------------
        ! Densities
        !
        local_label = 'electrolyte'
        !
        CALL this%density%init(cell, local_label)
        !
        local_label = 'gamma'
        !
        CALL this%gamma%init(cell, local_label)
        !
        local_label = 'dgamma'
        !
        CALL this%dgamma%init(cell, local_label)
        !
        IF (this%linearized) CALL this%de_dboundary_second_order%init(cell)
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
        INTEGER :: ityp
        !
        CHARACTER(LEN=80) :: sub_name = 'update_environ_electrolyte'
        !
        !--------------------------------------------------------------------------------
        !
        CALL env_start_clock(sub_name)
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
        CALL env_stop_clock(sub_name)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE update_environ_electrolyte
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE destroy_environ_electrolyte(this, lflag)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        LOGICAL, INTENT(IN) :: lflag
        !
        CLASS(environ_electrolyte), INTENT(INOUT) :: this
        !
        INTEGER :: ityp
        !
        CHARACTER(LEN=80) :: sub_name = 'destroy_environ_electrolyte'
        !
        !--------------------------------------------------------------------------------
        !
        CALL this%boundary%destroy(lflag)
        !
        DO ityp = 1, this%ntyp
            CALL this%ioncctype(ityp)%destroy()
        END DO
        !
        CALL this%gamma%destroy()
        !
        CALL this%dgamma%destroy()
        !
        CALL this%density%destroy()
        !
        IF (this%linearized) CALL this%de_dboundary_second_order%destroy()
        !
        IF (lflag) THEN
            !
            IF (.NOT. ALLOCATED(this%ioncctype)) &
                CALL env_errore(sub_name, 'Trying to destroy an empty object', 1)
            !
            DEALLOCATE (this%ioncctype)
        END IF
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
        CLASS(environ_electrolyte), TARGET, INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: sub_name = 'electrolyte_of_boundary'
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
        TYPE(environ_density), TARGET, INTENT(IN) :: potential
        !
        CLASS(environ_electrolyte), TARGET, INTENT(INOUT) :: this
        !
        REAL(DP), POINTER :: z, cbulk
        REAL(DP), DIMENSION(:), POINTER :: pot, rho, c, cfactor, gam
        !
        REAL(DP) :: kT, e, factor, sumcbulk, arg
        INTEGER :: ityp, ir
        TYPE(environ_density) :: denominator
        !
        REAL(DP), PARAMETER :: exp_arg_limit = 40.D0
        !
        CHARACTER(LEN=80) :: sub_name = 'calc_electrolyte_density'
        !
        !--------------------------------------------------------------------------------
        !
        gam => this%gamma%of_r
        pot => potential%of_r
        rho => this%density%of_r
        !
        rho = 0.D0
        kT = K_BOLTZMANN_RY * this%temperature
        !
        CALL denominator%init(potential%cell)
        !
        denominator%of_r = 1.D0
        !
        DO ityp = 1, this%ntyp
            cfactor => this%ioncctype(ityp)%cfactor%of_r
            cbulk => this%ioncctype(ityp)%cbulk
            z => this%ioncctype(ityp)%z
            cfactor = 1.D0
            !
            IF (this%linearized) THEN ! Linearized PB and modified PB
                !
                cfactor = 1.D0 - z * pot / kT
                !
                IF (this%cionmax > 0.D0) THEN
                    factor = cbulk / this%cionmax
                    !
                    IF (this%electrolyte_entropy == 'ions') &
                        denominator%of_r = denominator%of_r - factor * (1.D0 - gam)
                    !
                END IF
                !
            ELSE ! full PB
                !
                DO ir = 1, potential%cell%ir_end
                    !
                    !--------------------------------------------------------------------
                    ! Numerical problems arise when computing exp( -z*pot/kT )
                    ! in regions close to the nuclei (exponent is too large).
                    !
                    arg = -z * pot(ir) / kT
                    !
                    IF (arg > exp_arg_limit) THEN
                        cfactor(ir) = EXP(exp_arg_limit)
                    ELSE IF (arg < -exp_arg_limit) THEN
                        cfactor(ir) = EXP(-exp_arg_limit)
                    ELSE
                        cfactor(ir) = EXP(arg)
                    END IF
                    !
                END DO
                !
                IF (this%cionmax /= 0.D0) THEN ! full modified PB
                    !
                    factor = cbulk / this%cionmax
                    !
                    SELECT CASE (this%electrolyte_entropy)
                        !
                    CASE ('full')
                        denominator%of_r = denominator%of_r - factor * (1.D0 - cfactor)
                        !
                    CASE ('ions')
                        !
                        denominator%of_r = denominator%of_r - &
                                           factor * (1.D0 - gam * cfactor)
                        !
                    END SELECT
                    !
                END IF
                !
            END IF
            !
            NULLIFY (cfactor)
            NULLIFY (cbulk)
            NULLIFY (z)
            !
        END DO
        !
        DO ityp = 1, this%ntyp
            c => this%ioncctype(ityp)%c%of_r
            cfactor => this%ioncctype(ityp)%cfactor%of_r
            cbulk => this%ioncctype(ityp)%cbulk
            z => this%ioncctype(ityp)%z
            c = gam * cbulk * cfactor / denominator%of_r
            rho = rho + c * z
            !
            NULLIFY (c)
            NULLIFY (cfactor)
            NULLIFY (cbulk)
            NULLIFY (z)
        END DO
        !
        this%charge = this%density%integrate()
        !
        IF (this%linearized) THEN
            !
            !----------------------------------------------------------------------------
            ! Compute energy and de_dboundary terms that depend on the potential.
            ! These are the second order terms, first order terms are equal to zero.
            ! NOTE: the energy is equal to minus the electrostatic interaction of the
            ! electrolyte, and is identical for the various implementations of the entropy.
            !
            this%energy_second_order = 0.5D0 * this%density%scalar_product(potential)
            !
            IF (this%electrolyte_entropy == 'ions' .AND. &
                this%cionmax > 0.D0) THEN
                !
                sumcbulk = SUM(this%ioncctype(:)%cbulk)
                !
                this%de_dboundary_second_order%of_r = &
                    -0.5D0 * this%k2 / e2 / fpi * &
                    (1.D0 - sumcbulk / this%cionmax) * &
                    (pot / denominator%of_r)**2 * this%dgamma%of_r
                !
            ELSE
                !
                this%de_dboundary_second_order%of_r = -0.5D0 * this%k2 / e2 / fpi * &
                                                      pot * pot * this%dgamma%of_r
                !
            END IF
            !
        END IF
        !
        CALL denominator%destroy()
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
        REAL(DP) :: kT, sumcbulk, logterm, integral
        INTEGER :: ityp
        TYPE(environ_density) :: arg, f
        !
        !--------------------------------------------------------------------------------
        !
        energy = 0.D0
        !
        kT = K_BOLTZMANN_RY * this%temperature
        sumcbulk = SUM(this%ioncctype(:)%cbulk)
        !
        CALL arg%init(this%gamma%cell)
        !
        CALL f%init(this%gamma%cell)
        !
        IF (this%linearized) THEN
            !
            IF (this%cionmax == 0.D0) THEN ! Linearized PB
                !
                integral = this%gamma%integrate()
                energy = kT * sumcbulk * (this%gamma%cell%omega - integral)
                !
            ELSE ! Linearized modified PB
                !
                SELECT CASE (this%electrolyte_entropy)
                    !
                CASE ('full')
                    integral = this%gamma%integrate()
                    logterm = LOG(1.D0 - sumcbulk / this%cionmax)
                    !
                    energy = kT * this%cionmax * logterm * &
                             (integral - this%gamma%cell%omega)
                    !
                CASE ('ions')
                    !
                    arg%of_r = 1.D0 - &
                               sumcbulk / this%cionmax * &
                               (1.D0 - this%gamma%of_r)
                    !
                    f%of_r = LOG(arg%of_r)
                    integral = f%integrate()
                    !
                    energy = -kT * this%cionmax * integral
                    !
                END SELECT
                !
            END IF
            !
            energy = energy + this%energy_second_order
        ELSE
            arg%of_r = 0.D0
            !
            DO ityp = 1, this%ntyp
                !
                arg%of_r = arg%of_r + this%ioncctype(ityp)%cfactor%of_r * &
                           this%ioncctype(ityp)%cbulk
                !
            END DO
            !
            IF (this%cionmax == 0.D0) THEN ! Full PB
                !
                f%of_r = this%gamma%of_r * arg%of_r - sumcbulk
                integral = f%integrate()
                !
                energy = -kT * integral
                !
            ELSE ! full modified PB
                !
                SELECT CASE (this%electrolyte_entropy)
                    !
                CASE ('full')
                    arg%of_r = arg%of_r / (this%cionmax - sumcbulk)
                    arg%of_r = arg%of_r + 1.D0
                    f%of_r = this%gamma%of_r * LOG(arg%of_r)
                    integral = f%integrate()
                    logterm = LOG(1.D0 - sumcbulk / this%cionmax)
                    !
                    energy = -kT * this%cionmax * &
                             (integral + logterm * this%gamma%cell%omega)
                    !
                CASE ('ions')
                    arg%of_r = arg%of_r / this%cionmax * this%gamma%of_r
                    arg%of_r = arg%of_r + 1.D0 - sumcbulk / this%cionmax
                    f%of_r = LOG(arg%of_r)
                    integral = f%integrate()
                    energy = -kT * this%cionmax * integral
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
        CLASS(environ_electrolyte), TARGET, INTENT(IN) :: this
        !
        TYPE(environ_density), TARGET, INTENT(INOUT) :: de_dboundary
        !
        REAL(DP), POINTER :: gam(:)
        !
        REAL(DP) :: kT, sumcbulk
        INTEGER :: ityp
        TYPE(environ_density) :: arg
        !
        !--------------------------------------------------------------------------------
        !
        gam => this%gamma%of_r
        !
        kT = K_BOLTZMANN_RY * this%temperature
        sumcbulk = SUM(this%ioncctype(:)%cbulk)
        !
        IF (this%linearized) THEN
            !
            IF (this%cionmax == 0.D0) THEN ! linearized PB
                !
                de_dboundary%of_r = de_dboundary%of_r - &
                     & this%dgamma%of_r * kT * sumcbulk
                !
            ELSE ! linearized modified PB
                !
                SELECT CASE (this%electrolyte_entropy)
                    !
                CASE ('full')
                    !
                    de_dboundary%of_r = &
                        de_dboundary%of_r + &
                        this%dgamma%of_r * kT * this%cionmax * &
                        LOG(1.D0 - sumcbulk / this%cionmax)
                    !
                CASE ('ions')
                    !
                    de_dboundary%of_r = &
                        de_dboundary%of_r - &
                        this%dgamma%of_r * kT * sumcbulk / &
                        (1.D0 - sumcbulk / this%cionmax * (1.D0 - gam))
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
            DO ityp = 1, this%ntyp
                !
                arg%of_r = arg%of_r + &
                           this%ioncctype(ityp)%cfactor%of_r * &
                           this%ioncctype(ityp)%cbulk
                !
            END DO
            !
            IF (this%cionmax == 0.D0) THEN ! full PB
                !
                de_dboundary%of_r = de_dboundary%of_r - &
                                    this%dgamma%of_r * kT * arg%of_r
                !
            ELSE ! full modified PB
                !
                SELECT CASE (this%electrolyte_entropy)
                    !
                CASE ('full')
                    !
                    arg%of_r = arg%of_r / (this%cionmax - sumcbulk)
                    arg%of_r = arg%of_r + 1.D0
                    !
                    de_dboundary%of_r = &
                        de_dboundary%of_r - &
                        this%dgamma%of_r * kT * &
                        this%cionmax * LOG(arg%of_r)
                    !
                CASE ('ions')
                    !
                    de_dboundary%of_r = &
                        de_dboundary%of_r - &
                        this%dgamma%of_r * kT * arg%of_r / &
                        (1.D0 - (sumcbulk - arg%of_r * gam) / this%cionmax)
                    !
                END SELECT
                !
            END IF
            !
            CALL arg%destroy()
            !
        END IF
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
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE print_environ_electrolyte(this, local_verbose, local_depth)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_electrolyte), INTENT(IN) :: this
        INTEGER, INTENT(IN), OPTIONAL :: local_verbose
        INTEGER, INTENT(IN), OPTIONAL :: local_depth
        !
        INTEGER :: verbosity, passed_verbosity, passed_depth, ityp
        !
        CHARACTER(LEN=80) :: sub_name = 'print_environ_electrolyte'
        !
        !--------------------------------------------------------------------------------
        !
        IF (verbose == 0) RETURN
        !
        IF (PRESENT(local_verbose)) THEN
            verbosity = verbose + local_verbose
        ELSE
            verbosity = verbose
        END IF
        !
        IF (verbosity == 0) RETURN
        !
        IF (PRESENT(local_depth)) THEN
            passed_verbosity = verbosity - verbose - local_depth
            passed_depth = local_depth
        ELSE
            passed_verbosity = verbosity - verbose - depth
            passed_depth = depth
        END IF
        !
        IF (verbosity >= 1) THEN
            !
            IF (ionode) THEN
                !
                IF (verbosity >= verbose) THEN ! header
                    WRITE (environ_unit, 1000)
                ELSE
                    !
                    CALL env_block_divider(verbosity)
                    !
                    WRITE (environ_unit, 1001)
                END IF
                !
                WRITE (environ_unit, 1002) &
                    this%ntyp, this%temperature, 1.D0 / SQRT(this%k2)
                !
                IF (this%cionmax > 0.D0) WRITE (environ_unit, 1003) this%cionmax
                !
            END IF
            !
            DO ityp = 1, this%ntyp
                !
                IF (ionode) &
                    WRITE (environ_unit, 1004) &
                    this%ioncctype(ityp)%index, this%ioncctype(ityp)%cbulk, &
                    this%ioncctype(ityp)%cbulk * AMU_SI / BOHR_RADIUS_SI**3, &
                    this%ioncctype(ityp)%z
                !
                IF (verbosity >= 5) THEN
                    !
                    CALL this%ioncctype(ityp)%c%printout(passed_verbosity, passed_depth)
                    !
                    CALL this%ioncctype(ityp)%cfactor%printout(passed_verbosity, &
                                                               passed_depth)
                    !
                END IF
                !
            END DO
            !
            IF (verbosity >= 3) THEN
                !
                CALL this%density%printout(passed_verbosity, passed_depth)
                !
                CALL this%gamma%printout(passed_verbosity, passed_depth)
                !
            END IF
            !
            IF (verbosity >= 5) CALL this%dgamma%printout(passed_verbosity, passed_depth)
            !
            IF (ionode) THEN
                WRITE (environ_unit, 1005) this%linearized
                WRITE (environ_unit, 1006) this%charge
            END IF
            !
            IF (verbosity < verbose) CALL env_block_divider(verbosity)
            !
        END IF
        !
        FLUSH (environ_unit)
        !
        !--------------------------------------------------------------------------------
        !
1000    FORMAT(/, 4('%'), ' ELECTROLYTE ', 64('%'))
1001    FORMAT(/, ' ELECTROLYTE', /, ' ===========')
        !
1002    FORMAT(/, ' number electrol. species   = ', I4, /, &
                ' solvent temperature        = ', F7.1, /, &
                ' Debye length / sqrt(eps)   = ', F14.7)
        !
1003    FORMAT(/, ' modified Poisson-Boltzmann:', /, &
                ' maximum concentration      = ', F14.7)
        !
1004    FORMAT(/, ' electrolyte species:', I4, /, &
                ' bulk concentration  (a.u.) = ', E15.4, /, &
                ' bulk concentration (mol/L) = ', F14.7, /, &
                ' ionic charge               = ', F7.2)
        !
1005    FORMAT(/, ' electrolyte flags:', /, &
                ' linearized                 = ', L2)
        !
1006    FORMAT(/, ' total electrolyte charge   = ', F14.7)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE print_environ_electrolyte
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE class_electrolyte
!----------------------------------------------------------------------------------------
