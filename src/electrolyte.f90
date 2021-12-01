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
        LOGICAL :: lupdate = .FALSE.
        !
        TYPE(environ_electrolyte_base) :: base
        !
        !--------------------------------------------------------------------------------
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
        REAL(DP) :: energy_second_order = 0.D0
        !
        REAL(DP) :: charge = 0.D0
        !
        !--------------------------------------------------------------------------------
    CONTAINS
        !--------------------------------------------------------------------------------
        !
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
    SUBROUTINE init_environ_electrolyte(this, ntyp, mode, stype, rhomax, rhomin, &
                                        tbeta, const, alpha, softness, distance, &
                                        spread, solvent_radius, radial_scale, &
                                        radial_spread, filling_threshold, &
                                        filling_spread, field_awareness, &
                                        charge_asymmetry, field_max, field_min, &
                                        electrons, ions, system, temperature, cbulk, &
                                        cionmax, radius, z, electrolyte_entropy, &
                                        linearized, cores, deriv_method, cell)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        LOGICAL, INTENT(IN) :: linearized
        INTEGER, INTENT(IN) :: ntyp, stype
        CHARACTER(LEN=*), INTENT(IN) :: mode, electrolyte_entropy, deriv_method
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
        TYPE(core_container), TARGET, INTENT(IN) :: cores
        TYPE(environ_cell), INTENT(IN) :: cell
        !
        CLASS(environ_electrolyte), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: sub_name = 'init_environ_electrolyte'
        !
        !--------------------------------------------------------------------------------
        !
        CALL this%base%init(ntyp, const, distance, spread, temperature, cbulk, cionmax, &
                            radius, z, electrolyte_entropy, linearized, cell)
        !
        CALL this%boundary%init(.TRUE., .TRUE., .FALSE., mode, stype, rhomax, &
                                rhomin, tbeta, const, alpha, softness, distance, &
                                spread, solvent_radius, radial_scale, radial_spread, &
                                filling_threshold, filling_spread, field_awareness, &
                                charge_asymmetry, field_max, field_min, electrons, &
                                ions, system, cores, deriv_method, cell, 'electrolyte')
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
        CHARACTER(LEN=80) :: sub_name = 'destroy_environ_electrolyte'
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
        REAL(DP), POINTER :: z, cbulk, cionmax
        REAL(DP), DIMENSION(:), POINTER :: pot, rho, c, cfactor, gam
        TYPE(environ_electrolyte_base), POINTER :: base
        !
        REAL(DP) :: kT, factor, sumcbulk, arg
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
        base => this%base
        cionmax => base%cionmax
        !
        rho = 0.D0
        kT = K_BOLTZMANN_RY * base%temperature
        !
        CALL denominator%init(potential%cell)
        !
        denominator%of_r = 1.D0
        !
        DO ityp = 1, base%ntyp
            cfactor => base%ioncctype(ityp)%cfactor%of_r
            cbulk => base%ioncctype(ityp)%cbulk
            z => base%ioncctype(ityp)%z
            !
            cfactor = 1.D0
            !
            IF (base%linearized) THEN ! Linearized PB and modified PB
                cfactor = 1.D0 - z * pot / kT
                !
                IF (cionmax > 0.D0) THEN
                    factor = cbulk / cionmax
                    !
                    IF (base%electrolyte_entropy == 'ions') &
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
                IF (cionmax /= 0.D0) THEN ! full modified PB
                    factor = cbulk / cionmax
                    !
                    SELECT CASE (base%electrolyte_entropy)
                        !
                    CASE ('full')
                        denominator%of_r = denominator%of_r - factor * (1.D0 - cfactor)
                        !
                    CASE ('ions')
                        !
                        denominator%of_r = denominator%of_r - &
                                           factor * (1.D0 - gam * cfactor)
                        !
                    CASE DEFAULT
                        CALL io%error(sub_name, 'Unexpected electrolyte entropy', 1)
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
        END DO
        !
        DO ityp = 1, base%ntyp
            c => base%ioncctype(ityp)%c%of_r
            cfactor => base%ioncctype(ityp)%cfactor%of_r
            cbulk => base%ioncctype(ityp)%cbulk
            z => base%ioncctype(ityp)%z
            !
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
        IF (base%linearized) THEN
            !
            !----------------------------------------------------------------------------
            ! Compute energy and de_dboundary terms that depend on the potential.
            ! These are the second order terms, first order terms are equal to zero.
            ! NOTE: the energy is equal to minus the electrostatic interaction of the
            ! electrolyte, and is identical for the various implementations of the entropy.
            !
            this%energy_second_order = 0.5D0 * this%density%scalar_product(potential)
            !
            IF (base%electrolyte_entropy == 'ions' .AND. cionmax > 0.D0) THEN
                !
                sumcbulk = SUM(base%ioncctype%cbulk)
                !
                this%de_dboundary_second_order%of_r = &
                    -0.5D0 * base%k2 / e2 / fpi * (1.D0 - sumcbulk / cionmax) * &
                    (pot / denominator%of_r)**2 * this%dgamma%of_r
                !
            ELSE
                !
                this%de_dboundary_second_order%of_r = -0.5D0 * base%k2 / e2 / fpi * &
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
        CLASS(environ_electrolyte), TARGET, INTENT(IN) :: this
        !
        REAL(DP), INTENT(OUT) :: energy
        !
        INTEGER :: ityp
        REAL(DP) :: kT, sumcbulk, logterm, integral

        REAL(DP), POINTER :: cionmax
        TYPE(environ_electrolyte_base), POINTER :: base
        !
        TYPE(environ_density) :: arg, f
        !
        CHARACTER(LEN=80) :: sub_name = 'calc_eelectrolyte'
        !
        !--------------------------------------------------------------------------------
        !
        base => this%base
        cionmax => base%cionmax
        !
        energy = 0.D0
        !
        kT = K_BOLTZMANN_RY * base%temperature
        sumcbulk = SUM(base%ioncctype%cbulk)
        !
        CALL arg%init(this%gamma%cell)
        !
        CALL f%init(this%gamma%cell)
        !
        IF (base%linearized) THEN
            !
            IF (cionmax == 0.D0) THEN ! Linearized PB
                integral = this%gamma%integrate()
                energy = kT * sumcbulk * (this%gamma%cell%omega - integral)
            ELSE ! Linearized modified PB
                !
                SELECT CASE (base%electrolyte_entropy)
                    !
                CASE ('full')
                    integral = this%gamma%integrate()
                    logterm = LOG(1.D0 - sumcbulk / cionmax)
                    energy = kT * cionmax * logterm * (integral - this%gamma%cell%omega)
                    !
                CASE ('ions')
                    arg%of_r = 1.D0 - sumcbulk / cionmax * (1.D0 - this%gamma%of_r)
                    f%of_r = LOG(arg%of_r)
                    integral = f%integrate()
                    energy = -kT * cionmax * integral
                    !
                CASE DEFAULT
                    CALL io%error(sub_name, 'Unexpected electrolyte entropy', 1)
                    !
                END SELECT
                !
            END IF
            !
            energy = energy + this%energy_second_order
        ELSE
            arg%of_r = 0.D0
            !
            DO ityp = 1, base%ntyp
                !
                arg%of_r = arg%of_r + base%ioncctype(ityp)%cfactor%of_r * &
                           base%ioncctype(ityp)%cbulk
                !
            END DO
            !
            IF (cionmax == 0.D0) THEN ! Full PB
                f%of_r = this%gamma%of_r * arg%of_r - sumcbulk
                integral = f%integrate()
                energy = -kT * integral
            ELSE ! full modified PB
                !
                SELECT CASE (base%electrolyte_entropy)
                    !
                CASE ('full')
                    arg%of_r = arg%of_r / (cionmax - sumcbulk)
                    arg%of_r = arg%of_r + 1.D0
                    f%of_r = this%gamma%of_r * LOG(arg%of_r)
                    integral = f%integrate()
                    logterm = LOG(1.D0 - sumcbulk / cionmax)
                    energy = -kT * cionmax * (integral + logterm * this%gamma%cell%omega)
                    !
                CASE ('ions')
                    arg%of_r = arg%of_r / cionmax * this%gamma%of_r
                    arg%of_r = arg%of_r + 1.D0 - sumcbulk / cionmax
                    f%of_r = LOG(arg%of_r)
                    integral = f%integrate()
                    energy = -kT * cionmax * integral
                    !
                CASE DEFAULT
                    CALL io%error(sub_name, 'Unexpected electrolyte entropy', 1)
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
        REAL(DP), POINTER :: gam(:), cionmax
        TYPE(environ_electrolyte_base), POINTER :: base
        !
        INTEGER :: ityp
        REAL(DP) :: kT, sumcbulk
        !
        TYPE(environ_density) :: arg
        !
        CHARACTER(LEN=80) :: sub_name = 'calc_deelectrolyte_dboundary'
        !
        !--------------------------------------------------------------------------------
        !
        gam => this%gamma%of_r
        base => this%base
        cionmax => base%cionmax
        !
        kT = K_BOLTZMANN_RY * base%temperature
        sumcbulk = SUM(base%ioncctype%cbulk)
        !
        IF (base%linearized) THEN
            !
            IF (cionmax == 0.D0) THEN ! linearized PB
                de_dboundary%of_r = de_dboundary%of_r - this%dgamma%of_r * kT * sumcbulk
            ELSE ! linearized modified PB
                !
                SELECT CASE (base%electrolyte_entropy)
                    !
                CASE ('full')
                    !
                    de_dboundary%of_r = de_dboundary%of_r + &
                                        this%dgamma%of_r * kT * cionmax * &
                                        LOG(1.D0 - sumcbulk / cionmax)
                    !
                CASE ('ions')
                    !
                    de_dboundary%of_r = de_dboundary%of_r - &
                                        this%dgamma%of_r * kT * sumcbulk / &
                                        (1.D0 - sumcbulk / cionmax * (1.D0 - gam))
                    !
                CASE DEFAULT
                    CALL io%error(sub_name, 'Unexpected electrolyte entropy', 1)
                    !
                END SELECT
                !
            END IF
            !
            de_dboundary%of_r = de_dboundary%of_r + this%de_dboundary_second_order%of_r
        ELSE
            !
            CALL arg%init(de_dboundary%cell)
            !
            arg%of_r = 0.D0
            !
            DO ityp = 1, base%ntyp
                !
                arg%of_r = arg%of_r + &
                           base%ioncctype(ityp)%cfactor%of_r * &
                           base%ioncctype(ityp)%cbulk
                !
            END DO
            !
            IF (cionmax == 0.D0) THEN ! full PB
                de_dboundary%of_r = de_dboundary%of_r - this%dgamma%of_r * kT * arg%of_r
            ELSE ! full modified PB
                !
                SELECT CASE (base%electrolyte_entropy)
                    !
                CASE ('full')
                    arg%of_r = arg%of_r / (cionmax - sumcbulk)
                    arg%of_r = arg%of_r + 1.D0
                    !
                    de_dboundary%of_r = de_dboundary%of_r - &
                                        this%dgamma%of_r * kT * cionmax * LOG(arg%of_r)
                    !
                CASE ('ions')
                    !
                    de_dboundary%of_r = de_dboundary%of_r - &
                                        this%dgamma%of_r * kT * arg%of_r / &
                                        (1.D0 - (sumcbulk - arg%of_r * gam) / cionmax)
                    !
                CASE DEFAULT
                    CALL io%error(sub_name, 'Unexpected electrolyte entropy', 1)
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
        INTEGER, INTENT(IN), OPTIONAL :: verbose, debug_verbose, unit
        !
        INTEGER :: base_verbose, local_verbose, passed_verbose, local_unit, ityp
        !
        TYPE(environ_electrolyte_base), POINTER :: base
        !
        CHARACTER(LEN=80) :: sub_name = 'print_environ_electrolyte'
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
            DO ityp = 1, base%ntyp
                !
                IF (io%lnode) &
                    WRITE (local_unit, 1007) &
                    base%ioncctype(ityp)%index, base%ioncctype(ityp)%cbulk, &
                    base%ioncctype(ityp)%cbulk * AMU_SI / BOHR_RADIUS_SI**3, &
                    base%ioncctype(ityp)%z
                !
            END DO
            !
            IF (local_verbose >= 5) THEN
                !
                DO ityp = 1, base%ntyp
                    !
                    CALL base%ioncctype(ityp)%c%printout(passed_verbose, debug_verbose, &
                                                         local_unit)
                    !
                    CALL base%ioncctype(ityp)%cfactor%printout(passed_verbose, &
                                                               debug_verbose, &
                                                               local_unit)
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
1000    FORMAT(/, 4('%'), ' ELECTROLYTE ', 64('%'))
        !
1001    FORMAT(/, ' number electrolyte species = ', I14, /, &
                ' solvent temperature        = ', F14.1, /, &
                ' Debye length / sqrt(eps)   = ', F14.7)
        !
1002    FORMAT(/, ' modified Poisson-Boltzmann:', /, &
                ' maximum concentration      = ', F14.7)
        !
1003    FORMAT(/, ' electrolyte flags:', /, &
                ' linearized                 = ', L14)
        !
1004    FORMAT(/, ' total electrolyte charge   = ', F14.7)
        !
1005    FORMAT(/, ' species', /, 1X, 7('='))
        !
1006    FORMAT(/, '   i |  c_bulk (a.u.) | c_bulk (mol/L) | ionic charge', /, 1X, 52('-'))
        !
1007    FORMAT(1X, I3, ' | ', E14.4, ' | ', F14.7, ' | ', F12.2)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE print_environ_electrolyte
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE class_electrolyte
!----------------------------------------------------------------------------------------
