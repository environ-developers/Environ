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
! Authors: Oliviero Andreussi (Department of Physics, UNT)
!          Francesco Nattino  (THEOS and NCCR-MARVEL, EPFL)
!          Nicola Marzari     (THEOS and NCCR-MARVEL, EPFL)
!
!----------------------------------------------------------------------------------------
!> Module containing the main routines to handle environ_electrolyte
!! derived data types.
!!
!! Environ_electrolyte contains all the specifications and the details of
!! the electrolyte medium and of the ionic distribution in it.
!!
!----------------------------------------------------------------------------------------
MODULE tools_electrolyte
    !------------------------------------------------------------------------------------
    !
    USE modules_constants, ONLY: DP, e2, k_boltzmann_ry, fpi
    !
    USE physical_types, ONLY: environ_electrolyte, environ_electrons
    USE representation_types, ONLY: environ_density
    !
    USE utils_density, ONLY: init_environ_density, destroy_environ_density
    !
    USE tools_math, ONLY: scalar_product_environ_density, integrate_environ_density
    !
    !------------------------------------------------------------------------------------
    !
    PRIVATE
    !
    PUBLIC :: electrolyte_of_boundary, electrolyte_of_potential, calc_eelectrolyte, &
              calc_deelectrolyte_dboundary
    !
    !------------------------------------------------------------------------------------
CONTAINS
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE electrolyte_of_boundary(electrolyte)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_electrolyte), TARGET, INTENT(INOUT) :: electrolyte
        !
        CHARACTER(LEN=80) :: sub_name = 'electrolyte_of_boundary'
        !
        !--------------------------------------------------------------------------------
        ! Compute exclusion function gamma(r) and dgamma/ds(r)
        !
        electrolyte%gamma%of_r = 1.D0 - electrolyte%boundary%scaled%of_r
        electrolyte%dgamma%of_r = -1.D0
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE electrolyte_of_boundary
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE electrolyte_of_potential(potential, electrolyte)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_density), TARGET, INTENT(IN) :: potential
        !
        TYPE(environ_electrolyte), TARGET, INTENT(INOUT) :: electrolyte
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
        gam => electrolyte%gamma%of_r
        pot => potential%of_r
        rho => electrolyte%density%of_r
        !
        rho = 0.D0
        kT = k_boltzmann_ry * electrolyte%temperature
        !
        CALL init_environ_density(potential%cell, denominator)
        !
        denominator%of_r = 1.D0
        !
        DO ityp = 1, electrolyte%ntyp
            cfactor => electrolyte%ioncctype(ityp)%cfactor%of_r
            cbulk => electrolyte%ioncctype(ityp)%cbulk
            z => electrolyte%ioncctype(ityp)%z
            cfactor = 1.D0
            !
            IF (electrolyte%linearized) THEN
                !
                !------------------------------------------------------------------------
                ! Linearized PB and modified PB
                !
                cfactor = 1.D0 - z * pot / kT
                !
                IF (electrolyte%cionmax > 0.D0) THEN
                    factor = cbulk / electrolyte%cionmax
                    !
                    IF (electrolyte%electrolyte_entropy == 'ions') &
                        denominator%of_r = denominator%of_r - factor * (1.D0 - gam)
                    !
                END IF
                !
            ELSE
                !
                !------------------------------------------------------------------------
                ! Full PB
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
                IF (electrolyte%cionmax /= 0.D0) THEN
                    !
                    !--------------------------------------------------------------------
                    ! Full modified PB
                    !
                    factor = cbulk / electrolyte%cionmax
                    !
                    SELECT CASE (electrolyte%electrolyte_entropy)
                    CASE ('full')
                        !
                        denominator%of_r = denominator%of_r - &
                                           factor * (1.D0 - cfactor)
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
        DO ityp = 1, electrolyte%ntyp
            c => electrolyte%ioncctype(ityp)%c%of_r
            cfactor => electrolyte%ioncctype(ityp)%cfactor%of_r
            cbulk => electrolyte%ioncctype(ityp)%cbulk
            z => electrolyte%ioncctype(ityp)%z
            c = gam * cbulk * cfactor / denominator%of_r
            rho = rho + c * z
            !
            NULLIFY (c)
            NULLIFY (cfactor)
            NULLIFY (cbulk)
            NULLIFY (z)
        END DO
        !
        electrolyte%charge = integrate_environ_density(electrolyte%density)
        !
        IF (electrolyte%linearized) THEN
            !
            !----------------------------------------------------------------------------
            ! Compute energy and de_dboundary terms that depend on the potential.
            ! These are the second order terms, first order terms are equal to zero.
            ! NOTE: the energy is equal to minus the electrostatic interaction of the
            ! electrolyte, and is identical for the various implementations of the entropy.
            !
            electrolyte%energy_second_order = &
                0.5D0 * scalar_product_environ_density(electrolyte%density, potential)
            !
            IF (electrolyte%electrolyte_entropy == 'ions' .AND. &
                electrolyte%cionmax > 0.D0) THEN
                !
                sumcbulk = SUM(electrolyte%ioncctype(:)%cbulk)
                !
                electrolyte%de_dboundary_second_order%of_r = &
                    -0.5D0 * electrolyte%k2 / e2 / fpi * &
                    (1.D0 - sumcbulk / electrolyte%cionmax) * &
                    (pot / denominator%of_r)**2 * electrolyte%dgamma%of_r
                !
            ELSE
                !
                electrolyte%de_dboundary_second_order%of_r = &
                    -0.5D0 * electrolyte%k2 / e2 / fpi * &
                    pot * pot * electrolyte%dgamma%of_r
                !
            END IF
            !
        END IF
        !
        CALL destroy_environ_density(denominator)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE electrolyte_of_potential
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_eelectrolyte(electrolyte, energy)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_electrolyte), INTENT(IN) :: electrolyte
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
        kT = k_boltzmann_ry * electrolyte%temperature
        sumcbulk = SUM(electrolyte%ioncctype(:)%cbulk)
        !
        CALL init_environ_density(electrolyte%gamma%cell, arg)
        !
        CALL init_environ_density(electrolyte%gamma%cell, f)
        !
        IF (electrolyte%linearized) THEN
            !
            IF (electrolyte%cionmax == 0.D0) THEN
                !
                !------------------------------------------------------------------------
                ! Linearized PB
                !
                integral = integrate_environ_density(electrolyte%gamma)
                energy = kT * sumcbulk * (electrolyte%gamma%cell%omega - integral)
                !
            ELSE
                !
                !------------------------------------------------------------------------
                ! Linearized modified PB
                !
                SELECT CASE (electrolyte%electrolyte_entropy)
                CASE ('full')
                    integral = integrate_environ_density(electrolyte%gamma)
                    logterm = LOG(1.D0 - sumcbulk / electrolyte%cionmax)
                    !
                    energy = kT * electrolyte%cionmax * logterm * &
                             (integral - electrolyte%gamma%cell%omega)
                    !
                CASE ('ions')
                    !
                    arg%of_r = 1.D0 - &
                               sumcbulk / electrolyte%cionmax * &
                               (1.D0 - electrolyte%gamma%of_r)
                    !
                    f%of_r = LOG(arg%of_r)
                    integral = integrate_environ_density(f)
                    !
                    energy = -kT * electrolyte%cionmax * integral
                END SELECT
                !
            END IF
            !
            energy = energy + electrolyte%energy_second_order
        ELSE
            arg%of_r = 0.D0
            !
            DO ityp = 1, electrolyte%ntyp
                !
                arg%of_r = arg%of_r + electrolyte%ioncctype(ityp)%cfactor%of_r * &
                           electrolyte%ioncctype(ityp)%cbulk
                !
            END DO
            !
            IF (electrolyte%cionmax == 0.D0) THEN
                !
                !------------------------------------------------------------------------
                ! Full PB
                !
                f%of_r = electrolyte%gamma%of_r * arg%of_r - sumcbulk
                integral = integrate_environ_density(f)
                !
                energy = -kT * integral
            ELSE
                !
                !------------------------------------------------------------------------
                ! Full modified PB
                !
                SELECT CASE (electrolyte%electrolyte_entropy)
                CASE ('full')
                    arg%of_r = arg%of_r / (electrolyte%cionmax - sumcbulk)
                    arg%of_r = arg%of_r + 1.D0
                    f%of_r = electrolyte%gamma%of_r * LOG(arg%of_r)
                    integral = integrate_environ_density(f)
                    logterm = LOG(1.D0 - sumcbulk / electrolyte%cionmax)
                    !
                    energy = -kT * electrolyte%cionmax * &
                             (integral + logterm * electrolyte%gamma%cell%omega)
                    !
                CASE ('ions')
                    arg%of_r = arg%of_r / electrolyte%cionmax * electrolyte%gamma%of_r
                    arg%of_r = arg%of_r + 1.D0 - sumcbulk / electrolyte%cionmax
                    f%of_r = LOG(arg%of_r)
                    integral = integrate_environ_density(f)
                    energy = -kT * electrolyte%cionmax * integral
                END SELECT
                !
            END IF
            !
        END IF
        !
        CALL destroy_environ_density(arg)
        !
        CALL destroy_environ_density(f)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_eelectrolyte
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_deelectrolyte_dboundary(electrolyte, de_dboundary)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_electrolyte), TARGET, INTENT(IN) :: electrolyte
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
        gam => electrolyte%gamma%of_r
        !
        kT = k_boltzmann_ry * electrolyte%temperature
        sumcbulk = SUM(electrolyte%ioncctype(:)%cbulk)
        !
        IF (electrolyte%linearized) THEN
            !
            IF (electrolyte%cionmax == 0.D0) THEN
                !
                !------------------------------------------------------------------------
                ! Linearized PB
                !
                de_dboundary%of_r = de_dboundary%of_r - &
                     & electrolyte%dgamma%of_r * kT * sumcbulk
                !
            ELSE
                !
                !------------------------------------------------------------------------
                ! Linearized modified PB
                !
                SELECT CASE (electrolyte%electrolyte_entropy)
                CASE ('full')
                    !
                    de_dboundary%of_r = &
                        de_dboundary%of_r + &
                        electrolyte%dgamma%of_r * kT * electrolyte%cionmax * &
                        LOG(1.D0 - sumcbulk / electrolyte%cionmax)
                    !
                CASE ('ions')
                    !
                    de_dboundary%of_r = &
                        de_dboundary%of_r - &
                        electrolyte%dgamma%of_r * kT * sumcbulk / &
                        (1.D0 - sumcbulk / electrolyte%cionmax * (1.D0 - gam))
                    !
                END SELECT
                !
            END IF
            !
            de_dboundary%of_r = de_dboundary%of_r + &
                                electrolyte%de_dboundary_second_order%of_r
            !
        ELSE
            !
            CALL init_environ_density(de_dboundary%cell, arg)
            !
            arg%of_r = 0.D0
            !
            DO ityp = 1, electrolyte%ntyp
                !
                arg%of_r = arg%of_r + &
                           electrolyte%ioncctype(ityp)%cfactor%of_r * &
                           electrolyte%ioncctype(ityp)%cbulk
                !
            END DO
            !
            IF (electrolyte%cionmax == 0.D0) THEN
                !
                !------------------------------------------------------------------------
                ! Full PB
                !
                de_dboundary%of_r = de_dboundary%of_r - &
                                    electrolyte%dgamma%of_r * kT * arg%of_r
                !
            ELSE
                !
                !------------------------------------------------------------------------
                ! Full modified PB
                !
                SELECT CASE (electrolyte%electrolyte_entropy)
                CASE ('full')
                    !
                    arg%of_r = arg%of_r / (electrolyte%cionmax - sumcbulk)
                    arg%of_r = arg%of_r + 1.D0
                    !
                    de_dboundary%of_r = &
                        de_dboundary%of_r - &
                        electrolyte%dgamma%of_r * kT * &
                        electrolyte%cionmax * LOG(arg%of_r)
                    !
                CASE ('ions')
                    !
                    de_dboundary%of_r = &
                        de_dboundary%of_r - &
                        electrolyte%dgamma%of_r * kT * arg%of_r / &
                        (1.D0 - (sumcbulk - arg%of_r * gam) / electrolyte%cionmax)
                    !
                END SELECT
                !
            END IF
            !
            CALL destroy_environ_density(arg)
            !
        END IF
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_deelectrolyte_dboundary
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE tools_electrolyte
!----------------------------------------------------------------------------------------
