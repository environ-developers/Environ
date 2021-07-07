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
MODULE utils_electrolyte
    !------------------------------------------------------------------------------------
    !
    USE environ_param, ONLY: DP, e2, bohr_radius_si, amu_si, k_boltzmann_ry, fpi
    !
    USE types_core, ONLY: core_container
    USE types_cell, ONLY: environ_cell
    !
    USE types_physical, ONLY: environ_electrolyte, environ_electrons, environ_ions, &
                              environ_system
    !
    USE utils_boundary, ONLY: create_environ_boundary, init_environ_boundary_first, &
                              init_environ_boundary_second, destroy_environ_boundary
    !
    USE utils_density, ONLY: create_environ_density, init_environ_density, &
                             destroy_environ_density
    !
    USE tools_electrolyte, ONLY: electrolyte_of_boundary
    USE tools_functions, ONLY: density_of_functions
    !
    !------------------------------------------------------------------------------------
    !
    PRIVATE
    !
    PUBLIC :: create_environ_electrolyte, init_environ_electrolyte_first, &
              init_environ_electrolyte_second, update_environ_electrolyte, &
              destroy_environ_electrolyte
    !
    !------------------------------------------------------------------------------------
CONTAINS
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE create_environ_electrolyte(electrolyte)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_electrolyte), INTENT(INOUT) :: electrolyte
        !
        CHARACTER(LEN=80) :: label
        !
        CHARACTER(LEN=80) :: sub_name = 'create_environ_electrolyte'
        !
        !--------------------------------------------------------------------------------
        !
        label = 'electrolyte'
        !
        CALL create_environ_boundary(electrolyte%boundary, label)
        !
        CALL create_environ_density(electrolyte%density, label)
        !
        label = 'gamma'
        !
        CALL create_environ_density(electrolyte%gamma, label)
        !
        label = 'dgamma'
        !
        CALL create_environ_density(electrolyte%dgamma, label)
        !
        electrolyte%charge = 0.D0
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE create_environ_electrolyte
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_environ_electrolyte_first(ntyp, mode, stype, rhomax, rhomin, &
                                              tbeta, const, alpha, softness, distance, &
                                              spread, solvent_radius, radial_scale, &
                                              radial_spread, filling_threshold, &
                                              filling_spread, field_awareness, &
                                              charge_asymmetry, field_max, field_min, &
                                              electrons, ions, system, core, &
                                              temperature, cbulk, cionmax, radius, z, &
                                              electrolyte_entropy, &
                                              linearized, electrolyte)
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
        TYPE(core_container), TARGET, INTENT(IN) :: core
        TYPE(environ_electrolyte), INTENT(INOUT) :: electrolyte
        !
        INTEGER :: ityp
        REAL(DP) :: neutral, sumcbulk, amin, amax
        CHARACTER(LEN=80) :: ityps, label
        !
        CHARACTER(LEN=80) :: sub_name = 'init_environ_electrolyte_first'
        !
        !--------------------------------------------------------------------------------
        !
        CALL init_environ_boundary_first(.TRUE., .TRUE., .FALSE., mode, stype, &
                                         rhomax, rhomin, tbeta, const, alpha, &
                                         softness, distance, spread, solvent_radius, &
                                         radial_scale, radial_spread, &
                                         filling_threshold, filling_spread, &
                                         field_awareness, charge_asymmetry, field_max, &
                                         field_min, electrons, ions, system, core, &
                                         electrolyte%boundary)
        !
        !--------------------------------------------------------------------------------
        ! Setup all electrolyte parameters (with checks)
        !
        electrolyte%initialized = .FALSE.
        electrolyte%linearized = linearized
        electrolyte%ntyp = ntyp
        electrolyte%electrolyte_entropy = TRIM(electrolyte_entropy)
        electrolyte%temperature = temperature
        electrolyte%cionmax = 0.D0
        electrolyte%permittivity = const
        !
        IF (ALLOCATED(electrolyte%ioncctype)) &
            CALL env_errore(sub_name, &
                            'Trying to allocate an already allocated object', 1)
        !
        ALLOCATE (electrolyte%ioncctype(ntyp))
        !
        neutral = 0.D0
        !
        DO ityp = 1, ntyp
            electrolyte%ioncctype(ityp)%index = ityp
            !
            !----------------------------------------------------------------------------
            ! Convert bulk concentrations to atomic units
            !
            electrolyte%ioncctype(ityp)%cbulk = cbulk(ityp) * bohr_radius_si**3 / amu_si
            !
            electrolyte%ioncctype(ityp)%z = -z(ityp)
            neutral = neutral + cbulk(ityp) * z(ityp)
            !
            !----------------------------------------------------------------------------
            ! Create density for the local electrolyte concentration
            ! and related quantities
            !
            WRITE (ityps, '(I2.2)') ityp
            label = 'c_electrolyte_'//TRIM(ityps)
            !
            CALL create_environ_density(electrolyte%ioncctype(ityp)%c, label)
            !
            label = 'cfactor_electrolyte_'//TRIM(ityps)
            !
            CALL create_environ_density(electrolyte%ioncctype(ityp)%cfactor, label)
            !
        END DO
        !
        IF (neutral > 1.D-8) &
            CALL env_errore(sub_name, 'Bulk electrolyte is not neutral', 1)
        !
        !--------------------------------------------------------------------------------
        ! If cionmax is not provided in input but radius is, calculate cionmax
        !
        electrolyte%cionmax = cionmax * bohr_radius_si**3 / amu_si
        !
        IF (cionmax == 0.D0 .AND. radius > 0.D0) &
            electrolyte%cionmax = 0.64D0 * 3.D0 / fpi / radius**3
        !
        !--------------------------------------------------------------------------------
        ! Check suitability of cionmax value
        !
        sumcbulk = SUM(electrolyte%ioncctype(:)%cbulk)
        !
        IF (electrolyte%cionmax > 0.D0 .AND. electrolyte%cionmax <= sumcbulk) &
            CALL env_errore(sub_name, &
                            'cionmax should be larger than the sum of cbulks', 1)
        !
        electrolyte%energy_second_order = 0.D0
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_environ_electrolyte_first
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_environ_electrolyte_second(cell, electrolyte)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_cell), INTENT(IN) :: cell
        !
        TYPE(environ_electrolyte), INTENT(INOUT) :: electrolyte
        !
        INTEGER :: ityp
        REAL(DP) :: sum_cz2, arg, kT, e
        !
        !--------------------------------------------------------------------------------
        !
        sum_cz2 = 0.D0
        kT = k_boltzmann_ry * electrolyte%temperature
        !
        CALL init_environ_boundary_second(cell, electrolyte%boundary)
        !
        CALL init_environ_density(cell, electrolyte%gamma)
        !
        CALL init_environ_density(cell, electrolyte%dgamma)
        !
        CALL init_environ_density(cell, electrolyte%density)
        !
        DO ityp = 1, electrolyte%ntyp
            !
            CALL init_environ_density(cell, electrolyte%ioncctype(ityp)%c)
            !
            CALL init_environ_density(cell, electrolyte%ioncctype(ityp)%cfactor)
            !
            sum_cz2 = &
                sum_cz2 + &
                electrolyte%ioncctype(ityp)%cbulk * electrolyte%ioncctype(ityp)%z**2
            !
        END DO
        !
        !--------------------------------------------------------------------------------
        ! k^2 = eps / lambda_D^2
        !
        electrolyte%k2 = sum_cz2 / kT
        electrolyte%k2 = electrolyte%k2 * e2 * fpi
        !
        IF (electrolyte%linearized) &
            CALL init_environ_density(cell, electrolyte%de_dboundary_second_order)
        !
        electrolyte%initialized = .TRUE.
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_environ_electrolyte_second
    !------------------------------------------------------------------------------------
    !
    !>
    !------------------------------------------------------------------------------------
    SUBROUTINE update_environ_electrolyte(electrolyte)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_electrolyte), INTENT(INOUT) :: electrolyte
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
        IF (electrolyte%boundary%update_status > 0) electrolyte%update = .TRUE.
        !
        IF (electrolyte%update) THEN
            !
            !----------------------------------------------------------------------------
            ! Update the electrolyte in space if the boundary is ready
            !
            IF (electrolyte%boundary%update_status == 2) THEN
                !
                CALL electrolyte_of_boundary(electrolyte)
                !
                electrolyte%update = .FALSE.
                !
            END IF
            !
        END IF
        !
        CALL env_stop_clock(sub_name)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE update_environ_electrolyte
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE destroy_environ_electrolyte(lflag, electrolyte)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        LOGICAL, INTENT(IN) :: lflag
        !
        TYPE(environ_electrolyte), INTENT(INOUT) :: electrolyte
        !
        INTEGER :: ityp
        !
        CHARACTER(LEN=80) :: sub_name = 'destroy_environ_electrolyte'
        !
        !--------------------------------------------------------------------------------
        !
        CALL destroy_environ_boundary(lflag, electrolyte%boundary)
        !
        IF (electrolyte%initialized) THEN
            !
            DO ityp = 1, electrolyte%ntyp
                !
                CALL destroy_environ_density(electrolyte%ioncctype(ityp)%c)
                !
                CALL destroy_environ_density(electrolyte%ioncctype(ityp)%cfactor)
                !
            END DO
            !
            CALL destroy_environ_density(electrolyte%gamma)
            !
            CALL destroy_environ_density(electrolyte%dgamma)
            !
            CALL destroy_environ_density(electrolyte%density)
            !
            IF (electrolyte%linearized) &
                CALL destroy_environ_density(electrolyte%de_dboundary_second_order)
            !
            electrolyte%initialized = .FALSE.
            !
        END IF
        !
        IF (lflag) THEN
            !
            IF (.NOT. ALLOCATED(electrolyte%ioncctype)) &
                CALL env_errore(sub_name, 'Trying to destroy a non allocated object', 1)
            !
            DEALLOCATE (electrolyte%ioncctype)
        END IF
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE destroy_environ_electrolyte
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE utils_electrolyte
!----------------------------------------------------------------------------------------
