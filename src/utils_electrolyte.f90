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
!----------------------------------------------------------------------------------------
!  TYPE environ_electrolyte
!----------------------------------------------------------------------------------------
!         !
!         LOGICAL :: update = .FALSE.
!         LOGICAL :: initialized = .FALSE.
!         CHARACTER(LEN=80) :: electrolyte_entropy
!         CHARACTER(LEN=80) :: ion_adsorption
!         LOGICAL :: linearized = .FALSE.
!         INTEGER :: ntyp
!         TYPE(environ_ioncctype), ALLOCATABLE :: ioncctype(:)
!         !
!         REAL(DP) :: temperature
!         REAL(DP) :: k2
!         REAL(DP) :: cionmax
!         REAL(DP) :: permittivity
!         REAL(DP) :: adsorption_energy
!         !
!         TYPE(environ_boundary) :: boundary
!         TYPE(environ_density) :: density
!         !
!         !--------------------------------------------------------------------------------
!         ! The electrolyte switch function and related quantities
!         !
!         TYPE(environ_density) :: gamma
!         TYPE(environ_density) :: dgamma
!         !
!         TYPE(environ_functions) :: function_
!         !
!         TYPE(environ_density) :: de_dboundary_second_order
!         REAL(DP) :: energy_second_order
!         REAL(DP) :: charge = 0.0_DP
!         !
!----------------------------------------------------------------------------------------
!  END TYPE environ_electrolyte
!----------------------------------------------------------------------------------------
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
MODULE utils_electrolyte
    !------------------------------------------------------------------------------------
    !
    USE modules_constants, ONLY: DP, e2, bohr_radius_si, amu_si, k_boltzmann_ry, fpi
    !
    USE physical_types, ONLY: environ_electrolyte, environ_electrons, environ_ions, &
                              environ_system
    !
    USE core_types, ONLY: boundary_core
    USE cell_types, ONLY: environ_cell
    !
    USE utils_boundary, ONLY: create_environ_boundary, init_environ_boundary_first, &
                              init_environ_boundary_second, destroy_environ_boundary
    !
    USE utils_density, ONLY: create_environ_density, init_environ_density, &
                             destroy_environ_density
    !
    USE tools_functions, ONLY: density_of_functions
    !
    !------------------------------------------------------------------------------------
    !
    PRIVATE :: electrolyte_of_boundary
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
                                              electrolyte_entropy, ion_adsorption, &
                                              adsorption_energy, linearized, &
                                              electrolyte)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        LOGICAL, INTENT(IN) :: linearized
        INTEGER, INTENT(IN) :: ntyp, stype
        CHARACTER(LEN=80), INTENT(IN) :: mode, electrolyte_entropy, ion_adsorption
        !
        REAL(DP), INTENT(IN) :: rhomax, rhomin, tbeta, const, distance, spread, &
                                alpha, softness, temperature, solvent_radius, &
                                radial_scale, radial_spread, filling_threshold, &
                                filling_spread, field_awareness, charge_asymmetry, &
                                field_max, field_min, cionmax, radius, adsorption_energy
        !
        REAL(DP), DIMENSION(ntyp), INTENT(IN) :: cbulk, z
        !
        TYPE(environ_electrons), INTENT(IN) :: electrons
        TYPE(environ_ions), INTENT(IN) :: ions
        TYPE(environ_system), TARGET, INTENT(IN) :: system
        TYPE(boundary_core), TARGET, INTENT(IN) :: core
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
        ! electrolyte%ion_adsorption = TRIM(ion_adsorption) #TODO
        electrolyte%temperature = temperature
        electrolyte%cionmax = 0.D0
        electrolyte%permittivity = const
        !
        ! IF (electrolyte%ion_adsorption .NE. 'none') THEN #TODO
        !     !
        !     ! Setup exponential function
        !     !
        !     electrolyte%adsorption_energy = adsorption_energy
        !     electrolyte%function_%type_ = 3
        !     electrolyte%function_%pos => system%pos
        !     electrolyte%function_%dim = system%dim
        !     electrolyte%function_%axis = system%axis
        !     electrolyte%function_%width = distance
        !     electrolyte%function_%spread = spread
        !     !
        ! END IF
        !
        IF (ALLOCATED(electrolyte%ioncctype)) &
            CALL errore(sub_name, 'Trying to allocate an already allocated object', 1)
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
            ! IF (electrolyte%ion_adsorption .NE. 'none') THEN #TODO
            !     label = 'potential_'//TRIM(ityps)
            !     CALL create_environ_density(electrolyte%ioncctype(ityp)%potential, label)
            ! END IF
            !
        END DO
        !
        IF (neutral .GT. 1.D-8) &
            CALL errore(sub_name, 'Bulk electrolyte is not neutral', 1)
        !
        !--------------------------------------------------------------------------------
        ! If cionmax is not provided in input but radius is, calculate cionmax
        !
        electrolyte%cionmax = cionmax * bohr_radius_si**3 / amu_si
        !
        IF (cionmax .EQ. 0.D0 .AND. radius .GT. 0.D0) &
            electrolyte%cionmax = 0.64D0 * 3.D0 / fpi / radius**3
        !
        !--------------------------------------------------------------------------------
        ! Check suitability of cionmax value
        !
        sumcbulk = SUM(electrolyte%ioncctype(:)%cbulk)
        !
        IF (electrolyte%cionmax .GT. 0.D0 .AND. electrolyte%cionmax .LE. sumcbulk) &
            CALL errore(sub_name, 'cionmax should be larger than the sum of cbulks', 1)
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
            ! IF (electrolyte%ion_adsorption .NE. 'none') & #TODO
            !     CALL init_environ_density(cell, electrolyte%ioncctype(ityp)%potential)
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
        !--------------------------------------------------------------------------------
        !
        CALL start_clock('electrolyte')
        !
        !--------------------------------------------------------------------------------
        ! Check if the boundary is under update (status = 1) or
        ! has been fully updated (status = 2)
        !
        IF (electrolyte%boundary%update_status .GT. 0) electrolyte%update = .TRUE.
        !
        IF (electrolyte%update) THEN
            !
            !----------------------------------------------------------------------------
            ! Update the electrolyte in space if the boundary is ready
            !
            IF (electrolyte%boundary%update_status .EQ. 2) THEN
                !
                CALL electrolyte_of_boundary(electrolyte)
                !
                electrolyte%update = .FALSE.
                !
                ! SELECT CASE (electrolyte%ion_adsorption) #TODO
                ! CASE ('anion')
                !     !
                !     DO ityp = 1, electrolyte%ntyp
                !         !
                !         IF (electrolyte%ioncctype(ityp)%z .GT. 0.D0) THEN
                !             !
                !             CALL density_of_functions( &
                !                 electrolyte%function_, &
                !                 electrolyte%ioncctype(ityp)%potential, .TRUE.)
                !             !
                !             electrolyte%ioncctype(ityp)%potential%of_r(:) = &
                !                 electrolyte%adsorption_energy * &
                !                 (electrolyte%ioncctype(ityp)%potential%of_r(:)**2 - &
                !                  electrolyte%ioncctype(ityp)%potential%of_r(:) * 2.D0)
                !             !
                !         ELSE
                !             electrolyte%ioncctype(ityp)%potential%of_r(:) = 0.D0
                !         END IF
                !         !
                !     END DO
                !     !
                ! CASE ('cation')
                !     !
                !     DO ityp = 1, electrolyte%ntyp
                !         !
                !         IF (electrolyte%ioncctype(ityp)%z .LT. 0.D0) THEN
                !             !
                !             CALL density_of_functions( &
                !                 electrolyte%function_, &
                !                 electrolyte%ioncctype(ityp)%potential, .TRUE.)
                !             !
                !             electrolyte%ioncctype(ityp)%potential%of_r(:) = &
                !                 electrolyte%adsorption_energy * &
                !                 (electrolyte%ioncctype(ityp)%potential%of_r(:)**2 - &
                !                  electrolyte%ioncctype(ityp)%potential%of_r(:) * 2.D0)
                !             !
                !         ELSE
                !             electrolyte%ioncctype(ityp)%potential%of_r(:) = 0.D0
                !         END IF
                !         !
                !     END DO
                !     !
                ! CASE ('repulsion')
                !     !
                !     DO ityp = 1, electrolyte%ntyp
                !         !
                !         CALL density_of_functions( &
                !             electrolyte%function_, &
                !             electrolyte%ioncctype(ityp)%potential, .TRUE.)
                !         !
                !     END DO
                !     !
                ! END SELECT
                !
            END IF
            !
        END IF
        !
        CALL stop_clock('electrolyte')
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
                ! IF (electrolyte%ion_adsorption .NE. 'none') & ! #TODO
                !     CALL destroy_environ_density(electrolyte%ioncctype(ityp)%potential)
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
                CALL errore(sub_name, 'Trying to destroy a non allocated object', 1)
            !
            DEALLOCATE (electrolyte%ioncctype)
        END IF
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE destroy_environ_electrolyte
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
    !
    !------------------------------------------------------------------------------------
END MODULE utils_electrolyte
!----------------------------------------------------------------------------------------
