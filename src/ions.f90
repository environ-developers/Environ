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
!          Edan Bainglass     (Department of Physics, UNT)
!
!----------------------------------------------------------------------------------------
!>
!! Module containing the main routines to handle environ_ions
!! derived data types.
!!
!! Environ_ions contains all the details of the ions of the QM system,
!! including the atomic type, mass and charge, the solvation radii,
!! the sizes of the ionic cores and the spread to be used when treating
!! the ions as gaussians. The module also contains all the routines
!! to handle environ_ions and to generate smooth ionic density from the
!! ions specification.
!!
!----------------------------------------------------------------------------------------
MODULE class_ions
    !------------------------------------------------------------------------------------
    !
    USE class_io, ONLY: io
    USE env_char_ops, ONLY: env_lowercase

    USE environ_param, ONLY: DP, e2, pi, tpi, BOHR_RADIUS_ANGS
    !
    USE class_cell
    USE class_density
    USE class_function
    USE class_function_gaussian
    USE class_functions
    !
    USE class_iontype
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
    TYPE, PUBLIC :: environ_ions
        !--------------------------------------------------------------------------------
        !
        LOGICAL :: lupdate = .FALSE.
        !
        INTEGER :: number = 0
        REAL(DP) :: center(3) = 0.D0
        !
        !--------------------------------------------------------------------------------
        ! Specifications of point-like ions
        !
        INTEGER :: ntyp
        INTEGER, ALLOCATABLE :: ityp(:)
        REAL(DP), POINTER :: tau(:, :) => NULL()
        TYPE(environ_iontype), ALLOCATABLE :: iontype(:)
        !
        !--------------------------------------------------------------------------------
        ! Parameters of the fictitious gaussian ionic density
        ! needed by electrostatic calculations
        !
        LOGICAL :: use_smeared_ions = .FALSE.
        CLASS(environ_function), ALLOCATABLE :: smeared_ions(:)
        TYPE(environ_density) :: density
        !
        !--------------------------------------------------------------------------------
        ! Parameters of the density of core electrons
        !
        LOGICAL :: use_core_electrons = .FALSE.
        CLASS(environ_function), ALLOCATABLE :: core_electrons(:)
        TYPE(environ_density) :: core
        !
        REAL(DP) :: charge = 0.D0
        REAL(DP) :: dipole(3) = 0.D0
        REAL(DP) :: quadrupole_pc(3) = 0.D0
        REAL(DP) :: quadrupole_gauss(3) = 0.D0
        REAL(DP) :: quadrupole_correction = 0.D0
        REAL(DP) :: selfenergy_correction = 0.D0
        !
        REAL(DP) :: potential_shift ! due to Gaussian-spread description (if used)
        ! #TODO set to zero?
        !
        !--------------------------------------------------------------------------------
    CONTAINS
        !--------------------------------------------------------------------------------
        !
        PROCEDURE, PRIVATE :: create => create_environ_ions
        PROCEDURE :: init => init_environ_ions
        PROCEDURE :: update => update_environ_ions
        PROCEDURE :: destroy => destroy_environ_ions
        !
        PROCEDURE :: convert_iontype_to_ion_array_char
        PROCEDURE :: convert_iontype_to_ion_array_integer
        PROCEDURE :: convert_iontype_to_ion_array_real
        !
        GENERIC :: get_iontype_array => &
            convert_iontype_to_ion_array_char, &
            convert_iontype_to_ion_array_integer, &
            convert_iontype_to_ion_array_real
        !
        PROCEDURE :: printout => print_environ_ions
        PROCEDURE :: write_cube => write_cube_ions
        !
        !--------------------------------------------------------------------------------
    END TYPE environ_ions
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
    SUBROUTINE create_environ_ions(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_ions), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: sub_name = 'create_environ_ions'
        !
        !--------------------------------------------------------------------------------
        !
        IF (ASSOCIATED(this%tau)) CALL io%create_error(sub_name)
        !
        IF (ALLOCATED(this%ityp)) CALL io%create_error(sub_name)
        !
        IF (ALLOCATED(this%iontype)) CALL io%create_error(sub_name)
        !
        IF (ALLOCATED(this%smeared_ions)) CALL io%create_error(sub_name)
        !
        IF (ALLOCATED(this%core_electrons)) CALL io%create_error(sub_name)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE create_environ_ions
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_environ_ions(this, nat, ntyp, atom_label, ityp, zv, atomicspread, &
                                 corespread, solvationrad, radius_mode, lsoftcavity, &
                                 lsmearedions, lcoredensity, cell)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: nat, ntyp
        INTEGER, INTENT(IN) :: ityp(nat)
        REAL(DP), INTENT(IN) :: zv(ntyp)
        CHARACTER(LEN=3), INTENT(IN) :: atom_label(ntyp)
        REAL(DP), DIMENSION(ntyp), INTENT(IN) :: atomicspread, corespread, solvationrad
        CHARACTER(LEN=80), INTENT(IN) :: radius_mode
        LOGICAL, INTENT(IN) :: lsoftcavity, lsmearedions, lcoredensity
        TYPE(environ_cell), INTENT(IN) :: cell
        !
        CLASS(environ_ions), INTENT(INOUT) :: this
        !
        INTEGER :: i
        !
        CHARACTER(LEN=80) :: local_label
        !
        CHARACTER(LEN=80) :: sub_name = 'init_environ_ions'
        !
        !--------------------------------------------------------------------------------
        !
        CALL this%create()
        !
        this%number = nat
        this%ntyp = ntyp
        !
        ALLOCATE (this%ityp(nat))
        this%ityp = ityp
        !
        ALLOCATE (this%tau(3, nat))
        this%tau = 0.D0 ! ionic positions assigned at beginning of each ionic step
        !
        !--------------------------------------------------------------------------------
        ! Ion types
        !
        ALLOCATE (this%iontype(ntyp))
        !
        DO i = 1, ntyp
            !
            CALL this%iontype(i)%init(i, atom_label(i), zv(i), radius_mode, &
                                      atomicspread(i), corespread(i), solvationrad(i), &
                                      lsoftcavity, lsmearedions)
            !
        END DO
        !
        !--------------------------------------------------------------------------------
        ! Total ionic charge
        !
        DO i = 1, this%number
            this%charge = this%charge + this%iontype(this%ityp(i))%zv
        END DO
        !
        !--------------------------------------------------------------------------------
        ! Smeared ions
        !
        IF (lsmearedions) THEN
            local_label = 'smeared_ions'
            !
            CALL this%density%init(cell, local_label)
            !
            !----------------------------------------------------------------------------
            ! Build smeared ions from iontype data
            !
            ALLOCATE (environ_function_gaussian :: this%smeared_ions(this%number))
            !
            DO i = 1, this%number
                !
                CALL this%smeared_ions(i)%init( &
                    1, 1, 0, 0.0_DP, this%iontype(this%ityp(i))%atomicspread, &
                    this%iontype(this%ityp(i))%zv, this%tau(:, i))
                !
            END DO
            !
        END IF
        !
        this%use_smeared_ions = lsmearedions
        !
        !--------------------------------------------------------------------------------
        ! Core electrons
        !
        IF (lcoredensity) THEN
            local_label = 'core_electrons'
            !
            CALL this%core%init(cell, local_label)
            !
            !----------------------------------------------------------------------------
            ! Build core electrons from iontype data
            !
            ALLOCATE (environ_function_gaussian :: this%core_electrons(this%number))
            !
            DO i = 1, this%number
                !
                IF (TRIM(this%iontype(this%ityp(i))%label) == 'H') &
                    this%iontype(this%ityp(i))%corespread = 1.D-10
                !
                CALL this%core_electrons(i)%init( &
                    1, 1, 0, 0.0_DP, this%iontype(this%ityp(i))%corespread, &
                    -this%iontype(this%ityp(i))%zv, this%tau(:, i))
                !
            END DO
            !
        END IF
        !
        this%use_core_electrons = lcoredensity
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_environ_ions
    !------------------------------------------------------------------------------------
    !>
    !! Update ionic positions and compute derived quantities
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE update_environ_ions(this, nat, tau, center)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: nat
        REAL(DP), INTENT(IN) :: tau(3, nat)
        REAL(DP), INTENT(IN), OPTIONAL :: center(3)
        !
        CLASS(environ_ions), INTENT(INOUT) :: this
        !
        INTEGER :: i
        INTEGER :: dim, axis
        REAL(DP) :: charge, spread
        REAL(DP) :: pos(3)
        !
        CHARACTER(LEN=80) :: sub_name = 'update_environ_ions'
        !
        !--------------------------------------------------------------------------------
        !
        IF (this%number /= nat) CALL io%error(sub_name, 'Mismatch in number of atoms', 1)
        !
        this%tau = tau ! update positions
        !
        !--------------------------------------------------------------------------------
        ! Center of ionic charge used by three sub-modules
        !
        IF (PRESENT(center)) THEN
            this%center = center
        ELSE
            this%center = 0.D0
            !
            DO i = 1, this%number
                !
                this%center(:) = this%center(:) + &
                                 this%tau(:, i) * this%iontype(this%ityp(i))%zv
                !
            END DO
            !
            IF (ABS(this%charge) < 1.D-8) &
                CALL io%error(sub_name, 'Ionic charge equal to zero', 1)
            !
            this%center = this%center / this%charge
        END IF
        !
        !--------------------------------------------------------------------------------
        ! If needed, generate a fictitious ion density using gaussians
        !
        IF (this%use_smeared_ions) &
            CALL density_of_functions(this%smeared_ions, this%number, &
                                      this%density, .TRUE.)
        !
        !--------------------------------------------------------------------------------
        ! Compute quadrupole moment of point-like (and gaussian) nuclei
        !
        this%dipole = 0.D0 ! this is due to the choice of ionic center
        this%quadrupole_pc = 0.D0
        this%quadrupole_correction = 0.D0
        this%selfenergy_correction = 0.D0
        !
        DO i = 1, this%number
            !
            this%quadrupole_pc(:) = this%quadrupole_pc(:) + &
                                    this%iontype(this%ityp(i))%zv * &
                                    ((this%tau(:, i) - this%center(:)))**2
            !
            IF (this%use_smeared_ions) THEN
                !
                this%quadrupole_correction = &
                    this%quadrupole_correction + this%iontype(this%ityp(i))%zv * &
                    this%iontype(this%ityp(i))%atomicspread**2 * 0.5D0
                !
                this%selfenergy_correction = &
                    this%selfenergy_correction + this%iontype(this%ityp(i))%zv**2 / &
                    this%iontype(this%ityp(i))%atomicspread * SQRT(2.D0 / pi)
                !
            END IF
            !
        END DO
        !
        !--------------------------------------------------------------------------------
        ! Calculate potential shift due to Gaussian nuclei
        !
        IF (this%use_smeared_ions) THEN
            !
            this%potential_shift = this%quadrupole_correction * &
                                   tpi * e2 / this%density%cell%omega
            !
            this%quadrupole_gauss(:) = this%quadrupole_pc(:) + this%quadrupole_correction
        END IF
        !
        !--------------------------------------------------------------------------------
        ! Output current state
        !
        CALL this%printout()
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE update_environ_ions
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE destroy_environ_ions(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_ions), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: sub_name = 'destroy_environ_ions'
        !
        !--------------------------------------------------------------------------------
        !
        IF (this%use_smeared_ions) THEN
            !
            CALL this%density%destroy()
            !
            CALL destroy_environ_functions(this%smeared_ions, this%number)
            !
        END IF
        !
        IF (this%use_core_electrons) THEN
            !
            CALL this%core%destroy()
            !
            CALL destroy_environ_functions(this%core_electrons, this%number)
            !
        END IF
        !
        this%charge = 0.D0
        !
        IF (.NOT. ALLOCATED(this%ityp)) CALL io%destroy_error(sub_name)
        !
        IF (.NOT. ALLOCATED(this%iontype)) CALL io%destroy_error(sub_name)
        !
        IF (.NOT. ASSOCIATED(this%tau)) CALL io%destroy_error(sub_name)
        !
        DEALLOCATE (this%ityp)
        DEALLOCATE (this%iontype)
        DEALLOCATE (this%tau)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE destroy_environ_ions
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
    SUBROUTINE convert_iontype_to_ion_array_char(this, array)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_ions), INTENT(INOUT) :: this
        CHARACTER(LEN=3), ALLOCATABLE, INTENT(INOUT) :: array(:)
        !
        INTEGER :: i
        !
        CHARACTER(LEN=80) :: sub_name = 'convert_iontype_to_ion_array_char'
        !
        !--------------------------------------------------------------------------------
        !
        ALLOCATE (array(this%number))
        !
        DO i = 1, this%number
            array(i) = this%iontype(this%ityp(i))%label
        END DO
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE convert_iontype_to_ion_array_char
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE convert_iontype_to_ion_array_integer(this, array, item)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CHARACTER(LEN=20), INTENT(IN) :: item
        !
        CLASS(environ_ions), INTENT(INOUT) :: this
        INTEGER, ALLOCATABLE, INTENT(INOUT) :: array(:)
        !
        INTEGER :: i
        !
        CHARACTER(LEN=80) :: sub_name = 'convert_iontype_to_ion_array_integer'
        !
        !--------------------------------------------------------------------------------
        !
        ALLOCATE (array(this%number))
        !
        SELECT CASE (TRIM(item))
            !
        CASE ('index')
            !
            DO i = 1, this%number
                array(i) = this%iontype(this%ityp(i))%index
            END DO
            !
        CASE ('atnum')
            !
            DO i = 1, this%number
                array(i) = this%iontype(this%ityp(i))%atmnum
            END DO
            !
        CASE DEFAULT
            CALL io%error(sub_name, 'Unexpected keyword', 1)
            !
        END SELECT
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE convert_iontype_to_ion_array_integer
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE convert_iontype_to_ion_array_real(this, array, item)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CHARACTER(LEN=20), INTENT(IN) :: item
        !
        CLASS(environ_ions), INTENT(INOUT) :: this
        REAL(DP), ALLOCATABLE, INTENT(INOUT) :: array(:)
        !
        INTEGER :: i
        !
        CHARACTER(LEN=80) :: sub_name = 'convert_iontype_to_ion_array_real'
        !
        !--------------------------------------------------------------------------------
        !
        ALLOCATE (array(this%number))
        !
        SELECT CASE (TRIM(item))
            !
        CASE ('zv')
            !
            DO i = 1, this%number
                array(i) = this%iontype(this%ityp(i))%zv
            END DO
            !
        CASE ('atomicspread')
            !
            DO i = 1, this%number
                array(i) = this%iontype(this%ityp(i))%atomicspread
            END DO
            !
        CASE ('corespread')
            !
            DO i = 1, this%number
                array(i) = this%iontype(this%ityp(i))%corespread
            END DO
            !
        CASE ('solvationrad')
            !
            DO i = 1, this%number
                array(i) = this%iontype(this%ityp(i))%solvationrad
            END DO
            !
        CASE DEFAULT
            CALL io%error(sub_name, 'Unexpected keyword', 1)
            !
        END SELECT
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE convert_iontype_to_ion_array_real
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !                                   OUTPUT METHODS
    !
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !>
    !! Prints the details of the ions
    !!
    !! Nested objects receive a decremented passed verbose to trigger block printing
    !!
    !! @param verbose       : (INTEGER) adds verbosity to global verbose
    !! @param debug_verbose : (INTEGER) replaces global verbose for debugging
    !! @param unit          : (INTEGER) output target (default = io%debug_unit)
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE print_environ_ions(this, verbose, debug_verbose, unit)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_ions), INTENT(IN) :: this
        INTEGER, INTENT(IN), OPTIONAL :: verbose, debug_verbose, unit
        !
        INTEGER :: base_verbose, local_verbose, passed_verbose, local_unit, i
        !
        CHARACTER(LEN=80) :: sub_name = 'print_environ_ions'
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
        IF (local_verbose >= 1) THEN
            !
            IF (io%lnode) THEN
                WRITE (local_unit, 1000)
                !
                WRITE (local_unit, 1001) &
                    this%charge, this%center, this%dipole, this%quadrupole_pc
                !
                IF (this%use_smeared_ions) &
                    WRITE (local_unit, 1002) this%quadrupole_gauss
                !
                WRITE (local_unit, 1003)
                !
                DO i = 1, this%number
                    WRITE (local_unit, 1004) i, this%ityp(i), this%tau(:, i)
                END DO
                !
            END IF
            !
            IF (local_verbose >= 3) THEN
                !
                CALL print_environ_iontypes(this%iontype, this%ntyp, passed_verbose, &
                                            debug_verbose, local_unit)
                !
                IF (this%use_smeared_ions) THEN
                    !
                    CALL this%density%printout(passed_verbose, debug_verbose, local_unit)
                    !
                    IF (local_verbose >= 4) &
                        CALL print_environ_functions(this%smeared_ions, this%number, &
                                                     passed_verbose, debug_verbose, &
                                                     local_unit)
                    !
                END IF
                !
                IF (this%use_core_electrons) THEN
                    !
                    IF (local_verbose >= 4) &
                        CALL this%core%printout(passed_verbose, debug_verbose, &
                                                local_unit)
                    !
                    IF (local_verbose >= 5) &
                        CALL print_environ_functions(this%core_electrons, this%number, &
                                                     passed_verbose, debug_verbose, &
                                                     local_unit)
                    !
                END IF
                !
            END IF
            !
        END IF
        !
        FLUSH (local_unit)
        !
        !--------------------------------------------------------------------------------
        !
1000    FORMAT(/, 4('%'), ' IONS ', 70('%'))
        !
1001    FORMAT(/, ' total charge               = ', F14.7, /, &
                ' center of charge           = ', 3F14.7, /, &
                ' dipole                     = ', 3F14.7, /, &
                ' quadrupole (pc)            = ', 3F14.7)
        !
1002    FORMAT(' quadrupole (gauss)         = ', 3F14.7)
        !
1003    FORMAT(/, '   i | type | coordinates', /, 1X, 71('-'))
        !
1004    FORMAT(1X, I3, ' | ', I4, ' |                 ', 3F14.7)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE print_environ_ions
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE write_cube_ions(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_ions), TARGET, INTENT(IN) :: this
        !
        INTEGER :: iat
        !
        !--------------------------------------------------------------------------------
        ! Write cube cell data
        !
        DO iat = 1, this%number
            !
            WRITE (300, '(i5,4f12.6)') &
                this%iontype(this%ityp(iat))%atmnum, 0.D0, &
                this%tau(1, iat), this%tau(2, iat), this%tau(3, iat)
            !
        END DO
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE write_cube_ions
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE class_ions
!----------------------------------------------------------------------------------------
