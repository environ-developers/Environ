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
!! Module containing the main routines to handle environ_charges
!! derived data types.
!!
!! Environ_charges are used as wrappers for electrostatic charges to be
!! passed to electrostatic drivers. Charges may include source charges
!! from the QM system (electrons and ions), but also externally-defined
!! charges, dielectric polarization charges and electrolyte charges.
!!
!----------------------------------------------------------------------------------------
MODULE class_charges
    !------------------------------------------------------------------------------------
    !
    USE env_base_io, ONLY: ionode, environ_unit, verbose, depth
    !
    USE environ_param, ONLY: DP
    !
    USE class_cell
    USE class_density
    !
    USE class_dielectric
    USE class_electrons
    USE class_electrolyte
    USE class_externals
    USE class_ions
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
    TYPE, PUBLIC :: environ_charges
        !--------------------------------------------------------------------------------
        !
        LOGICAL :: include_ions
        TYPE(environ_ions), POINTER :: ions => NULL()
        !
        LOGICAL :: include_electrons
        TYPE(environ_electrons), POINTER :: electrons => NULL()
        !
        LOGICAL :: include_externals
        TYPE(environ_externals), POINTER :: externals => NULL()
        !
        LOGICAL :: include_dielectric
        TYPE(environ_dielectric), POINTER :: dielectric => NULL()
        !
        LOGICAL :: include_electrolyte
        TYPE(environ_electrolyte), POINTER :: electrolyte => NULL()
        !
        LOGICAL :: include_semiconductor
        TYPE(environ_semiconductor), POINTER :: semiconductor => NULL()
        !
        !--------------------------------------------------------------------------------
        ! Total smooth free charge
        !
        INTEGER :: number
        REAL(DP) :: charge
        TYPE(environ_density) :: density
        !
        !--------------------------------------------------------------------------------
    CONTAINS
        !--------------------------------------------------------------------------------
        !
        PROCEDURE, PRIVATE :: create => create_environ_charges
        PROCEDURE :: init => init_environ_charges
        PROCEDURE :: update => update_environ_charges
        PROCEDURE :: destroy => destroy_environ_charges
        !
        PROCEDURE :: add => add_environ_charges
        !
        PROCEDURE :: of_potential => charges_of_potential
        !
        PROCEDURE :: printout => print_environ_charges
        !
        !--------------------------------------------------------------------------------
    END TYPE environ_charges
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
    SUBROUTINE create_environ_charges(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_charges), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: sub_name = 'create_environ_charges'
        !
        !--------------------------------------------------------------------------------
        !
        IF (ASSOCIATED(this%ions)) &
            CALL env_errore(sub_name, 'Trying to create an existing object', 1)
        !
        IF (ASSOCIATED(this%electrons)) &
            CALL env_errore(sub_name, 'Trying to create an existing object', 1)
        !
        IF (ASSOCIATED(this%externals)) &
            CALL env_errore(sub_name, 'Trying to create an existing object', 1)
        !
        IF (ASSOCIATED(this%dielectric)) &
            CALL env_errore(sub_name, 'Trying to create an existing object', 1)
        !
        IF (ASSOCIATED(this%electrolyte)) &
            CALL env_errore(sub_name, 'Trying to create an existing object', 1)
        !
        IF (ASSOCIATED(this%semiconductor)) &
            CALL env_errore(sub_name, 'Trying to create an existing object', 1)
        !
        IF (ALLOCATED(this%density%of_r)) &
            CALL env_errore(sub_name, 'Trying to create an existing object', 1)
        !
        !--------------------------------------------------------------------------------
        !
        this%include_ions = .FALSE.
        NULLIFY (this%ions)
        !
        this%include_electrons = .FALSE.
        NULLIFY (this%electrons)
        !
        this%include_externals = .FALSE.
        NULLIFY (this%externals)
        !
        this%include_dielectric = .FALSE.
        NULLIFY (this%dielectric)
        !
        this%include_electrolyte = .FALSE.
        NULLIFY (this%electrolyte)
        !
        this%include_semiconductor = .FALSE.
        NULLIFY (this%semiconductor)
        !
        this%number = 0
        this%charge = 0.D0
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE create_environ_charges
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_environ_charges(this, cell)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_cell), INTENT(IN) :: cell
        !
        CLASS(environ_charges), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: local_label = 'charge'
        !
        !--------------------------------------------------------------------------------
        !
        CALL this%create()
        !
        CALL this%density%init(cell, local_label)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_environ_charges
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE update_environ_charges(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_charges), INTENT(INOUT) :: this
        !
        REAL(DP) :: local_charge
        !
        CHARACTER(LEN=80) :: sub_name = 'update_environ_charges'
        !
        !--------------------------------------------------------------------------------
        !
        this%number = 0
        this%charge = 0.D0
        this%density%of_r = 0.D0
        !
        IF (this%include_electrons) THEN
            !
            IF (.NOT. ASSOCIATED(this%electrons)) &
                CALL env_errore(sub_name, 'Missing expected charge component', 1)
            !
            this%number = this%number + this%electrons%number
            this%charge = this%charge + this%electrons%charge
            this%density%of_r = this%density%of_r + this%electrons%density%of_r
        END IF
        !
        IF (this%include_ions) THEN
            !
            IF (.NOT. ASSOCIATED(this%ions)) &
                CALL env_errore(sub_name, 'Missing expected charge component', 1)
            !
            this%number = this%number + this%ions%number
            this%charge = this%charge + this%ions%charge
            this%density%of_r = this%density%of_r + this%ions%density%of_r
        END IF
        !
        IF (this%include_externals) THEN
            !
            IF (.NOT. ASSOCIATED(this%externals)) &
                CALL env_errore(sub_name, 'Missing expected charge component', 1)
            !
            this%number = this%number + this%externals%number
            this%charge = this%charge + this%externals%charge
            this%density%of_r = this%density%of_r + this%externals%density%of_r
        END IF
        !
        CALL this%density%update()
        !
        local_charge = this%density%charge
        !
        IF (ABS(local_charge - this%charge) > 1.D-5) &
            CALL env_errore(sub_name, 'Inconsistent integral of total charge', 1)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE update_environ_charges
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE add_environ_charges(this, electrons, ions, externals, &
                                   dielectric, electrolyte, semiconductor)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_electrons), OPTIONAL, TARGET, INTENT(IN) :: electrons
        TYPE(environ_ions), OPTIONAL, TARGET, INTENT(IN) :: ions
        TYPE(environ_externals), OPTIONAL, TARGET, INTENT(IN) :: externals
        TYPE(environ_dielectric), OPTIONAL, TARGET, INTENT(IN) :: dielectric
        TYPE(environ_electrolyte), OPTIONAL, TARGET, INTENT(IN) :: electrolyte
        TYPE(environ_semiconductor), OPTIONAL, TARGET, INTENT(IN) :: semiconductor
        !
        CLASS(environ_charges), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: sub_name = 'add_environ_charges'
        !
        !--------------------------------------------------------------------------------
        !
        IF (PRESENT(ions)) THEN
            this%include_ions = .TRUE.
            this%ions => ions
        END IF
        !
        IF (PRESENT(electrons)) THEN
            this%include_electrons = .TRUE.
            this%electrons => electrons
        END IF
        !
        IF (PRESENT(externals)) THEN
            this%include_externals = .TRUE.
            this%externals => externals
        END IF
        !
        IF (PRESENT(dielectric)) THEN
            this%include_dielectric = .TRUE.
            this%dielectric => dielectric
        END IF
        !
        IF (PRESENT(electrolyte)) THEN
            this%include_electrolyte = .TRUE.
            this%electrolyte => electrolyte
        END IF
        !
        IF (PRESENT(semiconductor)) THEN
            this%include_semiconductor = .TRUE.
            this%semiconductor => semiconductor
        END IF
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE add_environ_charges
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE destroy_environ_charges(this, lflag)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        LOGICAL, INTENT(IN) :: lflag
        !
        CLASS(environ_charges), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: sub_name = 'destroy_environ_charges'
        !
        !--------------------------------------------------------------------------------
        !
        IF (lflag) THEN
            !
            IF (ASSOCIATED(this%ions)) NULLIFY (this%ions)
            !
            IF (ASSOCIATED(this%electrons)) NULLIFY (this%electrons)
            !
            IF (ASSOCIATED(this%externals)) NULLIFY (this%externals)
            !
            IF (ASSOCIATED(this%dielectric)) NULLIFY (this%dielectric)
            !
            IF (ASSOCIATED(this%electrolyte)) NULLIFY (this%electrolyte)
            !
        END IF
        !
        CALL this%density%destroy()
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE destroy_environ_charges
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
    SUBROUTINE charges_of_potential(this, potential)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_density), INTENT(IN) :: potential
        !
        CLASS(environ_charges), INTENT(INOUT) :: this
        !
        TYPE(environ_density) :: tot_charge_density
        !
        CHARACTER(LEN=80) :: sub_name = 'charges_of_potential'
        !
        !--------------------------------------------------------------------------------
        !
        CALL tot_charge_density%init(potential%cell)
        !
        tot_charge_density%of_r = this%density%of_r
        !
        IF (this%include_electrolyte) THEN
            !
            IF (.NOT. ASSOCIATED(this%electrolyte)) &
                CALL env_errore(sub_name, 'Missing expected charge component', 1)
            !
            CALL this%electrolyte%of_potential(potential)
            !
            !----------------------------------------------------------------------------
            ! The electrolyte charges are required in the total charge for the
            ! calculation of the dielectric polarization
            !
            tot_charge_density%of_r = tot_charge_density%of_r + &
                                      this%electrolyte%density%of_r
            !
        END IF
        !
        IF (this%include_dielectric) THEN
            !
            IF (.NOT. ASSOCIATED(this%dielectric)) &
                CALL env_errore(sub_name, 'Missing expected charge component', 1)
            !
            CALL this%dielectric%of_potential(tot_charge_density, potential)
            !
        END IF
        !
        CALL tot_charge_density%destroy()
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE charges_of_potential
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
    SUBROUTINE print_environ_charges(this, local_verbose, local_depth)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_charges), INTENT(IN) :: this
        INTEGER, INTENT(IN), OPTIONAL :: local_verbose
        INTEGER, INTENT(IN), OPTIONAL :: local_depth
        !
        INTEGER :: verbosity, passed_verbosity, passed_depth
        !
        CHARACTER(LEN=80) :: sub_name = 'print_environ_charges'
        !
        !--------------------------------------------------------------------------------
        !
        IF (.NOT. ionode .OR. verbose == 0) RETURN
        !
        IF (PRESENT(local_verbose)) THEN
            verbosity = verbose + local_verbose
        ELSE
            verbosity = verbose
        END IF
        !
        IF (verbosity == 0) RETURN
        !
        IF (verbosity >= 1) THEN
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
            WRITE (environ_unit, 1002) this%number
            WRITE (environ_unit, 1003) this%charge
            !
            IF (verbosity < verbose) CALL env_block_divider(verbosity)
            !
        END IF
        !
        FLUSH (environ_unit)
        !
        !--------------------------------------------------------------------------------
        !
1000    FORMAT(/, 4('%'), ' CHARGES ', 67('%'))
1001    FORMAT(/, ' CHARGES', /, ' =======')
        !
1002    FORMAT(/, ' total number of charges    = ', I10)
        !
1003    FORMAT(/, ' total charge               = ', F14.7)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE print_environ_charges
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE class_charges
!----------------------------------------------------------------------------------------
