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
!
! Authors: Oliviero Andreussi (Department of Physics, UNT)
!          Francesco Nattino  (THEOS and NCCR-MARVEL, EPFL)
!          Nicola Marzari     (THEOS and NCCR-MARVEL, EPFL)
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
MODULE utils_charges
    !------------------------------------------------------------------------------------
    !
    USE environ_param, ONLY: DP
    !
    USE types_representation, ONLY: environ_density
    USE types_cell, ONLY: environ_cell
    !
    USE types_physical, ONLY: environ_charges, environ_electrons, environ_ions, &
                              environ_externals, environ_dielectric, &
                              environ_electrolyte, environ_semiconductor
    !
    USE utils_density, ONLY: create_environ_density, init_environ_density, &
                             update_environ_density, destroy_environ_density
    !
    !------------------------------------------------------------------------------------
    !
    PRIVATE
    !
    PUBLIC :: create_environ_charges, init_environ_charges_first, &
              init_environ_charges_second, update_environ_charges, &
              destroy_environ_charges
    !
    !------------------------------------------------------------------------------------
CONTAINS
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE create_environ_charges(charges)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_charges) :: charges
        !
        CHARACTER(LEN=80) :: label = 'charges'
        !
        !--------------------------------------------------------------------------------
        !
        charges%include_ions = .FALSE.
        NULLIFY (charges%ions)
        !
        charges%include_electrons = .FALSE.
        NULLIFY (charges%electrons)
        !
        charges%include_externals = .FALSE.
        NULLIFY (charges%externals)
        !
        charges%include_dielectric = .FALSE.
        NULLIFY (charges%dielectric)
        !
        charges%include_electrolyte = .FALSE.
        NULLIFY (charges%electrolyte)
        !
        charges%include_semiconductor = .FALSE.
        NULLIFY (charges%semiconductor)
        !
        charges%number = 0
        charges%charge = 0.D0
        !
        CALL create_environ_density(charges%density, label)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE create_environ_charges
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_environ_charges_first(charges, electrons, ions, externals, &
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
        TYPE(environ_charges), INTENT(INOUT) :: charges
        !
        !--------------------------------------------------------------------------------
        !
        IF (PRESENT(ions)) THEN
            charges%include_ions = .TRUE.
            charges%ions => ions
        END IF
        !
        IF (PRESENT(electrons)) THEN
            charges%include_electrons = .TRUE.
            charges%electrons => electrons
        END IF
        !
        IF (PRESENT(externals)) THEN
            charges%include_externals = .TRUE.
            charges%externals => externals
        END IF
        !
        IF (PRESENT(dielectric)) THEN
            charges%include_dielectric = .TRUE.
            charges%dielectric => dielectric
        END IF
        !
        IF (PRESENT(electrolyte)) THEN
            charges%include_electrolyte = .TRUE.
            charges%electrolyte => electrolyte
        END IF
        !
        IF (PRESENT(semiconductor)) THEN
            charges%include_semiconductor = .TRUE.
            charges%semiconductor => semiconductor
        END IF
        !
        charges%initialized = .FALSE.
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_environ_charges_first
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_environ_charges_second(cell, charges)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_cell), INTENT(IN) :: cell
        !
        TYPE(environ_charges), INTENT(INOUT) :: charges
        !
        !--------------------------------------------------------------------------------
        !
        CALL init_environ_density(cell, charges%density)
        !
        charges%initialized = .TRUE.
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_environ_charges_second
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE update_environ_charges(charges)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_charges), INTENT(INOUT) :: charges
        !
        REAL(DP) :: local_charge
        !
        CHARACTER(LEN=80) :: sub_name = 'update_environ_charges'
        !
        !--------------------------------------------------------------------------------
        !
        charges%number = 0
        charges%charge = 0.D0
        charges%density%of_r = 0.D0
        !
        IF (charges%include_electrons) THEN
            !
            IF (.NOT. ASSOCIATED(charges%electrons)) &
                CALL env_errore(sub_name, 'Missing expected charge component', 1)
            !
            charges%number = charges%number + charges%electrons%number
            charges%charge = charges%charge + charges%electrons%charge
            !
            charges%density%of_r = charges%density%of_r + &
                                   charges%electrons%density%of_r
            !
        END IF
        !
        IF (charges%include_ions) THEN
            !
            IF (.NOT. ASSOCIATED(charges%ions)) &
                CALL env_errore(sub_name, 'Missing expected charge component', 1)
            !
            charges%number = charges%number + charges%ions%number
            charges%charge = charges%charge + charges%ions%charge
            charges%density%of_r = charges%density%of_r + charges%ions%density%of_r
        END IF
        !
        IF (charges%include_externals) THEN
            !
            IF (.NOT. ASSOCIATED(charges%externals)) &
                CALL env_errore(sub_name, 'Missing expected charge component', 1)
            !
            charges%number = charges%number + charges%externals%number
            charges%charge = charges%charge + charges%externals%charge
            !
            charges%density%of_r = charges%density%of_r + &
                                   charges%externals%density%of_r
            !
        END IF
        !
        CALL update_environ_density(charges%density)
        !
        local_charge = charges%density%charge
        !
        IF (ABS(local_charge - charges%charge) > 1.D-5) &
            CALL env_errore(sub_name, 'Inconsistent integral of total charge', 1)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE update_environ_charges
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE destroy_environ_charges(lflag, charges)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        LOGICAL, INTENT(IN) :: lflag
        !
        TYPE(environ_charges) :: charges
        !
        CHARACTER(LEN=80) :: sub_name = 'destroy_environ_charges'
        !
        !--------------------------------------------------------------------------------
        !
        IF (lflag) THEN
            !
            IF (ASSOCIATED(charges%ions)) NULLIFY (charges%ions)
            !
            IF (ASSOCIATED(charges%electrons)) NULLIFY (charges%electrons)
            !
            IF (ASSOCIATED(charges%externals)) NULLIFY (charges%externals)
            !
            IF (ASSOCIATED(charges%dielectric)) NULLIFY (charges%dielectric)
            !
            IF (ASSOCIATED(charges%electrolyte)) NULLIFY (charges%electrolyte)
            !
        END IF
        !
        IF (charges%initialized) THEN
            !
            CALL destroy_environ_density(charges%density)
            !
            charges%initialized = .FALSE.
        END IF
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE destroy_environ_charges
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE utils_charges
!----------------------------------------------------------------------------------------
