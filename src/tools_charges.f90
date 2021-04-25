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
MODULE tools_charges
    !------------------------------------------------------------------------------------
    !
    USE physical_types, ONLY: environ_charges
    USE representation_types, ONLY: environ_density
    !
    USE utils_density, ONLY: init_environ_density, destroy_environ_density
    !
    USE tools_dielectric, ONLY: dielectric_of_potential
    USE tools_electrolyte, ONLY: electrolyte_of_potential
    !
    !------------------------------------------------------------------------------------
CONTAINS
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE charges_of_potential(potential, charges)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_density), INTENT(IN) :: potential
        !
        TYPE(environ_charges), INTENT(INOUT) :: charges
        !
        TYPE(environ_density) :: tot_charge_density
        !
        CHARACTER(LEN=80) :: sub_name = 'charges_of_potential'
        !
        !--------------------------------------------------------------------------------
        !
        CALL init_environ_density(potential%cell, tot_charge_density)
        !
        tot_charge_density%of_r = charges%density%of_r
        !
        IF (charges%include_electrolyte) THEN
            !
            IF (.NOT. ASSOCIATED(charges%electrolyte)) &
                CALL errore(sub_name, 'Missing expected charge component', 1)
            !
            CALL electrolyte_of_potential(potential, charges%electrolyte)
            !
            !----------------------------------------------------------------------------
            ! The electrolyte charges are required in the total charge for the
            ! calculation of the dielectric polarization
            !
            tot_charge_density%of_r = tot_charge_density%of_r + &
                                      charges%electrolyte%density%of_r
            !
        END IF
        !
        IF (charges%include_dielectric) THEN
            !
            IF (.NOT. ASSOCIATED(charges%dielectric)) &
                CALL errore(sub_name, 'Missing expected charge component', 1)
            !
            CALL dielectric_of_potential(tot_charge_density, potential, &
                                         charges%dielectric)
            !
        END IF
        !
        CALL destroy_environ_density(tot_charge_density)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE charges_of_potential
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE tools_charges
!----------------------------------------------------------------------------------------
