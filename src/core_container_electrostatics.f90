!----------------------------------------------------------------------------------------
!
! Copyright (C) 2021 ENVIRON (www.quantum-environ.org)
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
! Authors: Edan Bainglass (Department of Physics, UNT)
!
!----------------------------------------------------------------------------------------
!>
!!
!----------------------------------------------------------------------------------------
MODULE class_core_container_electrostatics
    !------------------------------------------------------------------------------------
    !
    USE env_base_io, ONLY: io
    !
    USE environ_param, ONLY: DP
    !
    USE class_density
    USE class_gradient
    !
    USE class_core_container_corrections
    USE class_core_container_derivatives
    USE class_core_fft_electrostatics
    !
    USE class_ions
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
    TYPE, EXTENDS(container_derivatives), PUBLIC :: container_electrostatics
        !--------------------------------------------------------------------------------
        !
        CLASS(container_corrections), POINTER :: correction => NULL()
        !
        !--------------------------------------------------------------------------------
    CONTAINS
        !--------------------------------------------------------------------------------
        !
        PROCEDURE :: add_correction
        PROCEDURE :: destroy => destroy_electrostatics_container
        !
        PROCEDURE :: poisson => calc_poisson
        PROCEDURE :: gradpoisson => calc_gradpoisson
        PROCEDURE :: force => calc_force
        !
        !--------------------------------------------------------------------------------
    END TYPE container_electrostatics
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
    SUBROUTINE add_correction(this, correction)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(container_corrections), TARGET, INTENT(IN) :: correction
        !
        CLASS(container_electrostatics), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: sub_name = 'add_correction'
        !
        !--------------------------------------------------------------------------------
        !
        IF (ASSOCIATED(this%correction)) CALL io%create_error(sub_name)
        !
        !--------------------------------------------------------------------------------
        !
        this%correction => correction
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE add_correction
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE destroy_electrostatics_container(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(container_electrostatics), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: sub_name = 'destroy_electrostatics_container'
        !
        !--------------------------------------------------------------------------------
        !
        CALL this%core%destroy()
        !
        IF (.NOT. ASSOCIATED(this%core)) CALL io%destroy_error(sub_name)
        !
        NULLIFY (this%core)
        !
        IF (ASSOCIATED(this%correction)) THEN
            !
            CALL this%correction%destroy()
            !
            IF (.NOT. ASSOCIATED(this%correction)) CALL io%destroy_error(sub_name)
            !
            NULLIFY (this%correction)
        END IF
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE destroy_electrostatics_container
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !                               ELECTROSTATIC METHODS
    !
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_poisson(this, density, potential)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_density), INTENT(IN) :: density
        !
        CLASS(container_electrostatics), INTENT(INOUT) :: this
        TYPE(environ_density), INTENT(INOUT) :: potential
        !
        CHARACTER(LEN=80) :: sub_name = 'calc_poisson'
        !
        !--------------------------------------------------------------------------------
        !
        SELECT TYPE (core => this%core)
            !
        TYPE IS (core_fft_electrostatics)
            CALL core%poisson(density, potential)
            !
        CLASS DEFAULT
            CALL io%error(sub_name, 'Unexpected core', 1)
            !
        END SELECT
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_poisson
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_gradpoisson(this, density, gradient)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_density), INTENT(IN) :: density
        !
        CLASS(container_electrostatics), INTENT(INOUT) :: this
        TYPE(environ_gradient), INTENT(INOUT) :: gradient
        !
        CHARACTER(LEN=80) :: sub_name = 'calc_gradpoisson'
        !
        !--------------------------------------------------------------------------------
        !
        SELECT TYPE (core => this%core)
            !
        TYPE IS (core_fft_electrostatics)
            CALL core%gradpoisson(density, gradient)
            !
        CLASS DEFAULT
            CALL io%error(sub_name, 'Unexpected core', 1)
            !
        END SELECT
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_gradpoisson
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_force(this, natoms, density, ions, force)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: natoms
        TYPE(environ_density), INTENT(IN) :: density
        TYPE(environ_ions), TARGET, INTENT(IN) :: ions
        !
        CLASS(container_electrostatics), INTENT(INOUT) :: this
        REAL(DP), INTENT(INOUT) :: force(3, natoms)
        !
        CHARACTER(LEN=80) :: sub_name = 'calc_force'
        !
        !--------------------------------------------------------------------------------
        !
        SELECT TYPE (core => this%core)
            !
        TYPE IS (core_fft_electrostatics)
            CALL core%force(natoms, density, ions, force)
            !
        CLASS DEFAULT
            CALL io%error(sub_name, 'Unexpected core', 1)
            !
        END SELECT
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_force
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE class_core_container_electrostatics
!----------------------------------------------------------------------------------------
