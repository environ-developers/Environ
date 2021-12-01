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
! Authors: Edan Bainglass (Department of Physics, UNT)
!
!----------------------------------------------------------------------------------------
!>
!!
!----------------------------------------------------------------------------------------
MODULE class_container
    !------------------------------------------------------------------------------------
    !
    USE class_io, ONLY: io
    !
    USE class_core_container_derivatives
    USE class_core_container_electrostatics
    USE class_core_container_corrections
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
    TYPE, PUBLIC :: environ_container
        !--------------------------------------------------------------------------------
        !
        CHARACTER(LEN=80) :: label = ''
        !
        TYPE(container_derivatives), POINTER :: derivatives => NULL()
        TYPE(container_electrostatics), POINTER :: electrostatics => NULL()
        !
        !--------------------------------------------------------------------------------
    CONTAINS
        !--------------------------------------------------------------------------------
        !
        PROCEDURE, PRIVATE :: create => create_environ_container
        PROCEDURE :: init => init_environ_container
        PROCEDURE :: destroy => destroy_environ_container
        !
        !--------------------------------------------------------------------------------
    END TYPE environ_container
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
    SUBROUTINE create_environ_container(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_container), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: sub_name = 'create_environ_container'
        !
        !--------------------------------------------------------------------------------
        !
        IF (ASSOCIATED(this%derivatives)) CALL io%create_error(sub_name)
        !
        IF (ASSOCIATED(this%electrostatics)) CALL io%create_error(sub_name)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE create_environ_container
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_environ_container(this, label, derivatives, electrostatics)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CHARACTER(LEN=80), INTENT(IN) :: label
        TYPE(container_derivatives), TARGET, INTENT(IN), OPTIONAL :: derivatives
        TYPE(container_electrostatics), TARGET, INTENT(IN), OPTIONAL :: electrostatics
        !
        CLASS(environ_container), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: sub_name = 'init_environ_container'
        !
        !--------------------------------------------------------------------------------
        !
        CALL this%create()
        !
        this%label = label
        !
        IF (PRESENT(derivatives)) this%derivatives => derivatives
        !
        IF (PRESENT(electrostatics)) this%electrostatics => electrostatics
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_environ_container
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE destroy_environ_container(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_container), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: sub_name = 'destroy_environ_container'
        !
        !--------------------------------------------------------------------------------
        !
        IF (ASSOCIATED(this%derivatives)) CALL this%derivatives%destroy()
        !
        IF (ASSOCIATED(this%electrostatics)) CALL this%electrostatics%destroy()
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE destroy_environ_container
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE class_container
!----------------------------------------------------------------------------------------
