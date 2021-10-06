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
!! #TODO 1D_numeric, multigrid, bigdft
!!
!----------------------------------------------------------------------------------------
MODULE class_core_container
    !------------------------------------------------------------------------------------
    !
    USE class_io, ONLY: io
    !
    USE class_core_numerical
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
    TYPE, ABSTRACT, PUBLIC :: core_container
        !--------------------------------------------------------------------------------
        !
        CHARACTER(LEN=80) :: type_
        !
        CLASS(numerical_core), POINTER :: core => NULL()
        !
        !--------------------------------------------------------------------------------
    CONTAINS
        !--------------------------------------------------------------------------------
        !
        PROCEDURE, PRIVATE :: create => create_core_container
        PROCEDURE :: init => init_core_container
        PROCEDURE :: destroy => destroy_core_container
        !
        !--------------------------------------------------------------------------------
    END TYPE core_container
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
    SUBROUTINE create_core_container(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(core_container), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: sub_name = 'create_core_container'
        !
        !--------------------------------------------------------------------------------
        !
        IF (ASSOCIATED(this%core)) CALL io%create_error(sub_name)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE create_core_container
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_core_container(this, core, type_in)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(numerical_core), TARGET, INTENT(IN) :: core
        CHARACTER(LEN=80), INTENT(IN) :: type_in
        !
        CLASS(core_container), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: sub_name = 'init_core_container'
        !
        !--------------------------------------------------------------------------------
        !
        CALL this%create()
        !
        this%core => core
        this%type_ = type_in
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_core_container
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE destroy_core_container(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(core_container), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: sub_name = 'destroy_core_container'
        !
        !--------------------------------------------------------------------------------
        !
        IF (.NOT. ASSOCIATED(this%core)) CALL io%destroy_error(sub_name)
        !
        !--------------------------------------------------------------------------------
        !
        IF (ASSOCIATED(this%core%cell)) CALL this%core%destroy()
        !
        NULLIFY (this%core)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE destroy_core_container
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE class_core_container
!----------------------------------------------------------------------------------------
