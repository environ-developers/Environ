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
MODULE class_core_numerical
    !------------------------------------------------------------------------------------
    !
    USE class_cell
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
    TYPE, ABSTRACT, PUBLIC :: numerical_core
        !--------------------------------------------------------------------------------
        !
        CHARACTER(LEN=80) :: core_type
        !
        TYPE(environ_cell), POINTER :: cell => NULL()
        !
        !--------------------------------------------------------------------------------
    CONTAINS
        !--------------------------------------------------------------------------------
        !
        PROCEDURE(create_core), DEFERRED :: create
        PROCEDURE(destroy_core), DEFERRED :: destroy
        !
        !--------------------------------------------------------------------------------
    END TYPE numerical_core
    !------------------------------------------------------------------------------------
    !
    ABSTRACT INTERFACE
        SUBROUTINE create_core(this)
            IMPORT numerical_core
            CLASS(numerical_core), INTENT(INOUT) :: this
        END SUBROUTINE
        SUBROUTINE destroy_core(this)
            IMPORT numerical_core
            CLASS(numerical_core), INTENT(INOUT) :: this
        END SUBROUTINE
    END INTERFACE
    !
    !------------------------------------------------------------------------------------
END MODULE class_core_numerical
!----------------------------------------------------------------------------------------
