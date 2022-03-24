!----------------------------------------------------------------------------------------
!
! Copyright (C) 2021 ENVIRON (www.quantum-environ.org)
! Copyright (C) 2001-2008 Quantum ESPRESSO (www.quantum-espresso.org)
!
!----------------------------------------------------------------------------------------
!
!     This file is part of Environ version 2.1
!
!     Environ 2.1 is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation, either version 2 of the License, or
!     (at your option) any later version.
!
!     Environ 2.1 is distributed in the hope that it will be useful,
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
!! Array operations
!!
!----------------------------------------------------------------------------------------
MODULE env_array_ops
    !------------------------------------------------------------------------------------
    !
    USE class_io, ONLY: io
    !
    !------------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    PRIVATE
    !
    PUBLIC :: env_get_index
    !
    !------------------------------------------------------------------------------------
    !
    INTERFACE env_get_index
        MODULE PROCEDURE &
            env_get_index_integer, &
            env_get_index_char
    END INTERFACE env_get_index
    !
    !------------------------------------------------------------------------------------
CONTAINS
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    INTEGER FUNCTION env_get_index_integer(element, array) RESULT(index)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: element, array(:)
        !
        INTEGER :: i
        !
        CHARACTER(LEN=80) :: fun_name = 'env_get_index:integer'
        !
        !--------------------------------------------------------------------------------
        !
        DO i = 1, SIZE(array)
            !
            IF (array(i) == element) THEN
                index = i
                !
                RETURN
                !
            END IF
            !
        END DO
        !
        CALL io%error(fun_name, "Element not found", 1)
        !
        !--------------------------------------------------------------------------------
    END FUNCTION env_get_index_integer
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    INTEGER FUNCTION env_get_index_char(element, array) RESULT(index)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CHARACTER(LEN=*), INTENT(IN) :: element, array(:)
        !
        INTEGER :: i
        !
        CHARACTER(LEN=80) :: fun_name = 'env_get_index:char'
        !
        !--------------------------------------------------------------------------------
        !
        DO i = 1, SIZE(array)
            !
            IF (array(i) == element) THEN
                index = i
                !
                RETURN
                !
            END IF
            !
        END DO
        !
        CALL io%error(fun_name, "Element not found", 1)
        !
        !--------------------------------------------------------------------------------
    END FUNCTION env_get_index_char
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE env_array_ops
!----------------------------------------------------------------------------------------
