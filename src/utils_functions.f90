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
!
!----------------------------------------------------------------------------------------
!>
!! Module containing the main routines to handle environ_functions
!! derived data types.
!!
!! Environ_functions contains all the details of analytic functions needed
!! by Environ modules and defined on the three-dimensional real-space
!! domain, together with the routines to handle the derived data type and
!! to generate the functions from their parameters.
!!
!----------------------------------------------------------------------------------------
MODULE utils_functions
    !------------------------------------------------------------------------------------
    !
    USE modules_constants, ONLY: DP
    !
    USE representation_types, ONLY: environ_functions
    !
    !------------------------------------------------------------------------------------
    !
    INTERFACE create_environ_functions
        MODULE PROCEDURE create_environ_functions_scalar, create_environ_functions_array
    END INTERFACE create_environ_functions
    !
    !------------------------------------------------------------------------------------
    !
    PRIVATE
    !
    PUBLIC :: create_environ_functions, copy_environ_functions, &
              destroy_environ_functions
    !
    !------------------------------------------------------------------------------------
CONTAINS
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE create_environ_functions_scalar(type_, dimm, axis, pos, width, &
                                               spreadd, volume, f)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: type_, dimm, axis
        REAL(DP), INTENT(IN) :: width, spreadd, volume
        REAL(DP), TARGET, INTENT(IN) :: pos(3)
        !
        TYPE(environ_functions), INTENT(INOUT) :: f
        !
        INTEGER :: i
        !
        !--------------------------------------------------------------------------------
        !
        f%type_ = type_
        f%dim = dimm
        f%axis = axis
        f%spread = spreadd
        f%width = width
        f%volume = volume
        !
        f%pos => pos
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE create_environ_functions_scalar
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE create_environ_functions_array(n, type_, dimm, axis, pos, width, &
                                              spreadd, volume, f)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: n, type_
        INTEGER, DIMENSION(n), INTENT(IN) :: dimm, axis
        REAL(DP), DIMENSION(n), INTENT(IN) :: width, spreadd, volume
        REAL(DP), TARGET, INTENT(IN) :: pos(3, n)
        !
        TYPE(environ_functions), ALLOCATABLE, INTENT(INOUT) :: f(:)
        !
        INTEGER :: i
        !
        !--------------------------------------------------------------------------------
        !
        ALLOCATE (f(n))
        !
        DO i = 1, n
            !
            CALL create_environ_functions_scalar(type_, dimm(i), axis(i), pos(:, i), &
                                                 width(i), spreadd(i), volume(i), f(i))
            !
        END DO
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE create_environ_functions_array
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE copy_environ_functions(foriginal, fcopy)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_functions), INTENT(IN) :: foriginal
        !
        TYPE(environ_functions), INTENT(OUT) :: fcopy
        !
        !--------------------------------------------------------------------------------
        !
        fcopy%pos => foriginal%pos
        !
        fcopy%type_ = foriginal%type_
        fcopy%dim = foriginal%dim
        fcopy%axis = foriginal%axis
        fcopy%spread = foriginal%spread
        fcopy%width = foriginal%width
        fcopy%volume = foriginal%volume
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE copy_environ_functions
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE destroy_environ_functions(n, f)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: n
        TYPE(environ_functions), ALLOCATABLE, INTENT(INOUT) :: f(:)
        !
        INTEGER :: i
        !
        CHARACTER(LEN=80) :: sub_name = 'destroy_environ_functions'
        !
        !--------------------------------------------------------------------------------
        !
        IF (.NOT. ALLOCATED(f)) &
            CALL env_errore(sub_name, 'Trying to destroy a non allocated object', 1)
        !
        IF (SIZE(f) /= n) &
            CALL env_errore(sub_name, 'Inconsistent size of allocated object', 1)
        !
        DO i = 1, n
            IF (ASSOCIATED(f(i)%pos)) NULLIFY (f(i)%pos)
        END DO
        !
        DEALLOCATE (f)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE destroy_environ_functions
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE utils_functions
!----------------------------------------------------------------------------------------
