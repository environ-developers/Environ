!----------------------------------------------------------------------------------------
!
! Copyright (C) 2018 ENVIRON (www.quantum-environment.org)
! Copyright (C) 2011 Quantum ESPRESSO group
!
!----------------------------------------------------------------------------------------
!
! This file is part of Environ version 2.0
!
! Environ 2.0 is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 2 of the License, or
! (at your option) any later version.
!
! Environ 2.0 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more detail, either the file
! `License' in the root directory of the present distribution, or
! online at <http://www.gnu.org/licenses/>.
!
!----------------------------------------------------------------------------------------
!
! Authors:
!
!----------------------------------------------------------------------------------------
!>
!!
!----------------------------------------------------------------------------------------
MODULE env_char_ops
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
CONTAINS
    !------------------------------------------------------------------------------------
    !>
    !! Converts character to uppercase
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_uppercase(string, upper_str)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CHARACTER(LEN=*), INTENT(IN) :: string
        !
        CHARACTER(LEN=*), INTENT(OUT) :: upper_str
        !
        CHARACTER(LEN=1) :: c
        INTEGER :: i, ci
        !
        !--------------------------------------------------------------------------------
        !
        upper_str = TRIM(string)
        !
        DO i = 1, LEN(TRIM(string))
            c = string(i:i)
            ci = ICHAR(c)
            !
            IF (ci >= 97 .AND. ci <= 122) THEN
                upper_str(i:i) = CHAR(ci - 32)
            END IF
            !
        END DO
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_uppercase
    !------------------------------------------------------------------------------------
    !>
    !! Converts character to lowercase
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_lowercase(string, lower_str)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CHARACTER(LEN=*), INTENT(IN) :: string
        !
        CHARACTER(LEN=*), INTENT(OUT) :: lower_str
        !
        CHARACTER(LEN=1) :: c
        INTEGER :: i, ci
        !
        !--------------------------------------------------------------------------------
        !
        lower_str = TRIM(string)
        !
        DO i = 1, LEN(TRIM(string))
            c = string(i:i)
            ci = ICHAR(c)
            !
            IF (ci >= 65 .AND. ci <= 90) THEN
                lower_str(i:i) = CHAR(ci + 32)
            END IF
            !
        END DO
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_lowercase
    !------------------------------------------------------------------------------------
    !>
    !!
    !!
    !------------------------------------------------------------------------------------
    LOGICAL FUNCTION env_is_substring(string1, string2)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CHARACTER(LEN=*), INTENT(IN) :: string1, string2
        !
        INTEGER :: len1, len2, l
        !
        !--------------------------------------------------------------------------------
        !
        len1 = LEN_TRIM(string1)
        len2 = LEN_TRIM(string2)
        !
        DO l = 1, (len2 - len1 + 1)
            !
            IF (string1(1:len1) == string2(l:(l + len1 - 1))) THEN
                !
                env_is_substring = .TRUE.
                !
                RETURN
                !
            END IF
            !
        END DO
        !
        env_is_substring = .FALSE.
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END FUNCTION env_is_substring
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE env_char_ops
!----------------------------------------------------------------------------------------
