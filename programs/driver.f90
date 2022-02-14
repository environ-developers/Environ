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
PROGRAM driver
    !------------------------------------------------------------------------------------
    !
    USE cmdline_args
    !
    USE programs
    !
    !------------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=80) :: prog = ''
    !
    !------------------------------------------------------------------------------------
    !
    CALL parse_command_line()
    !
    CALL run_program()
    !
    CALL clean_up()
    !
    !------------------------------------------------------------------------------------
CONTAINS
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE parse_command_line()
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER :: i
        !
        CHARACTER(LEN=80) :: arg
        !
        !--------------------------------------------------------------------------------
        !
        DO i = 1, COMMAND_ARGUMENT_COUNT(), 2
            !
            CALL GET_COMMAND_ARGUMENT(i, arg)
            !
            SELECT CASE (arg)
                !
            CASE ('-n', '-name')
                CALL GET_COMMAND_ARGUMENT(i + 1, prog)
                !
            CASE ('-i', '-input')
                CALL GET_COMMAND_ARGUMENT(i + 1, inputfile)
                !
            CASE ('-c', '-cube')
                CALL GET_COMMAND_ARGUMENT(i + 1, cubefile)
                !
            CASE ('-with_pbc')
                use_pbc_corr = .TRUE.
                !
            END SELECT
            !
        END DO
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE parse_command_line
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE run_program()
        !--------------------------------------------------------------------------------
        !
        CALL general_setup()
        !
        SELECT CASE (prog)
            !
        CASE ('tester')
            CALL run_tester()
            !
        CASE ('from_cube')
            CALL run_environ_from_cube()
            !
        CASE DEFAULT
            !
            IF (prog /= '') THEN
                !
                PRINT '(/, 5X, A)', TRIM(prog)//" is not available"
                !
                CALL print_available_programs()
                !
            ELSE
                PRINT '(/, 5X, A, /)', "Missing calculation name"
            END IF
            !
        END SELECT
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE run_program
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END PROGRAM driver
!----------------------------------------------------------------------------------------
