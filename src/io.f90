!----------------------------------------------------------------------------------------
!
! Copyright (C) 2021 ENVIRON (www.quantum-environment.org)
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
! Authors: Compiled by Edan Bainglass (Department of Physics, UNT)
!
!----------------------------------------------------------------------------------------
!>
!!
!----------------------------------------------------------------------------------------
MODULE class_io
    !------------------------------------------------------------------------------------
    !
    USE env_char_ops, ONLY: env_uppercase
    !
    IMPLICIT NONE
    !
    PRIVATE
    !
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    TYPE environ_io
        !--------------------------------------------------------------------------------
        !
        LOGICAL :: lnode = .TRUE.
        INTEGER :: node = 0
        !
        LOGICAL :: lstdout = .FALSE. ! whether environ can print on standard output
        !
        INTEGER :: comm ! WE MAY NEED A SECOND COMMUNICATOR FOR IMAGE PARALLELIZATION
        !
        INTEGER :: unit
        INTEGER :: debug_unit
        !
        INTEGER :: verbosity = 0
        !
        !--------------------------------------------------------------------------------
    CONTAINS
        !--------------------------------------------------------------------------------
        !
        PROCEDURE, PRIVATE :: create => create_base_io
        PROCEDURE :: init => init_base_io
        PROCEDURE :: update_unit => update_output_program_unit
        !
        PROCEDURE, NOPASS :: find_free_unit
        !
        PROCEDURE, NOPASS :: error
        PROCEDURE, NOPASS :: create_error => env_create_error
        PROCEDURE, NOPASS :: destroy_error => env_destroy_error
        PROCEDURE, NOPASS :: invalid_opt => env_invalid_opt
        !
        PROCEDURE, NOPASS :: writer => env_write
        PROCEDURE, NOPASS :: divider => env_divider
        PROCEDURE, NOPASS :: header => env_header
        PROCEDURE, NOPASS :: warning => env_warning
        PROCEDURE, NOPASS :: default => env_default
        PROCEDURE, NOPASS :: block_divider => env_block_divider
        !
        !--------------------------------------------------------------------------------
    END TYPE environ_io
    !------------------------------------------------------------------------------------
    !
    TYPE(environ_io), PUBLIC, TARGET, SAVE :: io
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
    SUBROUTINE create_base_io(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_io), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: sub_name = 'create_base_io'
        !
        !--------------------------------------------------------------------------------
        !
        this%lnode = .TRUE.
        this%node = 0
        this%lstdout = .FALSE.
        this%comm = 0
        this%unit = 0
        this%debug_unit = 0
        this%verbosity = 0
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE create_base_io
    !------------------------------------------------------------------------------------
    !>
    !! Set global I/O constants
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_base_io(this, can_write, id, comm, prog_unit, lstdout)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_io), INTENT(INOUT) :: this
        !
        LOGICAL, INTENT(IN) :: can_write
        INTEGER, INTENT(IN) :: id
        INTEGER, INTENT(IN) :: comm
        INTEGER, INTENT(IN) :: prog_unit
        LOGICAL, INTENT(IN) :: lstdout
        !
        !--------------------------------------------------------------------------------
        !
        CALL this%create()
        !
        this%lnode = can_write
        this%node = id
        this%comm = comm
        this%unit = prog_unit
        this%lstdout = lstdout
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_base_io
    !------------------------------------------------------------------------------------
    !>
    !! Sets the output file target
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE update_output_program_unit(this, program_unit_in)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: program_unit_in
        !
        CLASS(environ_io), INTENT(INOUT) :: this
        !
        !--------------------------------------------------------------------------------
        !
        this%unit = program_unit_in
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE update_output_program_unit
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !                                  GENERAL METHODS
    !
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    INTEGER FUNCTION find_free_unit()
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, EXTERNAL :: env_find_free_unit
        !
        !--------------------------------------------------------------------------------
        !
        find_free_unit = env_find_free_unit()
        !
        !--------------------------------------------------------------------------------
    END FUNCTION find_free_unit
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !                                   ERROR METHODS
    !
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !>
    !! This is a simple routine which writes an error message to output:
    !! if ierr <= 0 it does nothing,
    !! if ierr  > 0 it stops.
    !!
    !! **** Important note for parallel execution ***
    !!
    !! In parallel execution unit 6 is written only by the first node;
    !! all other nodes have unit 6 redirected to nothing (/dev/null).
    !! As a consequence an error not occurring on the first node
    !! will be invisible. For T3E and ORIGIN machines, this problem
    !! is solved by writing an error message to unit * instead of 6.
    !! Whenever possible (IBM SP machines), we write to the standard
    !! error, unit 0 (the message will appear in the error files
    !! produced by loadleveler).
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE error(calling_routine, message, ierr)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CHARACTER(LEN=*), INTENT(IN) :: calling_routine
        CHARACTER(LEN=*), INTENT(IN) :: message
        INTEGER, INTENT(IN) :: ierr
        !
        !--------------------------------------------------------------------------------
        !
        CALL env_errore(calling_routine, message, ierr)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE error
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_create_error(routine)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CHARACTER(LEN=80), INTENT(IN) :: routine
        !
        !--------------------------------------------------------------------------------
        !
        IF (.NOT. io%lnode) RETURN
        !
        CALL io%error(routine, "Trying to create an existing object", 1)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_create_error
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_destroy_error(routine)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CHARACTER(LEN=80), INTENT(IN) :: routine
        !
        !--------------------------------------------------------------------------------
        !
        IF (.NOT. io%lnode) RETURN
        !
        CALL io%error(routine, "Trying to destroy an empty object", 1)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_destroy_error
    !------------------------------------------------------------------------------------
    !>
    !! Raises an error due to an invalid input option
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_invalid_opt(routine, param, input)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CHARACTER(LEN=*), INTENT(IN) :: routine ! the calling routine
        CHARACTER(LEN=*), INTENT(IN) :: param ! the input parameter
        CHARACTER(LEN=*), INTENT(IN) :: input ! the actual input
        !
        !--------------------------------------------------------------------------------
        !
        IF (.NOT. io%lnode) RETURN
        !
        CALL io%error(routine, "'"//TRIM(input)//"' is not a valid option for "//param, 1)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_invalid_opt
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !                                   OUTPUT METHODS
    !
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !>
    !! Writes an indented (5X) info message to output
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_write(message)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CHARACTER(LEN=*), INTENT(IN) :: message
        !
        !--------------------------------------------------------------------------------
        !
        IF (.NOT. io%lnode) RETURN
        !
        WRITE (io%unit, '(5X,A)') message
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_write
    !------------------------------------------------------------------------------------
    !>
    !! Writes an indented (5X) '=' divider (blank line above) to output
    !! Optional blank line below
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_divider(lblank)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        LOGICAL, OPTIONAL, INTENT(IN) :: lblank
        !
        LOGICAL :: local_blank = .FALSE.
        !
        !--------------------------------------------------------------------------------
        !
        IF (.NOT. io%lnode) RETURN
        !
        IF (PRESENT(lblank)) local_blank = lblank
        !
        WRITE (io%unit, FMT=1)
        !
        IF (local_blank) WRITE (io%unit, *) ! blank line
        !
1       FORMAT(/, 5X, 80('='))
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_divider
    !------------------------------------------------------------------------------------
    !>
    !! Writes an indented (5X) header (blank line above) to output
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_header(message)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CHARACTER(LEN=*), INTENT(IN) :: message
        !
        !--------------------------------------------------------------------------------
        !
        IF (.NOT. io%lnode) RETURN
        !
        WRITE (io%unit, *)
        !
        CALL io%writer(message)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_header
    !------------------------------------------------------------------------------------
    !>
    !! Writes a message message to output warning the user of a non-terminating issue
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_warning(message, ierr)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CHARACTER(LEN=*), INTENT(IN) :: message
        INTEGER, INTENT(IN) :: ierr
        !
        CHARACTER(LEN=6) :: cerr
        !
        !--------------------------------------------------------------------------------
        !
        IF (.NOT. io%lnode) RETURN
        !
        WRITE (cerr, '(I6)') ierr
        !
        WRITE (io%unit, 20) TRIM(ADJUSTL(cerr)), message
        !
20      FORMAT(/, "Environ Warning (", A, "): ", A,/)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_warning
    !------------------------------------------------------------------------------------
    !>
    !! Writes a commented message to output regarding a forced default setting
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_default(param, default, message)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CHARACTER(LEN=*), INTENT(IN) :: param ! input parameter
        CHARACTER(LEN=*), INTENT(IN) :: default ! default value
        CHARACTER(LEN=*), INTENT(IN) :: message ! additional comment
        !
        CHARACTER(LEN=80) :: comment
        !
        !--------------------------------------------------------------------------------
        !
        IF (.NOT. io%lnode) RETURN
        !
        IF (message /= '') THEN
            comment = '! '//message
        ELSE
            comment = ''
        END IF
        !
        WRITE (io%unit, '(5X,A," = ",A, A)') param, default, comment
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_default
    !------------------------------------------------------------------------------------
    !>
    !! Prints a block divider of variable length determnied by verbosity
    !!
    !! @param verbose      : (INTEGER) local verbose of the calling printer
    !! @param base_verbose : (INTEGER) verbosity to be measured against
    !! @param unit         : (INTEGER) output target
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_block_divider(verbose, base_verbose, unit)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: verbose, base_verbose, unit
        !
        CHARACTER(LEN=80) :: sub_name = 'env_block_divider'
        !
        !--------------------------------------------------------------------------------
        !
        IF (.NOT. io%lnode) RETURN
        !
        SELECT CASE (base_verbose - verbose - 1)
            !
        CASE (0)
            WRITE (unit, 10)
            !
        CASE (1)
            WRITE (unit, 11)
            !
        CASE (2)
            WRITE (unit, 12)
            !
        CASE (3)
            WRITE (unit, 13)
            !
        CASE (4)
            WRITE (unit, 14)
            !
        CASE DEFAULT
            CALL io%error(sub_name, "Unexpected verbose value", 1)
            !
        END SELECT
        !
        !--------------------------------------------------------------------------------
        !
10      FORMAT(1X, 79('_'))
11      FORMAT(1X, 75('_'))
12      FORMAT(1X, 71('_'))
13      FORMAT(1X, 67('_'))
14      FORMAT(1X, 63('_'))
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_block_divider
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE class_io
!----------------------------------------------------------------------------------------
