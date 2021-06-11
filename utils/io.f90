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
MODULE env_io
    !------------------------------------------------------------------------------------
    !
    USE env_utils_param
    !
    USE env_mp, ONLY: env_mp_bcast, env_mp_abort, env_mp_rank
    !
    USE env_char_ops, ONLY: env_uppercase
    !
#if defined(__PTRACE)&&defined(__INTEL_COMPILER)
    USE ifcore, ONLY: tracebackqq
#endif
    !
    !------------------------------------------------------------------------------------
    !
    PRIVATE
    !
    PUBLIC :: env_find_free_unit, env_read_line, env_field_count, env_get_field
    !
    !------------------------------------------------------------------------------------
CONTAINS
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    INTEGER FUNCTION env_find_free_unit()
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER :: iunit
        LOGICAL :: opnd
        !
        CHARACTER(LEN=80) :: fun_name = 'env_find_free_unit'
        !
        !--------------------------------------------------------------------------------
        !
        env_find_free_unit = -1
        !
        unit_loop: DO iunit = 99, 1, -1
            !
            INQUIRE (UNIT=iunit, OPENED=opnd)
            !
            IF (.NOT. opnd) THEN
                env_find_free_unit = iunit
                !
                RETURN
                !
            END IF
            !
        END DO unit_loop
        !
        CALL env_infomsg(fun_name, 'Free unit not found?!?')
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END FUNCTION env_find_free_unit
    !------------------------------------------------------------------------------------
    !>
    !! WE MAY WANT TO ADD A SECOND COMM ON IMAGES #TODO may be required for NEB
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_read_line(unit, line, nfield, field, end_of_file, error, ionode, &
                             ionode_id, comm)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: unit, ionode_id, comm
        LOGICAL, INTENT(IN) :: ionode
        CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: field
        INTEGER, OPTIONAL, INTENT(IN) :: nfield
        !
        CHARACTER(LEN=*), INTENT(OUT) :: line
        LOGICAL, OPTIONAL, INTENT(OUT) :: end_of_file, error
        !
        LOGICAL :: tend, terr
        !
        CHARACTER(LEN=80) :: sub_name = 'env_read_line'
        !
        !--------------------------------------------------------------------------------
        !
        IF (LEN(line) < 256) &
            CALL env_errore(sub_name, 'Input line too short', MAX(LEN(line), 1))
        !
        tend = .FALSE.
        terr = .FALSE.
        !
        IF (ionode) THEN
30          READ (unit, fmt='(A256)', ERR=15, END=10) line
            line = TRIM(line)
            !
            IF (line == ' ' .OR. (line(1:1) == '#' .OR. &
                                  line(1:1) == '!' .OR. &
                                  line(1:1) == '/')) &
                GOTO 30
            !
            IF (line(1:1) == '&') THEN
                line = env_uppercase(line(2:))
                !
                CALL env_warning('Skipping unnecessary '//TRIM(line)//' namelist')
                !
                DO WHILE (line(1:1) /= '/') ! consume namelist
                    READ (unit, fmt='(A256)', ERR=15, END=10) line
                END DO
                !
                GOTO 30
            END IF
            !
            GOTO 20
10          tend = .TRUE.
            GOTO 20
15          terr = .TRUE.
20          CONTINUE
        END IF
        !
        CALL env_mp_bcast(tend, ionode_id, comm)
        !
        CALL env_mp_bcast(terr, ionode_id, comm)
        !
        CALL env_mp_bcast(line, ionode_id, comm)
        !
        IF (PRESENT(end_of_file)) THEN
            end_of_file = tend
        ELSE IF (tend) THEN
            CALL env_infomsg(sub_name, 'End of file')
        END IF
        !
        IF (PRESENT(error)) THEN
            error = terr
        ELSE IF (terr) THEN
            CALL env_infomsg(sub_name, 'Read error')
        END IF
        !
        IF (PRESENT(field) .AND. .NOT. (tend .OR. terr)) &
            CALL env_field_compare(line, nfield, field)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_read_line
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_field_count(num, line, car)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CHARACTER(LEN=*), INTENT(IN) :: line
        CHARACTER(LEN=1), OPTIONAL, INTENT(IN) :: car
        !
        INTEGER, INTENT(OUT) :: num
        !
        CHARACTER(LEN=1) :: sep1, sep2
        INTEGER :: j
        !
        !--------------------------------------------------------------------------------
        !
        num = 0
        !
        IF (.NOT. PRESENT(car)) THEN
            sep1 = CHAR(32) ! blank character
            sep2 = CHAR(9) ! tab character
            !
            DO j = 2, MAX(LEN(line), 256)
                !
                IF (line(j:j) == '!' .OR. line(j:j) == CHAR(0)) THEN
                    !
                    IF (line(j - 1:j - 1) /= sep1 .AND. line(j - 1:j - 1) /= sep2) &
                        num = num + 1
                    !
                    EXIT
                    !
                END IF
                !
                IF ((line(j:j) == sep1 .OR. line(j:j) == sep2) .AND. &
                    (line(j - 1:j - 1) /= sep1 .AND. line(j - 1:j - 1) /= sep2)) &
                    num = num + 1
                !
            END DO
            !
        ELSE
            !
            sep1 = car
            !
            DO j = 2, MAX(LEN(line), 256)
                !
                IF (line(j:j) == '!' .OR. &
                    line(j:j) == CHAR(0) .OR. line(j:j) == CHAR(32)) THEN
                    !
                    IF (line(j - 1:j - 1) /= sep1) num = num + 1
                    !
                    EXIT
                    !
                END IF
                !
                IF (line(j:j) == sep1 .AND. line(j - 1:j - 1) /= sep1) num = num + 1
                !
            END DO
            !
        END IF
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_field_count
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_field_compare(str, nf, var)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: nf
        CHARACTER(LEN=*), INTENT(IN) :: str, var
        !
        INTEGER :: nc
        !
        CHARACTER(LEN=80) :: sub_name = 'env_field_compare'
        !
        !--------------------------------------------------------------------------------
        !
        CALL env_field_count(nc, str)
        !
        IF (nc < nf) CALL env_errore(sub_name, 'Wrong number of fields: '//TRIM(var), 1)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_field_compare
    !------------------------------------------------------------------------------------
    !>
    !! Extract whitespace-separated nth block from string
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_get_field(n, field, str, sep)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: n
        CHARACTER(LEN=*), INTENT(IN) :: str
        CHARACTER(LEN=1), OPTIONAL, INTENT(IN) :: sep
        !
        CHARACTER(LEN=*), INTENT(OUT) :: field
        !
        INTEGER :: i, j, z ! block start and end
        INTEGER :: k ! block counter
        CHARACTER(LEN=1) :: sep1, sep2
        !
        !--------------------------------------------------------------------------------
        !
        IF (PRESENT(sep)) THEN
            sep1 = sep
            sep2 = sep ! redundant, but easy
        ELSE
            sep1 = CHAR(32) ! blank character
            sep2 = CHAR(9) ! tab char
        END IF
        !
        k = 1 ! counter for the required block
        !
        DO i = 1, LEN(str)
            !
            z = MAX(i - 1, 1) ! look for the beginning of the required block
            !
            IF (k == n) EXIT
            !
            IF ((str(i:i) == sep1 .OR. str(i:i) == sep2) .AND. &
                (str(z:z) /= sep1 .AND. str(z:z) /= sep2)) &
                k = k + 1
            !
        END DO
        !
        DO j = i, LEN(str)
            !
            z = MAX(j - 1, 1) ! look for the beginning of the next block
            !
            IF ((str(j:j) == sep1 .OR. str(j:j) == sep2) .AND. &
                (str(z:z) /= sep1 .AND. str(z:z) /= sep2)) &
                k = k + 1
            !
            IF (k > n) EXIT
            !
        END DO
        !
        IF (j <= LEN(str)) THEN
            field = TRIM(ADJUSTL(str(i:j - 1)))
            ! if we are here, the reqired block was followed by a separator
            ! and another field. We have to trash one char (a separator)
        ELSE
            field = TRIM(ADJUSTL(str(i:LEN(str))))
            ! if we are here, it was the last block in str. We have to take
            ! all the remaining chars
        END IF
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_get_field
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE env_io
!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
!
! Message routines (outside of module for easy global use)
!
!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
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
!----------------------------------------------------------------------------------------
SUBROUTINE env_errore(calling_routine, message, ierr)
    !------------------------------------------------------------------------------------
    !
    USE env_mp, ONLY: env_mp_abort, env_mp_rank
    USE env_io, ONLY: env_find_free_unit
    USE env_utils_param
    !
    !------------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=*), INTENT(IN) :: calling_routine
    ! the name of the calling calling_routine
    !
    CHARACTER(LEN=*), INTENT(IN) :: message ! the output message
    INTEGER, INTENT(IN) :: ierr ! the error flag
    !
    INTEGER :: crashunit, mpime
    CHARACTER(LEN=6) :: cerr
    !
    !------------------------------------------------------------------------------------
    !
    IF (ierr <= 0) RETURN
    !
    !------------------------------------------------------------------------------------
    ! The error message is written in the "*" unit
    !
    WRITE (cerr, FMT='(I6)') ierr
    WRITE (UNIT=*, FMT='(/,1X,78("%"))')
    !
    WRITE (UNIT=*, FMT='(5X,"Error in routine ",A," (",A,"):")') &
        TRIM(calling_routine), TRIM(ADJUSTL(cerr))
    !
    WRITE (UNIT=*, FMT='(5X,A)') TRIM(message)
    WRITE (UNIT=*, FMT='(1X,78("%"),/)')
    !
#if defined(__MPI)&&defined(__AIX)
    !
    !------------------------------------------------------------------------------------
    ! In the case of ibm machines it is also written on the "0" unit
    ! which is automatically connected to stderr
    !
    WRITE (UNIT=0, FMT='(/,1X,78("%"))')
    !
    WRITE (UNIT=0, FMT='(5X,"Error in routine ",A," (",A,"):")') &
        TRIM(calling_routine), TRIM(ADJUSTL(cerr))
    !
    WRITE (UNIT=0, FMT='(5X,A)') TRIM(message)
    WRITE (UNIT=0, FMT='(1X,78("%"),/)')
    !
#endif
    !
    WRITE (*, '("     stopping ...")')
    !
    FLUSH (stdout)
    !
#if defined(__PTRACE)
#if defined(__INTEL_COMPILER)
    CALL tracebackqq(user_exit_code=-1)
#elif __GFORTRAN__
#if (__GNUC__>4)||((__GNUC__==4)&&(__GNUC_MINOR__>=8))
    CALL backtrace
#endif
#else
    WRITE (UNIT=0, FMT='(5X,A)') "Printing strace..."
    !
    CALL ptrace()
#endif
#endif
!
#if defined(__MPI)
    !
    mpime = env_mp_rank(MPI_COMM_WORLD)
    !
    !------------------------------------------------------------------------------------
    ! Write the message to a file and close it before exiting
    ! This will prevent loss of information on systems that
    ! do not flush the open streams
    ! Added by C.C.
    !
    crashunit = env_find_free_unit()
    OPEN (UNIT=crashunit, FILE=crash_file, POSITION='APPEND', STATUS='UNKNOWN')
    !
    WRITE (UNIT=crashunit, FMT='(/,1X,78("%"))')
    WRITE (UNIT=crashunit, FMT='(5X,"task #",I10)') mpime
    !
    WRITE (UNIT=crashunit, &
           FMT='(5X,"from ",A," : error #",I10)') calling_routine, ierr
    !
    WRITE (UNIT=crashunit, FMT='(5X,A)') message
    WRITE (UNIT=crashunit, FMT='(1X,78("%"),/)')
    !
    CLOSE (UNIT=crashunit)
    !
    CALL env_mp_abort(1, MPI_COMM_WORLD) ! try to exit in a smooth way
#endif
    !
    STOP 1
    !
    RETURN
    !
    !------------------------------------------------------------------------------------
END SUBROUTINE env_errore
!----------------------------------------------------------------------------------------
!>
!! This is a simple routine which writes an info message
!! from a given routine to output.
!!
!----------------------------------------------------------------------------------------
SUBROUTINE env_infomsg(routine, message)
    !------------------------------------------------------------------------------------
    !
    USE env_utils_param, ONLY: stdout
    !
    !------------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=*) :: routine ! the name of the calling routine
    CHARACTER(LEN=*) :: message ! the output message
    !
    !------------------------------------------------------------------------------------
    !
    WRITE (stdout, '(5X,"Message from routine ",A,":")') routine
    WRITE (stdout, '(5X,A)') message
    !
    RETURN
    !
    !------------------------------------------------------------------------------------
END SUBROUTINE env_infomsg
!----------------------------------------------------------------------------------------
!>
!! Writes a WARNING message to output
!!
!----------------------------------------------------------------------------------------
SUBROUTINE env_warning(message)
    !------------------------------------------------------------------------------------
    !
    USE env_utils_param, ONLY: stdout
    !
    !------------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=*) :: message ! the output message
    !
    !------------------------------------------------------------------------------------
    !
    WRITE (stdout, '("Warning: ", A)') message
    !
    RETURN
    !
    !------------------------------------------------------------------------------------
END SUBROUTINE env_warning
!----------------------------------------------------------------------------------------
!>
!! Raises an error due to an invalid input option
!!
!----------------------------------------------------------------------------------------
SUBROUTINE env_invalid_opt(routine, param, input)
    !------------------------------------------------------------------------------------
    !
    USE env_utils_param, ONLY: stdout
    !
    !------------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=*), INTENT(IN) :: routine ! the calling routine
    CHARACTER(LEN=*), INTENT(IN) :: param ! the input parameter
    CHARACTER(LEN=*), INTENT(IN) :: input ! the actual input
    !
    !------------------------------------------------------------------------------------
    !
    CALL env_errore(routine, &
                    "'"//TRIM(input)//"' is not a valid option for "//TRIM(param), 1)
    !
    RETURN
    !
    !------------------------------------------------------------------------------------
END SUBROUTINE env_invalid_opt
!----------------------------------------------------------------------------------------
