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
    USE env_util_param
    !
    USE env_mp, ONLY: env_mp_bcast, env_mp_abort, env_mp_rank
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
        CALL env_infomsg(fun_name, 'free unit not found ?!?')
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
            CALL env_errore(sub_name, ' input line too short ', MAX(LEN(line), 1))
        !
        tend = .FALSE.
        terr = .FALSE.
        !
        IF (ionode) THEN
30          READ (unit, fmt='(A256)', ERR=15, END=10) line
            !
            IF (line == ' ' .OR. line(1:1) == '#' .OR. line(1:1) == '!') GOTO 30
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
            CALL env_infomsg(sub_name, ' end of file ')
        END IF
        !
        IF (PRESENT(error)) THEN
            error = terr
        ELSE IF (terr) THEN
            CALL env_infomsg(sub_name, ' read error ')
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
        IF (nc < nf) &
            CALL env_errore(sub_name, ' wrong number of fields: '//TRIM(var), 1)
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
