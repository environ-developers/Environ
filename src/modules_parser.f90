!
! Copyright (C) 2001-2004 Carlo Cavazzoni and PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------------------
!>
!! SUBROUTINE field_count:   accepts two string (one of them is optional)
!!                           and one integer and count the number of fields
!!                           in the string separated by a blank or a tab
!!                           character. If the optional string is specified
!!                           (it has anyway len=1) it is assumed as the
!!                           separator character.
!!                           Ignores any character following the exclamation
!!                           mark (fortran comment)
!!
!! SUBROUTINE con_cam:       counts the number of fields in a string
!!                           separated by the optional character
!!
!! SUBROUTINE field_compare: accepts two strings and one integer. Counts the
!!                           fields contained in the first string and
!!                           compares it with the integer.
!!                           If they are less than the integer calls the
!!                           routine error and show by the second string the
!!                           name of the field where read-error occurred.
!!
!! SUBROUTINE version_parse: Determine the major, minor and patch numbers
!!                           from a version string with the fmt "i.j.k"
!!
!! FUNCTION version_compare: Compare two version strings; the result can be
!!                           "newer", "equal", "older", ""
!!
!! #TODO: several unused subroutines/functions
!!
!----------------------------------------------------------------------------------------
MODULE modules_parser
    !------------------------------------------------------------------------------------
    !
    USE environ_output, ONLY: program_unit
    USE modules_constants, ONLY: DP
    !
    PRIVATE
    !
    PUBLIC :: parse_unit, env_field_count, env_read_line, env_get_field
    PUBLIC :: env_version_parse, env_version_compare
    !
    INTEGER :: parse_unit = 5 ! normally 5, but can be set otherwise
    !
    !------------------------------------------------------------------------------------
CONTAINS
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
    SUBROUTINE env_read_line(line, nfield, field, end_of_file, error)
        !--------------------------------------------------------------------------------
        !
        USE mp, ONLY: mp_bcast
        USE environ_output, ONLY: comm ! WE MAY WANT TO ADD A SECOND COMM ON IMAGES #TODO may be required for NEB
        USE environ_output, ONLY: ionode, ionode_id
        !
        IMPLICIT NONE
        !
        CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: field
        INTEGER, OPTIONAL, INTENT(IN) :: nfield
        !
        CHARACTER(LEN=*), INTENT(OUT) :: line
        LOGICAL, OPTIONAL, INTENT(OUT) :: end_of_file, error
        !
        LOGICAL :: tend, terr
        !
        !--------------------------------------------------------------------------------
        !
        IF (LEN(line) < 256) &
            CALL errore(' read_line ', ' input line too short ', MAX(LEN(line), 1))
        !
        tend = .FALSE.
        terr = .FALSE.
        !
        IF (ionode) THEN
30          READ (parse_unit, fmt='(A256)', ERR=15, END=10) line
            !
            IF (line == ' ' .OR. line(1:1) == '#') GOTO 30
            !
            GOTO 20
10          tend = .TRUE.
            GOTO 20
15          terr = .TRUE.
20          CONTINUE
        END IF
        !
        CALL mp_bcast(tend, ionode_id, comm)
        !
        CALL mp_bcast(terr, ionode_id, comm)
        !
        CALL mp_bcast(line, ionode_id, comm)
        !
        IF (PRESENT(end_of_file)) THEN
            end_of_file = tend
        ELSE IF (tend) THEN
            CALL infomsg(' read_line ', ' end of file ')
        END IF
        !
        IF (PRESENT(error)) THEN
            error = terr
        ELSE IF (terr) THEN
            CALL infomsg(' read_line ', ' read error ')
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
    SUBROUTINE env_field_compare(str, nf, var)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CHARACTER(LEN=*), INTENT(IN) :: str
        INTEGER, INTENT(IN) :: nf
        CHARACTER(LEN=*), INTENT(IN) :: var
        !
        INTEGER :: nc
        !
        !--------------------------------------------------------------------------------
        !
        CALL env_field_count(nc, str)
        !
        IF (nc < nf) &
            CALL errore(' field_compare ', ' wrong number of fields: '//TRIM(var), 1)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_field_compare
    !------------------------------------------------------------------------------------
    !>
    !! #TODO: unused
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_con_cam(num, line, car)
        !--------------------------------------------------------------------------------
        !
        CHARACTER(LEN=*) :: line
        CHARACTER(LEN=1) :: sep
        CHARACTER(LEN=1), OPTIONAL :: car
        INTEGER :: num, j
        !
        !--------------------------------------------------------------------------------
        !
        num = 0
        !
        IF (LEN(line) > 256) THEN
            WRITE (program_unit, *) 'riga ', line
            WRITE (program_unit, *) 'lunga ', LEN(line)
            num = -1
            !
            RETURN
            !
        END IF
        !
        WRITE (program_unit, *) '1riga ', line
        WRITE (program_unit, *) '1lunga ', LEN(line)
        !
        IF (.NOT. PRESENT(car)) THEN
            sep = CHAR(32) ! char(32) is the blank character
        ELSE
            sep = car
        END IF
        !
        DO j = 2, MAX(LEN(line), 256)
            !
            IF (line(j:j) == '!' .OR. line(j:j) == CHAR(0)) RETURN
            !
            IF ((line(j:j) == sep) .AND. (line(j - 1:j - 1) /= sep)) num = num + 1
            !
        END DO
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_con_cam
    !------------------------------------------------------------------------------------
    !>
    !! Determine the major, minor and patch numbers from
    !! a version string with the fmt "i.j.k"
    !!
    !! The ierr variable assumes the following values
    !!
    !! ierr < 0     emtpy string
    !! ierr = 0     no problem
    !! ierr > 0     fatal error
    !!
    !! #TODO: unused
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_version_parse(str, major, minor, patch, ierr)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CHARACTER(*), INTENT(IN) :: str
        !
        INTEGER, INTENT(OUT) :: major, minor, patch, ierr
        !
        INTEGER :: i1, i2, length
        INTEGER :: ierrtot
        CHARACTER(10) :: num(3)
        !
        !--------------------------------------------------------------------------------
        !
        major = 0
        minor = 0
        patch = 0
        !
        length = LEN_TRIM(str)
        !
        IF (length == 0) THEN
            ierr = -1
            !
            RETURN
            !
        END IF
        !
        i1 = SCAN(str, ".")
        i2 = SCAN(str, ".", BACK=.TRUE.)
        !
        IF (i1 == 0 .OR. i2 == 0 .OR. i1 == i2) THEN
            ierr = 1
            !
            RETURN
            !
        END IF
        !
        num(1) = str(1:i1 - 1)
        num(2) = str(i1 + 1:i2 - 1)
        num(3) = str(i2 + 1:)
        !
        ierrtot = 0
        !
        READ (num(1), *, IOSTAT=ierr) major
        !
        IF (ierr /= 0) RETURN
        !
        READ (num(2), *, IOSTAT=ierr) minor
        !
        IF (ierr /= 0) RETURN
        !
        READ (num(3), *, IOSTAT=ierr) patch
        !
        IF (ierr /= 0) RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_version_parse
    !------------------------------------------------------------------------------------
    !>
    !! Compare two version strings; the result is
    !!
    !! "newer":   str1 is newer that str2
    !! "equal":   str1 is equal   to str2
    !! "older":   str1 is older than str2
    !! " ":       str1 or str2 has a wrong format
    !!
    !! #TODO: unused
    !!
    !------------------------------------------------------------------------------------
    FUNCTION env_version_compare(str1, str2)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        CHARACTER(*) :: str1, str2
        CHARACTER(10) :: env_version_compare
        !
        INTEGER :: version1(3), version2(3)
        INTEGER :: basis, icheck1, icheck2
        INTEGER :: ierr
        !
        !--------------------------------------------------------------------------------
        !
        env_version_compare = " "
        !
        CALL env_version_parse(str1, version1(1), version1(2), version1(3), ierr)
        !
        IF (ierr /= 0) RETURN
        !
        CALL env_version_parse(str2, version2(1), version2(2), version2(3), ierr)
        !
        IF (ierr /= 0) RETURN
        !
        basis = 1000
        !
        icheck1 = version1(1) * basis**2 + version1(2) * basis + version1(3)
        icheck2 = version2(1) * basis**2 + version2(2) * basis + version2(3)
        !
        IF (icheck1 > icheck2) THEN
            env_version_compare = 'newer'
        ELSE IF (icheck1 == icheck2) THEN
            env_version_compare = 'equal'
        ELSE
            env_version_compare = 'older'
        END IF
        !
        !--------------------------------------------------------------------------------
    END FUNCTION env_version_compare
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
        CHARACTER(len=*), INTENT(IN) :: str
        CHARACTER(len=1), OPTIONAL, INTENT(IN) :: sep
        !
        CHARACTER(len=*), INTENT(OUT) :: field
        !
        INTEGER :: i, j, z ! block start and end
        INTEGER :: k ! block counter
        CHARACTER(len=1) :: sep1, sep2
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
END MODULE modules_parser
!----------------------------------------------------------------------------------------