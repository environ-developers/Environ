!----------------------------------------------------------------------------------------
!
! Copyright (C) 2021 ENVIRON (www.quantum-environ.org)
! Copyright (C) 2002-2009 Quantum ESPRESSO (www.quantum-espresso.org)
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
! Authors: Compiled and modified by Edan Bainglass
!
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
    USE env_base_io, ONLY: io
    USE env_mp, ONLY: env_mp_abort, env_mp_rank
    USE env_parallel_include, ONLY: MPI_COMM_WORLD
    !
#if defined(__PTRACE)&&defined(__INTEL_COMPILER)
    USE ifcore, ONLY: tracebackqq
#endif
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
    ! #TODO figure this out
    !
    ! IF (.NOT. ionode) CALL env_mp_abort(1, MPI_COMM_WORLD) ! #TODO what if ionode != 0?
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
    FLUSH (io%unit) ! #TODO why?
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
    crashunit = io%find_free_unit()
    OPEN (UNIT=crashunit, FILE='CRASH', POSITION='APPEND', STATUS='UNKNOWN')
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
    !------------------------------------------------------------------------------------
END SUBROUTINE env_errore
!----------------------------------------------------------------------------------------
!>
!! Writes an indented (5X) info message to output
!!
!----------------------------------------------------------------------------------------
SUBROUTINE env_write(message)
    !------------------------------------------------------------------------------------
    !
    USE env_base_io, ONLY: io
    !
    !------------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=*), INTENT(IN) :: message
    !
    !------------------------------------------------------------------------------------
    !
    WRITE (io%unit, '(5X,A)') TRIM(message)
    !
    !------------------------------------------------------------------------------------
END SUBROUTINE env_write
!----------------------------------------------------------------------------------------
!>
!! Writes an indented (5X) '=' divider (blank line above) to output
!! Optional blank line below
!!
!----------------------------------------------------------------------------------------
SUBROUTINE env_divider(lblank)
    !------------------------------------------------------------------------------------
    !
    USE env_base_io, ONLY: io
    !
    !------------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    LOGICAL, INTENT(IN) :: lblank
    !
    !------------------------------------------------------------------------------------
    !
    WRITE (io%unit, FMT=1)
    !
    IF (lblank) WRITE (io%unit, *) ! blank line
    !
1   FORMAT(/, 5X, 80('='))
    !
    !------------------------------------------------------------------------------------
END SUBROUTINE env_divider
!----------------------------------------------------------------------------------------
!>
!! Writes an indented (5X) header (blank line above) to output
!!
!----------------------------------------------------------------------------------------
SUBROUTINE env_header(message)
    !------------------------------------------------------------------------------------
    !
    USE env_base_io, ONLY: io
    !
    !------------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=*), INTENT(IN) :: message
    !
    !------------------------------------------------------------------------------------
    !
    WRITE (io%unit, '(/,5X,A)') TRIM(message)
    !
    !------------------------------------------------------------------------------------
END SUBROUTINE env_header
!----------------------------------------------------------------------------------------
!>
!! Writes a message message to output warning the user of a non-terminating issue
!!
!----------------------------------------------------------------------------------------
SUBROUTINE env_warning(message)
    !------------------------------------------------------------------------------------
    !
    USE env_base_io, ONLY: io
    !
    !------------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=*), INTENT(IN) :: message ! the output message
    !
    !------------------------------------------------------------------------------------
    !
    WRITE (io%unit, '("Warning: ", A,/)') TRIM(message)
    !
    !------------------------------------------------------------------------------------
END SUBROUTINE env_warning
!----------------------------------------------------------------------------------------
!>
!! Writes a commented message to output regarding a forced default setting
!!
!----------------------------------------------------------------------------------------
SUBROUTINE env_default(param, default, message)
    !------------------------------------------------------------------------------------
    !
    USE env_base_io, ONLY: io
    !
    !------------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=*), INTENT(IN) :: param ! input parameter
    CHARACTER(LEN=*), INTENT(IN) :: default ! default value
    CHARACTER(LEN=*), INTENT(IN) :: message ! additional comment
    !
    CHARACTER(LEN=80) :: comment
    !
    CHARACTER(LEN=25) :: p
    CHARACTER(LEN=15) :: d
    !
    !------------------------------------------------------------------------------------
    !
    IF (message /= '') THEN
        comment = '! '//message
    ELSE
        comment = ''
    END IF
    !
    p = ADJUSTL(param)
    d = ADJUSTL(default)
    !
    WRITE (io%unit, '(5X,A," = ",A, A)') p, d, TRIM(ADJUSTL(comment))
    !
    !------------------------------------------------------------------------------------
END SUBROUTINE env_default
!----------------------------------------------------------------------------------------
!>
!! Raises an error due to an invalid input option
!!
!----------------------------------------------------------------------------------------
SUBROUTINE env_invalid_opt(routine, param, input)
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
    !------------------------------------------------------------------------------------
END SUBROUTINE env_invalid_opt
!----------------------------------------------------------------------------------------
!>
!! Prints a block divider of variable length determnied by verbosity
!!
!! @param verbose      : (INTEGER) local verbose of the calling printer
!! @param base_verbose : (INTEGER) verbosity to be measured against
!! @param unit         : (INTEGER) output target
!!
!----------------------------------------------------------------------------------------
SUBROUTINE env_block_divider(verbose, base_verbose, unit)
    !------------------------------------------------------------------------------------
    !
    USE env_base_io, ONLY: io
    !
    !------------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: verbose, base_verbose, unit
    !
    !------------------------------------------------------------------------------------
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
    END SELECT
    !
    !------------------------------------------------------------------------------------
    !
10  FORMAT(1X, 79('_'))
11  FORMAT(1X, 75('_'))
12  FORMAT(1X, 71('_'))
13  FORMAT(1X, 67('_'))
14  FORMAT(1X, 63('_'))
    !
    !------------------------------------------------------------------------------------
END SUBROUTINE env_block_divider
!----------------------------------------------------------------------------------------
!>
!!
!----------------------------------------------------------------------------------------
SUBROUTINE env_create_error(routine)
    !------------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=80), INTENT(IN) :: routine
    !
    !------------------------------------------------------------------------------------
    !
    CALL env_errore(routine, 'Trying to create an existing object', 1)
    !
    !------------------------------------------------------------------------------------
END SUBROUTINE env_create_error
!----------------------------------------------------------------------------------------
!>
!!
!----------------------------------------------------------------------------------------
SUBROUTINE env_destroy_error(routine)
    !------------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=80), INTENT(IN) :: routine
    !
    !------------------------------------------------------------------------------------
    !
    CALL env_errore(routine, 'Trying to destroy an empty object', 1)
    !
    !------------------------------------------------------------------------------------
END SUBROUTINE env_destroy_error
!----------------------------------------------------------------------------------------
