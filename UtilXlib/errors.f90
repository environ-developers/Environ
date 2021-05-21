!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
!
! Error routines (outside of module for easy global use)
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
    USE env_util_param
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
#if (__GNUC__>4) || ((__GNUC__==4) && (__GNUC_MINOR__>=8))
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
    USE env_util_param, ONLY: stdout
    !
    !------------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=*) :: routine ! the name of the calling routine
    CHARACTER(LEN=*) :: message ! the output message
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
!! This is a simple routine which writes an error message to output
!!
!----------------------------------------------------------------------------------------
SUBROUTINE env_fft_error(calling_routine, message, ierr)
    !------------------------------------------------------------------------------------
    !
    USE env_util_param, ONLY: MPI_COMM_WORLD
    !
    !------------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=*), INTENT(IN) :: calling_routine
    ! the name of the calling calling_routine

    CHARACTER(LEN=*), INTENT(IN) :: message ! the output message
    INTEGER, INTENT(IN) :: ierr
    !
    CHARACTER(LEN=6) :: cerr
    INTEGER :: info
    !
    !------------------------------------------------------------------------------------
    !
    IF (ierr <= 0) RETURN
    !
    !------------------------------------------------------------------------------------
    ! The error message is written un the "*" unit
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
    WRITE (*, '("     stopping ...")')
    !
#if defined(__MPI)
    CALL mpi_abort(MPI_COMM_WORLD, ierr, info)
#endif
    !
    STOP 1
    !
    RETURN
    !
    !------------------------------------------------------------------------------------
END SUBROUTINE env_fft_error
!----------------------------------------------------------------------------------------
