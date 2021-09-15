!----------------------------------------------------------------------------------------
!
! Copyright (C) 2021 ENVIRON (www.quantum-environ.org)
! Copyright (C) 2001-2007 Quantum ESPRESSO (www.quantum-espresso.org)
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
! Authors: Modified by Edan Bainglass
!
!----------------------------------------------------------------------------------------
!
! Clocks should be started, read, stopped either on all processors, or
! only on one, but not half and half! For parallel debugging, uncomment:
!
! #define __TRACE
!
! See also comments in subroutine print_this_clock about parallel case
!
!----------------------------------------------------------------------------------------
!>
!! Time-printing utilities - Contains the following subroutines:
!!
!! init_clocks( go )      initialization - must be called first
!!                        go = .TRUE. : up to "maxclock" clocks can be started
!!                        go = .FALSE.: only clock #1 can be started
!! env_start_clock( label )   starts clock "label" (max 12 characters)
!!                        if "label" has never been started, initializes it
!!                        issues warning if "label" already started
!! env_stop_clock( label )    stops  clock "label"
!!                        issues warning if "label" is either not running
!!                        or has never been started
!! env_print_clock( label )   print cpu and wall time measured by clock "label"
!!                        clock "label" may be running or stopped
!!                        and remains in the same state
!!                        issues warning if "label" has never been started
!! and the following function (real(kind=dp):
!! get_clock( label )     return wall time measured by clock "label"
!!                        returns -1 if "label" has never been started
!!
!! All output and warnings are written to stdout
!!
!----------------------------------------------------------------------------------------
MODULE env_clocks
    !------------------------------------------------------------------------------------
    !
    USE env_base_io, ONLY: io
    USE env_kinds, ONLY: DP
    !
#if defined(__CUDA)
    USE cudafor
#endif
    !
    !------------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    SAVE
    !
    INTEGER, PARAMETER :: maxclock = 128
    REAL(DP), PARAMETER :: notrunning = -1.0_DP
    !
    REAL(DP) :: cputime(maxclock), t0cpu(maxclock)
    REAL(DP) :: walltime(maxclock), t0wall(maxclock)
    REAL(DP) :: gputime(maxclock)
    CHARACTER(LEN=12) :: clock_label(maxclock)
    INTEGER :: called(maxclock)
    INTEGER :: gpu_called(maxclock)
    !
    REAL(DP) :: mpi_per_thread = 1.0_DP

    INTEGER :: nclock = 0
    LOGICAL :: no
    !
#if defined(__TRACE)
    INTEGER :: trace_depth = 0
    !
    INTEGER :: max_print_depth = maxclock
    ! used to gauge the ammount of output. default: a very deep depth
    !
    INTEGER :: mpime
#endif
    !
#if defined(__CUDA)
    TYPE(cudaEvent) :: gpu_starts(maxclock), gpu_stops(maxclock)
#else
    INTEGER :: gpu_starts(maxclock), gpu_stops(maxclock) ! dummy values never used
#endif
    !
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
    INTERFACE
        !--------------------------------------------------------------------------------
        !
        !--------------------------------------------------------------------------------
        !>
        !!
        !--------------------------------------------------------------------------------
        FUNCTION env_f_wall() BIND(C, name="env_cclock") RESULT(t)
            !----------------------------------------------------------------------------
            !
            USE ISO_C_BINDING, ONLY: C_DOUBLE
            REAL(kind=C_DOUBLE) :: t
            !
            !----------------------------------------------------------------------------
        END FUNCTION env_f_wall
        !--------------------------------------------------------------------------------
        !>
        !!
        !--------------------------------------------------------------------------------
        FUNCTION env_f_tcpu() BIND(C, name="env_scnds") RESULT(t)
            !----------------------------------------------------------------------------
            !
            USE ISO_C_BINDING, ONLY: C_DOUBLE
            REAL(kind=C_DOUBLE) :: t
            !
            !----------------------------------------------------------------------------
        END FUNCTION env_f_tcpu
        !--------------------------------------------------------------------------------
        !
        !--------------------------------------------------------------------------------
    END INTERFACE
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE env_clocks
!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
!
! Clock routines (outside of module for easy global use)
!
!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
!>
!!
!----------------------------------------------------------------------------------------
SUBROUTINE env_start_clock(label)
    !------------------------------------------------------------------------------------
    !
    USE env_base_io, ONLY: io
    USE env_kinds, ONLY: DP
    !
    USE env_clocks, ONLY: nclock, clock_label, notrunning, no, maxclock, &
                          t0cpu, t0wall, env_f_wall, env_f_tcpu
    !
#if defined(__TRACE)
    USE env_clocks, ONLY: trace_depth, mpime, max_print_depth
#endif
    !
    !------------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=*) :: label
    !
    CHARACTER(LEN=12) :: label_
    INTEGER :: n
    !
    !------------------------------------------------------------------------------------
    ! Used to gauge the ammount of output
    !
#if defined(__TRACE)
    IF (trace_depth <= max_print_depth) THEN
        !
        WRITE (io%unit, '(I3," depth=",I2," env_start_clock ",A )') &
            mpime, trace_depth, label
        !
        FLUSH (io%unit)
    END IF
    !
    trace_depth = trace_depth + 1
#endif
    !
    IF (no .AND. (nclock == 1)) RETURN
    !
    label_ = TRIM(label)
    !
    DO n = 1, nclock
        !
        IF (clock_label(n) == label_) THEN
            !
            !----------------------------------------------------------------------------
            ! Found previously defined clock: check if not already started,
            ! store in t0cpu the starting time
            !
            IF (t0cpu(n) /= notrunning) THEN
                !
                WRITE (io%unit, '("env_start_clock: clock # ",I2," for ",A12, &
                                &" already started")') n, label_
                !
            ELSE
                t0cpu(n) = env_f_tcpu()
                t0wall(n) = env_f_wall()
            END IF
            !
            RETURN
            !
        END IF
        !
    END DO
    !
    !------------------------------------------------------------------------------------
    ! Clock not found : add new clock for given label
    !
    IF (nclock == maxclock) THEN
        !
        WRITE (io%unit, '("env_start_clock(",A,"): &
                        &Too many clocks! call ignored")') label
        !
    ELSE
        nclock = nclock + 1
        clock_label(nclock) = label_
        t0cpu(nclock) = env_f_tcpu()
        t0wall(nclock) = env_f_wall()
    END IF
    !
    !------------------------------------------------------------------------------------
END SUBROUTINE env_start_clock
!----------------------------------------------------------------------------------------
!>
!!
!----------------------------------------------------------------------------------------
SUBROUTINE env_stop_clock(label)
    !------------------------------------------------------------------------------------
    !
    USE env_base_io, ONLY: io
    USE env_kinds, ONLY: DP
    !
    USE env_clocks, ONLY: no, nclock, clock_label, cputime, walltime, &
                          notrunning, called, t0cpu, t0wall, env_f_wall, env_f_tcpu
    !
#if defined(__TRACE)
    USE env_clocks, ONLY: trace_depth, mpime, max_print_depth
#endif
    !
    !------------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=*) :: label
    !
    CHARACTER(LEN=12) :: label_
    INTEGER :: n
    !
    !------------------------------------------------------------------------------------
    ! Used to gauge the ammount of output
    !
#if defined(__TRACE)
    trace_depth = trace_depth - 1
    !
    IF (trace_depth <= max_print_depth) THEN
        !
        WRITE (io%unit, '(I3," depth=",I2," env_stop_clock ",A )') &
            mpime, trace_depth, label
        !
        FLUSH (io%unit)
    END IF
#endif
    !
    IF (no) RETURN
    !
    label_ = TRIM(label)
    !
    DO n = 1, nclock
        !
        IF (clock_label(n) == label_) THEN
            !
            !----------------------------------------------------------------------------
            ! Found previously defined clock : check if properly initialised,
            ! add elapsed time, increase the counter of calls
            !
            IF (t0cpu(n) == notrunning) THEN
                !
                WRITE (io%unit, '("env_stop_clock: clock # ",I2," for ",A12, &
                                &" not running")') n, label
                !
            ELSE
                cputime(n) = cputime(n) + env_f_tcpu() - t0cpu(n)
                walltime(n) = walltime(n) + env_f_wall() - t0wall(n)
                t0cpu(n) = notrunning
                t0wall(n) = notrunning
                called(n) = called(n) + 1
            END IF
            !
            RETURN
            !
        END IF
        !
    END DO
    !
    WRITE (io%unit, '("env_stop_clock: no clock for ",A12," found !")') label
    ! clock not found
    !
    !------------------------------------------------------------------------------------
END SUBROUTINE env_stop_clock
!----------------------------------------------------------------------------------------
!>
!!
!----------------------------------------------------------------------------------------
SUBROUTINE env_start_clock_gpu(label)
    !------------------------------------------------------------------------------------
    !
    USE env_base_io, ONLY: io
    USE env_kinds, ONLY: DP
    !
#if defined(__TRACE)
    USE env_clocks, ONLY: trace_depth, mpime, max_print_depth
#endif
    !
    USE env_clocks, ONLY: nclock, clock_label, notrunning, no, maxclock, &
                          t0cpu, t0wall, env_f_wall, env_f_tcpu, &
                          gputime, gpu_starts, gpu_stops
    !
#if defined(__CUDA)
    USE cudafor
#endif
    !
    !------------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=*) :: label
    !
    CHARACTER(LEN=12) :: label_
    INTEGER :: n, ierr
    !
    !------------------------------------------------------------------------------------
    ! Used to gauge the ammount of output
    !
#if defined(__TRACE)
    IF (trace_depth <= max_print_depth) THEN
        !
        WRITE (io%unit, '(I3," depth=",I2," env_start_clock ",A )') &
            mpime, trace_depth, label
        !
        FLUSH (io%unit)
    END IF
    !
    trace_depth = trace_depth + 1
#endif
    !
    IF (no .AND. (nclock == 1)) RETURN
    !
    label_ = TRIM(label)
    !
    DO n = 1, nclock
        !
        IF (clock_label(n) == label_) THEN
            !
            !----------------------------------------------------------------------------
            ! Found previously defined clock: check if not already started,
            ! store in t0cpu the starting time
            !
            IF (t0cpu(n) /= notrunning) THEN
                WRITE (io%unit, '("env_start_clock: clock # ",I2," for ",A12, &
                                & " already started")') n, label_
            ELSE
                !
#if defined(__CUDA)
                ierr = cudaEventRecord(gpu_starts(n), 0)
#endif
                !
                t0cpu(n) = env_f_tcpu()
                t0wall(n) = env_f_wall()
            END IF
            !
            RETURN
            !
        END IF
        !
    END DO
    !
    !------------------------------------------------------------------------------------
    ! Clock not found : add new clock for given label
    !
    IF (nclock == maxclock) THEN
        WRITE (io%unit, '("env_start_clock(",A,"): Too many clocks! call ignored")') label
    ELSE
        !
        nclock = nclock + 1
        clock_label(nclock) = label_
#if defined(__CUDA)
        ierr = cudaEventRecord(gpu_starts(nclock), 0)
#endif
        t0cpu(nclock) = env_f_tcpu()
        t0wall(nclock) = env_f_wall()
        !
    END IF
    !
    !------------------------------------------------------------------------------------
END SUBROUTINE env_start_clock_gpu
!----------------------------------------------------------------------------------------
!>
!!
!----------------------------------------------------------------------------------------
SUBROUTINE env_stop_clock_gpu(label)
    !------------------------------------------------------------------------------------
    !
    USE env_base_io, ONLY: io
    USE env_kinds, ONLY: DP
    !
#if defined(__TRACE)
    USE env_clocks, ONLY: trace_depth, mpime, max_print_depth
#endif
    !
    USE env_clocks, ONLY: no, nclock, clock_label, cputime, walltime, &
                          notrunning, called, t0cpu, t0wall, env_f_wall, env_f_tcpu, &
                          gpu_called, gputime, gpu_starts, gpu_stops
    !
#if defined(__CUDA)
    USE cudafor
#endif
    !
    !------------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=*) :: label
    !
    CHARACTER(LEN=12) :: label_
    INTEGER :: n, ierr
    REAL :: time
    !
    !------------------------------------------------------------------------------------
    !
#if defined(__TRACE)
    trace_depth = trace_depth - 1
    !
    !------------------------------------------------------------------------------------
    ! Used to gauge the ammount of output
    !
    IF (trace_depth <= max_print_depth) THEN
        !
        WRITE (io%unit, '(I3," depth=",I2," env_stop_clock ",A )') &
            mpime, trace_depth, label
        !
        FLUSH (io%unit)
    END IF
#endif
    !
    IF (no) RETURN
    !
    !------------------------------------------------------------------------------------
    ! Initialize time used in CUDA APIs if __CUDA is present.
    !
    time = 0.0
    !
    label_ = TRIM(label)
    !
    DO n = 1, nclock
        !
        IF (clock_label(n) == label_) THEN
            !
            !----------------------------------------------------------------------------
            ! Found previously defined clock : check if properly initialised,
            ! add elapsed time, increase the counter of calls
            !
            IF (t0cpu(n) == notrunning) THEN
                !
                WRITE (io%unit, &
                       '("env_stop_clock: clock # ",I2," for ",A12, " not running")') &
                    n, label
                !
            ELSE
                cputime(n) = cputime(n) + env_f_tcpu() - t0cpu(n)
                !
#if defined(__CUDA)
                ierr = cudaEventRecord(gpu_stops(n), 0)
                ierr = cudaEventSynchronize(gpu_stops(n))
                ierr = cudaEventElapsedTime(time, gpu_starts(n), gpu_stops(n))
#endif
                !
                gputime(n) = gputime(n) + time
                gpu_called(n) = gpu_called(n) + 1
                !
                walltime(n) = walltime(n) + env_f_wall() - t0wall(n)
                t0cpu(n) = notrunning
                t0wall(n) = notrunning
                called(n) = called(n) + 1
                !
            END IF
            !
            RETURN
            !
        END IF
        !
    END DO
    !
    WRITE (io%unit, '("env_stop_clock_gpu: no clock for ",A12," found !")') label
    ! clock not found
    !
    !------------------------------------------------------------------------------------
END SUBROUTINE env_stop_clock_gpu
!----------------------------------------------------------------------------------------
!>
!!
!----------------------------------------------------------------------------------------
FUNCTION env_get_cpu_and_wall(n) RESULT(t)
    !------------------------------------------------------------------------------------
    !
    USE env_kinds, ONLY: DP
    !
    USE env_clocks, ONLY: cputime, walltime, mpi_per_thread, notrunning, t0cpu, t0wall, &
                          env_f_wall, env_f_tcpu
    !
    !------------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    INTEGER :: n
    REAL(DP) :: t(2)
    !
    !------------------------------------------------------------------------------------
    !
    IF (t0cpu(n) == notrunning) THEN
        t(1) = cputime(n)
        t(2) = walltime(n)
    ELSE
        t(1) = cputime(n) + env_f_tcpu() - t0cpu(n)
        t(2) = walltime(n) + env_f_wall() - t0wall(n)
    END IF
    !
#if defined(PRINT_AVG_CPU_TIME_PER_THREAD)
    t(1) = t(1) * mpi_per_thread ! rescale the elapsed cpu time on a per-thread basis
#endif
    !
    !------------------------------------------------------------------------------------
END FUNCTION env_get_cpu_and_wall
!----------------------------------------------------------------------------------------
!>
!!
!----------------------------------------------------------------------------------------
SUBROUTINE env_print_clock(label)
    !------------------------------------------------------------------------------------
    !
    USE env_base_io, ONLY: io
    USE env_clocks, ONLY: nclock, clock_label, gpu_called
    !
    !------------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=*) :: label
    !
    CHARACTER(LEN=12) :: label_
    INTEGER :: n
    LOGICAL :: print_gpu
    !
    !------------------------------------------------------------------------------------
    !
    print_gpu = ANY(gpu_called > 0)
    !
    IF (label == ' ') THEN
        WRITE (io%unit, *)
        !
        DO n = 1, nclock
            !
            CALL env_print_this_clock(n)
            !
            IF (print_gpu) CALL env_print_this_clock_gpu(n)
            !
        END DO
        !
    ELSE
        label_ = TRIM(label)
        !
        DO n = 1, nclock
            !
            IF (clock_label(n) == label_) THEN
                !
                CALL env_print_this_clock(n)
                !
                IF (print_gpu) CALL env_print_this_clock_gpu(n)
                !
                EXIT
                !
            END IF
            !
        END DO
        !
    END IF
    !
    !------------------------------------------------------------------------------------
END SUBROUTINE env_print_clock
!----------------------------------------------------------------------------------------
!>
!!
!----------------------------------------------------------------------------------------
SUBROUTINE env_print_this_clock(n)
    !------------------------------------------------------------------------------------
    !
    USE env_base_io, ONLY: io
    USE env_kinds, ONLY: DP
    !
    USE env_clocks, ONLY: clock_label, cputime, walltime, mpi_per_thread, &
                          notrunning, called, t0cpu, t0wall, env_f_wall, env_f_tcpu
    !
    !------------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    INTEGER :: n, nday, nhour, nmin, nmax, mday, mhour, mmin
    REAL(DP) :: elapsed_cpu_time, elapsed_wall_time, nsec, msec
    !
    !------------------------------------------------------------------------------------
    !
    IF (t0cpu(n) == notrunning) THEN
        !
        !--------------------------------------------------------------------------------
        ! Clock stopped, print the stored value for the cpu time
        !
        elapsed_cpu_time = cputime(n)
        elapsed_wall_time = walltime(n)
    ELSE
        !
        !--------------------------------------------------------------------------------
        ! clock not stopped, print the current value of the cpu time
        !
        elapsed_cpu_time = cputime(n) + env_f_tcpu() - t0cpu(n)
        elapsed_wall_time = walltime(n) + env_f_wall() - t0wall(n)
        called(n) = called(n) + 1
    END IF
    !
! #define PRINT_AVG_CPU_TIME_PER_THREAD
#if defined(PRINT_AVG_CPU_TIME_PER_THREAD)
    elapsed_cpu_time = elapsed_cpu_time * mpi_per_thread
    ! rescale the elapsed cpu time on a per-thread basis
#endif
    !
    nmax = called(n)
    !
    IF (n == 1) THEN
        !
        !--------------------------------------------------------------------------------
        ! The first clock is written as days/hour/min/sec
        !
#if defined(__CLOCK_SECONDS)
        WRITE (io%unit, '(5X,A12," : ",F9.2,"s CPU ",F9.2,"s WALL"/)') &
            clock_label(n), elapsed_cpu_time, elapsed_wall_time
#else
        !
        nday = elapsed_cpu_time / 86400
        nsec = elapsed_cpu_time - 86400 * nday
        nhour = nsec / 3600
        nsec = nsec - 3600 * nhour
        nmin = nsec / 60
        nsec = nsec - 60 * nmin
        !
        !--------------------------------------------------------------------------------
        ! The first clock writes elapsed (wall) time as well
        !
        mday = elapsed_wall_time / 86400
        msec = elapsed_wall_time - 86400 * mday
        mhour = msec / 3600
        msec = msec - 3600 * mhour
        mmin = msec / 60
        msec = msec - 60 * mmin
        !
        IF (nday > 0) THEN
            !
            WRITE (io%unit, ADVANCE='no', &
                   FMT='(5X,A12," : ",1X,I2,"d",I2,"h",I2,"m CPU ")') &
                clock_label(n), nday, nhour, nmin
            !
        ELSE IF (nhour > 0) THEN
            !
            WRITE (io%unit, ADVANCE='no', &
                   FMT='(5X,A12," : ",4X,I2,"h",I2,"m CPU ")') &
                clock_label(n), nhour, nmin
            !
        ELSE IF (nmin > 0) THEN
            !
            WRITE (io%unit, ADVANCE='no', &
                   FMT='(5X,A12," : ",1X,I2,"m",F5.2,"s CPU ")') &
                clock_label(n), nmin, nsec
            !
        ELSE
            !
            WRITE (io%unit, ADVANCE='no', &
                   FMT='(5X,A12," : ",4X,F5.2,"s CPU ")') &
                clock_label(n), nsec
            !
        END IF
        !
        IF (mday > 0) THEN
            WRITE (io%unit, '(1X,I2,"d",I2,"h",I2,"m WALL"/)') mday, mhour, mmin
        ELSE IF (mhour > 0) THEN
            WRITE (io%unit, '(4X,I2,"h",I2,"m WALL"/)') mhour, mmin
        ELSE IF (mmin > 0) THEN
            WRITE (io%unit, '(1X,I2,"m",F5.2,"s WALL"/)') mmin, msec
        ELSE
            WRITE (io%unit, '(4X,F5.2,"s WALL"/)') msec
        END IF
#endif
        !
    ELSE IF (nmax == 1 .OR. t0cpu(n) /= notrunning) THEN
        !
        ! for clocks that have been called only once
        WRITE (io%unit, &
               '(5X,A12," : ",F9.2,"s CPU ",F9.2,"s WALL (",I8," calls)")') &
            clock_label(n), elapsed_cpu_time, elapsed_wall_time, nmax
        !
    ELSE IF (nmax == 0) THEN
        !
        ! for clocks that have never been called
        WRITE (io%unit, &
               '("print_this: clock # ",I2," for ",A12," never called !"/)') &
            n, clock_label(n)
        !
    ELSE
        !
        ! for all other clocks
        WRITE (io%unit, &
               '(5X,A12," : ",F9.2,"s CPU ",F9.2,"s WALL (",I8," calls)")') &
            clock_label(n), elapsed_cpu_time, elapsed_wall_time, nmax
        !
    END IF
    !
    !------------------------------------------------------------------------------------
END SUBROUTINE env_print_this_clock
!----------------------------------------------------------------------------------------
!>
!!
!----------------------------------------------------------------------------------------
SUBROUTINE env_print_this_clock_gpu(n)
    !------------------------------------------------------------------------------------
    !
    USE env_base_io, ONLY: io
    USE env_kinds, ONLY: DP
    USE env_clocks, ONLY: clock_label, gputime, gpu_called
    !
    !------------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    INTEGER :: n
    REAL(DP) :: elapsed_gpu_time
    INTEGER :: nday, nhour, nmin, nmax, mday, mhour, mmin
    !
    !------------------------------------------------------------------------------------
    !
    nmax = gpu_called(n)
    elapsed_gpu_time = gputime(n) / 1000.D0 ! GPU times are stored in ms
    !
    IF (nmax == 0) RETURN
    !
    IF (n == 1) THEN
        !
        ! the first clock is written as days/hour/min/sec
        WRITE (io%unit, '(5X,A12," : ",F9.2,"s GPU "/)') &
            clock_label(n), elapsed_gpu_time
        !
    ELSE
        !
        WRITE (io%unit, '(35X,F9.2,"s GPU  (",I8," calls)")') elapsed_gpu_time, nmax
        ! for all other clocks
        !
    END IF
    !
    !------------------------------------------------------------------------------------
END SUBROUTINE env_print_this_clock_gpu
!----------------------------------------------------------------------------------------
!>
!!
!----------------------------------------------------------------------------------------
REAL(DP) FUNCTION env_get_clock(label)
    !------------------------------------------------------------------------------------
    !
    USE env_kinds, ONLY: DP
    !
    USE env_clocks, ONLY: no, nclock, clock_label, walltime, notrunning, &
                          t0wall, t0cpu, env_f_wall
    !
    !------------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=*) :: label
    INTEGER :: n
    !
    !------------------------------------------------------------------------------------
    !
    IF (no) THEN
        !
        IF (label == clock_label(1)) THEN
            env_get_clock = env_f_wall()
        ELSE
            env_get_clock = notrunning
        END IF
        !
        RETURN
        !
    END IF
    !
    DO n = 1, nclock
        !
        IF (label == clock_label(n)) THEN
            !
            IF (t0cpu(n) == notrunning) THEN
                env_get_clock = walltime(n)
            ELSE
                env_get_clock = walltime(n) + env_f_wall() - t0wall(n)
            END IF
            !
            RETURN
            !
        END IF
        !
    END DO
    !
    env_get_clock = notrunning ! clock not found
    !
    !------------------------------------------------------------------------------------
END FUNCTION env_get_clock
!----------------------------------------------------------------------------------------
