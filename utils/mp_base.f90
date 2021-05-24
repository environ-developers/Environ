!----------------------------------------------------------------------------------------
!
! Copyright (C) 2018 ENVIRON (www.quantum-environment.org)
! Copyright (C) 2011 Quantum ESPRESSO group
!
!----------------------------------------------------------------------------------------
!
! This file is part of Environ version 2.02
!
! Environ 2.02 is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 2 of the License, or
! (at your option) any later version.
!
! Environ 2.02 is distributed in the hope that it will be useful,
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
!
#define __MSGSIZ_MAX 2000000
#define __BCAST_MSGSIZ_MAX_GPU huge(n)
!
!  Some implementation of MPI (OpenMPI) if it is not well tuned for the given
!  network hardware (InfiniBand) tend to lose performance or get stuck inside
!  collective routines if processors are not well synchronized
!  A barrier fixes the problem
!
#define __USE_BARRIER
!
!----------------------------------------------------------------------------------------
!>
!! Wrapper for MPI implementations that have problems with large messages
!!
!! In some MPI implementation the communication subsystem
!! crashes when message exceeds a given size, so we need
!! to break down MPI communications in smaller pieces
!!
!----------------------------------------------------------------------------------------
MODULE env_mp_base
    !------------------------------------------------------------------------------------
    !
    USE env_util_param
    !
#if defined(__CUDA)
    USE env_data_buffer, ONLY: buff_i => mp_buff_i_d, &
                               buff_r => mp_buff_r_d
#else
    USE env_data_buffer, ONLY: buff_i => mp_buff_i, &
                               buff_r => mp_buff_r
#endif
    !
    !------------------------------------------------------------------------------------
    !
    PRIVATE
    !
    PUBLIC :: env_allocate_mp_buffers, env_deallocate_mp_buffers, env_mp_synchronize, &
              env_bcast_integer, env_bcast_integer8, env_bcast_real, &
              env_bcast_complex, env_bcast_logical, env_reduce_base_integer, &
              env_reduce_base_integer8, env_reduce_base_real, &
              env_reduce_base_complex, env_reduce_base_real_to, &
              env_reduce_base_complex_to
    !
    !------------------------------------------------------------------------------------
CONTAINS
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    ! These routines allocate buffer spaces used in env_reduce_base_real.
    ! These should be in data_buffer.f90 but need to be here becouse size is
    ! depends on the __MSGSIZ_MAX definition
    !
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_allocate_mp_buffers()
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, PARAMETER :: maxb = __MSGSIZ_MAX
        !
        !--------------------------------------------------------------------------------
        !
        IF (.NOT. ALLOCATED(buff_r)) ALLOCATE (buff_r(maxb))
        !
        IF (.NOT. ALLOCATED(buff_i)) ALLOCATE (buff_i(maxb))
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_allocate_mp_buffers
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_deallocate_mp_buffers()
        !--------------------------------------------------------------------------------
        !
        DEALLOCATE (buff_r, buff_i)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_deallocate_mp_buffers
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_mp_synchronize(gid)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: gid
        !
#if defined(__MPI)&&defined(__USE_BARRIER)
        INTEGER :: ierr
        !
        CHARACTER(LEN=80) :: sub_name = 'env_mp_synchronize'
        !
        !--------------------------------------------------------------------------------
        !
        CALL mpi_barrier(gid, ierr)
        !
        IF (ierr /= 0) CALL env_errore(sub_name, ' error in mpi_barrier ', ierr)
#endif
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_mp_synchronize
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_bcast_integer(array, n, root, gid)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: n, root, gid
        !
#if defined(__CUDA)
        INTEGER, DEVICE :: array(n)
#else
        INTEGER :: array(n)
#endif
        !
#if defined(__MPI)
        INTEGER :: msgsiz_max = __MSGSIZ_MAX
        INTEGER :: nblk, blksiz, iblk, istart, ierr
        !
        CHARACTER(LEN=80) :: sub_name = 'env_bcast_integer'
        !
        !--------------------------------------------------------------------------------
        !
#if defined(__TRACE)
        WRITE (*, *) TRIM(sub_name)//' IN'
#endif
        !
        IF (n <= 0) GO TO 1
        !
#if defined(__USE_BARRIER)
        CALL env_mp_synchronize(gid)
#endif
        !
        IF (n <= msgsiz_max) THEN
            !
            CALL MPI_BCAST(array, n, MPI_INTEGER, root, gid, ierr)
            !
            IF (ierr /= 0) CALL env_errore(sub_name, ' error in mpi_bcast 1 ', ierr)
            !
        ELSE
            nblk = n / msgsiz_max
            blksiz = msgsiz_max
            !
            DO iblk = 1, nblk
                istart = (iblk - 1) * msgsiz_max + 1
                !
                CALL MPI_BCAST(array(istart), blksiz, MPI_INTEGER, root, gid, ierr)
                !
                IF (ierr /= 0) CALL env_errore(sub_name, ' error in mpi_bcast 2 ', ierr)
                !
            END DO
            !
            blksiz = MOD(n, msgsiz_max)
            !
            IF (blksiz > 0) THEN
                istart = nblk * msgsiz_max + 1
                !
                CALL MPI_BCAST(array(istart), blksiz, MPI_INTEGER, root, gid, ierr)
                !
                IF (ierr /= 0) CALL env_errore(sub_name, ' error in mpi_bcast 3 ', ierr)
                !
            END IF
            !
        END IF
        !
        GO TO 2 ! skip sync, already done by MPI call
        !
1       CONTINUE
        !
#if defined(__CUDA)
        ierr = cudaDeviceSynchronize()
#endif
        !
2       CONTINUE
        !
#if defined(__TRACE)
        WRITE (*, *) TRIM(sub_name)//' OUT'
#endif
#endif
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_bcast_integer
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_bcast_integer8(array, n, root, gid)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: n, root, gid
        !
#if defined(__CUDA)
        INTEGER(i8b), DEVICE :: array(n)
#else
        INTEGER(i8b) :: array(n)
#endif
        !
#if defined(__MPI)
        INTEGER :: msgsiz_max = __MSGSIZ_MAX
        INTEGER :: nblk, blksiz, iblk, istart, ierr
        !
        CHARACTER(LEN=80) :: sub_name = 'env_bcast_integer8'
        !
        !--------------------------------------------------------------------------------
        !
#if defined(__TRACE)
        WRITE (*, *) TRIM(sub_name)//' IN'
#endif
        !
        IF (n <= 0) GOTO 1
        !
#if defined(__USE_BARRIER)
        CALL env_mp_synchronize(gid)
#endif
        !
        IF (n <= msgsiz_max) THEN
            !
            CALL MPI_BCAST(array, n, MPI_INTEGER8, root, gid, ierr)
            !
            IF (ierr /= 0) CALL env_errore(sub_name, ' error in mpi_bcast 1 ', ierr)
            !
        ELSE
            nblk = n / msgsiz_max
            blksiz = msgsiz_max
            !
            DO iblk = 1, nblk
                istart = (iblk - 1) * msgsiz_max + 1
                !
                CALL MPI_BCAST(array(istart), blksiz, MPI_INTEGER8, root, gid, ierr)
                !
                IF (ierr /= 0) CALL env_errore(sub_name, ' error in mpi_bcast 2 ', ierr)
                !
            END DO
            !
            blksiz = MOD(n, msgsiz_max)
            !
            IF (blksiz > 0) THEN
                istart = nblk * msgsiz_max + 1
                !
                CALL MPI_BCAST(array(istart), blksiz, MPI_INTEGER8, root, gid, ierr)
                !
                IF (ierr /= 0) CALL env_errore(sub_name, ' error in mpi_bcast 3 ', ierr)
                !
            END IF
            !
        END IF
        !
        GO TO 2 ! skip sync, already done by MPI call
        !
1       CONTINUE
        !
#if defined(__CUDA)
        ierr = cudaDeviceSynchronize()
#endif
        !
2       CONTINUE
        !
#if defined(__TRACE)
        WRITE (*, *) TRIM(sub_name)//' OUT'
#endif
#endif
        !
        RETURN
        !
        !------------------------------------------------------------------------------!
    END SUBROUTINE env_bcast_integer8
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_bcast_real(array, n, root, gid)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: n, root, gid
        !
#if defined(__CUDA)
        REAL(DP), DEVICE :: array(n)
#else
        REAL(DP) :: array(n)
#endif
        !
#if defined(__MPI)
        INTEGER :: msgsiz_max = __BCAST_MSGSIZ_MAX_GPU
        INTEGER :: nblk, blksiz, iblk, istart, ierr
        !
        CHARACTER(LEN=80) :: sub_name = 'env_bcast_real'
        !
        !--------------------------------------------------------------------------------
        !
#if defined(__TRACE)
        WRITE (*, *) TRIM(sub_name)//' IN'
#endif
        IF (n <= 0) GO TO 1
        !
#if defined(__USE_BARRIER)
        CALL env_mp_synchronize(gid)
#endif
        !
        IF (n <= msgsiz_max) THEN
            !
            CALL MPI_BCAST(array, n, MPI_DOUBLE_PRECISION, root, gid, ierr)
            !
            IF (ierr /= 0) CALL env_errore(sub_name, ' error in mpi_bcast 1 ', ierr)
            !
        ELSE
            nblk = n / msgsiz_max
            blksiz = msgsiz_max
            !
            DO iblk = 1, nblk
                istart = (iblk - 1) * msgsiz_max + 1
                !
                CALL MPI_BCAST(array(istart), blksiz, MPI_DOUBLE_PRECISION, root, &
                               gid, ierr)
                !
                IF (ierr /= 0) CALL env_errore(sub_name, ' error in mpi_bcast 2 ', ierr)
                !
            END DO
            !
            blksiz = MOD(n, msgsiz_max)
            !
            IF (blksiz > 0) THEN
                istart = nblk * msgsiz_max + 1
                !
                CALL MPI_BCAST(array(istart), blksiz, MPI_DOUBLE_PRECISION, root, &
                               gid, ierr)
                !
                IF (ierr /= 0) CALL env_errore(sub_name, ' error in mpi_bcast 3 ', ierr)
                !
            END IF
            !
        END IF
        !
        GO TO 2 ! skip sync, already done by MPI call
        !
1       CONTINUE
        !
#if defined(__CUDA)
        ierr = cudaDeviceSynchronize()
#endif
        !
2       CONTINUE
        !
#if defined(__TRACE)
        WRITE (*, *) TRIM(sub_name)//' OUT'
#endif
#endif
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_bcast_real
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_bcast_complex(array, n, root, gid)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: n, root, gid
        !
#if defined(__CUDA)
        COMPLEX(DP), DEVICE :: array(n)
#else
        COMPLEX(DP) :: array(n)
#endif
        !
#if defined(__MPI)
        INTEGER :: msgsiz_max = __BCAST_MSGSIZ_MAX_GPU
        INTEGER :: nblk, blksiz, iblk, istart, ierr
        !
        CHARACTER(LEN=80) :: sub_name = 'env_bcast_complex'
        !
        !--------------------------------------------------------------------------------
        !
#if defined(__TRACE)
        WRITE (*, *) TRIM(sub_name)//' IN'
#endif
        IF (n <= 0) GO TO 1
        !
#if defined(__USE_BARRIER)
        CALL env_mp_synchronize(gid)
#endif
        !
        IF (n <= msgsiz_max) THEN
            !
            CALL MPI_BCAST(array, n, MPI_DOUBLE_COMPLEX, root, gid, ierr)
            !
            IF (ierr /= 0) CALL env_errore(sub_name, ' error in mpi_bcast 1 ', ierr)
            !
        ELSE
            nblk = n / msgsiz_max
            blksiz = msgsiz_max
            !
            DO iblk = 1, nblk
                istart = (iblk - 1) * msgsiz_max + 1
                !
                CALL MPI_BCAST(array(istart), blksiz, MPI_DOUBLE_COMPLEX, root, &
                               gid, ierr)
                !
                IF (ierr /= 0) CALL env_errore(sub_name, ' error in mpi_bcast 2 ', ierr)
                !
            END DO
            !
            blksiz = MOD(n, msgsiz_max)
            !
            IF (blksiz > 0) THEN
                istart = nblk * msgsiz_max + 1
                !
                CALL MPI_BCAST(array(istart), blksiz, MPI_DOUBLE_COMPLEX, root, &
                               gid, ierr)
                !
                IF (ierr /= 0) CALL env_errore(sub_name, ' error in mpi_bcast 3 ', ierr)
                !
            END IF
            !
        END IF
        !
        GO TO 2 ! skip sync, already done by MPI call
        !
1       CONTINUE
        !
#if defined(__CUDA)
        ierr = cudaDeviceSynchronize()
#endif
        !
2       CONTINUE
        !
#if defined(__TRACE)
        WRITE (*, *) TRIM(sub_name)//' OUT'
#endif
#endif
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_bcast_complex
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_bcast_logical(array, n, root, gid)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: n, root, gid
        !
#if defined(__CUDA)
        LOGICAL, DEVICE :: array(n)
#else
        LOGICAL :: array(n)
#endif
        !
#if defined(__MPI)
        INTEGER :: msgsiz_max = __MSGSIZ_MAX
        INTEGER :: nblk, blksiz, iblk, istart, ierr
        !
        CHARACTER(LEN=80) :: sub_name = 'env_bcast_logical'
        !
        !--------------------------------------------------------------------------------
        !
#if defined(__TRACE)
        WRITE (*, *) TRIM(sub_name)//' IN'
#endif
        !
        IF (n <= 0) GO TO 1
        !
#if defined(__USE_BARRIER)
        CALL env_mp_synchronize(gid)
#endif
        !
        IF (n <= msgsiz_max) THEN
            !
            CALL MPI_BCAST(array, n, MPI_LOGICAL, root, gid, ierr)
            !
            IF (ierr /= 0) CALL env_errore(sub_name, ' error in mpi_bcast 1 ', ierr)
            !
        ELSE
            nblk = n / msgsiz_max
            blksiz = msgsiz_max
            !
            DO iblk = 1, nblk
                istart = (iblk - 1) * msgsiz_max + 1
                !
                CALL MPI_BCAST(array(istart), blksiz, MPI_LOGICAL, root, gid, ierr)
                !
                IF (ierr /= 0) CALL env_errore(sub_name, ' error in mpi_bcast 2 ', ierr)
                !
            END DO
            !
            blksiz = MOD(n, msgsiz_max)
            !
            IF (blksiz > 0) THEN
                istart = nblk * msgsiz_max + 1
                !
                CALL MPI_BCAST(array(istart), blksiz, MPI_LOGICAL, root, gid, ierr)
                !
                IF (ierr /= 0) CALL env_errore(sub_name, ' error in mpi_bcast 3 ', ierr)
                !
            END IF
            !
        END IF
        !
        GO TO 2 ! skip sync, already done by MPI call
        !
1       CONTINUE
        !
#if defined(__CUDA)
        ierr = cudaDeviceSynchronize()
#endif
        !
2       CONTINUE
        !
#if defined(__TRACE)
        WRITE (*, *) TRIM(sub_name)//' OUT'
#endif
#endif
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_bcast_logical
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
    !
    ! "reduce"-like subroutines
    !
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
#if defined(__USE_INPLACE_MPI)
    !>
    !! Sums a distributed variable ps(dim) over the processors.
    !! This version uses a fixed-length buffer of appropriate (?) dim
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_reduce_base_integer(dim, ps, comm, root)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: dim, comm
        INTEGER, INTENT(IN) :: root
        ! if root < 0 perform a reduction to all procs
        ! if root >= 0 perform a reduction only to root proc
        !
#if defined(__CUDA)
        INTEGER, DEVICE :: ps(dim) ! array whose elements have to be reduced
#else
        INTEGER :: ps(dim) ! array whose elements have to be reduced
#endif
        !
#if defined(__MPI)
        INTEGER :: info, myid
        !
        CHARACTER(LEN=80) :: sub_name = 'env_reduce_base_integer'
        !
        !--------------------------------------------------------------------------------
        !
#if defined(__TRACE)
        WRITE (*, *) TRIM(sub_name)//' IN'
#endif
        !
        IF (dim <= 0) GO TO 1 ! go to the end of the subroutine
        !
#if defined(__USE_BARRIER)
        CALL env_mp_synchronize(comm)
#endif
        !
        IF (root >= 0) THEN
            !
            CALL mpi_comm_rank(comm, myid, info)
            !
            IF (info /= 0) CALL errore(sub_name, 'error in mpi_comm_rank', info)
            !
            IF (myid == root) THEN
                !
                CALL MPI_REDUCE(MPI_IN_PLACE, ps, dim, MPI_INTEGER, MPI_SUM, &
                                root, comm, info)
                !
                IF (info /= 0) CALL errore(sub_name, 'error in mpi_reduce 1', info)
                !
            ELSE
                !
                CALL MPI_REDUCE(ps, ps, dim, MPI_INTEGER, MPI_SUM, &
                                root, comm, info)
                !
                IF (info /= 0) CALL errore(sub_name, 'error in mpi_reduce 1', info)
                !
            END IF
            !
        ELSE
            !
            CALL MPI_ALLREDUCE(MPI_IN_PLACE, ps, dim, MPI_INTEGER, MPI_SUM, &
                               comm, info)
            !
            IF (info /= 0) CALL errore(sub_name, 'error in mpi_allreduce 1', info)
            !
        END IF
        !
        GO TO 2 ! skip sync, already done by MPI call
        !
1       CONTINUE
        !
#if defined(__CUDA)
        ierr = cudaDeviceSynchronize()
#endif
        !
2       CONTINUE
        !
#if defined(__TRACE)
        WRITE (*, *) TRIM(sub_name)//' OUT'
#endif
#endif
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_reduce_base_integer
    !------------------------------------------------------------------------------------
    !>
    !! Sums a distributed variable ps(dim) over the processors.
    !! This version uses a fixed-length buffer of appropriate (?) dim
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_reduce_base_integer8(dim, ps, comm, root)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: dim, comm
        INTEGER, INTENT(IN) :: root
        ! if root < 0 perform a reduction to all procs
        ! if root >= 0 perform a reduction only to root proc
        !
#if defined(__CUDA)
        INTEGER(i8b), DEVICE :: ps(dim) ! array whose elements have to be reduced
#else
        INTEGER(i8b) :: ps(dim) ! array whose elements have to be reduced
#endif
        !
#if defined(__MPI)
        INTEGER :: info, myid
        !
        CHARACTER(LEN=80) :: sub_name = 'env_reduce_base_integer8'
        !
        !--------------------------------------------------------------------------------
        !
#if defined(__TRACE)
        WRITE (*, *) TRIM(sub_name)//' IN'
#endif
        !
        IF (dim <= 0) GO TO 1 ! go to the end of the subroutine
        !
#if defined(__USE_BARRIER)
        CALL env_mp_synchronize(comm)
#endif
        !
        IF (root >= 0) THEN
            !
            CALL mpi_comm_rank(comm, myid, info)
            !
            IF (info /= 0) CALL env_errore(sub_name, 'error in mpi_comm_rank', info)
            !
            IF (myid == root) THEN
                !
                CALL MPI_REDUCE(MPI_IN_PLACE, ps, dim, MPI_INTEGER8, MPI_SUM, root, &
                                comm, info)
                !
                IF (info /= 0) CALL env_errore(sub_name, 'error in mpi_reduce 1', info)
                !
            ELSE
                !
                CALL MPI_REDUCE(ps, ps, dim, MPI_INTEGER8, MPI_SUM, root, comm, info)
                !
                IF (info /= 0) CALL env_errore(sub_name, 'error in mpi_reduce 1', info)
                !
            END IF
            !
        ELSE
            !
            CALL MPI_ALLREDUCE(MPI_IN_PLACE, ps, dim, MPI_INTEGER8, MPI_SUM, comm, info)
            !
            IF (info /= 0) CALL env_errore(sub_name, 'error in mpi_allreduce 1', info)
            !
        END IF
        !
        GOTO 2 ! skip sync, already done by MPI call
        !
1       CONTINUE
        !
#if defined(__CUDA)
        ierr = cudaDeviceSynchronize()
#endif
        !
2       CONTINUE
        !
#if defined(__TRACE)
        WRITE (*, *) TRIM(sub_name)//' OUT'
#endif
#endif
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_reduce_base_integer8
    !------------------------------------------------------------------------------------
    !>
    !! Sums a distributed variable ps(dim) over the processors.
    !! This version uses a fixed-length buffer of appropriate (?) dim
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_reduce_base_real(dim, ps, comm, root)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: dim, comm
        INTEGER, INTENT(IN) :: root
        ! if root < 0 perform a reduction to all procs
        ! if root >= 0 perform a reduction only to root proc
        !
#if defined(__CUDA)
        REAL(DP), DEVICE :: ps(dim) ! array whose elements have to be reduced
#else
        REAL(DP) :: ps(dim) ! array whose elements have to be reduced
#endif
        !
#if defined(__MPI)
        INTEGER :: info, myid
        !
        CHARACTER(LEN=80) :: sub_name = 'env_reduce_base_real'
        !
        !--------------------------------------------------------------------------------
        !
#if defined(__TRACE)
        WRITE (*, *) TRIM(sub_name)//' IN'
#endif
        !
        IF (dim <= 0) GO TO 1 ! go to the end of the subroutine
        !
#if defined(__USE_BARRIER)
        CALL env_mp_synchronize(comm)
#endif
        !
        IF (root >= 0) THEN
            !
            CALL mpi_comm_rank(comm, myid, info)
            !
            IF (info /= 0) CALL errore(sub_name, 'error in mpi_comm_rank', info)
            !
            IF (myid == root) THEN
                !
                CALL MPI_REDUCE(MPI_IN_PLACE, ps, dim, MPI_DOUBLE_PRECISION, MPI_SUM, &
                                root, comm, info)
                !
                IF (info /= 0) CALL errore(sub_name, 'error in mpi_reduce 1', info)
                !
            ELSE
                !
                CALL MPI_REDUCE(ps, ps, dim, MPI_DOUBLE_PRECISION, MPI_SUM, &
                                root, comm, info)
                !
                IF (info /= 0) CALL errore(sub_name, 'error in mpi_reduce 1', info)
                !
            END IF
            !
        ELSE
            !
            CALL MPI_ALLREDUCE(MPI_IN_PLACE, ps, dim, MPI_DOUBLE_PRECISION, MPI_SUM, &
                               comm, info)
            !
            IF (info /= 0) CALL errore(sub_name, 'error in mpi_allreduce 1', info)
            !
        END IF
        !
        GO TO 2 ! skip sync, already done by MPI call
        !
1       CONTINUE
        !
#if defined(__CUDA)
        ierr = cudaDeviceSynchronize()
#endif
        !
2       CONTINUE
        !
#if defined(__TRACE)
        WRITE (*, *) TRIM(sub_name)//' OUT'
#endif
#endif
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_reduce_base_real
    !------------------------------------------------------------------------------------
    !>
    !! Sums a distributed variable ps(dim) over the processors.
    !! This version uses a fixed-length buffer of appropriate (?) dim
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_reduce_base_complex(dim, ps, comm, root)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: dim, comm
        INTEGER, INTENT(IN) :: root
        ! if root < 0 perform a reduction to all procs
        ! if root >= 0 perform a reduction only to root proc
        !
#if defined(__CUDA)
        COMPLEX(DP), DEVICE :: ps(dim) ! array whose elements have to be reduced
#else
        COMPLEX(DP) :: ps(dim) ! array whose elements have to be reduced
#endif
        !
#if defined(__MPI)
        INTEGER :: info, myid
        !
        CHARACTER(LEN=80) :: sub_name = 'env_reduce_base_real'
        !
        !--------------------------------------------------------------------------------
        !
#if defined(__TRACE)
        WRITE (*, *) TRIM(sub_name)//' IN'
#endif
        !
        IF (dim <= 0) GO TO 1 ! go to the end of the subroutine
        !
#if defined(__USE_BARRIER)
        CALL env_mp_synchronize(comm)
#endif
        !
        IF (root >= 0) THEN
            !
            CALL mpi_comm_rank(comm, myid, info)
            !
            IF (info /= 0) CALL errore(sub_name, 'error in mpi_comm_rank', info)
            !
            IF (myid == root) THEN
                !
                CALL MPI_REDUCE(MPI_IN_PLACE, ps, dim, MPI_DOUBLE_COMPLEX, MPI_SUM, &
                                root, comm, info)
                !
                IF (info /= 0) CALL errore(sub_name, 'error in mpi_reduce 1', info)
                !
            ELSE
                !
                CALL MPI_REDUCE(ps, ps, dim, MPI_DOUBLE_COMPLEX, MPI_SUM, &
                                root, comm, info)
                !
                IF (info /= 0) CALL errore(sub_name, 'error in mpi_reduce 1', info)
                !
            END IF
            !
        ELSE
            !
            CALL MPI_ALLREDUCE(MPI_IN_PLACE, ps, dim, MPI_DOUBLE_COMPLEX, MPI_SUM, &
                               comm, info)
            !
            IF (info /= 0) CALL errore(sub_name, 'error in mpi_allreduce 1', info)
            !
        END IF
        !
        GO TO 2 ! skip sync, already done by MPI call
        !
1       CONTINUE
        !
#if defined(__CUDA)
        ierr = cudaDeviceSynchronize()
#endif
        !
2       CONTINUE
        !
#if defined(__TRACE)
        WRITE (*, *) TRIM(sub_name)//' OUT'
#endif
#endif
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_reduce_base_complex
    !------------------------------------------------------------------------------------
#else
    !>
    !! Sums a distributed variable ps(dim) over the processors.
    !! This version uses a fixed-length buffer of appropriate (?) dim
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_reduce_base_integer(dim, ps, comm, root)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: dim, comm
        INTEGER, INTENT(IN) :: root
        ! if root < 0 perform a reduction to all procs
        ! if root >= 0 perform a reduction only to root proc
        !
#if defined(__CUDA)
        INTEGER, DEVICE :: ps(dim) ! array whose elements have to be reduced
#else
        INTEGER :: ps(dim) ! array whose elements have to be reduced
#endif
        !
#if defined(__MPI)
        INTEGER :: info, n, nbuf, nproc, myid
        INTEGER, PARAMETER :: maxb = __MSGSIZ_MAX
        !
        CHARACTER(LEN=80) :: sub_name = 'env_reduce_base_integer'
        !
        !--------------------------------------------------------------------------------
        !
#if defined(__TRACE)
        WRITE (*, *) TRIM(sub_name)//' IN'
#endif
        !
        CALL mpi_comm_size(comm, nproc, info)
        !
        IF (info /= 0) CALL env_errore(sub_name, 'error in mpi_comm_size', info)
        !
        CALL mpi_comm_rank(comm, myid, info)
        !
        IF (info /= 0) CALL env_errore(sub_name, 'error in mpi_comm_rank', info)
        !
        IF (dim <= 0 .OR. nproc <= 1) GO TO 1 ! go to the end of the subroutine
        !
#if defined(__USE_BARRIER)
        CALL env_mp_synchronize(comm)
#endif
        !
        nbuf = dim / maxb
        !
        DO n = 1, nbuf
            !
            IF (root >= 0) THEN
                !
                CALL MPI_REDUCE(ps(1 + (n - 1) * maxb), buff_i, maxb, &
                                MPI_INTEGER, MPI_SUM, root, comm, info)
                !
                IF (info /= 0) CALL env_errore(sub_name, 'error in mpi_reduce 1', info)
                !
            ELSE
                !
                CALL MPI_ALLREDUCE(ps(1 + (n - 1) * maxb), buff_i, maxb, &
                                   MPI_INTEGER, MPI_SUM, comm, info)
                !
                IF (info /= 0) &
                    CALL env_errore(sub_name, 'error in mpi_allreduce 1', info)
                !
            END IF
            !
            IF (root < 0 .OR. root == myid) THEN
                !
#if defined(__CUDA)
                info = cudaMemcpy(ps((1 + (n - 1) * maxb)), buff_i(1), maxb, &
                                  cudaMemcpyDeviceToDevice)
                !
                IF (info /= 0) CALL env_errore(sub_name, 'error in cudaMemcpy ', info)
#else
                !
                ps((1 + (n - 1) * maxb):(n * maxb)) = buff_i(1:maxb)
#endif
            END IF
            !
        END DO
        !
        !--------------------------------------------------------------------------------
        ! Possible remaining elements < maxb
        !
        IF ((dim - nbuf * maxb) > 0) THEN
            !
            IF (root >= 0) THEN
                !
                CALL MPI_REDUCE(ps(1 + nbuf * maxb), buff_i, &
                                (dim - nbuf * maxb), MPI_INTEGER, MPI_SUM, &
                                root, comm, info)
                !
                IF (info /= 0) CALL env_errore(sub_name, 'error in mpi_reduce 2', info)
                !
            ELSE
                !
                CALL MPI_ALLREDUCE(ps(1 + nbuf * maxb), buff_i, &
                                   (dim - nbuf * maxb), MPI_INTEGER, MPI_SUM, &
                                   comm, info)
                !
                IF (info /= 0) &
                    CALL env_errore(sub_name, 'error in mpi_allreduce 2', info)
                !
            END IF
            !
            IF (root < 0 .OR. root == myid) THEN
                !
#if defined(__CUDA)
                info = cudaMemcpy(ps((1 + nbuf * maxb)), buff_i(1), maxb, &
                                  cudaMemcpyDeviceToDevice)
                !
                IF (info /= 0) CALL env_errore(sub_name, 'error in cudaMemcpy ', info)
#else
                !
                ps((1 + nbuf * maxb):dim) = buff_i(1:(dim - nbuf * maxb))
#endif
            END IF
            !
        END IF
        !
        ! DO NOT skip sync, last step may not be an MPI Call
        !
1       CONTINUE
        !
#if defined(__CUDA)
        info = cudaDeviceSynchronize()
#endif
        !
2       CONTINUE
        !
#if defined(__TRACE)
        WRITE (*, *) TRIM(sub_name)//' OUT'
#endif
#endif
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_reduce_base_integer
    !------------------------------------------------------------------------------------
    !>
    !! Sums a distributed variable ps(dim) over the processors.
    !! This version uses a fixed-length buffer of appropriate (?) dim
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_reduce_base_integer8(dim, ps, comm, root)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: dim, comm
        INTEGER, INTENT(IN) :: root
        ! if root < 0 perform a reduction to all procs
        ! if root >= 0 perform a reduction only to root proc
        !
#if defined(__CUDA)
        INTEGER(i8b), DEVICE :: ps(dim) ! array whose elements have to be reduced
#else
        INTEGER(i8b) :: ps(dim) ! array whose elements have to be reduced
#endif
        !
#if defined(__MPI)
        INTEGER :: info, n, nbuf, nproc, myid
        INTEGER, PARAMETER :: maxb = __MSGSIZ_MAX
        !
        CHARACTER(LEN=80) :: sub_name = 'env_reduce_base_integer8'
        !
        !--------------------------------------------------------------------------------
        !
#if defined(__TRACE)
        WRITE (*, *) TRIM(sub_name)//' IN'
#endif
        !
        CALL mpi_comm_size(comm, nproc, info)
        !
        IF (info /= 0) CALL env_errore(sub_name, 'error in mpi_comm_size', info)
        !
        CALL mpi_comm_rank(comm, myid, info)
        !
        IF (info /= 0) CALL env_errore(sub_name, 'error in mpi_comm_rank', info)
        !
        IF (dim <= 0 .OR. nproc <= 1) GOTO 1 ! go to the end of the subroutine
        !
#if defined(__USE_BARRIER)
        CALL env_mp_synchronize(comm)
#endif
        !
        nbuf = dim / maxb
        !
        DO n = 1, nbuf
            !
            IF (root >= 0) THEN
                !
                CALL MPI_REDUCE(ps(1 + (n - 1) * maxb), buff_i, maxb, &
                                MPI_INTEGER8, MPI_SUM, root, comm, info)
                !
                IF (info /= 0) CALL env_errore(sub_name, 'error in mpi_reduce 1', info)
                !
            ELSE
                !
                CALL MPI_ALLREDUCE(ps(1 + (n - 1) * maxb), buff_i, maxb, &
                                   MPI_INTEGER8, MPI_SUM, comm, info)
                !
                IF (info /= 0) &
                    CALL env_errore(sub_name, 'error in mpi_allreduce 1', info)
                !
            END IF
            !
            IF (root < 0 .OR. root == myid) THEN
                !
#if defined(__CUDA)
                info = cudaMemcpy(ps((1 + (n - 1) * maxb)), buff_i(1), maxb, &
                                  cudaMemcpyDeviceToDevice)
                !
                IF (info /= 0) CALL env_errore(sub_name, 'error in cudaMemcpy ', info)
#else
                !
                ps((1 + (n - 1) * maxb):(n * maxb)) = buff_i(1:maxb)
#endif
            END IF
            !
        END DO
        !
        !--------------------------------------------------------------------------------
        ! Possible remaining elements < maxb
        !
        IF ((dim - nbuf * maxb) > 0) THEN
            !
            IF (root >= 0) THEN
                !
                CALL MPI_REDUCE(ps(1 + nbuf * maxb), buff_i, (dim - nbuf * maxb), &
                                MPI_INTEGER8, MPI_SUM, root, comm, info)
                !
                IF (info /= 0) CALL env_errore(sub_name, 'error in mpi_reduce 2', info)
                !
            ELSE
                !
                CALL MPI_ALLREDUCE(ps(1 + nbuf * maxb), buff_i, (dim - nbuf * maxb), &
                                   MPI_INTEGER8, MPI_SUM, comm, info)
                !
                IF (info /= 0) &
                    CALL env_errore(sub_name, 'error in mpi_allreduce 2', info)
                !
            END IF
            !
            IF (root < 0 .OR. root == myid) THEN
                !
#if defined(__CUDA)
                info = cudaMemcpy(ps((1 + nbuf * maxb)), buff_i(1), maxb, &
                                  cudaMemcpyDeviceToDevice)
                !
                IF (info /= 0) CALL env_errore(sub_name, 'error in cudaMemcpy ', info)
#else
                !
                ps((1 + nbuf * maxb):dim) = buff_i(1:(dim - nbuf * maxb))
#endif
            END IF
            !
        END IF
        !
        ! DO NOT skip sync, last step may not be an MPI Call
        !
1       CONTINUE
        !
#if defined(__CUDA)
        info = cudaDeviceSynchronize()
#endif
        !
2       CONTINUE
        !
#if defined(__TRACE)
        WRITE (*, *) TRIM(sub_name)//' OUT'
#endif
        !
#endif
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_reduce_base_integer8
    !------------------------------------------------------------------------------------
    !>
    !! Sums a distributed variable ps(dim) over the processors.
    !! This version uses a fixed-length buffer of appropriate (?) dim
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_reduce_base_real(dim, ps, comm, root)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: dim, comm
        INTEGER, INTENT(IN) :: root
        ! if root < 0 perform a reduction to all procs
        ! if root >= 0 perform a reduction only to root proc
        !
#if defined(__CUDA)
        REAL(DP), DEVICE :: ps(dim) ! array whose elements have to be reduced
#else
        REAL(DP) :: ps(dim) ! array whose elements have to be reduced
#endif
        !
#if defined(__MPI)
        INTEGER :: info, n, nbuf, nproc, myid
        INTEGER, PARAMETER :: maxb = __MSGSIZ_MAX
        !
        CHARACTER(LEN=80) :: sub_name = 'env_reduce_base_real'
        !
        !--------------------------------------------------------------------------------
        !
#if defined(__TRACE)
        WRITE (*, *) TRIM(sub_name)//' IN'
#endif
        !
        CALL mpi_comm_size(comm, nproc, info)
        !
        IF (info /= 0) CALL env_errore(sub_name, 'error in mpi_comm_size', info)
        !
        CALL mpi_comm_rank(comm, myid, info)
        !
        IF (info /= 0) CALL env_errore(sub_name, 'error in mpi_comm_rank', info)
        !
        IF (dim <= 0 .OR. nproc <= 1) GO TO 1 ! go to the end of the subroutine
        !
#if defined(__USE_BARRIER)
        CALL env_mp_synchronize(comm)
#endif
        !
        nbuf = dim / maxb
        !
        DO n = 1, nbuf
            !
            IF (root >= 0) THEN
                !
                CALL MPI_REDUCE(ps(1 + (n - 1) * maxb), buff_r, maxb, &
                                MPI_DOUBLE_PRECISION, MPI_SUM, root, comm, info)
                !
                IF (info /= 0) CALL env_errore(sub_name, 'error in mpi_reduce 1', info)
                !
            ELSE
                !
                CALL MPI_ALLREDUCE(ps(1 + (n - 1) * maxb), buff_r, maxb, &
                                   MPI_DOUBLE_PRECISION, MPI_SUM, comm, info)
                !
                IF (info /= 0) &
                    CALL env_errore(sub_name, 'error in mpi_allreduce 1', info)
                !
            END IF
            !
            IF (root < 0 .OR. root == myid) THEN
                !
#if defined(__CUDA)
                info = cudaMemcpy(ps((1 + (n - 1) * maxb)), buff_r(1), maxb, &
                                  cudaMemcpyDeviceToDevice)
                !
                IF (info /= 0) CALL env_errore(sub_name, 'error in cudaMemcpy ', info)
#else
                !
                ps((1 + (n - 1) * maxb):(n * maxb)) = buff_r(1:maxb)
#endif
            END IF
            !
        END DO
        !
        !--------------------------------------------------------------------------------
        ! Possible remaining elements < maxb
        !
        IF ((dim - nbuf * maxb) > 0) THEN
            !
            IF (root >= 0) THEN
                !
                CALL MPI_REDUCE(ps(1 + nbuf * maxb), buff_r, &
                                (dim - nbuf * maxb), MPI_DOUBLE_PRECISION, &
                                MPI_SUM, root, comm, info)
                !
                IF (info /= 0) CALL env_errore(sub_name, 'error in mpi_reduce 2', info)
                !
            ELSE
                !
                CALL MPI_ALLREDUCE(ps(1 + nbuf * maxb), buff_r, &
                                   (dim - nbuf * maxb), MPI_DOUBLE_PRECISION, &
                                   MPI_SUM, comm, info)
                !
                IF (info /= 0) &
                    CALL env_errore(sub_name, 'error in mpi_allreduce 2', info)
                !
            END IF
            !
            IF (root < 0 .OR. root == myid) THEN
                !
#if defined(__CUDA)
                info = cudaMemcpy(ps((1 + nbuf * maxb)), buff_r(1), maxb, &
                                  cudaMemcpyDeviceToDevice)
                !
                IF (info /= 0) CALL env_errore(sub_name, 'error in cudaMemcpy ', info)
#else
                !
                ps((1 + nbuf * maxb):dim) = buff_r(1:(dim - nbuf * maxb))
#endif
            END IF
            !
        END IF
        !
        ! DO NOT skip sync, last step may not be an MPI Call
        !
1       CONTINUE
        !
#if defined(__CUDA)
        info = cudaDeviceSynchronize()
#endif
        !
2       CONTINUE
        !
#if defined(__TRACE)
        WRITE (*, *) TRIM(sub_name)//' OUT'
#endif
        !
#endif
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_reduce_base_real
    !------------------------------------------------------------------------------------
    !>
    !! Sums a distributed variable ps(dim) over the processors.
    !! This version uses a fixed-length buffer of appropriate (?) dim
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_reduce_base_complex(dim, ps, comm, root)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: dim, comm
        INTEGER, INTENT(IN) :: root
        ! if root < 0 perform a reduction to all procs
        ! if root >= 0 perform a reduction only to root proc
        !
#if defined(__CUDA)
        COMPLEX(DP), DEVICE :: ps(dim) ! array whose elements have to be reduced
#else
        COMPLEX(DP) :: ps(dim) ! array whose elements have to be reduced
#endif
        !
#if defined(__MPI)
        INTEGER :: info, n, nbuf, nproc, myid
        INTEGER, PARAMETER :: maxb = __MSGSIZ_MAX
        !
        CHARACTER(LEN=80) :: sub_name = 'env_reduce_base_real'
        !
        !--------------------------------------------------------------------------------
        !
#if defined(__TRACE)
        WRITE (*, *) TRIM(sub_name)//' IN'
#endif
        !
        CALL mpi_comm_size(comm, nproc, info)
        !
        IF (info /= 0) CALL env_errore(sub_name, 'error in mpi_comm_size', info)
        !
        CALL mpi_comm_rank(comm, myid, info)
        !
        IF (info /= 0) CALL env_errore(sub_name, 'error in mpi_comm_rank', info)
        !
        IF (dim <= 0 .OR. nproc <= 1) GO TO 1 ! go to the end of the subroutine
        !
#if defined(__USE_BARRIER)
        CALL env_mp_synchronize(comm)
#endif
        !
        nbuf = dim / maxb
        !
        DO n = 1, nbuf
            !
            IF (root >= 0) THEN
                !
                CALL MPI_REDUCE(ps(1 + (n - 1) * maxb), buff_r, maxb, &
                                MPI_DOUBLE_COMPLEX, MPI_SUM, root, comm, info)
                !
                IF (info /= 0) CALL env_errore(sub_name, 'error in mpi_reduce 1', info)
                !
            ELSE
                !
                CALL MPI_ALLREDUCE(ps(1 + (n - 1) * maxb), buff_r, maxb, &
                                   MPI_DOUBLE_COMPLEX, MPI_SUM, comm, info)
                !
                IF (info /= 0) &
                    CALL env_errore(sub_name, 'error in mpi_allreduce 1', info)
                !
            END IF
            !
            IF (root < 0 .OR. root == myid) THEN
                !
#if defined(__CUDA)
                info = cudaMemcpy(ps((1 + (n - 1) * maxb)), buff_r(1), maxb, &
                                  cudaMemcpyDeviceToDevice)
                !
                IF (info /= 0) CALL env_errore(sub_name, 'error in cudaMemcpy ', info)
#else
                !
                ps((1 + (n - 1) * maxb):(n * maxb)) = buff_r(1:maxb)
#endif
            END IF
            !
        END DO
        !
        !--------------------------------------------------------------------------------
        ! Possible remaining elements < maxb
        !
        IF ((dim - nbuf * maxb) > 0) THEN
            !
            IF (root >= 0) THEN
                !
                CALL MPI_REDUCE(ps(1 + nbuf * maxb), buff_r, &
                                (dim - nbuf * maxb), MPI_DOUBLE_COMPLEX, &
                                MPI_SUM, root, comm, info)
                !
                IF (info /= 0) CALL env_errore(sub_name, 'error in mpi_reduce 2', info)
                !
            ELSE
                !
                CALL MPI_ALLREDUCE(ps(1 + nbuf * maxb), buff_r, &
                                   (dim - nbuf * maxb), MPI_DOUBLE_COMPLEX, &
                                   MPI_SUM, comm, info)
                !
                IF (info /= 0) &
                    CALL env_errore(sub_name, 'error in mpi_allreduce 2', info)
                !
            END IF
            !
            IF (root < 0 .OR. root == myid) THEN
                !
#if defined(__CUDA)
                info = cudaMemcpy(ps((1 + nbuf * maxb)), buff_r(1), maxb, &
                                  cudaMemcpyDeviceToDevice)
                !
                IF (info /= 0) CALL env_errore(sub_name, 'error in cudaMemcpy ', info)
#else
                !
                ps((1 + nbuf * maxb):dim) = buff_r(1:(dim - nbuf * maxb))
#endif
            END IF
            !
        END IF
        !
        ! DO NOT skip sync, last step may not be an MPI Call
        !
1       CONTINUE
        !
#if defined(__CUDA)
        info = cudaDeviceSynchronize()
#endif
        !
2       CONTINUE
        !
#if defined(__TRACE)
        WRITE (*, *) TRIM(sub_name)//' OUT'
#endif
        !
#endif
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_reduce_base_complex
    !------------------------------------------------------------------------------------
#endif
    !>
    !! Sums a distributed variable ps(dim) over the processors,
    !! and store the results in variable psout.
    !! This version uses a fixed-length buffer of appropriate (?) length
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_reduce_base_real_to(dim, ps, psout, comm, root)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: dim, comm
        INTEGER, INTENT(IN) :: root
        ! if root < 0 perform a reduction to all procs
        ! if root >= 0 perform a reduction only to root proc
        !
#if defined(__CUDA)
        REAL(DP), DEVICE, INTENT(IN) :: ps(dim)
        !
        REAL(DP), DEVICE, INTENT(OUT) :: psout(dim)
#else
        REAL(DP), INTENT(IN) :: ps(dim)
        !
        REAL(DP), INTENT(OUT) :: psout(dim)
#endif
        !
#if defined(__MPI)
        INTEGER :: info, n, nbuf, nproc, myid
        INTEGER, PARAMETER :: maxb = __MSGSIZ_MAX
        !
        CHARACTER(LEN=80) :: sub_name = 'env_reduce_base_real_to'
        !
        !--------------------------------------------------------------------------------
        !
#if defined(__TRACE)
        WRITE (*, *) TRIM(sub_name)//' IN'
#endif
        !
        CALL mpi_comm_size(comm, nproc, info)
        !
        IF (info /= 0) CALL env_errore(sub_name, 'error in mpi_comm_size', info)
        !
        CALL mpi_comm_rank(comm, myid, info)
        !
        IF (info /= 0) CALL env_errore(sub_name, 'error in mpi_comm_rank', info)
        !
        IF (dim > 0 .AND. nproc <= 1) psout = ps
        !
        IF (dim <= 0 .OR. nproc <= 1) GO TO 1
        ! go to the sync and later to the end of the subroutine
        !
#if defined(__USE_BARRIER)
        CALL env_mp_synchronize(comm)
#endif
        !
        nbuf = dim / maxb
        !
        DO n = 1, nbuf
            !
            IF (root >= 0) THEN
                !
                CALL MPI_REDUCE(ps(1 + (n - 1) * maxb), &
                                psout(1 + (n - 1) * maxb), maxb, &
                                MPI_DOUBLE_PRECISION, MPI_SUM, root, comm, info)
                !
                IF (info /= 0) CALL env_errore(sub_name, 'error in mpi_reduce 1', info)
                !
            ELSE
                !
                CALL MPI_ALLREDUCE(ps(1 + (n - 1) * maxb), &
                                   psout(1 + (n - 1) * maxb), maxb, &
                                   MPI_DOUBLE_PRECISION, MPI_SUM, comm, info)
                !
                IF (info /= 0) &
                    CALL env_errore(sub_name, 'error in mpi_allreduce 1', info)
                !
            END IF
            !
        END DO
        !
        !--------------------------------------------------------------------------------
        ! Possible remaining elements < maxb
        !
        IF ((dim - nbuf * maxb) > 0) THEN
            !
            IF (root >= 0) THEN
                !
                CALL MPI_REDUCE(ps(1 + nbuf * maxb), psout(1 + nbuf * maxb), &
                                (dim - nbuf * maxb), MPI_DOUBLE_PRECISION, MPI_SUM, &
                                root, comm, info)
                !
                IF (info /= 0) CALL env_errore(sub_name, 'error in mpi_reduce 2', info)
                !
            ELSE
                !
                CALL MPI_ALLREDUCE(ps(1 + nbuf * maxb), psout(1 + nbuf * maxb), &
                                   (dim - nbuf * maxb), MPI_DOUBLE_PRECISION, MPI_SUM, &
                                   comm, info)
                !
                IF (info /= 0) CALL &
                    env_errore(sub_name, 'error in mpi_allreduce 2', info)
                !
            END IF
            !
        END IF
        !
        GO TO 2 ! skip sync, already done by MPI
        !
1       CONTINUE
        !
#if defined(__CUDA)
        info = cudaDeviceSynchronize()
#endif
        !
2       CONTINUE
        !
#if defined(__TRACE)
        WRITE (*, *) TRIM(sub_name)//' OUT'
#endif
#endif
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_reduce_base_real_to
    !------------------------------------------------------------------------------------
    !>
    !! Sums a distributed variable ps(dim) over the processors,
    !! and store the results in variable psout.
    !! This version uses a fixed-length buffer of appropriate (?) length
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_reduce_base_complex_to(dim, ps, psout, comm, root)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: dim, comm
        INTEGER, INTENT(IN) :: root
        ! if root < 0 perform a reduction to all procs
        ! if root >= 0 perform a reduction only to root proc
        !
#if defined(__CUDA)
        COMPLEX(DP), DEVICE, INTENT(IN) :: ps(dim)
        !
        COMPLEX(DP), DEVICE, INTENT(OUT) :: psout(dim)
#else
        COMPLEX(DP), INTENT(IN) :: ps(dim)
        !
        COMPLEX(DP), INTENT(OUT) :: psout(dim)
#endif
        !
#if defined(__MPI)
        INTEGER :: info, n, nbuf, nproc, myid
        INTEGER, PARAMETER :: maxb = __MSGSIZ_MAX
        !
        CHARACTER(LEN=80) :: sub_name = 'env_reduce_base_real_to'
        !
        !--------------------------------------------------------------------------------
        !
#if defined(__TRACE)
        WRITE (*, *) TRIM(sub_name)//' IN'
#endif
        !
        CALL mpi_comm_size(comm, nproc, info)
        !
        IF (info /= 0) CALL env_errore(sub_name, 'error in mpi_comm_size', info)
        !
        CALL mpi_comm_rank(comm, myid, info)
        !
        IF (info /= 0) CALL env_errore(sub_name, 'error in mpi_comm_rank', info)
        !
        IF (dim > 0 .AND. nproc <= 1) psout = ps
        !
        IF (dim <= 0 .OR. nproc <= 1) GO TO 1
        ! go to the sync and later to the end of the subroutine
        !
#if defined(__USE_BARRIER)
        CALL env_mp_synchronize(comm)
#endif
        !
        nbuf = dim / maxb
        !
        DO n = 1, nbuf
            !
            IF (root >= 0) THEN
                !
                CALL MPI_REDUCE(ps(1 + (n - 1) * maxb), &
                                psout(1 + (n - 1) * maxb), maxb, &
                                MPI_DOUBLE_COMPLEX, MPI_SUM, root, comm, info)
                !
                IF (info /= 0) CALL env_errore(sub_name, 'error in mpi_reduce 1', info)
                !
            ELSE
                !
                CALL MPI_ALLREDUCE(ps(1 + (n - 1) * maxb), &
                                   psout(1 + (n - 1) * maxb), maxb, &
                                   MPI_DOUBLE_COMPLEX, MPI_SUM, comm, info)
                !
                IF (info /= 0) &
                    CALL env_errore(sub_name, 'error in mpi_allreduce 1', info)
                !
            END IF
            !
        END DO
        !
        !--------------------------------------------------------------------------------
        ! Possible remaining elements < maxb
        !
        IF ((dim - nbuf * maxb) > 0) THEN
            !
            IF (root >= 0) THEN
                !
                CALL MPI_REDUCE(ps(1 + nbuf * maxb), psout(1 + nbuf * maxb), &
                                (dim - nbuf * maxb), MPI_DOUBLE_COMPLEX, MPI_SUM, &
                                root, comm, info)
                !
                IF (info /= 0) CALL env_errore(sub_name, 'error in mpi_reduce 2', info)
                !
            ELSE
                !
                CALL MPI_ALLREDUCE(ps(1 + nbuf * maxb), psout(1 + nbuf * maxb), &
                                   (dim - nbuf * maxb), MPI_DOUBLE_COMPLEX, MPI_SUM, &
                                   comm, info)
                !
                IF (info /= 0) &
                    CALL env_errore(sub_name, 'error in mpi_allreduce 2', info)
                !
            END IF
            !
        END IF
        !
        GO TO 2 ! skip sync, already done by MPI
        !
1       CONTINUE
        !
#if defined(__CUDA)
        info = cudaDeviceSynchronize()
#endif
        !
2       CONTINUE
        !
#if defined(__TRACE)
        WRITE (*, *) TRIM(sub_name)//' OUT'
#endif
#endif
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_reduce_base_complex_to
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE env_mp_base
!----------------------------------------------------------------------------------------