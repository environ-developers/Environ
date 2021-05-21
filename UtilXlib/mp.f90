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
! Authors: Carlo Cavazzoni, Samuel Ponce, Pietro Bonfa'
!
!----------------------------------------------------------------------------------------
!>
!! This module contains interfaces to most low-level MPI operations:
!! initialization and stopping, broadcast, parallel sum, etc.
!!
!! Before hacking on the CUDA part remember that:
!!
!! 1. all mp_* interface should be blocking with respect to both MPI and CUDA.
!!    MPI will only wait for completion on the default stream therefore device
!!    synchronization must be enforced.
!!
!! 2. Host -> device memory copies of a memory block of 64 KB or less are
!!    asynchronous in the sense that they may return before the data is actually
!!    available on the GPU. However, the user is still free to change the buffer
!!    as soon as those calls return with no ill effects.
!!    (https://devtalk.nvidia.com/default/topic/471866/cuda-programming-and-performance/host-device-memory-copies-up-to-64-kb-are-asynchronous/)
!!
!! 3. For transfers from device to either pageable or pinned host memory,
!!    the function env_returns only once the copy has completed.
!!
!! 4. GPU synchronization is always enforced even if no communication takes place.
!!
!! #TODO COMPLEX mp routines differ from QE's implementation and need to be tested!!!
!!       Environ currently (2.0) does not broadcast any COMPLEX quantities
!!
!----------------------------------------------------------------------------------------
MODULE env_mp
    !------------------------------------------------------------------------------------
    !
    USE env_util_param
    !
    USE env_mp_base
    !
#if defined(__CUDA)
    USE cudafor
#endif
    !
    !------------------------------------------------------------------------------------
    !
    INTERFACE env_mp_bcast
        MODULE PROCEDURE &
            env_mp_bcast_i_0d, &
            env_mp_bcast_i_1d, &
            env_mp_bcast_i8_1d, &
            env_mp_bcast_i_2d, &
            env_mp_bcast_i_3d, &
            env_mp_bcast_i_4d, &
            env_mp_bcast_r_0d, &
            env_mp_bcast_r_1d, &
            env_mp_bcast_r_2d, &
            env_mp_bcast_r_3d, &
            env_mp_bcast_r_4d, &
            env_mp_bcast_r_5d, &
            env_mp_bcast_c_0d, &
            env_mp_bcast_c_1d, &
            env_mp_bcast_c_2d, &
            env_mp_bcast_c_3d, &
            env_mp_bcast_c_4d, &
            env_mp_bcast_c_5d, &
            env_mp_bcast_c_6d, &
            env_mp_bcast_z_0d, &
            env_mp_bcast_z_1d, &
            env_mp_bcast_l_0d, &
            env_mp_bcast_l_1d, &
            env_mp_bcast_l_2d
    END INTERFACE
    !
    INTERFACE env_mp_sum
        MODULE PROCEDURE &
            env_mp_sum_i_0d, &
            env_mp_sum_i_1d, &
            env_mp_sum_i8_1d, &
            env_mp_sum_i_2d, &
            env_mp_sum_i_3d, &
            env_mp_sum_i_4d, &
            env_mp_sum_i_5d, &
            env_mp_sum_r_0d, &
            env_mp_sum_r_1d, &
            env_mp_sum_r_2d, &
            env_mp_sum_r_2d_res, &
            env_mp_sum_r_3d, &
            env_mp_sum_r_4d, &
            env_mp_sum_r_5d, &
            env_mp_sum_r_6d, &
            env_mp_sum_c_0d, &
            env_mp_sum_c_1d, &
            env_mp_sum_c_2d, &
            env_mp_sum_c_2d_res, &
            env_mp_sum_c_3d, &
            env_mp_sum_c_4d, &
            env_mp_sum_c_5d, &
            env_mp_sum_c_6d
    END INTERFACE
    !
    !------------------------------------------------------------------------------------
    !
    PRIVATE
    !
    PUBLIC :: env_mp_rank, env_mp_size, env_mp_abort, env_mp_bcast, env_mp_sum
    !
    !------------------------------------------------------------------------------------
CONTAINS
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    INTEGER FUNCTION env_mp_rank(comm)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: comm
        !
        INTEGER :: ierr, taskid
        !
        !--------------------------------------------------------------------------------
        !
        ierr = 0
        taskid = 0
        !
#if defined(__MPI)
        CALL mpi_comm_rank(comm, taskid, ierr)
        !
        IF (ierr /= 0) CALL env_mp_stop(8001)
#endif
        !
        env_mp_rank = taskid
        !
        !--------------------------------------------------------------------------------
    END FUNCTION env_mp_rank
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    INTEGER FUNCTION env_mp_size(comm)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: comm
        !
        INTEGER :: ierr, numtask
        !
        !--------------------------------------------------------------------------------
        !
        ierr = 0
        numtask = 1
        !
#if defined(__MPI)
        CALL mpi_comm_size(comm, numtask, ierr)
        !
        IF (ierr /= 0) CALL env_mp_stop(8002)
#endif
        !
        env_mp_size = numtask
        !
        !--------------------------------------------------------------------------------
    END FUNCTION env_mp_size
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_mp_abort(errorcode, gid)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: errorcode, gid
        !
        INTEGER :: ierr
        !
        !--------------------------------------------------------------------------------
        !
#if defined(__MPI)
        CALL env_deallocate_mp_buffers()
        !
        CALL mpi_abort(gid, errorcode, ierr)
#endif
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_mp_abort
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_mp_stop(code)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: code
        !
        INTEGER :: ierr
        !
        !--------------------------------------------------------------------------------
        !
        WRITE (stdout, fmt='( "*** error in Message Passing (mp) module ***")')
        WRITE (stdout, fmt='( "*** error code: ",I5)') code
        !
#if defined(__MPI)
        CALL mpi_abort(MPI_COMM_WORLD, code, ierr)
        ! abort with extreme prejudice across the entire MPI set of tasks
#endif
        !
        STOP
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_mp_stop
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_mp_bcast_i_0d(msg, source, gid)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: source, gid
        !
#if defined(__CUDA)
        INTEGER, DEVICE :: msg
        INTEGER :: ierr
#else
        INTEGER :: msg
#endif
        !
#if defined(__MPI)
        !
        !--------------------------------------------------------------------------------
        !
#if defined(__CUDA)&&defined(__GPU_MPI)
        ierr = cudaDeviceSynchronize() ! this syncs __GPU_MPI case
#endif
        !
        CALL env_bcast_integer(msg, 1, source, gid)
        !
#if defined(__CUDA)
        ierr = cudaDeviceSynchronize() ! this syncs SERIAL, __MPI
#endif
        !
#endif
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_mp_bcast_i_0d
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_mp_bcast_i_1d(msg, source, gid)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: source, gid
        !
#if defined(__CUDA)
        INTEGER, DEVICE :: msg(:)
        INTEGER :: ierr
#else
        INTEGER :: msg(:)
#endif
        !
#if defined(__MPI)
        !
        !--------------------------------------------------------------------------------
        !
#if defined(__CUDA)&&defined(__GPU_MPI)
        ierr = cudaDeviceSynchronize() ! this syncs __GPU_MPI case
#endif
        !
        CALL env_bcast_integer(msg, SIZE(msg), source, gid)
        !
#if defined(__CUDA)
        ierr = cudaDeviceSynchronize() ! this syncs SERIAL, __MPI
#endif
        !
#endif
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_mp_bcast_i_1d
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_mp_bcast_i8_1d(msg, source, gid)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: source, gid
        !
#if defined(__CUDA)
        INTEGER(i8b), DEVICE :: msg(:)
        INTEGER :: ierr
#else
        INTEGER(i8b) :: msg(:)
#endif
        !
#if defined(__MPI)
        !
        !--------------------------------------------------------------------------------
        !
#if defined(__CUDA)&&defined(__GPU_MPI)
        ierr = cudaDeviceSynchronize() ! this syncs __GPU_MPI case
#endif
        !
        CALL env_bcast_integer8(msg, SIZE(msg), source, gid)
        !
#if defined(__CUDA)
        ierr = cudaDeviceSynchronize() ! this syncs SERIAL, __MPI
#endif
        !
#endif
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_mp_bcast_i8_1d
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_mp_bcast_i_2d(msg, source, gid)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: source, gid
        !
#if defined(__CUDA)
        INTEGER, DEVICE :: msg(:, :)
        INTEGER :: ierr
#else
        INTEGER :: msg(:, :)
#endif
        !
#if defined(__MPI)
        !
        !--------------------------------------------------------------------------------
        !
#if defined(__CUDA)&&defined(__GPU_MPI)
        ierr = cudaDeviceSynchronize() ! this syncs __GPU_MPI case
#endif
        !
        CALL env_bcast_integer(msg, SIZE(msg), source, gid)
        !
#if defined(__CUDA)
        ierr = cudaDeviceSynchronize() ! this syncs SERIAL, __MPI
#endif
        !
#endif
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_mp_bcast_i_2d
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_mp_bcast_i_3d(msg, source, gid)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: source, gid
        !
#if defined(__CUDA)
        INTEGER, DEVICE :: msg(:, :, :)
        INTEGER :: ierr
#else
        INTEGER :: msg(:, :, :)
#endif
        !
#if defined(__MPI)
        !
        !--------------------------------------------------------------------------------
        !
#if defined(__CUDA)&&defined(__GPU_MPI)
        ierr = cudaDeviceSynchronize() ! this syncs __GPU_MPI case
#endif
        !
        CALL env_bcast_integer(msg, SIZE(msg), source, gid)
        !
#if defined(__CUDA)
        ierr = cudaDeviceSynchronize() ! this syncs SERIAL, __MPI
#endif
        !
#endif
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_mp_bcast_i_3d
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_mp_bcast_i_4d(msg, source, gid)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: source, gid
        !
#if defined(__CUDA)
        INTEGER, DEVICE :: msg(:, :, :, :)
        INTEGER :: ierr
#else
        INTEGER :: msg(:, :, :, :)
#endif
        !
#if defined(__MPI)
        !
        !--------------------------------------------------------------------------------
        !
#if defined(__CUDA)&&defined(__GPU_MPI)
        ierr = cudaDeviceSynchronize() ! this syncs __GPU_MPI case
#endif
        !
        CALL env_bcast_integer(msg, SIZE(msg), source, gid)
        !
#if defined(__CUDA)
        ierr = cudaDeviceSynchronize() ! this syncs SERIAL, __MPI
#endif
        !
#endif
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_mp_bcast_i_4d
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_mp_bcast_r_0d(msg, source, gid)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: source, gid
        !
#if defined(__CUDA)
        REAL(DP), DEVICE :: msg
        INTEGER :: ierr
#else
        REAL(DP) :: msg
#endif
        !
#if defined(__MPI)
        !
        !--------------------------------------------------------------------------------
        !
#if defined(__CUDA)&&defined(__GPU_MPI)
        ierr = cudaDeviceSynchronize() ! this syncs __GPU_MPI case
#endif
        !
        CALL env_bcast_real(msg, 1, source, gid)
        !
#if defined(__CUDA)
        ierr = cudaDeviceSynchronize() ! this syncs SERIAL, __MPI
#endif
        !
#endif
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_mp_bcast_r_0d
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_mp_bcast_r_1d(msg, source, gid)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: source, gid
        !
#if defined(__CUDA)
        REAL(DP), DEVICE :: msg(:)
        INTEGER :: ierr
#else
        REAL(DP) :: msg(:)
#endif
        !
#if defined(__MPI)
        !
        !--------------------------------------------------------------------------------
        !
#if defined(__CUDA)&&defined(__GPU_MPI)
        ierr = cudaDeviceSynchronize() ! this syncs __GPU_MPI case
#endif
        !
        CALL env_bcast_real(msg, SIZE(msg), source, gid)
        !
#if defined(__CUDA)
        ierr = cudaDeviceSynchronize() ! this syncs SERIAL, __MPI
#endif
        !
#endif
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_mp_bcast_r_1d
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_mp_bcast_r_2d(msg, source, gid)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: source, gid
        !
#if defined(__CUDA)
        REAL(DP), DEVICE :: msg(:, :)
        INTEGER :: ierr
#else
        REAL(DP) :: msg(:, :)
#endif
        !
#if defined(__MPI)
        !
        !--------------------------------------------------------------------------------
        !
#if defined(__CUDA)&&defined(__GPU_MPI)
        ierr = cudaDeviceSynchronize() ! this syncs __GPU_MPI case
#endif
        !
        CALL env_bcast_real(msg, SIZE(msg), source, gid)
        !
#if defined(__CUDA)
        ierr = cudaDeviceSynchronize() ! this syncs SERIAL, __MPI
#endif
        !
#endif
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_mp_bcast_r_2d
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_mp_bcast_r_3d(msg, source, gid)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: source, gid
        !
#if defined(__CUDA)
        REAL(DP), DEVICE :: msg(:, :, :)
        INTEGER :: ierr
#else
        REAL(DP) :: msg(:, :, :)
#endif
        !
#if defined(__MPI)
        !
        !--------------------------------------------------------------------------------
        !
#if defined(__CUDA)&&defined(__GPU_MPI)
        ierr = cudaDeviceSynchronize() ! this syncs __GPU_MPI case
#endif
        !
        CALL env_bcast_real(msg, SIZE(msg), source, gid)
        !
#if defined(__CUDA)
        ierr = cudaDeviceSynchronize() ! this syncs SERIAL, __MPI
#endif
        !
#endif
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_mp_bcast_r_3d
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_mp_bcast_r_4d(msg, source, gid)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: source, gid
        !
#if defined(__CUDA)
        REAL(DP), DEVICE :: msg(:, :, :, :)
        INTEGER :: ierr
#else
        REAL(DP) :: msg(:, :, :, :)
#endif
        !
#if defined(__MPI)
        !
        !--------------------------------------------------------------------------------
        !
#if defined(__CUDA)&&defined(__GPU_MPI)
        ierr = cudaDeviceSynchronize() ! this syncs __GPU_MPI case
#endif
        !
        CALL env_bcast_real(msg, SIZE(msg), source, gid)
        !
#if defined(__CUDA)
        ierr = cudaDeviceSynchronize() ! this syncs SERIAL, __MPI
#endif
        !
#endif
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_mp_bcast_r_4d
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_mp_bcast_r_5d(msg, source, gid)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: source, gid
        !
#if defined(__CUDA)
        REAL(DP), DEVICE :: msg(:, :, :, :, :)
        INTEGER :: ierr
#else
        REAL(DP) :: msg(:, :, :, :, :)
#endif
        !
#if defined(__MPI)
        !
        !--------------------------------------------------------------------------------
        !
#if defined(__CUDA)&&defined(__GPU_MPI)
        ierr = cudaDeviceSynchronize() ! this syncs __GPU_MPI case
#endif
        !
        CALL env_bcast_real(msg, SIZE(msg), source, gid)
        !
#if defined(__CUDA)
        ierr = cudaDeviceSynchronize() ! this syncs SERIAL, __MPI
#endif
        !
#endif
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_mp_bcast_r_5d
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_mp_bcast_c_0d(msg, source, gid)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: source, gid
        !
#if defined(__CUDA)
        COMPLEX(DP), DEVICE :: msg
        INTEGER :: ierr
#else
        COMPLEX(DP) :: msg
#endif
        !
#if defined(__MPI)
        !
        !--------------------------------------------------------------------------------
        !
#if defined(__CUDA)&&defined(__GPU_MPI)
        ierr = cudaDeviceSynchronize() ! this syncs __GPU_MPI case
#endif
        !
        CALL env_bcast_complex(msg, 1, source, gid)
        !
#if defined(__CUDA)
        ierr = cudaDeviceSynchronize() ! this syncs SERIAL, __MPI
#endif
        !
#endif
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_mp_bcast_c_0d
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_mp_bcast_c_1d(msg, source, gid)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: source, gid
        !
#if defined(__CUDA)
        COMPLEX(DP), DEVICE :: msg(:)
        INTEGER :: ierr
#else
        COMPLEX(DP) :: msg(:)
#endif
        !
#if defined(__MPI)
        !
        !--------------------------------------------------------------------------------
        !
#if defined(__CUDA)&&defined(__GPU_MPI)
        ierr = cudaDeviceSynchronize() ! this syncs __GPU_MPI case
#endif
        !
        CALL env_bcast_complex(msg, SIZE(msg), source, gid)
        !
#if defined(__CUDA)
        ierr = cudaDeviceSynchronize() ! this syncs SERIAL, __MPI
#endif
        !
#endif
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_mp_bcast_c_1d
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_mp_bcast_c_2d(msg, source, gid)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: source, gid
        !
#if defined(__CUDA)
        COMPLEX(DP), DEVICE :: msg(:, :)
        INTEGER :: ierr
#else
        COMPLEX(DP) :: msg(:, :)
#endif
        !
#if defined(__MPI)
        !
        !--------------------------------------------------------------------------------
        !
#if defined(__CUDA)&&defined(__GPU_MPI)
        ierr = cudaDeviceSynchronize() ! this syncs __GPU_MPI case
#endif
        !
        CALL env_bcast_complex(msg, SIZE(msg), source, gid)
        !
#if defined(__CUDA)
        ierr = cudaDeviceSynchronize() ! this syncs SERIAL, __MPI
#endif
        !
#endif
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_mp_bcast_c_2d
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_mp_bcast_c_3d(msg, source, gid)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: source, gid
        !
#if defined(__CUDA)
        COMPLEX(DP), DEVICE :: msg(:, :, :)
        INTEGER :: ierr
#else
        COMPLEX(DP) :: msg(:, :, :)
#endif
        !
#if defined(__MPI)
        !
        !--------------------------------------------------------------------------------
        !
#if defined(__CUDA)&&defined(__GPU_MPI)
        ierr = cudaDeviceSynchronize() ! this syncs __GPU_MPI case
#endif
        !
        CALL env_bcast_complex(msg, SIZE(msg), source, gid)
        !
#if defined(__CUDA)
        ierr = cudaDeviceSynchronize() ! this syncs SERIAL, __MPI
#endif
        !
#endif
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_mp_bcast_c_3d
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_mp_bcast_c_4d(msg, source, gid)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: source, gid
        !
#if defined(__CUDA)
        COMPLEX(DP), DEVICE :: msg(:, :, :, :)
        INTEGER :: ierr
#else
        COMPLEX(DP) :: msg(:, :, :, :)
#endif
        !
#if defined(__MPI)
        !
        !--------------------------------------------------------------------------------
        !
#if defined(__CUDA)&&defined(__GPU_MPI)
        ierr = cudaDeviceSynchronize() ! this syncs __GPU_MPI case
#endif
        !
        CALL env_bcast_complex(msg, SIZE(msg), source, gid)
        !
#if defined(__CUDA)
        ierr = cudaDeviceSynchronize() ! this syncs SERIAL, __MPI
#endif
        !
#endif
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_mp_bcast_c_4d
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_mp_bcast_c_5d(msg, source, gid)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: source, gid
        !
#if defined(__CUDA)
        COMPLEX(DP), DEVICE :: msg(:, :, :, :, :)
        INTEGER :: ierr
#else
        COMPLEX(DP) :: msg(:, :, :, :, :)
#endif
        !
#if defined(__MPI)
        !
        !--------------------------------------------------------------------------------
        !
#if defined(__CUDA)&&defined(__GPU_MPI)
        ierr = cudaDeviceSynchronize() ! this syncs __GPU_MPI case
#endif
        !
        CALL env_bcast_complex(msg, SIZE(msg), source, gid)
        !
#if defined(__CUDA)
        ierr = cudaDeviceSynchronize() ! this syncs SERIAL, __MPI
#endif
        !
#endif
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_mp_bcast_c_5d
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_mp_bcast_c_6d(msg, source, gid)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: source, gid
        !
#if defined(__CUDA)
        COMPLEX(DP), DEVICE :: msg(:, :, :, :, :, :)
        INTEGER :: ierr
#else
        COMPLEX(DP) :: msg(:, :, :, :, :, :)
#endif
        !
#if defined(__MPI)
        !
        !--------------------------------------------------------------------------------
        !
#if defined(__CUDA)&&defined(__GPU_MPI)
        ierr = cudaDeviceSynchronize() ! this syncs __GPU_MPI case
#endif
        !
        CALL env_bcast_complex(msg, SIZE(msg), source, gid)
        !
#if defined(__CUDA)
        ierr = cudaDeviceSynchronize() ! this syncs SERIAL, __MPI
#endif
        !
#endif
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_mp_bcast_c_6d
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_mp_bcast_l_0d(msg, source, gid)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: source, gid
        !
#if defined(__CUDA)
        LOGICAL, DEVICE :: msg
        INTEGER :: ierr
#else
        LOGICAL :: msg
#endif
        !
#if defined(__MPI)
        !
        !--------------------------------------------------------------------------------
        !
#if defined(__CUDA)&&defined(__GPU_MPI)
        ierr = cudaDeviceSynchronize() ! this syncs __GPU_MPI case
#endif
        !
        CALL env_bcast_logical(msg, 1, source, gid)
        !
#if defined(__CUDA)
        ierr = cudaDeviceSynchronize() ! this syncs SERIAL, __MPI
#endif
        !
#endif
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_mp_bcast_l_0d
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_mp_bcast_l_1d(msg, source, gid)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: source, gid
        !
#if defined(__CUDA)
        LOGICAL, DEVICE :: msg(:)
        INTEGER :: ierr
#else
        LOGICAL :: msg(:)
#endif
        !
#if defined(__MPI)
        !
        !--------------------------------------------------------------------------------
        !
#if defined(__CUDA)&&defined(__GPU_MPI)
        ierr = cudaDeviceSynchronize() ! this syncs __GPU_MPI case
#endif
        !
        CALL env_bcast_logical(msg, SIZE(msg), source, gid)
        !
#if defined(__CUDA)
        ierr = cudaDeviceSynchronize() ! this syncs SERIAL, __MPI
#endif
        !
#endif
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_mp_bcast_l_1d
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_mp_bcast_l_2d(msg, source, gid)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: source, gid
        !
#if defined(__CUDA)
        LOGICAL, DEVICE :: msg(:, :)
        INTEGER :: ierr
#else
        LOGICAL :: msg(:, :)
#endif
        !
#if defined(__MPI)
        !
        !--------------------------------------------------------------------------------
        !
#if defined(__CUDA)&&defined(__GPU_MPI)
        ierr = cudaDeviceSynchronize() ! this syncs __GPU_MPI case
#endif
        !
        CALL env_bcast_logical(msg, SIZE(msg), source, gid)
        !
#if defined(__CUDA)
        ierr = cudaDeviceSynchronize() ! this syncs SERIAL, __MPI
#endif
        !
#endif
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_mp_bcast_l_2d
    !------------------------------------------------------------------------------------
    !>
    !! Converts a string into an array of integers for broadcasting
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_mp_bcast_z_0d(msg, source, gid)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: source, gid
        !
        CHARACTER(LEN=*) :: msg
        INTEGER, ALLOCATABLE :: imsg(:)
        !
        INTEGER :: msglen, ierr, i
        !
#if defined(__MPI)
        !
        !--------------------------------------------------------------------------------
        !
        ierr = 0
        msglen = LEN(msg)
        !
        ALLOCATE (imsg(1:msglen), STAT=ierr)
        !
        IF (ierr /= 0) CALL env_mp_stop(8003)
        !
        DO i = 1, msglen
            imsg(i) = ICHAR(msg(i:i))
        END DO
        !
        CALL env_bcast_integer(imsg, msglen, source, gid)
        !
        DO i = 1, msglen
            msg(i:i) = CHAR(imsg(i))
        END DO
        !
        DEALLOCATE (imsg, STAT=ierr)
        !
        IF (ierr /= 0) CALL env_mp_stop(8004)
        !
#endif
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_mp_bcast_z_0d
    !------------------------------------------------------------------------------------
    !>
    !! Converts an array of strings into an array of integer arrays for broadcasting
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_mp_bcast_z_1d(msg, source, gid)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: source, gid
        !
        CHARACTER(LEN=*) :: msg(:)
        INTEGER, ALLOCATABLE :: imsg(:, :)
        !
        INTEGER :: msglen, m1, m2, ierr, i, j
        !
#if defined(__MPI)
        !
        !--------------------------------------------------------------------------------
        !
        ierr = 0
        m1 = LEN(msg)
        m2 = SIZE(msg)
        msglen = LEN(msg) * SIZE(msg)
        !
        ALLOCATE (imsg(1:m1, 1:m2), STAT=ierr)
        !
        IF (ierr /= 0) CALL env_mp_stop(8005)
        !
        DO j = 1, m2
            !
            DO i = 1, m1
                imsg(i, j) = ICHAR(msg(j) (i:i))
            END DO
            !
        END DO
        !
        CALL env_bcast_integer(imsg, msglen, source, gid)
        !
        DO j = 1, m2
            !
            DO i = 1, m1
                msg(j) (i:i) = CHAR(imsg(i, j))
            END DO
            !
        END DO
        !
        DEALLOCATE (imsg, STAT=ierr)
        !
        IF (ierr /= 0) CALL env_mp_stop(8006)
        !
#endif
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_mp_bcast_z_1d
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_mp_sum_i_0d(msg, gid)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: gid
        !
#if defined(__CUDA)
        INTEGER, DEVICE :: msg
        INTEGER :: ierr
#else
        INTEGER :: msg
#endif
        !
#if defined(__MPI)
#if defined(__CUDA)
        !
        !--------------------------------------------------------------------------------
        ! Avoid unnecessary communications on __MPI and syncs SERIAL
        !
        IF (env_mp_size(gid) == 1) THEN
            ierr = cudaDeviceSynchronize()
            !
            RETURN
            !
        END IF
        !
#if defined(__GPU_MPI)
        ierr = cudaDeviceSynchronize() ! this syncs __GPU_MPI
#endif
#endif
        !
        CALL env_reduce_base_integer(1, msg, gid, -1)
        !
#if defined(__CUDA)
        ierr = cudaDeviceSynchronize() ! this syncs __MPI for small copies
#endif
        !
#endif
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_mp_sum_i_0d
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_mp_sum_i_1d(msg, gid)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: gid
        !
#if defined(__CUDA)
        INTEGER, DEVICE :: msg(:)
        INTEGER :: ierr
#else
        INTEGER :: msg(:)
#endif
        !
#if defined(__MPI)
#if defined(__CUDA)
        !
        !--------------------------------------------------------------------------------
        ! Avoid unnecessary communications on __MPI and syncs SERIAL
        !
        IF (env_mp_size(gid) == 1) THEN
            ierr = cudaDeviceSynchronize()
            !
            RETURN
            !
        END IF
        !
#if defined(__GPU_MPI)
        ierr = cudaDeviceSynchronize() ! this syncs __GPU_MPI
#endif
#endif
        !
        CALL env_reduce_base_integer(SIZE(msg), msg, gid, -1)
        !
#if defined(__CUDA)
        ierr = cudaDeviceSynchronize() ! this syncs __MPI for small copies
#endif
        !
#endif
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_mp_sum_i_1d
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_mp_sum_i8_1d(msg, gid)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: gid
        !
#if defined(__CUDA)
        INTEGER(i8b), DEVICE :: msg(:)
        INTEGER :: ierr
#else
        INTEGER(i8b) :: msg(:)
#endif
        !
#if defined(__MPI)
#if defined(__CUDA)
        !
        !--------------------------------------------------------------------------------
        ! Avoid unnecessary communications on __MPI and syncs SERIAL
        !
        IF (env_mp_size(gid) == 1) THEN
            ierr = cudaDeviceSynchronize()
            !
            RETURN
            !
        END IF
        !
#if defined(__GPU_MPI)
        ierr = cudaDeviceSynchronize() ! this syncs __GPU_MPI
#endif
#endif
        !
        CALL env_reduce_base_integer8(SIZE(msg), msg, gid, -1)
        !
#if defined(__CUDA)
        ierr = cudaDeviceSynchronize() ! this syncs __MPI for small copies
#endif
        !
#endif
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_mp_sum_i8_1d
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_mp_sum_i_2d(msg, gid)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: gid
        !
#if defined(__CUDA)
        INTEGER, DEVICE :: msg(:, :)
        INTEGER :: ierr
#else
        INTEGER :: msg(:, :)
#endif
        !
#if defined(__MPI)
#if defined(__CUDA)
        !
        !--------------------------------------------------------------------------------
        ! Avoid unnecessary communications on __MPI and syncs SERIAL
        !
        IF (env_mp_size(gid) == 1) THEN
            ierr = cudaDeviceSynchronize()
            !
            RETURN
            !
        END IF
        !
#if defined(__GPU_MPI)
        ierr = cudaDeviceSynchronize() ! this syncs __GPU_MPI
#endif
#endif
        !
        CALL env_reduce_base_integer(SIZE(msg), msg, gid, -1)
        !
#if defined(__CUDA)
        ierr = cudaDeviceSynchronize() ! this syncs __MPI for small copies
#endif
        !
#endif
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_mp_sum_i_2d
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_mp_sum_i_3d(msg, gid)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: gid
        !
#if defined(__CUDA)
        INTEGER, DEVICE :: msg(:, :, :)
        INTEGER :: ierr
#else
        INTEGER :: msg(:, :, :)
#endif
        !
#if defined(__MPI)
#if defined(__CUDA)
        !
        !--------------------------------------------------------------------------------
        ! Avoid unnecessary communications on __MPI and syncs SERIAL
        !
        IF (env_mp_size(gid) == 1) THEN
            ierr = cudaDeviceSynchronize()
            !
            RETURN
            !
        END IF
        !
#if defined(__GPU_MPI)
        ierr = cudaDeviceSynchronize() ! this syncs __GPU_MPI
#endif
#endif
        !
        CALL env_reduce_base_integer(SIZE(msg), msg, gid, -1)
        !
#if defined(__CUDA)
        ierr = cudaDeviceSynchronize() ! this syncs __MPI for small copies
#endif
        !
#endif
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_mp_sum_i_3d
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_mp_sum_i_4d(msg, gid)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: gid
        !
#if defined(__CUDA)
        INTEGER, DEVICE :: msg(:, :, :, :)
        INTEGER :: ierr
#else
        INTEGER :: msg(:, :, :, :)
#endif
        !
#if defined(__MPI)
#if defined(__CUDA)
        !
        !--------------------------------------------------------------------------------
        ! Avoid unnecessary communications on __MPI and syncs SERIAL
        !
        IF (env_mp_size(gid) == 1) THEN
            ierr = cudaDeviceSynchronize()
            !
            RETURN
            !
        END IF
        !
#if defined(__GPU_MPI)
        ierr = cudaDeviceSynchronize() ! this syncs __GPU_MPI
#endif
#endif
        !
        CALL env_reduce_base_integer(SIZE(msg), msg, gid, -1)
        !
#if defined(__CUDA)
        ierr = cudaDeviceSynchronize() ! this syncs __MPI for small copies
#endif
        !
#endif
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_mp_sum_i_4d
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_mp_sum_i_5d(msg, gid)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: gid
        !
#if defined(__CUDA)
        INTEGER, DEVICE :: msg(:, :, :, :, :)
        INTEGER :: ierr
#else
        INTEGER :: msg(:, :, :, :, :)
#endif
        !
#if defined(__MPI)
#if defined(__CUDA)
        !
        !--------------------------------------------------------------------------------
        ! Avoid unnecessary communications on __MPI and syncs SERIAL
        !
        IF (env_mp_size(gid) == 1) THEN
            ierr = cudaDeviceSynchronize()
            !
            RETURN
            !
        END IF
        !
#if defined(__GPU_MPI)
        ierr = cudaDeviceSynchronize() ! this syncs __GPU_MPI
#endif
#endif
        !
        CALL env_reduce_base_integer(SIZE(msg), msg, gid, -1)
        !
#if defined(__CUDA)
        ierr = cudaDeviceSynchronize() ! this syncs __MPI for small copies
#endif
        !
#endif
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_mp_sum_i_5d
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_mp_sum_r_0d(msg, gid)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: gid
        !
#if defined(__CUDA)
        REAL(DP), DEVICE :: msg
        INTEGER :: ierr
#else
        REAL(DP) :: msg
#endif
        !
#if defined(__MPI)
#if defined(__CUDA)
        !
        !--------------------------------------------------------------------------------
        ! Avoid unnecessary communications on __MPI and syncs SERIAL
        !
        IF (env_mp_size(gid) == 1) THEN
            ierr = cudaDeviceSynchronize()
            !
            RETURN
            !
        END IF
        !
#if defined(__GPU_MPI)
        ierr = cudaDeviceSynchronize() ! this syncs __GPU_MPI
#endif
#endif
        !
        CALL env_reduce_base_real(1, msg, gid, -1)
        !
#if defined(__CUDA)
        ierr = cudaDeviceSynchronize() ! this syncs __MPI for small copies
#endif
        !
#endif
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_mp_sum_r_0d
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_mp_sum_r_1d(msg, gid)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: gid
        !
#if defined(__CUDA)
        REAL(DP), DEVICE :: msg(:)
        INTEGER :: ierr
#else
        REAL(DP) :: msg(:)
#endif
        !
#if defined(__MPI)
#if defined(__CUDA)
        !
        !--------------------------------------------------------------------------------
        ! Avoid unnecessary communications on __MPI and syncs SERIAL
        !
        IF (env_mp_size(gid) == 1) THEN
            ierr = cudaDeviceSynchronize()
            !
            RETURN
            !
        END IF
        !
#if defined(__GPU_MPI)
        ierr = cudaDeviceSynchronize() ! this syncs __GPU_MPI
#endif
#endif
        !
        CALL env_reduce_base_real(SIZE(msg), msg, gid, -1)
        !
#if defined(__CUDA)
        ierr = cudaDeviceSynchronize() ! this syncs __MPI for small copies
#endif
        !
#endif
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_mp_sum_r_1d
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_mp_sum_r_2d(msg, gid)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: gid
        !
#if defined(__CUDA)
        REAL(DP), DEVICE :: msg(:, :)
        INTEGER :: ierr
#else
        REAL(DP) :: msg(:, :)
#endif
        !
#if defined(__MPI)
#if defined(__CUDA)
        !
        !--------------------------------------------------------------------------------
        ! Avoid unnecessary communications on __MPI and syncs SERIAL
        !
        IF (env_mp_size(gid) == 1) THEN
            ierr = cudaDeviceSynchronize()
            !
            RETURN
            !
        END IF
        !
#if defined(__GPU_MPI)
        ierr = cudaDeviceSynchronize() ! this syncs __GPU_MPI
#endif
#endif
        !
        CALL env_reduce_base_real(SIZE(msg), msg, gid, -1)
        !
#if defined(__CUDA)
        ierr = cudaDeviceSynchronize() ! this syncs __MPI for small copies
#endif
        !
#endif
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_mp_sum_r_2d
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_mp_sum_r_2d_res(msg, res, root, gid)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: root, gid
        !
#if defined(__CUDA)
        REAL(DP), INTENT(IN), DEVICE :: msg(:, :)
        !
        REAL(DP), INTENT(OUT), DEVICE :: res(:, :)
#else
        REAL(DP), INTENT(IN) :: msg(:, :)
        !
        REAL(DP), INTENT(OUT) :: res(:, :)
#endif
        !
        INTEGER :: msglen, taskid, ierr
        !
#if defined(__MPI)
        !
        !--------------------------------------------------------------------------------
        !
        msglen = SIZE(msg)
        !
        CALL mpi_comm_rank(gid, taskid, ierr)
        !
        IF (ierr /= 0) CALL env_mp_stop(8007)
        !
        IF (taskid == root .AND. msglen > SIZE(res)) CALL env_mp_stop(8008)
        !
#if defined(__CUDA)&&defined(__GPU_MPI)
        ierr = cudaDeviceSynchronize() ! this syncs __GPU_MPI
#endif
        !
        CALL env_reduce_base_real_to(msglen, msg, res, gid, root)
        !
#else
        res = msg
        !
#if defined(__CUDA)
        ierr = cudaDeviceSynchronize() ! this syncs SERIAL, __MPI
#endif
#endif
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_mp_sum_r_2d_res
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_mp_sum_r_3d(msg, gid)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: gid
        !
#if defined(__CUDA)
        REAL(DP), DEVICE :: msg(:, :, :)
        INTEGER :: ierr
#else
        REAL(DP) :: msg(:, :, :)
#endif
        !
#if defined(__MPI)
#if defined(__CUDA)
        !
        !--------------------------------------------------------------------------------
        ! Avoid unnecessary communications on __MPI and syncs SERIAL
        !
        IF (env_mp_size(gid) == 1) THEN
            ierr = cudaDeviceSynchronize()
            !
            RETURN
            !
        END IF
        !
#if defined(__GPU_MPI)
        ierr = cudaDeviceSynchronize() ! this syncs __GPU_MPI
#endif
#endif
        !
        CALL env_reduce_base_real(SIZE(msg), msg, gid, -1)
        !
#if defined(__CUDA)
        ierr = cudaDeviceSynchronize() ! this syncs __MPI for small copies
#endif
        !
#endif
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_mp_sum_r_3d
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_mp_sum_r_4d(msg, gid)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: gid
        !
#if defined(__CUDA)
        REAL(DP), DEVICE :: msg(:, :, :, :)
        INTEGER :: ierr
#else
        REAL(DP) :: msg(:, :, :, :)
#endif
        !
#if defined(__MPI)
#if defined(__CUDA)
        !
        !--------------------------------------------------------------------------------
        ! Avoid unnecessary communications on __MPI and syncs SERIAL
        !
        IF (env_mp_size(gid) == 1) THEN
            ierr = cudaDeviceSynchronize()
            !
            RETURN
            !
        END IF
        !
#if defined(__GPU_MPI)
        ierr = cudaDeviceSynchronize() ! this syncs __GPU_MPI
#endif
#endif
        !
        CALL env_reduce_base_real(SIZE(msg), msg, gid, -1)
        !
#if defined(__CUDA)
        ierr = cudaDeviceSynchronize() ! this syncs __MPI for small copies
#endif
        !
#endif
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_mp_sum_r_4d
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_mp_sum_r_5d(msg, gid)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: gid
        !
#if defined(__CUDA)
        REAL(DP), DEVICE :: msg(:, :, :, :, :)
        INTEGER :: ierr
#else
        REAL(DP) :: msg(:, :, :, :, :)
#endif
        !
#if defined(__MPI)
#if defined(__CUDA)
        !
        !--------------------------------------------------------------------------------
        ! Avoid unnecessary communications on __MPI and syncs SERIAL
        !
        IF (env_mp_size(gid) == 1) THEN
            ierr = cudaDeviceSynchronize()
            !
            RETURN
            !
        END IF
        !
#if defined(__GPU_MPI)
        ierr = cudaDeviceSynchronize() ! this syncs __GPU_MPI
#endif
#endif
        !
        CALL env_reduce_base_real(SIZE(msg), msg, gid, -1)
        !
#if defined(__CUDA)
        ierr = cudaDeviceSynchronize() ! this syncs __MPI for small copies
#endif
        !
#endif
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_mp_sum_r_5d
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_mp_sum_r_6d(msg, gid)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: gid
        !
#if defined(__CUDA)
        REAL(DP), DEVICE :: msg(:, :, :, :, :, :)
        INTEGER :: ierr
#else
        REAL(DP) :: msg(:, :, :, :, :, :)
#endif
        !
#if defined(__MPI)
#if defined(__CUDA)
        !
        !--------------------------------------------------------------------------------
        ! Avoid unnecessary communications on __MPI and syncs SERIAL
        !
        IF (env_mp_size(gid) == 1) THEN
            ierr = cudaDeviceSynchronize()
            !
            RETURN
            !
        END IF
        !
#if defined(__GPU_MPI)
        ierr = cudaDeviceSynchronize() ! this syncs __GPU_MPI
#endif
#endif
        !
        CALL env_reduce_base_real(SIZE(msg), msg, gid, -1)
        !
#if defined(__CUDA)
        ierr = cudaDeviceSynchronize() ! this syncs __MPI for small copies
#endif
        !
#endif
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_mp_sum_r_6d
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_mp_sum_c_0d(msg, gid)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: gid
        !
#if defined(__CUDA)
        COMPLEX(DP), DEVICE :: msg
        INTEGER :: ierr
#else
        COMPLEX(DP) :: msg
#endif
        !
#if defined(__MPI)
#if defined(__CUDA)
        !
        !--------------------------------------------------------------------------------
        ! Avoid unnecessary communications on __MPI and syncs SERIAL
        !
        IF (env_mp_size(gid) == 1) THEN
            ierr = cudaDeviceSynchronize()
            !
            RETURN
            !
        END IF
        !
#if defined(__GPU_MPI)
        ierr = cudaDeviceSynchronize() ! this syncs __GPU_MPI
#endif
#endif
        !
        CALL env_reduce_base_complex(1, msg, gid, -1)
        !
#if defined(__CUDA)
        ierr = cudaDeviceSynchronize() ! this syncs __MPI for small copies
#endif
        !
#endif
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_mp_sum_c_0d
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_mp_sum_c_1d(msg, gid)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: gid
        !
#if defined(__CUDA)
        COMPLEX(DP), DEVICE :: msg(:)
        INTEGER :: ierr
#else
        COMPLEX(DP) :: msg(:)
#endif
        !
#if defined(__MPI)
#if defined(__CUDA)
        !
        !--------------------------------------------------------------------------------
        ! Avoid unnecessary communications on __MPI and syncs SERIAL
        !
        IF (env_mp_size(gid) == 1) THEN
            ierr = cudaDeviceSynchronize()
            !
            RETURN
            !
        END IF
        !
#if defined(__GPU_MPI)
        ierr = cudaDeviceSynchronize() ! this syncs __GPU_MPI
#endif
#endif
        !
        CALL env_reduce_base_complex(SIZE(msg), msg, gid, -1)
        !
#if defined(__CUDA)
        ierr = cudaDeviceSynchronize() ! this syncs __MPI for small copies
#endif
        !
#endif
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_mp_sum_c_1d
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_mp_sum_c_2d(msg, gid)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: gid
        !
#if defined(__CUDA)
        COMPLEX(DP), DEVICE :: msg(:, :)
        INTEGER :: ierr
#else
        COMPLEX(DP) :: msg(:, :)
#endif
        !
#if defined(__MPI)
#if defined(__CUDA)
        !
        !--------------------------------------------------------------------------------
        ! Avoid unnecessary communications on __MPI and syncs SERIAL
        !
        IF (env_mp_size(gid) == 1) THEN
            ierr = cudaDeviceSynchronize()
            !
            RETURN
            !
        END IF
        !
#if defined(__GPU_MPI)
        ierr = cudaDeviceSynchronize() ! this syncs __GPU_MPI
#endif
#endif
        !
        CALL env_reduce_base_complex(SIZE(msg), msg, gid, -1)
        !
#if defined(__CUDA)
        ierr = cudaDeviceSynchronize() ! this syncs __MPI for small copies
#endif
        !
#endif
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_mp_sum_c_2d
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_mp_sum_c_2d_res(msg, res, gid)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: gid
        !
#if defined(__CUDA)
        COMPLEX(DP), INTENT(IN), DEVICE :: msg(:, :)
        !
        COMPLEX(DP), INTENT(OUT), DEVICE :: res(:, :)
        !
        INTEGER :: ierr
#else
        COMPLEX(DP), INTENT(IN) :: msg(:, :)
        !
        COMPLEX(DP), INTENT(OUT) :: res(:, :)
#endif
        !
#if defined(__MPI)
        !
        !--------------------------------------------------------------------------------
        !
#if defined(__CUDA)&&defined(__GPU_MPI)
        ierr = cudaDeviceSynchronize() ! this syncs __GPU_MPI
#endif
        !
        CALL env_reduce_base_complex_to(SIZE(msg), msg, res, gid, -1)
        !
#else
        res = msg
        !
#if defined(__CUDA)
        ierr = cudaDeviceSynchronize() ! this syncs SERIAL, __MPI
#endif
#endif
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_mp_sum_c_2d_res
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_mp_sum_c_3d(msg, gid)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: gid
        !
#if defined(__CUDA)
        COMPLEX(DP), DEVICE :: msg(:, :, :)
        INTEGER :: ierr
#else
        COMPLEX(DP) :: msg(:, :, :)
#endif
        !
#if defined(__MPI)
#if defined(__CUDA)
        !
        !--------------------------------------------------------------------------------
        ! Avoid unnecessary communications on __MPI and syncs SERIAL
        !
        IF (env_mp_size(gid) == 1) THEN
            ierr = cudaDeviceSynchronize()
            !
            RETURN
            !
        END IF
        !
#if defined(__GPU_MPI)
        ierr = cudaDeviceSynchronize() ! this syncs __GPU_MPI
#endif
#endif
        !
        CALL env_reduce_base_complex(SIZE(msg), msg, gid, -1)
        !
#if defined(__CUDA)
        ierr = cudaDeviceSynchronize() ! this syncs __MPI for small copies
#endif
        !
#endif
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_mp_sum_c_3d
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_mp_sum_c_4d(msg, gid)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: gid
        !
#if defined(__CUDA)
        COMPLEX(DP), DEVICE :: msg(:, :, :, :)
        INTEGER :: ierr
#else
        COMPLEX(DP) :: msg(:, :, :, :)
#endif
        !
#if defined(__MPI)
#if defined(__CUDA)
        !
        !--------------------------------------------------------------------------------
        ! Avoid unnecessary communications on __MPI and syncs SERIAL
        !
        IF (env_mp_size(gid) == 1) THEN
            ierr = cudaDeviceSynchronize()
            !
            RETURN
            !
        END IF
        !
#if defined(__GPU_MPI)
        ierr = cudaDeviceSynchronize() ! this syncs __GPU_MPI
#endif
#endif
        !
        CALL env_reduce_base_complex(SIZE(msg), msg, gid, -1)
        !
#if defined(__CUDA)
        ierr = cudaDeviceSynchronize() ! this syncs __MPI for small copies
#endif
        !
#endif
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_mp_sum_c_4d
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_mp_sum_c_5d(msg, gid)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: gid
        !
#if defined(__CUDA)
        COMPLEX(DP), DEVICE :: msg(:, :, :, :, :)
        INTEGER :: ierr
#else
        COMPLEX(DP) :: msg(:, :, :, :, :)
#endif
        !
#if defined(__MPI)
#if defined(__CUDA)
        !
        !--------------------------------------------------------------------------------
        ! Avoid unnecessary communications on __MPI and syncs SERIAL
        !
        IF (env_mp_size(gid) == 1) THEN
            ierr = cudaDeviceSynchronize()
            !
            RETURN
            !
        END IF
        !
#if defined(__GPU_MPI)
        ierr = cudaDeviceSynchronize() ! this syncs __GPU_MPI
#endif
#endif
        !
        CALL env_reduce_base_complex(SIZE(msg), msg, gid, -1)
        !
#if defined(__CUDA)
        ierr = cudaDeviceSynchronize() ! this syncs __MPI for small copies
#endif
        !
#endif
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_mp_sum_c_5d
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_mp_sum_c_6d(msg, gid)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: gid
        !
#if defined(__CUDA)
        COMPLEX(DP), DEVICE :: msg(:, :, :, :, :, :)
        INTEGER :: ierr
#else
        COMPLEX(DP) :: msg(:, :, :, :, :, :)
#endif
        !
#if defined(__MPI)
#if defined(__CUDA)
        !
        !--------------------------------------------------------------------------------
        ! Avoid unnecessary communications on __MPI and syncs SERIAL
        !
        IF (env_mp_size(gid) == 1) THEN
            ierr = cudaDeviceSynchronize()
            !
            RETURN
            !
        END IF
        !
#if defined(__GPU_MPI)
        ierr = cudaDeviceSynchronize() ! this syncs __GPU_MPI
#endif
#endif
        !
        CALL env_reduce_base_complex(SIZE(msg), msg, gid, -1)
        !
#if defined(__CUDA)
        ierr = cudaDeviceSynchronize() ! this syncs __MPI for small copies
#endif
        !
#endif
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_mp_sum_c_6d
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE env_mp
!----------------------------------------------------------------------------------------
