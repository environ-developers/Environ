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
! Authors: Carlo Cavazzoni and Stefano de Gironcoli, modified by PG
!
!----------------------------------------------------------------------------------------
!
#if defined(__CUDA)
!
#define __NON_BLOCKING_SCATTER
!
!----------------------------------------------------------------------------------------
!>
!! GPU version of env_fft_scatter.f90
!!
!----------------------------------------------------------------------------------------
MODULE env_fft_scatter_gpu
    !------------------------------------------------------------------------------------
    !
    USE cudafor
    !
    USE env_fft_param
    !
    USE env_types_fft, ONLY: env_fft_type_descriptor
    !
    USE env_fft_buffers, ONLY: env_check_fft_buffers_size, &
                               f_in => pin_space_scatter_in, &
                               f_aux => pin_space_scatter_out
    !
    !------------------------------------------------------------------------------------
    !
    PRIVATE
    !
    PUBLIC :: env_fft_scatter_xy_gpu, env_fft_scatter_yz_gpu, &
              env_fft_scatter_many_yz_gpu, env_fft_scatter_2d_gpu, &
              env_fft_scatter_many_columns_to_planes_store, &
              env_fft_scatter_many_columns_to_planes_send, &
              env_fft_scatter_many_planes_to_columns_store, &
              env_fft_scatter_many_planes_to_columns_send
    !
    !------------------------------------------------------------------------------------
CONTAINS
    !------------------------------------------------------------------------------------
    !>
    !! Transpose of the fft xy planes across the comm communicator.
    !! If the optional comm is not provided as input, the transpose is made
    !! across desc%comm2 communicator.
    !!
    !! a) From Y-oriented columns to X-oriented partial slices (isgn > 0)
    !!    Active columns along the Y direction corresponding to a subset of the
    !!    active X values and a range of Z values (in this order) are stored
    !!    consecutively for each processor and are such that the subgroup owns
    !!    all data for a range of Z values.
    !!
    !!    The Y pencil -> X-oriented partial slices transposition is performed
    !!    in the subgroup of processors (desc%comm2) owning this range of Z values.
    !!
    !!    The transpose takes place in two steps:
    !!    1) on each processor the columns are sliced into sections along Y
    !!       that are stored one after the other. On each processor, slices for
    !!       processor "iproc2" are desc%nr2p(iproc2)*desc%nr1p(me2)*desc%my_nr3p big.
    !!    2) all processors communicate to exchange slices (all sectin of columns with
    !!       Y in the slice belonging to "me" must be received, all the others
    !!       must be sent to "iproc2")
    !!
    !!    Finally one gets the "partial slice" representation: each processor has
    !!    all the X values of desc%my_nr2p Y and desc%my_nr3p Z values.
    !!    Data are organized with the X index running fastest, then Y, then Z.
    !!
    !!    f_in  contains the input Y columns, is destroyed on output
    !!    f_aux contains the output X-oriented partial slices.
    !!
    !! b) From planes to columns (isgn < 0)
    !!
    !!    Quite the same in the opposite direction
    !!    f_aux contains the input X-oriented partial slices, is destroyed on output
    !!    f_in  contains the output Y columns.
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_fft_scatter_xy_gpu(desc, f_in_d, f_aux_d, nxx_, isgn, stream)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(env_fft_type_descriptor), INTENT(IN) :: desc
        INTEGER, INTENT(IN) :: nxx_, isgn
        !
        INTEGER(kind=cuda_stream_kind), INTENT(IN) :: stream
        ! cuda stream for the execution
        !
        COMPLEX(DP), DEVICE, INTENT(INOUT) :: f_in_d(nxx_), f_aux_d(nxx_)
        !
#if defined(__MPI)
        !
        INTEGER :: ierr, me2, nproc2, iproc2, ncpx, nr1x, my_nr2p, nr2px, ip, ip0
        INTEGER :: i, it, j, k, kfrom, kdest, offset, ioff, mc, m1, m3, i1,
        INTEGER :: icompact, sendsize, aux
        INTEGER, ALLOCATABLE :: ncp_(:)
        INTEGER, POINTER, DEVICE :: nr1p__d(:), indx_d(:, :)
        !
#if defined(__NON_BLOCKING_SCATTER)
        INTEGER :: sh(desc%nproc2), rh(desc%nproc2)
#endif
        !
        CHARACTER(LEN=80) :: sub_name = 'env_fft_scatter_xy_gpu'
        !
        !--------------------------------------------------------------------------------
        !
        me2 = desc%mype2 + 1
        nproc2 = desc%nproc2
        !
        IF (ABS(isgn) == 3) nproc2 = 1
        !
        nr1x = desc%nr1x
        ! this mapping improves readability but, most importantly, it is needed
        ! IN the cuf kernels (as of PGI 18.10) since otherwise, when variables from
        ! `desc` appear in the body of the do loops, the compiler generates incorrect
        ! GPU code.
        !
        ALLOCATE (ncp_(nproc2)) ! allocate auxiliary array for columns distribution
        !
        ncp_ = desc%nr1p * desc%my_nr3p
        nr1p__d => desc%nr1p_d
        indx_d => desc%indp_d
        my_nr2p = desc%my_nr2p
        !
        CALL env_start_clock(sub_name)
        !
        !--------------------------------------------------------------------------------
        ! Calculate the message size
        !
        nr2px = MAXVAL(desc%nr2p) ! maximum number of Y values to be disributed
        !
        ncpx = MAXVAL(ncp_) ! maximum number of Y columns to be disributed
        !
        sendsize = ncpx * nr2px ! dimension of the scattered chunks (safe value)
        !
        CALL env_check_fft_buffers_size(desc)
        ! check host copy allocation of f_in and f_aux
        !
        ierr = 0
        !
        IF (isgn > 0) THEN
            !
            IF (nproc2 == 1) GO TO 10
            !
            !----------------------------------------------------------------------------
            ! "forward" scatter from columns to planes
            !
            ! step one: store contiguously the slices
            !
            offset = 0
            !
            DO iproc2 = 1, nproc2
                kdest = (iproc2 - 1) * sendsize
                kfrom = offset
                !
                ierr = cudaMemcpy2DAsync(f_aux(kdest + 1), nr2px, f_in_d(kfrom + 1), &
                                         desc%nr2x, desc%nr2p(iproc2), &
                                         ncp_(me2), cudaMemcpyDeviceToHost, stream)
                !
                offset = offset + desc%nr2p(iproc2)
            END DO
            !
            ! step two: communication across the nproc3 group
            !
#if defined(__NON_BLOCKING_SCATTER)
            ierr = cudaStreamSynchronize(stream)
            !
            DO iproc2 = 1, nproc2
                kdest = (iproc2 - 1) * sendsize
                !
                CALL mpi_irecv(f_in(kdest + 1), sendsize, &
                               MPI_DOUBLE_COMPLEX, iproc2 - 1, MPI_ANY_TAG, &
                               desc%comm2, rh(iproc2), ierr)
                !
                CALL mpi_isend(f_aux(kdest + 1), sendsize, &
                               MPI_DOUBLE_COMPLEX, iproc2 - 1, me2, desc%comm2, &
                               sh(iproc2), ierr)
                !
            END DO
#else
            !
            ierr = cudaStreamSynchronize(stream)
            !
            CALL mpi_alltoall(f_aux(1), sendsize, MPI_DOUBLE_COMPLEX, f_in(1), &
                              sendsize, MPI_DOUBLE_COMPLEX, desc%comm2, ierr)
            !
            IF (ABS(ierr) /= 0) CALL env_fft_error(sub_name, 'info<>0', ABS(ierr))
            !
            ierr = cudaMemcpyAsync(f_in_d, f_in, nxx_, cudaMemcpyHostToDevice, stream)
#endif
            !
10          CONTINUE
            !
            !$cuf kernel do (1) <<<*,*,0,stream>>>
            DO i = 1, nxx_
                f_aux_d(i) = (0.D0, 0.D0)
            END DO
            !
            DO iproc2 = 1, nproc2
#if defined(__NON_BLOCKING_SCATTER)
                !
                IF (nproc2 > 1) THEN
                    kdest = (iproc2 - 1) * sendsize
                    !
                    CALL mpi_wait(rh(iproc2), MPI_STATUSES_IGNORE, ierr)
                    !
                    CALL mpi_wait(sh(iproc2), MPI_STATUSES_IGNORE, ierr)
                    !
                    ierr = cudaMemcpyAsync(f_in_d(kdest + 1), f_in(kdest + 1), &
                                           sendsize, cudaMemcpyHostToDevice, stream)
                    !
                END IF
#endif
                !
                aux = ncp_(iproc2)
                !
                !$cuf kernel do(2) <<<*,*,0,stream>>>
                DO i = 1, aux
                    !
                    DO j = 1, my_nr2p
                        it = (iproc2 - 1) * sendsize + (i - 1) * nr2px
                        m3 = (i - 1) / nr1p__d(iproc2) + 1
                        i1 = MOD(i - 1, nr1p__d(iproc2)) + 1
                        m1 = indx_d(i1, iproc2)
                        icompact = m1 + (m3 - 1) * nr1x * my_nr2p + (j - 1) * nr1x
                        f_aux_d(icompact) = f_in_d(j + it)
                    END DO
                    !
                END DO
                !
            END DO
            !
        ELSE
            !
            !----------------------------------------------------------------------------
            ! "backward" scatter from planes to columns
            !
            DO iproc2 = 1, nproc2
                aux = ncp_(iproc2)
                !
                !$cuf kernel do (2) <<<*,*,0,stream>>>
                DO i = 1, aux
                    !
                    DO j = 1, my_nr2p
                        it = (iproc2 - 1) * sendsize + (i - 1) * nr2px
                        m3 = (i - 1) / nr1p__d(iproc2) + 1
                        i1 = MOD(i - 1, nr1p__d(iproc2)) + 1
                        m1 = indx_d(i1, iproc2)
                        icompact = m1 + (m3 - 1) * nr1x * my_nr2p + (j - 1) * nr1x
                        f_in_d(j + it) = f_aux_d(icompact)
                    END DO
                    !
                END DO
                !
#if defined(__NON_BLOCKING_SCATTER)
                !
                IF (nproc2 > 1) THEN
                    kdest = (iproc2 - 1) * sendsize
                    !
                    ierr = cudaMemcpyAsync(f_in(kdest + 1), f_in_d(kdest + 1), &
                                           sendsize, cudaMemcpyDeviceToHost, stream)
                    !
                END IF
#endif
                !
            END DO
            !
            IF (nproc2 == 1) GO TO 20
            !
            !  step two: communication
            !
#if defined(__NON_BLOCKING_SCATTER)
            DO iproc2 = 1, nproc2
                kdest = (iproc2 - 1) * sendsize
                !
                CALL mpi_irecv(f_aux((iproc2 - 1) * sendsize + 1), sendsize, &
                               MPI_DOUBLE_COMPLEX, iproc2 - 1, MPI_ANY_TAG, &
                               desc%comm2, rh(iproc2), ierr)
                !
                ierr = cudaStreamSynchronize(stream)
                !
                CALL mpi_isend(f_in(kdest + 1), &
                               sendsize, MPI_DOUBLE_COMPLEX, iproc2 - 1, me2, &
                               desc%comm2, sh(iproc2), ierr)
                !
            END DO
#else
            !
            ierr = cudaMemcpy(f_in, f_in_d, nxx_, cudaMemcpyDeviceToHost)
            !
            CALL mpi_alltoall(f_in(1), sendsize, MPI_DOUBLE_COMPLEX, f_aux(1), &
                              sendsize, MPI_DOUBLE_COMPLEX, desc%comm2, ierr)
            !
            IF (ABS(ierr) /= 0) CALL env_fft_error(sub_name, 'info<>0', ABS(ierr))
#endif
            !
            ! step one: store contiguously the columns
            !
            !----------------------------------------------------------------------------
            ! Ensures that no garbage is present in the output
            ! not useless ... clean the array to be returned from
            ! the garbage of previous A2A step
            !
            !$cuf kernel do (1) <<<*,*,0,stream>>>
            DO i = 1, nxx_
                f_in_d(i) = (0.D0, 0.D0)
            END DO
            !
            offset = 0
            !
            DO iproc2 = 1, nproc2
                kdest = (iproc2 - 1) * sendsize
                kfrom = offset
                !
#if defined(__NON_BLOCKING_SCATTER)
                CALL mpi_wait(rh(iproc2), MPI_STATUSES_IGNORE, ierr)
                !
                CALL mpi_wait(sh(iproc2), MPI_STATUSES_IGNORE, ierr)
#endif
                !
                ierr = cudaMemcpy2DAsync(f_in_d(kfrom + 1), &
                                         desc%nr2x, f_aux(kdest + 1), nr2px, &
                                         desc%nr2p(iproc2), ncp_(me2), &
                                         cudaMemcpyHostToDevice, stream)
                !
                offset = offset + desc%nr2p(iproc2)
            END DO
            !
20          CONTINUE
            !
        END IF
        !
        DEALLOCATE (ncp_)
        !
        CALL env_stop_clock(sub_name)
        !
#endif
        !
        RETURN
        !
99      FORMAT(20('(', 2F12.9, ')'))
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_fft_scatter_xy_gpu
    !------------------------------------------------------------------------------------
    !>
    !! Transpose of the fft yz planes across the desc%comm3 communicator
    !!
    !! a) From Z-oriented columns to Y-oriented colums (isgn > 0)
    !!    Active columns (or sticks or pencils) along the Z direction for each
    !!    processor are stored consecutively and are such that they correspond
    !!    to a subset of the active X values.
    !!
    !!    The pencil -> slices transposition is performed in the subgroup
    !!    of processors (desc%comm3) owning these X values.
    !!
    !!    The transpose takes place in two steps:
    !!    1) on each processor the columns are sliced into sections along Z
    !!       that are stored one after the other. On each processor, slices for
    !!       processor "iproc3" are desc%nr3p(iproc3)*desc%nsp/nsp(me) big.
    !!    2) all processors communicate to exchange slices (all columns with
    !!       Z in the slice belonging to "me" must be received, all the others
    !!       must be sent to "iproc3")
    !!
    !!    Finally one gets the "slice" representation: each processor has
    !!    desc%nr3p(mype3) Z values of all the active pencils along Y for the
    !!    X values of the current group. Data are organized with the Y index
    !!    running fastest, then the reordered X values, then Z.
    !!
    !!    f_in  contains the input Z columns, is destroyed on output
    !!    f_aux contains the output Y colums.
    !!
    !! b) From planes to columns (isgn < 0)
    !!
    !!    Quite the same in the opposite direction
    !!    f_aux contains the input Y columns, is destroyed on output
    !!    f_in  contains the output Z columns.
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_fft_scatter_yz_gpu(desc, f_in_d, f_aux_d, nxx_, isgn)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(env_fft_type_descriptor), INTENT(IN) :: desc
        INTEGER, INTENT(IN) :: nxx_, isgn
        !
        COMPLEX(DP), DEVICE, INTENT(INOUT) :: f_in_d(nxx_), f_aux_d(nxx_)
        !
#if defined(__MPI)
        INTEGER, DEVICE, POINTER :: desc_ismap_d(:)
        !
        INTEGER :: ierr, me, me2, me2_start, me2_end, me3
        INTEGER :: nproc3, iproc3, ncpx, nr3px, ip, ip0
        INTEGER :: nr1x, nr2x, nr3, nr3x
        INTEGER :: i, it, it0, k, kfrom, kdest, offset, ioff
        INTEGER :: mc, m1, m2, i1, sendsize, aux
        INTEGER, ALLOCATABLE :: ncp_(:)
        INTEGER, DEVICE, POINTER :: ir1p__d(:)
        INTEGER :: my_nr1p_
        !
#if defined(__NON_BLOCKING_SCATTER)
        INTEGER :: sh(desc%nproc3), rh(desc%nproc3)
#endif
        TYPE(cudaEvent) :: zero_event
        !
        CHARACTER(LEN=80) :: sub_name = 'env_fft_scatter_yz_gpu'
        !
        !--------------------------------------------------------------------------------
        !
        ierr = cudaEventCreate(zero_event)
        !
        me = desc%mype + 1
        me2 = desc%mype2 + 1
        me3 = desc%mype3 + 1
        nproc3 = desc%nproc3
        !
        nr1x = desc%nr1x
        nr2x = desc%nr2x
        nr3 = desc%nr3
        nr3x = desc%nr3x
        ! this mapping improves readability but, most importantly, it is needed
        ! in the cuf kernels (as of PGI 18.10) since otherwise, when variables from
        ! `desc` appear in the body of the do loops, the compiler generates incorrect
        ! GPU code.
        !
        ALLOCATE (ncp_(desc%nproc)) ! allocate auxiliary array for columns distribution

        me2_start = me2
        me2_end = me2
        !
        ncp_ = desc%nsp
        my_nr1p_ = COUNT(desc%ir1p > 0)
        ir1p__d => desc%ir1p_d
        !
        CALL env_start_clock(sub_name)
        !
        !--------------------------------------------------------------------------------
        ! Calculate the message size
        !
        nr3px = MAXVAL(desc%nr3p) ! maximum number of Z values to be exchanged
        ncpx = MAXVAL(ncp_) ! maximum number of Z columns to be exchanged
        !
        sendsize = ncpx * nr3px ! dimension of the scattered chunks
        !
        CALL env_check_fft_buffers_size(desc)
        ! check host copy allocation of f_in and f_aux
        !
        desc_ismap_d => desc%ismap_d
        !
        ierr = 0
        !
        IF (isgn > 0) THEN
            !
            IF (nproc3 == 1) GO TO 10
            !
            !----------------------------------------------------------------------------
            ! "forward" scatter from columns to planes
            !
            ! step one: store contiguously the slices
            !
            offset = 0
            !
            DO iproc3 = 1, nproc3
                kdest = (iproc3 - 1) * sendsize
                kfrom = offset
                !
                DO me2 = me2_start, me2_end
                    ip = desc%iproc(me2, me3)
                    !
                    ierr = cudaMemcpy2DAsync(f_aux(kdest + 1), nr3px, &
                                             f_in_d(kfrom + 1), nr3x, &
                                             desc%nr3p(iproc3), ncp_(ip), &
                                             cudaMemcpyDeviceToHost, &
                                             desc%stream_scatter_yz(iproc3))
                    !
                    kdest = kdest + nr3px * ncp_(ip)
                    kfrom = kfrom + nr3x * ncp_(ip)
                END DO
                !
                offset = offset + desc%nr3p(iproc3)
            END DO
            !
            !----------------------------------------------------------------------------
            ! Ensures that no garbage is present in the output
            ! useless; the later accessed elements are overwritten by the A2A step
            !
            ! step two: communication across the nproc3 group
            !
#if defined(__NON_BLOCKING_SCATTER)
            DO iproc3 = 1, nproc3
                ierr = cudaStreamSynchronize(desc%stream_scatter_yz(iproc3))
                !
                CALL mpi_irecv(f_in((iproc3 - 1) * sendsize + 1), sendsize, &
                               MPI_DOUBLE_COMPLEX, iproc3 - 1, MPI_ANY_TAG, &
                               desc%comm3, rh(iproc3), ierr)
                !
                CALL mpi_isend(f_aux((iproc3 - 1) * sendsize + 1), sendsize, &
                               MPI_DOUBLE_COMPLEX, iproc3 - 1, me3, desc%comm3, &
                               sh(iproc3), ierr)
                !
            END DO
#else
            !
            ierr = cudaDeviceSynchronize()
            !
            CALL mpi_alltoall(f_aux(1), sendsize, MPI_DOUBLE_COMPLEX, f_in(1), &
                              sendsize, MPI_DOUBLE_COMPLEX, desc%comm3, ierr)
            !
            IF (ABS(ierr) /= 0) CALL env_fft_error(sub_name, 'info<>0', ABS(ierr))
            !
            ierr = cudaMemcpy(f_in_d, f_in, nxx_, cudaMemcpyHostToDevice)
#endif
            !
10          CONTINUE
            !
            !$cuf kernel do (1) <<<*,*,0,desc%stream_scatter_yz(1)>>>
            DO i = 1, desc%my_nr3p * my_nr1p_ * nr2x
                f_aux_d(i) = (0.D0, 0.D0)
            END DO
            !
            ierr = cudaEventRecord(zero_event, desc%stream_scatter_yz(1))
            !
            DO iproc3 = 1, desc%nproc3
                it0 = (iproc3 - 1) * sendsize
                !
#if defined(__NON_BLOCKING_SCATTER)
                IF (nproc3 > 1) THEN
                    !
                    CALL mpi_wait(rh(iproc3), MPI_STATUSES_IGNORE, ierr)
                    !
                    CALL mpi_wait(sh(iproc3), MPI_STATUSES_IGNORE, ierr)
                    !
                    ierr = cudaMemcpyAsync(f_in_d(it0 + 1), f_in(it0 + 1), sendsize, &
                                           cudaMemcpyHostToDevice, &
                                           desc%stream_scatter_yz(iproc3))
                    !
                END IF
#endif
                !
                IF (iproc3 == 2) ierr = cudaEventSynchronize(zero_event)
                !
                DO me2 = me2_start, me2_end
                    ip = desc%iproc(me2, iproc3)
                    ioff = desc%iss(ip)
                    aux = ncp_(ip)
                    !
                    !$cuf kernel do(2) <<<*,*,0,desc%stream_scatter_yz(iproc3)>>>
                    DO i = 1, aux ! was ncp_(iproc3)
                        !
                        DO k = 1, desc%my_nr3p
                            it = it0 + (i - 1) * nr3px
                            !
                            mc = desc_ismap_d(i + ioff)
                            ! this is m1+(m2-1)*nr1x of the current pencil
                            !
                            m1 = MOD(mc - 1, nr1x) + 1
                            m2 = (mc - 1) / nr1x + 1
                            !
                            i1 = m2 + (ir1p__d(m1) - 1) * nr2x + &
                                 (k - 1) * nr2x * my_nr1p_
                            !
                            f_aux_d(i1) = f_in_d(k + it)
                            !
                        END DO
                        !
                    END DO
                    !
                    it0 = it0 + ncp_(ip) * nr3px
                END DO
                !
            END DO
            !
        ELSE
            !
            !----------------------------------------------------------------------------
            ! "backward" scatter from planes to columns
            !
            DO iproc3 = 1, nproc3
                it0 = (iproc3 - 1) * sendsize
                !
                DO me2 = me2_start, me2_end
                    ip = desc%iproc(me2, iproc3)
                    ioff = desc%iss(ip)
                    aux = ncp_(ip)
                    !
                    !$cuf kernel do(2) <<<*,*,0,desc%stream_scatter_yz(iproc3)>>>
                    DO i = 1, aux
                        !
                        DO k = 1, desc%my_nr3p
                            it = it0 + (i - 1) * nr3px
                            !
                            mc = desc_ismap_d(i + ioff)
                            ! this is m1+(m2-1)*nr1x of the current pencil
                            !
                            m1 = MOD(mc - 1, nr1x) + 1
                            m2 = (mc - 1) / nr1x + 1
                            !
                            i1 = m2 + (ir1p__d(m1) - 1) * nr2x + &
                                 (k - 1) * (nr2x * my_nr1p_)
                            !
                            f_in_d(k + it) = f_aux_d(i1)
                            !
                        END DO
                        !
                    END DO
                    !
                    it0 = it0 + ncp_(ip) * nr3px
                END DO
                !
#if defined(__NON_BLOCKING_SCATTER)
                IF (nproc3 > 1) THEN
                    kdest = (iproc3 - 1) * sendsize
                    !
                    ierr = cudaMemcpyAsync(f_in(kdest + 1), f_in_d(kdest + 1), &
                                           sendsize, cudaMemcpyDeviceToHost, &
                                           desc%stream_scatter_yz(iproc3))
                    !
                END IF
#endif
                !
            END DO
            !
            IF (nproc3 == 1) GO TO 20
            !
            ! step two: communication
            !
#if defined(__NON_BLOCKING_SCATTER)
            DO iproc3 = 1, nproc3
                !
                CALL mpi_irecv(f_aux((iproc3 - 1) * sendsize + 1), sendsize, &
                               MPI_DOUBLE_COMPLEX, iproc3 - 1, MPI_ANY_TAG, &
                               desc%comm3, rh(iproc3), ierr)
                !
                ierr = cudaStreamSynchronize(desc%stream_scatter_yz(iproc3))
                !
                CALL mpi_isend(f_in((iproc3 - 1) * sendsize + 1), &
                               sendsize, MPI_DOUBLE_COMPLEX, iproc3 - 1, &
                               me3, desc%comm3, sh(iproc3), ierr)
                !
            END DO
#else
            !
            ierr = cudaMemcpy(f_in, f_in_d, nxx_, cudaMemcpyDeviceToHost)
            !
            CALL mpi_alltoall(f_in(1), sendsize, MPI_DOUBLE_COMPLEX, f_aux(1), &
                              sendsize, MPI_DOUBLE_COMPLEX, desc%comm3, ierr)
            !
            IF (ABS(ierr) /= 0) CALL env_fft_error(sub_name, 'info<>0', ABS(ierr))
            !
#endif
            !
            !  step one: store contiguously the columns
            !
            offset = 0
            !
            DO iproc3 = 1, nproc3
                kdest = (iproc3 - 1) * sendsize
                !
#if defined(__NON_BLOCKING_SCATTER)
                IF (nproc3 > 1) THEN
                    !
                    CALL mpi_wait(rh(iproc3), MPI_STATUSES_IGNORE, ierr)
                    !
                    CALL mpi_wait(sh(iproc3), MPI_STATUSES_IGNORE, ierr)
                    !
                END IF
#endif
                !
                kfrom = offset
                !
                DO me2 = me2_start, me2_end
                    ip = desc%iproc(me2, me3)
                    !
                    ierr = cudaMemcpy2DAsync(f_in_d(kfrom + 1), nr3x, &
                                             f_aux(kdest + 1), nr3px, &
                                             desc%nr3p(iproc3), ncp_(ip), &
                                             cudaMemcpyHostToDevice, &
                                             desc%stream_scatter_yz(iproc3))
                    !
                    kdest = kdest + nr3px * ncp_(ip)
                    kfrom = kfrom + nr3x * ncp_(ip)
                END DO
                !
                offset = offset + desc%nr3p(iproc3)
            END DO
            !
            DO iproc3 = 1, nproc3
                ierr = cudaStreamSynchronize(desc%stream_scatter_yz(iproc3))
            END DO
            !
            !----------------------------------------------------------------------------
            ! Clean extra array elements in each stick
            !
            IF (nr3x /= nr3) THEN
                aux = ncp_(desc%mype + 1)
                !
                !$cuf kernel do(2) <<<*,*>>>
                DO k = 1, aux
                    !
                    DO i = nr3, nr3x
                        f_in_d((k - 1) * nr3x + i) = (0.D0, 0.D0)
                    END DO
                    !
                END DO
                !
            END IF
            !
20          CONTINUE
            !
        END IF
        !
        DEALLOCATE (ncp_)
        !
        CALL env_stop_clock(sub_name)
#endif
        !
        RETURN
        !
98      FORMAT(10('(', 2F12.9, ')'))
99      FORMAT(20('(', 2F12.9, ')'))
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_fft_scatter_yz_gpu
    !------------------------------------------------------------------------------------
    !>
    !! Transpose of the fft (many) yz planes across the desc%comm3 communicator
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_fft_scatter_many_yz_gpu(desc, f_in_d, f_aux_d, nxx_, isgn, howmany)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(env_fft_type_descriptor), INTENT(IN) :: desc
        INTEGER, INTENT(IN) :: nxx_, isgn, howmany
        !
        COMPLEX(DP), DEVICE, INTENT(INOUT) :: f_in_d(nxx_), f_aux_d(nxx_)

#if defined(__MPI)
        INTEGER, DEVICE, POINTER :: desc_ismap_d(:)
        !
        INTEGER :: ierr, me, me2, me3, nproc3, iproc3, ncpx, nr3px
        INTEGER :: ip, ip0, me2_start, me2_end
        INTEGER :: nr1x, nr2x, nr3, nr3x, nnr
        INTEGER :: i, j, it, it0, k, kfrom, kdest, offset, ioff
        INTEGER :: mc, m1, m2, i1, sendsize, aux
        INTEGER, ALLOCATABLE :: ncp_(:)
        INTEGER, DEVICE, POINTER :: ir1p__d(:)
        INTEGER :: my_nr1p_
        !
#if defined(__NON_BLOCKING_SCATTER)
        INTEGER :: sh(desc%nproc3), rh(desc%nproc3)
#endif
        !
        CHARACTER(LEN=80) :: sub_name = 'env_fft_scatter_many_yz_gpu'
        !
        !--------------------------------------------------------------------------------
        !
        me = desc%mype + 1
        me2 = desc%mype2 + 1
        me3 = desc%mype3 + 1
        nproc3 = desc%nproc3
        !
        nr1x = desc%nr1x
        nr2x = desc%nr2x
        nr3 = desc%nr3
        nr3x = desc%nr3x
        nnr = desc%nnr
        ! this mapping improves readability but, most importantly, it is needed
        ! in the cuf kernels (as of PGI 18.10) since otherwise, when variables from
        ! `desc` appear in the body of the do loops, the compiler generates incorrect
        ! GPU code.
        !
        ALLOCATE (ncp_(desc%nproc)) ! allocate auxiliary array for columns distribution
        !
        me2_start = me2
        me2_end = me2
        !
        ncp_ = desc%nsp
        my_nr1p_ = COUNT(desc%ir1p > 0)
        ir1p__d => desc%ir1p_d
        !
        CALL env_start_clock(sub_name)
        !
        !--------------------------------------------------------------------------------
        ! Calculate the message size
        !
        nr3px = MAXVAL(desc%nr3p) ! maximum number of Z values to be exchanged
        ncpx = MAXVAL(ncp_) ! maximum number of Z columns to be exchanged
        !
        sendsize = howmany * ncpx * nr3px ! dimension of the scattered chunks
        !
        CALL env_check_fft_buffers_size(desc, howmany)
        ! check dimensions of f_in and f_aux
        !
        desc_ismap_d => desc%ismap_d
        !
        ierr = 0
        IF (isgn > 0) THEN
            !
            IF (nproc3 == 1) GO TO 10
            !
            !----------------------------------------------------------------------------
            ! "forward" scatter from columns to planes
            !
            ! step one: store contiguously the slices
            !
            offset = 0
            !
            DO iproc3 = 1, nproc3
                kdest = (iproc3 - 1) * sendsize
                kfrom = offset
                !
                ierr = cudaMemcpy2D(f_aux(kdest + 1), nr3px, f_in_d(kfrom + 1), nr3x, &
                                    desc%nr3p(iproc3), howmany * ncpx, &
                                    cudaMemcpyDeviceToHost)
                !
                offset = offset + desc%nr3p(iproc3)
            END DO
            !
            !----------------------------------------------------------------------------
            ! Ensures that no garbage is present in the output
            ! useless; the later accessed elements are overwritten by the A2A step
            !
            ! step two: communication  across the    nproc3    group
            !
            CALL mpi_alltoall(f_aux(1), sendsize, MPI_DOUBLE_COMPLEX, f_in(1), &
                              sendsize, MPI_DOUBLE_COMPLEX, desc%comm3, ierr)
            !
            IF (ABS(ierr) /= 0) CALL env_fft_error(sub_name, 'info<>0', ABS(ierr))
            !
            ierr = cudaMemcpy(f_in_d, f_in, nxx_, cudaMemcpyHostToDevice)
            !
10          CONTINUE
            !
            !$cuf kernel do (1) <<<*,*>>>
            DO i = 1, nxx_
                f_aux_d(i) = (0.D0, 0.D0)
            END DO
            !
            DO iproc3 = 1, desc%nproc3
                it0 = (iproc3 - 1) * sendsize
                ioff = desc%iss(iproc3)
                aux = ncp_(iproc3)
                !
                !$cuf kernel do(3) <<<*,*>>>
                DO j = 0, howmany - 1
                    !
                    DO i = 1, aux ! was ncp_(iproc3)
                        !
                        DO k = 1, desc%my_nr3p
                            !
                            it = it0 + (i - 1) * nr3px + j * ncpx * nr3px
                            !desc%nnr !aux*nr3px
                            !
                            mc = desc_ismap_d(i + ioff)
                            ! this is m1+(m2-1)*nr1x of the current pencil
                            !
                            m1 = MOD(mc - 1, nr1x) + 1; m2 = (mc - 1) / nr1x + 1
                            !
                            i1 = m2 + (ir1p__d(m1) - 1) * nr2x + &
                                 (k - 1) * nr2x * my_nr1p_
                            !
                            f_aux_d(i1 + j * nnr) = f_in_d(k + it)
                            !
                        END DO
                        !
                    END DO
                    !
                END DO
                !
            END DO
            !
        ELSE
            !
            !----------------------------------------------------------------------------
            ! "backward" scatter from planes to columns
            !
            DO iproc3 = 1, nproc3
                it0 = (iproc3 - 1) * sendsize
                ip = desc%iproc(me2, iproc3)
                ioff = desc%iss(ip)
                aux = ncp_(ip)
                !
                !$cuf kernel do(3) <<<*,*>>>
                DO j = 0, howmany - 1
                    !
                    DO i = 1, aux
                        !
                        DO k = 1, desc%my_nr3p
                            it = it0 + (i - 1) * nr3px + j * ncpx * nr3px
                            !
                            mc = desc_ismap_d(i + ioff)
                            ! this is  m1+(m2-1)*nr1x  of the  current pencil
                            !
                            m1 = MOD(mc - 1, nr1x) + 1; m2 = (mc - 1) / nr1x + 1
                            !
                            i1 = m2 + (ir1p__d(m1) - 1) * nr2x + &
                                 (k - 1) * (nr2x * my_nr1p_)
                            !
                            f_in_d(k + it) = f_aux_d(i1 + j * nnr)
                        END DO
                        !
                    END DO
                    !
                END DO
                !
            END DO
            !
            IF (nproc3 == 1) GO TO 20
            !
            ! step two: communication
            !
            ierr = cudaMemcpy(f_in, f_in_d, nxx_, cudaMemcpyDeviceToHost)
            !
            CALL mpi_alltoall(f_in(1), sendsize, MPI_DOUBLE_COMPLEX, f_aux(1), &
                              sendsize, MPI_DOUBLE_COMPLEX, desc%comm3, ierr)
            !
            IF (ABS(ierr) /= 0) CALL env_fft_error(sub_name, 'info<>0', ABS(ierr))
            !
            ! step one: store contiguously the columns
            !
            offset = 0
            DO iproc3 = 1, nproc3
                kdest = (iproc3 - 1) * sendsize
                kfrom = offset
                ierr = cudaMemcpy2D(f_in_d(kfrom + 1), nr3x, f_aux(kdest + 1), nr3px, &
                                    desc%nr3p(iproc3), howmany * ncpx, &
                                    cudaMemcpyHostToDevice)
                !
                offset = offset + desc%nr3p(iproc3)
            END DO
            !
            !----------------------------------------------------------------------------
            ! Clean extra array elements in each stick
            !
            IF (nr3x /= nr3) THEN
                aux = ncp_(desc%mype + 1)
                !
                !$cuf kernel do(3) <<<*,*>>>
                DO j = 0, howmany - 1
                    !
                    DO k = 1, aux
                        !
                        DO i = nr3, nr3x
                            f_in_d(j * ncpx * nr3x + (k - 1) * nr3x + i) = 0.0D0
                        END DO
                        !
                    END DO
                    !
                END DO
                !
            END IF
            !
20          CONTINUE
            !
        END IF
        !
        DEALLOCATE (ncp_)
        !
        CALL env_stop_clock(sub_name)
#endif
        !
        RETURN
        !
98      FORMAT(10('(', 2F12.9, ')'))
99      FORMAT(20('(', 2F12.9, ')'))
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_fft_scatter_many_yz_gpu
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_fft_scatter_2d_gpu(dfft, f_in_d, f_in, nr3x, nxx_, f_aux_d, f_aux, &
                                      ncp_, npp_, isgn)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(env_fft_type_descriptor), TARGET, INTENT(IN) :: dfft
        INTEGER, INTENT(IN) :: nr3x, nxx_, isgn, ncp_(:), npp_(:)
        !
        COMPLEX(DP), DEVICE, INTENT(INOUT) :: f_in_d(nxx_), f_aux_d(nxx_)
        COMPLEX(DP), INTENT(INOUT) :: f_in(nxx_), f_aux(nxx_)
        !
        INTEGER :: cuf_i, cuf_j, nspip
        INTEGER :: istat
        INTEGER, POINTER, DEVICE :: p_ismap_d(:)
        !
#if defined(__MPI)
        INTEGER :: srh(2 * dfft%nproc)
        INTEGER :: k, offset, proc, ierr, me, nprocp, gproc, gcomm, i, kdest, kfrom
        INTEGER :: me_p, nppx, mc, j, npp, nnp, ii, it, ip, ioff, sendsiz
        INTEGER :: ncpx, ipp, nblk, nsiz
        INTEGER :: iter, dest, sorc
        INTEGER :: istatus(MPI_STATUS_SIZE)
        !
        CHARACTER(LEN=80) :: sub_name = 'env_fft_scatter_2d_gpu'
        !
        !--------------------------------------------------------------------------------
        !
        p_ismap_d => dfft%ismap_d
        me = dfft%mype + 1
        nprocp = dfft%nproc
        istat = cudaDeviceSynchronize()
        ncpx = MAXVAL(ncp_)
        nppx = MAXVAL(npp_)
        !
        IF (dfft%nproc == 1) nppx = dfft%nr3x
        ! this should never happend and should go away
        !
        sendsiz = ncpx * nppx
        ierr = 0
        !
        IF (isgn > 0) THEN
            !
            IF (nprocp == 1) GO TO 10
            !
            !----------------------------------------------------------------------------
            ! "forward" scatter from columns to planes
            !
            ! step one: store contiguously the slices
            !
            offset = 0
            !
            DO gproc = 1, nprocp
                kdest = (gproc - 1) * sendsiz
                kfrom = offset
                !
#ifdef __GPU_MPI
                !$cuf kernel do(2) <<<*,*>>>
                DO k = 1, ncp_(me)
                    !
                    DO i = 1, npp_(gproc)
                        !
                        f_aux_d(kdest + i + (k - 1) * nppx) = &
                            f_in_d(kfrom + i + (k - 1) * nr3x)
                        !
                    END DO
                    !
                END DO
#else
                !
                istat = cudaMemcpy2D(f_aux(kdest + 1), nppx, f_in_d(kfrom + 1), nr3x, &
                                     npp_(gproc), ncp_(me), cudaMemcpyDeviceToHost)
                !
                IF (istat) CALL env_fft_error(sub_bane, &
                                              "ERROR cudaMemcpy2D failed : ", istat)
#endif
                !
                offset = offset + npp_(gproc)
            END DO
            !
            ! step two: communication
            !
            gcomm = dfft%comm
            !
            CALL env_start_clock(sub_name)
            !
#ifdef __GPU_MPI
            istat = cudaDeviceSynchronize()
            !
            DO iter = 2, nprocp
                !
                IF (IAND(nprocp, nprocp - 1) == 0) THEN
                    sorc = IEOR(me - 1, iter - 1)
                ELSE
                    sorc = MOD(me - 1 - (iter - 1) + nprocp, nprocp)
                END IF
                !
                CALL MPI_IRECV(f_in_d((sorc) * sendsiz + 1), sendsiz, &
                               MPI_DOUBLE_COMPLEX, sorc, 0, gcomm, srh(iter - 1), ierr)
                !
            END DO
            !
            DO iter = 2, nprocp
                !
                IF (IAND(nprocp, nprocp - 1) == 0) THEN
                    dest = IEOR(me - 1, iter - 1)
                ELSE
                    dest = MOD(me - 1 + (iter - 1), nprocp)
                END IF
                !
                CALL MPI_ISEND(f_aux_d((dest) * sendsiz + 1), sendsiz, &
                               MPI_DOUBLE_COMPLEX, dest, 0, gcomm, &
                               srh(iter + nprocp - 2), ierr)
                !
            END DO
            !
            istat = cudaMemcpyAsync(f_in_d((me - 1) * sendsiz + 1), &
                                    f_aux_d((me - 1) * sendsiz + 1), sendsiz, &
                                    stream=dfft%a2a_comp)
            !
            CALL MPI_WAITALL(2 * nprocp - 2, srh, MPI_STATUSES_IGNORE, ierr)
            istat = cudaDeviceSynchronize()
#else
            !
            CALL mpi_alltoall(f_aux(1), sendsiz, MPI_DOUBLE_COMPLEX, f_in(1), sendsiz, &
                              MPI_DOUBLE_COMPLEX, gcomm, ierr)
#endif
            !
            CALL env_stop_clock(sub_name)
            !
            IF (ABS(ierr) /= 0) CALL env_fft_error(sub_name, 'info<>0', ABS(ierr))
            !
#ifndef __GPU_MPI
            f_in_d(1:sendsiz * dfft%nproc) = f_in(1:sendsiz * dfft%nproc)
#endif
            !
10          CONTINUE
            !
            !$cuf kernel do (1) <<<*,*>>>
            DO i = LBOUND(f_aux_d, 1), UBOUND(f_aux_d, 1)
                f_aux_d(i) = (0.D0, 0.D0)
            END DO
            !
            IF (isgn == 1) THEN
                npp = dfft%nr3p(me)
                nnp = dfft%nnp
                !
                DO ip = 1, dfft%nproc
                    ioff = dfft%iss(ip)
                    nspip = dfft%nsp(ip)
                    !
                    !$cuf kernel do(2) <<<*,*>>>
                    DO cuf_j = 1, npp
                        !
                        DO cuf_i = 1, nspip
                            it = (ip - 1) * sendsiz + (cuf_i - 1) * nppx
                            mc = p_ismap_d(cuf_i + ioff)
                            f_aux_d(mc + (cuf_j - 1) * nnp) = f_in_d(cuf_j + it)
                        END DO
                        !
                    END DO
                    !
                END DO
                !
            ELSE
                !
                npp = dfft%nr3p(me)
                nnp = dfft%nnp
                ip = 1
                !
                DO gproc = 1, dfft%nproc
                    ioff = dfft%iss(ip)
                    nspip = dfft%nsp(ip)
                    !
                    !$cuf kernel do(2) <<<*,*>>>
                    DO cuf_j = 1, npp
                        !
                        DO cuf_i = 1, nspip
                            mc = p_ismap_d(cuf_i + ioff)
                            it = (cuf_i - 1) * nppx + (gproc - 1) * sendsiz
                            f_aux_d(mc + (cuf_j - 1) * nnp) = f_in_d(cuf_j + it)
                        END DO
                        !
                    END DO
                    !
                    ip = ip + 1
                END DO
                !
            END IF
            !
        ELSE
            !
            !----------------------------------------------------------------------------
            ! "backward" scatter from planes to columns
            !
            IF (isgn == -1) THEN
                npp = dfft%nr3p(me)
                nnp = dfft%nnp
                !
                DO ip = 1, dfft%nproc
                    ioff = dfft%iss(ip)
                    nspip = dfft%nsp(ip)
                    !
                    !$cuf kernel do(2) <<<*,*>>>
                    DO cuf_j = 1, npp
                        !
                        DO cuf_i = 1, nspip
                            mc = p_ismap_d(cuf_i + ioff)
                            it = (ip - 1) * sendsiz + (cuf_i - 1) * nppx
                            f_in_d(cuf_j + it) = f_aux_d(mc + (cuf_j - 1) * nnp)
                        END DO
                        !
                    END DO
                    !
                END DO
                !
            ELSE
                npp = dfft%nr3p(me)
                nnp = dfft%nnp
                !
                DO gproc = 1, dfft%nproc
                    ioff = dfft%iss(gproc)
                    nspip = dfft%nsp(gproc)
                    !
                    !$cuf kernel do(2) <<<*,*>>>
                    DO cuf_j = 1, npp
                        !
                        DO cuf_i = 1, nspip
                            mc = p_ismap_d(cuf_i + ioff)
                            it = (cuf_i - 1) * nppx + (gproc - 1) * sendsiz
                            f_in_d(cuf_j + it) = f_aux_d(mc + (cuf_j - 1) * nnp)
                        END DO
                        !
                    END DO
                    !
                END DO
                !
            END IF
            !
            IF (nprocp == 1) GO TO 20
            !
            ! step two: communication
            !
            gcomm = dfft%comm
            !
#ifndef __GPU_MPI
            f_in(1:sendsiz * dfft%nproc) = f_in_d(1:sendsiz * dfft%nproc)
#endif
            !
            CALL env_start_clock(sub_name)
            !
#ifdef __GPU_MPI
            istat = cudaDeviceSynchronize()
            !
            DO iter = 2, nprocp
                !
                IF (IAND(nprocp, nprocp - 1) == 0) THEN
                    sorc = IEOR(me - 1, iter - 1)
                ELSE
                    sorc = MOD(me - 1 - (iter - 1) + nprocp, nprocp)
                END IF
                !
                CALL MPI_IRECV(f_aux_d((sorc) * sendsiz + 1), sendsiz, &
                               MPI_DOUBLE_COMPLEX, sorc, 0, gcomm, srh(iter - 1), ierr)
                !
            END DO
            !
            DO iter = 2, nprocp
                !
                IF (IAND(nprocp, nprocp - 1) == 0) THEN
                    dest = IEOR(me - 1, iter - 1)
                ELSE
                    dest = MOD(me - 1 + (iter - 1), nprocp)
                END IF
                !
                CALL MPI_ISEND(f_in_d((dest) * sendsiz + 1), sendsiz, &
                               MPI_DOUBLE_COMPLEX, dest, 0, &
                               gcomm, srh(iter + nprocp - 2), ierr)
                !
            END DO
            !
            istat = cudaMemcpyAsync(f_aux_d((me - 1) * sendsiz + 1), &
                                    f_in_d((me - 1) * sendsiz + 1), sendsiz, &
                                    stream=dfft%a2a_comp)
            !
            CALL MPI_WAITALL(2 * nprocp - 2, srh, MPI_STATUSES_IGNORE, ierr)
            !
            istat = cudaDeviceSynchronize()
#else
            !
            CALL mpi_alltoall(f_in(1), sendsiz, MPI_DOUBLE_COMPLEX, f_aux(1), sendsiz, &
                              MPI_DOUBLE_COMPLEX, gcomm, ierr)
#endif
            !
            CALL env_stop_clock(sub_name)
            !
            IF (ABS(ierr) /= 0) CALL env_fft_error(sub_name, 'info<>0', ABS(ierr))
            !
            ! step one: store contiguously the columns
            !
            offset = 0
            !
            DO gproc = 1, nprocp
                kdest = (gproc - 1) * sendsiz
                kfrom = offset
#ifdef __GPU_MPI
                !
                !$cuf kernel do(2) <<<*,*>>>
                DO k = 1, ncp_(me)
                    !
                    DO i = 1, npp_(gproc)
                        !
                        f_in_d(kfrom + i + (k - 1) * nr3x) = &
                            f_aux_d(kdest + i + (k - 1) * nppx)
                        !
                    END DO
                    !
                END DO
#else
                !
                istat = cudaMemcpy2D(f_in_d(kfrom + 1), nr3x, f_aux(kdest + 1), nppx, &
                                     npp_(gproc), ncp_(me), cudaMemcpyHostToDevice)
#endif
                !
                offset = offset + npp_(gproc)
            END DO
            !
20          CONTINUE
            !
        END IF
        !
        istat = cudaDeviceSynchronize()
#endif
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_fft_scatter_2d_gpu
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_fft_scatter_many_columns_to_planes_store(dfft, f_in_d, f_in, nr3x, &
                                                            nxx_, f_aux_d, f_aux, &
                                                            f_aux2_d, f_aux2, ncp_, &
                                                            npp_, isgn, batchsize, &
                                                            batch_id)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(env_fft_type_descriptor), TARGET, INTENT(IN) :: dfft
        INTEGER, INTENT(IN) :: nr3x, nxx_, isgn, ncp_(:), npp_(:)
        INTEGER, INTENT(IN) :: batchsize, batch_id
        !
        COMPLEX(DP), DEVICE, INTENT(INOUT) :: f_in_d(batchsize * nxx_), &
                                              f_aux_d(batchsize * nxx_), &
                                              f_aux2_d(batchsize * nxx_)
        !
        COMPLEX(DP), INTENT(INOUT) :: f_in(batchsize * nxx_), &
                                      f_aux(batchsize * nxx_), &
                                      f_aux2(batchsize * nxx_)
        !
        INTEGER :: istat
        INTEGER, POINTER, DEVICE :: p_ismap_d(:)
        !
#if defined(__MPI)
        INTEGER :: k, offset, proc, ierr, me, nprocp, gproc, gcomm, i, kdest, kfrom
        INTEGER :: me_p, nppx, mc, j, npp, nnp, nnr, ii, it, ip, ioff, sendsiz
        INTEGER :: ncpx, ipp, nblk, nsiz
        !
        INTEGER, ALLOCATABLE :: offset_proc(:)
        INTEGER :: iter, dest, sorc
        INTEGER :: istatus(MPI_STATUS_SIZE)
        !
        CHARACTER(LEN=80) :: sub_name = 'env_fft_scatter_many_columns_to_planes_store'
        !
        !--------------------------------------------------------------------------------
        !
        p_ismap_d => dfft%ismap_d
        me = dfft%mype + 1
        nprocp = dfft%nproc
        !
#ifdef __IPC
#ifndef __GPU_MPI
        CALL get_ipc_peers(dfft%IPC_PEER)
#endif
#endif
        !
        ncpx = MAXVAL(ncp_)
        nppx = MAXVAL(npp_)
        !
        IF (dfft%nproc == 1) nppx = dfft%nr3x
        !
        sendsiz = batchsize * ncpx * nppx
        nnr = dfft%nnr
        ierr = 0
        !
        IF (isgn < 0) CALL env_fft_error(sub_name, 'isign is wrong', isgn)
        !
        IF (nprocp == 1) GO TO 10
        !
        !--------------------------------------------------------------------------------
        ! "forward" scatter from columns to planes
        !
        ! step one: store contiguously the slices
        !
        ALLOCATE (offset_proc(nprocp))
        offset = 0
        !
        DO proc = 1, nprocp
            offset_proc(proc) = offset
            offset = offset + npp_(proc)
        END DO
        !
        DO iter = 2, nprocp
            !
            IF (IAND(nprocp, nprocp - 1) == 0) THEN
                dest = IEOR(me - 1, iter - 1)
            ELSE
                dest = MOD(me - 1 + (iter - 1), nprocp)
            END IF
            !
            proc = dest + 1
            kdest = (proc - 1) * sendsiz
            kfrom = offset_proc(proc)
            !
#ifdef __GPU_MPI
            istat = cudaMemcpy2DAsync(f_aux_d(kdest + 1), nppx, f_in_d(kfrom + 1), &
                                      nr3x, npp_(proc), batchsize * ncpx, &
                                      cudaMemcpyDeviceToDevice, dfft%bstreams(batch_id))
            !
            IF (istat /= cudaSuccess) &
                CALL env_fft_error(sub_name, 'cudaMemcpy2DAsync failed : ', istat)
#else
            !
#ifdef __IPC
            IF (dfft%IPC_PEER(dest + 1) .EQ. 1) THEN
                !
                istat = cudaMemcpy2DAsync(f_aux_d(kdest + 1), nppx, f_in_d(kfrom + 1), &
                                          nr3x, npp_(proc), batchsize * ncpx, &
                                          cudaMemcpyDeviceToDevice, &
                                          dfft%bstreams(batch_id))
                !
                IF (istat /= cudaSuccess) &
                    CALL env_fft_error(sub_name, 'cudaMemcpy2DAsync failed : ', istat)
                !
            ELSE
                !
                istat = cudaMemcpy2DAsync(f_aux(kdest + 1), nppx, f_in_d(kfrom + 1), &
                                          nr3x, npp_(proc), batchsize * ncpx, &
                                          cudaMemcpyDeviceToHost, &
                                          dfft%bstreams(batch_id))
                !
                IF (istat /= cudaSuccess) &
                    CALL env_fft_error(sub_name, 'cudaMemcpy2DAsync failed : ', istat)
                !
            END IF
#else
            !
            istat = cudaMemcpy2DAsync(f_aux(kdest + 1), nppx, f_in_d(kfrom + 1), nr3x, &
                                      npp_(proc), batchsize * ncpx, &
                                      cudaMemcpyDeviceToHost, dfft%bstreams(batch_id))
            !
            IF (istat /= cudaSuccess) &
                CALL env_fft_error(sub_name, 'cudaMemcpy2DAsync failed : ', istat)
#endif
            !
#endif
        END DO
        !
        istat = cudaEventRecord(dfft%bevents(batch_id), dfft%bstreams(batch_id))
        DEALLOCATE (offset_proc)
        !
10      CONTINUE
        !
#endif
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_fft_scatter_many_columns_to_planes_store
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_fft_scatter_many_columns_to_planes_send(dfft, f_in_d, f_in, nr3x, &
                                                           nxx_, f_aux_d, f_aux, &
                                                           f_aux2_d, f_aux2, ncp_, &
                                                           npp_, isgn, batchsize, &
                                                           batch_id)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(env_fft_type_descriptor), TARGET, INTENT(IN) :: dfft
        INTEGER, INTENT(IN) :: nr3x, nxx_, isgn, ncp_(:), npp_(:)
        INTEGER, INTENT(IN) :: batchsize, batch_id
        !
        COMPLEX(DP), DEVICE, INTENT(INOUT) :: f_in_d(batchsize * nxx_), &
                                              f_aux_d(batchsize * nxx_), &
                                              f_aux2_d(batchsize * nxx_)
        !
        COMPLEX(DP), INTENT(INOUT) :: f_in(batchsize * nxx_), &
                                      f_aux(batchsize * nxx_), &
                                      f_aux2(batchsize * nxx_)
        !
        INTEGER :: cuf_i, cuf_j, nspip
        INTEGER :: istat
        INTEGER, POINTER, DEVICE :: p_ismap_d(:)
        !
#if defined(__MPI)
        INTEGER :: k, offset, proc, ierr, me, nprocp, gproc, gcomm, i, kdest, kfrom
        INTEGER :: me_p, nppx, mc, j, npp, nnp, nnr, ii, it, ip, ioff, sendsiz
        INTEGER :: ncpx, ipp, nblk, nsiz
        !
        INTEGER :: iter, dest, sorc, req_cnt
        INTEGER :: istatus(MPI_STATUS_SIZE)
        !
        CHARACTER(LEN=80) :: sub_name = 'env_fft_scatter_many_columns_to_planes_send'
        !
        !--------------------------------------------------------------------------------
        !
        p_ismap_d => dfft%ismap_d
        me = dfft%mype + 1
        nprocp = dfft%nproc
        !
        ncpx = MAXVAL(ncp_)
        nppx = MAXVAL(npp_)
        !
        IF (dfft%nproc == 1) nppx = dfft%nr3x
        !
        sendsiz = batchsize * ncpx * nppx
        nnr = dfft%nnr
        !
#ifdef __IPC
        CALL get_ipc_peers(dfft%IPC_PEER)
#endif
        !
        ierr = 0
        !
        IF (isgn < 0) CALL env_fft_error(sub_name, 'isign is wrong', isgn)
        !
        IF (nprocp == 1) GO TO 10
        !
        ! step two: communication
        !
        gcomm = dfft%comm
        !
        istat = cudaEventSynchronize(dfft%bevents(batch_id))
        ! JR Note: Holding off staging receives until buffer is packed
        !
        CALL env_start_clock(sub_name)
        !
#ifdef __IPC
        CALL MPI_Barrier(gcomm, ierr)
#endif
        !
        req_cnt = 0
        !
        DO iter = 2, nprocp
            !
            IF (IAND(nprocp, nprocp - 1) == 0) THEN
                sorc = IEOR(me - 1, iter - 1)
            ELSE
                sorc = MOD(me - 1 - (iter - 1) + nprocp, nprocp)
            END IF
            !
#ifdef __IPC
            IF (dfft%IPC_PEER(sorc + 1) .EQ. 0) THEN
#endif
                !
#ifdef __GPU_MPI
                CALL MPI_IRECV(f_aux2_d((sorc) * sendsiz + 1), sendsiz, &
                               MPI_DOUBLE_COMPLEX, sorc, 0, gcomm, &
                               dfft%srh(req_cnt + 1, batch_id), ierr)
#else
                !
                CALL MPI_IRECV(f_aux2((sorc) * sendsiz + 1), sendsiz, &
                               MPI_DOUBLE_COMPLEX, sorc, 0, gcomm, &
                               dfft%srh(req_cnt + 1, batch_id), ierr)
#endif
                !
                req_cnt = req_cnt + 1
#ifdef __IPC
            END IF
#endif
            !
        END DO
        !
        DO iter = 2, nprocp
            !
            IF (IAND(nprocp, nprocp - 1) == 0) THEN
                dest = IEOR(me - 1, iter - 1)
            ELSE
                dest = MOD(me - 1 + (iter - 1), nprocp)
            END IF
            !
#ifdef __IPC
            IF (dfft%IPC_PEER(dest + 1) .EQ. 1) THEN
                !
                CALL ipc_send(f_aux_d((dest) * sendsiz + 1), sendsiz, &
                              f_aux2_d((me - 1) * sendsiz + 1), 1, dest, gcomm, ierr)
                !
            ELSE
#endif
                !
#ifdef __GPU_MPI
                CALL MPI_ISEND(f_aux_d((dest) * sendsiz + 1), sendsiz, &
                               MPI_DOUBLE_COMPLEX, dest, 0, gcomm, &
                               dfft%srh(req_cnt + 1, batch_id), ierr)
#else
                !
                CALL MPI_ISEND(f_aux((dest) * sendsiz + 1), sendsiz, &
                               MPI_DOUBLE_COMPLEX, dest, 0, gcomm, &
                               dfft%srh(req_cnt + 1, batch_id), ierr)
#endif
                !
                req_cnt = req_cnt + 1
#ifdef __IPC
            END IF
#endif
            !
        END DO
        !
        offset = 0
        !
        DO proc = 1, me - 1
            offset = offset + npp_(proc)
        END DO
        !
        istat = cudaMemcpy2DAsync(f_aux2_d((me - 1) * sendsiz + 1), nppx, &
                                  f_in_d(offset + 1), nr3x, npp_(me), batchsize * ncpx, &
                                  cudaMemcpyDeviceToDevice, dfft%bstreams(batch_id))
        !
        IF (istat /= cudaSuccess) &
            CALL env_fft_error(sub_name, 'cudaMemcpy2DAsync failed : ', istat)
        !
        IF (req_cnt > 0) &
            CALL MPI_WAITALL(req_cnt, dfft%srh(1:req_cnt, batch_id), &
                             MPI_STATUSES_IGNORE, ierr)
        !
#ifdef __IPC
        CALL sync_ipc_sends(gcomm)
        !
        CALL MPI_Barrier(gcomm, ierr)
#endif
        !
        CALL env_stop_clock(sub_name)
        !
        IF (ABS(ierr) /= 0) CALL env_fft_error(sub_name, 'info<>0', ABS(ierr))
        !
#ifndef __GPU_MPI
        DO proc = 1, nprocp
            !
            IF (proc /= me) THEN
#ifdef __IPC
                !
                IF (dfft%IPC_PEER(proc) .EQ. 0) THEN
                    kdest = (proc - 1) * sendsiz
                    !
                    istat = cudaMemcpyAsync(f_aux2_d(kdest + 1), f_aux2(kdest + 1), &
                                            sendsiz, stream=dfft%bstreams(batch_id))
                    !
                END IF
#else
                !
                kdest = (proc - 1) * sendsiz
                !
                istat = cudaMemcpyAsync(f_aux2_d(kdest + 1), f_aux2(kdest + 1), &
                                        sendsiz, stream=dfft%bstreams(batch_id))
#endif
                !
            END IF
            !
        END DO
#endif
        !
        i = cudaEventRecord(dfft%bevents(batch_id), dfft%bstreams(batch_id))
        i = cudaStreamWaitEvent(dfft%a2a_comp, dfft%bevents(batch_id), 0)
        !
10      CONTINUE
        !
        !$cuf kernel do (1) <<<*,*,0,dfft%a2a_comp>>>
        DO i = LBOUND(f_aux_d, 1), UBOUND(f_aux_d, 1)
            f_aux_d(i) = (0.D0, 0.D0)
        END DO
        !
        npp = dfft%nr3p(me)
        nnp = dfft%nnp
        !
        DO ip = 1, nprocp
            ioff = dfft%iss(ip)
            nspip = dfft%nsp(ip)
            !
            !$cuf kernel do(3) <<<*,*,0,dfft%a2a_comp>>>
            DO i = 0, batchsize - 1
                !
                DO cuf_j = 1, npp
                    !
                    DO cuf_i = 1, nspip
                        it = (ip - 1) * sendsiz + (cuf_i - 1) * nppx + i * nppx * ncpx
                        mc = p_ismap_d(cuf_i + ioff)
                        f_aux_d(mc + (cuf_j - 1) * nnp + i * nnr) = f_aux2_d(cuf_j + it)
                    END DO
                    !
                END DO
                !
            END DO
            !
        END DO
#endif
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_fft_scatter_many_columns_to_planes_send
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_fft_scatter_many_planes_to_columns_store(dfft, f_in_d, f_in, nr3x, &
                                                            nxx_, f_aux_d, f_aux, &
                                                            f_aux2_d, f_aux2, ncp_, &
                                                            npp_, isgn, batchsize, &
                                                            batch_id)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(env_fft_type_descriptor), TARGET, INTENT(IN) :: dfft
        INTEGER, INTENT(IN) :: nr3x, nxx_, isgn, ncp_(:), npp_(:)
        INTEGER, INTENT(IN) :: batchsize, batch_id
        !
        COMPLEX(DP), DEVICE, INTENT(INOUT) :: f_in_d(batchsize * nxx_), &
                                              f_aux_d(batchsize * nxx_), &
                                              f_aux2_d(batchsize * nxx_)
        !
        COMPLEX(DP), INTENT(INOUT) :: f_in(batchsize * nxx_), &
                                      f_aux(batchsize * nxx_), &
                                      f_aux2(batchsize * nxx_)
        !
        INTEGER :: cuf_i, cuf_j, nspip
        INTEGER :: istat
        INTEGER, POINTER, DEVICE :: p_ismap_d(:)
        !
#if defined(__MPI)
        INTEGER :: k, offset, proc, ierr, me, nprocp, gproc, gcomm, i, kdest, kfrom
        INTEGER :: me_p, nppx, mc, j, npp, nnp, nnr, ii, it, ip, ioff, sendsiz
        INTEGER :: ncpx, ipp, nblk, nsiz
        !
        LOGICAL :: use_tg
        INTEGER :: iter, dest, sorc
        INTEGER :: istatus(MPI_STATUS_SIZE)
        !
        CHARACTER(LEN=80) :: sub_name = 'env_fft_scatter_many_planes_to_columns_store'
        !
        !--------------------------------------------------------------------------------
        !
        p_ismap_d => dfft%ismap_d
        me = dfft%mype + 1
        nprocp = dfft%nproc
        !
#ifdef __IPC
#ifndef __GPU_MPI
        CALL get_ipc_peers(dfft%IPC_PEER)
#endif
#endif
        !
        ncpx = MAXVAL(ncp_)
        nppx = MAXVAL(npp_)
        !
        IF (dfft%nproc == 1) nppx = dfft%nr3x
        !
        sendsiz = batchsize * ncpx * nppx
        nnr = dfft%nnr
        ierr = 0
        !
        IF (isgn > 0) CALL env_fft_error(sub_name, 'isign is wrong', isgn)
        !
        !--------------------------------------------------------------------------------
        ! "backward" scatter from planes to columns
        !
        npp = dfft%nr3p(me)
        nnp = dfft%nnp
        !
        DO iter = 1, nprocp
            !
            IF (IAND(nprocp, nprocp - 1) == 0) THEN
                dest = IEOR(me - 1, iter - 1)
            ELSE
                dest = MOD(me - 1 + (iter - 1), nprocp)
            END IF
            !
            ip = dest + 1
            ioff = dfft%iss(ip)
            nspip = dfft%nsp(ip)
            !
            !$cuf kernel do(3) <<<*,*,0,dfft%a2a_comp>>>
            DO i = 0, batchsize - 1
                !
                DO cuf_j = 1, npp
                    !
                    DO cuf_i = 1, nspip
                        mc = p_ismap_d(cuf_i + ioff)
                        it = (ip - 1) * sendsiz + (cuf_i - 1) * nppx + i * nppx * ncpx
                        f_aux2_d(cuf_j + it) = f_aux_d(mc + (cuf_j - 1) * nnp + i * nnr)
                    END DO
                    !
                END DO
                !
            END DO
            !
        END DO
        !
#ifndef __GPU_MPI
        i = cudaEventRecord(dfft%bevents(batch_id), dfft%a2a_comp)
        i = cudaStreamWaitEvent(dfft%bstreams(batch_id), dfft%bevents(batch_id), 0)
        !
        DO proc = 1, nprocp
            !
            IF (proc /= me) THEN
                !
#ifdef __IPC
                IF (dfft%IPC_PEER(proc) .EQ. 0) THEN
                    kdest = (proc - 1) * sendsiz
                    !
                    istat = cudaMemcpyAsync(f_aux2(kdest + 1), f_aux2_d(kdest + 1), &
                                            sendsiz, stream=dfft%bstreams(batch_id))
                    !
                END IF
#else
                kdest = (proc - 1) * sendsiz
                !
                istat = cudaMemcpyAsync(f_aux2(kdest + 1), f_aux2_d(kdest + 1), &
                                        sendsiz, stream=dfft%bstreams(batch_id))
#endif
                !
            END IF
            !
        END DO
#endif
        !
#ifdef __GPU_MPI
        istat = cudaEventRecord(dfft%bevents(batch_id), dfft%a2a_comp)
#else
        istat = cudaEventRecord(dfft%bevents(batch_id), dfft%bstreams(batch_id))
#endif
#endif
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_fft_scatter_many_planes_to_columns_store
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_fft_scatter_many_planes_to_columns_send(dfft, f_in_d, f_in, nr3x, &
                                                           nxx_, f_aux_d, f_aux, &
                                                           f_aux2_d, f_aux2, ncp_, &
                                                           npp_, isgn, batchsize, &
                                                           batch_id)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(env_fft_type_descriptor), TARGET, INTENT(IN) :: dfft
        INTEGER, INTENT(IN) :: nr3x, nxx_, isgn, ncp_(:), npp_(:)
        INTEGER, INTENT(IN) :: batchsize, batch_id
        !
        COMPLEX(DP), DEVICE, INTENT(INOUT) :: f_in_d(batchsize * nxx_), &
                                              f_aux_d(batchsize * nxx_), &
                                              f_aux2_d(batchsize * nxx_)
        !
        COMPLEX(DP), INTENT(INOUT) :: f_in(batchsize * nxx_), &
                                      f_aux(batchsize * nxx_), &
                                      f_aux2(batchsize * nxx_)
        !
        INTEGER :: istat
        INTEGER, POINTER, DEVICE :: p_ismap_d(:)
        !
#if defined(__MPI)
        INTEGER :: k, offset, proc, ierr, me, nprocp, gproc, gcomm, i, kdest, kfrom
        INTEGER :: me_p, nppx, mc, j, npp, nnp, nnr, ii, it, ip, ioff, sendsiz
        INTEGER :: ncpx, ipp, nblk, nsiz
        !
        INTEGER :: iter, dest, sorc, req_cnt
        INTEGER :: istatus(MPI_STATUS_SIZE)
        !
        CHARACTER(LEN=80) :: sub_name = 'env_fft_scatter_many_planes_to_columns_send'
        !
        !--------------------------------------------------------------------------------
        !
        p_ismap_d => dfft%ismap_d
        me = dfft%mype + 1
        nprocp = dfft%nproc
        !
        ncpx = MAXVAL(ncp_)
        ! max number of sticks among processors ( should be of wave func )
        !
        nppx = MAXVAL(npp_)
        ! max size of the "Z" section of each processor in the nproc3 group along Z
        !
        IF (dfft%nproc == 1) nppx = dfft%nr3x
        !
        sendsiz = batchsize * ncpx * nppx
        nnr = dfft%nnr
        !
#ifdef __IPC
        CALL get_ipc_peers(dfft%IPC_PEER)
#endif
        !
        ierr = 0
        !
        IF (isgn > 0) CALL env_fft_error(sub_name, 'isign is wrong', isgn)
        !
        !--------------------------------------------------------------------------------
        ! "backward" scatter from planes to columns
        !
        IF (nprocp == 1) GO TO 20
        !
        !--------------------------------------------------------------------------------
        ! Communication takes place here:
        ! fractions of sticks will be moved to form complete sticks
        !
        gcomm = dfft%comm
        !
        istat = cudaEventSynchronize(dfft%bevents(batch_id))
        ! JR Note: Holding off staging receives until buffer is packed.
        !
        CALL env_start_clock(sub_name)
        !
#ifdef __IPC
        CALL MPI_Barrier(gcomm, ierr)
#endif
        !
        req_cnt = 0
        !
        DO iter = 2, nprocp
            !
            IF (IAND(nprocp, nprocp - 1) == 0) THEN
                sorc = IEOR(me - 1, iter - 1)
            ELSE
                sorc = MOD(me - 1 - (iter - 1) + nprocp, nprocp)
            END IF
            !
#ifdef __IPC
            IF (dfft%IPC_PEER(sorc + 1) .EQ. 0) THEN
#endif
                !
#ifdef __GPU_MPI
                CALL MPI_IRECV(f_aux_d((sorc) * sendsiz + 1), sendsiz, &
                               MPI_DOUBLE_COMPLEX, sorc, 0, gcomm, &
                               dfft%srh(req_cnt + 1, batch_id), ierr)
#else
                !
                CALL MPI_IRECV(f_aux((sorc) * sendsiz + 1), sendsiz, &
                               MPI_DOUBLE_COMPLEX, sorc, 0, gcomm, &
                               dfft%srh(req_cnt + 1, batch_id), ierr)
#endif
                !
                req_cnt = req_cnt + 1
#ifdef __IPC
            END IF
#endif
            !
        END DO
        !
        DO iter = 2, nprocp
            !
            IF (IAND(nprocp, nprocp - 1) == 0) THEN
                dest = IEOR(me - 1, iter - 1)
            ELSE
                dest = MOD(me - 1 + (iter - 1), nprocp)
            END IF
            !
#ifdef __IPC
            IF (dfft%IPC_PEER(dest + 1) .EQ. 1) THEN
                !
                CALL ipc_send(f_aux2_d((dest) * sendsiz + 1), sendsiz, &
                              f_aux_d((me - 1) * sendsiz + 1), 0, dest, gcomm, ierr)
                !
            ELSE
#endif
                !
#ifdef __GPU_MPI
                CALL MPI_ISEND(f_aux2_d((dest) * sendsiz + 1), sendsiz, &
                               MPI_DOUBLE_COMPLEX, dest, 0, gcomm, &
                               dfft%srh(req_cnt + 1, batch_id), ierr)
#else
                !
                CALL MPI_ISEND(f_aux2((dest) * sendsiz + 1), sendsiz, &
                               MPI_DOUBLE_COMPLEX, dest, 0, gcomm, &
                               dfft%srh(req_cnt + 1, batch_id), ierr)
#endif
                !
                req_cnt = req_cnt + 1
#ifdef __IPC
            END IF
#endif
            !
        END DO
        !
        !--------------------------------------------------------------------------------
        ! Move the data that we already have (and therefore doesn't to pass 
        ! through MPI avove) directly from f_aux_2 to f_in. The rest will be 
        ! done below
        !
        offset = 0
        !
        DO proc = 1, me - 1
            offset = offset + npp_(proc)
        END DO
        !
        istat = cudaMemcpy2DAsync(f_in_d(offset + 1), nr3x, &
                                  f_aux2_d((me - 1) * sendsiz + 1), nppx, npp_(me), &
                                  batchsize * ncpx, cudaMemcpyDeviceToDevice, &
                                  dfft%bstreams(batch_id))
        !
        IF (req_cnt > 0) &
            CALL MPI_WAITALL(req_cnt, dfft%srh(1:req_cnt, batch_id), &
                             MPI_STATUSES_IGNORE, ierr)
        !
#ifdef __IPC
        CALL sync_ipc_sends(gcomm)
        !
        CALL MPI_Barrier(gcomm, ierr)
#endif
        !
        CALL env_stop_clock(sub_name)
        !
        IF (ABS(ierr) /= 0) CALL env_fft_error(sub_name, 'info<>0', ABS(ierr))
        !
        !--------------------------------------------------------------------------------
        ! Store contiguously the (remaining) columns (one already done above).
        !
        offset = 0
        !
        DO gproc = 1, nprocp
            kdest = (gproc - 1) * sendsiz
            kfrom = offset
            !
            IF (gproc /= me) THEN ! (me already done above)
                !
#ifdef __GPU_MPI
                istat = cudaMemcpy2DAsync(f_in_d(kfrom + 1), nr3x, f_aux_d(kdest + 1), &
                                          nppx, npp_(gproc), batchsize * ncpx, &
                                          cudaMemcpyDeviceToDevice, &
                                          dfft%bstreams(batch_id))
                !
                IF (istat /= cudaSuccess) &
                    CALL env_fft_error(sub_name, 'cudaMemcpy2DAsync failed : ', istat)
#else
                !
#ifdef __IPC
                IF (dfft%IPC_PEER(gproc) .EQ. 1) THEN
                    !
                    istat = cudaMemcpy2DAsync(f_in_d(kfrom + 1), nr3x, &
                                              f_aux_d(kdest + 1), nppx, npp_(gproc), &
                                              batchsize * ncpx, &
                                              cudaMemcpyDeviceToDevice, &
                                              dfft%bstreams(batch_id))
                    !
                ELSE
                    !
                    istat = cudaMemcpy2DAsync(f_in_d(kfrom + 1), nr3x, &
                                              f_aux(kdest + 1), nppx, npp_(gproc), &
                                              batchsize * ncpx, &
                                              cudaMemcpyHostToDevice, &
                                              dfft%bstreams(batch_id))
                    !
                END IF
                !
                IF (istat /= cudaSuccess) &
                    CALL env_fft_error(sub_name, 'cudaMemcpy2DAsync failed : ', istat)
#else
                !
                istat = cudaMemcpy2DAsync(f_in_d(kfrom + 1), nr3x, f_aux(kdest + 1), &
                                          nppx, npp_(gproc), batchsize * ncpx, &
                                          cudaMemcpyHostToDevice, &
                                          dfft%bstreams(batch_id))
                !
                IF (istat /= cudaSuccess) &
                    CALL env_fft_error(sub_name, 'cudaMemcpy2DAsync failed : ', istat)
#endif
#endif
                !
            END IF
            !
            offset = offset + npp_(gproc)
        END DO
        !
20      CONTINUE
#endif
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_fft_scatter_many_planes_to_columns_send
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE env_fft_scatter_gpu
!----------------------------------------------------------------------------------------
#endif
