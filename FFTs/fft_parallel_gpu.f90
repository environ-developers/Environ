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
! Authors:  Carlo Cavazzoni (Apr. 2009)
!           Stefano de Gironcoli (Sep-Nov 2016)
!
!----------------------------------------------------------------------------------------
!
#if defined(__CUDA)
!
!----------------------------------------------------------------------------------------
!>
!! GPU version of env_fft_parallel.f90
!!
!----------------------------------------------------------------------------------------
MODULE env_fft_parallel_gpu
    !------------------------------------------------------------------------------------
    !
    USE cudafor
    !
    USE env_fft_param
    !
    USE env_fft_types, ONLY: env_fft_type_descriptor
    !
    USE env_fft_buffers, ONLY: env_check_fft_buffers_size, &
                               f_h => pin_space_scatter_in, &
                               aux_h => pin_space_scatter_out, &
                               aux_d => dev_space_fftparallel, &
                               aux2_h => pin_space_scatter_dblbuffer, &
                               aux2_d => dev_space_scatter_dblbuffer
    !
    USE env_fft_scalar, ONLY: env_cft_1z_gpu, env_cft_2xy_gpu
    !
    USE env_fft_scatter_gpu
    !
    !------------------------------------------------------------------------------------
    !
    PRIVATE
    !
    PUBLIC :: env_tg_cft3s_gpu, env_many_cft3s_gpu, env_tg_cft3s_2d_gpu, &
              env_many_cft3s_2d_gpu
    !
    !------------------------------------------------------------------------------------
CONTAINS
    !------------------------------------------------------------------------------------
    !>
    !! General purpose driver, GPU version
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_tg_cft3s_gpu(f_d, dfft, isgn)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(env_fft_type_descriptor), INTENT(IN) :: dfft ! descriptor of fft data layout
        INTEGER, INTENT(IN) :: isgn ! fft direction
        !
        COMPLEX(DP), DEVICE, INTENT(INOUT) :: f_d(:)
        ! array containing data to be transformed
        !
        INTEGER :: n1, n2, n3, nx1, nx2, nx3
        INTEGER :: nnr_
        INTEGER :: nsticks_x, nsticks_y, nsticks_z
        INTEGER :: ierr, i
        INTEGER(kind=cuda_stream_kind) :: stream = 0
        !
        CHARACTER(LEN=80) :: sub_name = 'env_tg_cft3s_gpu'
        !
        !--------------------------------------------------------------------------------
        !
        n1 = dfft%nr1
        n2 = dfft%nr2
        n3 = dfft%nr3
        nx1 = dfft%nr1x
        nx2 = dfft%nr2x
        nx3 = dfft%nr3x
        !
        nnr_ = dfft%nnr
        nsticks_x = dfft%my_nr2p * dfft%my_nr3p
        nsticks_y = dfft%nr1p(dfft%mype2 + 1) * dfft%my_nr3p
        nsticks_z = dfft%nsp(dfft%mype + 1)
        !
        CALL env_check_fft_buffers_size(dfft)
        !
        IF (isgn > 0) THEN ! G -> R
            !
            CALL env_cft_1z_gpu(f_d, nsticks_z, n3, nx3, isgn, aux_d, stream, &
                                in_place=.TRUE.)
            !
            CALL env_fft_scatter_yz_gpu(dfft, f_d, aux_d, nnr_, isgn)
            !
            CALL env_cft_1z_gpu(aux_d, nsticks_y, n2, nx2, isgn, f_d, stream)
            !
            CALL env_fft_scatter_xy_gpu(dfft, f_d, aux_d, nnr_, isgn, stream)
            !
            CALL env_cft_1z_gpu(aux_d, nsticks_x, n1, nx1, isgn, f_d, stream)
            !
            !----------------------------------------------------------------------------
            ! Clean garbage beyond the intended dimension
            ! Should not be needed but apparently it is
            !
            IF (nsticks_x * nx1 < nnr_) THEN
                !
                !$cuf kernel do(1)<<<*,*,0,stream>>>
                DO i = nsticks_x * nx1 + 1, nnr_
                    f_d(i) = (0.0_DP, 0.0_DP)
                END DO
                !
            END IF
            !
        ELSE ! R -> G
            !
            CALL env_cft_1z_gpu(f_d, nsticks_x, n1, nx1, isgn, aux_d, stream)
            !
            CALL env_fft_scatter_xy_gpu(dfft, f_d, aux_d, nnr_, isgn, stream)
            !
            CALL env_cft_1z_gpu(f_d, nsticks_y, n2, nx2, isgn, aux_d, stream)
            !
            CALL env_fft_scatter_yz_gpu(dfft, f_d, aux_d, nnr_, isgn)
            !
            CALL env_cft_1z_gpu(f_d, nsticks_z, n3, nx3, isgn, aux_d, stream, &
                                in_place=.TRUE.)
            !
        END IF
        !
        RETURN
        !
99      FORMAT(20('(', 2F12.9, ')'))
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_tg_cft3s_gpu
    !------------------------------------------------------------------------------------
    !>
    !! Specific driver for the new 'many' call, GPU version
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_many_cft3s_gpu(f_d, dfft, isgn, howmany)
        !----------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(env_fft_type_descriptor), INTENT(IN) :: dfft ! descriptor of fft data layout
        INTEGER, INTENT(IN) :: howmany ! number of FFTs grouped together
        INTEGER, INTENT(IN) :: isgn ! fft direction
        !
        COMPLEX(DP), DEVICE, INTENT(INOUT) :: f_d(:)
        ! array containing data to be transformed
        !
        INTEGER :: n1, n2, n3, nx1, nx2, nx3
        INTEGER :: nnr_
        INTEGER :: nsticks_x, nsticks_y, nsticks_z
        INTEGER :: nsticks_zx
        INTEGER :: ierr, i, j
        INTEGER :: nstreams
        INTEGER(kind=cuda_stream_kind) :: stream = 0 ! cuda_default_stream
        !
        CHARACTER(LEN=80) :: sub_name = 'env_many_cft3s_gpu'
        !
        !--------------------------------------------------------------------------------
        !
        n1 = dfft%nr1
        n2 = dfft%nr2
        n3 = dfft%nr3
        nx1 = dfft%nr1x
        nx2 = dfft%nr2x
        nx3 = dfft%nr3x
        !
        nnr_ = dfft%nnr
        nsticks_x = dfft%my_nr2p * dfft%my_nr3p
        nsticks_y = dfft%nr1p(dfft%mype2 + 1) * dfft%my_nr3p
        nsticks_z = dfft%nsp(dfft%mype + 1)
        nsticks_zx = MAXVAL(dfft%nsp)
        !
        CALL env_check_fft_buffers_size(dfft, howmany)
        !
        nstreams = dfft%nstream_many
        ierr = cudaDeviceSynchronize()
        !
        IF (isgn > 0) THEN ! G -> R
            !
            DO i = 0, howmany - 1
                !
                CALL env_cft_1z_gpu(f_d(i * nnr_ + 1:), nsticks_z, n3, nx3, isgn, &
                                    aux_d(nx3 * nsticks_zx * i + 1:), stream)
                !
            END DO
            !
            ! this brings back data from packed to unpacked
            CALL env_fft_scatter_many_yz_gpu(dfft, aux_d(1), f_d(1), howmany * nnr_, &
                                             isgn, howmany)
            !
            DO i = 0, howmany - 1
                !
                CALL env_cft_1z_gpu(f_d(i * nnr_ + 1:), nsticks_y, n2, nx2, isgn, &
                                    aux_d(i * nnr_ + 1:), stream, in_place=.TRUE.)
                !
            END DO
            !
            DO i = 0, howmany - 1
                !
                CALL env_fft_scatter_xy_gpu(dfft, f_d(i * nnr_ + 1:), &
                                            aux_d(i * nnr_ + 1:), nnr_, isgn, stream)
                !
            END DO
            !
            DO i = 0, howmany - 1
                !
                CALL env_cft_1z_gpu(aux_d(i * nnr_ + 1:), nsticks_x, n1, nx1, isgn, &
                                    _d(i * nnr_ + 1:), stream)
                !
            END DO
            !
            !----------------------------------------------------------------------------
            ! Clean garbage beyond the intended dimension
            ! Should not be needed but apparently it is
            !
            IF (nsticks_x * nx1 < nnr_) THEN
                !
                !$cuf kernel do(2)<<<*,*>>>
                DO i = 0, howmany - 1
                    !
                    DO j = nsticks_x * nx1 + 1, nnr_
                        f_d(j + i * nnr_) = (0.0_DP, 0.0_DP)
                    END DO
                    !
                END DO
                !
            END IF
            !
        ELSE ! R -> G
            !
            DO i = 0, howmany - 1
                !
                CALL env_cft_1z_gpu(f_d(i * nnr_ + 1:), nsticks_x, n1, nx1, isgn, &
                                    aux_d(i * nnr_ + 1:), stream)
                !
            END DO
            !
            DO i = 0, howmany - 1
                !
                CALL env_fft_scatter_xy_gpu(dfft, f_d(i * nnr_ + 1), &
                                            aux_d(i * nnr_ + 1), nnr_, isgn, stream)
                !
            END DO
            !
            DO i = 0, howmany - 1
                !
                CALL env_cft_1z_gpu(f_d(i * nnr_ + 1:), nsticks_y, n2, nx2, isgn, &
                                    aux_d(i * nnr_ + 1:), stream, in_place=.TRUE.)
                !
            END DO
            !
            CALL env_fft_scatter_many_yz_gpu(dfft, aux_d, f_d, howmany * nnr_, isgn, &
                                             howmany)
            !
            DO i = 0, howmany - 1
                !
                CALL env_cft_1z_gpu(aux_d(nx3 * nsticks_zx * i + 1:), nsticks_z, &
                                    n3, nx3, isgn, f_d(i * nnr_ + 1:), stream)
                !
            END DO
            !
            !----------------------------------------------------------------------------
            ! Clean garbage beyond the intended dimension
            ! Should not be needed but apparently it is
            !
            IF (nsticks_z * nx3 < nnr_) THEN
                !
                !$cuf kernel do(2)<<<*,*>>>
                DO i = 0, howmany - 1
                    !
                    DO j = nsticks_z * nx3 + 1, nnr_
                        f_d(j + i * nnr_) = (0.0_DP, 0.0_DP)
                    END DO
                    !
                END DO
                !
            END IF
            !
        END IF
        !
        RETURN
        !
99      FORMAT(20('(', 2F12.9, ')'))
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_many_cft3s_gpu
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_tg_cft3s_2d_gpu(f_d, dfft, isgn)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(env_fft_type_descriptor), INTENT(IN) :: dfft
        INTEGER, INTENT(IN) :: isgn ! fft direction
        ! descriptor of fft data layout
        !
        COMPLEX(DP), DEVICE, INTENT(INOUT) :: f_d(dfft%nnr)
        ! array containing data to be transformed
        !
        INTEGER :: me_p, istat
        INTEGER :: n1, n2, n3, nx1, nx2, nx3
        COMPLEX(DP), ALLOCATABLE :: yf(:)
        INTEGER :: planes(dfft%nr1x)
        INTEGER(kind=cuda_stream_kind) :: stream = 0
        !
        CHARACTER(LEN=80) :: sub_name = 'env_tg_cft3s_2d_gpu'
        !
        !--------------------------------------------------------------------------------
        !
        n1 = dfft%nr1
        n2 = dfft%nr2
        n3 = dfft%nr3
        nx1 = dfft%nr1x
        nx2 = dfft%nr2x
        nx3 = dfft%nr3x
        !
        IF (dfft%has_task_groups) &
            CALL env_fft_error(sub_name, &
                               'task groups in 2D + 1D decomposition not implemented', 1)
        !
        CALL env_check_fft_buffers_size(dfft)
        !
        me_p = dfft%mype + 1
        !
        !--------------------------------------------------------------------------------
        ! Transpose data for the 2-D FFT on the x-y plane
        !
        ! NOGRP*dfft%nnr: The length of aux and f
        ! nr3x          : The length of each Z-stick
        ! aux           : input - output
        ! f             : working space
        ! isgn          : type of scatter
        ! dfft%nsw(me)  : holds the number of Z-sticks proc. me has.
        ! dfft%nr3p     : number of planes per processor
        !
        IF (isgn > 0) THEN
            !
            CALL env_cft_1z_gpu(f_d, dfft%nsp(me_p), n3, nx3, isgn, aux_d, stream)
            !
            planes = dfft%iplp
            !
            ! forward scatter from stick to planes
            CALL env_fft_scatter_2d_gpu(dfft, aux_d, aux_h, nx3, dfft%nnr, f_d, f_h, &
                                        dfft%nsp, dfft%nr3p, isgn)
            !
            CALL env_cft_2xy_gpu(f_d, aux_d, dfft%my_nr3p, n1, n2, nx1, nx2, isgn, &
                                 stream, planes)
            !
        ELSE
            planes = dfft%iplp
            !
            CALL env_cft_2xy_gpu(f_d, aux_d, dfft%my_nr3p, n1, n2, nx1, nx2, isgn, &
                                 stream, planes)
            !
            ! backward scatter from stick to planes
            CALL env_fft_scatter_2d_gpu(dfft, aux_d, aux_h, nx3, dfft%nnr, f_d, f_h, &
                                        dfft%nsp, dfft%nr3p, isgn)
            !
            CALL env_cft_1z_gpu(aux_d, dfft%nsp(me_p), n3, nx3, isgn, f_d, stream)
            !
        END IF
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_tg_cft3s_2d_gpu
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_many_cft3s_2d_gpu(f_d, dfft, isgn, batchsize)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(env_fft_type_descriptor), INTENT(IN) :: dfft ! descriptor of fft data layout
        INTEGER, INTENT(IN) :: isgn ! fft direction
        INTEGER, INTENT(IN) :: batchsize
        !
        COMPLEX(DP), DEVICE, INTENT(INOUT) :: f_d(batchsize * dfft%nnr)
        ! array containing data to be transformed
        !
        INTEGER :: me_p, istat, i, j, currsize
        INTEGER :: n1, n2, n3, nx1, nx2, nx3, ncpx, nppx, proc
        COMPLEX(DP), ALLOCATABLE :: yf(:)
        INTEGER :: planes(dfft%nr1x)
        INTEGER :: sticks(dfft%nproc)
        INTEGER(kind=cuda_stream_kind) :: stream = 0
        !
        CHARACTER(LEN=80) :: sub_name = 'env_many_cft3s_2d_gpu'
        !
        !--------------------------------------------------------------------------------
        !
        n1 = dfft%nr1
        n2 = dfft%nr2
        n3 = dfft%nr3
        nx1 = dfft%nr1x
        nx2 = dfft%nr2x
        nx3 = dfft%nr3x
        !
        CALL env_check_fft_buffers_size(dfft, batchsize)
        !
        me_p = dfft%mype + 1
        !
        ncpx = 0
        nppx = 0
        !
        DO proc = 1, dfft%nproc
            !
            IF (ABS(isgn) == 1) ncpx = MAX(ncpx, dfft%nsp(proc))
            !
            nppx = MAX(nppx, dfft%nr3p(proc))
        END DO
        !
        sticks = dfft%nsp
        !
        IF (dfft%nproc <= 1) &
            CALL env_fft_error(sub_name, &
                               'this subroutine should never be called with nproc= ', &
                               dfft%nproc)
        !
        IF (isgn > 0) THEN
            DO j = 0, batchsize - 1, dfft%subbatchsize
                currsize = MIN(dfft%subbatchsize, batchsize - j)
                !
                planes = dfft%iplp
                !
                DO i = 0, currsize - 1
                    !
                    CALL env_cft_1z_gpu(f_d((j + i) * dfft%nnr + 1:), &
                                        sticks(me_p), n3, nx3, isgn, &
                                        aux_d(j * dfft%nnr + i * ncpx * nx3 + 1:), &
                                        dfft%a2a_comp)
                    !
                END DO
                !
                i = cudaEventRecord(dfft%bevents(j / dfft%subbatchsize + 1), &
                                    dfft%a2a_comp)
                !
                i = cudaStreamWaitEvent(dfft%bstreams(j / dfft%subbatchsize + 1), &
                                        dfft%bevents(j / dfft%subbatchsize + 1), 0)
                !
                IF (j > 0) &
                    i = cudaStreamWaitEvent(dfft%bstreams(j / dfft%subbatchsize + 1), &
                                            dfft%bevents(j / dfft%subbatchsize), 0)

                CALL env_fft_scatter_many_columns_to_planes_store( &
                    dfft, aux_d(j * dfft%nnr + 1:), aux_h(j * dfft%nnr + 1:), nx3, &
                    dfft%nnr, f_d(j * dfft%nnr + 1:), f_h(j * dfft%nnr + 1:), &
                    aux2_d(j * dfft%nnr + 1:), aux2_h(j * dfft%nnr + 1:), &
                    sticks, dfft%nr3p, isgn, currsize, j / dfft%subbatchsize + 1)
                !
            END DO
            !
            DO j = 0, batchsize - 1, dfft%subbatchsize
                currsize = MIN(dfft%subbatchsize, batchsize - j)
                !
                CALL env_fft_scatter_many_columns_to_planes_send( &
                    dfft, aux_d(j * dfft%nnr + 1:), aux_h(j * dfft%nnr + 1:), nx3, &
                    dfft%nnr, f_d(j * dfft%nnr + 1:), f_h(j * dfft%nnr + 1:), &
                    aux2_d(j * dfft%nnr + 1:), aux2_h(j * dfft%nnr + 1:), sticks, &
                    dfft%nr3p, isgn, currsize, j / dfft%subbatchsize + 1)
                !
                IF (currsize == dfft%subbatchsize) THEN
                    !
                    CALL env_cft_2xy_gpu(f_d(j * dfft%nnr + 1:), &
                                         aux_d(j * dfft%nnr + 1:), currsize * nppx, &
                                         n1, n2, nx1, nx2, isgn, dfft%a2a_comp, planes)
                    !
                ELSE
                    !
                    DO i = 0, currsize - 1
                        !
                        CALL env_cft_2xy_gpu(f_d((j + i) * dfft%nnr + 1:), &
                                             aux_d((j + i) * dfft%nnr + 1:), &
                                             dfft%nr3p(me_p), n1, n2, nx1, nx2, isgn, &
                                             dfft%a2a_comp, planes)
                        !
                    END DO
                    !
                END IF
                !
            END DO
            !
        ELSE
            !
            DO j = 0, batchsize - 1, dfft%subbatchsize
                currsize = MIN(dfft%subbatchsize, batchsize - j)
                !
                planes = dfft%iplp
                !
                IF (currsize == dfft%subbatchsize) THEN
                    !
                    CALL env_cft_2xy_gpu(f_d(j * dfft%nnr + 1:), &
                                         aux_d(j * dfft%nnr + 1:), currsize * nppx, &
                                         n1, n2, nx1, nx2, isgn, dfft%a2a_comp, planes)
                    !
                ELSE
                    !
                    DO i = 0, currsize - 1
                        !
                        CALL env_cft_2xy_gpu(f_d((j + i) * dfft%nnr + 1:), &
                                             aux_d((j + i) * dfft%nnr + 1:), &
                                             dfft%nr3p(me_p), n1, n2, nx1, nx2, isgn, &
                                             dfft%a2a_comp, planes)
                        !
                    END DO
                    !
                END IF
                !
                IF (j > 0) &
                    i = cudaStreamWaitEvent(dfft%bstreams(j / dfft%subbatchsize + 1), &
                                            dfft%bevents(j / dfft%subbatchsize), 0)
                !
                CALL env_fft_scatter_many_planes_to_columns_store( &
                    dfft, aux_d(j * dfft%nnr + 1:), aux_h(j * dfft%nnr + 1:), nx3, &
                    dfft%nnr, f_d(j * dfft%nnr + 1:), f_h(j * dfft%nnr + 1:), &
                    aux2_d(j * dfft%nnr + 1:), aux2_h(j * dfft%nnr + 1:), sticks, &
                    dfft%nr3p, isgn, currsize, j / dfft%subbatchsize + 1)
                !
            END DO
            !
            DO j = 0, batchsize - 1, dfft%subbatchsize
                currsize = MIN(dfft%subbatchsize, batchsize - j)
                !
                CALL env_fft_scatter_many_planes_to_columns_send( &
                    dfft, aux_d(j * dfft%nnr + 1:), aux_h(j * dfft%nnr + 1:), nx3, &
                    dfft%nnr, f_d(j * dfft%nnr + 1:), f_h(j * dfft%nnr + 1:), &
                    aux2_d(j * dfft%nnr + 1:), aux2_h(j * dfft%nnr + 1:), sticks, &
                    dfft%nr3p, isgn, currsize, j / dfft%subbatchsize + 1)
                !
                i = cudaEventRecord(dfft%bevents(j / dfft%subbatchsize + 1), &
                                    dfft%bstreams(j / dfft%subbatchsize + 1))
                !
                i = cudaStreamWaitEvent(dfft%a2a_comp, &
                                        dfft%bevents(j / dfft%subbatchsize + 1), 0)
                !
                DO i = 0, currsize - 1
                    !
                    CALL env_cft_1z_gpu(aux_d(j * dfft%nnr + i * ncpx * nx3 + 1:), &
                                        sticks(me_p), n3, nx3, isgn, &
                                        f_d((j + i) * dfft%nnr + 1:), dfft%a2a_comp)
                    !
                END DO
                !
            END DO
            !
        END IF
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_many_cft3s_2d_gpu
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE env_fft_parallel_gpu
!----------------------------------------------------------------------------------------
#endif
