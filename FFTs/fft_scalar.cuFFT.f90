!----------------------------------------------------------------------------------------
!
! Copyright (C) 2021 ENVIRON (www.quantum-environ.org)
! Copyright (C) Quantum ESPRESSO (www.quantum-espresso.org)
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
! Authors: Carlo Cavazzoni, modified by P. Giannozzi, contributions
!          by Martin Hilgemans, Guido Roma, Pascal Thibaudeau, Stephane Lefranc,
!          Nicolas Lacorne, Filippo Spiga, Nicola Varini, Jason Wood
!          Modified by Edan Bainglass
!
!----------------------------------------------------------------------------------------
!
#if defined(__CUDA)
!
#define __CUFFT_ALL_XY_PLANES
!
!----------------------------------------------------------------------------------------
!>
!! Cuda FFT routines
!!
!----------------------------------------------------------------------------------------
MODULE env_fft_scalar_cuFFT
    !------------------------------------------------------------------------------------
    !
    USE cudafor
    !
    USE cufft
    !
    USE env_fft_param
    !
    USE, INTRINSIC :: ISO_C_BINDING
    ! iso_c_binding provides C_PTR, C_NULL_PTR, C_ASSOCIATED
    !
#if defined(TRACK_FLOPS)
    USE flops_tracker, ONLY: fft_ops
#endif
    !
    !------------------------------------------------------------------------------------
    !
    PRIVATE
    !
    PUBLIC :: env_cft_1z_gpu, env_cft_2xy_gpu, env_cfft3d_gpu
    !
    !------------------------------------------------------------------------------------
CONTAINS
    !------------------------------------------------------------------------------------
    !>
    !! FFT along "x" and "y" direction
    !!
    !! Driver routine for nsl 1d complex fft's of length nz
    !! ldz >= nz is the distance between sequences to be transformed
    !! (ldz>nz is used on some architectures to reduce memory conflicts)
    !! input  :  c_d(ldz*nsl)   (complex)
    !! ### GPU VERION IN PLACE!!!
    !! ### output : cout_d(ldz*nsl) (complex - NOTA BENE: transform is not in-place!)
    !! isign > 0 : forward (f(G)=>f(R)), isign <0 backward (f(R) => f(G))
    !! Up to "ndims" initializations (for different combinations of input
    !! parameters nz, nsl, ldz) are stored and re-used if available
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_cft_1z_gpu(c_d, nsl, nz, ldz, isign, cout_d, stream, in_place)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: isign, nsl, nz, ldz
        LOGICAL, INTENT(IN), OPTIONAL :: in_place
        !
        INTEGER(kind=cuda_stream_kind) :: stream
        !
        COMPLEX(DP), DEVICE :: c_d(:), cout_d(:)
        !
        REAL(DP) :: tscale
        INTEGER :: i, err, idir, ip, void, istat
        !
#if defined(TRACK_FLOPS)
        REAL(DP), SAVE :: zflops(ndims) = 0.D0
#endif
        !
        INTEGER, SAVE :: zdims(3, ndims) = -1
        INTEGER, SAVE :: icurrent = 1
        LOGICAL :: found
        LOGICAL :: is_inplace
        !
        INTEGER :: tid
        !
        INTEGER, SAVE :: cufft_planz(ndims) = 0 ! contains FFT factors ( PLAN )
        !
        CHARACTER(LEN=80) :: sub_name = 'env_fft_scalar: env_cft_1z_gpu'
        !
        !--------------------------------------------------------------------------------
        !
        IF (nsl < 0) CALL env_errore(sub_name, 'nsl out of range', nsl)
        !
        !--------------------------------------------------------------------------------
        ! Here initialize table only if necessary
        !
        CALL lookup()
        !
        IF (.NOT. found) CALL init_plan()
        ! no table exist for these parameters; initialize a new one
        !
        !--------------------------------------------------------------------------------
        ! Now perform the FFTs using machine specific drivers
        !
        IF (PRESENT(in_place)) THEN
            is_inplace = in_place
        ELSE
            is_inplace = .FALSE.
        END IF
        !
        istat = cufftSetStream(cufft_planz(ip), stream)
        !
#if defined(__FFT_CLOCKS)
        CALL env_start_clock(sub_name)
#endif
        !
        IF (isign < 0) THEN
            istat = cufftExecZ2Z(cufft_planz(ip), c_d(1), c_d(1), CUFFT_FORWARD)
            tscale = 1.0_DP / nz
            !
            IF (is_inplace) THEN
                !
                !$cuf kernel do(1) <<<*,*,0,stream>>>
                DO i = 1, ldz * nsl
                    c_d(i) = c_d(i) * tscale
                END DO
                !
            ELSE
                !
                !$cuf kernel do(1) <<<*,*,0,stream>>>
                DO i = 1, ldz * nsl
                    cout_d(i) = c_d(i) * tscale
                END DO
                !
            END IF
            !
        ELSE IF (isign > 0) THEN
            !
            IF (is_inplace) THEN
                istat = cufftExecZ2Z(cufft_planz(ip), c_d(1), c_d(1), CUFFT_INVERSE)
            ELSE
                istat = cufftExecZ2Z(cufft_planz(ip), c_d(1), cout_d(1), CUFFT_INVERSE)
            END IF
            !
        END IF
        !
#if defined(__FFT_CLOCKS)
        CALL env_stop_clock(sub_name)
#endif
        !
#if defined(TRACK_FLOPS)
        fft_ops = fft_ops + zflops(ip)
#endif
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    CONTAINS
        !--------------------------------------------------------------------------------
        !>
        !! Check if there is already a table initialized
        !! for this combination of parameters
        !!
        !--------------------------------------------------------------------------------
        SUBROUTINE lookup()
            !----------------------------------------------------------------------------
            !
            DO ip = 1, ndims
                found = (nz == zdims(1, ip)) .AND. &
                        (nsl == zdims(2, ip)) .AND. &
                        (ldz == zdims(3, ip))
                !
                IF (found) EXIT
            END DO
            !
            !----------------------------------------------------------------------------
        END SUBROUTINE lookup
        !--------------------------------------------------------------------------------
        !>
        !!
        !--------------------------------------------------------------------------------
        SUBROUTINE init_plan()
            !----------------------------------------------------------------------------
            !
            IMPLICIT NONE
            !
            INTEGER, PARAMETER :: RANK = 1
            INTEGER :: FFT_DIM(RANK), DATA_DIM(RANK)
            INTEGER :: STRIDE, DIST, BATCH
            !
            !----------------------------------------------------------------------------
            !
            FFT_DIM(1) = nz
            DATA_DIM(1) = ldz
            STRIDE = 1
            DIST = ldz
            BATCH = nsl
            !
            IF (cufft_planz(icurrent) /= 0) istat = cufftDestroy(cufft_planz(icurrent))
            !
            istat = cufftPlanMany(cufft_planz(icurrent), RANK, FFT_DIM, &
                                  DATA_DIM, STRIDE, DIST, &
                                  DATA_DIM, STRIDE, DIST, &
                                  CUFFT_Z2Z, BATCH)
            !
#if defined(__CUDA_DEBUG)
            PRINT *, "INIT CUFFT Z PLAN: ", nz, "x", nsl, "x", ldz
#endif
            !
#if defined(TRACK_FLOPS)
            zflops(icurrent) = 5.0D0 * REAL(nz) * LOG(REAL(nz)) / LOG(2.D0)
#endif
            !
            zdims(1, icurrent) = nz
            zdims(2, icurrent) = nsl
            zdims(3, icurrent) = ldz
            ip = icurrent
            icurrent = MOD(icurrent, ndims) + 1
            !
            !----------------------------------------------------------------------------
        END SUBROUTINE init_plan
        !--------------------------------------------------------------------------------
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_cft_1z_gpu
    !------------------------------------------------------------------------------------
    !>
    !! 3D scalar FFTs
    !!
    !! Driver routine for nzl 2d complex fft's of lengths nx and ny
    !! input : r_d(ldx*ldy)  complex, transform is in-place
    !! ldx >= nx, ldy >= ny are the physical dimensions of the equivalent
    !! 2d array: r2d(ldx, ldy) (x first dimension, y second dimension)
    !! (ldx>nx, ldy>ny used on some architectures to reduce memory conflicts)
    !! pl2ix(nx) (optional) is 1 for columns along y to be transformed
    !! isign > 0 : forward (f(G)=>f(R)), isign <0 backward (f(R) => f(G))
    !! Up to "ndims" initializations (for different combinations of input
    !! parameters nx, ny, nzl, ldx) are stored and re-used if available
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_cft_2xy_gpu(r_d, temp_d, nzl, nx, ny, ldx, ldy, isign, stream, pl2ix)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: isign, ldx, ldy, nx, ny, nzl
        INTEGER, OPTIONAL, INTENT(IN) :: pl2ix(:)
        !
        INTEGER(kind=cuda_stream_kind) :: stream
        !
        COMPLEX(DP), DEVICE :: r_d(ldx, ldy, nzl), temp_d(ldy, nzl, ldx)
        INTEGER :: i, k, j, err, idir, ip, kk, void, istat
        REAL(DP) :: tscale
        INTEGER, SAVE :: icurrent = 1
        INTEGER, SAVE :: dims(6, ndims) = -1
        LOGICAL :: dofft(nfftx), found
        !
#if defined(TRACK_FLOPS)
        REAL(DP), SAVE :: xyflops(ndims) = 0.D0
#endif
        !
#if defined(__CUFFT_ALL_XY_PLANES)
        INTEGER, SAVE :: cufft_plan_2d(ndims) = 0
#else
        INTEGER, SAVE :: cufft_plan_x(ndims) = 0
        INTEGER, SAVE :: cufft_plan_y(2, ndims) = 0
#endif
        INTEGER :: batch_1, batch_2
        !
        CHARACTER(LEN=80) :: sub_name = 'env_fft_scalar: env_cft_2xy_gpu'
        !
        !--------------------------------------------------------------------------------
        !
        dofft(1:nx) = .TRUE.
        batch_1 = nx
        batch_2 = 0
        !
        IF (PRESENT(pl2ix)) THEN
            !
            IF (SIZE(pl2ix) < nx) &
                CALL env_errore(sub_name, 'Wrong dimension for arg no. 8', 1)
            !
            DO i = 1, nx
                !
                IF (pl2ix(i) < 1) dofft(i) = .FALSE.
                !
            END DO
            !
            i = 1
            !
            DO WHILE (pl2ix(i) >= 1 .AND. i <= nx); i = i + 1; END DO
            !
            batch_1 = i - 1
            !
            DO WHILE (pl2ix(i) < 1 .AND. i <= nx); i = i + 1; END DO
            !
            batch_2 = nx - i + 1
            !
#if 0
            DO WHILE (i <= nx)
                !
                DO WHILE (pl2ix(i) < 1 .AND. i <= nx); i = i + 1; END DO
                !
                batch_start = i
                !
                DO WHILE (pl2ix(i) >= 1 .AND. i <= nx); i = i + 1; END DO
                !
                batch_end = i - 1
                batch_count = batch_end - batch_start + 1
            END DO
#endif
            !
        END IF
        !
        !--------------------------------------------------------------------------------
        ! Here initialize table only if necessary
        !
        CALL lookup()
        !
        IF (.NOT. found) CALL init_plan()
        ! no table exist for these parameters; initialize a new one
        !
        !--------------------------------------------------------------------------------
        ! Now perform the FFTs using machine specific drivers
        !
#if defined(__CUFFT_ALL_XY_PLANES)
        istat = cufftSetStream(cufft_plan_2d(ip), stream)
#else
        istat = cufftSetStream(cufft_plan_x(ip), stream)
        istat = cufftSetStream(cufft_plan_y(1, ip), stream)
        istat = cufftSetStream(cufft_plan_y(2, ip), stream)
#endif
        !
#if defined(__FFT_CLOCKS)
        CALL env_start_clock(sub_name)
#endif
        !
        IF (isign < 0) THEN
            tscale = 1.0_DP / (nx * ny)
            !
#if defined(__CUFFT_ALL_XY_PLANES)
            !
            istat = cufftExecZ2Z(cufft_plan_2d(ip), r_d(1, 1, 1), r_d(1, 1, 1), &
                                 CUFFT_FORWARD)
            !
            !$cuf kernel do(3) <<<*,(16,16,1), 0, stream>>>
            DO k = 1, nzl
                !
                DO j = 1, ldy
                    !
                    DO i = 1, ldx
                        r_d(i, j, k) = r_d(i, j, k) * tscale
                    END DO
                    !
                END DO
                !
            END DO
#else
            !
            istat = cufftExecZ2Z(cufft_plan_x(ip), r_d(1, 1, 1), r_d(1, 1, 1), &
                                 CUFFT_FORWARD)
            !
            IF (istat) PRINT *, 'error in fftxy fftx istat = ', istat
            !
            !$cuf kernel do(3) <<<*,(16,16,1), 0, stream>>>
            DO k = 1, nzl
                !
                DO i = 1, ldx
                    !
                    DO j = 1, ldy
                        temp_d(j, k, i) = r_d(i, j, k)
                    END DO
                    !
                END DO
                !
            END DO
            !
            IF (batch_1 > 0) THEN
                !
                istat = cufftExecZ2Z(cufft_plan_y(1, ip), temp_d(1, 1, 1), &
                                     temp_d(1, 1, 1), CUFFT_FORWARD)
                !
                IF (istat) PRINT *, 'error in fftxy ffty batch_1 istat = ', istat
                !
            END IF
            !
            IF (batch_2 > 0) THEN
                !
                istat = cufftExecZ2Z(cufft_plan_y(2, ip), &
                                     temp_d(1, 1, nx - batch_2 + 1), &
                                     temp_d(1, 1, nx - batch_2 + 1), CUFFT_FORWARD)
                !
                IF (istat) PRINT *, 'error in fftxy ffty batch_2 istat = ', istat
                !
            END IF
            !
            !$cuf kernel do(3) <<<*,(16,16,1), 0, stream>>>
            DO k = 1, nzl
                !
                DO j = 1, ldy
                    !
                    DO i = 1, ldx
                        r_d(i, j, k) = temp_d(j, k, i) * tscale
                    END DO
                    !
                END DO
                !
            END DO
#endif
            !
        ELSE IF (isign > 0) THEN
            !
#if defined(__CUFFT_ALL_XY_PLANES)
            istat = cufftExecZ2Z(cufft_plan_2d(ip), r_d(1, 1, 1), r_d(1, 1, 1), &
                                 CUFFT_INVERSE)
#else
            !
            !$cuf kernel do(3) <<<*,(16,16,1), 0, stream>>>
            DO k = 1, nzl
                !
                DO i = 1, ldx
                    !
                    DO j = 1, ldy
                        temp_d(j, k, i) = r_d(i, j, k)
                    END DO
                    !
                END DO
                !
            END DO
            !
            IF (batch_1 > 0) THEN
                !
                istat = cufftExecZ2Z(cufft_plan_y(1, ip), temp_d(1, 1, 1), &
                                     temp_d(1, 1, 1), CUFFT_INVERSE)
                !
                IF (istat) PRINT *, 'error in fftxy ffty batch_1 istat = ', istat
                !
            END IF
            !
            IF (batch_2 > 0) THEN
                !
                istat = cufftExecZ2Z(cufft_plan_y(2, ip), &
                                     temp_d(1, 1, nx - batch_2 + 1), &
                                     temp_d(1, 1, nx - batch_2 + 1), CUFFT_INVERSE)
                !
                IF (istat) PRINT *, 'error in fftxy ffty batch_2 istat = ', istat
                !
            END IF
            !
            !$cuf kernel do(3) <<<*,(16,16,1), 0, stream>>>
            DO k = 1, nzl
                !
                DO j = 1, ldy
                    !
                    DO i = 1, ldx
                        r_d(i, j, k) = temp_d(j, k, i)
                    END DO
                    !
                END DO
                !
            END DO
            !
            istat = cufftExecZ2Z(cufft_plan_x(ip), r_d(1, 1, 1), r_d(1, 1, 1), &
                                 CUFFT_INVERSE)
            !
            IF (istat) PRINT *, 'error in fftxy fftx istat = ', istat
#endif
            !
        END IF
        !
#if defined(__FFT_CLOCKS)
        CALL env_stop_clock(sub_name)
#endif
        !
#if defined(TRACK_FLOPS)
        fft_ops = fft_ops + xyflops(ip)
#endif
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    CONTAINS
        !--------------------------------------------------------------------------------
        !>
        !! Check if there is already a table initialized
        !! for this combination of parameters
        !!
        !--------------------------------------------------------------------------------
        SUBROUTINE lookup()
            !----------------------------------------------------------------------------
            !
            DO ip = 1, ndims
                found = (ny == dims(1, ip)) .AND. (nx == dims(3, ip))
                found = found .AND. (ldx == dims(2, ip)) .AND. (nzl == dims(4, ip))
                !
#if ! defined(__CUFFT_ALL_XY_PLANES)
                found = found .AND. (batch_1 == dims(5, ip)) .AND. &
                        (batch_2 == dims(6, ip))
#endif
                !
                IF (found) EXIT
                !
            END DO
            !
            !----------------------------------------------------------------------------
        END SUBROUTINE lookup
        !--------------------------------------------------------------------------------
        !>
        !!
        !--------------------------------------------------------------------------------
        SUBROUTINE init_plan()
            !----------------------------------------------------------------------------
            !
            IMPLICIT NONE
            !
#if defined(__CUFFT_ALL_XY_PLANES)
            INTEGER, PARAMETER :: RANK = 2
            INTEGER :: FFT_DIM(RANK), DATA_DIM(RANK)
            INTEGER :: STRIDE, DIST, BATCH
            !
            !----------------------------------------------------------------------------
            !
            FFT_DIM(1) = ny
            FFT_DIM(2) = nx
            DATA_DIM(1) = ldy
            DATA_DIM(2) = ldx
            STRIDE = 1
            DIST = ldx * ldy
            BATCH = nzl
            !
            IF (cufft_plan_2d(icurrent) /= 0) &
                istat = cufftDestroy(cufft_plan_2d(icurrent))
            !
            istat = cufftPlanMany(cufft_plan_2d(icurrent), RANK, FFT_DIM, &
                                  DATA_DIM, STRIDE, DIST, &
                                  DATA_DIM, STRIDE, DIST, &
                                  CUFFT_Z2Z, BATCH)
            !
#if defined(__CUDA_DEBUG)
            PRINT *, 'INIT CUFFT ALL_XY PLAN: ', nx, 'x', ny, 'x', nzl, 'ldx:', ldx, &
                'batch:', batch_1, batch_2
#endif
            !
#else
            INTEGER, PARAMETER :: RANK = 1
            INTEGER :: FFT_DIM_X(RANK), DATA_DIM_X(RANK)
            INTEGER :: FFT_DIM_Y(RANK), DATA_DIM_Y(RANK)
            INTEGER :: STRIDE_X, STRIDE_Y, DIST_X, DIST_Y, BATCH_X, BATCH_Y1, BATCH_Y2
            !
            !----------------------------------------------------------------------------
            !
            FFT_DIM_X(1) = nx
            DATA_DIM_X(1) = ldx
            STRIDE_X = 1
            DIST_X = ldx
            BATCH_X = ny * nzl
            !
            FFT_DIM_Y(1) = ny
            DATA_DIM_Y(1) = ldy
            STRIDE_Y = 1
            DIST_Y = ldy
            BATCH_Y1 = nzl * BATCH_1
            BATCH_Y2 = nzl * BATCH_2
            !
            IF (cufft_plan_x(icurrent) /= 0) &
                istat = cufftDestroy(cufft_plan_x(icurrent))
            !
            IF (cufft_plan_y(1, icurrent) /= 0) &
                istat = cufftDestroy(cufft_plan_y(1, icurrent))
            !
            IF (cufft_plan_y(2, icurrent) /= 0) &
                istat = cufftDestroy(cufft_plan_y(2, icurrent))
            !
#if defined(__CUDA_DEBUG)
            PRINT *, 'INIT CUFFT XY PLAN: ', nx, 'x', ny, 'x', nzl, 'ldx:', ldx, &
                'batch:', batch_1, batch_2
#endif
            !
            istat = cufftPlanMany(cufft_plan_x(icurrent), RANK, FFT_DIM_X, &
                                  DATA_DIM_X, STRIDE_X, DIST_X, &
                                  DATA_DIM_X, STRIDE_X, DIST_X, &
                                  CUFFT_Z2Z, BATCH_X)
            !
            istat = cufftPlanMany(cufft_plan_y(1, icurrent), RANK, FFT_DIM_Y, &
                                  DATA_DIM_Y, STRIDE_Y, DIST_Y, &
                                  DATA_DIM_Y, STRIDE_Y, DIST_Y, &
                                  CUFFT_Z2Z, BATCH_Y1)
            !
            istat = cufftPlanMany(cufft_plan_y(2, icurrent), RANK, FFT_DIM_Y, &
                                  DATA_DIM_Y, STRIDE_Y, DIST_Y, &
                                  DATA_DIM_Y, STRIDE_Y, DIST_Y, &
                                  CUFFT_Z2Z, BATCH_Y2)
#endif
            !
#if defined(TRACK_FLOPS)
            xyflops(icurrent) = &
                REAL(ny * nzl) * 5.0D0 * REAL(nx) * LOG(REAL(nx)) / LOG(2.D0) + &
                REAL(nzl * BATCH_1 + nzl * BATCH_2) * 5.0D0 * REAL(ny) * LOG(REAL(ny)) &
                / LOG(2.D0)
#endif
            !
            dims(1, icurrent) = ny
            dims(2, icurrent) = ldx
            dims(3, icurrent) = nx
            dims(4, icurrent) = nzl
            dims(5, icurrent) = BATCH_1
            dims(6, icurrent) = BATCH_2
            ip = icurrent
            icurrent = MOD(icurrent, ndims) + 1
            !
            !----------------------------------------------------------------------------
        END SUBROUTINE init_plan
        !--------------------------------------------------------------------------------
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_cft_2xy_gpu
    !------------------------------------------------------------------------------------
    !>
    !! 3D scalar FFTs
    !!
    !! Driver routine for 3d complex fft of lengths nx, ny, nz
    !! input  :  f(ldx*ldy*ldz)  complex, transform is in-place
    !! ldx >= nx, ldy >= ny, ldz >= nz are the physical dimensions
    !! of the equivalent 3d array: f3d(ldx,ldy,ldz)
    !! (ldx>nx, ldy>ny, ldz>nz may be used on some architectures
    !! to reduce memory conflicts - not implemented for FFTW)
    !! isign > 0 : f(G) => f(R)   ; isign < 0 : f(R) => f(G)
    !!
    !! Up to "ndims" initializations (for different combinations of input
    !! parameters nx, ny, nz) are stored and re-used if available
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_cfft3d_gpu(f_d, nx, ny, nz, ldx, ldy, ldz, howmany, isign, stream)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: nx, ny, nz, ldx, ldy, ldz, howmany, isign
        !
        COMPLEX(DP), device :: f_d(:)
        INTEGER(kind=cuda_stream_kind) :: stream
        INTEGER :: i, k, j, err, idir, ip, istat
        REAL(DP) :: tscale
        INTEGER, SAVE :: icurrent = 1
        INTEGER, SAVE :: dims(4, ndims) = -1
        !
        INTEGER, SAVE :: cufft_plan_3d(ndims) = 0
        !
        CHARACTER(LEN=80) :: sub_name = 'env_fft_scalar: env_cfft3d_gpu'
        !
        !--------------------------------------------------------------------------------
        !
        IF (nx < 1) CALL env_errore(sub_name, ' nx is less than 1 ', 1)
        !
        IF (ny < 1) CALL env_errore(sub_name, ' ny is less than 1 ', 1)
        !
        IF (nz < 1) CALL env_errore(sub_name, ' nz is less than 1 ', 1)
        !
        IF (nx /= ldx .OR. ny /= ldy .OR. nz /= ldz) &
            CALL env_errore(sub_name, 'Leading dimensions must match data dimension', 1)
        !
        !--------------------------------------------------------------------------------
        ! Here initialize table only if necessary
        !
        CALL lookup()
        !
        IF (ip == -1) CALL init_plan()
        ! no table exist for these parameters; initialize a new one
        !
        !--------------------------------------------------------------------------------
        ! Now perform the 3D FFT using the machine specific driver
        !
        istat = cufftSetStream(cufft_plan_3d(ip), stream)
        !
        IF (isign < 0) THEN
            istat = cufftExecZ2Z(cufft_plan_3d(ip), f_d(1), f_d(1), CUFFT_FORWARD)
            !
            tscale = 1.0_DP / DBLE(nx * ny * nz)
            !
            !$cuf kernel do(1) <<<*,(16,16,1),0,stream>>>
            DO i = 1, ldx * ldy * ldz * howmany
                f_d(i) = f_d(i) * tscale
            END DO
            !
        ELSE IF (isign > 0) THEN
            istat = cufftExecZ2Z(cufft_plan_3d(ip), f_d(1), f_d(1), CUFFT_INVERSE)
        END IF
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    CONTAINS
        !--------------------------------------------------------------------------------
        !>
        !! Check if there is already a table initialized
        !! for this combination of parameters
        !!
        !--------------------------------------------------------------------------------
        SUBROUTINE lookup()
            !----------------------------------------------------------------------------
            !
            ip = -1
            !
            DO i = 1, ndims
                !
                IF ((nx == dims(1, i)) .AND. &
                    (ny == dims(2, i)) .AND. &
                    (nz == dims(3, i)) .AND. &
                    (howmany == dims(4, i))) THEN
                    !
                    ip = i
                    !
                    EXIT
                    !
                END IF
                !
            END DO
            !
            !----------------------------------------------------------------------------
        END SUBROUTINE lookup
        !--------------------------------------------------------------------------------
        !>
        !!
        !--------------------------------------------------------------------------------
        SUBROUTINE init_plan()
            !----------------------------------------------------------------------------
            !
            IMPLICIT NONE
            !
            INTEGER, PARAMETER :: RANK = 3
            INTEGER :: FFT_DIM(RANK), DATA_DIM(RANK)
            INTEGER :: STRIDE, DIST, BATCH
            !
            !----------------------------------------------------------------------------
            !
            FFT_DIM(1) = nz
            FFT_DIM(2) = ny
            FFT_DIM(3) = nx
            DATA_DIM(1) = ldz
            DATA_DIM(2) = ldy
            DATA_DIM(3) = ldx
            STRIDE = 1
            DIST = ldx * ldy * ldz
            BATCH = howmany
            !
            IF (cufft_plan_3d(icurrent) /= 0) &
                istat = cufftDestroy(cufft_plan_3d(icurrent))
            !
            istat = cufftPlanMany(cufft_plan_3d(icurrent), RANK, FFT_DIM, &
                                  DATA_DIM, STRIDE, DIST, &
                                  DATA_DIM, STRIDE, DIST, &
                                  CUFFT_Z2Z, BATCH)
            !
            dims(1, icurrent) = nx
            dims(2, icurrent) = ny
            dims(3, icurrent) = nz
            dims(4, icurrent) = howmany
            ip = icurrent
            icurrent = MOD(icurrent, ndims) + 1
            !
            !----------------------------------------------------------------------------
        END SUBROUTINE init_plan
        !--------------------------------------------------------------------------------
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_cfft3d_gpu
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE env_fft_scalar_cuFFT
!----------------------------------------------------------------------------------------
#endif
