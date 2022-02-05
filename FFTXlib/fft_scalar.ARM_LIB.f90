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
#if defined(__ARM_LIB)
#if defined(_OPENMP)&&defined(__FFT_SCALAR_THREAD_SAFE)
! thread safety guard
#error ARM_LIB is not compatiable with __FFT_SCALAR_THREAD_SAFE
#endif
!
!----------------------------------------------------------------------------------------
!>
!! ARMlib routines (both 3d for serial execution and 1d+2d FFTs for parallel execution)
!!
!----------------------------------------------------------------------------------------
MODULE env_fft_scalar_arm
    !------------------------------------------------------------------------------------
    !
    USE class_io, ONLY: io
    !
    USE env_fft_param
    !
    USE, INTRINSIC :: ISO_C_BINDING
    ! iso_c_binding provides C_PTR, C_NULL_PTR, C_ASSOCIATED
    !
    USE env_fftw_interfaces
    !
    !------------------------------------------------------------------------------------
    !
    PRIVATE
    !
    PUBLIC :: env_cft_1z, env_cft_2xy, env_cfft3d
    !
    !------------------------------------------------------------------------------------
CONTAINS
    !------------------------------------------------------------------------------------
    !>
    !! FFT along "z"
    !!
    !! Driver routine for nsl 1d complex fft's of length nz
    !! ldz >= nz is the distance between sequences to be transformed
    !! (ldz>nz is used on some architectures to reduce memory conflicts)
    !! input  :  c(ldz*nsl)   (complex)
    !! output : cout(ldz*nsl) (complex - NOTA BENE: transform is not in-place)
    !! isign > 0 : forward (f(G)=>f(R)), isign <0 backward (f(R) => f(G))
    !! Up to "ndims" initializations (for different combinations of input
    !! parameters nz, nsl, ldz) are stored and re-used if available
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_cft_1z(c, nsl, nz, ldz, isign, cout)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: isign, nsl, nz, ldz
        !
        COMPLEX(DP) :: c(:), cout(:)
        !
        REAL(DP) :: tscale
        INTEGER :: i, err, idir, ip, void
        INTEGER, SAVE :: zdims(3, ndims) = -1
        INTEGER, SAVE :: icurrent = 1
        LOGICAL :: found
        !
        INTEGER :: tid
        !
        INTEGER, PARAMETER :: ltabl = 3 * nfftx + 100
        INTEGER :: INFO
        !
        COMPLEX(DP), SAVE :: fw_tablez(ltabl, ndims)
        COMPLEX(DP), SAVE :: bw_tablez(ltabl, ndims)
        !
        CHARACTER(LEN=80) :: sub_name = 'env_fft_scalar: env_cft_1z'
        !
        !--------------------------------------------------------------------------------
        !
        IF (nsl < 0) CALL io%error(sub_name, 'nsl out of range', nsl)
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
#if defined(__FFT_CLOCKS)
        CALL env_start_clock(sub_name)
#endif
        !
        IF (isign < 0) THEN
            tscale = 1.0_DP / nz
            !
            CALL ZFFT1MX(-1, tscale, .FALSE., nsl, ldz, c(1), 1, ldz, cout(1), 1, ldz, &
                         fw_tablez(1, ip), INFO)
            !
        ELSE IF (isign > 0) THEN
            !
            CALL ZFFT1MX(1, 1.0_DP, .FALSE., nsl, ldz, c(1), 1, ldz, cout(1), 1, ldz, &
                         bw_tablez(1, ip), INFO)
            !
        END IF

#if defined(__FFT_CLOCKS)
        CALL env_stop_clock(sub_name)
#endif
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
                found = (nz == zdims(1, ip))
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
            tscale = 1.0_DP / nz
            !
            CALL ZFFT1MX(0, tscale, .FALSE., nsl, nz, c(1), 1, ldz, cout(1), 1, ldz, &
                         fw_tablez(1, icurrent), INFO)
            !
            CALL ZFFT1MX(0, 1.0_DP, .FALSE., nsl, nz, c(1), 1, ldz, cout(1), 1, ldz, &
                         bw_tablez(1, icurrent), INFO)
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
    END SUBROUTINE env_cft_1z
    !------------------------------------------------------------------------------------
    !>
    !! FFT along "x" and "y" direction
    !!
    !! Driver routine for nzl 2d complex fft's of lengths nx and ny
    !! input : r(ldx*ldy)  complex, transform is in-place
    !! ldx >= nx, ldy >= ny are the physical dimensions of the equivalent
    !! 2d array: r2d(ldx, ldy) (x first dimension, y second dimension)
    !! (ldx>nx, ldy>ny used on some architectures to reduce memory conflicts)
    !! pl2ix(nx) (optional) is 1 for columns along y to be transformed
    !! isign > 0 : forward (f(G)=>f(R)), isign <0 backward (f(R) => f(G))
    !! Up to "ndims" initializations (for different combinations of input
    !! parameters nx, ny, nzl, ldx) are stored and re-used if available
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_cft_2xy(r, nzl, nx, ny, ldx, ldy, isign, pl2ix)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: isign, ldx, ldy, nx, ny, nzl
        INTEGER, OPTIONAL, INTENT(IN) :: pl2ix(:)
        !
        COMPLEX(DP) :: r(:)
        INTEGER :: i, k, j, err, idir, ip, kk, void
        REAL(DP) :: tscale
        INTEGER, SAVE :: icurrent = 1
        INTEGER, SAVE :: dims(4, ndims) = -1
        LOGICAL :: dofft(nfftx), found
        !
        INTEGER, PARAMETER :: ltabl = 3 * nfftx + 100
        INTEGER :: INFO
        !
        COMPLEX(DP), SAVE :: fw_tablex(ltabl, ndims), bw_tablex(ltabl, ndims)
        COMPLEX(DP), SAVE :: fw_tabley(ltabl, ndims), bw_tabley(ltabl, ndims)
        !
        CHARACTER(LEN=80) :: sub_name = 'env_fft_scalar: env_cft_2xy'
        !
        !--------------------------------------------------------------------------------
        dofft(1:nx) = .TRUE.
        !
        IF (PRESENT(pl2ix)) THEN
            !
            IF (SIZE(pl2ix) < nx) &
                CALL io%error(sub_name, 'Wrong dimension for arg no. 8', 1)
            !
            DO i = 1, nx
                IF (pl2ix(i) < 1) dofft(i) = .FALSE.
            END DO
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
#if defined(__FFT_CLOCKS)
        CALL env_start_clock(sub_name)
#endif
        !
#if defined(_OPENMP)
        IF (isign < 0) THEN
            tscale = 1.0_DP / (nx * ny)
            !
            DO k = 1, nzl
                kk = 1 + (k - 1) * ldx * ldy
                !
                CALL ZFFT1MX(-1, tscale, .TRUE., ny, nx, r(kk), 1, ldx, r(kk), 1, ldx, &
                             fw_tablex(1, ip), INFO)
                !
                CALL ZFFT1MX(-1, 1.0_DP, .TRUE., nx, ny, r(kk), ldx, 1, r(kk), ldx, 1, &
                             fw_tabley(1, ip), INFO)
                !
            END DO
            !
        ELSE IF (isign > 0) THEN
            !
            DO k = 1, nzl
                kk = 1 + (k - 1) * ldx * ldy
                !
                CALL ZFFT1MX(1, 1.0_DP, .TRUE., nx, ny, r(kk), ldx, 1, r(kk), ldx, 1, &
                             bw_tabley(1, ip), INFO)
                !
                CALL ZFFT1MX(1, 1.0_DP, .TRUE., ny, nx, r(kk), 1, ldx, r(kk), 1, ldx, &
                             bw_tablex(1, ip), INFO)
                !
            END DO
            !
        END IF
#else
        !
        IF (isign < 0) THEN
            tscale = 1.0_DP / (nx * ny)
            !
            DO k = 1, nzl
                kk = 1 + (k - 1) * ldx * ldy
                !
                CALL ZFFT1MX(-1, tscale, .TRUE., ny, nx, r(kk), 1, ldx, r(kk), 1, ldx, &
                             fw_tablex(1, ip), INFO)
                !
                DO i = 1, nx
                    !
                    IF (dofft(i)) THEN
                        kk = i + (k - 1) * ldx * ldy
                        !
                        CALL ZFFT1MX(-1, 1.0_DP, .TRUE., 1, ny, r(kk), ldx, 1, r(kk), &
                                     ldx, 1, fw_tabley(1, ip), INFO)
                        !
                    END IF
                    !
                END DO
                !
            END DO
            !
        ELSE IF (isign > 0) THEN
            !
            DO k = 1, nzl
                !
                DO i = 1, nx
                    !
                    IF (dofft(i)) THEN
                        kk = i + (k - 1) * ldx * ldy
                        !
                        CALL ZFFT1MX(1, 1.0_DP, .TRUE., 1, ny, r(kk), ldx, 1, r(kk), &
                                     ldx, 1, bw_tabley(1, ip), INFO)
                        !
                    END IF
                    !
                END DO
                !
                kk = 1 + (k - 1) * ldx * ldy
                !
                CALL ZFFT1MX(1, 1.0_DP, .TRUE., ny, nx, r(kk), 1, ldx, r(kk), 1, ldx, &
                             bw_tablex(1, ip), INFO)
                !
            END DO
            !
        END IF
#endif
        !
#if defined(__FFT_CLOCKS)
        CALL env_stop_clock(sub_name)
#endif
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
                !
                found = (ny == dims(1, ip)) .AND. (nx == dims(3, ip))
                found = found .AND. (ldx == dims(2, ip)) .AND. (nzl == dims(4, ip))
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
            tscale = 1.0_DP / (nx * ny)
            !
#if defined(_OPENMP)
            !
            CALL ZFFT1MX(0, 1.0_DP, .TRUE., nx, ny, r(1), ldx, 1, r(1), ldx, 1, &
                         fw_tabley(1, icurrent), INFO)
            !
            CALL ZFFT1MX(0, 1.0_DP, .TRUE., nx, ny, r(1), ldx, 1, r(1), ldx, 1, &
                         bw_tabley(1, icurrent), INFO)
            !
            CALL ZFFT1MX(0, tscale, .TRUE., ny, nx, r(1), 1, ldx, r(1), 1, ldx, &
                         fw_tablex(1, icurrent), INFO)
            !
            CALL ZFFT1MX(0, 1.0_DP, .TRUE., ny, nx, r(1), 1, ldx, r(1), 1, ldx, &
                         bw_tablex(1, icurrent), INFO)
            !
#else
            !
            CALL ZFFT1MX(0, 1.0_DP, .TRUE., 1, ny, r(1), ldx, 1, r(1), ldx, 1, &
                         fw_tabley(1, icurrent), INFO)
            !
            CALL ZFFT1MX(0, 1.0_DP, .TRUE., 1, ny, r(1), ldx, 1, r(1), ldx, 1, &
                         bw_tabley(1, icurrent), INFO)
            !
            CALL ZFFT1MX(0, tscale, .TRUE., ny, nx, r(1), 1, ldx, r(1), 1, ldx, &
                         fw_tablex(1, icurrent), INFO)
            !
            CALL ZFFT1MX(0, 1.0_DP, .TRUE., ny, nx, r(1), 1, ldx, r(1), 1, ldx, &
                         bw_tablex(1, icurrent), INFO)
            !
#endif
            !
            dims(1, icurrent) = ny
            dims(2, icurrent) = ldx
            dims(3, icurrent) = nx
            dims(4, icurrent) = nzl
            !
            ip = icurrent
            icurrent = MOD(icurrent, ndims) + 1
            !
            !----------------------------------------------------------------------------
        END SUBROUTINE init_plan
        !--------------------------------------------------------------------------------
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_cft_2xy
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
    SUBROUTINE env_cfft3d(f, nx, ny, nz, ldx, ldy, ldz, howmany, isign)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: nx, ny, nz, ldx, ldy, ldz, howmany, isign
        !
        COMPLEX(DP) :: f(:)
        INTEGER :: i, k, j, err, idir, ip
        REAL(DP) :: tscale
        INTEGER, SAVE :: icurrent = 1
        INTEGER, SAVE :: dims(3, ndims) = -1
        !
#if defined(__ARM_LIB_BUG)
        INTEGER, PARAMETER :: ltabl = 4 * nfftx + 300
        INTEGER :: INFO
        !
        COMPLEX(DP), SAVE :: fw_table(ltabl, ndims)
        COMPLEX(DP), SAVE :: bw_table(ltabl, ndims)
#else
        TYPE(C_PTR), SAVE :: fw_plan(ndims) = C_NULL_PTR
        TYPE(C_PTR), SAVE :: bw_plan(ndims) = C_NULL_PTR
#endif
        !
        CHARACTER(LEN=80) :: sub_name = 'fft scalar: env_cfft3d'
        !
        !--------------------------------------------------------------------------------
        !
        IF (nx < 1) CALL io%error(sub_name, 'nx is less than 1', 1)
        !
        IF (ny < 1) CALL io%error(sub_name, 'ny is less than 1', 1)
        !
        IF (nz < 1) CALL io%error(sub_name, 'nz is less than 1', 1)
        !
        IF (howmany /= 1) &
            CALL io%error(sub_name, &
                          'howmany different from 1 not yet implemented for ARM', 1)
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
#if defined(__ARM_LIB_BUG)
        IF (isign < 0) THEN
            tscale = 1.0_DP / DBLE(nx * ny * nz)
            !
            CALL ZFFT3DY(-1, tscale, .TRUE., nx, ny, nz, f(1), 1, ldx, ldx * ldy, &
                         f(1), 1, ldx, ldx * ldy, fw_table(1, ip), ltabl, INFO)
            !
        ELSE IF (isign > 0) THEN
            !
            CALL ZFFT3DY(1, 1.0_DP, .TRUE., nx, ny, nz, f(1), 1, ldx, ldx * ldy, &
                         f(1), 1, ldx, ldx * ldy, bw_table(1, ip), ltabl, INFO)
            !
        END IF
#else
        !
        IF (isign < 0) THEN
            !
            CALL env_fftw_inplace_drv_3d(fw_plan(ip), 1, f(1), 1, 1)
            !
            tscale = 1.0_DP / DBLE(nx * ny * nz)
            !
            CALL ZDSCAL(nx * ny * nz, tscale, f(1), 1)
            !
        ELSE IF (isign > 0) THEN
            CALL env_fftw_inplace_drv_3d(bw_plan(ip), 1, f(1), 1, 1)
        END IF
#endif
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
                    (nz == dims(3, i))) THEN
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
            IF (nx /= ldx .OR. ny /= ldy .OR. nz /= ldz) &
                CALL io%error(sub_name, 'Not implemented', 1)
            !
#if defined(__ARM_LIB_BUG)
            tscale = 1.0_DP / DBLE(nx * ny * nz)
            !
            CALL ZFFT3DY(0, tscale, .TRUE., nx, ny, nz, f(1), 1, ldx, ldx * ldy, &
                         f(1), 1, ldx, ldx * ldy, fw_table(1, icurrent), ltabl, INFO)
            !
            CALL ZFFT3DY(0, 1.0_DP, .TRUE., nx, ny, nz, f(1), 1, ldx, ldx * ldy, &
                         f(1), 1, ldx, ldx * ldy, bw_table(1, icurrent), ltabl, INFO)
            !
#else
            !
            IF (C_ASSOCIATED(fw_plan(icurrent))) &
                CALL env_destroy_plan_3d(fw_plan(icurrent))
            !
            IF (C_ASSOCIATED(bw_plan(icurrent))) &
                CALL env_destroy_plan_3d(bw_plan(icurrent))
            !
            idir = -1
            !
            CALL env_create_plan_3d(fw_plan(icurrent), nx, ny, nz, idir)
            !
            idir = 1
            !
            CALL env_create_plan_3d(bw_plan(icurrent), nx, ny, nz, idir)
#endif
            !
            dims(1, icurrent) = nx
            dims(2, icurrent) = ny
            dims(3, icurrent) = nz
            ip = icurrent
            icurrent = MOD(icurrent, ndims) + 1
            !
            !----------------------------------------------------------------------------
        END SUBROUTINE init_plan
        !--------------------------------------------------------------------------------
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_cfft3d
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE env_fft_scalar_arm
!----------------------------------------------------------------------------------------
#endif
