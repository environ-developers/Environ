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
#if defined(__SX6)
#if defined(_OPENMP)&&defined(__FFT_SCALAR_THREAD_SAFE)
! thread safety guard
#error SX6 is not compatiable with __FFT_SCALAR_THREAD_SAFE
#endif
!
!----------------------------------------------------------------------------------------
!>
!! Uses legacy NEC ASL libraries (3d only, no parallel execution)
!!
!----------------------------------------------------------------------------------------
MODULE env_fft_scalar_sx6
    !------------------------------------------------------------------------------------
    !
    USE class_io, ONLY: io
    !
    USE env_fft_param
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
        INTEGER, PARAMETER :: ltabl = 2 * nfftx + 64
        ! dimension of the tables of factors calculated at the initialization stage
        !
        REAL(DP), SAVE :: tablez(ltabl, ndims)
        REAL(DP) :: work(4 * nz * nsl)
        COMPLEX(DP) :: DUMMY
        INTEGER, SAVE :: isys = 1
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
            idir = -1
            tscale = 1.0_DP / nz
        ELSE IF (isign > 0) THEN
            idir = 1
            tscale = 1.0_DP
        END IF
        !
        IF (isign /= 0) &
            CALL ZZFFTM(idir, nz, nsl, tscale, c(1), ldz, cout(1), ldz, tablez(1, ip), &
                        work, isys)
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
            CALL ZZFFTM(0, nz, 1, 1.0_DP, DUMMY, ldz, DUMMY, ldz, tablez(1, icurrent), &
                        work, isys)
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
    !! parameters nx,ny,nzl,ldx) are stored and re-used if available
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
        INTEGER, PARAMETER :: ltabl = 2 * nfftx + 64
        REAL(DP), SAVE :: tablex(ltabl, ndims), tabley(ltabl, ndims)
        REAL(DP) :: work(4 * nx * ny)
        COMPLEX(DP) :: XY(ldx * ny)
        COMPLEX(DP) :: DUMMY
        INTEGER, SAVE :: isys = 1
        !
        CHARACTER(LEN=80) :: sub_name = 'env_fft_scalar: env_cft_2xy'
        !
        !--------------------------------------------------------------------------------
        !
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
        IF (isign < 0) THEN
            idir = -1
            tscale = 1.0_DP / (nx * ny)
            !
            DO k = 0, nzl - 1
                kk = k * ldx * ldy
                !
                !------------------------------------------------------------------------
                ! FORWARD: ny FFTs in the X direction
                !
                CALL ZZFFTM(idir, nx, ny, tscale, r(kk + 1), ldx, r(kk + 1), ldx, &
                            tablex(1, ip), work(1), isys)
                !
                !------------------------------------------------------------------------
                ! FORWARD: nx FFTs in the Y direction
                !
                DO i = 1, nx
                    !
                    IF (dofft(i)) THEN
                        !
                        DO j = 0, ny - 1
                            XY(j + 1) = r(i + (j) * ldx + kk)
                        END DO
                        !
                        CALL ZZFFT(idir, ny, 1.0_DP, XY, XY, tabley(1, ip), work(1), &
                                   isys)
                        !
                        DO j = 0, ny - 1
                            r(i + (j) * ldx + kk) = XY(j + 1)
                        END DO
                        !
                    END IF
                    !
                END DO
                !
            END DO
            !
        ELSE IF (isign > 0) THEN
            !
            idir = 1
            tscale = 1.0_DP
            !
            DO k = 0, nzl - 1
                !
                !------------------------------------------------------------------------
                ! BACKWARD: nx FFTs in the Y direction
                !
                kk = (k) * ldx * ldy
                !
                DO i = 1, nx
                    !
                    IF (dofft(i)) THEN
                        !
                        DO j = 0, ny - 1
                            XY(j + 1) = r(i + (j) * ldx + kk)
                        END DO
                        !
                        CALL ZZFFT(idir, ny, 1.0_DP, XY, XY, tabley(1, ip), work(1), &
                                   isys)
                        !
                        DO j = 0, ny - 1
                            r(i + (j) * ldx + kk) = XY(j + 1)
                        END DO
                        !
                    END IF
                    !
                END DO
                !
                !------------------------------------------------------------------------
                ! BACKWARD: ny FFTs in the X direction
                !
                CALL ZZFFTM(idir, nx, ny, tscale, r(kk + 1), ldx, r(kk + 1), ldx, &
                            tablex(1, ip), work(1), isys)
                !
            END DO
            !
        END IF
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
            CALL ZZFFT(0, ny, 1.0_DP, DUMMY, DUMMY, tabley(1, icurrent), work, isys)
            !
            CALL ZZFFTM(0, nx, 1, 1.0_DP, DUMMY, ldx, DUMMY, ldx, tablex(1, icurrent), &
                        work, isys)
            !
            dims(1, icurrent) = ny
            dims(2, icurrent) = ldx
            dims(3, icurrent) = nx
            dims(4, icurrent) = nzl
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
    !! parameters nx,ny,nz) are stored and re-used if available
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
        INTEGER, PARAMETER :: ltabl = 60
        INTEGER, PARAMETER :: lwork = 195 + 6 * nfftx
        INTEGER, SAVE :: iw0(ltabl, ndims)
        INTEGER :: k_off, kj_offset
        REAL(DP), SAVE :: auxp(lwork, ndims)
        COMPLEX(DP), ALLOCATABLE :: cw2(:)
        COMPLEX(DP) :: f_out(SIZE(f))
        !
#if defined(ASL)&&defined(MICRO)
        INTEGER :: nbtasks
        COMMON / NEC_ASL_PARA / nbtasks
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
                          'homany different from 1 not yet implemented for SX6', 1)
        !
#if defined(ASL)
        ALLOCATE (cw2(ldx * ldy * ldz))
        !
        CALL zfc3cl(f(1), nx, ny, nz, ldx, ldy, ldz, err)
        !
#else
        ALLOCATE (cw2(6 * ldx * ldy * ldz))
#endif
        !
        !--------------------------------------------------------------------------------
        ! Here initialize table only if necessary
        !
        ip = -1
        !
        CALL lookup()
        !
        IF (ip == -1) CALL init_plan()
        ! no table exist for these parameters; initialize a new one
        !
        !--------------------------------------------------------------------------------
        ! Now perform the 3D FFT using the machine specific driver
        !
        IF (ip == -1) CALL init_plan()
        ! no table exist for these parameters; initialize a new one
        !
        !--------------------------------------------------------------------------------
        ! Now perform the 3D FFT using the machine specific driver
        !
#if defined(ASL)
#if defined(MICRO)
        CALL hfc3bf(nx, ny, nz, f(1), ldx, ldy, ldz, -isign, iw0(1, ip), auxp(1, ip), &
                    cw2(1), nbtasks, err)
#else
        CALL zfc3bf(nx, ny, nz, f(1), ldx, ldy, ldz, -isign, iw0(1, ip), auxp(1, ip), &
                    cw2(1), err)
#endif
        IF (isign < 0) THEN
            tscale = 1.0_DP / DBLE(nx * ny * nz)
            !
            CALL ZDSCAL(ldx * ldy * ldz, tscale, f(1), 1)
            !
        END IF
#else
        err = 0
        !
        tscale = 1.0_DP
        !
        IF (isign < 0) tscale = tscale / DBLE(nx * ny * nz)
        !
        CALL ZZFFT3D(isign, nx, ny, nz, tscale, f(1), ldx, ldy, f_out(1), ldx, ldy, &
                     auxp(1, ip), cw2(1), err)
        !
        !$omp parallel do private(j, i, k_off, kj_offset)
        DO k = 1, nz
            k_off = (k - 1) * ldx * ldy
            !
            DO j = 1, ny
                kj_offset = (j - 1) * ldx + k_off
                !
                DO i = 1, nx
                    f(i + kj_offset) = f_out(i + kj_offset)
                END DO
                !
            END DO
            !
        END DO
        !$omp end parallel do
#endif
        !
        IF (err /= 0) CALL io%error(sub_name, 'FFT returned an error', err)
        !
        DEALLOCATE (cw2)
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
#if defined(ASL)
#if defined(MICRO)
            CALL hfc3fb(nx, ny, nz, f(1), ldx, ldy, ldz, 0, iw0(1, icurrent), &
                        auxp(1, icurrent), cw2(1), nbtasks, err)
#else
            CALL zfc3fb(nx, ny, nz, f(1), ldx, ldy, ldz, 0, iw0(1, icurrent), &
                        auxp(1, icurrent), cw2(1), err)
#endif
#else
            !
            err = 0
            !
            CALL ZZFFT3D(0, nx, ny, nz, 1.0_DP, f(1), ldx, ldy, f(1), ldx, ldy, &
                         auxp(1, icurrent), cw2(1), err)
#endif
            !
            IF (err /= 0) CALL io%error(sub_name, 'FFT init returned an error', err)
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
END MODULE env_fft_scalar_sx6
!----------------------------------------------------------------------------------------
#endif
