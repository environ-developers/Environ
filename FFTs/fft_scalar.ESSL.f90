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
#if defined(__LINUX_ESSL)
#if defined(_OPENMP)&&defined(__FFT_SCALAR_THREAD_SAFE)
! thread safety guard
#error ESSL is not compatiable with __FFT_SCALAR_THREAD_SAFE
#endif
!
!----------------------------------------------------------------------------------------
!>
!! IBM ESSL routines
!!
!----------------------------------------------------------------------------------------
MODULE env_fft_scalar_essl
    !------------------------------------------------------------------------------------
    !
    USE env_fft_param
    !
    !------------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    SAVE
    !
    ! Local Parameter
    !
    ! Workspace that is statically allocated is defined here
    ! in order to avoid multiple copies of the same workspace
    ! lwork:   Dimension of the work space array (if any)
    !
    ! ESSL IBM library: see the ESSL manual for DCFT
    !
    INTEGER, PARAMETER :: lwork = 20000 + (2 * nfftx + 256) * 64 + 3 * nfftx
    REAL(DP) :: work(lwork)
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
        INTEGER, PARAMETER :: ltabl = 20000 + 3 * nfftx
        ! dimension of the tables of factors calculated at the initialization stage
        !
        REAL(DP), SAVE :: fw_tablez(ltabl, ndims)
        REAL(DP), SAVE :: bw_tablez(ltabl, ndims)
        !
        CHARACTER(LEN=80) :: sub_name = 'env_fft_scalar: env_cft_1z'
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
#if defined(__FFT_CLOCKS)
        CALL env_start_clock(sub_name)
#endif
        !
        !--------------------------------------------------------------------------------
        ! ESSL uses a different convention for forward/backward transforms
        ! w.r.t most other implementations: notice the sign of "idir"
        !
        IF (isign < 0) THEN
            idir = +1
            tscale = 1.0_DP / nz
            !
            CALL DCFT(0, c(1), 1, ldz, cout(1), 1, ldz, nz, nsl, idir, tscale, &
                      fw_tablez(1, ip), ltabl, work, lwork)
            !
        ELSE IF (isign > 0) THEN
            idir = -1
            tscale = 1.0_DP
            !
            CALL DCFT(0, c(1), 1, ldz, cout(1), 1, ldz, nz, nsl, idir, tscale, &
                      bw_tablez(1, ip), ltabl, work, lwork)
            !
        END IF
        !
#if defined(__FFT_CLOCKS)
        CALL env_stop_clock(sub_name)
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
                !
                found = (nz == zdims(1, ip)) .AND. &
                        (nsl == zdims(2, ip)) .AND. &
                        (ldz == zdims(3, ip))
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
            CALL DCFT(1, c(1), 1, ldz, cout(1), 1, ldz, nz, nsl, 1, tscale, &
                      fw_tablez(1, icurrent), ltabl, work(1), lwork)
            !
            CALL DCFT(1, c(1), 1, ldz, cout(1), 1, ldz, nz, nsl, -1, 1.0_DP, &
                      bw_tablez(1, icurrent), ltabl, work(1), lwork)
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
        INTEGER, PARAMETER :: ltabl = 20000 + 3 * nfftx
        REAL(DP), SAVE :: fw_tablex(ltabl, ndims), fw_tabley(ltabl, ndims)
        REAL(DP), SAVE :: bw_tablex(ltabl, ndims), bw_tabley(ltabl, ndims)
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
                CALL env_errore(sub_name, 'Wrong dimension for arg no. 8', 1)
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
                CALL DCFT(0, r(kk), 1, ldx, r(kk), 1, ldx, nx, ny, &
                          1, tscale, fw_tablex(1, ip), ltabl, work(1), lwork)
                !
                CALL DCFT(0, r(kk), ldx, 1, r(kk), ldx, 1, ny, nx, &
                          1, 1.0_DP, fw_tabley(1, ip), ltabl, work(1), lwork)
                !
            END DO
            !
        ELSE IF (isign > 0) THEN
            !
            DO k = 1, nzl
                kk = 1 + (k - 1) * ldx * ldy
                !
                CALL DCFT(0, r(kk), ldx, 1, r(kk), ldx, 1, ny, nx, &
                          -1, 1.0_DP, bw_tabley(1, ip), ltabl, work(1), lwork)
                !
                CALL DCFT(0, r(kk), 1, ldx, r(kk), 1, ldx, nx, ny, &
                          -1, 1.0_DP, bw_tablex(1, ip), ltabl, work(1), lwork)
                !
            END DO
            !
        END IF
#else
        !
        IF (isign < 0) THEN
            idir = 1
            tscale = 1.0_DP / (nx * ny)
            !
            DO k = 1, nzl
                kk = 1 + (k - 1) * ldx * ldy
                !
                CALL DCFT(0, r(kk), 1, ldx, r(kk), 1, ldx, nx, ny, idir, &
                          tscale, fw_tablex(1, ip), ltabl, work(1), lwork)
                !
                DO i = 1, nx
                    !
                    IF (dofft(i)) THEN
                        kk = i + (k - 1) * ldx * ldy
                        !
                        CALL DCFT(0, r(kk), ldx, 1, r(kk), ldx, 1, ny, 1, &
                                  idir, 1.0_DP, fw_tabley(1, ip), ltabl, work(1), lwork)
                        !
                    END IF
                    !
                END DO
                !
            END DO
            !
        ELSE IF (isign > 0) THEN
            idir = -1
            !
            DO k = 1, nzl
                !
                DO i = 1, nx
                    !
                    IF (dofft(i)) THEN
                        kk = i + (k - 1) * ldx * ldy
                        !
                        CALL DCFT(0, r(kk), ldx, 1, r(kk), ldx, 1, ny, 1, &
                                  idir, 1.0_DP, bw_tabley(1, ip), ltabl, work(1), lwork)
                        !
                    END IF
                    !
                END DO
                !
                kk = 1 + (k - 1) * ldx * ldy
                !
                CALL DCFT(0, r(kk), 1, ldx, r(kk), 1, ldx, nx, ny, idir, &
                          1.0_DP, bw_tablex(1, ip), ltabl, work(1), lwork)
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
#if defined(_OPENMP)
            tscale = 1.0_DP / (nx * ny)
            !
            CALL DCFT(1, r(1), ldx, 1, r(1), ldx, 1, ny, nx, 1, 1.0_DP, &
                      fw_tabley(1, icurrent), ltabl, work(1), lwork)
            !
            CALL DCFT(1, r(1), ldx, 1, r(1), ldx, 1, ny, nx, -1, 1.0_DP, &
                      bw_tabley(1, icurrent), ltabl, work(1), lwork)
            !
            CALL DCFT(1, r(1), 1, ldx, r(1), 1, ldx, nx, ny, 1, &
                      tscale, fw_tablex(1, icurrent), ltabl, work(1), lwork)
            !
            CALL DCFT(1, r(1), 1, ldx, r(1), 1, ldx, nx, ny, -1, &
                      1.0_DP, bw_tablex(1, icurrent), ltabl, work(1), lwork)
            !
#else
            !
            tscale = 1.0_DP / (nx * ny)
            !
            CALL DCFT(1, r(1), ldx, 1, r(1), ldx, 1, ny, 1, 1, 1.0_DP, &
                      fw_tabley(1, icurrent), ltabl, work(1), lwork)
            !
            CALL DCFT(1, r(1), ldx, 1, r(1), ldx, 1, ny, 1, -1, 1.0_DP, &
                      bw_tabley(1, icurrent), ltabl, work(1), lwork)
            !
            CALL DCFT(1, r(1), 1, ldx, r(1), 1, ldx, nx, ny, 1, &
                      tscale, fw_tablex(1, icurrent), ltabl, work(1), lwork)
            !
            CALL DCFT(1, r(1), 1, ldx, r(1), 1, ldx, nx, ny, -1, &
                      1.0_DP, bw_tablex(1, icurrent), ltabl, work(1), lwork)
            !
#endif
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
        CHARACTER(LEN=80) :: sub_name = 'env_fft_scalar: env_cfft3d'
        !
        !--------------------------------------------------------------------------------
        !
        IF (nx < 1) CALL env_errore(sub_name, 'nx is less than 1', 1)
        !
        IF (ny < 1) CALL env_errore(sub_name, 'ny is less than 1', 1)
        !
        IF (nz < 1) CALL env_errore(sub_name, 'nz is less than 1', 1)
        !
        IF (howmany /= 1) &
            CALL env_errore(sub_name, &
                            'howmany different from 1 not yet implemented for ESSL', 1)
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
        IF (isign < 0) THEN
            tscale = 1.0_DP / (nx * ny * nz)
            idir = +1
        ELSE IF (isign > 0) THEN
            tscale = 1.0_DP
            idir = -1
        END IF
        !
        IF (isign /= 0) &
            CALL dcft3(f(1), ldx, ldx * ldy, f(1), ldx, ldx * ldy, &
                       nx, ny, nz, idir, tscale, work(1), lwork)
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
        !! No initialization for 3d FFT's from ESSL
        !!
        !--------------------------------------------------------------------------------
        SUBROUTINE init_plan()
            !----------------------------------------------------------------------------
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
END MODULE env_fft_scalar_essl
!----------------------------------------------------------------------------------------
#endif
