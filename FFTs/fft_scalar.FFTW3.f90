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
!          stick and plane revision - Stefano de Gironcoli - September 2016
!          Modified by Edan Bainglass
!
!----------------------------------------------------------------------------------------
!
#if defined(__FFTW3)
#if defined(_OPENMP)&&defined(__FFT_SCALAR_THREAD_SAFE)
! thread safety guard
#error FFTW3 is not compatible with __FFT_SCALAR_THREAD_SAFE
#endif
!
!----------------------------------------------------------------------------------------
!>
!! FFTW v.3 routines
!!
!----------------------------------------------------------------------------------------
MODULE env_fft_scalar_fftw3
    !------------------------------------------------------------------------------------
    !
    USE env_base_io, ONLY: io
    !
    USE env_fft_param
    !
    USE, INTRINSIC :: ISO_C_BINDING
    ! iso_c_binding provides C_PTR, C_NULL_PTR, C_ASSOCIATED
    !
#if defined(_OPENMP)
#include "fftw3.f03"
#else
#include "fftw3.f"
#endif
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
        INTEGER :: i, err, idir, ip
        INTEGER, SAVE :: zdims(3, ndims) = -1
        INTEGER, SAVE :: icurrent = 1
        LOGICAL :: found
        !
        TYPE(C_PTR), SAVE :: fw_planz(ndims) = C_NULL_PTR
        TYPE(C_PTR), SAVE :: bw_planz(ndims) = C_NULL_PTR
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
            !
            CALL dfftw_execute_dft(fw_planz(ip), c, cout)
            !
            tscale = 1.0_DP / nz
            cout(1:ldz * nsl) = cout(1:ldz * nsl) * tscale
        ELSE IF (isign > 0) THEN
            CALL dfftw_execute_dft(bw_planz(ip), c, cout)
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
            IMPLICIT NONE
            !
#if defined(_OPENMP)
            INTEGER :: void
            INTEGER, EXTERNAL :: omp_get_max_threads
            !
            !----------------------------------------------------------------------------
            !
            CALL dfftw_cleanup_threads()
            !
            void = fftw_init_threads()
            !
            CALL dfftw_plan_with_nthreads(omp_get_max_threads())
#endif
            !
            IF (C_ASSOCIATED(fw_planz(icurrent))) &
                CALL dfftw_destroy_plan(fw_planz(icurrent))
            !
            IF (C_ASSOCIATED(bw_planz(icurrent))) &
                CALL dfftw_destroy_plan(bw_planz(icurrent))
            !
            idir = -1
            !
            CALL dfftw_plan_many_dft(fw_planz(icurrent), 1, nz, nsl, c, &
                                     (/SIZE(c)/), 1, ldz, cout, &
                                     (/SIZE(cout)/), 1, ldz, idir, FFTW_ESTIMATE)
            !
            idir = 1
            !
            CALL dfftw_plan_many_dft(bw_planz(icurrent), 1, nz, nsl, c, &
                                     (/SIZE(c)/), 1, ldz, cout, &
                                     (/SIZE(cout)/), 1, ldz, idir, FFTW_ESTIMATE)
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
        INTEGER :: i, k, j, err, idir, ip, kk
        REAL(DP) :: tscale
        INTEGER, SAVE :: icurrent = 1
        INTEGER, SAVE :: dims(4, ndims) = -1
        LOGICAL :: dofft(nfftx), found
        !
        TYPE(C_PTR), SAVE :: fw_plan(2, ndims) = C_NULL_PTR
        TYPE(C_PTR), SAVE :: bw_plan(2, ndims) = C_NULL_PTR
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
        IF (ldx /= nx .OR. ldy /= ny) THEN
            !
            IF (isign < 0) THEN
                !
                DO j = 0, nzl - 1
                    !
                    CALL dfftw_execute_dft(fw_plan(1, ip), &
                                           r(1 + j * ldx * ldy:), r(1 + j * ldx * ldy:))
                    !
                END DO
                !
                DO i = 1, nx
                    !
                    DO k = 1, nzl
                        !
                        IF (dofft(i)) THEN
                            j = i + ldx * ldy * (k - 1)
                            !
                            CALL dfftw_execute_dft(fw_plan(2, ip), r(j:), r(j:))
                            !
                        END IF
                        !
                    END DO
                    !
                END DO
                !
                tscale = 1.0_DP / (nx * ny)
                !
                CALL ZDSCAL(ldx * ldy * nzl, tscale, r(1), 1)
                !
            ELSE IF (isign > 0) THEN
                !
                DO i = 1, nx
                    !
                    DO k = 1, nzl
                        !
                        IF (dofft(i)) THEN
                            j = i + ldx * ldy * (k - 1)
                            !
                            CALL dfftw_execute_dft(bw_plan(2, ip), r(j:), r(j:))
                            !
                        END IF
                        !
                    END DO
                    !
                END DO
                !
                DO j = 0, nzl - 1
                    !
                    CALL dfftw_execute_dft(bw_plan(1, ip), &
                                           r(1 + j * ldx * ldy:), r(1 + j * ldx * ldy:))
                    !
                END DO
                !
            END IF
            !
        ELSE
            !
            IF (isign < 0) THEN
                !
                CALL dfftw_execute_dft(fw_plan(1, ip), r(1:), r(1:))
                !
                tscale = 1.0_DP / (nx * ny)
                !
                CALL ZDSCAL(ldx * ldy * nzl, tscale, r(1), 1)
                !
            ELSE IF (isign > 0) THEN
                CALL dfftw_execute_dft(bw_plan(1, ip), r(1:), r(1:))
            END IF
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
            IMPLICIT NONE
            !
#if defined(_OPENMP)
            INTEGER :: void
            INTEGER, EXTERNAL :: omp_get_max_threads
            !
            !----------------------------------------------------------------------------
            !
            CALL dfftw_cleanup_threads()
            !
            void = fftw_init_threads()
            !
            CALL dfftw_plan_with_nthreads(omp_get_max_threads())
#endif
            !
            IF (ldx /= nx .OR. ldy /= ny) THEN
                !
                IF (C_ASSOCIATED(fw_plan(2, icurrent))) &
                    CALL dfftw_destroy_plan(fw_plan(2, icurrent))
                !
                IF (C_ASSOCIATED(bw_plan(2, icurrent))) &
                    CALL dfftw_destroy_plan(bw_plan(2, icurrent))
                !
                idir = -1
                !
                CALL dfftw_plan_many_dft(fw_plan(2, icurrent), 1, ny, 1, r(1:), &
                                         (/ldx * ldy/), ldx, 1, r(1:), (/ldx * ldy/), &
                                         ldx, 1, idir, FFTW_ESTIMATE)
                !
                idir = 1
                !
                CALL dfftw_plan_many_dft(bw_plan(2, icurrent), 1, ny, 1, r(1:), &
                                         (/ldx * ldy/), ldx, 1, r(1:), (/ldx * ldy/), &
                                         ldx, 1, idir, FFTW_ESTIMATE)
                !
                IF (C_ASSOCIATED(fw_plan(1, icurrent))) &
                    CALL dfftw_destroy_plan(fw_plan(1, icurrent))
                !
                IF (C_ASSOCIATED(bw_plan(1, icurrent))) &
                    CALL dfftw_destroy_plan(bw_plan(1, icurrent))
                !
                idir = -1
                !
                CALL dfftw_plan_many_dft(fw_plan(1, icurrent), 1, nx, ny, r(1:), &
                                         (/ldx * ldy/), 1, ldx, r(1:), (/ldx * ldy/), &
                                         1, ldx, idir, FFTW_ESTIMATE)
                !
                idir = 1
                !
                CALL dfftw_plan_many_dft(bw_plan(1, icurrent), 1, nx, ny, r(1:), &
                                         (/ldx * ldy/), 1, ldx, r(1:), (/ldx * ldy/), &
                                         1, ldx, idir, FFTW_ESTIMATE)
                !
            ELSE
                !
                IF (C_ASSOCIATED(fw_plan(1, icurrent))) &
                    CALL dfftw_destroy_plan(fw_plan(1, icurrent))
                !
                IF (C_ASSOCIATED(bw_plan(1, icurrent))) &
                    CALL dfftw_destroy_plan(bw_plan(1, icurrent))
                !
                idir = -1
                !
                CALL dfftw_plan_many_dft(fw_plan(1, icurrent), 2, (/nx, ny/), nzl, &
                                         r(1:), (/nx, ny/), 1, nx * ny, r(1:), &
                                         (/nx, ny/), 1, nx * ny, idir, FFTW_ESTIMATE)
                !
                idir = 1
                !
                CALL dfftw_plan_many_dft(bw_plan(1, icurrent), 2, (/nx, ny/), nzl, &
                                         r(1:), (/nx, ny/), 1, nx * ny, r(1:), &
                                         (/nx, ny/), 1, nx * ny, idir, FFTW_ESTIMATE)
                !
            END IF
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
        TYPE(C_PTR), SAVE :: fw_plan(ndims) = C_NULL_PTR
        TYPE(C_PTR), SAVE :: bw_plan(ndims) = C_NULL_PTR
        !
        CHARACTER(LEN=80) :: sub_name = 'env_fft_scalar: env_cfft3d'
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
                          'howmany different from 1 not yet implemented for FFTW3', 1)
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
            !
            CALL dfftw_execute_dft(fw_plan(ip), f(1:), f(1:))
            !
            tscale = 1.0_DP / DBLE(nx * ny * nz)
            !
            CALL ZDSCAL(nx * ny * nz, tscale, f(1), 1)
            !
        ELSE IF (isign > 0) THEN
            CALL dfftw_execute_dft(bw_plan(ip), f(1:), f(1:))
        END IF
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
            IMPLICIT NONE
            !
#if defined(_OPENMP)
            INTEGER :: void
            INTEGER, EXTERNAL :: omp_get_max_threads
            !
            !----------------------------------------------------------------------------
            !
            CALL dfftw_cleanup_threads()
            !
            void = fftw_init_threads()
            !
            CALL dfftw_plan_with_nthreads(omp_get_max_threads())
#endif
            !
            IF (nx /= ldx .OR. ny /= ldy .OR. nz /= ldz) &
                CALL io%error(sub_name, 'Not implemented', 3)
            !
            IF (C_ASSOCIATED(fw_plan(icurrent))) &
                CALL dfftw_destroy_plan(fw_plan(icurrent))
            !
            IF (C_ASSOCIATED(bw_plan(icurrent))) &
                CALL dfftw_destroy_plan(bw_plan(icurrent))
            !
            idir = -1
            !
            CALL dfftw_plan_dft_3d(fw_plan(icurrent), nx, ny, nz, f(1:), &
                                   f(1:), idir, FFTW_ESTIMATE)
            !
            idir = 1
            !
            CALL dfftw_plan_dft_3d(bw_plan(icurrent), nx, ny, nz, f(1:), &
                                   f(1:), idir, FFTW_ESTIMATE)
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
END MODULE env_fft_scalar_fftw3
!----------------------------------------------------------------------------------------
#endif
