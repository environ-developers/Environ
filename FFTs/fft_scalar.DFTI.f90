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
#if defined(__DFTI)
#if defined(_OPENMP)&&defined(__FFT_SCALAR_THREAD_SAFE)
! compiler safeguard for thread-safe eligibility
#if defined(__PGI)
#error PGI compiler breaks the use of __FFT_SCALAR_THREAD_SAFE in DFTI
#endif
#endif
!
#include "mkl_dfti.f90"
!
!----------------------------------------------------------------------------------------
!>
!! Intel DFTI routines
!!
!----------------------------------------------------------------------------------------
MODULE env_fft_scalar_dfti
    !------------------------------------------------------------------------------------
    !
    USE class_io, ONLY: io
    !
    USE env_fft_param
    !
    USE MKL_DFTI ! this can be found in the MKL include directory
    !
    !------------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    TYPE env_dfti_descriptor_array
        TYPE(DFTI_DESCRIPTOR), POINTER :: desc
    END TYPE
    !------------------------------------------------------------------------------------
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
    SUBROUTINE env_cft_1z(c, nsl, nz, ldz, isign, cout, in_place)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: isign, nsl, nz, ldz
        LOGICAL, INTENT(IN), OPTIONAL :: in_place
        !
        COMPLEX(DP) :: c(:), cout(:)
        !
        REAL(DP) :: tscale
        INTEGER :: i, ip
        INTEGER, SAVE :: zdims(3, ndims) = -1
        INTEGER, SAVE :: icurrent = 1
        LOGICAL :: found
        !
        TYPE(env_dfti_descriptor_array), SAVE :: hand(ndims)
        ! Intel MKL native FFT driver
        !
        LOGICAL, SAVE :: dfti_first = .TRUE.
        LOGICAL, SAVE :: is_inplace
        INTEGER :: dfti_status = 0
        INTEGER :: placement
        !
        CHARACTER(LEN=80) :: sub_name = 'env_fft_scalar: env_cft_1z'
        !
        !--------------------------------------------------------------------------------
        !
#if defined(__FFT_SCALAR_THREAD_SAFE)
        !$omp threadprivate(hand, dfti_first, zdims, icurrent, is_inplace)
#endif
        !
        IF (PRESENT(in_place)) THEN
            is_inplace = in_place
        ELSE
            is_inplace = .FALSE.
        END IF
        !
        !--------------------------------------------------------------------------------
        ! Check dimensions and corner cases
        !
        ! Starting from MKL 2019, it is no longer possible to define "empty" plans,
        ! i.e. plans with 0 FFTs. Just return immediately in this case
        !
        IF (nsl <= 0) THEN
            !
            IF (nsl < 0) CALL io%error(sub_name, 'nsl out of range', nsl)
            !
            RETURN
            !
        END IF
        !
        !--------------------------------------------------------------------------------
        ! Here initialize table only if necessary
        !
        CALL lookup()
        !
        IF (.NOT. found) CALL init_dfti()
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
            IF (is_inplace) THEN
                dfti_status = DftiComputeForward(hand(ip)%desc, c)
            ELSE
                dfti_status = DftiComputeForward(hand(ip)%desc, c, cout)
            END IF
            !
            IF (dfti_status /= 0) &
                !
                CALL io%error(sub_name, &
                              'Stopped in DftiComputeForward '// &
                              DftiErrorMessage(dfti_status), dfti_status)
            !
        ELSE IF (isign > 0) THEN
            !
            IF (is_inplace) THEN
                dfti_status = DftiComputeBackward(hand(ip)%desc, c)
            ELSE
                dfti_status = DftiComputeBackward(hand(ip)%desc, c, cout)
            END IF
            !
            IF (dfti_status /= 0) &
                CALL io%error(sub_name, &
                              'Stopped in DftiComputeBackward '// &
                              DftiErrorMessage(dfti_status), dfti_status)
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
            IF (dfti_first) THEN
                !
                DO ip = 1, ndims
                    hand(ip)%desc => NULL()
                END DO
                !
                dfti_first = .FALSE.
            END IF
            !
            DO ip = 1, ndims
                !
                found = (nz == zdims(1, ip)) .AND. &
                        (nsl == zdims(2, ip)) .AND. &
                        (ldz == zdims(3, ip))
                !
                dfti_status = DftiGetValue(hand(ip)%desc, DFTI_PLACEMENT, placement)
                !
                found = found .AND. is_inplace .AND. (placement == DFTI_INPLACE)
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
        SUBROUTINE init_dfti()
            !----------------------------------------------------------------------------
            !
            IF (ASSOCIATED(hand(icurrent)%desc)) THEN
                dfti_status = DftiFreeDescriptor(hand(icurrent)%desc)
                !
                IF (dfti_status /= 0) THEN
                    WRITE (*, *) 'Stopped in DftiFreeDescriptor', dfti_status
                    !
                    STOP
                    !
                END IF
                !
            END IF
            !
            dfti_status = DftiCreateDescriptor(hand(icurrent)%desc, DFTI_DOUBLE, &
                                               DFTI_COMPLEX, 1, nz)
            !
            IF (dfti_status /= 0) &
                CALL io%error(sub_name, 'Stopped in DftiCreateDescriptor', dfti_status)
            !
            dfti_status = DftiSetValue(hand(icurrent)%desc, &
                                       DFTI_NUMBER_OF_TRANSFORMS, nsl)
            !
            IF (dfti_status /= 0) &
                CALL io%error(sub_name, &
                              'Stopped in DFTI_NUMBER_OF_TRANSFORMS', dfti_status)
            !
            dfti_status = DftiSetValue(hand(icurrent)%desc, DFTI_INPUT_DISTANCE, ldz)
            !
            IF (dfti_status /= 0) &
                CALL io%error(sub_name, 'Stopped in DFTI_INPUT_DISTANCE', dfti_status)
            !
            IF (is_inplace) THEN
                !
                dfti_status = DftiSetValue(hand(icurrent)%desc, DFTI_PLACEMENT, &
                                           DFTI_INPLACE)
                !
            ELSE
                !
                dfti_status = DftiSetValue(hand(icurrent)%desc, DFTI_PLACEMENT, &
                                           DFTI_NOT_INPLACE)
                !
            END IF
            !
            IF (dfti_status /= 0) &
                CALL io%error(sub_name, 'Stopped in DFTI_PLACEMENT', dfti_status)
            !
            dfti_status = DftiSetValue(hand(icurrent)%desc, DFTI_OUTPUT_DISTANCE, ldz)
            !
            IF (dfti_status /= 0) &
                CALL io%error(sub_name, 'Stopped in DFTI_OUTPUT_DISTANCE', dfti_status)
            !
            tscale = 1.0_DP / nz
            dfti_status = DftiSetValue(hand(icurrent)%desc, DFTI_FORWARD_SCALE, tscale)
            !
            IF (dfti_status /= 0) &
                CALL io%error(sub_name, 'Stopped in DFTI_FORWARD_SCALE', dfti_status)
            !
            dfti_status = DftiSetValue(hand(icurrent)%desc, DFTI_BACKWARD_SCALE, DBLE(1))
            !
            IF (dfti_status /= 0) &
                CALL io%error(sub_name, 'Stopped in DFTI_BACKWARD_SCALE', dfti_status)
            !
            dfti_status = DftiCommitDescriptor(hand(icurrent)%desc)
            !
            IF (dfti_status /= 0) &
                CALL io%error(sub_name, 'Stopped in DftiCommitDescriptor', dfti_status)
            !
            zdims(1, icurrent) = nz
            zdims(2, icurrent) = nsl
            zdims(3, icurrent) = ldz
            ip = icurrent
            icurrent = MOD(icurrent, ndims) + 1
            !
            !----------------------------------------------------------------------------
        END SUBROUTINE init_dfti
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
        INTEGER :: i, k, j, ip
        REAL(DP) :: tscale
        INTEGER, SAVE :: icurrent = 1
        INTEGER, SAVE :: dims(4, ndims) = -1
        LOGICAL :: dofft(nfftx), found
        !
        TYPE(env_dfti_descriptor_array), SAVE :: hand(ndims)
        LOGICAL, SAVE :: dfti_first = .TRUE.
        INTEGER :: dfti_status = 0
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
        IF (.NOT. found) CALL init_dfti()
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
            dfti_status = DftiComputeForward(hand(ip)%desc, r(:))
            !
            IF (dfti_status /= 0) THEN
                WRITE (*, *) 'Stopped in DftiComputeForward', dfti_status
                !
                STOP
                !
            END IF
            !
        ELSE IF (isign > 0) THEN
            dfti_status = DftiComputeBackward(hand(ip)%desc, r(:))
            !
            IF (dfti_status /= 0) THEN
                WRITE (*, *) 'Stopped in DftiComputeBackward', dfti_status
                !
                STOP
                !
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
            IF (dfti_first) THEN
                !
                DO ip = 1, ndims
                    hand(ip)%desc => NULL()
                END DO
                !
                dfti_first = .FALSE.
            END IF
            !
            DO ip = 1, ndims
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
        SUBROUTINE init_dfti()
            !----------------------------------------------------------------------------
            !
            IF (ASSOCIATED(hand(icurrent)%desc)) THEN
                dfti_status = DftiFreeDescriptor(hand(icurrent)%desc)
                !
                IF (dfti_status /= 0) THEN
                    WRITE (*, *) 'Stopped in DftiFreeDescriptor', dfti_status
                    !
                    STOP
                    !
                END IF
                !
            END IF
            !
            dfti_status = DftiCreateDescriptor(hand(icurrent)%desc, DFTI_DOUBLE, &
                                               DFTI_COMPLEX, 2, (/nx, ny/))
            !
            IF (dfti_status /= 0) THEN
                WRITE (*, *) 'Stopped in DftiCreateDescriptor', dfti_status
                !
                STOP
                !
            END IF
            !
            dfti_status = DftiSetValue(hand(icurrent)%desc, &
                                       DFTI_NUMBER_OF_TRANSFORMS, nzl)
            !
            IF (dfti_status /= 0) THEN
                WRITE (*, *) 'Stopped in DFTI_NUMBER_OF_TRANSFORMS', dfti_status
                !
                STOP
                !
            END IF
            !
            dfti_status = DftiSetValue(hand(icurrent)%desc, &
                                       DFTI_INPUT_DISTANCE, ldx * ldy)
            !
            IF (dfti_status /= 0) THEN
                WRITE (*, *) 'Stopped in DFTI_INPUT_DISTANCE', dfti_status
                !
                STOP
                !
            END IF
            !
            dfti_status = DftiSetValue(hand(icurrent)%desc, DFTI_PLACEMENT, &
                                       DFTI_INPLACE)
            !
            IF (dfti_status /= 0) THEN
                WRITE (*, *) 'Stopped in DFTI_PLACEMENT', dfti_status
                !
                STOP
                !
            END IF
            !
            tscale = 1.0_DP / (nx * ny)
            dfti_status = DftiSetValue(hand(icurrent)%desc, &
                                       DFTI_FORWARD_SCALE, tscale)
            !
            IF (dfti_status /= 0) THEN
                WRITE (*, *) 'Stopped in DFTI_FORWARD_SCALE', dfti_status
                !
                STOP
                !
            END IF
            !
            dfti_status = DftiSetValue(hand(icurrent)%desc, &
                                       DFTI_BACKWARD_SCALE, DBLE(1))
            !
            IF (dfti_status /= 0) THEN
                WRITE (*, *) 'Stopped in DFTI_BACKWARD_SCALE', dfti_status
                !
                STOP
                !
            END IF
            !
            dfti_status = DftiCommitDescriptor(hand(icurrent)%desc)
            !
            IF (dfti_status /= 0) THEN
                WRITE (*, *) 'Stopped in DftiCommitDescriptor', dfti_status
                !
                STOP
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
        END SUBROUTINE init_dfti
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
        INTEGER, SAVE :: dims(4, ndims) = -1
        !
        TYPE(env_dfti_descriptor_array), SAVE :: hand(ndims)
        ! Intel MKL native FFT driver
        !
        LOGICAL, SAVE :: dfti_first = .TRUE.
        INTEGER :: dfti_status = 0
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
        IF (howmany < 1) CALL io%error(sub_name, 'howmany is less than 1', 1)
        !
        !--------------------------------------------------------------------------------
        ! Here initialize table only if necessary
        !
        CALL lookup()
        !
        IF (ip == -1) CALL init_dfti()
        ! no table exist for these parameters; initialize a new one
        !
        !--------------------------------------------------------------------------------
        ! Now perform the 3D FFT using the machine specific driver
        !
        IF (isign < 0) THEN
            dfti_status = DftiComputeForward(hand(ip)%desc, f(1:))
            !
            IF (dfti_status /= 0) THEN
                WRITE (*, *) 'Stopped in env_cfft3d, DftiComputeForward', dfti_status
                !
                STOP
                !
            END IF
            !
        ELSE IF (isign > 0) THEN
            dfti_status = DftiComputeBackward(hand(ip)%desc, f(1:))
            !
            IF (dfti_status /= 0) THEN
                WRITE (*, *) 'Stopped in env_cfft3d, DftiComputeBackward', dfti_status
                !
                STOP
                !
            END IF
            !
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
            IF (dfti_first) THEN
                !
                DO ip = 1, ndims
                    hand(ip)%desc => NULL()
                END DO
                !
                dfti_first = .FALSE.
            END IF
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
        SUBROUTINE init_dfti()
            !----------------------------------------------------------------------------
            !
            IF (ASSOCIATED(hand(icurrent)%desc)) THEN
                dfti_status = DftiFreeDescriptor(hand(icurrent)%desc)
                !
                IF (dfti_status /= 0) THEN
                    WRITE (*, *) 'Stopped in env_cfft3d, DftiFreeDescriptor', dfti_status
                    !
                    STOP
                    !
                END IF
                !
            END IF
            !
            dfti_status = DftiCreateDescriptor(hand(icurrent)%desc, DFTI_DOUBLE, &
                                               DFTI_COMPLEX, 3, (/nx, ny, nz/))
            !
            IF (dfti_status /= 0) THEN
                WRITE (*, *) 'Stopped in env_cfft3d, DftiCreateDescriptor', dfti_status
                !
                STOP
                !
            END IF
            !
            dfti_status = DftiSetValue(hand(icurrent)%desc, &
                                       DFTI_NUMBER_OF_TRANSFORMS, howmany)
            !
            IF (dfti_status /= 0) THEN
                !
                WRITE (*, *) 'Stopped in env_cfft3d, DFTI_NUMBER_OF_TRANSFORMS', &
                    dfti_status
                !
                STOP
                !
            END IF
            !
            dfti_status = DftiSetValue(hand(icurrent)%desc, &
                                       DFTI_INPUT_DISTANCE, ldx * ldy * ldz)
            !
            IF (dfti_status /= 0) THEN
                WRITE (*, *) 'Stopped in cfft3dm, DFTI_INPUT_DISTANCE', dfti_status
                !
                STOP
                !
            END IF
            !
            dfti_status = DftiSetValue(hand(icurrent)%desc, &
                                       DFTI_PLACEMENT, DFTI_INPLACE)
            !
            IF (dfti_status /= 0) THEN
                WRITE (*, *) 'Stopped in env_cfft3d, DFTI_PLACEMENT', dfti_status
                !
                STOP
                !
            END IF
            !
            tscale = 1.0_DP / (nx * ny * nz)
            dfti_status = DftiSetValue(hand(icurrent)%desc, &
                                       DFTI_FORWARD_SCALE, tscale)
            !
            IF (dfti_status /= 0) THEN
                WRITE (*, *) 'Stopped in env_cfft3d, DFTI_FORWARD_SCALE', dfti_status
                !
                STOP
                !
            END IF
            !
            tscale = 1.0_DP
            dfti_status = DftiSetValue(hand(icurrent)%desc, &
                                       DFTI_BACKWARD_SCALE, tscale)
            !
            IF (dfti_status /= 0) THEN
                WRITE (*, *) 'Stopped in env_cfft3d, DFTI_BACKWARD_SCALE', dfti_status
                !
                STOP
                !
            END IF
            !
            dfti_status = DftiCommitDescriptor(hand(icurrent)%desc)
            !
            IF (dfti_status /= 0) THEN
                WRITE (*, *) 'Stopped in env_cfft3d, DftiCreateDescriptor', dfti_status
                !
                STOP
                !
            END IF
            !
            dims(1, icurrent) = nx
            dims(2, icurrent) = ny
            dims(3, icurrent) = nz
            dims(4, icurrent) = howmany
            ip = icurrent
            icurrent = MOD(icurrent, ndims) + 1
            !
            !----------------------------------------------------------------------------
        END SUBROUTINE init_dfti
        !--------------------------------------------------------------------------------
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_cfft3d
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE env_fft_scalar_dfti
!----------------------------------------------------------------------------------------
#endif