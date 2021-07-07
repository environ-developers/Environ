!----------------------------------------------------------------------------------------
!
! Copyright (C) 2021 ENVIRON (www.quantum-environ.org)
! Copyright (C) 2002-2018 Quantum ESPRESSO (www.quantum-espresso.org)
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
! Authors: Modified by Edan Bainglass
!
!----------------------------------------------------------------------------------------
!>
!! Utility functions to perform memcpy and memset on the device with CUDA Fortran
!! cuf_memXXX contains a CUF KERNEL to perform the selected operation
!!
!----------------------------------------------------------------------------------------
MODULE env_utils_cuda
    !------------------------------------------------------------------------------------
    !
    USE env_utils_param, ONLY: DP
    !
    !------------------------------------------------------------------------------------
    !
    INTERFACE env_cuf_memcpy
        MODULE PROCEDURE &
            env_cuf_memcpy_r_0dd, &
            env_cuf_memcpy_r2d, &
            env_cuf_memcpy_r3d, &
            env_cuf_memcpy_c_0dd, &
            env_cuf_memcpy_c2d, &
            env_cuf_memcpy_c3d
    END INTERFACE
    !
    INTERFACE env_cuf_memset
        MODULE PROCEDURE &
            env_cuf_memset_r_0dd, &
            env_cuf_memset_r2d, &
            env_cuf_memset_r3d, &
            env_cuf_memset_c_0dd, &
            env_cuf_memset_c2d, &
            env_cuf_memset_c3d
    END INTERFACE
    !
    !------------------------------------------------------------------------------------
    !
    PRIVATE
    !
    PUBLIC :: env_cuf_memcpy, env_cuf_memset
    !
    !------------------------------------------------------------------------------------
CONTAINS
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_cuf_memcpy_r_0dd(array_out, array_in, range1)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP), INTENT(IN) :: array_in(:)
        INTEGER, INTENT(IN) :: range1(2)
        !
        REAL(DP), INTENT(INOUT) :: array_out(:)
        !
#if defined(__CUDA)
        attributes(DEVICE) :: array_out, array_in
#endif
        !
        INTEGER :: i1, d1s, d1e
        !
        !--------------------------------------------------------------------------------
        !
        d1s = range1(1)
        d1e = range1(2)
        !
        !$cuf kernel do(1)
        DO i1 = d1s, d1e
            array_out(i1) = array_in(i1)
        END DO
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_cuf_memcpy_r_0dd
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_cuf_memcpy_r2d(array_out, array_in, range1, range2)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP), INTENT(IN) :: array_in(:, :)
        INTEGER, DIMENSION(2), INTENT(IN) :: range1, range2
        !
        REAL(DP), INTENT(INOUT) :: array_out(:, :)
        !
#if defined(__CUDA)
        attributes(DEVICE) :: array_out, array_in
#endif
        !
        INTEGER :: i1, d1s, d1e
        INTEGER :: i2, d2s, d2e
        !
        !--------------------------------------------------------------------------------
        !
        d1s = range1(1)
        d1e = range1(2)
        d2s = range2(1)
        d2e = range2(2)
        !
        !$cuf kernel do(2)
        DO i2 = d2s, d2e
            !
            DO i1 = d1s, d1e
                array_out(i1, i2) = array_in(i1, i2)
            END DO
            !
        END DO
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_cuf_memcpy_r2d
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_cuf_memcpy_r3d(array_out, array_in, range1, range2, range3)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP), INTENT(IN) :: array_in(:, :, :)
        INTEGER, DIMENSION(2), INTENT(IN) :: range1, range2, range3
        !
        REAL(DP), INTENT(INOUT) :: array_out(:, :, :)
        !
#if defined(__CUDA)
        attributes(DEVICE) :: array_out, array_in
#endif
        !
        INTEGER :: i1, d1s, d1e
        INTEGER :: i2, d2s, d2e
        INTEGER :: i3, d3s, d3e
        !
        !--------------------------------------------------------------------------------
        !
        d1s = range1(1)
        d1e = range1(2)
        d2s = range2(1)
        d2e = range2(2)
        d3s = range3(1)
        d3e = range3(2)
        !
        !$cuf kernel do(3)
        DO i3 = d3s, d3e
            !
            DO i2 = d2s, d2e
                !
                DO i1 = d1s, d1e
                    array_out(i1, i2, i3) = array_in(i1, i2, i3)
                END DO
                !
            END DO
            !
        END DO
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_cuf_memcpy_r3d
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_cuf_memcpy_c_0dd(array_out, array_in, range1)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        COMPLEX(DP), INTENT(IN) :: array_in(:)
        INTEGER, INTENT(IN) :: range1(2)
        !
        COMPLEX(DP), INTENT(INOUT) :: array_out(:)
        !
#if defined(__CUDA)
        attributes(DEVICE) :: array_out, array_in
#endif
        !
        INTEGER :: i1, d1s, d1e
        !
        !--------------------------------------------------------------------------------
        !
        d1s = range1(1)
        d1e = range1(2)
        !
        !$cuf kernel do(1)
        DO i1 = d1s, d1e
            array_out(i1) = array_in(i1)
        END DO
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_cuf_memcpy_c_0dd
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_cuf_memcpy_c2d(array_out, array_in, range1, range2)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        COMPLEX(DP), INTENT(IN) :: array_in(:, :)
        INTEGER, DIMENSION(2), INTENT(IN) :: range1, range2
        !
        COMPLEX(DP), INTENT(INOUT) :: array_out(:, :)
        !
#if defined(__CUDA)
        attributes(DEVICE) :: array_out, array_in
#endif
        !
        INTEGER :: i1, d1s, d1e
        INTEGER :: i2, d2s, d2e
        !
        !--------------------------------------------------------------------------------
        !
        d1s = range1(1)
        d1e = range1(2)
        d2s = range2(1)
        d2e = range2(2)
        !
        !$cuf kernel do(2)
        DO i2 = d2s, d2e
            !
            DO i1 = d1s, d1e
                array_out(i1, i2) = array_in(i1, i2)
            END DO
            !
        END DO
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_cuf_memcpy_c2d
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_cuf_memcpy_c3d(array_out, array_in, range1, range2, range3)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        COMPLEX(DP), INTENT(IN) :: array_in(:, :, :)
        INTEGER, DIMENSION(2), INTENT(IN) :: range1, range2, range3
        !
        COMPLEX(DP), INTENT(INOUT) :: array_out(:, :, :)
        !
#if defined(__CUDA)
        attributes(DEVICE) :: array_out, array_in
#endif
        !
        INTEGER :: i1, d1s, d1e
        INTEGER :: i2, d2s, d2e
        INTEGER :: i3, d3s, d3e
        !
        !--------------------------------------------------------------------------------
        !
        d1s = range1(1)
        d1e = range1(2)
        d2s = range2(1)
        d2e = range2(2)
        d3s = range3(1)
        d3e = range3(2)
        !
        !$cuf kernel do(3)
        DO i3 = d3s, d3e
            !
            DO i2 = d2s, d2e
                !
                DO i1 = d1s, d1e
                    array_out(i1, i2, i3) = array_in(i1, i2, i3)
                END DO
                !
            END DO
            !
        END DO
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_cuf_memcpy_c3d
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_cuf_memset_r_0dd(array_out, val, range1)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP), INTENT(IN) :: val
        INTEGER, INTENT(IN) :: range1(2)
        !
        REAL(DP), INTENT(INOUT) :: array_out(:)
        !
#if defined(__CUDA)
        attributes(DEVICE) :: array_out
#endif
        !
        INTEGER :: i1, d1s, d1e
        !
        !--------------------------------------------------------------------------------
        !
        d1s = range1(1)
        d1e = range1(2)
        !
        !$cuf kernel do(1)
        DO i1 = d1s, d1e
            array_out(i1) = val
        END DO
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_cuf_memset_r_0dd
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_cuf_memset_r2d(array_out, val, range1, range2)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP), INTENT(IN) :: val
        INTEGER, DIMENSION(2), INTENT(IN) :: range1, range2
        !
        REAL(DP), INTENT(INOUT) :: array_out(:, :)
        !
#if defined(__CUDA)
        attributes(DEVICE) :: array_out
#endif
        !
        INTEGER :: i1, d1s, d1e
        INTEGER :: i2, d2s, d2e
        !
        !--------------------------------------------------------------------------------
        !
        d1s = range1(1)
        d1e = range1(2)
        d2s = range2(1)
        d2e = range2(2)
        !
        !$cuf kernel do(2)
        DO i2 = d2s, d2e
            !
            DO i1 = d1s, d1e
                array_out(i1, i2) = val
            END DO
            !
        END DO
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_cuf_memset_r2d
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_cuf_memset_r3d(array_out, val, range1, range2, range3)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP), INTENT(IN) :: val
        INTEGER, DIMENSION(2), INTENT(IN) :: range1, range2, range3
        !
        REAL(DP), INTENT(INOUT) :: array_out(:, :, :)
        !
#if defined(__CUDA)
        attributes(DEVICE) :: array_out
#endif
        !
        INTEGER :: i1, d1s, d1e
        INTEGER :: i2, d2s, d2e
        INTEGER :: i3, d3s, d3e
        !
        !--------------------------------------------------------------------------------
        !
        d1s = range1(1)
        d1e = range1(2)
        d2s = range2(1)
        d2e = range2(2)
        d3s = range3(1)
        d3e = range3(2)
        !
        !$cuf kernel do(3)
        DO i3 = d3s, d3e
            !
            DO i2 = d2s, d2e
                !
                DO i1 = d1s, d1e
                    array_out(i1, i2, i3) = val
                END DO
                !
            END DO
            !
        END DO
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_cuf_memset_r3d
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_cuf_memset_c_0dd(array_out, val, range1)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        COMPLEX(DP), INTENT(IN) :: val
        INTEGER, INTENT(IN) :: range1(2)
        !
        COMPLEX(DP), INTENT(INOUT) :: array_out(:)
        !
#if defined(__CUDA)
        attributes(DEVICE) :: array_out
#endif
        !
        INTEGER :: i1, d1s, d1e
        !
        !--------------------------------------------------------------------------------
        !
        d1s = range1(1)
        d1e = range1(2)
        !
        !$cuf kernel do(1)
        DO i1 = d1s, d1e
            array_out(i1) = val
        END DO
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_cuf_memset_c_0dd
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_cuf_memset_c2d(array_out, val, range1, range2)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        COMPLEX(DP), INTENT(IN) :: val
        INTEGER, DIMENSION(2), INTENT(IN) :: range1, range2
        !
        COMPLEX(DP), INTENT(INOUT) :: array_out(:, :)
        !
#if defined(__CUDA)
        attributes(DEVICE) :: array_out
#endif
        !
        INTEGER :: i1, d1s, d1e
        INTEGER :: i2, d2s, d2e
        !
        !--------------------------------------------------------------------------------
        !
        d1s = range1(1)
        d1e = range1(2)
        d2s = range2(1)
        d2e = range2(2)
        !
        !$cuf kernel do(2)
        DO i2 = d2s, d2e
            !
            DO i1 = d1s, d1e
                array_out(i1, i2) = val
            END DO
            !
        END DO
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_cuf_memset_c2d
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_cuf_memset_c3d(array_out, val, range1, range2, range3)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        COMPLEX(DP), INTENT(IN) :: val
        INTEGER, DIMENSION(2), INTENT(IN) :: range1, range2, range3
        !
        COMPLEX(DP), INTENT(INOUT) :: array_out(:, :, :)
        !
#if defined(__CUDA)
        attributes(DEVICE) :: array_out
#endif
        !
        INTEGER :: i1, d1s, d1e
        INTEGER :: i2, d2s, d2e
        INTEGER :: i3, d3s, d3e
        !
        !--------------------------------------------------------------------------------
        !
        d1s = range1(1)
        d1e = range1(2)
        d2s = range2(1)
        d2e = range2(2)
        d3s = range3(1)
        d3e = range3(2)
        !
        !$cuf kernel do(3)
        DO i3 = d3s, d3e
            !
            DO i2 = d2s, d2e
                !
                DO i1 = d1s, d1e
                    array_out(i1, i2, i3) = val
                END DO
                !
            END DO
            !
        END DO
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_cuf_memset_c3d
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE env_utils_cuda
!----------------------------------------------------------------------------------------
