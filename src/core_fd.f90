! Copyright (C) 2018 ENVIRON (www.quantum-environment.org)
! Copyright (C) 2006-2010 Quantum ESPRESSO group
!
!    This file is part of Environ version 1.1
!
!    Environ 1.1 is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 2 of the License, or
!    (at your option) any later version.
!
!    Environ 1.1 is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more detail, either the file
!    `License' in the root directory of the present distribution, or
!    online at <http://www.gnu.org/licenses/>.
!
!----------------------------------------------------------------------------------------
!
! Authors: Oliviero Andreussi (Department of Physics, UNT)
!
!----------------------------------------------------------------------------------------
!>
!! Module to compute finite-differences gradients on dense real space grid
!!
!----------------------------------------------------------------------------------------
MODULE core_fd
    !------------------------------------------------------------------------------------
    !
    USE environ_types
    USE core_types
    USE cell_types
    !
    PRIVATE
    !
    PUBLIC :: gradient_fd
    !
    !------------------------------------------------------------------------------------
CONTAINS
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE gradient_fd(fd, f, grad)
        !--------------------------------------------------------------------------------
        !
        USE scatter_mod, ONLY: scatter_grid
        USE mp, ONLY: mp_sum
        !
        IMPLICIT NONE
        !
        TYPE(fd_core), TARGET, INTENT(IN) :: fd
        TYPE(environ_density), INTENT(IN) :: f
        !
        TYPE(environ_gradient), INTENT(INOUT) :: grad
        !
        INTEGER :: i, ir, ipol, in
        INTEGER, DIMENSION(:), ALLOCATABLE :: ix, iy, iz
        REAL(DP), DIMENSION(:, :), ALLOCATABLE :: gradtmp, gradaux
        !
        LOGICAL :: physical
        INTEGER, POINTER :: nnr, nfdpoint
        TYPE(environ_cell), POINTER :: cell
        TYPE(fft_type_descriptor), POINTER :: dfft
        !
        !--------------------------------------------------------------------------------
        !
        dfft => fd%cell%dfft
        cell => fd%cell
        nnr => fd%cell%nnr
        nfdpoint => fd%nfdpoint
        !
        ALLOCATE (ix(-nfdpoint:nfdpoint), iy(-nfdpoint:nfdpoint), &
                  iz(-nfdpoint:nfdpoint))
        !
        ALLOCATE (gradtmp(dfft%nr1x * dfft%nr2x * dfft%nr3x, 3))
        gradtmp = 0.D0
        !
        DO ir = 1, cell%ir_end
            !
            CALL ir2ijk(cell, ir, ix(0), iy(0), iz(0), physical)
            !
            IF (.NOT. physical) CYCLE ! do not include points outside the physical range
            !
            DO in = 1, fd%nfdpoint
                ix(in) = ix(in - 1) + 1
                !
                IF (ix(in) > dfft%nr1x - 1) ix(in) = 0
                !
                ix(-in) = ix(-in + 1) - 1
                !
                IF (ix(-in) < 0) ix(-in) = dfft%nr1x - 1
                !
                iy(in) = iy(in - 1) + 1
                !
                IF (iy(in) > dfft%nr2x - 1) iy(in) = 0
                !
                iy(-in) = iy(-in + 1) - 1
                !
                IF (iy(-in) < 0) iy(-in) = dfft%nr2x - 1
                !
                iz(in) = iz(in - 1) + 1
                !
                IF (iz(in) > dfft%nr3x - 1) iz(in) = 0
                !
                iz(-in) = iz(-in + 1) - 1
                !
                IF (iz(-in) < 0) iz(-in) = dfft%nr3x - 1
                !
            END DO
            !
            DO in = -fd%nfdpoint, fd%nfdpoint
                i = ix(in) + iy(0) * dfft%nr1x + iz(0) * dfft%nr1x * dfft%nr2x + 1
                gradtmp(i, 1) = gradtmp(i, 1) - fd%icfd(in) * f%of_r(ir) * dfft%nr1
                i = ix(0) + iy(in) * dfft%nr1x + iz(0) * dfft%nr1x * dfft%nr2x + 1
                gradtmp(i, 2) = gradtmp(i, 2) - fd%icfd(in) * f%of_r(ir) * dfft%nr2
                i = ix(0) + iy(0) * dfft%nr1x + iz(in) * dfft%nr1x * dfft%nr2x + 1
                gradtmp(i, 3) = gradtmp(i, 3) - fd%icfd(in) * f%of_r(ir) * dfft%nr3
            END DO
            !
        END DO
        !
        DEALLOCATE (ix, iy, iz)
        !
        ALLOCATE (gradaux(nnr, 3))
        !
#if defined (__MPI)
        DO ipol = 1, 3
            !
            CALL mp_sum(gradtmp(:, ipol), cell%dfft%comm)
            !
            CALL scatter_grid(dfft, gradtmp(:, ipol), gradaux(:, ipol))
            !
        END DO
        !
#else
        gradaux(1:nnr, :) = gradtmp(1:nnr, :)
#endif
        !
        DEALLOCATE (gradtmp)
        !
        DO ir = 1, nnr
            grad%of_r(:, ir) = MATMUL(cell%bg, gradaux(ir, :))
        END DO
        !
        DEALLOCATE (gradaux)
        !
        grad%of_r = grad%of_r / DBLE(fd%ncfd) / cell%alat
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE gradient_fd
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE core_fd
!----------------------------------------------------------------------------------------
