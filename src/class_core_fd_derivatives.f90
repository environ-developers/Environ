!----------------------------------------------------------------------------------------
!
! Copyright (C) 2018-2021 ENVIRON (www.quantum-environ.org)
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
! Authors: Oliviero Andreussi (Department of Physics, UNT)
!          Francesco Nattino  (THEOS and NCCR-MARVEL, EPFL)
!          Nicola Marzari     (THEOS and NCCR-MARVEL, EPFL)
!          Edan Bainglass     (Department of Physics, UNT)
!
!----------------------------------------------------------------------------------------
!>
!!
!----------------------------------------------------------------------------------------
MODULE class_core_fd_derivatives
    !------------------------------------------------------------------------------------
    !
    USE env_mp, ONLY: env_mp_sum
    !
    USE env_types_fft, ONLY: env_fft_type_descriptor
    USE env_base_scatter, ONLY: env_scatter_grid
    !
    USE environ_param, ONLY: DP
    !
    USE class_cell
    USE class_density
    USE class_gradient
    !
    USE class_core_fd
    !
    !------------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    PRIVATE
    !
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    TYPE, EXTENDS(core_fd), PUBLIC :: core_fd_derivatives
        !--------------------------------------------------------------------------------
        !
        !--------------------------------------------------------------------------------
    CONTAINS
        !--------------------------------------------------------------------------------
        !
        PROCEDURE :: gradient => gradient_fd
        !
        !--------------------------------------------------------------------------------
    END TYPE core_fd_derivatives
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
CONTAINS
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !                                  GENERAL METHODS
    !
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE gradient_fd(this, f, grad)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(core_fd_derivatives), TARGET, INTENT(IN) :: this
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
        TYPE(env_fft_type_descriptor), POINTER :: dfft
        !
        !--------------------------------------------------------------------------------
        !
        dfft => this%cell%dfft
        cell => this%cell
        nnr => this%cell%nnr
        nfdpoint => this%nfdpoint
        !
        ALLOCATE (ix(-nfdpoint:nfdpoint), iy(-nfdpoint:nfdpoint), &
                  iz(-nfdpoint:nfdpoint))
        !
        ALLOCATE (gradtmp(dfft%nr1x * dfft%nr2x * dfft%nr3x, 3))
        gradtmp = 0.D0
        !
        DO ir = 1, cell%ir_end
            !
            CALL cell%ir2ijk(ir, ix(0), iy(0), iz(0), physical)
            !
            IF (.NOT. physical) CYCLE ! do not include points outside the physical range
            !
            DO in = 1, this%nfdpoint
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
            DO in = -this%nfdpoint, this%nfdpoint
                i = ix(in) + iy(0) * dfft%nr1x + iz(0) * dfft%nr1x * dfft%nr2x + 1
                gradtmp(i, 1) = gradtmp(i, 1) - this%icfd(in) * f%of_r(ir) * dfft%nr1
                i = ix(0) + iy(in) * dfft%nr1x + iz(0) * dfft%nr1x * dfft%nr2x + 1
                gradtmp(i, 2) = gradtmp(i, 2) - this%icfd(in) * f%of_r(ir) * dfft%nr2
                i = ix(0) + iy(0) * dfft%nr1x + iz(in) * dfft%nr1x * dfft%nr2x + 1
                gradtmp(i, 3) = gradtmp(i, 3) - this%icfd(in) * f%of_r(ir) * dfft%nr3
            END DO
            !
        END DO
        !
        DEALLOCATE (ix, iy, iz)
        !
        ALLOCATE (gradaux(nnr, 3))
        !
#if defined(__MPI)
        DO ipol = 1, 3
            !
            CALL env_mp_sum(gradtmp(:, ipol), cell%dfft%comm)
            !
            CALL env_scatter_grid(dfft, gradtmp(:, ipol), gradaux(:, ipol))
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
        grad%of_r = grad%of_r / DBLE(this%ncfd)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE gradient_fd
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE class_core_fd_derivatives
!----------------------------------------------------------------------------------------
