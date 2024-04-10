!----------------------------------------------------------------------------------------
!
! Copyright (C) 2018-2022 ENVIRON (www.quantum-environ.org)
!
!----------------------------------------------------------------------------------------
!
!     This file is part of Environ version 3.0
!
!     Environ 3.0 is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation, either version 2 of the License, or
!     (at your option) any later version.
!
!     Environ 3.0 is distributed in the hope that it will be useful,
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
!! Environ_boundary contains all the specifications and the details of
!! the smooth interface between the QM and the continuum regions of the
!! simulation cell. The main interface function is stored in the %scaled
!! component, the type also stores boundary real-space derivatives (gradient,
!! laplacian, dsurface, hessian) and other quantities needed by Environ
!! modules.
!!
!----------------------------------------------------------------------------------------
MODULE boundary_tools
    !------------------------------------------------------------------------------------
    !
    USE class_io, ONLY: io
    !
    USE environ_param, ONLY: DP, sqrtpi, tpi
    !
    USE class_cell
    USE class_density
    USE class_gradient
    USE class_hessian
    !
    USE tools_math, ONLY: environ_erfc
    !
    !------------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    PRIVATE
    !
    PUBLIC :: calc_dsurface_no_pre
    !
    PUBLIC :: calc_partial_of_boundary
    !
    PUBLIC :: gradient_of_boundary, laplacian_of_boundary, dsurface_of_boundary
    !
    PUBLIC :: sfunct0, dsfunct0, &
              sfunct1, dsfunct1, d2sfunct1, &
              sfunct2, dsfunct2, d2sfunct2
    !
    PUBLIC :: reduced2nnr
    !
    !------------------------------------------------------------------------------------
    !
    INTERFACE gradient_of_boundary
        MODULE PROCEDURE &
            calc_gradient_of_boundary_highmem, &
            calc_gradient_of_boundary_lowmem
    END INTERFACE gradient_of_boundary
    !
    INTERFACE laplacian_of_boundary
        MODULE PROCEDURE &
            calc_laplacian_of_boundary_highmem, &
            calc_laplacian_of_boundary_lowmem
    END INTERFACE laplacian_of_boundary
    !
    INTERFACE dsurface_of_boundary
        MODULE PROCEDURE &
            calc_dsurface_of_boundary_highmem, &
            calc_dsurface_of_boundary_lowmem
    END INTERFACE dsurface_of_boundary
    !
    !------------------------------------------------------------------------------------
    !
    REAL(DP), PARAMETER :: bound_tol = 1.D-60
    !
    !------------------------------------------------------------------------------------
CONTAINS
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !                                BOUNDARY DERIVATIVES
    !
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_dsurface_no_pre(cell, grad, hess, dsurf)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_cell), INTENT(IN) :: cell
        REAL(DP), INTENT(IN) :: grad(3, cell%nnr)
        REAL(DP), INTENT(IN) :: hess(3, 3, cell%nnr)
        !
        REAL(DP), INTENT(OUT) :: dsurf(cell%nnr)
        !
        REAL(DP), PARAMETER :: toldsurface = 1.D-50
        !
        INTEGER :: j, k, i
        REAL(DP) :: gmod
        !
        !--------------------------------------------------------------------------------
        !
        DO i = 1, cell%ir_end
            dsurf(i) = 0.D0
            gmod = SUM(grad(:, i)**2)
            !
            IF (gmod < toldsurface) CYCLE
            !
            DO j = 1, 3
                !
                DO k = 1, 3
                    !
                    IF (j == k) CYCLE
                    !
                    dsurf(i) = dsurf(i) + &
                               grad(j, i) * grad(k, i) * hess(j, k, i) - &
                               grad(j, i) * grad(j, i) * hess(k, k, i)
                    !
                END DO
                !
            END DO
            !
            dsurf(i) = dsurf(i) / gmod / SQRT(gmod)
        END DO
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_dsurface_no_pre
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_partial_of_boundary(n, i, local, gradlocal, partial)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: n, i
        TYPE(environ_density), INTENT(IN) :: local(n)
        TYPE(environ_gradient), INTENT(IN) :: gradlocal(n)
        !
        TYPE(environ_gradient), INTENT(INOUT) :: partial
        !
        INTEGER :: j, k
        !
        TYPE(environ_density) :: denlocal
        !
        CHARACTER(LEN=80) :: routine = 'calc_partial_of_boundary'
        !
        !--------------------------------------------------------------------------------
        !
        IF (i > n) CALL io%error(routine, "Index out of bound", 1)
        !
        DO j = 1, 3
            partial%of_r(j, :) = gradlocal(i)%of_r(j, :)
            !
            DO k = 1, n
                !
                IF (k == i) CYCLE
                !
                partial%of_r(j, :) = partial%of_r(j, :) * local(k)%of_r(:)
            END DO
            !
        END DO
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_partial_of_boundary
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_gradient_of_boundary_highmem(n, local, gradlocal, grad)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: n
        TYPE(environ_density), INTENT(IN) :: local(n)
        TYPE(environ_gradient), INTENT(IN) :: gradlocal(n)
        !
        TYPE(environ_gradient), INTENT(INOUT) :: grad
        !
        INTEGER :: i
        TYPE(environ_cell), POINTER :: cell
        TYPE(environ_gradient) :: partial
        !
        !--------------------------------------------------------------------------------
        !
        cell => grad%cell
        !
        CALL partial%init(cell)
        !
        grad%of_r = 0.D0
        !
        DO i = 1, n
            !
            CALL calc_partial_of_boundary(n, i, local, gradlocal, partial)
            !
            grad%of_r = grad%of_r + partial%of_r
        END DO
        !
        CALL partial%destroy()
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_gradient_of_boundary_highmem
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_laplacian_of_boundary_highmem(n, local, gradloc, laplloc, lapl)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: n
        TYPE(environ_density), INTENT(IN) :: local(n)
        TYPE(environ_gradient), INTENT(IN) :: gradloc(n)
        TYPE(environ_density), INTENT(IN) :: laplloc(n)
        !
        TYPE(environ_density), INTENT(INOUT) :: lapl
        !
        INTEGER :: i, j, k
        TYPE(environ_cell), POINTER :: cell
        TYPE(environ_density) :: tmp
        !
        !--------------------------------------------------------------------------------
        !
        cell => lapl%cell
        !
        CALL tmp%init(cell)
        !
        lapl%of_r = 0.D0
        !
        DO i = 1, n
            !
            DO j = 1, n
                !
                IF (j == i) THEN
                    tmp%of_r = laplloc(i)%of_r
                ELSE
                    CALL gradloc(i)%scalar_product(gradloc(j), tmp)
                END IF
                !
                DO k = 1, n
                    !
                    IF (k == j .OR. k == i) CYCLE
                    !
                    tmp%of_r = tmp%of_r * local(k)%of_r
                END DO
                !
                lapl%of_r = lapl%of_r + tmp%of_r
            END DO
            !
        END DO
        !
        CALL tmp%destroy()
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_laplacian_of_boundary_highmem
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_dsurface_of_boundary_highmem(n, local, gradloc, hessloc, grad, hess, dsurf)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: n
        TYPE(environ_density), INTENT(IN) :: local(n)
        TYPE(environ_gradient), INTENT(IN) :: gradloc(n)
        TYPE(environ_hessian), INTENT(IN) :: hessloc(n)
        !
        TYPE(environ_gradient), INTENT(INOUT) :: grad
        TYPE(environ_density), INTENT(INOUT) :: dsurf
        TYPE(environ_hessian), INTENT(INOUT) :: hess
        !
        INTEGER :: i, j, k, l, m
        TYPE(environ_cell), POINTER :: cell
        TYPE(environ_density) :: dens
        TYPE(environ_gradient) :: partial
        !
        !--------------------------------------------------------------------------------
        !
        cell => grad%cell
        !
        CALL dens%init(cell)
        !
        CALL partial%init(cell)
        !
        grad%of_r = 0.D0
        !
        DO i = 1, n
            !
            CALL calc_partial_of_boundary(n, i, local, gradloc, partial)
            !
            grad%of_r = grad%of_r + partial%of_r
            !
            DO j = 1, n
                !
                DO k = 1, 3
                    !
                    DO l = 1, 3
                        !
                        IF (j == i) THEN
                            dens%of_r = hessloc(i)%of_r(k, l, :)
                        ELSE
                            dens%of_r(:) = gradloc(i)%of_r(k, :) * &
                                           gradloc(j)%of_r(l, :)
                        END IF
                        !
                        DO m = 1, n
                            !
                            IF (m == j .OR. m == i) CYCLE
                            !
                            dens%of_r = dens%of_r * local(m)%of_r
                        END DO
                        !
                        hess%of_r(k, l, :) = hess%of_r(k, l, :) + dens%of_r
                        !
                    END DO
                    !
                END DO
                !
            END DO
            !
        END DO
        !
        !--------------------------------------------------------------------------------
        ! Final operations
        !
        CALL calc_dsurface_no_pre(cell, grad%of_r, hess%of_r, dsurf%of_r)
        !
        CALL dens%destroy()
        !
        CALL partial%destroy()
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_dsurface_of_boundary_highmem
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_gradient_of_boundary_lowmem(n, local, gradlocal, scal, grad)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: n
        TYPE(environ_density), INTENT(IN) :: scal ! soft sphere interface function
        TYPE(environ_density), INTENT(IN) :: local(n)
        TYPE(environ_gradient), INTENT(IN) :: gradlocal(n)
        !
        TYPE(environ_gradient), INTENT(INOUT) :: grad
        !
        INTEGER :: i, j, k
        TYPE(environ_cell), POINTER :: cell
        !
        !--------------------------------------------------------------------------------
        !
        cell => grad%cell
        !
        grad%of_r = 0.D0
        !
        !--------------------------------------------------------------------------------
        ! Temporary quotient
        !
        DO i = 1, n
            !
            DO j = 1, cell%nnr
                !
                IF (ABS(local(i)%of_r(j)) <= bound_tol) CYCLE
                !
                DO k = 1, 3
                    !
                    grad%of_r(k, j) = &
                        grad%of_r(k, j) + &
                        (gradlocal(i)%of_r(k, j) / &
                         local(i)%of_r(j) * scal%of_r(j))
                    !
                END DO
                !
            END DO
            !
        END DO
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_gradient_of_boundary_lowmem
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_laplacian_of_boundary_lowmem(n, local, gradlocal, laploc, scal, grad, lapl)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: n
        TYPE(environ_density), INTENT(IN) :: scal ! soft sphere interface function
        TYPE(environ_density), INTENT(IN) :: local(n)
        TYPE(environ_gradient), INTENT(IN) :: gradlocal(n)
        TYPE(environ_density), INTENT(IN) :: laploc(n)
        TYPE(environ_gradient), INTENT(IN) :: grad
        !
        TYPE(environ_density), INTENT(INOUT) :: lapl
        !
        INTEGER :: i, j, k
        TYPE(environ_cell), POINTER :: cell
        !
        !--------------------------------------------------------------------------------
        !
        cell => lapl%cell
        !
        DO i = 1, n
            !
            DO j = 1, cell%nnr
                !
                IF (ABS(local(i)%of_r(j)) <= bound_tol) CYCLE
                !
                lapl%of_r(j) = lapl%of_r(j) + &
                               (laploc(i)%of_r(j) / &
                               local(i)%of_r(j) * scal%of_r(j))
                !
                DO k = 1, 3
                    !
                    lapl%of_r(j) = &
                        lapl%of_r(j) - &
                        ((gradlocal(i)%of_r(k, j)**2 / &
                        local(i)%of_r(j)**2) * scal%of_r(j))
                    !
                    lapl%of_r(j) = &
                        lapl%of_r(j) + &
                        (grad%of_r(k, j) * &
                        gradlocal(i)%of_r(k, j) / local(i)%of_r(j))
                    !
                END DO
                !
            END DO
            !
        END DO
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_laplacian_of_boundary_lowmem
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_dsurface_of_boundary_lowmem(n, local, gradloc, hessloc, grad, hess, scal, dsurf)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: n
        TYPE(environ_density), INTENT(IN) :: scal
        TYPE(environ_density), INTENT(IN) :: local(n)
        TYPE(environ_gradient), INTENT(IN) :: gradloc(n)
        TYPE(environ_hessian), INTENT(IN) :: hessloc(n)
        TYPE(environ_gradient), INTENT(IN) :: grad
        !
        TYPE(environ_density), INTENT(INOUT) :: dsurf
        TYPE(environ_hessian), INTENT(INOUT) :: hess
        !
        INTEGER :: i, j, k, l
        TYPE(environ_cell), POINTER :: cell
        !
        !--------------------------------------------------------------------------------
        !
        cell => grad%cell
        !
        DO i = 1, n
            !
            DO j = 1, cell%nnr
                !
                IF (ABS(local(i)%of_r(j)) <= bound_tol) CYCLE
                !
                DO k = 1, 3
                    !
                    DO l = 1, 3
                        !
                        hess%of_r(k, l, j) = &
                            hess%of_r(k, l, j) + &
                            (hessloc(i)%of_r(k, l, j) / local(i)%of_r(j) * scal%of_r(j))
                        !
                        hess%of_r(k, l, j) = &
                            hess%of_r(k, l, j) - &
                            ((gradloc(i)%of_r(k, j) * gradloc(i)%of_r(l, j) / &
                              local(i)%of_r(j)**2) * scal%of_r(j))
                        !
                        hess%of_r(k, l, j) = &
                            hess%of_r(k, l, j) + &
                            (grad%of_r(k, j) * gradloc(i)%of_r(l, j) / &
                            local(i)%of_r(j))
                        !
                    END DO
                    !
                END DO
                !
            END DO
            !
        END DO
        !
        CALL calc_dsurface_no_pre(cell, grad%of_r, hess%of_r, dsurf%of_r)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_dsurface_of_boundary_lowmem
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !                                BOUNDARY GENERATORS
    !
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !>
    !! Switching function 0: goes from 1 to 0 when passing through the
    !! threshold
    !!
    !! \f[
    !!    1 + \frac{1 - (x/x_t)^k}{1 + (x/x_t)^k}
    !! \f]
    !! where \f$x_t\f$ is the threshold
    !!
    !------------------------------------------------------------------------------------
    REAL(DP) FUNCTION sfunct0(x, xthr, fact)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP) :: x, xthr, fact
        !
        REAL(DP) :: arg
        !
        !--------------------------------------------------------------------------------
        !
        arg = (ABS(x) / xthr)**fact
        sfunct0 = 0.5D0 * (1.D0 + (1.D0 - arg) / (1.D0 + arg))
        !
        !--------------------------------------------------------------------------------
    END FUNCTION sfunct0
    !------------------------------------------------------------------------------------
    !>
    !! Derivative of switching function 0
    !!
    !------------------------------------------------------------------------------------
    REAL(DP) FUNCTION dsfunct0(x, xthr, fact)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP) :: x, xthr, fact
        !
        REAL(DP) :: arg
        !
        !--------------------------------------------------------------------------------
        !
        arg = (ABS(x) / xthr)**fact
        dsfunct0 = -fact * ABS(x)**(fact - 1.D0) / xthr**fact / (1.D0 + arg)**2
        !
        !--------------------------------------------------------------------------------
    END FUNCTION dsfunct0
    !------------------------------------------------------------------------------------
    !>
    !! Switching function 1 that goes from 1 to 0 when passing from
    !! xmin to xmax.
    !!
    !! NOTE: fact should be equal to LOG(xmax/xmin) but is
    !! passed in input to save time
    !!
    !! \f[
    !!    x - \sin(x)
    !! \f]
    !!
    !------------------------------------------------------------------------------------
    REAL(DP) FUNCTION sfunct1(x, xmax, xmin, fact)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP) :: x, xmax, xmin, fact
        !
        REAL(DP) :: arg
        !
        !--------------------------------------------------------------------------------
        !
        IF (x <= xmin) THEN
            sfunct1 = 1.D0
        ELSE IF (x < xmax) THEN
            arg = tpi * LOG(xmax / ABS(x)) / fact
            sfunct1 = (arg - SIN(arg)) / tpi
        ELSE
            sfunct1 = 0.D0
        END IF
        !
        !--------------------------------------------------------------------------------
    END FUNCTION sfunct1
    !------------------------------------------------------------------------------------
    !>
    !! @brief Derivative of switching function 1
    !!
    !! NOTE: fact should be equal to LOG(xmax/xmin) but is passed in
    !! input to save time.
    !!
    !------------------------------------------------------------------------------------
    REAL(DP) FUNCTION dsfunct1(x, xmax, xmin, fact)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP) :: x, xmax, xmin, fact
        !
        REAL(DP) :: arg
        !
        !--------------------------------------------------------------------------------
        !
        IF (x <= xmin) THEN
            dsfunct1 = 0.D0
        ELSE IF (x < xmax) THEN
            arg = tpi * LOG(xmax / ABS(x)) / fact
            dsfunct1 = (COS(arg) - 1.D0) / ABS(x) / fact ! #TODO in fact should not use ABS(x)
        ELSE
            dsfunct1 = 0.D0
        END IF
        !
        !--------------------------------------------------------------------------------
    END FUNCTION dsfunct1
    !------------------------------------------------------------------------------------
    !>
    !! @brief Second derivative of switching function 1
    !!
    !! Note: fact should be equal to LOG(xmax/xmin) but is passed in
    !! input to save time
    !!
    !------------------------------------------------------------------------------------
    REAL(DP) FUNCTION d2sfunct1(x, xmax, xmin, fact)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP) :: x, xmax, xmin, fact
        !
        REAL(DP) :: arg
        !
        !--------------------------------------------------------------------------------
        !
        IF (x <= xmin) THEN
            d2sfunct1 = 0.D0
        ELSE IF (x < xmax) THEN
            arg = tpi * LOG(xmax / ABS(x)) / fact
            d2sfunct1 = (tpi * SIN(arg) + fact * (1.D0 - COS(arg))) / (x * fact)**2
        ELSE
            d2sfunct1 = 0.D0
        END IF
        !
        !--------------------------------------------------------------------------------
    END FUNCTION d2sfunct1
    !------------------------------------------------------------------------------------
    !>
    !! Switching function 2, erfc() that goes from 1 to 0 when passing
    !! through xthr.
    !!
    !------------------------------------------------------------------------------------
    REAL(DP) FUNCTION sfunct2(x, xthr, spread)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP) :: x, xthr, spread
        !
        REAL(DP) :: arg
        !
        !--------------------------------------------------------------------------------
        !
        arg = (x - xthr) / spread
        sfunct2 = 0.5D0 * environ_erfc(arg)
        !
        !--------------------------------------------------------------------------------
    END FUNCTION sfunct2
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    REAL(DP) FUNCTION dsfunct2(x, xthr, spread)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP) :: x, xthr, spread
        !
        REAL(DP) :: arg
        !
        !--------------------------------------------------------------------------------
        !
        arg = (x - xthr) / spread
        !
        IF (ABS(arg) > 6.D0) THEN ! 6.D0 is the threshold of environ_erfc(x)
            dsfunct2 = 0.D0
        ELSE
            dsfunct2 = -EXP(-arg**2) / sqrtpi / spread
        END IF
        !
        !--------------------------------------------------------------------------------
    END FUNCTION dsfunct2
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    REAL(DP) FUNCTION d2sfunct2(x, xthr, spread)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP) :: x, xthr, spread
        !
        REAL(DP) :: arg
        !
        !--------------------------------------------------------------------------------
        !
        arg = (x - xthr) / spread
        IF (ABS(arg) > 6.D0) THEN
            d2sfunct2 = 0.D0
        ELSE
            d2sfunct2 = EXP(-arg**2) / sqrtpi / spread**2 * 2.D0 * arg
        END IF
        !
        !--------------------------------------------------------------------------------
    END FUNCTION d2sfunct2
    !------------------------------------------------------------------------------------
    !>
    !! @brief Maps reduced arrays to size of nnr
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE reduced2nnr(ir, nnr, init_val, den_of_r, den_vals, grad_of_r, grad_vals)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: ir(:), nnr
        REAL(DP), INTENT(IN) :: init_val
        REAL(DP), OPTIONAL, INTENT(IN) :: den_vals(:), grad_vals(:, :)
        !
        REAL(DP), OPTIONAL, INTENT(INOUT) :: den_of_r(:), grad_of_r(:, :)
        !
        INTEGER :: i, idx
        !
        CHARACTER(LEN=80) :: routine = 'reduced2nnr'
        !
        !--------------------------------------------------------------------------------
        !
        idx = 1
        !
        IF (PRESENT(den_of_r)) THEN
            den_of_r = init_val
            !
            IF (.NOT. PRESENT(den_vals)) &
                CALL io%error(routine, 'Missing stored density values', 1)
            !
        ELSE IF (PRESENT(grad_of_r)) THEN
            grad_of_r = init_val
            !
            IF (.NOT. PRESENT(grad_vals)) &
                CALL io%error(routine, 'Missing stored gradient values', 1)
            !
        ELSE
            CALL io%error(routine, 'Missing of_r array', 1)
        END IF
        !
        !--------------------------------------------------------------------------------
        ! Map array to nnr
        !
        DO i = 1, nnr
            !
            IF (.NOT. ir(idx) == i) CYCLE
            !
            IF (PRESENT(den_of_r)) THEN
                den_of_r(i) = den_vals(idx)
            ELSE IF (PRESENT(grad_of_r)) THEN
                grad_of_r(:, i) = grad_vals(idx, :)
            END IF
            !
            idx = idx + 1
        END DO
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE reduced2nnr
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE boundary_tools
!----------------------------------------------------------------------------------------
