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
              sfunct2, dsfunct2, d2sfunct2, &
              boundfunct, dboundfunct, d2boundfunct
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
    INTEGER, PARAMETER :: bound_tol = 1.D-60
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
    SUBROUTINE calc_dsurface_no_pre(cell, grad, hess, dsurface)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_cell) :: cell
        REAL(DP), INTENT(IN) :: grad(3, cell%nnr)
        REAL(DP), INTENT(IN) :: hess(3, 3, cell%nnr)
        !
        REAL(DP), INTENT(OUT) :: dsurface(cell%nnr)
        !
        REAL(DP), PARAMETER :: toldsurface = 1.D-50
        !
        INTEGER :: ipol, jpol, i
        REAL(DP) :: gmod
        !
        !
        !--------------------------------------------------------------------------------
        !
        DO i = 1, cell%ir_end
            dsurface(i) = 0.D0
            gmod = SUM(grad(:, i)**2)
            !
            IF (gmod < toldsurface) CYCLE
            !
            DO ipol = 1, 3
                !
                DO jpol = 1, 3
                    !
                    IF (ipol == jpol) CYCLE
                    !
                    dsurface(i) = dsurface(i) + &
                                  grad(ipol, i) * grad(jpol, i) * hess(ipol, jpol, i) - &
                                  grad(ipol, i) * grad(ipol, i) * hess(jpol, jpol, i)
                    !
                END DO
                !
            END DO
            !
            dsurface(i) = dsurface(i) / gmod / SQRT(gmod)
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
        INTEGER :: j, ipol
        !
        CHARACTER(LEN=80) :: sub_name = 'calc_partial_of_boundary'
        !
        !--------------------------------------------------------------------------------
        !
        IF (i > n) CALL io%error(sub_name, 'Index out of bound', 1)
        !
        DO ipol = 1, 3
            partial%of_r(ipol, :) = gradlocal(i)%of_r(ipol, :)
            !
            DO j = 1, n
                !
                IF (j == i) CYCLE
                !
                partial%of_r(ipol, :) = partial%of_r(ipol, :) * local(j)%of_r(:)
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
    SUBROUTINE calc_gradient_of_boundary_highmem(n, local, gradlocal, gradient)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: n
        TYPE(environ_density), INTENT(IN) :: local(n)
        TYPE(environ_gradient), INTENT(IN) :: gradlocal(n)
        !
        TYPE(environ_gradient), INTENT(INOUT) :: gradient
        !
        INTEGER :: i, j, ipol
        TYPE(environ_cell), POINTER :: cell
        TYPE(environ_gradient) :: partial
        !
        !--------------------------------------------------------------------------------
        !
        cell => gradient%cell
        !
        CALL partial%init(cell)
        !
        gradient%of_r = 0.D0
        !
        DO i = 1, n
            !
            CALL calc_partial_of_boundary(n, i, local, gradlocal, partial)
            !
            gradient%of_r = gradient%of_r + partial%of_r
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
    SUBROUTINE calc_laplacian_of_boundary_highmem(n, local, gradlocal, lapllocal, &
                                                  laplacian)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: n
        TYPE(environ_density), INTENT(IN) :: local(n)
        TYPE(environ_gradient), INTENT(IN) :: gradlocal(n)
        TYPE(environ_density), INTENT(IN) :: lapllocal(n)
        !
        TYPE(environ_density), INTENT(INOUT) :: laplacian
        !
        INTEGER :: i, j, k, ipol
        TYPE(environ_cell), POINTER :: cell
        TYPE(environ_density) :: tmp
        !
        !--------------------------------------------------------------------------------
        !
        cell => laplacian%cell
        !
        CALL tmp%init(cell)
        !
        laplacian%of_r = 0.D0
        !
        DO i = 1, n
            !
            DO j = 1, n
                !
                IF (j == i) THEN
                    tmp%of_r = lapllocal(i)%of_r
                ELSE
                    CALL gradlocal(i)%scalar_product(gradlocal(j), tmp)
                END IF
                !
                DO k = 1, n
                    !
                    IF (k == j .OR. k == i) CYCLE
                    !
                    tmp%of_r = tmp%of_r * local(k)%of_r
                END DO
                !
                laplacian%of_r = laplacian%of_r + tmp%of_r
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
    SUBROUTINE calc_dsurface_of_boundary_highmem(n, local, gradlocal, hesslocal, &
                                                 gradient, laplacian, hessian, dsurface)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: n
        TYPE(environ_density), INTENT(IN) :: local(n)
        TYPE(environ_gradient), INTENT(IN) :: gradlocal(n)
        TYPE(environ_hessian), INTENT(IN) :: hesslocal(n)
        !
        TYPE(environ_gradient), INTENT(INOUT) :: gradient
        TYPE(environ_density), INTENT(INOUT) :: laplacian, dsurface
        TYPE(environ_hessian), INTENT(INOUT) :: hessian
        !
        INTEGER :: i, j, k, ipol, jpol
        TYPE(environ_cell), POINTER :: cell
        TYPE(environ_density) :: dens
        TYPE(environ_gradient) :: partial
        !
        !--------------------------------------------------------------------------------
        !
        cell => laplacian%cell
        !
        CALL dens%init(cell)
        !
        CALL partial%init(cell)
        !
        gradient%of_r = 0.D0
        !
        DO i = 1, n
            !
            CALL calc_partial_of_boundary(n, i, local, gradlocal, partial)
            !
            gradient%of_r = gradient%of_r + partial%of_r
            !
            DO j = 1, n
                !
                DO ipol = 1, 3
                    !
                    DO jpol = 1, 3
                        !
                        IF (j == i) THEN
                            dens%of_r(:) = hesslocal(i)%of_r(ipol, jpol, :)
                        ELSE
                            !
                            dens%of_r(:) = gradlocal(i)%of_r(ipol, :) * &
                                           gradlocal(j)%of_r(jpol, :)
                            !
                        END IF
                        !
                        DO k = 1, n
                            !
                            IF (k == j .OR. k == i) CYCLE
                            !
                            dens%of_r = dens%of_r * local(k)%of_r
                        END DO
                        !
                        hessian%of_r(ipol, jpol, :) = hessian%of_r(ipol, jpol, :) + &
                                                      dens%of_r(:)
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
        laplacian%of_r = hessian%of_r(1, 1, :) + hessian%of_r(2, 2, :) + &
                         hessian%of_r(3, 3, :)
        !
        CALL calc_dsurface_no_pre(cell, gradient%of_r, hessian%of_r, dsurface%of_r)
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
    SUBROUTINE calc_gradient_of_boundary_lowmem(n, local, gradlocal, scaled, gradient)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: n
        TYPE(environ_density), INTENT(IN) :: scaled ! soft sphere interface function
        TYPE(environ_density), INTENT(IN) :: local(n)
        TYPE(environ_gradient), INTENT(IN) :: gradlocal(n)
        !
        TYPE(environ_gradient), INTENT(INOUT) :: gradient
        !
        INTEGER :: i, j, ipol
        TYPE(environ_cell), POINTER :: cell
        !
        !--------------------------------------------------------------------------------
        !
        cell => gradient%cell
        !
        gradient%of_r = 0.D0
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
                DO ipol = 1, 3
                    gradient%of_r(ipol, j) = gradient%of_r(ipol, j) + &
                                             (gradlocal(i)%of_r(ipol, j) / &
                                              local(i)%of_r(j) * scaled%of_r(j))
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
    SUBROUTINE calc_laplacian_of_boundary_lowmem(n, local, gradlocal, lapllocal, &
                                                 scaled, gradient, laplacian)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: n
        TYPE(environ_density), INTENT(IN) :: scaled ! soft sphere interface function
        TYPE(environ_density), INTENT(IN) :: local(n)
        TYPE(environ_gradient), INTENT(IN) :: gradlocal(n)
        TYPE(environ_density), INTENT(IN) :: lapllocal(n)
        TYPE(environ_gradient), INTENT(IN) :: gradient
        !
        TYPE(environ_density), INTENT(INOUT) :: laplacian
        !
        INTEGER :: i, j, k, ipol
        TYPE(environ_cell), POINTER :: cell
        !
        !--------------------------------------------------------------------------------
        !
        cell => laplacian%cell
        !
        DO i = 1, n
            !
            DO j = 1, cell%nnr
                !
                IF (ABS(local(i)%of_r(j)) <= bound_tol) CYCLE
                !
                laplacian%of_r(j) = laplacian%of_r(j) + &
                                    (lapllocal(i)%of_r(j) / &
                                     local(i)%of_r(j) * scaled%of_r(j))
                !
                DO ipol = 1, 3
                    !
                    laplacian%of_r(j) = laplacian%of_r(j) - &
                                        ((gradlocal(i)%of_r(ipol, j)**2 / &
                                          local(i)%of_r(j)**2) * scaled%of_r(j))
                    !
                    laplacian%of_r(j) = laplacian%of_r(j) + &
                                        (gradient%of_r(ipol, j) * &
                                         gradlocal(i)%of_r(ipol, j) / local(i)%of_r(j))
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
    SUBROUTINE calc_dsurface_of_boundary_lowmem(n, local, gradlocal, hesslocal, &
                                                gradient, laplacian, hessian, &
                                                scaled, dsurface)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: n
        TYPE(environ_density), INTENT(IN) :: scaled
        TYPE(environ_density), INTENT(IN) :: local(n)
        TYPE(environ_gradient), INTENT(IN) :: gradlocal(n)
        TYPE(environ_hessian), INTENT(IN) :: hesslocal(n)
        TYPE(environ_gradient), INTENT(IN) :: gradient
        !
        TYPE(environ_density), INTENT(INOUT) :: laplacian
        TYPE(environ_density), INTENT(INOUT) :: dsurface
        TYPE(environ_hessian), INTENT(INOUT) :: hessian
        !
        INTEGER :: i, j, k, ipol, jpol
        TYPE(environ_cell), POINTER :: cell
        !
        !--------------------------------------------------------------------------------
        !
        cell => laplacian%cell
        !
        DO i = 1, n
            !
            DO j = 1, cell%nnr
                !
                IF (ABS(local(i)%of_r(j)) <= bound_tol) CYCLE
                !
                DO ipol = 1, 3
                    !
                    DO jpol = 1, 3
                        !
                        hessian%of_r(ipol, jpol, j) = &
                            hessian%of_r(ipol, jpol, j) + &
                            (hesslocal(i)%of_r(ipol, jpol, j) / &
                             local(i)%of_r(j) * scaled%of_r(j))
                        !
                        hessian%of_r(ipol, jpol, j) = &
                            hessian%of_r(ipol, jpol, j) - &
                            ((gradlocal(i)%of_r(ipol, j) * gradlocal(i)%of_r(jpol, j) / &
                              local(i)%of_r(j)**2) * scaled%of_r(j))
                        !
                        hessian%of_r(ipol, jpol, j) = &
                            hessian%of_r(ipol, jpol, j) + &
                            (gradient%of_r(ipol, j) * gradlocal(i)%of_r(jpol, j) / &
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
        laplacian%of_r = hessian%of_r(1, 1, :) + hessian%of_r(2, 2, :) + &
                         hessian%of_r(3, 3, :)
        !
        CALL calc_dsurface_no_pre(cell, gradient%of_r, hessian%of_r, dsurface%of_r)
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
    !! Calculates the density-dependent dielectric constant
    !!
    !! ifunct = 0 => original Fattebert and Gygi function
    !!
    !------------------------------------------------------------------------------------
    REAL(DP) FUNCTION boundfunct(rho, rhomax, rhomin, tbeta, const, ifunct)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP) :: rho
        REAL(DP) :: rhomax
        REAL(DP) :: rhomin
        REAL(DP) :: tbeta
        REAL(DP) :: const
        !
        INTEGER :: ifunct
        !
        REAL(DP) :: arg
        !
        CHARACTER(LEN=80) :: fun_name = 'boundfunct'
        !
        !--------------------------------------------------------------------------------
        !
        SELECT CASE (ifunct)
            !
        CASE (0)
            boundfunct = 1.D0 - sfunct0(rho, rhomax, tbeta)
            !
        CASE (1)
            boundfunct = 1.D0 - sfunct1(rho, rhomax, rhomin, tbeta)
            !
        CASE (2)
            !
            boundfunct = &
                (const - EXP(LOG(const) * sfunct1(rho, rhomax, rhomin, tbeta))) / &
                (const - 1.D0)
            !
        CASE DEFAULT
            CALL io%error(fun_name, 'Unknown boundary type', 1)
            !
        END SELECT
        !
        !--------------------------------------------------------------------------------
    END FUNCTION boundfunct
    !------------------------------------------------------------------------------------
    !>
    !! @brief Calculates the derivative of the density-dependent dielectric
    !! constant
    !!
    !! ifunct = 0 => original Fattebert and Gygi function
    !!
    !! @param[in]    rho      electrostatic density
    !! @param[in]    rhomax   maximum density cutoff
    !! @param[in]    rhomin   minimum density cutoff
    !! @param[in]    tbeta
    !! @param[in]    const
    !! @param[in]    ifunct
    !! @return       the second derivative of the boundary function
    !!
    !------------------------------------------------------------------------------------
    REAL(DP) FUNCTION dboundfunct(rho, rhomax, rhomin, tbeta, const, ifunct)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP) :: rho
        REAL(DP) :: rhomax
        REAL(DP) :: rhomin
        REAL(DP) :: tbeta
        REAL(DP) :: const
        !
        INTEGER :: ifunct
        !
        REAL(DP) :: arg
        !
        CHARACTER(LEN=80) :: fun_name = 'dboundfunct'
        !
        !--------------------------------------------------------------------------------
        !
        SELECT CASE (ifunct)
            !
        CASE (0)
            dboundfunct = -dsfunct0(rho, rhomax, tbeta)
            !
        CASE (1)
            dboundfunct = -dsfunct1(rho, rhomax, rhomin, tbeta)
            !
        CASE (2)
            !
            dboundfunct = -EXP(LOG(const) * sfunct1(rho, rhomax, rhomin, tbeta)) / &
                          (const - 1.D0) * LOG(const) * &
                          dsfunct1(rho, rhomax, rhomin, tbeta)
            !
        CASE DEFAULT
            CALL io%error(fun_name, 'Unknown boundary type', 1)
            !
        END SELECT
        !
        !--------------------------------------------------------------------------------
    END FUNCTION dboundfunct
    !------------------------------------------------------------------------------------
    !>
    !! @brief Calculates the second derivative of the density-dependent
    !! dielectric constant
    !!
    !! ifunct = 0 => original Fattebery and Gygi function
    !!
    !! @param[in]    rho      electrostatic density
    !! @param[in]    rhomax   maximum density cutoff
    !! @param[in]    rhomin   minimum density cutoff
    !! @param[in]    tbeta
    !! @param[in]    const
    !! @param[in]    ifunct
    !! @return       the second derivative of the boundary function
    !!
    !------------------------------------------------------------------------------------
    REAL(DP) FUNCTION d2boundfunct(rho, rhomax, rhomin, tbeta, const, ifunct)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP) :: rho
        REAL(DP) :: rhomax
        REAL(DP) :: rhomin
        REAL(DP) :: tbeta
        REAL(DP) :: const
        !
        INTEGER :: ifunct
        !
        REAL(DP) :: arg, arg2
        !
        CHARACTER(LEN=80) :: fun_name = 'd2boundfunct'
        !
        !--------------------------------------------------------------------------------
        !
        SELECT CASE (ifunct)
            !
        CASE (0)
            CALL io%error(fun_name, 'Option not yet implemented', 1)
            !
        CASE (1)
            d2boundfunct = -d2sfunct1(rho, rhomax, rhomin, tbeta)
            !
        CASE (2)
            !
            d2boundfunct = -EXP(LOG(const) * sfunct1(rho, rhomax, rhomin, tbeta)) / &
                           (const - 1.D0) * LOG(const) * &
                           (LOG(const) * dsfunct1(rho, rhomax, rhomin, tbeta)**2 + &
                            d2sfunct1(rho, rhomax, rhomin, tbeta))
            !
        CASE DEFAULT
            CALL io%error(fun_name, 'Unknown boundary type', 1)
            !
        END SELECT
        !
        !--------------------------------------------------------------------------------
    END FUNCTION d2boundfunct
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE boundary_tools
!----------------------------------------------------------------------------------------
