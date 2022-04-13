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
! Authors: Gabriel Medrano    (Department of Physics, UNT)
!
!----------------------------------------------------------------------------------------
!>
!!
!----------------------------------------------------------------------------------------
MODULE class_function_bspline
    !------------------------------------------------------------------------------------
    !
    USE class_io, ONLY: io
    USE env_mp, ONLY: env_mp_sum
    !
    USE environ_param, ONLY: DP, sqrtpi, pi, fpi
    !
    USE class_cell
    USE class_density
    USE class_function
    USE class_gradient
    USE class_hessian
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
    TYPE, PUBLIC :: knot_span
        !--------------------------------------------------------------------------------
        !
        INTEGER, ALLOCATABLE :: powers(:, :, :)
        REAL(DP), ALLOCATABLE :: coeff(:, :, :)
        !
        !--------------------------------------------------------------------------------
    END TYPE knot_span
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    TYPE, EXTENDS(environ_function), PUBLIC :: environ_function_bspline
        !--------------------------------------------------------------------------------
        !
        TYPE(knot_span), ALLOCATABLE :: spans(:, :)
        REAL(DP), ALLOCATABLE :: u(:, :)
        !
        INTEGER :: span_num, degree, knot_num
        REAL(DP) :: m_spread, norm, xyz_norm(3)
        !
        !--------------------------------------------------------------------------------
    CONTAINS
        !--------------------------------------------------------------------------------
        !
        PROCEDURE :: density => density_of_function
        PROCEDURE :: gradient => gradient_of_function
        PROCEDURE :: setup => setup_of_function
        PROCEDURE :: se_quad_corr => self_energy_quadrapole_corrections
        !
        PROCEDURE, PRIVATE :: get_u, calc_val, calc_grad_val, bsplinevolume
        !
        !--------------------------------------------------------------------------------
    END TYPE environ_function_bspline
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
    SUBROUTINE density_of_function(this, density, zero, ir_vals, vals, r_vals, dist_vals)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_function_bspline), INTENT(INOUT) :: this
        INTEGER, OPTIONAL, INTENT(OUT) :: ir_vals(:)
        REAL(DP), OPTIONAL, INTENT(OUT) :: vals(:), r_vals(:, :), dist_vals(:)
        LOGICAL, OPTIONAL, INTENT(IN) :: zero
        !
        TYPE(environ_density), INTENT(INOUT) :: density
        !
        INTEGER :: i, uidx(3)
        LOGICAL :: physical
        REAL(DP) :: r(3), length
        REAL(DP), ALLOCATABLE :: local(:)
        !
        CHARACTER(LEN=80) :: sub_name = 'density_of_function'
        !
        !--------------------------------------------------------------------------------
        ! If called directly and not through a functions object, initialize the register
        !
        IF (PRESENT(zero)) THEN
            IF (zero) density%of_r = 0.D0
        END IF
        !
        !--------------------------------------------------------------------------------
        !
        ASSOCIATE (cell => density%cell, &
                   pos => this%pos, &
                   dim => this%dim, &
                   charge => this%volume, &
                   u => this%u, &
                   axis => this%axis)
            !
            !----------------------------------------------------------------------------
            ! Set local parameters
            !
            CALL this%setup(cell%at, cell%nr)
            !
            SELECT CASE (dim)
                !
            CASE (0)
                this%norm = charge / this%bsplinevolume()
                !
            CASE (1)
                length = ABS(cell%at(axis, axis))
                this%norm = charge / length / this%bsplinevolume()
                !
            CASE (2)
                length = ABS(cell%at(axis, axis))
                this%norm = charge * length / cell%omega / this%bsplinevolume()
                !
            CASE DEFAULT
                CALL io%error(sub_name, "Unexpected system dimensions", 1)
                !
            END SELECT
            !
            !----------------------------------------------------------------------------
            !
            ALLOCATE (local(cell%nnr))
            local = 0.D0
            !
            DO i = 1, cell%ir_end
                !
                CALL cell%ir2r(i, r, physical)
                ! compute r vector
                !
                IF (.NOT. physical) CYCLE
                !
                IF (SQRT(SUM((r - this%pos)**2.D0)) > 1.6D0) CYCLE
                !
                uidx = this%get_u(r)
                !
                ! Calculate the bspline value at a given point
                !
                local(i) = this%calc_val(r, uidx)
                !
            END DO
            !
            density%of_r = density%of_r + local * this%norm
            DEALLOCATE (local)
            !
        END ASSOCIATE
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE density_of_function
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE gradient_of_function(this, gradient, zero, ir_vals, vals, grid_pts, &
                                    r_vals, dist_vals)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_function_bspline), INTENT(IN) :: this
        INTEGER, OPTIONAL, INTENT(IN) :: ir_vals(:), grid_pts
        REAL(DP), OPTIONAL, INTENT(OUT) :: vals(:, :), r_vals(:, :), dist_vals(:)
        LOGICAL, OPTIONAL, INTENT(IN) :: zero
        !
        TYPE(environ_gradient), INTENT(INOUT) :: gradient
        !
        INTEGER :: i, uidx(3)
        LOGICAL :: physical
        REAL(DP) :: r(3)
        !
        CHARACTER(LEN=80) :: sub_name = 'gradient_of_function'
        !
        !--------------------------------------------------------------------------------
        !
        IF (.NOT. ALLOCATED(this%spans)) &
            CALL io%error(sub_name, "Powers and coefficients not calculated", 1)
        !
        !--------------------------------------------------------------------------------
        ! If called directly and not through a functions object, initialize the register
        !
        IF (PRESENT(zero)) THEN
            IF (zero) gradient%of_r = 0.D0
        END IF
        !
        !--------------------------------------------------------------------------------
        !
        ASSOCIATE (cell => gradient%cell, &
                   pos => this%pos, &
                   dim => this%dim, &
                   u => this%u, &
                   axis => this%axis)
            !
            DO i = 1, cell%ir_end
                !
                CALL cell%ir2r(i, r, physical)
                ! compute r vector
                !
                IF (.NOT. physical) CYCLE
                !
                IF (SQRT(SUM((r - this%pos)**2.D0)) > 1.6D0) CYCLE
                !
                uidx = this%get_u(r)
                !
                ! Calculate gradient of bspline function at a given point
                gradient%of_r(:, i) = gradient%of_r(:, i) + this%calc_grad_val(r, uidx)
                !
            END DO
            !
        END ASSOCIATE
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE gradient_of_function
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    REAL(DP) FUNCTION self_energy_quadrapole_corrections(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_function_bspline), INTENT(IN) :: this
        !
        REAL(DP) :: term1, term2, x2_vals(3), c1, c2
        INTEGER :: i, j
        !
        CHARACTER(LEN=80) :: sub_name = 'self_energy_quadrapole_corrections'
        !
        !--------------------------------------------------------------------------------
        !
        x2_vals = 0.D0
        !
        DO i = 1, 3
            !
            c1 = this%u(i, this%knot_num) + this%u(i, 1)
            c2 = (c1 / 2.D0)**2.D0
            !
            DO j = 1, this%span_num
                !
                ASSOCIATE (pows => this%spans(i, 1)%powers(j, this%degree, :), &
                           coeffs => this%spans(i, 1)%coeff(j, this%degree, :))
                    !
                    term1 = SUM(coeffs * this%u(i, j)**(pows + 3) / (pows + 3))
                    term2 = SUM(coeffs * this%u(i, j + 1)**(pows + 3) / (pows + 3))
                    !
                    x2_vals(i) = x2_vals(i) + term2 - term1
                    !
                    term1 = SUM(coeffs * this%u(i, j)**(pows + 2) / (pows + 2))
                    term2 = SUM(coeffs * this%u(i, j + 1)**(pows + 2) / (pows + 2))
                    !
                    x2_vals(i) = x2_vals(i) - c1 * (term2 - term1)
                    !
                END ASSOCIATE
                !
            END DO
            !
            x2_vals(i) = x2_vals(i) + c2 * this%xyz_norm(i)
            !
            x2_vals(i) = x2_vals(i) * this%volume / this%xyz_norm(i)
            !
        END DO
        !
        self_energy_quadrapole_corrections = SUM(x2_vals) / 3.D0
        !
        !--------------------------------------------------------------------------------
    END FUNCTION self_energy_quadrapole_corrections
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !                               PRIVATE HELPER METHODS
    !
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    FUNCTION get_u(this, u_in) RESULT(u_out)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_function_bspline), INTENT(IN) :: this
        REAL(DP), INTENT(IN) :: u_in(3)
        INTEGER :: u_out(3)
        !
        CHARACTER(LEN=80) :: sub_name = 'get_u'
        !
        INTEGER :: i, a
        !
        !--------------------------------------------------------------------------------
        !
        u_out = this%span_num
        !
        DO a = 1, 3
            !
            DO i = 1, this%span_num
                !
                IF (u_in(a) >= this%u(a, i) .AND. u_in(a) < this%u(a, i + 1)) THEN
                    !
                    u_out(a) = i
                    !
                END IF
                !
            END DO
            !
        END DO
        !
        !--------------------------------------------------------------------------------
    END FUNCTION get_u
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE setup_of_function(this, at, nr)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_function_bspline), INTENT(INOUT) :: this
        REAL(DP), INTENT(IN) :: at(3, 3)
        INTEGER, INTENT(IN) :: nr(3)
        !
        CHARACTER(LEN=80) :: sub_name = 'setup_of_function'
        !
        INTEGER :: i, j, k, a, idx
        INTEGER, ALLOCATABLE :: pows(:)
        REAL(DP) :: cvals(4), fluff, dx
        !
        !--------------------------------------------------------------------------------
        !
        this%knot_num = 6
        ALLOCATE (this%u(3, this%knot_num))
        !
        DO i = 1, 3
            !
            dx = at(i, i) / nr(i)
            !
            DO j = 1, this%knot_num
                !
                this%u(i, j) = this%pos(i) - dx * this%knot_num / 2.D0 + (j - 1 + 0.5D0) * dx
                !
            END DO
            !
            IF (this%u(i, this%knot_num) - this%u(i, 1) > 1.6D0) THEN
                PRINT *, this%u(i, this%knot_num) - this%u(i, 1)
                CALL io%error(sub_name, 'Spread of Bspline is larger than 1.6 Bohr.', 1)
            END IF
            !
        END DO
        !
        this%span_num = this%knot_num - 1
        this%degree = this%span_num - 1
        !
        ALLOCATE (pows(0:this%degree))
        ALLOCATE (this%spans(3, this%span_num))
        !
        DO a = 1, 3
            !
            pows = -1
            !
            DO i = 0, this%degree
                !
                pows(i) = i
                !
                DO j = 1, this%span_num
                    !
                    IF (i == 0) THEN
                        !
                        ALLOCATE (this%spans(a, j)%coeff(this%span_num, 0:this%degree, 0:this%degree))
                        ALLOCATE (this%spans(a, j)%powers(this%span_num, 0:this%degree, 0:this%degree))
                        !
                        this%spans(a, j)%coeff = 0.D0
                        this%spans(a, j)%powers = -1
                        !
                        this%spans(a, j)%coeff(j, 0, 0) = 1.D0
                        this%spans(a, j)%powers(j, 0, 0) = 0
                        !
                    ELSE IF (j + i <= this%span_num) THEN
                        !
                        cvals(1) = 1.D0 / (this%u(a, j + i) - this%u(a, j))
                        cvals(2) = -this%u(a, j) * cvals(1)
                        cvals(4) = -1.D0 / (this%u(a, j + i + 1) - this%u(a, j + 1))
                        cvals(3) = -this%u(a, j + i + 1) * cvals(4)
                        !
                        DO k = 1, i
                            !
                            ! Updating variable powers
                            this%spans(a, j)%powers(k + j - 1, i, :) = pows
                            this%spans(a, j)%powers(k + j, i, :) = pows
                            !
                            ! First term in B-spline equation
                            this%spans(a, j)%coeff(k + j - 1, i, :) = this%spans(a, j)%coeff(k + j - 1, i, :) + &
                                                                      this%spans(a, j)%coeff(k + j - 1, i - 1, :) * cvals(2)
                            this%spans(a, j)%coeff(k + j - 1, i, 1:i) = this%spans(a, j)%coeff(k + j - 1, i, 1:i) + &
                                                                        this%spans(a, j)%coeff(k + j - 1, i - 1, 0:i - 1) * cvals(1)
                            !
                            ! Second term in B-spline equation
                            this%spans(a, j)%coeff(k + j, i, :) = this%spans(a, j)%coeff(k + j, i, :) + &
                                                                  this%spans(a, j + 1)%coeff(k + j, i - 1, :) * cvals(3)
                            this%spans(a, j)%coeff(k + j, i, 1:i) = this%spans(a, j)%coeff(k + j, i, 1:i) + &
                                                                    this%spans(a, j + 1)%coeff(k + j, i - 1, 0:i - 1) * cvals(4)
                            !
                        END DO
                        !
                    END IF
                    !
                END DO
                !
            END DO
            !
        END DO
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE setup_of_function
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    REAL(DP) FUNCTION calc_val(this, u_in, idx)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_function_bspline), INTENT(IN) :: this
        REAL(DP), INTENT(IN) :: u_in(3)
        INTEGER, INTENT(IN) :: idx(3)
        !
        CHARACTER(LEN=80) :: sub_name = 'calc_val'
        !
        INTEGER :: i, degree, span
        REAL(DP) :: p, cons
        !
        !--------------------------------------------------------------------------------
        !
        degree = this%degree
        span = 1
        !
        calc_val = 1.D0
        !
        DO i = 1, 3
            !
            IF (u_in(i) > MAXVAL(this%u(i, :)) .OR. u_in(i) < MINVAL(this%u(i, :))) calc_val = 0.D0
            !
            ASSOCIATE (pows => this%spans(i, span)%powers(idx(i), degree, :), &
                       coeffs => this%spans(i, span)%coeff(idx(i), degree, :))
                !
                calc_val = calc_val * SUM(coeffs * u_in(i)**pows)
                !
            END ASSOCIATE
            !
        END DO
        !
        !--------------------------------------------------------------------------------
    END FUNCTION calc_val
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    FUNCTION calc_grad_val(this, u_in, idx)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_function_bspline), INTENT(IN) :: this
        REAL(DP), INTENT(IN) :: u_in(3)
        INTEGER, INTENT(IN) :: idx(3)
        REAL(DP) :: calc_grad_val(3)
        !
        CHARACTER(LEN=80) :: sub_name = 'calc_grad_val'
        !
        INTEGER :: i, j, k, degree, span
        REAL(DP) :: local_vals(3), grad_val
        !
        !--------------------------------------------------------------------------------
        !
        degree = this%degree
        span = 1
        !
        DO i = 1, 3
            !
            IF (u_in(i) > MAXVAL(this%u(i, :)) .OR. u_in(i) < MINVAL(this%u(i, :))) THEN
                !
                calc_grad_val = 0.D0
                RETURN
                !
            END IF
            !
        END DO
        !
        calc_grad_val = 1.D0
        !
        DO i = 1, 3
            !
            ASSOCIATE (pows => this%spans(i, span)%powers(idx(i), degree, :), &
                       coeffs => this%spans(i, span)%coeff(idx(i), degree, :))
                !
                local_vals(i) = SUM(coeffs * u_in(i)**REAL(pows, DP))
                !
            END ASSOCIATE
            !
        END DO
        !
        DO i = 1, 3
            !
            DO j = 1, 3
                !
                IF (i == j) THEN
                    !
                    grad_val = 0.D0
                    DO k = 2, this%degree + 1
                        !
                        ASSOCIATE (pows => this%spans(i, span)%powers(idx(i), degree, :), &
                                   coeffs => this%spans(i, span)%coeff(idx(i), degree, :))
                            !
                            grad_val = grad_val + coeffs(k) * pows(k) * u_in(i)**REAL(pows(k) - 1, DP)
                            !
                        END ASSOCIATE
                        !
                    END DO
                    !
                    calc_grad_val(i) = calc_grad_val(i) * grad_val
                    !
                ELSE
                    !
                    calc_grad_val(i) = calc_grad_val(i) * local_vals(j)
                    !
                END IF
                !
            END DO
            !
        END DO
        !
        calc_grad_val = calc_grad_val * this%norm
        !
        !--------------------------------------------------------------------------------
    END FUNCTION calc_grad_val
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    REAL(DP) FUNCTION bsplinevolume(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_function_bspline), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: sub_name = 'bsplinevolume'
        !
        INTEGER :: i, j
        REAL(DP) :: term1, term2
        !
        !--------------------------------------------------------------------------------
        !
        bsplinevolume = 1.D0
        this%xyz_norm = 0.D0
        !
        DO i = 1, 3
            !
            DO j = 1, this%span_num
                !
                ASSOCIATE (pows => this%spans(i, 1)%powers(j, this%degree, :), &
                           coeffs => this%spans(i, 1)%coeff(j, this%degree, :))
                    !
                    term1 = SUM(coeffs * this%u(i, j)**(pows + 1) / (pows + 1))
                    term2 = SUM(coeffs * this%u(i, j + 1)**(pows + 1) / (pows + 1))
                    !
                    this%xyz_norm(i) = this%xyz_norm(i) + term2 - term1
                    !
                END ASSOCIATE
                !
            END DO
            !
        END DO
        !
        bsplinevolume = this%xyz_norm(1) * this%xyz_norm(2) * this%xyz_norm(3)
        !
        !--------------------------------------------------------------------------------
    END FUNCTION bsplinevolume
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE class_function_bspline
!----------------------------------------------------------------------------------------
