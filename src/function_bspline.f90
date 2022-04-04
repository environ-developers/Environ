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
        INTEGER, ALLOCATABLE :: powers(:,:,:)
        REAL(DP), ALLOCATABLE :: coeff(:,:,:)
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
        TYPE( knot_span ), ALLOCATABLE :: spans(:,:)
        REAL(DP), ALLOCATABLE :: u(:,:)
        !
        INTEGER :: span_num, degree, knot_num
        REAL(DP) :: m_spread, norm
        !
        !--------------------------------------------------------------------------------
    CONTAINS
        !--------------------------------------------------------------------------------
        !
        PROCEDURE :: density => density_of_function
        PROCEDURE :: gradient => gradient_of_function
        PROCEDURE :: setup => setup_of_function
        !
        PROCEDURE :: get_u, calc_val, pts_in_span
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
    SUBROUTINE density_of_function(this, density, zero, ir_vals, vals)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_function_bspline), INTENT(INOUT) :: this
        INTEGER, OPTIONAL, INTENT(OUT) :: ir_vals(:)
        REAL(DP), OPTIONAL, INTENT(OUT) :: vals(:)
        LOGICAL, OPTIONAL, INTENT(IN) :: zero
        !
        TYPE(environ_density), INTENT(INOUT) :: density
        !
        INTEGER :: i, uidx(3)
        LOGICAL :: physical
        REAL(DP) :: r(3), r2, length
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
        CALL this%pts_in_span(density%cell%at, density%cell%nr)
        CALL this%setup(this%pos)
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
            SELECT CASE (dim)
                !
            CASE (0)
                this%norm = charge
                !
            CASE (1)
                length = ABS(cell%at(axis, axis))
                this%norm = charge / length
                !
            CASE (2)
                length = ABS(cell%at(axis, axis))
                this%norm = charge * length / cell%omega
                !
            CASE DEFAULT
                CALL io%error(sub_name, "Unexpected system dimensions", 1)
                !
            END SELECT
            !
            !----------------------------------------------------------------------------
            !
            DO i = 1, cell%ir_end
                !
                CALL cell%get_min_distance(i, dim, axis, pos, r, r2, physical)
                ! compute minimum distance using minimum image convention
                !
                IF (.NOT. physical) CYCLE
                !
                uidx = this%get_u(r)
                !
                IF (ANY(uidx >= this%span_num)) CYCLE
                !
                ! Calculate the bspline value at a given point
                !
                density%of_r(i) = density%of_r(i) + this%calc_val(r, uidx)
                !
            END DO
            !
            this%norm = this%norm / density%integrate()
            density%of_r = density%of_r * this%norm
            !
        END ASSOCIATE
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE density_of_function
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE gradient_of_function(this, gradient, zero, ir_vals, vals, grid_pts)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_function_bspline), INTENT(IN) :: this
        INTEGER, OPTIONAL, INTENT(IN) :: ir_vals(:), grid_pts
        REAL(DP), OPTIONAL, INTENT(OUT) :: vals(:, :)
        LOGICAL, OPTIONAL, INTENT(IN) :: zero
        !
        TYPE(environ_gradient), INTENT(INOUT) :: gradient
        !
        INTEGER :: i, a, uidx(3)
        LOGICAL :: physical
        REAL(DP) :: r(3), r2, val, uval(3)
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
            !----------------------------------------------------------------------------
            ! Set local parameters
            !
            !----------------------------------------------------------------------------
            !
            DO i = 1, cell%ir_end
                !
                CALL cell%get_min_distance(i, dim, axis, pos, r, r2, physical)
                ! compute minimum distance using minimum image convention
                !
                IF (.NOT. physical) CYCLE
                !
                DO a = 1, 3
                    uval(a) = r(a) + cell%at(a,a)*0.5D0
                END DO
                !
                IF (ANY(uval > MAXVAL(this%u)) .OR. ANY(uval < MINVAL(this%u))) THEN
                    val = 0.D0
                ELSE
                    !
                    uidx = this%get_u(uval)
                    !
                    !
                    ! Calculate gradient of bspline function at a given point
                    ! #TODO Need to work on calculating gradient correctly for 3D bspline
                    val = this%calc_val(uval, uidx, this%degree - 1)
                    val = val - this%calc_val(uval, uidx, this%degree - 1, .TRUE.)
                    !
                END IF
                !
                print *, i, uval(1), val
                !
            END DO
            !
        END ASSOCIATE
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE gradient_of_function
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
                IF (u_in(a) >= this%u(a,i) .AND. u_in(a) < this%u(a,i+1)) THEN
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
    SUBROUTINE setup_of_function(this, pos)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_function_bspline), INTENT(INOUT) :: this
        REAL(DP), INTENT(IN) :: pos(3)
        !
        CHARACTER(LEN=80) :: sub_name = 'setup_of_function'
        !
        INTEGER :: i, j, k, a, idx
        INTEGER, ALLOCATABLE :: pows(:)
        REAL(DP) :: cvals(4), fluff, dx
        !
        !--------------------------------------------------------------------------------
        !
        print *, 'Number of knots: ', this%knot_num
        ALLOCATE(this%u(3,this%knot_num))
        dx = this%m_spread * this%spread / REAL(this%knot_num - 1, DP)
        !
        DO i = 1, 3
            !
            DO j = 1, this%knot_num
                !
                this%u(i,j) = pos(i) - this%m_spread * this%spread / 2.D0 + (j - 1) * dx
                !
            END DO
            !
            print *, 'Knot values span: ', this%u(i,1), this%u(i,this%knot_num)
            !
        END DO
        !
        this%span_num = this%knot_num - 1
        this%degree = this%span_num - 1
        !
        ALLOCATE(pows(0:this%degree))
        ALLOCATE(this%spans(3,this%span_num))
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
                        ALLOCATE(this%spans(a,j)%coeff(this%span_num,0:this%degree,0:this%degree))
                        ALLOCATE(this%spans(a,j)%powers(this%span_num,0:this%degree,0:this%degree))
                        !
                        this%spans(a,j)%coeff = 0.D0
                        this%spans(a,j)%powers = -1
                        !
                        this%spans(a,j)%coeff(j,0,0) = 1.D0
                        this%spans(a,j)%powers(j,0,0) = 0
                        !
                    ELSE IF (j + i <= this%span_num) THEN
                        !
                        cvals(1) = 1.D0 / (this%u(a,j+i) - this%u(a,j))
                        cvals(2) = -this%u(a,j) * cvals(1)
                        cvals(4) = -1.D0 / (this%u(a,j+i+1) - this%u(a,j+1))
                        cvals(3) = -this%u(a,j+i+1) * cvals(4)
                        !
                        DO k = 1, i
                            !
                            ! Updating variable powers
                            this%spans(a,j)%powers(k+j-1,i,:) = pows
                            this%spans(a,j)%powers(k+j,i,:) = pows
                            !
                            ! First term in B-spline equation
                            this%spans(a,j)%coeff(k+j-1,i,:) = this%spans(a,j)%coeff(k+j-1,i,:) + &
                                                        this%spans(a,j)%coeff(k+j-1,i-1,:) * cvals(2)
                            this%spans(a,j)%coeff(k+j-1,i,1:i) = this%spans(a,j)%coeff(k+j-1,i,1:i) + &
                                                    this%spans(a,j)%coeff(k+j-1,i-1,0:i-1) * cvals(1)
                            !
                            ! Second term in B-spline equation
                            this%spans(a,j)%coeff(k+j,i,:) = this%spans(a,j)%coeff(k+j,i,:) + &
                                                        this%spans(a,j+1)%coeff(k+j,i-1,:) * cvals(3)
                            this%spans(a,j)%coeff(k+j,i,1:i) = this%spans(a,j)%coeff(k+j,i,1:i) + &
                                                    this%spans(a,j+1)%coeff(k+j,i-1,0:i-1) * cvals(4)
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
        OPEN(unit=300,file='setup-info',status='unknown')
        DO a=1,3
           WRITE(300,"(A,I4)") 'Dimension: ', a
           DO i=0,this%degree
               WRITE(300,"(5X,A,I4)") 'Degree: ', i
               DO j=1,this%span_num
                   WRITE(300,"(10X,A,I4)") 'Span: ', j
                   DO k=1,this%span_num
                       WRITE(300,"(15X,A,10I4)") 'Linear Powers: ', this%spans(a,j)%powers(k,i,:)
                       WRITE(300,"(15X,A,10F17.8)") 'Coefficients: ', this%spans(a,j)%coeff(k,i,:)
                   END DO
                   WRITE(300,"(/)")
               END DO
               WRITE(300,"(/)")
           END DO
           WRITE(300,"(/)")
        END DO
        !flush(6)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE setup_of_function
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    REAL(DP) FUNCTION calc_val(this, u_in, idx, deg, next)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_function_bspline), INTENT(IN) :: this
        REAL(DP), INTENT(IN) :: u_in(3)
        INTEGER, INTENT(IN) :: idx(3)
        INTEGER, OPTIONAL, INTENT(IN) :: deg
        LOGICAL, OPTIONAL, INTENT(IN) :: next
        !
        CHARACTER(LEN=80) :: sub_name = 'calc_val'
        !
        INTEGER :: i, degree, span
        REAL(DP) :: p, cons
        !
        !--------------------------------------------------------------------------------
        !
        degree = this%degree
        IF (PRESENT(deg)) degree = deg
        !
        span = 1
        IF (PRESENT(next)) THEN
            IF (next) span = 2
        END IF
        !
        calc_val = 1.D0
        !
        DO i = 1, 3
            !
            ASSOCIATE (pows => this%spans(i,span)%powers, &
                       coeffs => this%spans(i,span)%coeff)
                !
                calc_val = calc_val * SUM( coeffs(idx(i),degree,:)*u_in(i)**REAL(pows(idx(i),degree,:),DP))
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
    SUBROUTINE pts_in_span(this, at, nr)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_function_bspline), INTENT(INOUT) :: this
        REAL(DP), INTENT(IN) :: at(3,3)
        INTEGER, INTENT(IN) :: nr(3)
        !
        INTEGER :: i
        REAL(DP) :: length, frac
        !
        CHARACTER(LEN=80) :: sub_name = 'pts_in_span'
        !
        !--------------------------------------------------------------------------------
        !
        this%knot_num = 6
        this%m_spread = 5.D0
        !
        DO i = 1, 3
            !
            length = SQRT(SUM(at(i,:)*at(i,:)))
            frac = this%m_spread * this%spread / length
            this%knot_num = CEILING( nr(i) * frac ) - 1
            !
        END DO
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE pts_in_span
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE class_function_bspline
!----------------------------------------------------------------------------------------
