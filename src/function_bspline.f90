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
        TYPE( knot_span ), ALLOCATABLE :: spans(:)
        !
        INTEGER :: span_num
        !
        !--------------------------------------------------------------------------------
    CONTAINS
        !--------------------------------------------------------------------------------
        !
        PROCEDURE :: density => density_of_function
        PROCEDURE :: gradient => gradient_of_function
        PROCEDURE :: setup => setup_of_function
        !
        PROCEDURE, PRIVATE :: get_u
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
        INTEGER :: i, j, k, pow, uidx
        LOGICAL :: physical
        REAL(DP) :: r(3), r2, dist, uval, val, coef
        !
        CHARACTER(LEN=80) :: sub_name = 'density_of_function'
        !
        !--------------------------------------------------------------------------------
        !
        IF (this%degree < 0) &
            CALL io%error(sub_name, 'Wrong value for target degree', 1)
        !
        IF (this%degree > SIZE(this%u) - 1) &
            CALL io%error(sub_name, 'Degree cant be larger than number of knot spans', 1)
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
                   u => this%u, &
                   axis => this%axis)
            !
            CALL this%setup()
            !
            DO i = 1, cell%ir_end
                !
                CALL cell%get_min_distance(i, dim, axis, pos, r, r2, physical)
                ! compute minimum distance using minimum image convention
                !
                IF (.NOT. physical) CYCLE
                !
                dist = SQRT(r2)
                IF (r(1)+cell%at(1,1)*0.5D0 == cell%at(1,1) .AND. i /= 1 ) EXIT
                !
                uval = r(1) + cell%at(1,1)*0.5D0
                uidx = this%get_u(uval)
                !
                ! Calculate the bspline value at a given point
                !
                val = 0.D0
                !
                print *, i, uval, val
                !
            END DO
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
        CHARACTER(LEN=80) :: sub_name = 'gradient_of_function'
        !
        !--------------------------------------------------------------------------------
        !
        CALL io%error(sub_name, "Not implemented", 1)
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
    INTEGER FUNCTION get_u(this, u_in)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_function_bspline), INTENT(IN) :: this
        REAL(DP), INTENT(IN) :: u_in
        !
        CHARACTER(LEN=80) :: sub_name = 'get_u'
        !
        INTEGER :: i
        !
        !--------------------------------------------------------------------------------
        !
        get_u = this%span_num
        !
        DO i = 1, this%span_num
            !
            IF (u_in >= this%u(i) .AND. u_in < this%u(i+1)) THEN
                !
                get_u = i
                !
            END IF
            !
        END DO
        !
        !--------------------------------------------------------------------------------
    END FUNCTION get_u
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE setup_of_function(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_function_bspline), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: sub_name = 'setup_of_function'
        !
        INTEGER :: i, j, k, pows(0:this%degree), idx
        REAL(DP) :: cvals(4)
        !
        !--------------------------------------------------------------------------------
        !
        this%span_num = SIZE(this%u) - 1
        ALLOCATE(this%spans(this%span_num))
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
                    ALLOCATE(this%spans(j)%coeff(this%span_num,0:this%degree,0:this%degree))
                    ALLOCATE(this%spans(j)%powers(this%span_num,0:this%degree,0:this%degree))
                    !
                    this%spans(j)%coeff = 0.D0
                    this%spans(j)%powers = -1
                    !
                    this%spans(j)%coeff(j,0,0) = 1.D0
                    this%spans(j)%powers(j,0,0) = 0
                    !
                ELSE IF (j + i <= this%span_num) THEN
                    !
                    cvals(1) = 1.D0 / (this%u(j+i) - this%u(j))
                    cvals(2) = -this%u(j) * cvals(1)
                    cvals(4) = -1.D0 / (this%u(j+i+1) - this%u(j+1))
                    cvals(3) = -this%u(j+i+1) * cvals(4)
                    !
                    DO k = 1, i
                        !
                        ! Updating variable powers
                        this%spans(j)%powers(k+j-1,i,:) = pows
                        this%spans(j)%powers(k+j,i,:) = pows
                        !
                        ! First term in B-spline equation
                        this%spans(j)%coeff(k+j-1,i,:) = this%spans(j)%coeff(k+j-1,i,:) + &
                                                    this%spans(j)%coeff(k+j-1,i-1,:) * cvals(2)
                        this%spans(j)%coeff(k+j-1,i,1:i) = this%spans(j)%coeff(k+j-1,i,1:i) + &
                                                this%spans(j)%coeff(k+j-1,i-1,0:i-1) * cvals(1)
                        !
                        ! Second term in B-spline equation
                        this%spans(j)%coeff(k+j,i,:) = this%spans(j)%coeff(k+j,i,:) + &
                                                    this%spans(j+1)%coeff(k+j,i-1,:) * cvals(3)
                        this%spans(j)%coeff(k+j,i,1:i) = this%spans(j)%coeff(k+j,i,1:i) + &
                                                this%spans(j+1)%coeff(k+j,i-1,0:i-1) * cvals(4)
                        !
                    END DO
                    !
                END IF
                !
            END DO
            !
        END DO
        !
        !DO i=0,this%degree
        !    WRITE(*,"(A,I4)") 'Degree: ', i
        !    DO j=1,this%span_num
        !        WRITE(*,"(5X,A,I4)") 'Span: ', j
        !        DO k=1,this%span_num
        !            WRITE(*,"(10X,A,10I4)") 'Linear Powers: ', this%spans(j)%powers(k,i,:)
        !            WRITE(*,"(10X,A,10F17.8)") 'Coefficients: ', this%spans(j)%coeff(k,i,:)
        !        END DO
        !        WRITE(*,"(/)")
        !    END DO
        !    WRITE(*,"(/)")
        !END DO
        flush(6)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE setup_of_function
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE class_function_bspline
!----------------------------------------------------------------------------------------
