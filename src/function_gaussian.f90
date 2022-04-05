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
!          Edan Bainglass     (Department of Physics, UNT)
!
!----------------------------------------------------------------------------------------
!>
!!
!----------------------------------------------------------------------------------------
MODULE class_function_gaussian
    !------------------------------------------------------------------------------------
    !
    USE class_io, ONLY: io
    !
    USE environ_param, ONLY: DP, sqrtpi
    !
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
    TYPE, EXTENDS(environ_function), PUBLIC :: environ_function_gaussian
        !--------------------------------------------------------------------------------
        !
        !--------------------------------------------------------------------------------
    CONTAINS
        !--------------------------------------------------------------------------------
        !
        PROCEDURE :: density => density_of_function
        PROCEDURE :: gradient => gradient_of_function
        !
        !--------------------------------------------------------------------------------
    END TYPE environ_function_gaussian
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
        CLASS(environ_function_gaussian), INTENT(INOUT) :: this
        INTEGER, OPTIONAL, INTENT(OUT) :: ir_vals(:)
        REAL(DP), OPTIONAL, INTENT(OUT) :: vals(:), r_vals(:, :), dist_vals(:)
        LOGICAL, OPTIONAL, INTENT(IN) :: zero
        !
        TYPE(environ_density), INTENT(INOUT) :: density
        !
        INTEGER :: i
        LOGICAL :: physical
        REAL(DP) :: r(3), r2, scale, length
        REAL(DP), ALLOCATABLE :: local(:)
        !
        CHARACTER(LEN=80) :: sub_name = 'density_of_function'
        !
        !--------------------------------------------------------------------------------
        !
        IF (ABS(this%volume) < func_tol) RETURN
        !
        IF (ABS(this%spread) < func_tol) &
            CALL io%error(sub_name, "Wrong spread for Gaussian function", 1)
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
                   spread => this%spread, &
                   charge => this%volume, &
                   dim => this%dim, &
                   axis => this%axis)
            !
            !----------------------------------------------------------------------------
            ! Set local parameters
            !
            SELECT CASE (dim)
                !
            CASE (0)
                scale = charge / (sqrtpi * spread)**3
                !
            CASE (1)
                length = ABS(cell%at(axis, axis))
                scale = charge / length / (sqrtpi * spread)**2
                !
            CASE (2)
                length = ABS(cell%at(axis, axis))
                scale = charge * length / cell%omega / (sqrtpi * spread)
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
                CALL cell%get_min_distance(i, dim, axis, pos, r, r2, physical)
                ! compute minimum distance using minimum image convention
                !
                IF (.NOT. physical) CYCLE
                !
                r2 = r2 / spread**2
                !
                IF (r2 <= exp_tol) local(i) = EXP(-r2) ! compute Gaussian function
                !
            END DO
            !
            density%of_r = density%of_r + scale * local
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
    SUBROUTINE gradient_of_function(this, gradient, zero, ir_vals, vals, grid_pts)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_function_gaussian), INTENT(IN) :: this
        INTEGER, OPTIONAL, INTENT(IN) :: ir_vals( : ), grid_pts
        REAL(DP), OPTIONAL, INTENT(OUT) :: vals(:,:)
        LOGICAL, OPTIONAL, INTENT(IN) :: zero
        !
        TYPE(environ_gradient), INTENT(INOUT) :: gradient
        !
        INTEGER :: i
        LOGICAL :: physical
        REAL(DP) :: r(3), r2, scale, length
        REAL(DP), ALLOCATABLE :: gradlocal(:, :)
        !
        CHARACTER(LEN=80) :: sub_name = 'gradient_of_function'
        !
        !--------------------------------------------------------------------------------
        !
        IF (ABS(this%volume) < func_tol) RETURN
        !
        IF (ABS(this%spread) < func_tol) &
            CALL io%error(sub_name, "Wrong spread for Gaussian function", 1)
        !
        IF (this%axis < 1 .OR. this%axis > 3) &
            CALL io%error(sub_name, "Wrong value of axis", 1)
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
                   spread => this%spread, &
                   charge => this%volume, &
                   dim => this%dim, &
                   axis => this%axis)
            !
            !----------------------------------------------------------------------------
            ! Set local parameters
            !
            SELECT CASE (dim)
                !
            CASE (0)
                scale = charge / (sqrtpi * spread)**3
                !
            CASE (1)
                length = ABS(cell%at(axis, axis))
                scale = charge / length / (sqrtpi * spread)**2
                !
            CASE (2)
                length = ABS(cell%at(axis, axis))
                scale = charge * length / cell%omega / (sqrtpi * spread)
                !
            CASE DEFAULT
                CALL io%error(sub_name, "Unexpected system dimensions", 1)
                !
            END SELECT
            !
            scale = scale * 2.D0 / spread**2
            !
            !----------------------------------------------------------------------------
            !
            ALLOCATE (gradlocal(3, cell%nnr))
            gradlocal = 0.D0
            !
            DO i = 1, cell%ir_end
                !
                CALL cell%get_min_distance(i, dim, axis, pos, r, r2, physical)
                ! compute minimum distance using minimum image convention
                !
                IF (.NOT. physical) CYCLE
                !
                r2 = r2 / spread**2
                !
                IF (r2 <= exp_tol) gradlocal(:, i) = -EXP(-r2) * r
                ! compute gradient of Gaussian function
                !
            END DO
            !
            gradient%of_r = gradient%of_r + scale * gradlocal
            DEALLOCATE (gradlocal)
            !
        END ASSOCIATE
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE gradient_of_function
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE class_function_gaussian
!----------------------------------------------------------------------------------------
