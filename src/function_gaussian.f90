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
        PROCEDURE :: density => density_of_function_gaussian
        PROCEDURE :: gradient => gradient_of_function_gaussian
        PROCEDURE :: laplacian => laplacian_of_function_gaussian
        PROCEDURE :: hessian => hessian_of_function_gaussian
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
    SUBROUTINE density_of_function_gaussian(this, density, zero)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_function_gaussian), TARGET, INTENT(IN) :: this
        LOGICAL, INTENT(IN), OPTIONAL :: zero
        !
        TYPE(environ_density), TARGET, INTENT(INOUT) :: density
        !
        LOGICAL :: physical
        INTEGER :: ir
        REAL(DP) :: r(3), r2, scale, spr2, length
        REAL(DP), ALLOCATABLE :: local(:)
        !
        CHARACTER(LEN=80) :: sub_name = 'density_of_function_gaussian'
        !
        !--------------------------------------------------------------------------------
        !
        IF (PRESENT(zero)) THEN
            IF (zero) density%of_r = 0.D0
        END IF
        !
        ASSOCIATE (cell => density%cell, &
                   pos => this%pos, &
                   spread => this%spread, &
                   charge => this%volume, &
                   dim => this%dim, &
                   axis => this%axis)
            !
            IF (ABS(charge) < func_tol) RETURN
            !
            IF (ABS(spread) < func_tol) &
                CALL io%error(sub_name, 'Wrong spread for Gaussian function', 1)
            !
            IF (axis < 1 .OR. axis > 3) CALL io%error(sub_name, 'Wrong value of axis', 1)
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
                CALL io%error(sub_name, 'Wrong value of dim', 1)
                !
            END SELECT
            !
            spr2 = spread**2
            !
            ALLOCATE (local(cell%nnr))
            local = 0.D0
            !
            DO ir = 1, cell%ir_end
                !
                CALL cell%get_min_distance(ir, dim, axis, pos, r, r2, physical)
                ! compute minimum distance using minimum image convention
                !
                IF (.NOT. physical) CYCLE
                !
                r2 = r2 / spr2
                !
                IF (r2 <= exp_tol) local(ir) = EXP(-r2) ! compute Gaussian function
                !
            END DO
            !
            density%of_r = density%of_r + scale * local
            DEALLOCATE (local)
            !
        END ASSOCIATE
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE density_of_function_gaussian
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE gradient_of_function_gaussian(this, gradient, zero)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_function_gaussian), TARGET, INTENT(IN) :: this
        LOGICAL, INTENT(IN), OPTIONAL :: zero
        !
        TYPE(environ_gradient), TARGET, INTENT(INOUT) :: gradient
        !
        LOGICAL :: physical
        INTEGER :: ir
        REAL(DP) :: r(3), r2, scale, spr2, length
        REAL(DP), ALLOCATABLE :: gradlocal(:, :)
        !
        CHARACTER(LEN=80) :: sub_name = 'gradient_of_function_gaussian'
        !
        !--------------------------------------------------------------------------------
        !
        IF (PRESENT(zero)) THEN
            IF (zero) gradient%of_r = 0.D0
        END IF
        !
        ASSOCIATE (cell => gradient%cell, &
                   pos => this%pos, &
                   spread => this%spread, &
                   charge => this%volume, &
                   dim => this%dim, &
                   axis => this%axis)
            !
            IF (ABS(charge) < func_tol) RETURN
            !
            IF (ABS(spread) < func_tol) &
                CALL io%error(sub_name, 'Wrong spread for Gaussian function', 1)
            !
            IF (axis < 1 .OR. axis > 3) CALL io%error(sub_name, 'Wrong value of axis', 1)
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
                CALL io%error(sub_name, 'Wrong value of dim', 1)
                !
            END SELECT
            !
            scale = scale * 2.D0 / spread**2
            !
            spr2 = spread**2
            !
            ALLOCATE (gradlocal(3, cell%nnr))
            gradlocal = 0.D0
            !
            DO ir = 1, cell%ir_end
                !
                CALL cell%get_min_distance(ir, dim, axis, pos, r, r2, physical)
                ! compute minimum distance using minimum image convention
                !
                IF (.NOT. physical) CYCLE
                !
                r2 = r2 / spr2
                !
                IF (r2 <= exp_tol) gradlocal(:, ir) = -EXP(-r2) * r
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
    END SUBROUTINE gradient_of_function_gaussian
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE laplacian_of_function_gaussian(this, laplacian, zero)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_function_gaussian), TARGET, INTENT(IN) :: this
        LOGICAL, INTENT(IN), OPTIONAL :: zero
        !
        TYPE(environ_density), TARGET, INTENT(INOUT) :: laplacian
        !
        CHARACTER(LEN=80) :: sub_name = 'laplacian_of_function_gaussian'
        !
        !--------------------------------------------------------------------------------
        !
        CALL io%error(sub_name, 'Options not yet implemented', 1)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE laplacian_of_function_gaussian
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE hessian_of_function_gaussian(this, hessian, zero)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_function_gaussian), TARGET, INTENT(IN) :: this
        LOGICAL, INTENT(IN), OPTIONAL :: zero
        !
        TYPE(environ_hessian), TARGET, INTENT(INOUT) :: hessian
        !
        CHARACTER(LEN=80) :: sub_name = 'hessian_of_function_gaussian'
        !
        !--------------------------------------------------------------------------------
        !
        CALL io%error(sub_name, 'Options not yet implemented', 1)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE hessian_of_function_gaussian
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE class_function_gaussian
!----------------------------------------------------------------------------------------
