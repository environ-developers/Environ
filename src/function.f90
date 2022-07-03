!----------------------------------------------------------------------------------------
!
! Copyright (C) 2018-2021 ENVIRON (www.quantum-environ.org)
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
!          Edan Bainglass     (Department of Physics, UNT)
!
!----------------------------------------------------------------------------------------
!>
!! Module containing the main routines to handle environ_functions
!! derived data types.
!!
!! Environ_functions contains all the details of analytic functions needed
!! by Environ modules and defined on the three-dimensional real-space
!! domain, together with the routines to handle the derived data type and
!! to generate the functions from their parameters.
!!
!! Common method parameters:
!!
!! - zero      (LOGICAL) -> zero out register
!! - ir       *(INTEGER) -> indices of points of interest
!! - vals        *(REAL) -> values of points of interest
!! - r_vals      *(REAL) -> displacements of points of interest
!! - dist_vals   *(REAL) -> distances of points of interest
!!
!! * OPTIONAL
!!
!----------------------------------------------------------------------------------------
MODULE class_function
    !------------------------------------------------------------------------------------
    !
    USE class_io, ONLY: io
    !
    USE environ_param, ONLY: DP
    !
    USE class_cell
    USE class_density
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
    TYPE, ABSTRACT, PUBLIC :: environ_function
        !--------------------------------------------------------------------------------
        !
        INTEGER :: f_type, axis, dim
        REAL(DP) :: width, spread, volume
        !
        REAL(DP), POINTER :: pos(:) => NULL()
        ! environ_functions are not designed to be mobile, thus position
        ! can be included in the definition of the type
        !
        !--------------------------------------------------------------------------------
    CONTAINS
        !--------------------------------------------------------------------------------
        !
        PROCEDURE, PRIVATE :: create => create_environ_function
        PROCEDURE :: init => init_environ_function
        PROCEDURE :: destroy => destroy_environ_function
        !
        PROCEDURE :: density => density_of_function
        PROCEDURE :: gradient => gradient_of_function
        PROCEDURE :: laplacian => laplacian_of_function
        PROCEDURE :: hessian => hessian_of_function
        PROCEDURE :: derivative => derivative_of_function
        !
        PROCEDURE :: quad_corr => quadrapole_corrections
        !
        !--------------------------------------------------------------------------------
    END TYPE environ_function
    !------------------------------------------------------------------------------------
    !
    REAL(DP), PUBLIC, PARAMETER :: func_tol = 1.D-10
    REAL(DP), PUBLIC, PARAMETER :: exp_tol = 4.D1
    !
    !------------------------------------------------------------------------------------
CONTAINS
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !                                   ADMIN METHODS
    !
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE create_environ_function(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_function), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: routine = 'create_environ_function'
        !
        !--------------------------------------------------------------------------------
        !
        IF (ASSOCIATED(this%pos)) CALL io%create_error(routine)
        !
        !--------------------------------------------------------------------------------
        !
        this%f_type = 0
        this%axis = 0
        this%dim = 0
        this%width = 0.D0
        this%spread = 0.D0
        this%volume = 0.D0
        !
        NULLIFY (this%pos)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE create_environ_function
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_environ_function(this, f_type, f_axis, f_dim, f_width, f_spread, &
                                     f_volume, f_pos)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: f_type, f_dim, f_axis
        REAL(DP), INTENT(IN) :: f_width, f_spread, f_volume
        REAL(DP), OPTIONAL, TARGET, INTENT(IN) :: f_pos(:)
        !
        CLASS(environ_function), INTENT(INOUT) :: this
        !
        !--------------------------------------------------------------------------------
        !
        CALL this%create()
        !
        this%f_type = f_type
        this%dim = f_dim
        this%axis = f_axis
        this%spread = f_spread
        this%width = f_width
        this%volume = f_volume
        !
        IF (PRESENT(f_pos)) THEN
            this%pos => f_pos
        ELSE
            ALLOCATE (this%pos(3))
            this%pos = 0.D0
        END IF
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_environ_function
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE destroy_environ_function(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_function), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: routine = 'destroy_environ_function'
        !
        !--------------------------------------------------------------------------------
        !
        IF (.NOT. ASSOCIATED(this%pos)) CALL io%destroy_error(routine)
        !
        !--------------------------------------------------------------------------------
        !
        NULLIFY (this%pos)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE destroy_environ_function
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !                                  FUNCTION METHODS
    !
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE density_of_function(this, density, zero, ir, vals, r_vals, dist_vals)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        LOGICAL, OPTIONAL, INTENT(IN) :: zero
        !
        CLASS(environ_function), INTENT(INOUT) :: this
        TYPE(environ_density), INTENT(INOUT) :: density
        !
        INTEGER, OPTIONAL, INTENT(OUT) :: ir(:)
        REAL(DP), OPTIONAL, INTENT(OUT) :: vals(:), r_vals(:, :), dist_vals(:)
        !
        CHARACTER(LEN=80) :: routine = 'density_of_function'
        !
        !--------------------------------------------------------------------------------
        !
        CALL io%error(routine, "Not implemented", 1)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE density_of_function
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE gradient_of_function(this, gradient, zero, ir, vals, r_vals, dist_vals)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_function), INTENT(IN) :: this
        !
        LOGICAL, OPTIONAL, INTENT(IN) :: zero
        INTEGER, OPTIONAL, INTENT(IN) :: ir(:)
        REAL(DP), OPTIONAL, INTENT(IN) :: r_vals(:, :), dist_vals(:)
        !
        TYPE(environ_gradient), INTENT(INOUT) :: gradient
        !
        REAL(DP), OPTIONAL, INTENT(OUT) :: vals(:, :)
        !
        CHARACTER(LEN=80) :: routine = 'gradient_of_function'
        !
        !--------------------------------------------------------------------------------
        !
        CALL io%error(routine, "Not implemented", 1)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE gradient_of_function
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE laplacian_of_function(this, laplacian, zero, ir, r_vals, dist_vals)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_function), INTENT(IN) :: this
        !
        LOGICAL, OPTIONAL, INTENT(IN) :: zero
        INTEGER, OPTIONAL, INTENT(IN) :: ir(:)
        REAL(DP), OPTIONAL, INTENT(IN) :: r_vals(:, :), dist_vals(:)
        !
        TYPE(environ_density), INTENT(INOUT) :: laplacian
        !
        CHARACTER(LEN=80) :: routine = 'laplacian_of_function'
        !
        !--------------------------------------------------------------------------------
        !
        CALL io%error(routine, "Not implemented", 1)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE laplacian_of_function
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE hessian_of_function(this, hessian, zero, ir, r_vals, dist_vals)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_function), INTENT(IN) :: this
        !
        LOGICAL, OPTIONAL, INTENT(IN) :: zero
        INTEGER, OPTIONAL, INTENT(IN) :: ir(:)
        REAL(DP), OPTIONAL, INTENT(IN) :: r_vals(:, :), dist_vals(:)
        !
        TYPE(environ_hessian), INTENT(INOUT) :: hessian
        !
        CHARACTER(LEN=80) :: routine = 'hessian_of_function'
        !
        !--------------------------------------------------------------------------------
        !
        CALL io%error(routine, "Not implemented", 1)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE hessian_of_function
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE derivative_of_function(this, derivative, zero, ir, r_vals, dist_vals)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_function), INTENT(IN) :: this
        !
        LOGICAL, OPTIONAL, INTENT(IN) :: zero
        INTEGER, OPTIONAL, INTENT(IN) :: ir(:)
        REAL(DP), OPTIONAL, INTENT(IN) :: r_vals(:, :), dist_vals(:)
        !
        TYPE(environ_density), INTENT(INOUT) :: derivative
        !
        CHARACTER(LEN=80) :: routine = 'derivative_of_function'
        !
        !--------------------------------------------------------------------------------
        !
        CALL io%error(routine, "Not implemented", 1)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE derivative_of_function
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    REAL(DP) FUNCTION quadrapole_corrections(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_function), INTENT(IN) :: this
        !
        CHARACTER(LEN=80) :: routine = 'quadrapole_corrections'
        !
        !--------------------------------------------------------------------------------
        !
        quadrapole_corrections = 0.D0
        !
        CALL io%error(routine, "Not implemented", 1)
        !
        !--------------------------------------------------------------------------------
    END FUNCTION quadrapole_corrections
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE class_function
!----------------------------------------------------------------------------------------
