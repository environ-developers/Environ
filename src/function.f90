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
!! Module containing the main routines to handle environ_functions
!! derived data types.
!!
!! Environ_functions contains all the details of analytic functions needed
!! by Environ modules and defined on the three-dimensional real-space
!! domain, together with the routines to handle the derived data type and
!! to generate the functions from their parameters.
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
        PROCEDURE :: copy => copy_environ_function
        PROCEDURE :: destroy => destroy_environ_function
        !
        PROCEDURE(get_density), DEFERRED :: density
        PROCEDURE(get_gradient), DEFERRED :: gradient
        PROCEDURE(get_laplacian), DEFERRED :: laplacian
        PROCEDURE(get_hessian), DEFERRED :: hessian
        PROCEDURE(get_derivative), DEFERRED :: derivative
        !
        !--------------------------------------------------------------------------------
    END TYPE environ_function
    !------------------------------------------------------------------------------------
    !
    ABSTRACT INTERFACE
        SUBROUTINE get_density(this, density, zero)
            IMPORT environ_function, environ_density
            CLASS(environ_function), TARGET, INTENT(IN) :: this
            LOGICAL, INTENT(IN), OPTIONAL :: zero
            TYPE(environ_density), TARGET, INTENT(INOUT) :: density
        END SUBROUTINE
        SUBROUTINE get_gradient(this, gradient, zero)
            IMPORT environ_function, environ_gradient
            CLASS(environ_function), TARGET, INTENT(IN) :: this
            LOGICAL, INTENT(IN), OPTIONAL :: zero
            TYPE(environ_gradient), TARGET, INTENT(INOUT) :: gradient
        END SUBROUTINE
        SUBROUTINE get_laplacian(this, laplacian, zero)
            IMPORT environ_function, environ_density
            CLASS(environ_function), TARGET, INTENT(IN) :: this
            LOGICAL, INTENT(IN), OPTIONAL :: zero
            TYPE(environ_density), TARGET, INTENT(INOUT) :: laplacian
        END SUBROUTINE
        SUBROUTINE get_hessian(this, hessian, zero)
            IMPORT environ_function, environ_hessian
            CLASS(environ_function), TARGET, INTENT(IN) :: this
            LOGICAL, INTENT(IN), OPTIONAL :: zero
            TYPE(environ_hessian), TARGET, INTENT(INOUT) :: hessian
        END SUBROUTINE
        SUBROUTINE get_derivative(this, derivative, zero)
            IMPORT environ_function, environ_density
            CLASS(environ_function), TARGET, INTENT(IN) :: this
            LOGICAL, INTENT(IN), OPTIONAL :: zero
            TYPE(environ_density), TARGET, INTENT(INOUT) :: derivative
        END SUBROUTINE
    END INTERFACE
    !
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
        CHARACTER(LEN=80) :: sub_name = 'create_environ_function'
        !
        !--------------------------------------------------------------------------------
        !
        IF (ASSOCIATED(this%pos)) CALL io%create_error(sub_name)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE create_environ_function
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_environ_function(this, type_in, axis, dim_in, width, spread_in, &
                                     volume_in, pos)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: type_in, dim_in, axis
        REAL(DP), INTENT(IN) :: width, spread_in, volume_in
        REAL(DP), TARGET, INTENT(IN), OPTIONAL :: pos(:)
        !
        CLASS(environ_function), INTENT(INOUT) :: this
        !
        !--------------------------------------------------------------------------------
        !
        CALL this%create()
        !
        this%f_type = type_in
        this%dim = dim_in
        this%axis = axis
        this%spread = spread_in
        this%width = width
        this%volume = volume_in
        !
        IF (PRESENT(pos)) THEN
            this%pos => pos
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
    SUBROUTINE copy_environ_function(this, copy)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_function), INTENT(IN) :: this
        !
        CLASS(environ_function), INTENT(OUT) :: copy
        !
        CHARACTER(LEN=80) :: sub_name = 'copy_environ_function'
        !
        !--------------------------------------------------------------------------------
        !
        CALL copy%init(this%f_type, this%axis, this%dim, this%width, this%spread, &
                       this%volume, this%pos)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE copy_environ_function
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
        CHARACTER(LEN=80) :: sub_name = 'destroy_environ_function'
        !
        !--------------------------------------------------------------------------------
        !
        IF (.NOT. ASSOCIATED(this%pos)) CALL io%destroy_error(sub_name)
        !
        !--------------------------------------------------------------------------------
        !
        NULLIFY (this%pos)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE destroy_environ_function
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE class_function
!----------------------------------------------------------------------------------------
