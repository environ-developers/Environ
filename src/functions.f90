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
MODULE class_functions
    !------------------------------------------------------------------------------------
    !
    USE class_io, ONLY: io
    !
    USE environ_param, ONLY: DP
    !
    USE class_density
    USE class_function
    USE class_function_erfc
    USE class_function_gaussian
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
    TYPE, PUBLIC :: environ_functions
        !--------------------------------------------------------------------------------
        !
        INTEGER :: number = 0
        INTEGER :: f_type
        !
        CLASS(environ_function), ALLOCATABLE :: array(:)
        !
        !--------------------------------------------------------------------------------
    CONTAINS
        !--------------------------------------------------------------------------------
        !
        PROCEDURE, PRIVATE :: create => create_environ_functions
        PROCEDURE :: init => init_environ_functions
        PROCEDURE :: destroy => destroy_environ_functions
        !
        PROCEDURE :: density => density_of_functions
        PROCEDURE :: gradient => gradient_of_functions
        PROCEDURE :: laplacian => laplacian_of_functions
        PROCEDURE :: hessian => hessian_of_functions
        PROCEDURE :: derivative => derivative_of_functions
        !
        PROCEDURE :: printout => print_environ_functions
        !
        !--------------------------------------------------------------------------------
    END TYPE environ_functions
    !------------------------------------------------------------------------------------
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
    SUBROUTINE create_environ_functions(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_functions), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: sub_name = 'create_environ_functions'
        !
        !--------------------------------------------------------------------------------
        !
        IF (ALLOCATED(this%array)) CALL io%create_error(sub_name)
        !
        !--------------------------------------------------------------------------------
        !
        this%number = 0
        this%f_type = 0
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE create_environ_functions
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_environ_functions(this, n, f_type, f_axis, f_dim, f_width, &
                                      f_spread, f_volume, f_pos)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: n
        INTEGER, INTENT(IN) :: f_type
        INTEGER, DIMENSION(n), INTENT(IN) :: f_dim, f_axis
        REAL(DP), DIMENSION(n), INTENT(IN) :: f_width, f_spread, f_volume
        REAL(DP), OPTIONAL, TARGET, INTENT(IN) :: f_pos(3, n)
        !
        CLASS(environ_functions), INTENT(INOUT) :: this
        !
        INTEGER :: i
        !
        CHARACTER(LEN=80) :: sub_name = 'init_environ_functions'
        !
        !--------------------------------------------------------------------------------
        !
        CALL this%create()
        !
        !--------------------------------------------------------------------------------
        ! Cast function as concrete type
        !
        SELECT CASE (f_type)
            !
        CASE (1)
            ALLOCATE (environ_function_gaussian :: this%array(n))
            !
        CASE (2, 3, 4)
            ALLOCATE (environ_function_erfc :: this%array(n))
            !
        CASE DEFAULT
            CALL io%error(sub_name, "Unexpected function type", 1)
            !
        END SELECT
        !
        !--------------------------------------------------------------------------------
        !
        this%number = n
        this%f_type = f_type
        !
        DO i = 1, this%number
            !
            CALL this%array(i)%init(f_type, f_axis(i), f_dim(i), f_width(i), &
                                    f_spread(i), f_volume(i), f_pos(:, i))
            !
        END DO
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_environ_functions
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE destroy_environ_functions(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_functions), INTENT(INOUT) :: this
        !
        INTEGER :: i
        !
        CHARACTER(LEN=80) :: sub_name = 'destroy_environ_functions'
        !
        !--------------------------------------------------------------------------------
        !
        IF (.NOT. ALLOCATED(this%array)) CALL io%destroy_error(sub_name)
        !
        !--------------------------------------------------------------------------------
        !
        DO i = 1, this%number
            CALL this%array(i)%destroy()
        END DO
        !
        DEALLOCATE (this%array)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE destroy_environ_functions
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
    SUBROUTINE density_of_functions(this, density, zero)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_functions), INTENT(INOUT) :: this
        LOGICAL, OPTIONAL, INTENT(IN) :: zero
        !
        TYPE(environ_density), INTENT(INOUT) :: density
        !
        INTEGER :: i
        !
        !--------------------------------------------------------------------------------
        !
        IF (PRESENT(zero)) THEN
            IF (zero) density%of_r = 0.D0
        END IF
        !
        DO i = 1, this%number
            CALL this%array(i)%density(density)
        END DO
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE density_of_functions
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE gradient_of_functions(this, gradient, zero)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_functions), INTENT(IN) :: this
        LOGICAL, OPTIONAL, INTENT(IN) :: zero
        !
        TYPE(environ_gradient), INTENT(INOUT) :: gradient
        !
        INTEGER :: i
        !
        !--------------------------------------------------------------------------------
        !
        IF (PRESENT(zero)) THEN
            IF (zero) gradient%of_r = 0.D0
        END IF
        !
        DO i = 1, this%number
            CALL this%array(i)%gradient(gradient)
        END DO
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE gradient_of_functions
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE laplacian_of_functions(this, laplacian, zero)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_functions), INTENT(IN) :: this
        LOGICAL, OPTIONAL, INTENT(IN) :: zero
        !
        TYPE(environ_density), INTENT(INOUT) :: laplacian
        !
        INTEGER :: i
        !
        !--------------------------------------------------------------------------------
        !
        IF (PRESENT(zero)) THEN
            IF (zero) laplacian%of_r = 0.D0
        END IF
        !
        DO i = 1, this%number
            CALL this%array(i)%laplacian(laplacian)
        END DO
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE laplacian_of_functions
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE hessian_of_functions(this, hessian, zero)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_functions), INTENT(IN) :: this
        LOGICAL, OPTIONAL, INTENT(IN) :: zero
        !
        TYPE(environ_hessian), INTENT(INOUT) :: hessian
        !
        INTEGER :: i
        !
        !--------------------------------------------------------------------------------
        !
        IF (PRESENT(zero)) THEN
            IF (zero) hessian%of_r = 0.D0
        END IF
        !
        DO i = 1, this%number
            CALL this%array(i)%hessian(hessian)
        END DO
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE hessian_of_functions
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE derivative_of_functions(this, derivative, zero)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_functions), INTENT(IN) :: this
        LOGICAL, OPTIONAL, INTENT(IN) :: zero
        !
        TYPE(environ_density), INTENT(INOUT) :: derivative
        !
        INTEGER :: i
        !
        !--------------------------------------------------------------------------------
        !
        IF (PRESENT(zero)) THEN
            IF (zero) derivative%of_r = 0.D0
        END IF
        !
        DO i = 1, this%number
            CALL this%array(i)%derivative(derivative)
        END DO
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE derivative_of_functions
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !                                   OUTPUT METHODS
    !
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !>
    !! Prints the details of the functions
    !!
    !! If called by a parent object, prints details in block format
    !!
    !! @param verbose       : (INTEGER) adds verbosity to global verbose
    !! @param debug_verbose : (INTEGER) replaces global verbose for debugging
    !! @param unit          : (INTEGER) output target (default = io%debug_unit)
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE print_environ_functions(this, verbose, debug_verbose, unit)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_functions), INTENT(IN) :: this
        INTEGER, OPTIONAL, INTENT(IN) :: verbose, debug_verbose, unit
        !
        INTEGER :: base_verbose, local_verbose, local_unit, i
        !
        CHARACTER(LEN=80) :: sub_name = 'print_environ_functions'
        !
        !--------------------------------------------------------------------------------
        !
        IF (.NOT. io%lnode) RETURN
        !
        IF (PRESENT(debug_verbose)) THEN
            base_verbose = debug_verbose
            !
            IF (PRESENT(verbose)) THEN
                local_verbose = verbose
            ELSE
                local_verbose = debug_verbose
            END IF
            !
        ELSE IF (io%verbosity > 0) THEN
            base_verbose = io%verbosity
            !
            IF (PRESENT(verbose)) THEN
                local_verbose = base_verbose + verbose
            ELSE
                local_verbose = base_verbose
            END IF
            !
        ELSE
            RETURN
        END IF
        !
        IF (PRESENT(unit)) THEN
            local_unit = unit
        ELSE
            local_unit = io%debug_unit
        END IF
        !
        IF (local_verbose >= 1) THEN
            !
            IF (local_verbose >= base_verbose) THEN ! header
                WRITE (local_unit, 1000)
            ELSE
                !
                CALL io%block_divider(local_verbose, base_verbose, local_unit)
                !
                WRITE (local_unit, 1001)
            END IF
            !
            WRITE (local_unit, 1002) ! legend
            WRITE (local_unit, 1003) ! table headers
            !
            DO i = 1, this%number
                !
                ASSOCIATE (array => this%array)
                    !
                    WRITE (local_unit, 1004) &
                        i, array(i)%f_type, array(i)%dim, array(i)%axis, &
                        array(i)%width, array(i)%spread, array(i)%volume, array(i)%pos
                    !
                END ASSOCIATE
                !
            END DO
            !
            IF (local_verbose < base_verbose) &
                CALL io%block_divider(local_verbose, base_verbose, local_unit)
            !
        END IF
        !
        FLUSH (local_unit)
        !
        !--------------------------------------------------------------------------------
        !
1000    FORMAT(/, 4('%'), " FUNCTIONS ", 65('%'))
1001    FORMAT(/, " FUNCTIONS", /, " =========")
        !
1002    FORMAT(/, " 1 - Gaussian", /, &
                " 2 - Complementary error function", /, &
                " 3 - Scaled complementary error function", /, &
                " 4 - Scaled error function")
        !
1003    FORMAT(/, "   i | type | dim | axis | width | spread |  volume  | position", /, &
                1X, 84('-'))
!
1004    FORMAT(1X, I3, " | ", I4, " | ", I3, " | ", I4, " | ", F5.3, " | ", &
               F6.3, " | ", F8.3, " |", 3F10.5)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE print_environ_functions
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE class_functions
!----------------------------------------------------------------------------------------
