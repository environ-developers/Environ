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
    USE env_base_io, ONLY: ionode, environ_unit, global_verbose
    !
    USE environ_param, ONLY: DP
    !
    USE class_cell
    USE class_density
    USE class_function
    USE class_function_erfc
    USE class_function_exponential
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
    PUBLIC :: init_environ_functions, copy_environ_functions, &
              destroy_environ_functions, density_of_functions, gradient_of_functions, &
              laplacian_of_functions, hessian_of_functions, print_environ_functions
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
    SUBROUTINE init_environ_functions(f, fsrc, n, type_in, axis_in, dim_in, width_in, &
                                      spread_in, volume_in, pos_in)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_function), INTENT(IN) :: fsrc
        INTEGER, INTENT(IN) :: n
        INTEGER, INTENT(IN) :: type_in
        INTEGER, DIMENSION(n), INTENT(IN) :: dim_in, axis_in
        REAL(DP), DIMENSION(n), INTENT(IN) :: width_in, spread_in, volume_in
        REAL(DP), TARGET, INTENT(IN) :: pos_in(3, n)
        !
        CLASS(environ_function), ALLOCATABLE, INTENT(INOUT) :: f(:)
        !
        INTEGER :: i
        !
        CHARACTER(LEN=80) :: sub_name = 'init_environ_functions'
        !
        !--------------------------------------------------------------------------------
        !
        IF (ALLOCATED(f)) &
            CALL env_errore(sub_name, 'Trying to create an existing object', 1)
        !
        !--------------------------------------------------------------------------------
        !
        SELECT TYPE (fsrc)
            !
        TYPE IS (environ_function_gaussian)
            ALLOCATE (environ_function_gaussian :: f(n))
            !
        TYPE IS (environ_function_exponential)
            ALLOCATE (environ_function_exponential :: f(n))
            !
        TYPE IS (environ_function_erfc)
            ALLOCATE (environ_function_erfc :: f(n))
            !
        END SELECT
        !
        DO i = 1, n
            !
            CALL f(i)%init(type_in, axis_in(i), dim_in(i), width_in(i), spread_in(i), &
                           volume_in(i), pos_in(:, i))
            !
        END DO
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_environ_functions
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE copy_environ_functions(f, n, copy)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: n
        CLASS(environ_function), ALLOCATABLE, INTENT(IN) :: f(:)
        !
        CLASS(environ_function), ALLOCATABLE, INTENT(OUT) :: copy(:)
        !
        INTEGER :: i
        !
        CHARACTER(LEN=80) :: sub_name = 'copy_environ_functions'
        !
        !--------------------------------------------------------------------------------
        !
        IF (ALLOCATED(copy)) &
            CALL env_errore(sub_name, 'Trying to create an existing array', 1)
        !
        IF (.NOT. ALLOCATED(f)) &
            CALL env_errore(sub_name, 'Trying to copy an empty array', 1)
        !
        !--------------------------------------------------------------------------------
        !
        SELECT TYPE (f)
            !
        TYPE IS (environ_function_gaussian)
            ALLOCATE (environ_function_gaussian :: copy(n))
            !
        TYPE IS (environ_function_exponential)
            ALLOCATE (environ_function_exponential :: copy(n))
            !
        TYPE IS (environ_function_erfc)
            ALLOCATE (environ_function_erfc :: copy(n))
            !
        END SELECT
        !
        DO i = 1, n
            CALL f(i)%copy(copy(i))
        END DO
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE copy_environ_functions
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE destroy_environ_functions(f, n)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: n
        !
        CLASS(environ_function), ALLOCATABLE, INTENT(INOUT) :: f(:)
        !
        INTEGER :: i
        !
        CHARACTER(LEN=80) :: sub_name = 'destroy_environ_functions'
        !
        !--------------------------------------------------------------------------------
        !
        IF (.NOT. ALLOCATED(f)) &
            CALL env_errore(sub_name, 'Trying to destroy an empty array', 1)
        !
        IF (SIZE(f) /= n) &
            CALL env_errore(sub_name, 'Inconsistent size of allocated object', 1)
        !
        !--------------------------------------------------------------------------------
        !
        DO i = 1, n
            CALL f(i)%destroy()
        END DO
        !
        DEALLOCATE (f)
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
    SUBROUTINE density_of_functions(f, n, density, zero)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: n
        LOGICAL, INTENT(IN), OPTIONAL :: zero
        CLASS(environ_function), INTENT(IN) :: f(:)
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
        DO i = 1, n
            CALL f(i)%density(density)
        END DO
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE density_of_functions
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE gradient_of_functions(f, n, gradient, zero)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: n
        LOGICAL, INTENT(IN), OPTIONAL :: zero
        CLASS(environ_function), INTENT(IN) :: f(:)
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
        DO i = 1, n
            CALL f(i)%gradient(gradient)
        END DO
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE gradient_of_functions
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE laplacian_of_functions(f, n, laplacian, zero)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: n
        LOGICAL, INTENT(IN), OPTIONAL :: zero
        CLASS(environ_function), INTENT(IN) :: f(:)
        !
        TYPE(environ_density), INTENT(INOUT) :: laplacian
        !
        INTEGER :: i
        !
        CHARACTER(LEN=80) :: sub_name = 'laplacian_of_functions'
        !
        !--------------------------------------------------------------------------------
        !
        IF (PRESENT(zero)) THEN
            IF (zero) laplacian%of_r = 0.D0
        END IF
        !
        DO i = 1, n
            CALL f(i)%laplacian(laplacian)
        END DO
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE laplacian_of_functions
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE hessian_of_functions(f, n, hessian, zero)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: n
        LOGICAL, INTENT(IN), OPTIONAL :: zero
        CLASS(environ_function), INTENT(IN) :: f(:)
        !
        TYPE(environ_hessian), INTENT(INOUT) :: hessian
        !
        INTEGER :: i
        !
        CHARACTER(LEN=80) :: sub_name = 'hessian_of_functions'
        !
        !--------------------------------------------------------------------------------
        !
        IF (PRESENT(zero)) THEN
            IF (zero) hessian%of_r = 0.D0
        END IF
        !
        DO i = 1, n
            CALL f(i)%hessian(hessian)
        END DO
        !
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE hessian_of_functions
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
    !! @param unit          : (INTEGER) output target (default = environ_unit)
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE print_environ_functions(f, n, verbose, debug_verbose, unit)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: n
        CLASS(environ_function), INTENT(IN) :: f(:)
        INTEGER, INTENT(IN), OPTIONAL :: verbose, debug_verbose, unit
        !
        INTEGER :: base_verbose, local_verbose, local_unit, i
        !
        CHARACTER(LEN=80) :: sub_name = 'print_environ_functions'
        !
        !--------------------------------------------------------------------------------
        !
        IF (.NOT. ionode) RETURN
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
        ELSE IF (global_verbose > 0) THEN
            base_verbose = global_verbose
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
            local_unit = environ_unit
        END IF
        !
        IF (local_verbose >= 1) THEN
            !
            IF (local_verbose >= base_verbose) THEN ! header
                WRITE (local_unit, 1000)
            ELSE
                !
                CALL env_block_divider(ionode, local_verbose, base_verbose, local_unit)
                !
                WRITE (local_unit, 1001)
            END IF
            !
            WRITE (local_unit, 1002) ! legend
            WRITE (local_unit, 1003) ! table headers
            !
            DO i = 1, n
                !
                WRITE (local_unit, 1004) &
                    i, f(i)%f_type, f(i)%dim, f(i)%axis, f(i)%width, f(i)%spread, &
                    f(i)%volume, f(i)%pos
                !
            END DO
            !
            IF (local_verbose < base_verbose) &
                CALL env_block_divider(ionode, local_verbose, base_verbose, local_unit)
            !
        END IF
        !
        FLUSH (local_unit)
        !
        !--------------------------------------------------------------------------------
        !
1000    FORMAT(/, 4('%'), ' FUNCTIONS ', 65('%'))
1001    FORMAT(/, ' FUNCTIONS', /, ' =========')
        !
1002    FORMAT(/, ' 1 - Gaussian', /, &
                ' 2 - Complementary error function', /, &
                ' 3 - Exponential', /, &
                ' 4 - Scaled complementary error function', /, &
                ' 5 - Scaled error function')
        !
1003    FORMAT(/, '   i | type | dim | axis | width | spread |  volume  | position', /, &
                1X, 84('-'))
!
1004    FORMAT(1X, I3, ' | ', I4, ' | ', I3, ' | ', I4, ' | ', F5.3, ' | ', &
               F6.3, ' | ', F8.3, ' |', 3F10.5)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE print_environ_functions
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE class_functions
!----------------------------------------------------------------------------------------
