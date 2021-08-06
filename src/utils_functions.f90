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
    USE env_base_io, ONLY: ionode, environ_unit, verbose, depth
    !
    USE environ_param, ONLY: DP
    !
    USE class_cell
    USE class_density
    USE class_gradient
    USE class_hessian
    !
    USE generate_functions, ONLY: generate_gaussian, generate_erfc, &
                                  generate_exponential, generate_gradgaussian, &
                                  generate_graderfc, generate_gradexponential, &
                                  generate_laplerfc, generate_hesserfc, erfcvolume
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
    !>
    !!
    !------------------------------------------------------------------------------------
    TYPE, PUBLIC :: environ_function
        !--------------------------------------------------------------------------------
        !
        INTEGER :: f_type, axis, dim
        REAL(DP) :: width, spread, volume
        REAL(DP), POINTER :: pos(:)
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
        PROCEDURE :: density => density_of_function
        PROCEDURE :: gradient => gradient_of_function
        PROCEDURE :: laplacian => laplacian_of_function
        PROCEDURE :: hessian => hessian_of_function
        !
        !--------------------------------------------------------------------------------
    END TYPE environ_function
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
    SUBROUTINE create_environ_function(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_function), INTENT(INOUT) :: this
        !
        !--------------------------------------------------------------------------------
        !
        NULLIFY (this%pos)
        !
        RETURN
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
        REAL(DP), TARGET, INTENT(IN) :: pos(3)
        !
        CLASS(environ_function), INTENT(INOUT) :: this
        !
        INTEGER :: i
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
        this%pos => pos
        !
        RETURN
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
        TYPE(environ_function), INTENT(OUT) :: copy
        !
        !--------------------------------------------------------------------------------
        !
        copy%pos => this%pos
        !
        copy%f_type = this%f_type
        copy%dim = this%dim
        copy%axis = this%axis
        copy%spread = this%spread
        copy%width = this%width
        copy%volume = this%volume
        !
        RETURN
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
        IF (ASSOCIATED(this%pos)) NULLIFY (this%pos)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE destroy_environ_function
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
    SUBROUTINE density_of_function(this, density, zero)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_function), TARGET, INTENT(IN) :: this
        LOGICAL, INTENT(IN), OPTIONAL :: zero
        !
        TYPE(environ_density), TARGET, INTENT(INOUT) :: density
        !
        INTEGER :: i
        REAL(DP) :: local_charge
        !
        INTEGER, POINTER :: type_, dim, axis
        TYPE(environ_cell), POINTER :: cell
        REAL(DP), POINTER :: charge, spread, width
        REAL(DP), POINTER :: pos(:)
        !
        CHARACTER(LEN=80) :: sub_name = 'density_of_function'
        !
        !--------------------------------------------------------------------------------
        !
        IF (PRESENT(zero)) THEN
            IF (zero) density%of_r = 0.D0
        END IF
        !
        cell => density%cell
        !
        type_ => this%f_type
        pos => this%pos
        spread => this%spread
        charge => this%volume
        width => this%width
        dim => this%dim
        axis => this%axis
        !
        SELECT CASE (type_)
            !
        CASE (1)
            CALL generate_gaussian(dim, axis, charge, spread, pos, density) ! gaussian
            !
        CASE (2)
            !
            CALL generate_erfc(dim, axis, charge, width, spread, pos, density)
            ! CHARGE * NORMALIZED_ERFC_HALF(X); integrates to charge
            !
        CASE (3)
            !
            CALL generate_exponential(dim, axis, width, spread, pos, density)
            ! exponential
            !
        CASE (4)
            !
            !----------------------------------------------------------------------------
            ! CHARGE * NORMALIZED_ERFC_HALF(X) * VOLUME_NORMALIZED_ERFC_HALF
            ! Goes from charge to 0
            !
            local_charge = erfcvolume(dim, axis, width, spread, cell) * charge
            !
            CALL generate_erfc(dim, axis, local_charge, width, spread, pos, density)
            !
        CASE (5)
            !
            !----------------------------------------------------------------------------
            ! CHARGE * ( 1 - NORMALIZED_ERFC_HALF(x) * VOLUME_NORMALIZED_ERFC_HALF )
            ! goes from 0 to charge
            !
            local_charge = -erfcvolume(dim, axis, width, spread, cell) * charge
            !
            CALL generate_erfc(dim, axis, local_charge, width, spread, pos, density)
            !
            density%of_r = density%of_r + charge
            !
        END SELECT
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE density_of_function
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE gradient_of_function(this, gradient, zero)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        LOGICAL, INTENT(IN), OPTIONAL :: zero
        CLASS(environ_function), TARGET, INTENT(IN) :: this
        !
        TYPE(environ_gradient), TARGET, INTENT(INOUT) :: gradient
        !
        INTEGER :: i
        REAL(DP) :: local_charge
        !
        INTEGER, POINTER :: type_, dim, axis
        TYPE(environ_cell), POINTER :: cell
        REAL(DP), POINTER :: charge, spread, width
        REAL(DP), POINTER :: pos(:)
        !
        CHARACTER(LEN=80) :: sub_name = 'gradient_of_function'
        !
        !--------------------------------------------------------------------------------
        !
        cell => gradient%cell
        !
        IF (PRESENT(zero)) THEN
            IF (zero) gradient%of_r = 0.D0
        END IF
        !
        type_ => this%f_type
        pos => this%pos
        spread => this%spread
        charge => this%volume
        width => this%width
        dim => this%dim
        axis => this%axis
        !
        SELECT CASE (type_)
            !
        CASE (1)
            !
            CALL generate_gradgaussian(dim, axis, charge, spread, pos, gradient)
            ! gaussian
            !
        CASE (2)
            !
            CALL generate_graderfc(dim, axis, charge, width, spread, pos, gradient)
            ! CHARGE * NORMALIZED_ERFC_HALF(X) ! integrates to charge
            !
        CASE (3)
            !
            CALL generate_gradexponential(dim, axis, width, spread, pos, gradient)
            ! exponential
            !
        CASE (4)
            !
            !----------------------------------------------------------------------------
            ! CHARGE * NORMALIZED_ERFC_HALF(X) * VOLUME_NORMALIZED_ERFC_HALF
            ! goes from charge to 0
            !
            local_charge = erfcvolume(dim, axis, width, spread, cell) * charge
            !
            CALL generate_graderfc(dim, axis, local_charge, width, spread, pos, gradient)
            !
        CASE (5)
            !
            !----------------------------------------------------------------------------
            ! CHARGE * ( 1 - NORMALIZED_ERFC_HALF(x) * VOLUME_NORMALIZED_ERFC_HALF )
            ! goes from 0 to charge
            !
            local_charge = -erfcvolume(dim, axis, width, spread, cell) * charge
            !
            CALL generate_graderfc(dim, axis, local_charge, width, spread, pos, gradient)
            !
        END SELECT
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE gradient_of_function
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE laplacian_of_function(this, laplacian, zero)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        LOGICAL, INTENT(IN), OPTIONAL :: zero
        CLASS(environ_function), TARGET, INTENT(IN) :: this
        !
        TYPE(environ_density), TARGET, INTENT(INOUT) :: laplacian
        !
        INTEGER :: i
        REAL(DP) :: local_charge
        !
        INTEGER, POINTER :: type_, dim, axis
        TYPE(environ_cell), POINTER :: cell
        REAL(DP), POINTER :: charge, spread, width
        REAL(DP), POINTER :: pos(:)
        !
        CHARACTER(LEN=80) :: sub_name = 'laplacian_of_function'
        !
        !--------------------------------------------------------------------------------
        !
        cell => laplacian%cell
        !
        IF (PRESENT(zero)) THEN
            IF (zero) laplacian%of_r = 0.D0
        END IF
        !
        type_ => this%f_type
        pos => this%pos
        spread => this%spread
        charge => this%volume
        width => this%width
        dim => this%dim
        axis => this%axis
        !
        SELECT CASE (type_)
            !
        CASE (1)
            CALL env_errore(sub_name, 'Options not yet implemented', 1) ! gaussian
            !
        CASE (2)
            !
            CALL generate_laplerfc(dim, axis, charge, width, spread, pos, laplacian)
            ! CHARGE * NORMALIZED_ERFC_HALF(X) ! integrates to charge
            !
        CASE (3)
            CALL env_errore(sub_name, 'Options not yet implemented', 1) ! exponential
            !
        CASE (4)
            !
            !----------------------------------------------------------------------------
            ! CHARGE * NORMALIZED_ERFC_HALF(X) * VOLUME_NORMALIZED_ERFC_HALF
            ! goes from charge to 0
            !
            local_charge = erfcvolume(dim, axis, width, spread, cell) * charge
            !
            CALL generate_laplerfc(dim, axis, local_charge, width, spread, &
                                   pos, laplacian)
            !
        CASE (5)
            !
            !----------------------------------------------------------------------------
            ! CHARGE * ( 1 - NORMALIZED_ERFC_HALF(x) * VOLUME_NORMALIZED_ERFC_HALF )
            ! goes from 0 to charge
            !
            local_charge = -erfcvolume(dim, axis, width, spread, cell) * charge
            !
            CALL generate_laplerfc(dim, axis, local_charge, width, spread, &
                                   pos, laplacian)
            !
        END SELECT
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE laplacian_of_function
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE hessian_of_function(this, hessian, zero)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        LOGICAL, INTENT(IN), OPTIONAL :: zero
        CLASS(environ_function), TARGET, INTENT(IN) :: this
        !
        TYPE(environ_hessian), TARGET, INTENT(INOUT) :: hessian
        !
        INTEGER :: i
        REAL(DP) :: local_charge
        !
        INTEGER, POINTER :: type_, dim, axis
        TYPE(environ_cell), POINTER :: cell
        REAL(DP), POINTER :: charge, spread, width
        REAL(DP), POINTER :: pos(:)
        !
        CHARACTER(LEN=80) :: sub_name = 'hessian_of_function'
        !
        !--------------------------------------------------------------------------------
        !
        cell => hessian%cell
        !
        IF (PRESENT(zero)) THEN
            IF (zero) hessian%of_r = 0.D0
        END IF
        !
        type_ => this%f_type
        pos => this%pos
        spread => this%spread
        charge => this%volume
        width => this%width
        dim => this%dim
        axis => this%axis
        !
        SELECT CASE (type_)
            !
        CASE (1)
            CALL env_errore(sub_name, 'Options not yet implemented', 1) ! gaussian
            !
        CASE (2)
            !
            CALL generate_hesserfc(dim, axis, charge, width, spread, pos, hessian)
            ! CHARGE * NORMALIZED_ERFC_HALF(X) ! integrates to charge
            !
        CASE (3)
            CALL env_errore(sub_name, 'Options not yet implemented', 1) ! exponential
            !
        CASE (4)
            !
            !----------------------------------------------------------------------------
            ! CHARGE * NORMALIZED_ERFC_HALF(X) * VOLUME_NORMALIZED_ERFC_HALF
            ! goes from charge to 0
            !
            local_charge = erfcvolume(dim, axis, width, spread, cell) * charge
            !
            CALL generate_hesserfc(dim, axis, local_charge, width, spread, pos, hessian)
            !
        CASE (5)
            !
            !----------------------------------------------------------------------------
            ! CHARGE * ( 1 - NORMALIZED_ERFC_HALF(x) * VOLUME_NORMALIZED_ERFC_HALF )
            ! goes from 0 to charge
            !
            local_charge = -erfcvolume(dim, axis, width, spread, cell) * charge
            !
            CALL generate_hesserfc(dim, axis, local_charge, width, spread, pos, hessian)
            !
        END SELECT
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE hessian_of_function
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !                                  ARRAY ROUTINES
    !
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_environ_functions(f, n, type_in, axis, dim_in, width, spread_in, &
                                      volume_in, pos)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: n, type_in
        INTEGER, DIMENSION(n), INTENT(IN) :: dim_in, axis
        REAL(DP), DIMENSION(n), INTENT(IN) :: width, spread_in, volume_in
        REAL(DP), TARGET, INTENT(IN) :: pos(3, n)
        !
        TYPE(environ_function), ALLOCATABLE, INTENT(INOUT) :: f(:)
        !
        INTEGER :: i
        !
        CHARACTER(LEN=80) :: sub_name = 'init_environ_functions'
        !
        !--------------------------------------------------------------------------------
        !
        IF (ALLOCATED(f)) RETURN
        !
        ALLOCATE (f(n))
        !
        DO i = 1, n
            !
            CALL f(i)%init(type_in, axis(i), dim_in(i), width(i), spread_in(i), &
                           volume_in(i), pos(:, i))
            !
        END DO
        !
        RETURN
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
        TYPE(environ_function), INTENT(IN) :: f(n)
        !
        TYPE(environ_function), ALLOCATABLE, INTENT(OUT) :: copy(:)
        !
        INTEGER :: i
        !
        !--------------------------------------------------------------------------------
        !
        ALLOCATE (copy(n))
        !
        DO i = 1, n
            CALL f(i)%copy(copy(i))
        END DO
        !
        RETURN
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
        TYPE(environ_function), ALLOCATABLE, INTENT(INOUT) :: f(:)
        !
        INTEGER :: i
        !
        CHARACTER(LEN=80) :: sub_name = 'destroy_environ_functions'
        !
        !--------------------------------------------------------------------------------
        !
        IF (SIZE(f) /= n) &
            CALL env_errore(sub_name, 'Inconsistent size of allocated object', 1)
        !
        DO i = 1, n
            CALL f(i)%destroy()
        END DO
        !
        DEALLOCATE (f)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE destroy_environ_functions
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
        TYPE(environ_function), INTENT(IN) :: f(n)
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
        RETURN
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
        TYPE(environ_function), INTENT(IN) :: f(n)
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
        RETURN
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
        TYPE(environ_function), INTENT(IN) :: f(n)
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
        DO i = 1, n
            CALL f(i)%laplacian(laplacian)
        END DO
        !
        RETURN
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
        TYPE(environ_function), INTENT(IN) :: f(n)
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
        DO i = 1, n
            CALL f(i)%hessian(hessian)
        END DO
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE hessian_of_functions
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE print_environ_functions(f, n, local_verbose, local_depth)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: n
        TYPE(environ_function), INTENT(IN) :: f(n)
        INTEGER, INTENT(IN), OPTIONAL :: local_verbose, local_depth
        !
        INTEGER :: i, verbosity, passed_verbosity, passed_depth
        !
        CHARACTER(LEN=80) :: sub_name = 'print_environ_functions'
        !
        !--------------------------------------------------------------------------------
        !
        IF (.NOT. ionode .OR. verbose == 0) RETURN
        !
        IF (PRESENT(local_verbose)) THEN
            verbosity = verbose + local_verbose
        ELSE
            verbosity = verbose
        END IF
        !
        IF (verbosity == 0) RETURN
        !
        IF (verbosity >= 1) THEN
            !
            IF (verbosity >= verbose) THEN ! header
                WRITE (environ_unit, 1000)
            ELSE
                !
                CALL env_block_divider(verbosity)
                !
                WRITE (environ_unit, 1001)
            END IF
            !
            WRITE (environ_unit, 1002) ! legend
            WRITE (environ_unit, 1003) ! table headers
            !
            DO i = 1, n
                !
                WRITE (environ_unit, 1004) &
                    i, f(i)%f_type, f(i)%dim, f(i)%axis, f(i)%width, f(i)%spread, &
                    f(i)%volume, f(i)%pos
                !
            END DO
            !
            IF (verbosity < verbose) CALL env_block_divider(verbosity)
            !
        END IF
        !
        FLUSH (environ_unit)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
        !
1000    FORMAT(/, 4('%'), ' FUNCTIONS ', 65('%'))
1001    FORMAT(/, ' FUNCTIONS', /, ' =========')
        !
1002    FORMAT(/, ' 1 - Gaussian', /, &
                ' 2 - Complimentory error function', /, &
                ' 3 - Exponential', /, &
                ' 4 - Scaled complimentory error function', /, &
                ' 5 - Scaled error function')
        !
1003    FORMAT(/, '   i | type | dim | axis | width | spread |  volume  | position', /, &
                1X, 83('-'))
!
1004    FORMAT(1X, I3, ' | ', I4, ' | ', I3, ' | ', I4, ' | ', F5.3, ' | ', &
               F6.3, ' | ', F8.3, ' |', 3F10.7)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE print_environ_functions
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE class_functions
!----------------------------------------------------------------------------------------
