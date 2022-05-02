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
!          Francesco Nattino  (THEOS and NCCR-MARVEL, EPFL)
!          Nicola Marzari     (THEOS and NCCR-MARVEL, EPFL)
!          Edan Bainglass     (Department of Physics, UNT)
!
!----------------------------------------------------------------------------------------
!>
!! Environ_boundary contains all the specifications and the details of
!! the smooth interface between the QM and the continuum regions of the
!! simulation cell. The main interface function is stored in the %scaled
!! component, the type also stores boundary real-space derivatives (gradient,
!! laplacian, dsurface, hessian) and other quantities needed by Environ
!! modules.
!!
!----------------------------------------------------------------------------------------
MODULE class_boundary_system
    !------------------------------------------------------------------------------------
    !
    USE class_io, ONLY: io
    !
    USE environ_param, ONLY: DP
    !
    USE class_density
    USE class_gradient
    USE class_hessian
    !
    USE class_boundary
    USE class_function_erfc
    USE class_system
    !
    USE boundary_tools
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
    TYPE, EXTENDS(environ_boundary), PUBLIC :: environ_boundary_system
        !--------------------------------------------------------------------------------
        !
        TYPE(environ_system), POINTER :: system => NULL()
        !
        TYPE(environ_function_erfc) :: simple ! components needed for boundary of system
        !
        !--------------------------------------------------------------------------------
    CONTAINS
        !--------------------------------------------------------------------------------
        !
        PROCEDURE :: create => create_environ_boundary
        PROCEDURE :: init => init_environ_boundary
        PROCEDURE :: update => update_environ_boundary
        PROCEDURE :: destroy => destroy_environ_boundary
        PROCEDURE :: build => boundary_of_system
        !
        PROCEDURE :: dboundary_dions => calc_dboundary_dions
        !
        PROCEDURE :: printout => print_environ_boundary
        !
        !--------------------------------------------------------------------------------
    END TYPE environ_boundary_system
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
    SUBROUTINE create_environ_boundary(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_boundary_system), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: routine = 'create_environ_boundary'
        !
        !--------------------------------------------------------------------------------
        !
        IF (ASSOCIATED(this%system)) CALL io%create_error(routine)
        !
        !--------------------------------------------------------------------------------
        !
        NULLIFY (this%system)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE create_environ_boundary
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_environ_boundary(this, system_distance, system_spread, system)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_system), TARGET, INTENT(IN) :: system
        !
        REAL(DP), INTENT(IN) :: system_distance
        REAL(DP), INTENT(IN) :: system_spread
        !
        CLASS(environ_boundary_system), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: routine = 'init_environ_boundary'
        !
        !--------------------------------------------------------------------------------
        !
        CALL this%create()
        !
        this%system => system
        !
        CALL this%simple%init(3, system%axis, system%dim, system_distance, &
                              system_spread, 1.D0, system%com)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_environ_boundary
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE update_environ_boundary(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_boundary_system), INTENT(INOUT) :: this
        !
        LOGICAL :: update_anything
        !
        CHARACTER(LEN=80) :: routine = 'update_environ_boundary'
        !
        !--------------------------------------------------------------------------------
        !
        update_anything = .FALSE.
        !
        update_anything = update_anything .OR. this%system%lupdate
        !
        IF (.NOT. update_anything) THEN
            !
            IF (this%update_status == 2) this%update_status = 0
            ! nothing is under update, change update_status and exit
            !
            RETURN
            !
        END IF
        !
        IF (this%system%lupdate) THEN
            !
            !----------------------------------------------------------------------------
            ! Only ions are needed, fully update the boundary
            !
            CALL this%build()
            !
            ! TO DEBUG SOLVENT-AWARE
            ! !
            ! CALL this%invert()
            ! !
            ! CALL test_de_dboundary(this)
            !
            this%update_status = 2 ! boundary has changed and is ready
        ELSE
            !
            IF (this%update_status == 2) this%update_status = 0
            ! boundary has not changed
            !
            RETURN
            !
        END IF
        !
        CALL this%update_solvent_aware() ! update solvent aware if applicable
        !
        !--------------------------------------------------------------------------------
        ! Output current state
        !
        IF (this%update_status == 2) CALL this%printout()
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE update_environ_boundary
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE destroy_environ_boundary(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_boundary_system), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: routine = 'destroy_environ_boundary'
        !
        !--------------------------------------------------------------------------------
        !
        CALL this%pre_destroy()
        !
        CALL this%simple%destroy()
        !
        NULLIFY (this%system)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE destroy_environ_boundary
    !------------------------------------------------------------------------------------
    !>
    !! Updates the boundary using a function
    !!
    !! Calculates the dielectric constant as a function of the charge
    !! density, and the derivatives of the dielectric constant with
    !! respect to the charge density. Also updates the volume and surface
    !! components. Expects an explicity defined system density function.
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE boundary_of_system(this, density)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_density), OPTIONAL, TARGET, INTENT(IN) :: density
        !
        CLASS(environ_boundary_system), TARGET, INTENT(INOUT) :: this
        !
        TYPE(environ_hessian), POINTER :: hessloc
        !
        CHARACTER(LEN=80) :: routine = 'boundary_of_system'
        !
        !--------------------------------------------------------------------------------
        !
        ASSOCIATE (cell => this%scaled%cell, &
                   deriv => this%deriv, &
                   derivatives => this%cores%derivatives, &
                   simple => this%simple, &
                   scal => this%scaled, &
                   grad => this%gradient, &
                   lapl => this%laplacian, &
                   dsurf => this%dsurface)
            !
            !----------------------------------------------------------------------------
            !
            CALL simple%density(scal, .TRUE.)
            ! compute soft spheres and generate boundary
            !
            !----------------------------------------------------------------------------
            ! Generate boundary derivatives, if needed
            !
            IF (deriv >= 3) THEN
                !
                IF (this%solvent_aware) THEN
                    hessloc => this%hessian
                ELSE
                    ALLOCATE (hessloc)
                    !
                    CALL hessloc%init(cell)
                    !
                END IF
                !
            END IF
            !
            SELECT CASE (this%derivatives_method)
                !
            CASE ('fft')
                !
                IF (deriv == 1 .OR. deriv == 2) CALL derivatives%gradient(scal, grad)
                !
                IF (deriv == 2) CALL derivatives%laplacian(scal, lapl)
                !
                IF (deriv == 3) CALL this%calc_dsurface(scal, grad, lapl, hessloc, dsurf)
                !
            CASE ('chain')
                !
                IF (deriv >= 1) CALL simple%gradient(grad, .TRUE.)
                !
                IF (deriv >= 2) CALL simple%laplacian(lapl, .TRUE.)
                !
                IF (deriv >= 3) THEN
                    !
                    CALL simple%hessian(hessloc, .TRUE.)
                    !
                    CALL calc_dsurface_no_pre(cell, grad%of_r, hessloc%of_r, dsurf%of_r)
                    !
                END IF
                !
            CASE DEFAULT
                CALL io%error(routine, "Unexpected derivatives method", 1)
                !
            END SELECT
            !
            IF (deriv >= 3) THEN
                !
                IF (.NOT. this%solvent_aware) THEN
                    !
                    CALL hessloc%destroy()
                    !
                    DEALLOCATE (hessloc)
                END IF
                !
            END IF
            !
            this%volume = scal%integrate()
            !
            IF (deriv >= 1) THEN
                !
                CALL grad%update_modulus()
                !
                this%surface = grad%modulus%integrate()
            END IF
            !
        END ASSOCIATE
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE boundary_of_system
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
    SUBROUTINE calc_dboundary_dions(this, index, partial)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_boundary_system), TARGET, INTENT(IN) :: this
        INTEGER, INTENT(IN) :: index
        !
        TYPE(environ_gradient), INTENT(INOUT) :: partial
        !
        CHARACTER(LEN=80) :: routine = 'calc_dboundary_dions'
        !
        !--------------------------------------------------------------------------------
        !
        IF (index > this%system%ions%number) &
            CALL io%error(routine, "Index greater than number of ions", 1)
        !
        IF (index <= 0) CALL io%error(routine, "Index of ion is zero or lower", 1)
        !
        ! PROBABLY THERE IS A UNIFORM CONTRIBUTION TO THE FORCES
        ! WHICH SHOULD ONLY AFFECT THE COM OF THE SYSTEM, POSSIBLY NEED TO ADD
        ! A CHECK ON ATOMS THAT BELONG TO THE SYSTEM
        !
        partial%of_r = 0.D0
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_dboundary_dions
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !                                   OUTPUT METHODS
    !
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !>
    !! Prints the details of the boundary
    !!
    !! Nested objects receive a decremented passed verbose to trigger block printing
    !!
    !! @param verbose       : (INTEGER) adds verbosity to global verbose
    !! @param debug_verbose : (INTEGER) replaces global verbose for debugging
    !! @param unit          : (INTEGER) output target (default = io%debug_unit)
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE print_environ_boundary(this, verbose, debug_verbose, unit)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_boundary_system), INTENT(IN) :: this
        INTEGER, OPTIONAL, INTENT(IN) :: verbose, debug_verbose, unit
        !
        INTEGER :: base_verbose, local_verbose, passed_verbose, local_unit
        !
        CHARACTER(LEN=80) :: routine = 'print_environ_boundary'
        !
        !--------------------------------------------------------------------------------
        !
        IF (.NOT. PRESENT(debug_verbose) .AND. io%verbosity <= 0) RETURN
        !
        CALL this%pre_printout(verbose, debug_verbose, unit)
        !
        CALL this%print_setup(base_verbose, local_verbose, passed_verbose, local_unit, &
                              verbose, debug_verbose, unit)
        !
        !--------------------------------------------------------------------------------
        !
        IF (local_verbose >= 1) THEN
            !
            IF (io%lnode) THEN
                !
                WRITE (local_unit, 1100) &
                    this%simple%pos, this%simple%width, this%simple%spread, &
                    this%simple%dim, this%simple%axis
                !
            END IF
            !
        END IF
        !
        FLUSH (local_unit)
        !
        !--------------------------------------------------------------------------------
        !
1100    FORMAT(/, " boundary is built as an analytic function centered on system position:", /, &
                " center of the boundary     = ", 3F14.7, /, &
                " distance from the center   = ", F14.7, /, &
                " spread of the interface    = ", F14.7, /, &
                " dimensionality             = ", I14, /, &
                " axis                       = ", I14)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE print_environ_boundary
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE class_boundary_system
!----------------------------------------------------------------------------------------
