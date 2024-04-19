!----------------------------------------------------------------------------------------
!
! Copyright (C) 2018-2022 ENVIRON (www.quantum-environ.org)
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
MODULE class_boundary_electronic
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
    USE class_electrons
    USE class_ions
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
    TYPE, EXTENDS(environ_boundary), PUBLIC :: environ_boundary_electronic
        !--------------------------------------------------------------------------------
        !
        TYPE(environ_electrons), POINTER :: electrons => NULL()
        !
        TYPE(environ_ions), POINTER :: ions => NULL() ! if including core correction
        !
        REAL(DP) :: rhomax
        REAL(DP) :: rhomin
        REAL(DP) :: fact
        !
        TYPE(environ_density) :: density
        !
        TYPE(environ_density) :: dscaled
        TYPE(environ_density) :: d2scaled
        !
        !--------------------------------------------------------------------------------
    CONTAINS
        !--------------------------------------------------------------------------------
        !
        PROCEDURE :: create => create_environ_boundary
        PROCEDURE :: init => init_environ_boundary
        PROCEDURE :: update => update_environ_boundary
        PROCEDURE :: destroy => destroy_environ_boundary
        PROCEDURE :: build => boundary_of_density
        !
        PROCEDURE :: dboundary_dions => calc_dboundary_dions
        !
        PROCEDURE :: printout => print_environ_boundary
        !
        !--------------------------------------------------------------------------------
    END TYPE environ_boundary_electronic
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
        CLASS(environ_boundary_electronic), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: routine = 'create_environ_boundary'
        !
        !--------------------------------------------------------------------------------
        !
        IF (ASSOCIATED(this%electrons)) CALL io%create_error(routine)
        !
        IF (ASSOCIATED(this%ions)) CALL io%create_error(routine)
        !
        !--------------------------------------------------------------------------------
        !
        this%rhomax = 0.D0
        this%rhomin = 0.D0
        this%fact = 0.D0
        !
        NULLIFY (this%electrons)
        NULLIFY (this%ions)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE create_environ_boundary
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_environ_boundary(this, rhomax, rhomin, electrons, ions)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        !
        TYPE(environ_electrons), TARGET, INTENT(IN) :: electrons
        TYPE(environ_ions), OPTIONAL, TARGET, INTENT(IN) :: ions
        !
        REAL(DP), INTENT(IN) :: rhomax
        REAL(DP), INTENT(IN) :: rhomin
        !
        CLASS(environ_boundary_electronic), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: routine = 'init_environ_boundary'
        !
        !--------------------------------------------------------------------------------
        !
        CALL this%create()
        !
        this%rhomax = rhomax
        this%rhomin = rhomin
        this%fact = LOG(rhomax / rhomin)
        !
        this%electrons => electrons
        !
        IF (this%mode == 'full') THEN
            !
            IF (.NOT. PRESENT(ions)) CALL io%error(routine, "Missing ions", 1)

            this%ions => ions
        END IF
        !
        CALL this%density%init(this%cell, 'boundary_density_'//this%label)
        !
        CALL this%dscaled%init(this%cell, 'dboundary_'//this%label)
        !
        CALL this%d2scaled%init(this%cell, 'd2boundary_'//this%label)
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
        CLASS(environ_boundary_electronic), INTENT(INOUT) :: this
        !
        LOGICAL :: update_anything
        !
        CHARACTER(LEN=80) :: routine = 'update_environ_boundary'
        !
        !--------------------------------------------------------------------------------
        !
        update_anything = .FALSE.
        !
        IF (ASSOCIATED(this%ions)) update_anything = this%ions%lupdate
        !
        update_anything = update_anything .OR. this%electrons%lupdate
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
        SELECT CASE (this%mode)
            !
        CASE ('full')
            !
            IF (this%ions%lupdate) THEN
                !
                !------------------------------------------------------------------------
                ! Compute the ionic part
                !
                CALL this%ions%core_electrons%density(this%ions%core, .TRUE.)
                !
                this%update_status = 1 ! waiting to finish update
            END IF
            !
            IF (this%electrons%lupdate) THEN
                !
                !------------------------------------------------------------------------
                ! Check if the ionic part has been updated
                !
                IF (this%update_status == 0) &
                    CALL io%error(routine, &
                                  "Wrong update status, possibly missing ionic update", 1)
                !
                this%density%of_r = this%electrons%density%of_r + this%ions%core%of_r
                !
                CALL this%build()
                !
                this%update_status = 2 ! boundary has changed and is ready
            END IF
            !
        CASE ('electronic')
            !
            IF (this%electrons%lupdate) THEN
                this%density%of_r = this%electrons%density%of_r
                !
                IF ( this%electrons%use_local_gradient ) THEN
                    this%derivatives_method = 'chain-local'
                    this%gradient%of_r = this%electrons%gradient%of_r
                ENDIF
                !
                CALL this%build()
                !
                this%update_status = 2 ! boundary has changes and is ready
                !
                ! CALL test_energy_derivatives(1, this) ! DEBUGGING
                !
            ELSE
                !
                IF (this%update_status == 2) this%update_status = 0
                ! boundary has not changed
                !
                RETURN
                !
            END IF
            !
        CASE DEFAULT
            CALL io%error(routine, "Unrecognized boundary mode", 1)
            !
        END SELECT
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
        CLASS(environ_boundary_electronic), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: routine = 'destroy_environ_boundary'
        !
        !--------------------------------------------------------------------------------
        !
        CALL this%pre_destroy()
        !
        CALL this%density%destroy()
        !
        CALL this%dscaled%destroy()
        !
        CALL this%d2scaled%destroy()
        !
        IF (ASSOCIATED(this%ions)) NULLIFY (this%ions)
        !
        NULLIFY (this%electrons)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE destroy_environ_boundary
    !------------------------------------------------------------------------------------
    !>
    !! Calculates the dielectric constant as a function of the charge
    !! density, and the derivatives of the the dielectric constant
    !! with respect to the charge density. Additionally calculates the
    !! volume and surface components.
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE boundary_of_density(this, density)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_density), OPTIONAL, TARGET, INTENT(IN) :: density
        !
        CLASS(environ_boundary_electronic), TARGET, INTENT(INOUT) :: this
        !
        REAL(DP), POINTER :: rhomax, rhomin, fact
        !
        TYPE(environ_density), POINTER :: denloc
        TYPE(environ_gradient), POINTER :: gradloc
        TYPE(environ_hessian), POINTER :: hessloc
        !
        INTEGER :: i, j
        !
        CHARACTER(LEN=80) :: routine = 'boundary_of_density'
        !
        !--------------------------------------------------------------------------------
        !
        IF (PRESENT(density)) THEN
            denloc => density
        ELSE
            denloc => this%density
        END IF
        !
        IF (.NOT. ASSOCIATED(denloc%cell, this%scaled%cell)) &
            CALL io%error(routine, "Inconsistent domains", 1)
        !
        rhomax => this%rhomax
        rhomin => this%rhomin
        fact => this%fact
        !
        !--------------------------------------------------------------------------------
        !
        ASSOCIATE (cell => denloc%cell, &
                   derivatives => this%cores%derivatives, &
                   rho => denloc%of_r, &
                   scal => this%scaled, &
                   dscal => this%dscaled, &
                   d2scal => this%d2scaled, &
                   grad => this%gradient, &
                   lapl => this%laplacian, &
                   dsurf => this%dsurface)
            !
            !----------------------------------------------------------------------------
            !
            DO i = 1, cell%ir_end
                scal%of_r(i) = 1.D0 - sfunct1(rho(i), rhomax, rhomin, fact)
                dscal%of_r(i) = -dsfunct1(rho(i), rhomax, rhomin, fact)
                d2scal%of_r(i) = -d2sfunct1(rho(i), rhomax, rhomin, fact)
            END DO
            !
            !----------------------------------------------------------------------------
            ! Compute boundary derivatives, if needed
            !
            IF (this%need_hessian) THEN
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

            call io%writer(this%derivatives_method)

            SELECT CASE (this%derivatives_method)
                !
            CASE ('fft')
                CALL this%compute_boundary_derivatives_fft(scal, hessloc)
                !
            CASE ('chain')
                CALL this%compute_boundary_derivatives_fft(denloc, hessloc)
                !
                !------------------------------------------------------------------------
                ! Apply chain rule
                !
                IF (this%need_hessian) THEN
                    !
                    IF (this%solvent_aware) THEN
                        !
                        DO i = 1, 3
                            !
                            DO j = 1, 3
                                !
                                hessloc%of_r(i, j, :) = &
                                    hessloc%of_r(i, j, :) * dscal%of_r + &
                                    grad%of_r(i, :) * grad%of_r(j, :) * d2scal%of_r
                                !
                            END DO
                            !
                        END DO
                        !
                    END IF
                    !
                END IF
                !
                IF (this%need_laplacian) &
                    lapl%of_r = lapl%of_r * dscal%of_r + &
                                (grad%of_r(1, :)**2 + &
                                 grad%of_r(2, :)**2 + &
                                 grad%of_r(3, :)**2) * d2scal%of_r
                !
                IF (this%need_gradient) THEN
                    !
                    DO i = 1, 3
                        grad%of_r(i, :) = grad%of_r(i, :) * dscal%of_r
                    END DO
                    !
                END IF
                !
            CASE ('chain-local')
                CALL io%writer('switched boundary derivative to chain-local')
                ! Save local gradient of electronic density
                ALLOCATE (gradloc)
                CALL gradloc%init(cell)
                gradloc%of_r = grad%of_r
                !
                CALL this%compute_boundary_derivatives_fft(denloc, hessloc)
                grad%of_r = gradloc%of_r
                DEALLOCATE(gradloc)
                !
                !------------------------------------------------------------------------
                ! Apply chain rule
                !
                IF (this%need_hessian) THEN
                    !
                    IF (this%solvent_aware) THEN
                        !
                        DO i = 1, 3
                            !
                            DO j = 1, 3
                                !
                                hessloc%of_r(i, j, :) = &
                                    hessloc%of_r(i, j, :) * dscal%of_r + &
                                    grad%of_r(i, :) * grad%of_r(j, :) * d2scal%of_r
                                !
                            END DO
                            !
                        END DO
                        !
                    END IF
                    !
                END IF
                !
                IF (this%need_laplacian) &
                    lapl%of_r = lapl%of_r * dscal%of_r + &
                                (grad%of_r(1, :)**2 + &
                                 grad%of_r(2, :)**2 + &
                                 grad%of_r(3, :)**2) * d2scal%of_r
                !
                IF (this%need_gradient) THEN
                    !
                    DO i = 1, 3
                        grad%of_r(i, :) = grad%of_r(i, :) * dscal%of_r
                    END DO
                    !
                END IF
                !
            CASE DEFAULT
                CALL io%error(routine, "Unexpected derivatives method", 1)
                !
            END SELECT
            !
            !----------------------------------------------------------------------------
            ! Final updates
            !
            this%volume = scal%integrate()
            !
            IF (this%need_gradient) THEN
                !
                CALL grad%update_modulus()
                !
                this%surface = grad%modulus%integrate()
            END IF
            !
            IF (this%need_hessian .AND. .NOT. this%solvent_aware) CALL hessloc%destroy()
            !
        END ASSOCIATE
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE boundary_of_density
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
        CLASS(environ_boundary_electronic), TARGET, INTENT(IN) :: this
        INTEGER, INTENT(IN) :: index
        !
        TYPE(environ_gradient), INTENT(INOUT) :: partial
        !
        INTEGER :: i
        REAL(DP) :: spurious_force
        !
        CHARACTER(LEN=80) :: routine = 'calc_dboundary_dions'
        !
        !--------------------------------------------------------------------------------
        !
        IF (this%mode == 'electronic') RETURN
        ! exit if boundary is only defined on electronic density
        !
        IF (index > this%ions%number) &
            CALL io%error(routine, "Index greater than number of ions", 1)
        !
        IF (index <= 0) CALL io%error(routine, "Index of ion is zero or lower", 1)
        !
        IF (this%ions%core_electrons%number == 0) &
            CALL io%error(routine, "Missing details of core electrons", 1)
        !
        IF (.NOT. ASSOCIATED(this%dscaled%cell, partial%cell)) &
            CALL io%error(routine, "Mismatch or unassociated boundary derivative", 1)
        !
        CALL this%ions%core_electrons%array(index)%gradient(partial, .TRUE.)
        !
        DO i = 1, 3
            partial%of_r(i, :) = -partial%of_r(i, :) * this%dscaled%of_r
        END DO
        !
        CALL partial%update_modulus()
        !
        spurious_force = partial%modulus%integrate()
        !
        IF (io%lnode .AND. spurious_force > 1.D-5) &
            WRITE (io%unit, 1000) index, spurious_force
        !
        !--------------------------------------------------------------------------------
        !
1000    FORMAT(" WARNING: Unphysical forces due to core electrons are non-negligible ", /, &
               " atom type ", I3, " is subject to a spurious force of ", F12.6)
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
        CLASS(environ_boundary_electronic), INTENT(IN) :: this
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
        IF (io%lnode) THEN
            WRITE (local_unit, 1100) this%rhomax, this%rhomin
            !
            IF (local_verbose >= 3) WRITE (local_unit, 1101) this%fact
            !
        END IF
        !
        IF (local_verbose >= 4) THEN
            !
            CALL this%density%printout(passed_verbose, debug_verbose, local_unit)
            !
            IF (io%lnode .AND. ASSOCIATED(this%ions)) WRITE (local_unit, 1102)
            !
        END IF
        !
        IF (local_verbose >= 5) THEN
            !
            CALL this%dscaled%printout(passed_verbose, debug_verbose, local_unit)
            !
            CALL this%d2scaled%printout(passed_verbose, debug_verbose, local_unit)
            !
        END IF
        !
        FLUSH (local_unit)
        !
        !--------------------------------------------------------------------------------
        !
1100    FORMAT(/, " using the optimal SCCS function:", /, &
                " rhomax                     = ", F14.7, /, &
                " rhomin                     = ", F14.7)
        !
1101    FORMAT(" log(rhomax/rhomin)         = ", F14.7)
        !
1102    FORMAT(/, " adding fictitious core-electrons")
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE print_environ_boundary
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE class_boundary_electronic
!----------------------------------------------------------------------------------------
