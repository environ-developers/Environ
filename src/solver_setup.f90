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
! Authors: Francesco Nattino  (THEOS and NCCR-MARVEL, EPFL)
!          Oliviero Andreussi (Department of Physics, UNT)
!          Nicola Marzari     (THEOS and NCCR-MARVEL, EPFL)
!          Edan Bainglass     (Department of Physics, UNT)
!
!----------------------------------------------------------------------------------------
!>
!!
!----------------------------------------------------------------------------------------
MODULE class_solver_setup
    !------------------------------------------------------------------------------------
    !
    USE class_io, ONLY: io
    !
    USE environ_param, ONLY: DP, e2, tpi
    !
    USE class_cell
    USE class_density
    !
    USE class_core_fft_electrostatics
    USE class_core_1da
    !
    USE class_solver
    USE class_solver_direct
    USE class_solver_fixedpoint
    USE class_solver_gradient
    USE class_solver_iterative
    USE class_solver_newton
    !
    USE class_charges
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
    TYPE, PUBLIC :: electrostatic_setup
        !--------------------------------------------------------------------------------
        !
        CHARACTER(LEN=80) :: problem
        !
        CLASS(electrostatic_solver), POINTER :: solver => NULL()
        !
        CLASS(electrostatic_setup), POINTER :: inner => NULL()
        !
        !--------------------------------------------------------------------------------
    CONTAINS
        !--------------------------------------------------------------------------------
        !
        PROCEDURE, PRIVATE :: create => create_electrostatic_setup
        PROCEDURE :: init => init_electrostatic_setup
        PROCEDURE :: set_flags => set_electrostatic_flags
        PROCEDURE :: destroy => destroy_electrostatic_setup
        !
        PROCEDURE :: calc_v => calc_velectrostatic
        PROCEDURE :: calc_e => calc_eelectrostatic
        PROCEDURE :: calc_f => calc_felectrostatic
        !
        !--------------------------------------------------------------------------------
    END TYPE electrostatic_setup
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
    SUBROUTINE create_electrostatic_setup(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(electrostatic_setup), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: sub_name = 'create_electrostatic_setup'
        !
        !--------------------------------------------------------------------------------
        !
        IF (ASSOCIATED(this%solver)) CALL io%create_error(sub_name)
        !
        IF (ASSOCIATED(this%inner)) CALL io%create_error(sub_name)
        !
        !--------------------------------------------------------------------------------
        !
        this%problem = 'poisson'
        !
        NULLIFY (this%solver)
        NULLIFY (this%inner)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE create_electrostatic_setup
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_electrostatic_setup(this, problem, solver)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CHARACTER(LEN=80), INTENT(IN) :: problem
        CLASS(electrostatic_solver), TARGET, INTENT(IN) :: solver
        !
        CLASS(electrostatic_setup), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: sub_name = 'init_electrostatic_setup'
        !
        !--------------------------------------------------------------------------------
        ! Sanity check on the global setup
        !
        CALL this%create()
        !
        this%problem = problem
        this%solver => solver
        !
        SELECT CASE (TRIM(ADJUSTL(this%problem)))
            !
        CASE ('poisson')
            !
        CASE ('generalized', 'gpe')
            !
            SELECT TYPE (solver)
                !
            TYPE IS (solver_direct)
                !
                CALL io%error(sub_name, &
                              'Cannot use a direct solver for &
                              &the Generalized Poisson eq.', 1)
            END SELECT
            !
        CASE ('linpb', 'linmodpb', 'linearized-pb')
            !
            SELECT TYPE (solver)
                !
            TYPE IS (solver_direct)
                !
                CALL io%error(sub_name, &
                              'Only gradient-based solver for &
                              &the linearized Poisson-Boltzmann eq.', 1)
                !
            TYPE IS (solver_fixedpoint)
                !
                CALL io%error(sub_name, &
                              'Only gradient-based solver for &
                              &the linearized Poisson-Boltzmann eq.', 1)
                !
            END SELECT
            !
            IF (ASSOCIATED(solver%cores%correction) .AND. &
                solver%cores%correction%type_ /= '1da') THEN
                !
                CALL io%error(sub_name, &
                              'Linearized-PB problem requires &
                              &parabolic pbc correction.', 1)
                !
            END IF
            !
        CASE ('pb', 'modpb', 'poisson-boltzmann')
            !
            SELECT TYPE (solver)
                !
            TYPE IS (solver_direct)
                !
                CALL io%error(sub_name, &
                              'No direct or gradient-based solver for &
                              &the full Poisson-Boltzmann eq.', 1)
                !
            TYPE IS (solver_gradient)
                !
                CALL io%error(sub_name, &
                              'No direct or gradient-based solver for &
                              &the full Poisson-Boltzmann eq.', 1)
                !
            END SELECT
            !
            IF (ASSOCIATED(solver%cores%correction) .AND. &
                solver%cores%correction%type_ /= '1da') THEN
                !
                CALL io%error(sub_name, &
                              'Full-PB problem requires parabolic pbc correction.', 1)
                !
            END IF
            !
        CASE DEFAULT
            CALL io%error(sub_name, 'Unexpected keyword for electrostatic problem', 1)
            !
        END SELECT
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_electrostatic_setup
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE set_electrostatic_flags(this, need_auxiliary, need_gradient, &
                                       need_factsqrt)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(electrostatic_setup), INTENT(INOUT) :: this
        !
        LOGICAL, INTENT(INOUT) :: need_auxiliary, need_gradient, need_factsqrt
        !
        !--------------------------------------------------------------------------------
        !
        SELECT CASE (this%problem)
            !
        CASE ('generalized', 'linpb', 'linmodpb', 'pb', 'modpb')
            !
            SELECT TYPE (solver => this%solver)
                !
            TYPE IS (solver_gradient)
                !
                SELECT CASE (solver%preconditioner)
                    !
                CASE ('sqrt')
                    need_factsqrt = .TRUE.
                    !
                CASE ('left', 'none')
                    need_gradient = .TRUE.
                    !
                END SELECT
                !
            END SELECT
            !
            SELECT TYPE (solver => this%solver)
                !
            CLASS IS (solver_iterative)
                !
                IF (solver%auxiliary /= 'none') need_auxiliary = .TRUE.
                !
            END SELECT
            !
        END SELECT
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE set_electrostatic_flags
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    RECURSIVE SUBROUTINE destroy_electrostatic_setup(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(electrostatic_setup), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: sub_name = 'destroy_electrostatic_setup'
        !
        !--------------------------------------------------------------------------------
        !
        IF (.NOT. ASSOCIATED(this%solver)) CALL io%destroy_error(sub_name)
        !
        !--------------------------------------------------------------------------------
        !
        SELECT TYPE (solver => this%solver)
            !
        CLASS IS (solver_direct)
            CALL solver%destroy()
            !
        END SELECT
        !
        NULLIFY (this%solver)
        !
        IF (ASSOCIATED(this%inner)) THEN
            !
            CALL this%inner%destroy()
            !
            NULLIFY (this%inner)
        END IF
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE destroy_electrostatic_setup
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !                                  GENERAL METHODS
    !
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !>
    !! Calculates the electrostatic embedding contribution to the Kohn-Sham potential
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_velectrostatic(this, charges, potential)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(electrostatic_setup), INTENT(INOUT) :: this
        TYPE(environ_charges), INTENT(INOUT) :: charges
        TYPE(environ_density), INTENT(INOUT) :: potential
        !
        CHARACTER(LEN=80) :: sub_name = 'calc_velectrostatic'
        !
        !--------------------------------------------------------------------------------
        !
        CALL env_start_clock(sub_name)
        !
        potential%of_r = 0.D0
        !
        !--------------------------------------------------------------------------------
        ! Select the appropriate combination of problem and solver
        !
        SELECT CASE (this%problem)
            !
        CASE ('poisson')
            !
            SELECT TYPE (solver => this%solver)
                !
            TYPE IS (solver_direct)
                CALL solver%poisson(charges, potential)
                !
            CLASS DEFAULT
                CALL io%error(sub_name, 'Unexpected solver', 1)
                !
            END SELECT
            !
        CASE ('generalized')
            !
            IF (.NOT. ASSOCIATED(charges%dielectric)) &
                CALL io%error(sub_name, 'Missing details of dielectric medium', 1)
            !
            SELECT TYPE (solver => this%solver)
                !
            TYPE IS (solver_direct)
                !
                CALL io%error(sub_name, 'Option not yet implemented', 1)
                !
                ! CALL generalized_direct() #TODO future work
                !
            TYPE IS (solver_gradient)
                CALL solver%generalized(charges, potential)
                !
            TYPE IS (solver_fixedpoint)
                CALL solver%generalized(charges, potential)
                !
            CLASS DEFAULT
                CALL io%error(sub_name, 'Unexpected solver', 1)
                !
            END SELECT
            !
        CASE ('linpb', 'linmodpb')
            !
            IF (.NOT. ASSOCIATED(charges%electrolyte)) &
                CALL io%error(sub_name, 'Missing details of electrolyte ions', 1)
            !
            SELECT TYPE (solver => this%solver)
                !
            TYPE IS (solver_direct)
                !
                CALL io%error(sub_name, 'Option not yet implemented', 1)
                !
                ! CALL linpb_direct() #TODO future work
                !
            TYPE IS (solver_gradient)
                CALL solver%linearized_pb(charges, potential)
                !
            CLASS DEFAULT
                CALL io%error(sub_name, 'Unexpected solver', 1)
                !
            END SELECT
            !
        CASE ('pb', 'modpb')
            !
            IF (.NOT. ASSOCIATED(charges%electrolyte)) &
                CALL io%error(sub_name, 'Missing details of electrolyte ions', 1)
            !
            SELECT TYPE (solver => this%solver)
                !
            TYPE IS (solver_direct)
                !
                CALL io%error(sub_name, 'Option not yet implemented', 1)
                !
            TYPE IS (solver_fixedpoint)
                !
                IF (ASSOCIATED(this%inner)) THEN
                    !
                    IF (.NOT. ASSOCIATED(charges%dielectric)) &
                        CALL io%error(sub_name, &
                                      'Missing details of dielectric medium', 1)
                    !
                    CALL solver%pb_nested(charges, potential, this%inner%solver)
                    !
                ELSE
                    CALL solver%pb_nested(charges, potential)
                END IF
                !
            TYPE IS (solver_newton)
                !
                IF (.NOT. ASSOCIATED(this%inner)) &
                    CALL io%error(sub_name, &
                                  'Missing details of inner electrostatic setup', 1)
                !
                CALL solver%pb_nested(this%inner%solver, charges, potential)
                !
            CLASS DEFAULT
                CALL io%error(sub_name, 'Unexpected solver', 1)
                !
            END SELECT
            !
        CASE DEFAULT
            CALL io%error(sub_name, 'Unexpected problem keyword', 1)
            !
        END SELECT
        !
        CALL env_stop_clock(sub_name)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_velectrostatic
    !------------------------------------------------------------------------------------
    !>
    !! Calculates the electrostatic embedding contribution to the energy
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_eelectrostatic(this, charges, potential, energy)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_density), INTENT(IN) :: potential
        !
        CLASS(electrostatic_setup), INTENT(INOUT) :: this
        TYPE(environ_charges), INTENT(INOUT) :: charges
        REAL(DP), INTENT(OUT) :: energy
        !
        REAL(DP) :: degauss, eself
        !
        CHARACTER(LEN=80) :: sub_name = 'calc_eelectrostatic'
        !
        !--------------------------------------------------------------------------------
        !
        CALL env_start_clock(sub_name)
        !
        IF (.NOT. ASSOCIATED(charges%density%cell, potential%cell)) &
            CALL io%error(sub_name, 'Mismatch in charges and potential domains', 1)
        !
        energy = 0.D0
        eself = 0.D0
        degauss = 0.D0
        !
        !--------------------------------------------------------------------------------
        ! Electrons and nuclei and external charges
        !
        energy = energy + 0.5D0 * charges%density%scalar_product(potential)
        degauss = degauss + charges%charge
        !
        !--------------------------------------------------------------------------------
        ! Include environment contributions
        !
        IF (charges%include_dielectric) &
            degauss = degauss + charges%dielectric%charge * 0.5D0 ! polarization charge
        !
        IF (charges%include_electrolyte) THEN
            !
            ! note: electrolyte electrostatic interaction should be negative
            energy = energy - 0.5D0 * &
                     charges%electrolyte%density%scalar_product(potential)
            !
            degauss = degauss + charges%electrolyte%charge ! electrolyte charge
        END IF
        !
        !--------------------------------------------------------------------------------
        ! Adding correction for point-like nuclei: only affects simulations of charged
        ! systems, it does not affect forces, but shift the energy depending on the
        ! fictitious Gaussian spread of the nuclei
        !
        ! Compute spurious self-polarization energy
        !
        eself = charges%ions%selfenergy_correction * e2
        !
        SELECT TYPE (core => this%solver%cores%core)
            !
        TYPE IS (core_fft_electrostatics)
            !
            IF (core%use_internal_pbc_corr .OR. &
                ASSOCIATED(this%solver%cores%correction)) THEN
                degauss = 0.D0
            ELSE
                !
                degauss = -degauss * charges%ions%quadrupole_correction * &
                          e2 * tpi / charges%density%cell%omega
                !
            END IF
            !
        END SELECT
        !
        energy = energy + eself + degauss
        !
        CALL env_stop_clock(sub_name)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_eelectrostatic
    !------------------------------------------------------------------------------------
    !>
    !! Calculates the electrostatic contribution to the interatomic forces
    !!
    !! If called by the reference solver, aux will only include ions and electrons.
    !! If called by the outer solver, aux will include the full charge density,
    !! including any environment-specific contributions
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_felectrostatic(this, natoms, charges, force, ldoublecell)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: natoms
        TYPE(environ_charges), INTENT(IN) :: charges
        LOGICAL, INTENT(IN) :: ldoublecell
        !
        CLASS(electrostatic_setup), INTENT(INOUT) :: this
        REAL(DP), INTENT(INOUT) :: force(3, natoms)
        !
        TYPE(environ_density) :: aux
        !
        CHARACTER(LEN=80) :: sub_name = 'calc_felectrostatic'
        !
        !--------------------------------------------------------------------------------
        !
        CALL env_start_clock(sub_name)
        !
        !--------------------------------------------------------------------------------
        ! Electrostatic contributions to forces are formulated in terms
        ! of the auxiliary charges times the derivatives of the pseudopotentials
        !
        CALL aux%init(charges%density%cell)
        !
        IF (ldoublecell) THEN
            aux%of_r = charges%density%of_r ! ions, electrons, and (if active) externals
            !
        ELSE IF (charges%include_externals) THEN
            !
            IF (.NOT. ASSOCIATED(charges%externals)) &
                CALL io%error(sub_name, 'Missing expected charge component', 1)
            !
            aux%of_r = aux%of_r + charges%externals%density%of_r
        END IF
        !
        IF (charges%include_dielectric) THEN
            !
            IF (.NOT. ASSOCIATED(charges%dielectric)) &
                CALL io%error(sub_name, 'Missing expected charge component', 1)
            !
            aux%of_r = aux%of_r + charges%dielectric%density%of_r
        END IF
        !
        IF (charges%include_electrolyte) THEN
            !
            IF (.NOT. ASSOCIATED(charges%electrolyte)) &
                CALL io%error(sub_name, 'Missing expected charge component', 1)
            !
            aux%of_r = aux%of_r + charges%electrolyte%density%of_r
        END IF
        !
        SELECT TYPE (core => this%solver%cores%core)
            !
        TYPE IS (core_fft_electrostatics)
            CALL this%solver%cores%force(natoms, aux, charges%ions, force)
            !
        END SELECT
        !
        IF (ASSOCIATED(this%solver%cores%correction)) THEN
            !
            IF (.NOT. ldoublecell) aux%of_r = aux%of_r + charges%density%of_r
            !
            CALL this%solver%cores%correction%calc_f(natoms, charges, aux, force)
            !
        END IF
        !
        CALL aux%destroy()
        !
        CALL env_stop_clock(sub_name)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_felectrostatic
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE class_solver_setup
!----------------------------------------------------------------------------------------
