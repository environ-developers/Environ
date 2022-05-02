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
    USE class_core_fft
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
        CHARACTER(LEN=80) :: routine = 'create_electrostatic_setup'
        !
        !--------------------------------------------------------------------------------
        !
        IF (ASSOCIATED(this%solver)) CALL io%create_error(routine)
        !
        IF (ASSOCIATED(this%inner)) CALL io%create_error(routine)
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
        CHARACTER(LEN=*), INTENT(IN) :: problem
        CLASS(electrostatic_solver), TARGET, INTENT(IN) :: solver
        !
        CLASS(electrostatic_setup), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: routine = 'init_electrostatic_setup'
        !
        !--------------------------------------------------------------------------------
        !
        CALL this%create()
        !
        this%problem = problem
        this%solver => solver
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_electrostatic_setup
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
        CHARACTER(LEN=80) :: routine = 'destroy_electrostatic_setup'
        !
        !--------------------------------------------------------------------------------
        !
        IF (.NOT. ASSOCIATED(this%solver)) CALL io%destroy_error(routine)
        !
        !--------------------------------------------------------------------------------
        !
        CALL this%solver%destroy()
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
    SUBROUTINE calc_velectrostatic(this, charges, v)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(electrostatic_setup), INTENT(IN) :: this
        !
        TYPE(environ_density), INTENT(INOUT) :: v
        TYPE(environ_charges), INTENT(INOUT) :: charges
        !
        CHARACTER(LEN=80) :: routine = 'calc_velectrostatic'
        !
        !--------------------------------------------------------------------------------
        !
        CALL env_start_clock(routine)
        !
        v%of_r = 0.D0
        !
        !--------------------------------------------------------------------------------
        ! Select the appropriate combination of problem and solver
        !
        SELECT CASE (this%problem)
            !
        CASE ('poisson')
            CALL this%solver%poisson(charges, v)
            !
        CASE ('generalized')
            !
            IF (.NOT. ASSOCIATED(charges%dielectric)) &
                CALL io%error(routine, "Missing details of dielectric medium", 1)
            !
            CALL this%solver%generalized(charges, v)
            !
        CASE ('linpb', 'linmodpb')
            !
            IF (.NOT. ASSOCIATED(charges%electrolyte)) &
                CALL io%error(routine, "Missing details of electrolyte ions", 1)
            !
            CALL this%solver%linearized_pb(charges, v)
            !
        CASE ('pb', 'modpb')
            !
            IF (.NOT. ASSOCIATED(charges%electrolyte)) &
                CALL io%error(routine, "Missing details of electrolyte ions", 1)
            !
            IF (ASSOCIATED(this%inner)) THEN
                !
                IF (.NOT. ASSOCIATED(charges%dielectric)) &
                    CALL io%error(routine, "Missing details of dielectric medium", 1)
                !
                CALL this%solver%pb_nested(charges, v, this%inner%solver)
                !
            ELSE
                CALL this%solver%pb_nested(charges, v)
            END IF
            !
        CASE DEFAULT
            CALL io%error(routine, "Unexpected 'problem'", 1)
            !
        END SELECT
        !
        CALL env_stop_clock(routine)
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
        CHARACTER(LEN=80) :: routine = 'calc_eelectrostatic'
        !
        !--------------------------------------------------------------------------------
        !
        CALL env_start_clock(routine)
        !
        IF (.NOT. ASSOCIATED(charges%density%cell, potential%cell)) &
            CALL io%error(routine, "Mismatch in charges and potential domains", 1)
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
        IF (this%solver%cores%internal_correction .OR. &
            this%solver%cores%has_corrections) THEN
            degauss = 0.D0
        ELSE
            !
            degauss = -degauss * charges%ions%quadrupole_correction * &
                      e2 * tpi / charges%density%cell%omega
            !
        END IF
        !
        energy = energy + eself + degauss
        !
        CALL env_stop_clock(routine)
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
    SUBROUTINE calc_felectrostatic(this, nat, charges, force, ldoublecell)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: nat
        TYPE(environ_charges), INTENT(IN) :: charges
        LOGICAL, INTENT(IN) :: ldoublecell
        !
        CLASS(electrostatic_setup), INTENT(INOUT) :: this
        REAL(DP), INTENT(INOUT) :: force(3, nat)
        !
        TYPE(environ_density) :: aux
        !
        CHARACTER(LEN=80) :: routine = 'calc_felectrostatic'
        !
        !--------------------------------------------------------------------------------
        !
        CALL env_start_clock(routine)
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
                CALL io%error(routine, "Missing expected charge component", 1)
            !
            aux%of_r = aux%of_r + charges%externals%density%of_r
        END IF
        !
        IF (charges%include_dielectric) THEN
            !
            IF (.NOT. ASSOCIATED(charges%dielectric)) &
                CALL io%error(routine, "Missing expected charge component", 1)
            !
            aux%of_r = aux%of_r + charges%dielectric%density%of_r
        END IF
        !
        IF (charges%include_electrolyte) THEN
            !
            IF (.NOT. ASSOCIATED(charges%electrolyte)) &
                CALL io%error(routine, "Missing expected charge component", 1)
            !
            aux%of_r = aux%of_r + charges%electrolyte%density%of_r
        END IF
        !
        ASSOCIATE (cores => this%solver%cores)
            !
            CALL cores%electrostatics%force(nat, aux, charges%ions%smeared_ions, force)
            !
            IF (cores%has_corrections) THEN
                !
                IF (.NOT. ldoublecell) aux%of_r = aux%of_r + charges%density%of_r
                !
                CALL cores%corrections%force_periodic(nat, charges%ions%smeared_ions, &
                                                      aux, force)
                !
            END IF
            !
        END ASSOCIATE
        !
        CALL aux%destroy()
        !
        CALL env_stop_clock(routine)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_felectrostatic
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE class_solver_setup
!----------------------------------------------------------------------------------------
