!----------------------------------------------------------------------------------------
!
! Copyright (C) 2021 ENVIRON (www.quantum-environ.org)
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
!          Ismaila Dabo       (DMSE, Penn State)
!          Nicola Marzari     (THEOS and NCCR-MARVEL, EPFL)
!          Edan Bainglass     (Department of Physics, UNT)
!
!----------------------------------------------------------------------------------------
!>
!!
!----------------------------------------------------------------------------------------
MODULE class_environ
    !------------------------------------------------------------------------------------
    !
    USE class_io, ONLY: io
    !
    USE environ_param, ONLY: DP, RYTOEV
    !
    USE class_setup
    !
    USE env_base_input
    !
    USE class_cell
    USE class_density
    USE class_gradient
    !
    USE class_boundary
    USE class_boundary_electronic
    USE class_boundary_ionic
    USE class_boundary_system
    USE class_charges
    USE class_dielectric
    USE class_electrolyte
    USE class_electrons
    USE class_externals
    USE class_ions
    USE class_semiconductor
    USE class_system
    !
    USE env_write_cube
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
    TYPE, PUBLIC :: environ_main
        !--------------------------------------------------------------------------------
        !
        LOGICAL :: initialized = .FALSE.
        !
        !--------------------------------------------------------------------------------
        !
        TYPE(environ_setup), POINTER :: setup => NULL()
        !
        !--------------------------------------------------------------------------------
        ! Mapped quantities
        !
        TYPE(environ_ions) :: system_ions, environment_ions
        TYPE(environ_system) :: system_system, environment_system
        TYPE(environ_electrons) :: system_electrons, environment_electrons
        TYPE(environ_charges) :: system_charges, environment_charges
        !
        !--------------------------------------------------------------------------------
        ! Details of the continuum interface
        !
        CLASS(environ_boundary), ALLOCATABLE :: solvent
        !
        !--------------------------------------------------------------------------------
        ! Response properties
        !
        TYPE(environ_electrons) :: system_response_electrons, &
                                   environment_response_electrons
        !
        TYPE(environ_charges) :: system_response_charges, &
                                 environment_response_charges
        !
        !--------------------------------------------------------------------------------
        ! Physical objects
        !
        TYPE(environ_externals) :: externals
        TYPE(environ_dielectric) :: static, optical
        TYPE(environ_electrolyte) :: electrolyte
        TYPE(environ_semiconductor) :: semiconductor
        !
        !--------------------------------------------------------------------------------
        ! Computed physical variables
        !
        TYPE(environ_density) :: &
            vzero, velectrostatic, vreference, dvtot, vconfine, vsoftcavity
        !
        REAL(DP) :: evolume = 0.0_DP
        REAL(DP) :: esurface = 0.0_DP
        REAL(DP) :: econfine = 0.0_DP
        REAL(DP) :: deenviron = 0.0_DP
        REAL(DP) :: eelectrolyte = 0.0_DP
        REAL(DP) :: eelectrostatic = 0.0_DP
        !
        !--------------------------------------------------------------------------------
        ! MBPOL charge density
        !
        TYPE(environ_density) :: additional_charges
        !
        !--------------------------------------------------------------------------------
    CONTAINS
        !--------------------------------------------------------------------------------
        !
        PROCEDURE, PRIVATE :: create => create_environ_base
        PROCEDURE :: init => init_environ_base
        PROCEDURE :: add_charges => environ_add_charges
        !
        PROCEDURE :: get_vzero
        PROCEDURE :: get_dvtot
        !
        PROCEDURE :: update_electrons => environ_update_electrons
        PROCEDURE :: update_ions => environ_update_ions
        PROCEDURE :: update_potential => environ_update_potential
        PROCEDURE :: update_response => environ_update_response
        PROCEDURE :: update_cell_dependent_quantities
        !
        PROCEDURE, PRIVATE :: init_physical => environ_init_physical
        PROCEDURE, PRIVATE :: init_potential => environ_init_potential
        !
        PROCEDURE :: print_energies => print_environ_energies
        PROCEDURE :: print_potential_shift => print_environ_potential_shift
        !
        !--------------------------------------------------------------------------------
    END TYPE environ_main
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
    SUBROUTINE create_environ_base(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_main), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: routine = 'create_environ_base'
        !
        !--------------------------------------------------------------------------------
        !
        this%initialized = .FALSE.
        this%evolume = 0.0_DP
        this%esurface = 0.0_DP
        this%econfine = 0.0_DP
        this%deenviron = 0.0_DP
        this%eelectrolyte = 0.0_DP
        this%eelectrostatic = 0.0_DP
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE create_environ_base
    !------------------------------------------------------------------------------------
    !>
    !! Subroutine to initialize fundamental quantities needed by the
    !! environ modules. This subroutine is called by init_run.f90, thus
    !! only once per pw.x execution.
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_environ_base(this, nat, ntyp, ityp, zv, label, number, weight)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: nat, ntyp
        INTEGER, INTENT(IN) :: ityp(nat)
        REAL(DP), INTENT(IN) :: zv(ntyp)
        CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: label(ntyp)
        INTEGER, OPTIONAL, INTENT(IN) :: number(ntyp)
        REAL(DP), OPTIONAL, INTENT(IN) :: weight(ntyp)
        !
        CLASS(environ_main), INTENT(INOUT) :: this
        !
        !--------------------------------------------------------------------------------
        !
        CALL this%create()
        !
        CALL this%init_potential()
        !
        CALL this%init_physical(nat, ntyp, ityp, zv, label, number, weight)
        !
        this%initialized = .TRUE.
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_environ_base
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE environ_add_charges(this, nnr, density, label)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: nnr
        REAL(DP), INTENT(IN) :: density(nnr)
        CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: label
        !
        CLASS(environ_main), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: local_label = 'additional_charges'
        !
        !--------------------------------------------------------------------------------
        !
        IF (PRESENT(label)) local_label = label
        !
        CALL this%additional_charges%init(this%setup%system_cell, local_label)
        !
        this%additional_charges%of_r = density
        !
        CALL this%system_charges%add(additional_charges=this%additional_charges)
        !
        this%setup%laddcharges = .TRUE.
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE environ_add_charges
    !------------------------------------------------------------------------------------
    !>
    !! Save local potential that will be overwritten by environ
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE environ_update_potential(this, nnr, vltot)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: nnr
        REAL(DP), INTENT(IN) :: vltot(nnr)
        !
        CLASS(environ_main), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: routine = 'environ_update_potential'
        !
        !--------------------------------------------------------------------------------
        !
        this%vzero%lupdate = .TRUE.
        !
        IF (.NOT. ASSOCIATED(this%vzero%cell)) RETURN
        !
        IF (this%vzero%cell%nnr /= nnr) &
            CALL io%error(routine, "Inconsistent size in input potential", 1)
        !
        this%vzero%of_r = vltot
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE environ_update_potential
    !------------------------------------------------------------------------------------
    !>
    !! Update fixed quantities defined inside the cell
    !! NOTE: updating depends on cell%lupdate
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE update_cell_dependent_quantities(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_main), TARGET, INTENT(INOUT) :: this
        !
        TYPE(environ_setup), POINTER :: setup
        !
        CHARACTER(LEN=80) :: routine = 'update_cell_dependent_quantities'
        !
        !--------------------------------------------------------------------------------
        !
        setup => this%setup
        !
        IF (setup%lstatic) CALL this%static%update()
        !
        IF (setup%loptical) CALL this%optical%update()
        !
        IF (setup%lelectrolyte) CALL this%electrolyte%update()
        !
        IF (setup%lexternals) CALL this%externals%update()
        !
        IF (setup%lsemiconductor) CALL this%semiconductor%update()
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE update_cell_dependent_quantities
    !------------------------------------------------------------------------------------
    !>
    !! Initialize the ions-related quantities to be used in the Environ
    !! modules. This initialization is called by electrons.f90, thus it
    !! is performed at every step of electronic optimization. It may not
    !! be the most efficient choice, but it is a safe choice.
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE environ_update_ions(this, nat, tau, com)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: nat
        REAL(DP), INTENT(IN) :: tau(3, nat)
        REAL(DP), OPTIONAL, INTENT(IN) :: com(3)
        !
        CLASS(environ_main), TARGET, INTENT(INOUT) :: this
        !
        TYPE(environ_setup), POINTER :: setup
        !
        !--------------------------------------------------------------------------------
        !
        setup => this%setup
        !
        this%system_ions%lupdate = .TRUE.
        this%environment_ions%lupdate = .TRUE.
        !
        setup%niter_ionic = setup%niter_ionic + 1
        !
        !--------------------------------------------------------------------------------
        ! Update system ions parameters
        !
        CALL this%system_ions%update(nat, tau, com)
        !
        !--------------------------------------------------------------------------------
        ! Update system system parameters
        !
        this%system_system%lupdate = .TRUE.
        !
        IF (environ_debug) THEN
            !
            CALL this%system_system%update(system_pos)
            ! using fixed system_pos from input for debugging with finite-differences
            !
        ELSE
            CALL this%system_system%update(com)
        END IF
        !
        !--------------------------------------------------------------------------------
        ! Update cell mapping
        !
        CALL setup%update_mapping(this%system_system%com)
        !
        !--------------------------------------------------------------------------------
        ! Update environment ions parameters
        !
        CALL this%environment_ions%update(nat, tau, com)
        !
        !--------------------------------------------------------------------------------
        ! Update environment system parameters
        !
        this%environment_system%lupdate = .TRUE.
        !
        CALL this%environment_system%update(this%system_system%com)
        !
        !--------------------------------------------------------------------------------
        ! Update cores
        !
        IF (setup%l1da) CALL setup%env_1da%update_origin(this%environment_system%com)
        !
        !--------------------------------------------------------------------------------
        ! Update rigid environ properties, defined on ions
        !
        IF (setup%lrigidcavity) THEN
            !
            IF (setup%lsolvent) THEN
                !
                CALL this%solvent%update()
                !
                !------------------------------------------------------------------------
                ! Update quantities that depend on the solvent boundary
                !
                IF (setup%lstatic) CALL this%static%update()
                !
                IF (setup%loptical) CALL this%optical%update()
                !
            END IF
            !
            IF (setup%lelectrolyte) THEN
                !
                CALL this%electrolyte%boundary%update()
                !
                CALL this%electrolyte%update()
                !
            END IF
            !
        END IF
        !
        !--------------------------------------------------------------------------------
        ! Update externals positions
        ! Note: if using ms-gcs, set positions of Helmholtz planes
        !
        IF (setup%lexternals) THEN
            !
            IF (setup%niter_ionic == 1) THEN
                !
                IF (setup%lmsgcs) THEN
                    extcharge_pos = 0.D0
                    extcharge_pos(3, 1) = MINVAL(this%system_ions%tau(3, :)) - 15.1178D0
                    extcharge_pos(3, 2) = MAXVAL(this%system_ions%tau(3, :)) + 15.1178D0
                END IF
                !
                CALL this%externals%functions%update(env_external_charges, &
                                                     extcharge_pos)
                !
            END IF
            !
            CALL this%externals%update()
            !
        END IF
        !
        IF (setup%lelectrostatic .OR. setup%lconfine) THEN
            !
            CALL this%system_charges%update()
            !
            CALL this%environment_charges%update()
            !
        END IF
        !
        this%system_system%lupdate = .FALSE.
        this%system_ions%lupdate = .FALSE.
        this%environment_system%lupdate = .FALSE.
        this%environment_ions%lupdate = .FALSE.
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE environ_update_ions
    !------------------------------------------------------------------------------------
    !>
    !! Initialize the electrons-related quantities to be used in the Environ
    !! modules. This initialization is called by electrons.f90, thus it
    !! is performed at every step of electronic optimization.
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE environ_update_electrons(this, nnr, rho, nelec)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: nnr
        REAL(DP), INTENT(IN) :: rho(nnr)
        REAL(DP), OPTIONAL, INTENT(IN) :: nelec
        !
        CLASS(environ_main), TARGET, INTENT(INOUT) :: this
        !
        REAL(DP), ALLOCATABLE :: aux(:)
        !
        TYPE(environ_setup), POINTER :: setup
        TYPE(environ_cell), POINTER :: environment_cell
        !
        !--------------------------------------------------------------------------------
        !
        setup => this%setup
        environment_cell => setup%environment_cell
        !
        !--------------------------------------------------------------------------------
        !
        this%system_electrons%lupdate = .TRUE.
        this%environment_electrons%lupdate = .TRUE.
        !
        !--------------------------------------------------------------------------------
        ! Update electrons parameters
        !
        CALL this%system_electrons%update(nnr, rho, nelec)
        !
        this%system_electrons%density%label = 'small_electrons'
        !
        IF (setup%ldoublecell) THEN
            ALLOCATE (aux(environment_cell%nnr))
            !
            CALL setup%mapping%to_large(nnr, environment_cell%nnr, rho, aux)
            !
            CALL this%environment_electrons%update(environment_cell%nnr, aux, nelec)
            !
        ELSE
            CALL this%environment_electrons%update(nnr, rho, nelec)
        END IF
        !
        this%environment_electrons%density%label = 'large_electrons'
        !
        IF (setup%lelectrostatic .OR. setup%lconfine) THEN
            !
            CALL this%system_charges%update()
            !
            CALL this%environment_charges%update()
            !
        END IF
        !
        !--------------------------------------------------------------------------------
        ! Update soft environ properties, defined on electrons
        !
        IF (setup%lsoftcavity) THEN
            !
            IF (setup%lsoftsolvent) THEN
                !
                CALL this%solvent%update()
                !
                !------------------------------------------------------------------------
                ! Update quantities that depend on the solvent boundary
                !
                IF (setup%lstatic) CALL this%static%update()
                !
                IF (setup%loptical) CALL this%optical%update()
                !
            END IF
            !
            IF (setup%lsoftelectrolyte) THEN
                !
                CALL this%electrolyte%boundary%update()
                !
                CALL this%electrolyte%update()
                !
            END IF
            !
        END IF
        !
        this%system_electrons%lupdate = .FALSE.
        this%environment_electrons%lupdate = .FALSE.
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE environ_update_electrons
    !------------------------------------------------------------------------------------
    !>
    !! Initialize the response charges to be used in the TDDFPT + Environ
    !! modules. This initialization is called by plugin_tddfpt_potential.f90
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE environ_update_response(this, nnr, drho)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: nnr
        REAL(DP), INTENT(IN) :: drho(nnr)
        !
        CLASS(environ_main), TARGET, INTENT(INOUT) :: this
        !
        REAL(DP), ALLOCATABLE :: aux(:)
        !
        TYPE(environ_setup), POINTER :: setup
        TYPE(environ_cell), POINTER :: environment_cell
        !
        !--------------------------------------------------------------------------------
        !
        setup => this%setup
        environment_cell => setup%environment_cell
        !
        IF (setup%optical_permittivity == 1.D0) RETURN
        !
        !--------------------------------------------------------------------------------
        !
        this%system_response_electrons%lupdate = .TRUE.
        this%environment_response_electrons%lupdate = .TRUE.
        !
        !--------------------------------------------------------------------------------
        ! Update response charges in system cell
        !
        CALL this%system_response_electrons%update(nnr, drho, 0.D0)
        !
        CALL this%system_response_charges%update()
        !
        !--------------------------------------------------------------------------------
        ! Update response charges in environment cell
        !
        IF (setup%ldoublecell) THEN
            !
            ALLOCATE (aux(environment_cell%nnr))
            !
            CALL setup%mapping%to_large(nnr, environment_cell%nnr, drho, aux)
            !
            CALL this%environment_response_electrons%update(environment_cell%nnr, &
                                                            aux, 0.D0)
            !
        ELSE
            CALL this%environment_response_electrons%update(nnr, drho, 0.D0)
        END IF
        !
        CALL this%environment_response_charges%update()
        !
        this%system_response_electrons%lupdate = .FALSE.
        this%environment_response_electrons%lupdate = .FALSE.
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE environ_update_response
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !                                   ACCESS METHODS
    !
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    FUNCTION get_vzero(this, nnr) RESULT(vzero)
        !--------------------------------------------------------------------------------
        !
        INTEGER, INTENT(IN) :: nnr
        CLASS(environ_main), INTENT(IN) :: this
        !
        REAL(DP) :: vzero(nnr)
        !
        CHARACTER(LEN=80) :: routine = 'get_vzero'
        !
        !--------------------------------------------------------------------------------
        !
        IF (nnr /= this%vzero%cell%nnr) &
            CALL io%error(routine, "Mismatch in grid size", 1)
        !
        !--------------------------------------------------------------------------------
        !
        vzero = this%vzero%of_r
        !
        !--------------------------------------------------------------------------------
    END FUNCTION get_vzero
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    FUNCTION get_dvtot(this, nnr) RESULT(dvtot)
        !--------------------------------------------------------------------------------
        !
        INTEGER, INTENT(IN) :: nnr
        CLASS(environ_main), INTENT(IN) :: this
        !
        REAL(DP) :: dvtot(nnr)
        !
        CHARACTER(LEN=80) :: routine = 'get_dvtot'
        !
        !--------------------------------------------------------------------------------
        !
        IF (nnr /= this%dvtot%cell%nnr) &
            CALL io%error(routine, "Mismatch in grid size", 1)
        !
        !--------------------------------------------------------------------------------
        !
        dvtot = this%dvtot%of_r
        !
        !--------------------------------------------------------------------------------
    END FUNCTION get_dvtot
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !                               PRIVATE HELPER METHODS
    !
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE environ_init_potential(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_main), TARGET, INTENT(INOUT) :: this
        !
        TYPE(environ_setup), POINTER :: setup
        TYPE(environ_cell), POINTER :: system_cell, environment_cell
        !
        CHARACTER(LEN=80) :: routine = 'environ_init_potential'
        !
        !--------------------------------------------------------------------------------
        !
        setup => this%setup
        system_cell => setup%system_cell
        environment_cell => setup%environment_cell
        !
        !--------------------------------------------------------------------------------
        ! Create local storage for base potential, that needs to be modified
        !
        CALL this%vzero%init(system_cell, 'vzero')
        !
        CALL this%dvtot%init(system_cell, 'dvtot')
        !
        !--------------------------------------------------------------------------------
        ! Electrostatic contribution
        !
        IF (setup%lelectrostatic) THEN
            !
            CALL this%velectrostatic%init(environment_cell, 'velectrostatic')
            !
            CALL this%vreference%init(system_cell, 'vreference')
            !
        END IF
        !
        !--------------------------------------------------------------------------------
        ! Contribution to the potential due to boundary
        !
        IF (setup%lsoftcavity) &
            CALL this%vsoftcavity%init(environment_cell, 'vsoftcavity')
        !
        !--------------------------------------------------------------------------------
        ! Confinement contribution
        !
        IF (setup%lconfine) CALL this%vconfine%init(environment_cell, 'vconfine')
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE environ_init_potential
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE environ_init_physical(this, nat, ntyp, ityp, zv, label, number, weight)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: nat, ntyp
        INTEGER, INTENT(IN) :: ityp(nat)
        REAL(DP), INTENT(IN) :: zv(ntyp)
        CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: label(ntyp)
        INTEGER, OPTIONAL, INTENT(IN) :: number(ntyp)
        REAL(DP), OPTIONAL, INTENT(IN) :: weight(ntyp)
        !
        CLASS(environ_main), TARGET, INTENT(INOUT) :: this
        !
        TYPE(environ_setup), POINTER :: setup
        TYPE(environ_cell), POINTER :: system_cell, environment_cell
        !
        CHARACTER(LEN=80) :: routine = 'environ_init_physical'
        !
        !--------------------------------------------------------------------------------
        !
        setup => this%setup
        system_cell => setup%system_cell
        environment_cell => setup%environment_cell
        !
        !--------------------------------------------------------------------------------
        ! Response
        !
        IF (setup%loptical) THEN
            !
            CALL this%system_response_electrons%init(system_cell)
            !
            CALL this%system_response_charges%init(system_cell)
            !
            CALL this%system_response_charges%add( &
                electrons=this%system_response_electrons)
            !
            CALL this%environment_response_electrons%init(environment_cell)
            !
            CALL this%environment_response_charges%init(environment_cell)
            !
            CALL this%environment_response_charges%add( &
                electrons=this%environment_response_electrons, &
                dielectric=this%optical)
            !
        END IF
        !
        !--------------------------------------------------------------------------------
        ! Ions
        !
        CALL this%system_ions%init(nat, ntyp, ityp, zv, atomicspread, corespread, &
                                   solvationrad, radius_mode, setup%lsoftcavity, &
                                   setup%lsmearedions, setup%lcoredensity, system_cell, &
                                   label, number, weight)
        !
        CALL this%environment_ions%init(nat, ntyp, ityp, zv, atomicspread, corespread, &
                                        solvationrad, radius_mode, setup%lsoftcavity, &
                                        setup%lsmearedions, setup%lcoredensity, &
                                        environment_cell, label, number, weight)
        !
        !--------------------------------------------------------------------------------
        ! Electrons
        !
        CALL this%system_electrons%init(system_cell)
        !
        CALL this%environment_electrons%init(environment_cell)
        !
        !--------------------------------------------------------------------------------
        ! System
        !
        CALL this%system_system%init(system_ntyp, system_dim, system_axis, &
                                     this%system_ions)
        !
        CALL this%environment_system%init(system_ntyp, system_dim, system_axis, &
                                          this%environment_ions)
        !
        !--------------------------------------------------------------------------------
        ! Free charges (if computing electrostatics or confinement)
        !
        CALL this%system_charges%init(system_cell)
        !
        CALL this%environment_charges%init(environment_cell)
        !
        IF (setup%lelectrostatic .OR. setup%lconfine) THEN
            !
            CALL this%system_charges%add(electrons=this%system_electrons)
            !
            CALL this%environment_charges%add(electrons=this%environment_electrons)
            !
        END IF
        !
        IF (setup%lelectrostatic) THEN
            !
            CALL this%system_charges%add(ions=this%system_ions)
            !
            CALL this%environment_charges%add(ions=this%environment_ions)
            !
        END IF
        !
        !--------------------------------------------------------------------------------
        ! External charges
        !
        IF (setup%lmsgcs) THEN
            !
            !----------------------------------------------------------------------------
            ! If ms-gcs calculation, add helmholtz planes
            !
            setup%lexternals = .TRUE.
            env_external_charges = 2
            !
            IF (ALLOCATED(extcharge_dim) .OR. &
                ALLOCATED(extcharge_axis) .OR. &
                ALLOCATED(extcharge_charge) .OR. &
                ALLOCATED(extcharge_spread) .OR. &
                ALLOCATED(extcharge_pos)) &
                CALL io%error(routine, "ms-gcs does not support user-defined external charges", 1)
            !
            ALLOCATE (extcharge_dim(env_external_charges))
            ALLOCATE (extcharge_axis(env_external_charges))
            ALLOCATE (extcharge_charge(env_external_charges))
            ALLOCATE (extcharge_spread(env_external_charges))
            ALLOCATE (extcharge_pos(3, env_external_charges))
            !
            extcharge_dim(1) = 2
            extcharge_axis(1) = 3
            extcharge_spread(1) = 0.25
            extcharge_charge(1) = 0.0
            !
            extcharge_dim(2) = 2
            extcharge_axis(2) = 3
            extcharge_spread(2) = 0.25
            extcharge_charge(2) = 0.0
        END IF
        !
        IF (setup%lexternals) THEN
            !
            CALL this%externals%init(env_external_charges, extcharge_dim, &
                                     extcharge_axis, extcharge_spread, &
                                     extcharge_charge, environment_cell)
            !
            CALL this%environment_charges%add(externals=this%externals)
            !
            DEALLOCATE (extcharge_axis)
            DEALLOCATE (extcharge_dim)
            DEALLOCATE (extcharge_spread)
            DEALLOCATE (extcharge_charge)
        END IF
        !
        !--------------------------------------------------------------------------------
        ! Solvent boundary
        !
        IF (setup%lsolvent) THEN
            !
            !----------------------------------------------------------------------------
            ! Casting and general setup
            !
            SELECT CASE (solvent_mode)
                !
            CASE ('electronic', 'full')
                ALLOCATE (environ_boundary_electronic :: this%solvent)
                !
            CASE ('ionic')
                ALLOCATE (environ_boundary_ionic :: this%solvent)
                !
            CASE ('system')
                ALLOCATE (environ_boundary_system :: this%solvent)
                !
            CASE DEFAULT
                CALL io%error(routine, "Unrecognized boundary mode", 1)
                !
            END SELECT
            !
            CALL this%solvent%pre_init( &
                solvent_mode, setup%lgradient, setup%need_factsqrt, setup%lsurface, &
                setup%outer_container, deriv_method, environment_cell, 'solvent')
            !
            !----------------------------------------------------------------------------
            ! Boundary awareness
            !
            IF (solvent_radius > 0.D0) &
                CALL this%solvent%init_solvent_aware(solvent_radius, radial_scale, &
                                                     radial_spread, filling_threshold, &
                                                     filling_spread)
            !
            IF (field_aware) &
                CALL this%solvent%init_field_aware(field_factor, field_asymmetry, &
                                                   field_max, field_min)
            !
            !----------------------------------------------------------------------------
            ! Specific setup
            !
            SELECT TYPE (solvent => this%solvent)
                !
            TYPE IS (environ_boundary_electronic)
                !
                CALL solvent%init(rhomax, rhomin, this%environment_electrons, &
                                  this%environment_ions)
                !
            TYPE IS (environ_boundary_ionic)
                !
                CALL solvent%init(alpha, softness, this%environment_ions, &
                                  this%environment_electrons)
                !
            TYPE IS (environ_boundary_system)
                !
                CALL solvent%init(solvent_distance, solvent_spread, &
                                  this%environment_system)
                !
            CLASS DEFAULT
                CALL io%error(routine, "Unrecognized boundary mode", 1)
                !
            END SELECT
            !
        END IF
        !
        !--------------------------------------------------------------------------------
        ! Electrolyte
        !
        IF (setup%lelectrolyte) THEN
            !
            CALL this%electrolyte%init( &
                env_electrolyte_ntyp, electrolyte_mode, electrolyte_rhomax, &
                electrolyte_rhomin, env_static_permittivity, &
                electrolyte_alpha, electrolyte_softness, electrolyte_distance, &
                electrolyte_spread, solvent_radius, radial_scale, radial_spread, &
                filling_threshold, filling_spread, field_aware, field_factor, &
                field_asymmetry, field_max, field_min, this%environment_electrons, &
                this%environment_ions, this%environment_system, temperature, cion, &
                cionmax, rion, zion, electrolyte_entropy, electrolyte_linearized, &
                setup%outer_container, electrolyte_deriv_method, environment_cell)
            !
            CALL this%environment_charges%add(electrolyte=this%electrolyte)
            !
        END IF
        !
        !--------------------------------------------------------------------------------
        ! Semiconductor
        !
        IF (setup%lsemiconductor) THEN
            !
            CALL this%semiconductor%init(temperature, sc_permittivity, &
                                         sc_carrier_density, sc_electrode_chg, &
                                         sc_distance, sc_spread, sc_chg_thr, &
                                         setup%lmsgcs, this%environment_system, &
                                         environment_cell)
            !
            CALL this%environment_charges%add(semiconductor=this%semiconductor)
            !
        END IF
        !
        !--------------------------------------------------------------------------------
        ! Dielectric
        !
        IF (setup%lstatic) THEN
            !
            CALL this%static%init(env_static_permittivity, this%solvent, &
                                  setup%need_gradient, setup%need_factsqrt, &
                                  setup%need_auxiliary, env_dielectric_regions, &
                                  environment_cell)
            !
            IF (this%static%nregions > 0) &
                CALL this%static%set_regions(env_dielectric_regions, epsregion_dim, &
                                             epsregion_axis, epsregion_pos, &
                                             epsregion_width, epsregion_spread, &
                                             epsregion_eps(1, :))
            !
            CALL this%environment_charges%add(dielectric=this%static)
            !
        END IF
        !
        IF (setup%loptical) THEN
            !
            CALL this%optical%init(env_optical_permittivity, this%solvent, &
                                   setup%need_gradient, setup%need_factsqrt, &
                                   setup%need_auxiliary, env_dielectric_regions, &
                                   environment_cell)
            !
            IF (this%optical%nregions > 0) &
                CALL this%optical%set_regions(env_dielectric_regions, epsregion_dim, &
                                              epsregion_axis, epsregion_pos, &
                                              epsregion_width, epsregion_spread, &
                                              epsregion_eps(2, :))
            !
        END IF
        !
        IF (ALLOCATED(epsregion_axis)) DEALLOCATE (epsregion_axis)
        !
        IF (ALLOCATED(epsregion_dim)) DEALLOCATE (epsregion_dim)
        !
        IF (ALLOCATED(epsregion_width)) DEALLOCATE (epsregion_width)
        !
        IF (ALLOCATED(epsregion_spread)) DEALLOCATE (epsregion_spread)
        !
        IF (ALLOCATED(epsregion_eps)) DEALLOCATE (epsregion_eps)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE environ_init_physical
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !                                   OUTPUT METHODS
    !
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !>
    !! Write out the different Environ contributions to the energy.
    !! Called by electrons.f90
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE print_environ_energies(this, prog, de_flag)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_main), TARGET, INTENT(IN) :: this
        CHARACTER(LEN=*), INTENT(IN) :: prog
        LOGICAL, OPTIONAL, INTENT(IN) :: de_flag
        !
        INTEGER, POINTER :: unit
        TYPE(environ_setup), POINTER :: setup
        !
        LOGICAL :: print_de
        !
        CHARACTER(LEN=80) :: routine = 'print_environ_energies'
        !
        !--------------------------------------------------------------------------------
        !
        IF (.NOT. io%lnode) RETURN
        !
        IF (PRESENT(de_flag)) THEN
            print_de = de_flag
        ELSE
            print_de = .TRUE.
        END IF
        !
        unit => io%unit
        setup => this%setup
        !
        SELECT CASE (prog)
            !
        CASE ('PW')
            !
            IF (setup%lelectrostatic) WRITE (unit, 1000) this%eelectrostatic
            !
            IF (setup%lsurface) WRITE (unit, 1001) this%esurface
            !
            IF (setup%lvolume) WRITE (unit, 1002) this%evolume
            !
            IF (setup%lconfine) WRITE (unit, 1003) this%econfine
            !
            IF (setup%lelectrolyte) WRITE (unit, 1004) this%eelectrolyte
            !
            IF (print_de) WRITE (unit, 1005) this%deenviron
            !
        CASE ('CP') ! converted to Hartree
            !
            IF (setup%lelectrostatic) WRITE (unit, 1006) this%eelectrostatic * 0.5D0
            !
            IF (setup%lsurface) WRITE (unit, 1007) this%esurface * 0.5D0
            !
            IF (setup%lvolume) WRITE (unit, 1008) this%evolume * 0.5D0
            !
            IF (setup%lconfine) WRITE (unit, 1009) this%econfine * 0.5D0
            !
            IF (setup%lelectrolyte) WRITE (unit, 1010) this%eelectrolyte * 0.5D0
            !
        CASE DEFAULT
            CALL io%error(routine, "Unexpected calling program", 1)
            !
        END SELECT
        !
        IF (lvolume .OR. lsurface) THEN
            WRITE(unit, *) ''
            !
            IF (lvolume) WRITE(unit, 1011) this%solvent%volume
            !
            IF (lsurface) WRITE(unit, 1012) this%solvent%surface
            !
            WRITE(unit, *) ''
        END IF
        !
        !--------------------------------------------------------------------------------
        !
1000    FORMAT("     electrostatic embedding   =", F17.8, " Ry")
1001    FORMAT("     cavitation energy         =", F17.8, " Ry")
1002    FORMAT("     PV energy                 =", F17.8, " Ry")
1003    FORMAT("     confinement energy        =", F17.8, " Ry")
1004    FORMAT("     electrolyte free energy   =", F17.8, " Ry")
1005    FORMAT("     correction to one-el term =", F17.8, " Ry")
        !
1006    FORMAT("     electrostatic embedding = ", F14.5, " Hartree a.u.")
1007    FORMAT("           cavitation energy = ", F14.5, " Hartree a.u.")
1008    FORMAT("                   PV energy = ", F14.5, " Hartree a.u.")
1009    FORMAT("     electrolyte free energy = ", F14.5, " Hartree a.u.")
1010    FORMAT("          confinement energy = ", F14.5, " Hartree a.u.")
        !
1011    FORMAT(5X, "Total volume of the QM region = ", F18.8, " Bohr^3")
1012    FORMAT(5X, "Total surface of the QM region = ", F17.8, " Bohr^2")
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE print_environ_energies
    !------------------------------------------------------------------------------------
    !>
    !! If Gaussian nuclei are used, write out the corresponding potential shift
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE print_environ_potential_shift(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_main), INTENT(IN) :: this
        !
        !--------------------------------------------------------------------------------
        !
        IF (.NOT. io%lnode) RETURN
        !
        IF (this%setup%lsmearedions) &
            WRITE (io%unit, 1100) this%environment_ions%potential_shift * RYTOEV
        !
1100    FORMAT(/, 5(' '), &
                "the potential shift due to the Gaussian-smeared nuclei is ", &
                F10.4, " ev")
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE print_environ_potential_shift
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE class_environ
!----------------------------------------------------------------------------------------
