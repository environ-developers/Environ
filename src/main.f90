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
    USE env_base_input
    USE class_setup
    !
    USE class_cell
    USE class_density
    USE class_gradient
    !
    USE class_boundary
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
    TYPE, PUBLIC :: environ_obj
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
        TYPE(environ_boundary) :: solvent
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
    CONTAINS
        !--------------------------------------------------------------------------------
        !
        PROCEDURE, PRIVATE :: create => create_environ_base
        PROCEDURE :: init => init_environ_base
        !
        PROCEDURE, PRIVATE :: init_physical => environ_init_physical
        PROCEDURE, PRIVATE :: init_potential => environ_init_potential
        !
        PROCEDURE :: update_electrons => environ_update_electrons
        PROCEDURE :: update_ions => environ_update_ions
        PROCEDURE :: update_potential => environ_update_potential
        PROCEDURE :: update_response => environ_update_response
        PROCEDURE :: update_cell_dependent_quantities
        !
        PROCEDURE :: print_energies => print_environ_energies
        PROCEDURE :: print_potential_shift => print_environ_potential_shift
        !
        !--------------------------------------------------------------------------------
    END TYPE environ_obj
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
        CLASS(environ_obj), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: sub_name = 'create_environ_base'
        !
        !--------------------------------------------------------------------------------
        !
        IF (ASSOCIATED(this%setup)) &
            CALL io%error(sub_name, 'Trying to create an existing object', 1)
        !
        !--------------------------------------------------------------------------------
        !
        NULLIFY (this%setup)
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
    SUBROUTINE init_environ_base(this, setup, nelec, nat, ntyp, atom_label, ityp, zv)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_setup), TARGET, INTENT(IN) :: setup
        !
        INTEGER, INTENT(IN) :: nat, ntyp
        INTEGER, INTENT(IN) :: nelec, ityp(nat)
        REAL(DP), INTENT(IN) :: zv(ntyp)
        CHARACTER(LEN=3), INTENT(IN) :: atom_label(:)
        !
        CLASS(environ_obj), TARGET, INTENT(INOUT) :: this
        !
        !--------------------------------------------------------------------------------
        ! Set verbosity and open debug file
        !
        io%verbosity = verbose ! set internal verbosity from input
        !
        IF (io%verbosity >= 1) THEN
            io%debug_unit = io%find_free_unit()
            !
            OPEN (unit=io%debug_unit, file='environ.debug', status='unknown')
            !
        END IF
        !
        !--------------------------------------------------------------------------------
        !
        this%setup => setup
        !
        CALL this%init_potential()
        !
        CALL this%init_physical(nelec, nat, ntyp, atom_label, ityp, zv)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_environ_base
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
        CLASS(environ_obj), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: sub_name = 'environ_update_potential'
        !
        !--------------------------------------------------------------------------------
        !
        this%vzero%lupdate = .TRUE.
        !
        IF (.NOT. ASSOCIATED(this%vzero%cell)) RETURN
        !
        IF (this%vzero%cell%nnr /= nnr) &
            CALL io%error(sub_name, 'Inconsistent size in input potential', 1)
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
        CLASS(environ_obj), TARGET, INTENT(INOUT) :: this
        !
        TYPE(environ_setup), POINTER :: setup
        !
        CHARACTER(LEN=80) :: sub_name = 'update_cell_dependent_quantities'
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
    SUBROUTINE environ_update_ions(this, nat, tau)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: nat
        REAL(DP), INTENT(IN) :: tau(3, nat)
        !
        CLASS(environ_obj), TARGET, INTENT(INOUT) :: this
        !
        REAL(DP) :: local_pos(3)
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
        !--------------------------------------------------------------------------------
        ! Update system ions parameters
        !
        CALL this%system_ions%update(nat, tau)
        !
        !--------------------------------------------------------------------------------
        ! Update system system parameters
        !
        this%system_system%lupdate = .TRUE.
        !
        CALL this%system_system%update()
        !
        !--------------------------------------------------------------------------------
        ! Update cell mapping
        !
        IF (.NOT. setup%mapping%initialized) THEN
            !
            IF (setup%ldoublecell) THEN
                !
                IF (environ_debug) THEN
                    local_pos = mapping_pos ! debugging with finite-differences
                ELSE
                    local_pos = this%system_system%pos ! center of charge
                END IF
                !
                CALL setup%update_mapping(local_pos)
                !
            ELSE
                setup%mapping%initialized = .TRUE. ! one-to-one mapping
            END IF
            !
        END IF
        !
        !--------------------------------------------------------------------------------
        ! Update environment ions parameters
        !
        CALL this%environment_ions%update(nat, tau)
        !
        !--------------------------------------------------------------------------------
        ! Update environment system parameters
        !
        this%environment_system%lupdate = .TRUE.
        !
        CALL this%environment_system%update()
        !
        !--------------------------------------------------------------------------------
        ! Set soft-sphere parameters
        !
        IF (setup%lsolvent) CALL this%solvent%set_soft_spheres()
        !
        IF (setup%lelectrolyte) CALL this%electrolyte%boundary%set_soft_spheres()
        !
        !--------------------------------------------------------------------------------
        ! Update cores
        !
        IF (setup%l1da) &
            CALL setup%core_1da_elect%update_origin(this%environment_system%pos)
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
        ! External charges rely on the environment cell, which is defined
        ! with respect to the system origin
        !
        IF (setup%lexternals) CALL this%externals%update()
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
        REAL(DP), INTENT(IN), OPTIONAL :: nelec
        !
        CLASS(environ_obj), TARGET, INTENT(INOUT) :: this
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
        CLASS(environ_obj), TARGET, INTENT(INOUT) :: this
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
        CLASS(environ_obj), TARGET, INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: local_label
        !
        TYPE(environ_setup), POINTER :: setup
        TYPE(environ_cell), POINTER :: system_cell, environment_cell
        !
        CHARACTER(LEN=80) :: sub_name = 'environ_init_potential'
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
        local_label = 'vzero'
        !
        CALL this%vzero%init(system_cell, local_label)
        !
        local_label = 'dvtot'
        !
        CALL this%dvtot%init(system_cell, local_label)
        !
        !--------------------------------------------------------------------------------
        ! Electrostatic contribution
        !
        IF (setup%lelectrostatic) THEN
            local_label = 'velectrostatic'
            !
            CALL this%velectrostatic%init(environment_cell, local_label)
            !
            local_label = 'vreference'
            !
            CALL this%vreference%init(system_cell, local_label)
            !
        END IF
        !
        !--------------------------------------------------------------------------------
        ! Contribution to the potential due to boundary
        !
        IF (setup%lsoftcavity) THEN
            local_label = 'vsoftcavity'
            !
            CALL this%vsoftcavity%init(environment_cell, local_label)
            !
        END IF
        !
        !--------------------------------------------------------------------------------
        ! Confinement contribution
        !
        IF (setup%lconfine) THEN
            local_label = 'vconfine'
            !
            CALL this%vconfine%init(environment_cell, local_label)
            !
        END IF
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE environ_init_potential
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE environ_init_physical(this, nelec, nat, ntyp, atom_label, ityp, zv)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: nelec, nat, ntyp
        INTEGER, INTENT(IN) :: ityp(nat)
        REAL(DP), INTENT(IN) :: zv(ntyp)
        CHARACTER(LEN=3), INTENT(IN) :: atom_label(:)
        !
        CLASS(environ_obj), TARGET, INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: local_label
        !
        TYPE(environ_setup), POINTER :: setup
        TYPE(environ_cell), POINTER :: system_cell, environment_cell
        !
        CHARACTER(LEN=80) :: sub_name = 'environ_init_physical'
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
            CALL this%system_response_electrons%init(0, system_cell)
            !
            CALL this%system_response_charges%init(system_cell)
            !
            CALL this%system_response_charges%add( &
                electrons=this%system_response_electrons)
            !
            CALL this%environment_response_electrons%init(0, environment_cell)
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
        CALL this%system_ions%init(nat, ntyp, atom_label, ityp, zv, atomicspread, &
                                   corespread, solvationrad, radius_mode, &
                                   setup%lsoftcavity, setup%lsmearedions, &
                                   setup%lcoredensity, system_cell)
        !
        CALL this%environment_ions%init(nat, ntyp, atom_label, ityp, zv, atomicspread, &
                                        corespread, solvationrad, radius_mode, &
                                        setup%lsoftcavity, setup%lsmearedions, &
                                        setup%lcoredensity, environment_cell)
        !
        !--------------------------------------------------------------------------------
        ! Electrons
        !
        CALL this%system_electrons%init(nelec, system_cell)
        !
        CALL this%environment_electrons%init(nelec, environment_cell)
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
        IF (setup%lexternals) THEN
            !
            CALL this%externals%init(env_external_charges, extcharge_dim, &
                                     extcharge_axis, extcharge_pos, &
                                     extcharge_spread, extcharge_charge, &
                                     environment_cell)
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
            local_label = 'solvent'
            !
            CALL this%solvent%init( &
                setup%lgradient, setup%need_factsqrt, setup%lsurface, &
                solvent_mode, stype, rhomax, rhomin, tbeta, env_static_permittivity, &
                alpha, softness, solvent_distance, solvent_spread, solvent_radius, &
                radial_scale, radial_spread, filling_threshold, filling_spread, &
                field_awareness, charge_asymmetry, field_max, field_min, &
                this%environment_electrons, this%environment_ions, &
                this%environment_system, setup%derivatives, environment_cell, &
                local_label)
            !
        END IF
        !
        !--------------------------------------------------------------------------------
        ! Electrolyte
        !
        IF (setup%lelectrolyte) THEN
            !
            CALL this%electrolyte%init( &
                env_electrolyte_ntyp, electrolyte_mode, stype, electrolyte_rhomax, &
                electrolyte_rhomin, electrolyte_tbeta, env_static_permittivity, &
                electrolyte_alpha, electrolyte_softness, electrolyte_distance, &
                electrolyte_spread, solvent_radius, radial_scale, radial_spread, &
                filling_threshold, filling_spread, field_awareness, charge_asymmetry, &
                field_max, field_min, this%environment_electrons, this%environment_ions, &
                this%environment_system, setup%derivatives, temperature, cion, &
                cionmax, rion, zion, electrolyte_entropy, electrolyte_linearized, &
                environment_cell)
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
                                         this%environment_system, environment_cell)
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
    SUBROUTINE print_environ_energies(this, prog)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CHARACTER(LEN=2), INTENT(IN) :: prog
        !
        CLASS(environ_obj), TARGET, INTENT(INOUT) :: this
        !
        INTEGER, POINTER :: unit
        TYPE(environ_setup), POINTER :: setup
        !
        CHARACTER(LEN=80) :: sub_name = 'print_environ_energies'
        !
        !--------------------------------------------------------------------------------
        !
        unit => io%unit
        setup => this%setup
        !
        IF (io%lnode) THEN
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
                WRITE (unit, 1005) this%deenviron
                !
            CASE ('CP')
                !
                IF (setup%lelectrostatic) WRITE (unit, 1006) this%eelectrostatic
                !
                IF (setup%lsurface) WRITE (unit, 1007) this%esurface
                !
                IF (setup%lvolume) WRITE (unit, 1008) this%evolume
                !
                IF (setup%lconfine) WRITE (unit, 1009) this%econfine
                !
                IF (setup%lelectrolyte) WRITE (unit, 1010) this%eelectrolyte
                !
            CASE DEFAULT
                CALL io%error(sub_name, 'Wrong program calling Environ', 1)
                !
            END SELECT
            !
        END IF
        !
        !--------------------------------------------------------------------------------
        !
1000    FORMAT('     electrostatic embedding   =', F17.8, ' Ry')
1001    FORMAT('     cavitation energy         =', F17.8, ' Ry')
1002    FORMAT('     PV energy                 =', F17.8, ' Ry')
1003    FORMAT('     confinement energy        =', F17.8, ' Ry')
1004    FORMAT('     electrolyte free energy   =', F17.8, ' Ry')
1005    FORMAT('     correction to one-el term =', F17.8, ' Ry')
        !
1006    FORMAT('     electrostatic embedding = ', F14.5, ' Hartree a.u.')
1007    FORMAT('           cavitation energy = ', F14.5, ' Hartree a.u.')
1008    FORMAT('                   PV energy = ', F14.5, ' Hartree a.u.')
1009    FORMAT('     electrolyte free energy = ', F14.5, ' Hartree a.u.')
1010    FORMAT('          confinement energy = ', F14.5, ' Hartree a.u.')
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
        CLASS(environ_obj), INTENT(INOUT) :: this
        !
        !--------------------------------------------------------------------------------
        !
        IF (this%setup%lsmearedions) &
            WRITE (io%unit, 1100) this%environment_ions%potential_shift * RYTOEV
        !
1100    FORMAT(/, 5(' '), &
                'the potential shift due to the Gaussian-smeared nuclei is ', &
                F10.4, ' ev')
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE print_environ_potential_shift
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE class_environ
!----------------------------------------------------------------------------------------
