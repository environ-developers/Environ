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
    USE env_base_io
    !
    USE environ_param, ONLY: DP, BOHR_RADIUS_SI, RYDBERG_SI
    !
    USE class_cell
    USE class_density
    USE class_gradient
    USE class_mapping
    !
    USE class_core_container_corrections
    USE class_core_container_derivatives
    USE class_core_container_electrostatics
    USE class_core_fd_derivatives
    USE class_core_fft_derivatives
    USE class_core_fft_electrostatics
    USE class_core_1da_electrostatics
    !
    USE class_solver_setup
    USE class_solver
    USE class_solver_direct
    USE class_solver_fixedpoint
    USE class_solver_gradient
    USE class_solver_iterative
    USE class_solver_newton
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
    USE env_base_input
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
    TYPE environ_obj
        !--------------------------------------------------------------------------------
        !
        !--------------------------------------------------------------------------------
        ! Environ execution parameters
        !
        LOGICAL :: environ_restart
        REAL(DP) :: environ_thr
        INTEGER :: environ_nskip
        !
        !--------------------------------------------------------------------------------
        ! Control flags
        !
        LOGICAL :: ltddfpt, lsolvent, lelectrostatic, lsoftsolvent, lsoftelectrolyte, &
                   lsoftcavity, lrigidsolvent, lrigidelectrolyte, lrigidcavity, &
                   lcoredensity, lsmearedions, lgradient
        !
        !--------------------------------------------------------------------------------
        ! Internal parameters for mapping between environment and system cell
        !
        LOGICAL :: ldoublecell
        TYPE(environ_mapping) :: mapping
        TYPE(environ_cell) :: system_cell
        TYPE(environ_cell), POINTER :: environment_cell
        TYPE(environ_ions) :: system_ions, environment_ions
        TYPE(environ_electrons) :: system_electrons, environment_electrons
        !
        !--------------------------------------------------------------------------------
        !
        LOGICAL :: lexternals
        TYPE(environ_externals) :: externals ! external charges
        TYPE(environ_charges) :: system_charges, environment_charges
        TYPE(environ_system) :: system_system, environment_system
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
        ! Details of the continuum interface
        !
        LOGICAL :: lboundary
        TYPE(environ_boundary) :: solvent
        TYPE(container_derivatives) :: derivatives
        !
        !--------------------------------------------------------------------------------
        ! Dielectric parameters (solvent)
        !
        LOGICAL :: lstatic, loptical, ldielectric
        REAL(DP) :: env_static_permittivity, env_optical_permittivity
        TYPE(environ_dielectric) :: static, optical
        !
        !--------------------------------------------------------------------------------
        ! Ionic countercharge parameters
        !
        LOGICAL :: lelectrolyte
        INTEGER :: env_electrolyte_ntyp
        TYPE(environ_electrolyte) :: electrolyte
        !
        !--------------------------------------------------------------------------------
        ! Semiconductor parameters
        !
        LOGICAL :: lsemiconductor, louterloop
        TYPE(environ_semiconductor) :: semiconductor
        !
        !--------------------------------------------------------------------------------
        ! Cavitation energy parameters
        !
        LOGICAL :: lsurface
        REAL(DP) :: env_surface_tension
        !
        !--------------------------------------------------------------------------------
        ! PV term parameters
        !
        LOGICAL :: lvolume
        REAL(DP) :: env_pressure
        !
        !--------------------------------------------------------------------------------
        ! Confinement potential parameters
        !
        LOGICAL :: lconfine
        REAL(DP) :: env_confine
        !
        !--------------------------------------------------------------------------------
        ! Periodicity correction parameters
        !
        LOGICAL :: lperiodic
        !
        !--------------------------------------------------------------------------------
        ! Temporary parameters
        !
        INTEGER :: nrep
        INTEGER :: niter ! stores the iteration of environ for debugging purposes
        !
        !--------------------------------------------------------------------------------
        ! Computed physical variables
        !
        REAL(DP) :: deenviron, eelectrostatic, esurface, evolume, eelectrolyte, econfine
        !
        TYPE(environ_density) :: &
            vzero, velectrostatic, vreference, dvtot, vconfine, vsoftcavity
        !
        !--------------------------------------------------------------------------------
        ! Electrostatic
        !
        TYPE(electrostatic_setup) :: reference, outer, inner
        TYPE(container_electrostatics) :: reference_cores, outer_cores, inner_cores
        TYPE(container_corrections) :: pbc_core
        !
        !--------------------------------------------------------------------------------
        ! Internal setup of numerical solvers
        !
        TYPE(solver_direct) :: reference_direct, direct
        TYPE(solver_gradient) :: gradient, inner_gradient
        TYPE(solver_fixedpoint) :: fixedpoint, inner_fixedpoint
        TYPE(solver_newton) :: newton
        !
        !--------------------------------------------------------------------------------
        ! General setup of periodic boundary conditions
        !
        LOGICAL :: need_pbc_correction, need_electrolyte, need_semiconductor, &
                   need_outer_loop
        !
        !--------------------------------------------------------------------------------
        !
        LOGICAL :: need_gradient, need_factsqrt, need_auxiliary
        ! logical flags that need to be used outside
        !
        !--------------------------------------------------------------------------------
        ! Numerical cores
        !
        LOGICAL :: lfd
        TYPE(core_fd_derivatives) :: core_fd
        !
        LOGICAL :: lfft_system
        TYPE(core_fft_electrostatics) :: core_fft_sys
        !
        LOGICAL :: lfft_environment
        TYPE(core_fft_derivatives) :: core_fft_deriv
        TYPE(core_fft_electrostatics) :: core_fft_elect
        !
        LOGICAL :: l1da
        TYPE(core_1da_electrostatics) :: core_1da_elect
        !
        !--------------------------------------------------------------------------------
    CONTAINS
        !--------------------------------------------------------------------------------
        !
        PROCEDURE :: init_first => set_environ_bases
        PROCEDURE :: init_second => environ_initbase
        PROCEDURE :: init_potential => environ_initpotential
        PROCEDURE :: init_cell => environ_initcell
        PROCEDURE :: init_ions => environ_initions
        PROCEDURE :: init_electrons => environ_initelectrons
        PROCEDURE :: init_response => environ_initresponse
        !
        PROCEDURE, PRIVATE :: set_core_base, set_electrostatic_base, set_environ_base
        !
        PROCEDURE :: potential => calc_venviron
        PROCEDURE :: energy => calc_eenviron
        PROCEDURE :: force => calc_fenviron
        PROCEDURE :: dpotential => calc_dvenviron
        PROCEDURE :: denergy => calc_deenviron
        !
        PROCEDURE :: summary => environ_summary
        !
        !--------------------------------------------------------------------------------
    END TYPE environ_obj
    !------------------------------------------------------------------------------------
    !
    CHARACTER(LEN=256) :: bibliography(4)
    !
    DATA bibliography/ &
        "O. Andreussi, I. Dabo and N. Marzari, &
        &J. Chem. Phys. 136, 064102 (2012)", &
        "I. Timrov, O. Andreussi, A. Biancardi, N. Marzari, and S. Baroni, &
        &J. Chem. Phys. 142, 034111 (2015)", &
        "O. Andreussi, N.G. Hoermann, F. Nattino, G. Fisicaro, S. Goedecker, and N. Marzari, &
        &J. Chem. Theory Comput. 15, 1996 (2019)", &
        "F. Nattino, M. Truscott, N. Marzari, and O. Andreussi, &
        &J. Chem. Phys. 150, 041722 (2019)"/
    !
    TYPE(environ_obj), PUBLIC, SAVE :: env
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
    SUBROUTINE set_environ_bases(this, nelec, nat, ntyp, atom_label, &
                                 use_internal_pbc_corr)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        LOGICAL, INTENT(IN) :: use_internal_pbc_corr
        INTEGER, INTENT(IN) :: nelec, nat, ntyp
        CHARACTER(LEN=3), INTENT(IN) :: atom_label(:)
        !
        CLASS(environ_obj), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: sub_name = 'set_environ_bases'
        !
        !--------------------------------------------------------------------------------
        ! Set attributes (DO NOT CHANGE THE CALL ORDER)
        !
        CALL this%set_electrostatic_base()
        !
        CALL this%set_environ_base(nelec, nat, ntyp, atom_label)
        !
        CALL this%set_core_base(use_internal_pbc_corr)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE set_environ_bases
    !------------------------------------------------------------------------------------
    !>
    !! Subroutine to initialize fundamental quantities needed by the
    !! environ modules. This subroutine is called by init_run.f90, thus
    !! only once per pw.x execution.
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE environ_initbase(this, at, comm_in, me, root, gcutm)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: comm_in, me, root
        REAL(DP), INTENT(IN) :: at(3, 3)
        REAL(DP), INTENT(IN) :: gcutm
        !
        CLASS(environ_obj), TARGET, INTENT(INOUT) :: this
        !
        INTEGER :: m(3)
        !
        INTEGER :: ipol
        INTEGER :: environment_nr(3)
        REAL(DP) :: environment_at(3, 3)
        !
        CHARACTER(LEN=80) :: label = ' '
        !
        !--------------------------------------------------------------------------------
        ! Common initialization for simulations with Environ
        !
        CALL this%system_cell%init(gcutm, comm_in, at)
        !
        !--------------------------------------------------------------------------------
        ! Double cell and mapping
        !
        IF (this%ldoublecell) THEN
            !
            !----------------------------------------------------------------------------
            ! Scale environment lattice (and corresponding ffts) by 2 * nrep(i) + 1
            !
            DO ipol = 1, 3
                !
                environment_at(:, ipol) = at(:, ipol) * &
                                          (2.D0 * this%mapping%nrep(ipol) + 1.D0)
                !
            END DO
            !
            environment_nr(1) = this%system_cell%dfft%nr1 * (2 * this%mapping%nrep(1) + 1)
            environment_nr(2) = this%system_cell%dfft%nr2 * (2 * this%mapping%nrep(2) + 1)
            environment_nr(3) = this%system_cell%dfft%nr3 * (2 * this%mapping%nrep(3) + 1)
            !
            CALL this%environment_cell%init(gcutm, comm_in, environment_at, &
                                            environment_nr)
            !
        ELSE
            this%environment_cell => this%system_cell
        END IF
        !
        CALL this%mapping%init_second(this%system_cell, this%environment_cell)
        !
        !--------------------------------------------------------------------------------
        ! Initialize numerical cores
        !
        IF (this%lfd) CALL this%core_fd%init_second(this%environment_cell)
        !
        IF (this%lfft_system) CALL this%core_fft_sys%init_second(gcutm, this%system_cell)
        !
        IF (this%lfft_environment) THEN
            !
            CALL this%core_fft_deriv%init_second(gcutm, this%environment_cell)
            !
            IF (this%lelectrostatic) &
                CALL this%core_fft_elect%init_second(gcutm, this%environment_cell)
            !
        END IF
        !
        IF (this%l1da) &
            CALL this%core_1da_elect%init_second(this%environment_cell)
        !
        !--------------------------------------------------------------------------------
        ! Create local storage for base potential, that needs to be modified
        !
        label = 'vzero'
        !
        CALL this%vzero%create(label)
        !
        CALL this%vzero%init(this%system_cell)
        !
        this%deenviron = 0.0_DP
        !
        label = 'dvtot'
        !
        CALL this%dvtot%create(label)
        !
        CALL this%dvtot%init(this%system_cell)
        !
        !--------------------------------------------------------------------------------
        ! Electrostatic contribution
        !
        this%eelectrostatic = 0.0_DP
        !
        IF (this%lelectrostatic) THEN
            !
            label = 'velectrostatic'
            !
            CALL this%velectrostatic%create(label)
            !
            CALL this%velectrostatic%init(this%environment_cell)
            !
            label = 'vreference'
            !
            CALL this%vreference%create(label)
            !
            CALL this%vreference%init(this%system_cell)
            !
        END IF
        !
        !--------------------------------------------------------------------------------
        ! Contribution to the potential due to boundary
        !
        IF (this%lsoftcavity) THEN
            !
            label = 'vsoftcavity'
            !
            CALL this%vsoftcavity%create(label)
            !
            CALL this%vsoftcavity%init(this%environment_cell)
            !
        END IF
        !
        !--------------------------------------------------------------------------------
        ! Surface & volume contributions
        !
        this%esurface = 0.0_DP ! cavity contribution
        this%evolume = 0.0_DP ! pressure contribution
        !
        !--------------------------------------------------------------------------------
        ! Confinement contribution
        !
        this%econfine = 0.0_DP
        !
        IF (this%lconfine) THEN
            label = 'vconfine'
            !
            CALL this%vconfine%create(label)
            !
            CALL this%vconfine%init(this%environment_cell)
            !
        END IF
        !
        !--------------------------------------------------------------------------------
        !
        this%eelectrolyte = 0.0_DP ! non-electrostatic electrolyte contribution
        !
        !--------------------------------------------------------------------------------
        ! Second step of initialization of some environ derived type
        !
        CALL this%system_electrons%init_second(this%system_cell)
        !
        CALL this%environment_electrons%init_second(this%environment_cell)
        !
        IF (this%lsolvent) CALL this%solvent%init_second(this%environment_cell)
        !
        IF (this%lstatic) CALL this%static%init_second(this%environment_cell)
        !
        IF (this%loptical) THEN
            !
            CALL this%system_response_electrons%init_second(this%system_cell)
            !
            CALL this%system_response_charges%init_second(this%system_cell)
            !
            CALL this%environment_response_electrons%init_second(this%environment_cell)
            !
            CALL this%environment_response_charges%init_second(this%environment_cell)
            !
            CALL this%optical%init_second(this%environment_cell)
            !
        END IF
        !
        IF (this%lelectrolyte) CALL this%electrolyte%init_second(this%environment_cell)
        !
        IF (this%lexternals) CALL this%externals%init_second(this%environment_cell)
        !
        IF (this%lelectrostatic .OR. this%lconfine) THEN
            !
            CALL this%system_charges%init_second(this%system_cell)
            !
            CALL this%environment_charges%init_second(this%environment_cell)
            !
        END IF
        !
        IF (this%lsemiconductor) THEN
            !
            CALL this%semiconductor%init_second(this%environment_cell)
            !
        END IF
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE environ_initbase
    !------------------------------------------------------------------------------------
    !>
    !! Save local potential that will be overwritten by environ
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE environ_initpotential(this, nnr, vltot)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: nnr
        REAL(DP), INTENT(IN) :: vltot(nnr)
        !
        CLASS(environ_obj), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: sub_name = 'environ_initpotential'
        !
        !--------------------------------------------------------------------------------
        !
        this%vzero%lupdate = .TRUE.
        !
        IF (.NOT. ASSOCIATED(this%vzero%cell)) RETURN
        !
        IF (this%vzero%cell%nnr /= nnr) &
            CALL env_errore(sub_name, 'Inconsistent size in input potential', 1)
        !
        this%vzero%of_r = vltot
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE environ_initpotential
    !------------------------------------------------------------------------------------
    !>
    !! Initialize the cell-related quantities to be used in the Environ
    !! modules. This initialization is called by electrons.f90, thus it
    !! is performed at every step of electronic optimization.
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE environ_initcell(this, at)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP), INTENT(IN) :: at(3, 3)
        !
        CLASS(environ_obj), INTENT(INOUT) :: this
        !
        INTEGER :: ipol
        REAL(DP) :: environment_at(3, 3)
        !
        !--------------------------------------------------------------------------------
        !
        this%system_cell%lupdate = .TRUE.
        !
        CALL this%system_cell%update(at)
        !
        IF (this%ldoublecell) THEN
            this%environment_cell%lupdate = .TRUE.
            !
            DO ipol = 1, 3
                !
                environment_at(:, ipol) = at(:, ipol) * &
                                          (2.D0 * this%mapping%nrep(ipol) + 1.D0)
                !
            END DO
            !
            CALL this%environment_cell%update(environment_at)
            !
        END IF
        !
        !--------------------------------------------------------------------------------
        ! Update cores
        !
        IF (this%lfft_system) CALL this%core_fft_sys%update_cell(this%system_cell)
        !
        IF (this%lfft_environment) THEN
            !
            CALL this%core_fft_deriv%update_cell(this%environment_cell)
            !
            IF (this%lelectrostatic) &
                CALL this%core_fft_elect%update_cell(this%environment_cell)
            !
        END IF
        !
        IF (this%l1da) &
            CALL this%core_1da_elect%update_cell(this%environment_cell)
        !
        !--------------------------------------------------------------------------------
        ! Update fixed quantities defined inside the cell
        !
        IF (this%lstatic) CALL this%static%update()
        !
        IF (this%loptical) CALL this%optical%update()
        !
        IF (this%lelectrolyte) CALL this%electrolyte%update()
        !
        IF (this%lexternals) CALL this%externals%update()
        !
        IF (this%lsemiconductor) CALL this%semiconductor%update()
        !
        this%system_cell%lupdate = .FALSE.
        !
        IF (this%ldoublecell) this%environment_cell%lupdate = .FALSE.
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE environ_initcell
    !------------------------------------------------------------------------------------
    !>
    !! Initialize the ions-related quantities to be used in the Environ
    !! modules. This initialization is called by electrons.f90, thus it
    !! is performed at every step of electronic optimization. It may not
    !! be the most efficient choice, but it is a safe choice.
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE environ_initions(this, nat, ntyp, ityp, zv, tau)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: nat, ntyp
        INTEGER, INTENT(IN) :: ityp(nat)
        REAL(DP), INTENT(IN) :: zv(ntyp), tau(3, nat)
        !
        CLASS(environ_obj), INTENT(INOUT) :: this
        !
        !--------------------------------------------------------------------------------
        !
        this%system_ions%lupdate = .TRUE.
        this%environment_ions%lupdate = .TRUE.
        !
        !--------------------------------------------------------------------------------
        !
        CALL this%system_ions%init_second(nat, ntyp, ityp, zv, this%system_cell)
        ! second step of initialization for system ions
        !
        !--------------------------------------------------------------------------------
        ! Update system ions parameters
        !
        CALL this%system_ions%update(nat, tau)
        !
        CALL this%system_ions%printout()
        !
        !--------------------------------------------------------------------------------
        ! Update system system parameters
        !
        this%system_system%lupdate = .TRUE.
        !
        CALL this%system_system%update()
        !
        CALL this%system_system%printout()
        !
        !--------------------------------------------------------------------------------
        !
        CALL this%mapping%update(this%system_system%pos)
        ! update mapping with correct shift of environment cell
        !
        !--------------------------------------------------------------------------------
        ! Second step of initialization for environment ions
        !
        CALL this%environment_ions%init_second(nat, ntyp, ityp, zv, &
                                               this%environment_cell)
        !
        !--------------------------------------------------------------------------------
        ! Update environment ions parameters
        !
        CALL this%environment_ions%update(nat, tau)
        !
        CALL this%environment_ions%printout()
        !
        !--------------------------------------------------------------------------------
        ! Update environment system parameters
        !
        this%environment_system%lupdate = .TRUE.
        !
        CALL this%environment_system%update()
        !
        CALL this%environment_system%printout()
        !
        !--------------------------------------------------------------------------------
        ! Set soft-sphere parameters
        !
        IF (this%lsolvent) CALL this%solvent%set_soft_spheres()
        !
        IF (this%lelectrolyte) CALL this%electrolyte%boundary%set_soft_spheres()
        !
        !--------------------------------------------------------------------------------
        ! Update cores
        !
        IF (this%l1da) &
            CALL this%core_1da_elect%update_origin(this%environment_system%pos)
        !
        !--------------------------------------------------------------------------------
        ! Update rigid environ properties, defined on ions
        !
        IF (this%lrigidcavity) THEN
            !
            IF (this%lsolvent) THEN
                !
                CALL this%solvent%update()
                !
                IF (this%solvent%update_status == 2) CALL this%solvent%printout()
                !
                !------------------------------------------------------------------------
                ! Update quantities that depend on the solvent boundary
                !
                IF (this%lstatic) THEN
                    !
                    CALL this%static%update()
                    !
                    IF (.NOT. this%static%lupdate) CALL this%static%printout()
                    !
                END IF
                !
                IF (this%loptical) THEN
                    !
                    CALL this%optical%update()
                    !
                    IF (.NOT. this%optical%lupdate) CALL this%optical%printout()
                    !
                END IF
                !
            END IF
            !
            IF (this%lelectrolyte) THEN
                !
                CALL this%electrolyte%boundary%update()
                !
                IF (this%electrolyte%boundary%update_status == 2) &
                    CALL this%electrolyte%boundary%printout()
                !
                CALL this%electrolyte%update()
                !
                IF (.NOT. this%electrolyte%lupdate) CALL this%electrolyte%printout()
                !
            END IF
            !
        END IF
        !
        !--------------------------------------------------------------------------------
        ! External charges rely on the environment cell, which is defined
        ! with respect to the system origin
        !
        IF (this%lexternals) CALL this%externals%update()
        !
        IF (this%lelectrostatic .OR. this%lconfine) THEN
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
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE environ_initions
    !------------------------------------------------------------------------------------
    !>
    !! Initialize the electrons-related quantities to be used in the Environ
    !! modules. This initialization is called by electrons.f90, thus it
    !! is performed at every step of electronic optimization.
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE environ_initelectrons(this, nnr, rho, nelec)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: nnr
        REAL(DP), INTENT(IN) :: rho(nnr)
        REAL(DP), INTENT(IN), OPTIONAL :: nelec
        !
        CLASS(environ_obj), INTENT(INOUT) :: this
        !
        REAL(DP), ALLOCATABLE :: aux(:)
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
        ! CALL this%system_electrons%density%write_cube() ! DEBUGGING
        !
        IF (this%ldoublecell) THEN
            ALLOCATE (aux(this%environment_cell%nnr))
            !
            CALL this%mapping%to_large(nnr, this%environment_cell%nnr, rho, aux)
            !
            CALL this%environment_electrons%update(this%environment_cell%nnr, aux, nelec)
            !
        ELSE
            CALL this%environment_electrons%update(nnr, rho, nelec)
        END IF
        !
        this%environment_electrons%density%label = 'large_electrons'
        !
        ! CALL this%environment_electrons%density%write_cube() ! DEBUGGING
        !
        ! STOP ! DEBUGGING
        !
        CALL this%system_electrons%printout()
        !
        CALL this%environment_electrons%printout()
        !
        IF (this%lelectrostatic .OR. this%lconfine) THEN
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
        IF (this%lsoftcavity) THEN
            !
            IF (this%lsoftsolvent) THEN
                !
                CALL this%solvent%update()
                !
                IF (this%solvent%update_status == 2) CALL this%solvent%printout()
                !
                !------------------------------------------------------------------------
                ! Update quantities that depend on the solvent boundary
                !
                IF (this%lstatic) THEN
                    !
                    CALL this%static%update()
                    !
                    IF (.NOT. this%static%lupdate) CALL this%static%printout()
                    !
                END IF
                !
                IF (this%loptical) THEN
                    !
                    CALL this%optical%update()
                    !
                    IF (.NOT. this%optical%lupdate) CALL this%optical%printout()
                    !
                END IF
                !
            END IF
            !
            IF (this%lsoftelectrolyte) THEN
                !
                CALL this%electrolyte%boundary%update()
                !
                IF (this%electrolyte%boundary%update_status == 2) &
                    CALL this%electrolyte%boundary%printout()
                !
                CALL this%electrolyte%update()
                !
                IF (.NOT. this%electrolyte%lupdate) CALL this%electrolyte%printout()
                !
            END IF
            !
        END IF
        !
        this%system_electrons%lupdate = .FALSE.
        this%environment_electrons%lupdate = .FALSE.
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE environ_initelectrons
    !------------------------------------------------------------------------------------
    !>
    !! Initialize the response charges to be used in the TDDFPT + Environ
    !! modules. This initialization is called by plugin_tddfpt_potential.f90
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE environ_initresponse(this, nnr, drho)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: nnr
        REAL(DP), INTENT(IN) :: drho(nnr)
        !
        CLASS(environ_obj), INTENT(INOUT) :: this
        !
        REAL(DP), ALLOCATABLE :: aux(:)
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
        IF (this%ldoublecell) THEN
            !
            ALLOCATE (aux(this%environment_cell%nnr))
            !
            CALL this%mapping%to_large(nnr, this%environment_cell%nnr, drho, aux)
            !
            CALL this%environment_response_electrons%update(this%environment_cell%nnr, &
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
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE environ_initresponse
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !                                  GENERAL METHODS
    !
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !>
    !! Calculates the Environ contribution to the local potential. All
    !! the Environ modules need to be called here. The potentials are
    !! all computed on the dense real-space grid and added to vtot.
    !
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_venviron(this, update, local_verbose)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        LOGICAL, INTENT(IN) :: update
        INTEGER, INTENT(IN), OPTIONAL :: local_verbose
        !
        CLASS(environ_obj), INTENT(INOUT) :: this
        !
        TYPE(environ_density) :: aux
        TYPE(environ_density) :: de_dboundary
        !
        CHARACTER(LEN=80) :: sub_name = 'calc_venviron'
        !
        !--------------------------------------------------------------------------------
        ! If not updating the potentials, add old potentials and exit
        !
        IF (.NOT. update) THEN
            !
            IF (PRESENT(local_verbose)) THEN
                !
                CALL write_cube(this%dvtot, this%system_ions, local_verbose)
                !
                IF (this%lelectrostatic) THEN
                    !
                    CALL write_cube(this%vreference, this%system_ions, local_verbose)
                    !
                    CALL write_cube(this%velectrostatic, this%system_ions, local_verbose)
                    !
                    CALL this%system_charges%printout(local_verbose, 0)
                    !
                END IF
                !
                IF (this%lconfine) &
                    CALL write_cube(this%vconfine, this%system_ions, local_verbose)
                !
                IF (this%lsoftcavity) &
                    CALL write_cube(this%vsoftcavity, this%system_ions, local_verbose)
                !
            END IF
            !
            RETURN
            !
        END IF
        !
        this%dvtot%of_r = 0.D0
        !
        CALL aux%init(this%system_cell)
        !
        !--------------------------------------------------------------------------------
        ! If any form of electrostatic embedding is present, calculate its contribution
        !
        IF (this%lelectrostatic) THEN
            !
            !----------------------------------------------------------------------------
            ! Electrostatics is also computed inside the calling program,
            ! need to remove the reference #TODO to-be-decided
            !
            CALL this%reference%calc_v(this%system_charges, this%vreference)
            !
            CALL write_cube(this%vreference, this%system_ions)
            !
            CALL this%outer%calc_v(this%environment_charges, this%velectrostatic)
            !
            ! IF (this%lexternals) CALL this%environment_charges%update(this%lexternals)
            ! #TODO keep for now until external tests are fully debugged
            !
            CALL write_cube(this%velectrostatic, this%system_ions)
            !
            CALL this%mapping%to_small(this%velectrostatic, aux)
            !
            this%dvtot%of_r = aux%of_r - this%vreference%of_r
            !
            CALL this%environment_charges%of_potential(this%velectrostatic)
            !
            CALL this%environment_charges%printout()
            !
        END IF
        !
        !--------------------------------------------------------------------------------
        ! #TODO add brief description of confine
        !
        IF (this%lconfine) THEN
            !
            CALL this%solvent%vconfine(this%env_confine, this%vconfine)
            !
            CALL write_cube(this%vconfine, this%system_ions)
            !
            CALL this%mapping%to_small(this%vconfine, aux)
            !
            this%dvtot%of_r = this%dvtot%of_r + aux%of_r
        END IF
        !
        !--------------------------------------------------------------------------------
        ! Compute the total potential depending on the boundary
        !
        IF (this%lsoftcavity) THEN
            this%vsoftcavity%of_r = 0.D0
            !
            CALL de_dboundary%init(this%environment_cell)
            !
            IF (this%lsoftsolvent) THEN
                de_dboundary%of_r = 0.D0
                !
                ! if surface tension greater than zero, calculate cavity contribution
                IF (this%lsurface) &
                    CALL this%solvent%desurface_dboundary(this%env_surface_tension, &
                                                          de_dboundary)
                !
                ! if external pressure different from zero, calculate PV contribution
                IF (this%lvolume) &
                    CALL this%solvent%devolume_dboundary(this%env_pressure, de_dboundary)
                !
                ! if confinement potential not zero, calculate confine contribution
                IF (this%lconfine) &
                    CALL this%solvent%deconfine_dboundary(this%env_confine, &
                                                          this%environment_charges%electrons%density, &
                                                          de_dboundary)
                !
                ! if dielectric embedding, calculate dielectric contribution
                IF (this%lstatic) &
                    CALL this%static%de_dboundary(this%velectrostatic, de_dboundary)
                !
                ! if solvent-aware interface correct the potential
                IF (this%solvent%solvent_aware) &
                    CALL this%solvent%sa_de_dboundary(de_dboundary)
                !
                IF (this%solvent%field_aware) THEN
                    CALL env_errore(sub_name, 'field-aware not yet implimented ', 1)
                ELSE
                    !
                    this%vsoftcavity%of_r = de_dboundary%of_r * this%solvent%dscaled%of_r
                    ! multiply by derivative of the boundary w.r.t electronic density
                    !
                END IF
                !
            END IF
            !
            IF (this%lsoftelectrolyte) THEN
                de_dboundary%of_r = 0.D0
                !
                CALL this%electrolyte%de_dboundary(de_dboundary)
                ! if electrolyte is present add its non-electrostatic contribution
                !
                ! if solvent-aware interface correct the potential
                IF (this%electrolyte%boundary%solvent_aware) &
                    CALL this%electrolyte%boundary%sa_de_dboundary(de_dboundary)
                !
                IF (this%electrolyte%boundary%field_aware) THEN
                    CALL env_errore(sub_name, 'field-aware not yet implimented ', 1)
                ELSE
                    !
                    ! multiply for the derivative of the boundary w.r.t electronic density
                    this%vsoftcavity%of_r = &
                        this%vsoftcavity%of_r + &
                        de_dboundary%of_r * this%electrolyte%boundary%dscaled%of_r
                    !
                END IF
                !
            END IF
            !
            CALL write_cube(this%vsoftcavity, this%system_ions)
            !
            CALL this%mapping%to_small(this%vsoftcavity, aux)
            !
            this%dvtot%of_r = this%dvtot%of_r + aux%of_r
            !
            CALL de_dboundary%destroy()
            !
        END IF
        !
        CALL aux%destroy()
        !
        CALL write_cube(this%dvtot, this%system_ions, local_verbose)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_venviron
    !------------------------------------------------------------------------------------
    !>
    !! Calculates the Environ contribution to the energy. We must remove
    !! int v_environ * rhoelec that is automatically included in the
    !! energy computed as the sum of Kohn-Sham eigenvalues.
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_eenviron(this, total_energy)
        !--------------------------------------------------------------------------------
        !
        ! USE embedding_confine, ONLY: calc_econfine #TODO to-be-decided
        !
        IMPLICIT NONE
        !
        CLASS(environ_obj), INTENT(INOUT) :: this
        REAL(DP), INTENT(INOUT) :: total_energy
        !
        REAL(DP) :: ereference
        !
        !--------------------------------------------------------------------------------
        !
        this%eelectrostatic = 0.D0
        this%esurface = 0.D0
        this%evolume = 0.D0
        this%econfine = 0.D0
        this%eelectrolyte = 0.D0
        !
        this%niter = this%niter + 1
        !
        !--------------------------------------------------------------------------------
        ! If electrostatic is on, compute electrostatic energy
        !
        IF (this%lelectrostatic) THEN
            !
            CALL this%reference%calc_e(this%system_charges, this%vreference, ereference)
            !
            CALL this%outer%calc_e(this%environment_charges, this%velectrostatic, &
                                   this%eelectrostatic)
            !
            this%eelectrostatic = this%eelectrostatic - ereference
        END IF
        !
        !--------------------------------------------------------------------------------
        !
        ! if surface tension not zero, compute cavitation energy
        IF (this%lsurface) &
            CALL this%solvent%esurface(this%env_surface_tension, this%esurface)
        !
        IF (this%lvolume) CALL this%solvent%evolume(this%env_pressure, this%evolume)
        ! if pressure not zero, compute PV energy
        !
        ! if confinement potential not zero compute confine energy
        IF (this%lconfine) &
            this%econfine = this%environment_electrons%density%scalar_product(this%vconfine)
        !
        ! if electrolyte is present, calculate its non-electrostatic contribution
        IF (this%lelectrolyte) CALL this%electrolyte%energy(this%eelectrolyte)
        !
        total_energy = total_energy + this%eelectrostatic + this%esurface + &
                       this%evolume + this%econfine + this%eelectrolyte
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_eenviron
    !------------------------------------------------------------------------------------
    !>
    !! Calculates the Environ contribution to the forces. Due to
    !! Hellman-Feynman only a few of the Environ modules have an
    !! effect on the atomic forces.
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_fenviron(this, nat, force_environ)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: nat
        !
        CLASS(environ_obj), INTENT(INOUT) :: this
        REAL(DP), INTENT(INOUT) :: force_environ(3, nat)
        !
        INTEGER :: i
        TYPE(environ_density) :: de_dboundary
        TYPE(environ_gradient) :: partial
        !
        !--------------------------------------------------------------------------------
        !
        force_environ = 0.D0
        !
        ! compute the electrostatic embedding contribution to the interatomic forces
        IF (this%lelectrostatic) &
            CALL this%outer%calc_f(nat, this%environment_charges, force_environ)
        !
        !--------------------------------------------------------------------------------
        ! Compute the total forces depending on the boundary
        !
        IF (this%lrigidcavity) THEN
            !
            CALL de_dboundary%init(this%environment_cell)
            !
            CALL partial%init(this%environment_cell)
            !
            IF (this%lrigidsolvent) THEN
                !
                de_dboundary%of_r = 0.D0
                !
                ! if surface tension greater than zero, calculate cavity contribution
                IF (this%lsurface) &
                    CALL this%solvent%desurface_dboundary(this%env_surface_tension, &
                                                          de_dboundary)
                !
                ! if external pressure not zero, calculate PV contribution
                IF (this%lvolume) &
                    CALL this%solvent%devolume_dboundary(this%env_pressure, de_dboundary)
                !
                ! if confinement potential not zero, calculate confine contribution
                IF (this%lconfine) &
                    CALL this%solvent%deconfine_dboundary(this%env_confine, &
                                                          this%environment_charges%electrons%density, &
                                                          de_dboundary)
                !
                ! if dielectric embedding, calculate dielectric contribution
                IF (this%lstatic) &
                    CALL this%static%de_dboundary(this%velectrostatic, de_dboundary)
                !
                ! if solvent-aware, correct the potential
                IF (this%solvent%solvent_aware) &
                    CALL this%solvent%sa_de_dboundary(de_dboundary)
                !
                !------------------------------------------------------------------------
                ! Multiply by derivative of the boundary w.r.t ionic positions
                !
                DO i = 1, nat
                    !
                    CALL this%solvent%dboundary_dions(i, partial)
                    !
                    force_environ(:, i) = force_environ(:, i) - &
                                          partial%scalar_product_density(de_dboundary)
                    !
                END DO
                !
            END IF
            !
            IF (this%lrigidelectrolyte) THEN
                !
                de_dboundary%of_r = 0.D0
                !
                ! if electrolyte is present, add its non-electrostatic contribution
                CALL this%electrolyte%de_dboundary(de_dboundary)
                !
                ! if solvent-aware, correct the potential
                IF (this%electrolyte%boundary%solvent_aware) &
                    CALL this%electrolyte%boundary%sa_de_dboundary(de_dboundary)
                !
                !------------------------------------------------------------------------
                ! Multiply by derivative of the boundary w.r.t ionic positions
                !
                DO i = 1, nat
                    !
                    CALL this%electrolyte%boundary%dboundary_dions(i, partial)
                    !
                    force_environ(:, i) = force_environ(:, i) - &
                                          partial%scalar_product_density(de_dboundary)
                    !
                END DO
                !
            END IF
            !
            CALL partial%destroy()
            !
            CALL de_dboundary%destroy()
            !
        END IF
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_fenviron
    !------------------------------------------------------------------------------------
    !>
    !! Calculates the Environ contribution to the response potential in TD calculations
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_dvenviron(this, nnr, dv)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: nnr
        !
        CLASS(environ_obj), INTENT(INOUT) :: this
        REAL(DP), INTENT(INOUT) :: dv(nnr)
        !
        TYPE(environ_density) :: aux
        TYPE(environ_density) :: dvreference
        TYPE(environ_density) :: dvelectrostatic
        TYPE(environ_density) :: dvsoftcavity
        TYPE(environ_density) :: dv_dboundary
        !
        !--------------------------------------------------------------------------------
        !
        CALL aux%init(this%system_cell)
        !
        !--------------------------------------------------------------------------------
        ! If any form of electrostatic embedding is present, calculate its contribution
        !
        IF (this%lelectrostatic) THEN
            !
            !----------------------------------------------------------------------------
            ! Electrostatics is also computed inside the calling program,
            ! need to remove the reference
            !
            CALL dvreference%init(this%system_cell)
            !
            CALL this%reference%calc_v(this%system_response_charges, dvreference)
            !
            CALL dvelectrostatic%init(this%environment_cell)
            !
            CALL this%outer%calc_v(this%environment_response_charges, dvelectrostatic)
            !
            CALL this%mapping%to_small(dvelectrostatic, aux)
            !
            dv(:) = dv(:) + aux%of_r(:) - dvreference%of_r(:)
            !
            CALL dvreference%destroy()
            !
        END IF
        !
        !--------------------------------------------------------------------------------
        ! Compute the response potential depending on the boundary
        !
        IF (this%lsoftcavity) THEN
            !
            CALL dvsoftcavity%init(this%environment_cell)
            !
            CALL dv_dboundary%init(this%environment_cell)
            !
            IF (this%lsoftsolvent) THEN
                dv_dboundary%of_r = 0.D0
                !
                ! if dielectric embedding, calcultes dielectric contribution
                IF (this%loptical) &
                    CALL this%optical%dv_dboundary(this%velectrostatic, &
                                                   dvelectrostatic, dv_dboundary)
                !
                dvsoftcavity%of_r = dv_dboundary%of_r * this%solvent%dscaled%of_r
            END IF
            !
            CALL this%mapping%to_small(dvsoftcavity, aux)
            !
            dv(:) = dv(:) + aux%of_r(:)
            !
            CALL dv_dboundary%destroy()
            !
            CALL dvsoftcavity%destroy()
            !
        END IF
        !
        IF (this%lelectrostatic) CALL dvelectrostatic%destroy()
        !
        CALL aux%destroy()
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_dvenviron
    !------------------------------------------------------------------------------------
    !>
    !! Calculates the energy corrections in PW calculations
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_deenviron(this, total_energy)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_obj), INTENT(INOUT) :: this
        REAL(DP), INTENT(INOUT) :: total_energy
        !
        !--------------------------------------------------------------------------------
        !
        this%deenviron = -this%system_electrons%density%scalar_product(this%dvtot)
        total_energy = total_energy + this%deenviron
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_deenviron
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
    SUBROUTINE set_electrostatic_base(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_obj), TARGET, INTENT(INOUT) :: this
        !
        CLASS(electrostatic_solver), POINTER :: local_outer_solver, local_inner_solver
        !
        LOGICAL :: lconjugate, lnested
        !
        CHARACTER(LEN=80) :: local_type, local_auxiliary, local_problem, inner_problem
        !
        CHARACTER(LEN=20) :: sub_name = 'set_electrostatic_base'
        !
        !--------------------------------------------------------------------------------
        ! Initial setup of core flags
        !
        this%lfd = .FALSE.
        this%l1da = .FALSE.
        this%lfft_system = .FALSE.
        this%lfft_environment = .FALSE.
        !
        !--------------------------------------------------------------------------------
        ! Setup nested scheme if required
        !
        IF (inner_solver /= 'none') THEN
            lnested = .TRUE.
        ELSE
            lnested = .FALSE.
        END IF
        !
        !--------------------------------------------------------------------------------
        ! Set reference core according to calling program
        !
        SELECT CASE (prog)
            !
        CASE ('PW', 'CP', 'TD', 'XS')
            this%lfft_system = .TRUE.
            local_type = 'fft'
            !
            CALL this%reference_cores%init(this%core_fft_sys, local_type)
            !
        CASE DEFAULT
            CALL env_errore(sub_name, 'Unexpected name of host code', 1)
            !
        END SELECT
        !
        !--------------------------------------------------------------------------------
        ! Numerical core for periodic boundary or electrolyte corrections
        !
        this%need_pbc_correction = .FALSE.
        this%need_electrolyte = .FALSE.
        this%need_semiconductor = .FALSE.
        this%need_outer_loop = .FALSE.
        !
        !--------------------------------------------------------------------------------
        ! First check keywords specified in input
        !
        IF (pbc_dim >= 0) THEN
            !
            SELECT CASE (TRIM(ADJUSTL(pbc_correction)))
                !
            CASE ('none')
                this%need_pbc_correction = .FALSE.
                !
            CASE ('parabolic')
                this%need_pbc_correction = .TRUE.
                this%l1da = .TRUE.
                local_type = '1da'
                !
            CASE ('gcs', 'gouy-chapman', 'gouy-chapman-stern')
                this%need_pbc_correction = .TRUE.
                this%need_electrolyte = .TRUE.
                this%l1da = .TRUE.
                local_type = 'gcs'
                !
            CASE ('ms', 'mott-schottky')
                this%need_pbc_correction = .TRUE.
                this%need_semiconductor = .TRUE.
                this%l1da = .TRUE.
                local_type = 'ms'
                !
            CASE ('ms-gcs', 'mott-schottky-gouy-chapman-stern')
                this%need_pbc_correction = .TRUE.
                this%need_semiconductor = .TRUE.
                this%need_outer_loop = .TRUE.
                this%need_electrolyte = .TRUE.
                this%l1da = .TRUE.
                local_type = 'ms-gcs'
                !
            CASE DEFAULT
                CALL env_errore(sub_name, 'Option not yet implemented', 1)
                !
            END SELECT
            !
        END IF
        !
        IF (this%need_pbc_correction .AND. this%l1da) &
            CALL this%pbc_core%init(this%core_1da_elect, local_type)
        !
        !--------------------------------------------------------------------------------
        ! Set up main (outer) core and inner core (if nested scheme)
        !
        SELECT CASE (core)
            !
        CASE ('fft')
            this%lfft_environment = .TRUE.
            local_type = 'fft'
            !
            CALL this%outer_cores%init(this%core_fft_elect, local_type)
            !
            IF (lnested) CALL this%inner_cores%init(this%core_fft_elect, local_type)
            !
        CASE ('1d-analytic', '1da')
            this%l1da = .TRUE.
            local_type = '1d-analytic'
            !
            CALL this%outer_cores%init(this%core_1da_elect, local_type)
            !
            IF (lnested) CALL this%inner_cores%init(this%core_1da_elect, local_type)
            !
        CASE DEFAULT
            !
            CALL env_errore(sub_name, &
                            'Unexpected value for core container type keyword', 1)
            !
        END SELECT
        !
        !--------------------------------------------------------------------------------
        ! Setup corrections if necessary
        !
        IF (this%need_pbc_correction) CALL this%outer_cores%add_correction(this%pbc_core)
        !
        IF (lnested .AND. this%need_pbc_correction) &
            CALL this%inner_cores%add_correction(this%pbc_core)
        !
        !--------------------------------------------------------------------------------
        ! Set reference solver according to calling program
        !
        SELECT CASE (prog)
            !
        CASE ('PW', 'CP', 'TD', 'XS')
            CALL this%reference_direct%set_core(this%reference_cores)
            !
        CASE DEFAULT
            CALL env_errore(sub_name, 'Unexpected name of host code', 1)
            !
        END SELECT
        !
        !--------------------------------------------------------------------------------
        ! Set up main (outer) solver
        !
        lconjugate = .FALSE.
        !
        SELECT CASE (solver)
            !
        CASE ('direct')
            CALL this%direct%set_core(this%outer_cores)
            !
            local_outer_solver => this%direct
            !
        CASE ('cg', 'sd')
            !
            IF (TRIM(ADJUSTL(solver)) == 'cg') lconjugate = .TRUE.
            !
            CALL this%gradient%create(this%outer_cores, maxstep, tol, auxiliary)
            !
            CALL this%gradient%init(solver, lconjugate, step_type, step, &
                                    preconditioner, screening_type, screening)
            !
            local_outer_solver => this%gradient
            !
        CASE ('fp')
            !
            CALL this%fixedpoint%create(this%outer_cores, maxstep, tol, auxiliary)
            !
            CALL this%fixedpoint%init(mix_type, mix, ndiis)
            !
            local_outer_solver => this%fixedpoint
            !
        CASE ('newton')
            CALL this%newton%create(this%outer_cores, maxstep, tol, auxiliary)
            !
            CALL this%newton%init()
            !
            local_outer_solver => this%newton
            !
        CASE DEFAULT
            !
            CALL env_errore(sub_name, &
                            'Unexpected value for electrostatic solver keyword', 1)
            !
        END SELECT
        !
        !--------------------------------------------------------------------------------
        ! If nested scheme, set up inner solver
        !
        IF (lnested) THEN
            lconjugate = .FALSE.
            !
            SELECT CASE (solver)
                !
            CASE ('fp')
                !
                IF (auxiliary == 'ioncc') THEN
                    inner_problem = 'generalized'
                    !
                    SELECT CASE (inner_solver)
                        !
                    CASE ('cg', 'sd')
                        !
                        IF (TRIM(ADJUSTL(inner_solver)) == 'cg') lconjugate = .TRUE.
                        !
                        CALL this%inner_gradient%create(this%inner_cores, &
                                                        inner_maxstep, inner_tol, &
                                                        auxiliary)
                        !
                        CALL this%inner_gradient%init(inner_solver, lconjugate, &
                                                      step_type, step, preconditioner, &
                                                      screening_type, screening)
                        !
                        local_inner_solver => this%inner_gradient
                        !
                    CASE ('fp')
                        local_auxiliary = 'full'
                        !
                        CALL this%inner_fixedpoint%create(this%inner_cores, &
                                                          inner_maxstep, inner_tol, &
                                                          local_auxiliary)
                        !
                        CALL this%inner_fixedpoint%init(mix_type, inner_mix, ndiis)
                        !
                        local_inner_solver => this%inner_fixedpoint
                        !
                    END SELECT
                    !
                ELSE
                    !
                    CALL env_errore(sub_name, &
                                    'Unexpected value for auxiliary charge &
                                    &in nested solver', 1)
                    !
                END IF
                !
            CASE ('newton')
                inner_problem = 'linpb'
                lconjugate = .TRUE.
                !
                CALL this%inner_gradient%create(this%inner_cores, inner_maxstep, &
                                                inner_tol, auxiliary)
                !
                CALL this%inner_gradient%init(inner_solver, lconjugate, step_type, &
                                              step, preconditioner, screening_type, &
                                              screening)
                !
                local_inner_solver => this%inner_gradient
                !
            END SELECT
            !
        END IF
        !
        !--------------------------------------------------------------------------------
        ! Set reference setup according to calling program
        !
        CALL this%reference%create()
        !
        SELECT CASE (prog)
            !
        CASE ('PW', 'CP', 'TD', 'XS')
            local_problem = 'poisson'
            !
        CASE DEFAULT
            CALL env_errore(sub_name, 'Unexpected name of host code', 1)
            !
        END SELECT
        !
        CALL this%reference%init(local_problem, this%reference_direct)
        !
        !--------------------------------------------------------------------------------
        ! Create outer setup
        !
        CALL this%outer%create()
        !
        CALL this%outer%init(problem, local_outer_solver)
        !
        !--------------------------------------------------------------------------------
        ! If nested scheme, create inner setup
        !
        IF (lnested) THEN
            !
            CALL this%inner%create()
            !
            CALL this%inner%init(inner_problem, local_inner_solver)
            !
            this%outer%inner => this%inner
        END IF
        !
        !--------------------------------------------------------------------------------
        ! Set logical flags according to electrostatic setup
        !
        this%need_gradient = .FALSE.
        this%need_factsqrt = .FALSE.
        this%need_auxiliary = .FALSE.
        !
        CALL this%reference%set_flags(this%need_auxiliary, this%need_gradient, &
                                      this%need_factsqrt)
        !
        CALL this%outer%set_flags(this%need_auxiliary, this%need_gradient, &
                                  this%need_factsqrt)
        !
        IF (lnested) &
            CALL this%inner%set_flags(this%need_auxiliary, this%need_gradient, &
                                      this%need_factsqrt)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE set_electrostatic_base
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE set_environ_base(this, nelec, nat, ntyp, atom_label)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: nelec, nat, ntyp
        CHARACTER(LEN=3), INTENT(IN) :: atom_label(:)
        !
        CLASS(environ_obj), INTENT(INOUT) :: this
        !
        INTEGER :: i
        !
        CHARACTER(LEN=80) :: label
        !
        CHARACTER(LEN=20) :: sub_name = 'set_environ_base'
        !
        !--------------------------------------------------------------------------------
        ! TDDFPT flag
        !
        SELECT CASE (prog)
            !
        CASE ('TD')
            this%ltddfpt = .TRUE.
            !
        CASE DEFAULT
            this%ltddfpt = .FALSE.
            !
        END SELECT
        !
        !--------------------------------------------------------------------------------
        ! Create necessary local types
        !
        CALL this%system_electrons%create()
        !
        CALL this%system_ions%create()
        !
        CALL this%system_system%create()
        !
        CALL this%environment_electrons%create()
        !
        CALL this%environment_ions%create()
        !
        CALL this%environment_system%create()
        !
        !--------------------------------------------------------------------------------
        ! General flags
        !
        this%environ_restart = environ_restart
        this%environ_thr = environ_thr
        this%environ_nskip = environ_nskip
        !
        !--------------------------------------------------------------------------------
        ! Set main environment flags, convert to internal units
        !
        this%env_static_permittivity = env_static_permittivity
        this%env_optical_permittivity = env_optical_permittivity
        !
        this%env_surface_tension = &
            env_surface_tension * 1.D-3 * BOHR_RADIUS_SI**2 / RYDBERG_SI
        !
        this%env_pressure = env_pressure * 1.D9 / RYDBERG_SI * BOHR_RADIUS_SI**3
        this%env_confine = env_confine
        this%env_electrolyte_ntyp = env_electrolyte_ntyp
        !
        !--------------------------------------------------------------------------------
        ! Set basic logical flags
        !
        this%lstatic = env_static_permittivity > 1.D0
        this%loptical = env_optical_permittivity > 1.D0
        !
        IF (env_dielectric_regions > 0) THEN
            !
            DO i = 1, env_dielectric_regions
                this%lstatic = this%lstatic .OR. (epsregion_eps(1, i) > 1.D0)
                this%loptical = this%loptical .OR. (epsregion_eps(2, i) > 1.D0)
            END DO
            !
        END IF
        !
        this%lsurface = this%env_surface_tension > 0.D0
        this%lvolume = this%env_pressure /= 0.D0
        this%lconfine = this%env_confine /= 0.D0
        this%lexternals = env_external_charges > 0
        this%lelectrolyte = this%env_electrolyte_ntyp > 0 .OR. this%need_electrolyte
        this%lsemiconductor = this%need_semiconductor
        this%lperiodic = this%need_pbc_correction
        this%ldoublecell = SUM(env_nrep) > 0
        this%louterloop = this%need_outer_loop
        !
        !--------------------------------------------------------------------------------
        ! Derived flags
        !
        this%ldielectric = this%lstatic .OR. this%loptical
        !
        this%lsolvent = this%ldielectric .OR. this%lsurface .OR. this%lvolume .OR. &
                        this%lconfine
        !
        this%lelectrostatic = this%ldielectric .OR. this%lelectrolyte .OR. &
                              this%lexternals .OR. this%lperiodic
        !
        this%lsoftsolvent = this%lsolvent .AND. (solvent_mode == 'electronic' .OR. &
                                                 solvent_mode == 'full' .OR. &
                                                 solvent_mode(1:2) == 'fa')
        !
        this%lsoftelectrolyte = this%lelectrolyte .AND. &
                                (electrolyte_mode == 'electronic' .OR. &
                                 electrolyte_mode == 'full' .OR. &
                                 electrolyte_mode(1:2) == 'fa') ! field-aware
        !
        this%lsoftcavity = this%lsoftsolvent .OR. this%lsoftelectrolyte
        this%lrigidsolvent = this%lsolvent .AND. solvent_mode /= 'electronic'
        this%lrigidelectrolyte = this%lelectrolyte .AND. electrolyte_mode /= 'electronic'
        this%lrigidcavity = this%lrigidsolvent .OR. this%lrigidelectrolyte
        !
        this%lcoredensity = (this%lsolvent .AND. solvent_mode == 'full') .OR. &
                            (this%lelectrolyte .AND. electrolyte_mode == 'full')
        !
        this%lsmearedions = this%lelectrostatic
        this%lboundary = this%lsolvent .OR. this%lelectrolyte
        this%lgradient = this%ldielectric .OR. (solvent_mode(1:2) == 'fa') ! field-aware
        !
        !--------------------------------------------------------------------------------
        ! Create optional types
        !
        IF (this%lexternals) CALL this%externals%create()
        !
        IF (this%lsolvent) THEN
            label = 'solvent'
            !
            CALL this%solvent%create(label)
            !
        END IF
        !
        IF (this%lelectrolyte) CALL this%electrolyte%create()
        !
        IF (this%lsemiconductor) CALL this%semiconductor%create()
        !
        IF (this%lstatic) CALL this%static%create()
        !
        !--------------------------------------------------------------------------------
        ! Set response properties for TD calculations
        !
        IF (this%loptical) THEN
            !
            CALL this%optical%create()
            !
            CALL this%system_response_electrons%create()
            !
            CALL this%system_response_electrons%init_first(0)
            !
            CALL this%system_response_charges%create()
            !
            CALL this%system_response_charges%init_first( &
                electrons=this%system_response_electrons)
            !
            CALL this%environment_response_electrons%create()
            !
            CALL this%environment_response_electrons%init_first(0)
            !
            CALL this%environment_response_charges%create()
            !
            CALL this%environment_response_charges%init_first( &
                electrons=this%environment_response_electrons, &
                dielectric=this%optical)
            !
        END IF
        !
        IF (this%lelectrostatic .OR. this%lconfine) THEN
            !
            CALL this%system_charges%create()
            !
            CALL this%environment_charges%create()
            !
        END IF
        !
        !--------------------------------------------------------------------------------
        ! Allocate and set basic properties of ions
        !
        CALL this%system_ions%init_first(nat, ntyp, this%lsoftcavity, &
                                         this%lcoredensity, this%lsmearedions, &
                                         radius_mode, atom_label, atomicspread, &
                                         corespread, solvationrad)
        !
        CALL this%environment_ions%init_first(nat, ntyp, this%lsoftcavity, &
                                              this%lcoredensity, this%lsmearedions, &
                                              radius_mode, atom_label, atomicspread, &
                                              corespread, solvationrad)
        !
        !--------------------------------------------------------------------------------
        ! Set basic properties of electrons
        !
        CALL this%system_electrons%init_first(nelec)
        !
        CALL this%environment_electrons%init_first(nelec)
        !
        !--------------------------------------------------------------------------------
        ! Set basic properties of the selected system
        !
        CALL this%system_system%init(system_ntyp, system_dim, system_axis, &
                                     this%system_ions)
        !
        CALL this%environment_system%init(system_ntyp, system_dim, system_axis, &
                                          this%environment_ions)
        !
        !--------------------------------------------------------------------------------
        ! Collect free charges if computing electrostatics or confinement
        !
        IF (this%lelectrostatic .OR. this%lconfine) THEN
            !
            CALL this%system_charges%init_first(this%system_electrons)
            !
            CALL this%environment_charges%init_first(this%environment_electrons)
            !
        END IF
        !
        IF (this%lelectrostatic) THEN
            !
            CALL this%system_charges%init_first(ions=this%system_ions)
            !
            CALL this%environment_charges%init_first(ions=this%environment_ions)
            !
        END IF
        !
        !--------------------------------------------------------------------------------
        ! Allocate and set basic properties of external charges
        !
        IF (this%lexternals) THEN
            !
            CALL this%externals%init_first(env_external_charges, extcharge_dim, &
                                           extcharge_axis, extcharge_pos, &
                                           extcharge_spread, extcharge_charge)
            !
            CALL this%environment_charges%init_first(externals=this%externals)
            !
            DEALLOCATE (extcharge_axis)
            DEALLOCATE (extcharge_dim)
            DEALLOCATE (extcharge_spread)
            DEALLOCATE (extcharge_charge)
        END IF
        !
        !--------------------------------------------------------------------------------
        ! Setup cores needed for derivatives of boundaries
        !
        IF (this%lboundary) THEN
            this%lfft_environment = .TRUE.
            !
            IF (derivatives == 'fd') this%lfd = .TRUE.
            !
            CALL this%derivatives%init(this%core_fft_deriv, derivatives)
            !
        END IF
        !
        !--------------------------------------------------------------------------------
        ! Set the parameters of the solvent boundary
        !
        IF (this%lsolvent) THEN
            !
            CALL this%solvent%init_first( &
                this%lgradient, this%need_factsqrt, this%lsurface, solvent_mode, stype, &
                rhomax, rhomin, tbeta, this%env_static_permittivity, alpha, softness, &
                solvent_distance, solvent_spread, solvent_radius, radial_scale, &
                radial_spread, filling_threshold, filling_spread, field_awareness, &
                charge_asymmetry, field_max, field_min, this%environment_electrons, &
                this%environment_ions, this%environment_system, this%derivatives)
            !
        END IF
        !
        !--------------------------------------------------------------------------------
        ! Set the parameters of the electrolyte and of its boundary
        !
        IF (this%lelectrolyte) THEN
            !
            CALL this%electrolyte%init_first( &
                env_electrolyte_ntyp, electrolyte_mode, stype, electrolyte_rhomax, &
                electrolyte_rhomin, electrolyte_tbeta, this%env_static_permittivity, &
                electrolyte_alpha, electrolyte_softness, electrolyte_distance, &
                electrolyte_spread, solvent_radius, radial_scale, radial_spread, &
                filling_threshold, filling_spread, field_awareness, charge_asymmetry, &
                field_max, field_min, this%environment_electrons, this%environment_ions, &
                this%environment_system, this%derivatives, temperature, cion, cionmax, &
                rion, zion, electrolyte_entropy, electrolyte_linearized)
            !
            CALL this%environment_charges%init_first(electrolyte=this%electrolyte)
            !
        END IF
        !
        !--------------------------------------------------------------------------------
        ! Set the parameters of the semiconductor
        !
        IF (this%lsemiconductor) THEN
            !
            CALL this%semiconductor%init_first(temperature, sc_permittivity, &
                                               sc_carrier_density, sc_electrode_chg, &
                                               sc_distance, sc_spread, sc_chg_thr, &
                                               this%environment_system)
            !
            CALL this%environment_charges%init_first(semiconductor=this%semiconductor)
            !
        END IF
        !
        !--------------------------------------------------------------------------------
        ! Set the parameters of the dielectric
        !
        IF (this%lstatic) THEN
            !
            CALL this%static%init_first(this%env_static_permittivity, this%solvent, &
                                        this%need_gradient, this%need_factsqrt, &
                                        this%need_auxiliary)
            !
            IF (env_dielectric_regions > 0) &
                CALL this%static%set_regions(env_dielectric_regions, epsregion_dim, &
                                             epsregion_axis, epsregion_pos, &
                                             epsregion_width, epsregion_spread, &
                                             epsregion_eps(1, :))
            !
            CALL this%environment_charges%init_first(dielectric=this%static)
            !
        END IF
        !
        IF (this%loptical) THEN
            !
            CALL this%optical%init_first(this%env_optical_permittivity, this%solvent, &
                                         this%need_gradient, this%need_factsqrt, &
                                         this%need_auxiliary)
            !
            IF (env_dielectric_regions > 0) &
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
        ! Set the parameters for double cell mapping
        !
        IF (this%ldoublecell) ALLOCATE (this%environment_cell)
        !
        CALL this%mapping%init_first(env_nrep) ! #TODO should this be inside the above IF?
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE set_environ_base
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE set_core_base(this, use_internal_pbc_corr)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        LOGICAL, INTENT(IN) :: use_internal_pbc_corr
        !
        CLASS(environ_obj), INTENT(INOUT) :: this
        !
        !--------------------------------------------------------------------------------
        ! Set up active numerical cores
        !
        IF (this%lfd) CALL this%core_fd%init_first(ifdtype, nfdpoint)
        !
        IF (this%lfft_system) CALL this%core_fft_sys%init_first(use_internal_pbc_corr)
        !
        IF (this%lfft_environment) THEN
            !
            CALL this%core_fft_deriv%init_first()
            !
            IF (this%lelectrostatic) &
                CALL this%core_fft_elect%init_first(use_internal_pbc_corr)
            !
        END IF
        !
        IF (this%l1da) CALL this%core_1da_elect%init_first(pbc_dim, pbc_axis)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE set_core_base
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !                                   OUTPUT METHODS
    !
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !!
    !> Write out the main parameters of Environ calculations, summarizing
    !! the input keywords (some info also on internal vs input units).
    !! Called by summary.f90
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE environ_summary(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_obj), INTENT(INOUT) :: this
        !
        !--------------------------------------------------------------------------------
        !
        WRITE (program_unit, *)
        WRITE (program_unit, 1000)
        WRITE (program_unit, 1001) bibliography(1)
        !
        !--------------------------------------------------------------------------------
        ! Environ Summary
        !
        WRITE (program_unit, 1002) this%environ_thr
        !
        IF (this%lsolvent) THEN
            !
            IF (this%solvent%b_type == 0) THEN
                WRITE (program_unit, 1003) 'Fatteber-Gygi'
                WRITE (program_unit, 1004) this%solvent%rhozero, this%solvent%tbeta
            ELSE IF (this%solvent%b_type == 1) THEN
                WRITE (program_unit, 1003) 'SCCS'
                WRITE (program_unit, 1005) this%solvent%rhomax, this%solvent%rhomin
            END IF
            !
            IF (this%solvent%solvent_aware) WRITE (program_unit, 1006)
            !
            IF (this%solvent%field_aware) THEN
                WRITE (program_unit, 1007)
                !
                WRITE (program_unit, 1008) &
                    this%solvent%field_factor, this%solvent%charge_asymmetry
                !
                WRITE (program_unit, 1009) this%solvent%field_min, this%solvent%field_max
                !
            END IF
            !
        END IF
        !
        IF (this%env_static_permittivity > 1.D0) THEN
            WRITE (program_unit, 1010) this%env_static_permittivity
            !
            IF (this%ltddfpt) WRITE (program_unit, 1011) this%env_optical_permittivity
            !
            WRITE (program_unit, 1012) TRIM(this%solvent%mode)
        END IF
        !
        IF (this%env_surface_tension > 0.D0) &
            WRITE (program_unit, 1013) &
            this%env_surface_tension / 1.D-3 / BOHR_RADIUS_SI**2 * RYDBERG_SI, &
            this%env_surface_tension
        !
        IF (this%env_pressure /= 0.D0) &
            WRITE (program_unit, 1014) &
            this%env_pressure * RYDBERG_SI / BOHR_RADIUS_SI**3 * 1.D-9, &
            this%env_pressure
        !
        !------------------------------------------------------------------------
        ! Electrostatic Summary
        !
        IF (this%lelectrostatic) THEN
            WRITE (program_unit, 1015)
            WRITE (program_unit, 1016) this%outer%problem, this%outer%solver%solver_type
            !
            SELECT TYPE (solver => this%outer%solver)
                !
            CLASS IS (solver_iterative)
                WRITE (program_unit, 1017) solver%auxiliary
                !
            END SELECT
            !
            WRITE (program_unit, 1018) this%outer%solver%cores%core%core_type
            !
            IF (ASSOCIATED(this%outer%solver%cores%correction)) THEN
                !
                WRITE (program_unit, 1019) &
                    this%outer%solver%cores%correction%core%core_type
                !
            END IF
            !
            IF (this%lfd) THEN
                !
                IF (this%core_fd%ifdtype == 1) THEN
                    WRITE (program_unit, 1020) 'central diff.', this%core_fd%nfdpoint
                ELSE IF (this%core_fd%ifdtype == 2 .OR. this%core_fd%ifdtype == 3) THEN
                    WRITE (program_unit, 1020) 'lanczos diff.', this%core_fd%nfdpoint
                ELSE IF (this%core_fd%ifdtype == 4 .OR. this%core_fd%ifdtype == 5) THEN
                    WRITE (program_unit, 1020) 'noise-robust diff.', this%core_fd%nfdpoint
                END IF
                !
            END IF
            !
        END IF
        !
        WRITE (program_unit, 1021)
        !
        !--------------------------------------------------------------------------------
        !
1000    FORMAT(/, 5X, 'Environ Module', /, 5X, '==============')
!
1001    FORMAT(/, 5X, 'Please cite', /, 9X, A80, &
                /, 5X, 'in publications or presentations arising from this work.',/)
        !
1002    FORMAT('     compensation onset threshold      = ', E24.4, ' ')
1003    FORMAT('     switching function adopted        = ', A24, ' ')
        !
1004    FORMAT('     solvation density threshold       = ', E24.4, ' ' &
               /'     smoothness exponent (2 x beta)    = ', F24.2, ' ')
        !
1005    FORMAT('     density limit for vacuum region   = ', E24.4, ' ' &
               /'     density limit for bulk solvent    = ', E24.4, ' ')
        !
1006    FORMAT('     interface is solvent aware            ')
1007    FORMAT('     interface is field aware            ')
        !
1008    FORMAT('     field aware factor                = ', F24.2, ' ' &
               /'     asymmetry of field-awareness      = ', F24.2, ' ')
        !
1009    FORMAT('     field limit for no correction     = ', F24.2, ' ' &
               /'     field limit for full correction   = ', F24.2, ' ')
        !
1010    FORMAT('     static permittivity               = ', F24.2, ' ')
1011    FORMAT('     optical permittivity              = ', F24.4, ' ')
1012    FORMAT('     epsilon calculation mode          = ', A24, ' ')
        !
1013    FORMAT('     surface tension in input (dyn/cm) = ', F24.2, ' ' &
               /'     surface tension in internal units = ', E24.4, ' ')
        !
1014    FORMAT('     external pressure in input (GPa)  = ', F24.2, ' ' &
               /'     external pressure in inter. units = ', E24.4, ' ')
        !
1015    FORMAT(/, 5X, 'Electrostatic Setup', /, 5X, '-------------------')
        !
1016    FORMAT('     electrostatic problem to solve    = ', A24, ' ' &
               /'     numerical solver adopted          = ', A24, ' ')
        !
1017    FORMAT('     type of auxiliary density adopted = ', A24, ' ')
1018    FORMAT('     type of core tool for poisson     = ', A24, ' ')
1019    FORMAT('     type of core tool for correction  = ', A24, ' ')
        !
1020    FORMAT('     type of numerical differentiator  = ', A24, ' ' &
               /'     number of points in num. diff.    = ', I24, ' ')
        !
1021    FORMAT(/)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE environ_summary
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE class_environ
!----------------------------------------------------------------------------------------
