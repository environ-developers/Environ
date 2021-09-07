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
! Authors: Edan Bainglass (Department of Physics, UNT)
!
!----------------------------------------------------------------------------------------
!>
!!
!----------------------------------------------------------------------------------------
MODULE env_setup
    !------------------------------------------------------------------------------------
    !
    USE env_base_io, ONLY: program_unit, prog, global_verbose
    !
    USE environ_param, ONLY: DP, BOHR_RADIUS_SI, RYDBERG_SI
    !
    USE env_base_input
    USE environ_input, ONLY: env_read_input
    !
    USE class_cell
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
    TYPE, PUBLIC :: environ_setup
        !--------------------------------------------------------------------------------
        !
        !--------------------------------------------------------------------------------
        ! Main flags
        !
        LOGICAL :: restart = .FALSE.
        REAL(DP) :: threshold = 0.0_DP
        INTEGER :: nskip = 0
        INTEGER :: niter = 0
        INTEGER :: nrep = 1
        !
        !--------------------------------------------------------------------------------
        ! Environment flags
        !
        REAL(DP) :: static_permittivity = 0.0_DP
        REAL(DP) :: optical_permittivity = 0.0_DP
        !
        REAL(DP) :: surface_tension = 0.0_DP
        REAL(DP) :: pressure = 0.0_DP
        REAL(DP) :: confine = 0.0_DP
        !
        !--------------------------------------------------------------------------------
        ! Set basic logical flags = .FALSE.
        !
        LOGICAL :: lstatic = .FALSE.
        LOGICAL :: loptical = .FALSE.
        LOGICAL :: lsurface = .FALSE.
        LOGICAL :: lvolume = .FALSE.
        LOGICAL :: lconfine = .FALSE.
        LOGICAL :: lexternals = .FALSE.
        LOGICAL :: lregions = .FALSE.
        LOGICAL :: lelectrolyte = .FALSE.
        LOGICAL :: lsemiconductor = .FALSE.
        LOGICAL :: lperiodic = .FALSE.
        LOGICAL :: ldoublecell = .FALSE.
        LOGICAL :: louterloop = .FALSE.
        LOGICAL :: ltddfpt = .FALSE.
        !
        !--------------------------------------------------------------------------------
        ! Derived flags
        !
        LOGICAL :: ldielectric = .FALSE.
        LOGICAL :: lsolvent = .FALSE.
        LOGICAL :: lelectrostatic = .FALSE.
        LOGICAL :: lsoftsolvent = .FALSE.
        LOGICAL :: lsoftelectrolyte = .FALSE.
        LOGICAL :: lsoftcavity = .FALSE.
        LOGICAL :: lrigidsolvent = .FALSE.
        LOGICAL :: lrigidelectrolyte = .FALSE.
        LOGICAL :: lrigidcavity = .FALSE.
        LOGICAL :: lcoredensity = .FALSE.
        LOGICAL :: lsmearedions = .FALSE.
        LOGICAL :: lboundary = .FALSE.
        LOGICAL :: lgradient = .FALSE.
        !
        !--------------------------------------------------------------------------------
        ! Core flags
        !
        LOGICAL :: lfd = .FALSE.
        LOGICAL :: l1da = .FALSE.
        LOGICAL :: lfft_system = .FALSE.
        LOGICAL :: lfft_environment = .FALSE.
        !
        !--------------------------------------------------------------------------------
        ! Solver flags
        !
        LOGICAL :: need_gradient = .FALSE.
        LOGICAL :: need_factsqrt = .FALSE.
        LOGICAL :: need_auxiliary = .FALSE.
        !
        !--------------------------------------------------------------------------------
        ! Corrections flags
        !
        LOGICAL :: need_pbc_correction = .FALSE.
        LOGICAL :: need_electrolyte = .FALSE.
        LOGICAL :: need_semiconductor = .FALSE.
        LOGICAL :: need_outer_loop = .FALSE.
        !
        !--------------------------------------------------------------------------------
        ! Simulation space
        !
        TYPE(environ_cell) :: system_cell
        TYPE(environ_cell), POINTER :: environment_cell => NULL()
        TYPE(environ_mapping) :: mapping
        !
        !--------------------------------------------------------------------------------
        ! Derivatives
        !
        TYPE(container_derivatives) :: derivatives
        !
        !--------------------------------------------------------------------------------
        ! Electrostatic
        !
        TYPE(electrostatic_setup) :: reference, outer, inner
        TYPE(container_electrostatics) :: reference_cores, outer_cores, inner_cores
        TYPE(container_corrections) :: pbc_core
        !
        !--------------------------------------------------------------------------------
        ! Numerical solvers
        !
        TYPE(solver_direct) :: reference_direct, direct
        TYPE(solver_gradient) :: gradient, inner_gradient
        TYPE(solver_fixedpoint) :: fixedpoint, inner_fixedpoint
        TYPE(solver_newton) :: newton
        !
        !--------------------------------------------------------------------------------
        ! Numerical cores
        !
        TYPE(core_fd_derivatives) :: core_fd
        TYPE(core_fft_electrostatics) :: core_fft_sys
        TYPE(core_fft_derivatives) :: core_fft_deriv
        TYPE(core_fft_electrostatics) :: core_fft_elect
        TYPE(core_1da_electrostatics) :: core_1da_elect
        !
        !--------------------------------------------------------------------------------
    CONTAINS
        !--------------------------------------------------------------------------------
        !
        PROCEDURE :: read_input => environ_read_input
        !
        PROCEDURE :: set_flags => set_environ_flags
        PROCEDURE :: set_numerical_base => set_environ_numerical_base
        !
        PROCEDURE :: init_cell => environ_init_cell
        PROCEDURE :: update_cell => environ_update_cell
        PROCEDURE :: end_cell_update => environ_end_cell_update
        !
        PROCEDURE :: update_mapping => environ_update_mapping
        !
        PROCEDURE :: init_cores => init_environ_numerical_cores
        PROCEDURE :: update_cores => update_environ_numerical_cores
        !
        !
        !
        PROCEDURE, PRIVATE :: set_core_containers => set_environ_core_containers
        PROCEDURE, PRIVATE :: set_electrostatics => set_environ_electrostatic
        !
        PROCEDURE, PRIVATE :: &
            set_execution_flags, &
            set_environment_flags, &
            set_dielectric_flags, &
            set_derived_flags
        !
        PROCEDURE :: printout => environ_setup_summary
        !
        !--------------------------------------------------------------------------------
    END TYPE environ_setup
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
    !! Read Environ input
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE environ_read_input(this, filename)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CHARACTER(LEN=80), INTENT(IN), OPTIONAL :: filename
        !
        CLASS(environ_setup), INTENT(INOUT) :: this
        !
        INTEGER :: i
        !
        CHARACTER(LEN=80) :: local_filename = 'environ.in'
        !
        CHARACTER(LEN=80) :: sub_name = 'environ_read_input'
        !
        !--------------------------------------------------------------------------------
        !
        IF (PRESENT(filename)) local_filename = filename
        !
        CALL env_read_input(local_filename)
        !
        global_verbose = verbose ! set internal verbosity from input
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE environ_read_input
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE set_environ_flags(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_setup), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: sub_name = 'set_environ_flags'
        !
        !--------------------------------------------------------------------------------
        !
        IF (pbc_dim >= 0) THEN
            !
            SELECT CASE (TRIM(ADJUSTL(pbc_correction)))
                !
            CASE ('none')
                !
            CASE ('parabolic')
                this%need_pbc_correction = .TRUE.
                !
            CASE ('gcs', 'gouy-chapman', 'gouy-chapman-stern')
                this%need_pbc_correction = .TRUE.
                this%need_electrolyte = .TRUE.
                !
            CASE ('ms', 'mott-schottky')
                this%need_pbc_correction = .TRUE.
                this%need_semiconductor = .TRUE.
                !
            CASE ('ms-gcs', 'mott-schottky-gouy-chapman-stern')
                this%need_pbc_correction = .TRUE.
                this%need_semiconductor = .TRUE.
                this%need_outer_loop = .TRUE.
                this%need_electrolyte = .TRUE.
                !
            CASE DEFAULT
                CALL env_errore(sub_name, 'Option not yet implemented', 1)
                !
            END SELECT
            !
        END IF
        !
        !--------------------------------------------------------------------------------
        !
        CALL this%set_execution_flags()
        !
        CALL this%set_environment_flags()
        !
        CALL this%set_dielectric_flags()
        !
        CALL this%set_derived_flags()
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE set_environ_flags
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE set_environ_numerical_base(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_setup), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: sub_name = 'init_environ_numerical_base'
        !
        !--------------------------------------------------------------------------------
        !
        CALL this%set_core_containers()
        !
        CALL this%set_electrostatics()
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE set_environ_numerical_base
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE environ_init_cell(this, gcutm, comm_in, at)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: comm_in
        REAL(DP), INTENT(IN) :: at(3, 3)
        REAL(DP), INTENT(IN) :: gcutm
        !
        CLASS(environ_setup), TARGET, INTENT(INOUT) :: this
        !
        INTEGER :: ipol
        INTEGER :: environment_nr(3)
        REAL(DP) :: environment_at(3, 3)
        !
        CHARACTER(LEN=80) :: local_label 
        !
        CHARACTER(LEN=80) :: sub_name = 'environ_init_cell'
        !
        !--------------------------------------------------------------------------------
        !
        CALL this%system_cell%init(gcutm, comm_in, at)
        !
        !--------------------------------------------------------------------------------
        ! Double cell and mapping
        !
        IF (this%ldoublecell) THEN
            ALLOCATE (this%environment_cell)
            !
            !----------------------------------------------------------------------------
            ! Scale environment lattice (and corresponding ffts) by 2 * nrep(i) + 1
            !
            DO ipol = 1, 3
                environment_at(:, ipol) = at(:, ipol) * (2.D0 * env_nrep(ipol) + 1.D0)
            END DO
            !
            environment_nr(1) = this%system_cell%dfft%nr1 * (2 * env_nrep(1) + 1)
            environment_nr(2) = this%system_cell%dfft%nr2 * (2 * env_nrep(2) + 1)
            environment_nr(3) = this%system_cell%dfft%nr3 * (2 * env_nrep(3) + 1)
            !
            local_label = 'environment'
            !
            CALL this%environment_cell%init(gcutm, comm_in, environment_at, &
                                            environment_nr, local_label)
            !
        ELSE
            this%environment_cell => this%system_cell
        END IF
        !
        CALL this%mapping%init(env_nrep, this%system_cell, this%environment_cell)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE environ_init_cell
    !------------------------------------------------------------------------------------
    !>
    !! Set up active numerical cores
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_environ_numerical_cores(this, gcutm, use_internal_pbc_corr)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP), INTENT(IN) :: gcutm
        LOGICAL, INTENT(IN), OPTIONAL :: use_internal_pbc_corr
        !
        CLASS(environ_setup), INTENT(INOUT) :: this
        !
        !--------------------------------------------------------------------------------
        !
        IF (this%lfd) CALL this%core_fd%init(ifdtype, nfdpoint, this%environment_cell)
        !
        IF (this%lfft_system) &
            CALL this%core_fft_sys%init(gcutm, this%system_cell, use_internal_pbc_corr)
        !
        IF (this%lfft_environment) THEN
            !
            CALL this%core_fft_deriv%init(gcutm, this%environment_cell)
            !
            IF (this%lelectrostatic) &
                CALL this%core_fft_elect%init(gcutm, this%environment_cell, &
                                              use_internal_pbc_corr)
            !
        END IF
        !
        IF (this%l1da) CALL this%core_1da_elect%init(pbc_dim, pbc_axis, &
                                                     this%environment_cell)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_environ_numerical_cores
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !                                   UPDATE METHODS
    !
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !>
    !! Initialize the cell-related quantities to be used in the Environ
    !! modules. This initialization is called by electrons.f90, thus it
    !! is performed at every step of electronic optimization.
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE environ_update_cell(this, at)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP), INTENT(IN) :: at(3, 3)
        !
        CLASS(environ_setup), INTENT(INOUT) :: this
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
        CALL this%system_cell%printout()
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
            CALL this%environment_cell%printout()
            !
        END IF
        !
        CALL this%update_cores()
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE environ_update_cell
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE environ_end_cell_update(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_setup), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: sub_name = 'environ_end_cell_update'
        !
        !--------------------------------------------------------------------------------
        !
        this%system_cell%lupdate = .FALSE.
        !
        IF (this%ldoublecell) this%environment_cell%lupdate = .FALSE.
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE environ_end_cell_update
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE update_environ_numerical_cores(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_setup), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: sub_name = 'update_environ_numerical_cores'
        !
        !--------------------------------------------------------------------------------
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
        IF (this%l1da) CALL this%core_1da_elect%update_cell(this%environment_cell)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE update_environ_numerical_cores
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE environ_update_mapping(this, pos)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP), INTENT(IN) :: pos(3)
        !
        CLASS(environ_setup), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: sub_name = 'environ_update_mapping'
        !
        !--------------------------------------------------------------------------------
        !
        CALL this%mapping%update(pos)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE environ_update_mapping
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
    SUBROUTINE set_execution_flags(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_setup), INTENT(INOUT) :: this
        !
        !--------------------------------------------------------------------------------
        !
        this%restart = environ_restart
        this%threshold = environ_thr
        this%nskip = environ_nskip
        this%ldoublecell = SUM(env_nrep) > 0
        this%lperiodic = this%need_pbc_correction
        this%louterloop = this%need_outer_loop
        this%ltddfpt = prog == 'TD'
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE set_execution_flags
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE set_environment_flags(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_setup), INTENT(INOUT) :: this
        !
        !--------------------------------------------------------------------------------
        !
        this%surface_tension = &
            env_surface_tension * 1.D-3 * BOHR_RADIUS_SI**2 / RYDBERG_SI
        !
        this%lsurface = this%surface_tension > 0.D0
        !
        this%pressure = env_pressure * 1.D9 / RYDBERG_SI * BOHR_RADIUS_SI**3
        this%lvolume = this%pressure /= 0.D0
        !
        this%confine = env_confine
        this%lconfine = this%confine /= 0.D0
        !
        this%lexternals = env_external_charges > 0
        this%lelectrolyte = env_electrolyte_ntyp > 0 .OR. this%need_electrolyte
        this%lsemiconductor = this%need_semiconductor
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE set_environment_flags
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE set_dielectric_flags(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_setup), INTENT(INOUT) :: this
        !
        INTEGER :: i
        !
        !--------------------------------------------------------------------------------
        !
        this%static_permittivity = env_static_permittivity
        this%optical_permittivity = env_optical_permittivity
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
        this%ldielectric = this%lstatic .OR. this%loptical
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE set_dielectric_flags
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE set_derived_flags(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_setup), INTENT(INOUT) :: this
        !
        INTEGER :: i
        !
        !--------------------------------------------------------------------------------
        ! Derived flags
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
    END SUBROUTINE set_derived_flags
    !------------------------------------------------------------------------------------
    !>
    !! Initialize derivative and electrostatic core containers
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE set_environ_core_containers(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_setup), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: local_type
        !
        CHARACTER(LEN=80) :: sub_name = 'init_environ_core_containers'
        !
        !--------------------------------------------------------------------------------
        ! Calling program reference core
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
        ! Outer/inner cores
        !
        SELECT CASE (core)
            !
        CASE ('fft')
            this%lfft_environment = .TRUE.
            local_type = 'fft'
            !
            CALL this%outer_cores%init(this%core_fft_elect, local_type)
            !
            IF (inner_solver /= 'none') &
                CALL this%inner_cores%init(this%core_fft_elect, local_type)
            !
        CASE ('1d-analytic', '1da')
            this%l1da = .TRUE.
            local_type = '1d-analytic'
            !
            CALL this%outer_cores%init(this%core_1da_elect, local_type)
            !
            IF (inner_solver /= 'none') &
                CALL this%inner_cores%init(this%core_1da_elect, local_type)
            !
        CASE DEFAULT
            !
            CALL env_errore(sub_name, &
                            'Unexpected value for core container type keyword', 1)
            !
        END SELECT
        !
        !--------------------------------------------------------------------------------
        ! Correction cores
        !
        IF (pbc_dim >= 0) THEN
            !
            SELECT CASE (TRIM(ADJUSTL(pbc_correction)))
                !
            CASE ('none')
                !
            CASE ('parabolic')
                this%l1da = .TRUE.
                local_type = '1da'
                !
            CASE ('gcs', 'gouy-chapman', 'gouy-chapman-stern')
                this%l1da = .TRUE.
                local_type = 'gcs'
                !
            CASE ('ms', 'mott-schottky')
                this%l1da = .TRUE.
                local_type = 'ms'
                !
            CASE ('ms-gcs', 'mott-schottky-gouy-chapman-stern')
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
        IF (this%need_pbc_correction) CALL this%outer_cores%add_correction(this%pbc_core)
        !
        IF (inner_solver /= 'none' .AND. this%need_pbc_correction) &
            CALL this%inner_cores%add_correction(this%pbc_core)
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
    END SUBROUTINE set_environ_core_containers
    !------------------------------------------------------------------------------------
    !>
    !! Initialize electrostatic solvers and setups
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE set_environ_electrostatic(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_setup), TARGET, INTENT(INOUT) :: this
        !
        CLASS(electrostatic_solver), POINTER :: local_outer_solver, local_inner_solver
        !
        LOGICAL :: lconjugate
        !
        CHARACTER(LEN=80) :: local_auxiliary, local_problem, inner_problem
        !
        CHARACTER(LEN=80) :: sub_name = 'init_environ_electrostatic'
        !
        !--------------------------------------------------------------------------------
        ! Calling program reference solver
        !
        SELECT CASE (prog)
            !
        CASE ('PW', 'CP', 'TD', 'XS')
            CALL this%reference_direct%init_direct(this%reference_cores)
            !
        CASE DEFAULT
            CALL env_errore(sub_name, 'Unexpected name of host code', 1)
            !
        END SELECT
        !
        !--------------------------------------------------------------------------------
        ! Outer solver
        !
        SELECT CASE (solver)
            !
        CASE ('direct')
            CALL this%direct%init_direct(this%outer_cores)
            !
            local_outer_solver => this%direct
            !
        CASE ('cg', 'sd')
            !
            IF (TRIM(ADJUSTL(solver)) == 'cg') lconjugate = .TRUE.
            !
            CALL this%gradient%init(lconjugate, step_type, step, preconditioner, &
                                    screening_type, screening, this%outer_cores, &
                                    maxstep, tol, auxiliary)
            !
            local_outer_solver => this%gradient
            !
        CASE ('fp')
            !
            CALL this%fixedpoint%init(mix_type, mix, ndiis, this%outer_cores, maxstep, &
                                      tol, auxiliary)
            !
            local_outer_solver => this%fixedpoint
            !
        CASE ('newton')
            CALL this%newton%init(this%outer_cores, maxstep, tol, auxiliary)
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
        ! Inner solver
        !
        IF (inner_solver /= 'none') THEN
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
                        CALL this%inner_gradient%init( &
                            lconjugate, step_type, step, preconditioner, &
                            screening_type, screening, this%inner_cores, inner_maxstep, &
                            inner_tol, auxiliary)
                        !
                        local_inner_solver => this%inner_gradient
                        !
                    CASE ('fp')
                        local_auxiliary = 'full'
                        !
                        CALL this%inner_fixedpoint%init( &
                            mix_type, inner_mix, ndiis, this%inner_cores, &
                            inner_maxstep, inner_tol, local_auxiliary)
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
                CALL this%inner_gradient%init( &
                    lconjugate, step_type, step, preconditioner, screening_type, &
                    screening, this%inner_cores, inner_maxstep, inner_tol, auxiliary)
                !
                local_inner_solver => this%inner_gradient
                !
            END SELECT
            !
        END IF
        !
        !--------------------------------------------------------------------------------
        ! Calling program reference setup
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
        ! Outer setup
        !
        CALL this%outer%init(problem, local_outer_solver)
        !
        !--------------------------------------------------------------------------------
        ! Inner setup
        !
        IF (inner_solver /= 'none') THEN
            !
            CALL this%inner%init(inner_problem, local_inner_solver)
            !
            this%outer%inner => this%inner
        END IF
        !
        !--------------------------------------------------------------------------------
        ! Set electrostatic flags
        !
        CALL this%reference%set_flags(this%need_auxiliary, this%need_gradient, &
                                      this%need_factsqrt)
        !
        CALL this%outer%set_flags(this%need_auxiliary, this%need_gradient, &
                                  this%need_factsqrt)
        !
        IF (inner_solver /= 'none') &
            CALL this%inner%set_flags(this%need_auxiliary, this%need_gradient, &
                                      this%need_factsqrt)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE set_environ_electrostatic
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
    SUBROUTINE environ_setup_summary(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_setup), INTENT(INOUT) :: this
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
        WRITE (program_unit, 1002) environ_thr
        !
        IF (this%lsolvent) THEN
            !
            IF (stype == 0) THEN
                WRITE (program_unit, 1003) 'Fatteber-Gygi'
                WRITE (program_unit, 1004) (rhomax + rhomin) * 0.5_DP, tbeta
            ELSE IF (stype == 1) THEN
                WRITE (program_unit, 1003) 'SCCS'
                WRITE (program_unit, 1005) rhomax, rhomin
            END IF
            !
            IF (solvent_radius > 0.D0) WRITE (program_unit, 1006)
            !
            IF (field_awareness > 0.D0) THEN
                WRITE (program_unit, 1007)
                WRITE (program_unit, 1008) field_awareness, charge_asymmetry
                WRITE (program_unit, 1009) field_min, field_max
            END IF
            !
        END IF
        !
        IF (env_static_permittivity > 1.D0) THEN
            WRITE (program_unit, 1010) env_static_permittivity
            !
            IF (this%ltddfpt) WRITE (program_unit, 1011) env_optical_permittivity
            !
            WRITE (program_unit, 1012) TRIM(solvent_mode)
        END IF
        !
        IF (this%surface_tension > 0.D0) &
            WRITE (program_unit, 1013) env_surface_tension, this%surface_tension
        !
        IF (this%pressure /= 0.D0) &
            WRITE (program_unit, 1014) env_pressure, this%pressure
        !
        !--------------------------------------------------------------------------------
        ! Electrostatic Summary
        !
        IF (this%lelectrostatic) THEN
            WRITE (program_unit, 1015)
            WRITE (program_unit, 1016) problem, solver
            WRITE (program_unit, 1017) auxiliary
            WRITE (program_unit, 1018) core
            !
            IF (this%need_pbc_correction) WRITE (program_unit, 1019) '1d-analytic'
            !
            IF (derivatives == 'fd') THEN
                !
                IF (ifdtype == 1) THEN
                    WRITE (program_unit, 1020) 'central diff.', nfdpoint
                ELSE IF (ifdtype == 2 .OR. ifdtype == 3) THEN
                    WRITE (program_unit, 1020) 'lanczos diff.', nfdpoint
                ELSE IF (ifdtype == 4 .OR. ifdtype == 5) THEN
                    WRITE (program_unit, 1020) 'noise-robust diff.', nfdpoint
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
1001    FORMAT(/, 5X, 'Please cite', /, 9X, A80, /, &
                5X, 'in publications or presentations arising from this work.',/)
        !
1002    FORMAT('     compensation onset threshold      = ', E24.4)
1003    FORMAT('     switching function adopted        = ', A24)
        !
1004    FORMAT('     solvation density threshold       = ', E24.4, /, &
               '     smoothness exponent (2 x beta)    = ', F24.2)
        !
1005    FORMAT('     density limit for vacuum region   = ', E24.4, /, &
               '     density limit for bulk solvent    = ', E24.4)
        !
1006    FORMAT('     interface is solvent aware            ')
1007    FORMAT('     interface is field aware            ')
        !
1008    FORMAT('     field aware factor                = ', F24.2, /, &
               '     asymmetry of field-awareness      = ', F24.2)
        !
1009    FORMAT('     field limit for no correction     = ', F24.2, /, &
               '     field limit for full correction   = ', F24.2)
        !
1010    FORMAT('     static permittivity               = ', F24.2)
1011    FORMAT('     optical permittivity              = ', F24.4)
1012    FORMAT('     epsilon calculation mode          = ', A24)
        !
1013    FORMAT('     surface tension in input (dyn/cm) = ', F24.2, /, &
               '     surface tension in internal units = ', E24.4)
        !
1014    FORMAT('     external pressure in input (GPa)  = ', F24.2, /, &
               '     external pressure in inter. units = ', E24.4)
        !
1015    FORMAT(/, 5X, 'Electrostatic Setup', /, 5X, '-------------------')
        !
1016    FORMAT('     electrostatic problem to solve    = ', A24, /, &
               '     numerical solver adopted          = ', A24)
        !
1017    FORMAT('     type of auxiliary density adopted = ', A24)
1018    FORMAT('     type of core tool for poisson     = ', A24)
1019    FORMAT('     type of core tool for correction  = ', A24)
        !
1020    FORMAT('     type of numerical differentiator  = ', A24, /, &
               '     number of points in num. diff.    = ', I24)
        !
1021    FORMAT(/)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE environ_setup_summary
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE env_setup
!----------------------------------------------------------------------------------------
