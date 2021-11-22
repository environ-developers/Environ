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
MODULE class_setup
    !------------------------------------------------------------------------------------
    !
    USE class_io, ONLY: io
    !
    USE environ_param, ONLY: DP, BOHR_RADIUS_SI, RYDBERG_SI
    !
    USE env_base_input
    !
    USE class_cell
    USE class_mapping
    !
    USE class_core_container
    !
    USE class_core
    USE class_core_fft
    USE class_core_fft
    USE class_core_1da
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
        ! Set basic logical flags
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
        LOGICAL :: ltddfpt = .FALSE.
        LOGICAL :: laddcharges = .FALSE.
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
        ! Correction flags
        !
        LOGICAL :: use_inter_corr = .FALSE.
        !
        !--------------------------------------------------------------------------------
        ! Core flags
        !
        LOGICAL :: l1da = .FALSE.
        LOGICAL :: lfft_system = .FALSE.
        LOGICAL :: lfft_environment = .FALSE.
        !
        !--------------------------------------------------------------------------------
        ! Solver flags
        !
        LOGICAL :: need_inner = .FALSE.
        LOGICAL :: need_gradient = .FALSE.
        LOGICAL :: need_factsqrt = .FALSE.
        LOGICAL :: need_auxiliary = .FALSE.
        !
        !--------------------------------------------------------------------------------
        ! Simulation space
        !
        TYPE(environ_cell) :: system_cell
        TYPE(environ_cell), POINTER :: environment_cell => NULL()
        TYPE(environ_mapping) :: mapping
        !
        !--------------------------------------------------------------------------------
        ! Electrostatic
        !
        TYPE(electrostatic_setup) :: reference, outer, inner
        !
        !--------------------------------------------------------------------------------
        ! Containers
        !
        TYPE(core_container) :: reference_container
        TYPE(core_container) :: outer_container, inner_container
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
        TYPE(core_fft) :: env_fft, ref_fft
        TYPE(core_1da) :: env_1da
        !
        !--------------------------------------------------------------------------------
    CONTAINS
        !--------------------------------------------------------------------------------
        !
        PROCEDURE :: init => init_environ_setup
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
        PROCEDURE :: set_tddfpt
        PROCEDURE :: set_restart
        PROCEDURE :: get_threshold
        PROCEDURE :: get_nskip
        PROCEDURE :: is_tddfpt
        PROCEDURE :: is_restart
        !
        PROCEDURE, PRIVATE :: set_flags => set_environ_flags
        PROCEDURE, PRIVATE :: set_numerical_base => set_environ_numerical_base
        !
        PROCEDURE, PRIVATE :: set_core_containers => set_environ_core_containers
        PROCEDURE, PRIVATE :: set_electrostatics => set_environ_electrostatic
        !
        PROCEDURE, PRIVATE :: set_execution_flags
        PROCEDURE, PRIVATE :: set_simulation_flags
        PROCEDURE, PRIVATE :: set_environment_flags
        PROCEDURE, PRIVATE :: set_derived_flags
        !
        PROCEDURE :: print_summary => environ_setup_summary
        PROCEDURE :: print_potential_warning => print_environ_potential_warning
        PROCEDURE :: print_clocks => print_environ_clocks
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
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_environ_setup(this, use_internal_pbc_corr)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        LOGICAL, INTENT(IN), OPTIONAL :: use_internal_pbc_corr
        !
        CLASS(environ_setup), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: sub_name = 'init_environ_setup'
        !
        !--------------------------------------------------------------------------------
        ! Internal pbc correction applied directly to FFT cores
        !
        IF (PRESENT(use_internal_pbc_corr)) this%use_inter_corr = use_internal_pbc_corr
        !
        !--------------------------------------------------------------------------------
        !
        CALL this%set_flags()
        !
        CALL this%set_numerical_base()
        !
        !--------------------------------------------------------------------------------
        ! Open Environ output file
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
    END SUBROUTINE init_environ_setup
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE environ_init_cell(this, gcutm, comm_in, at, nr)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: comm_in
        REAL(DP), INTENT(IN) :: at(3, 3)
        REAL(DP), INTENT(IN) :: gcutm
        INTEGER, INTENT(IN), OPTIONAL :: nr(3)
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
        CALL this%system_cell%init(gcutm, comm_in, at, nr)
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
    SUBROUTINE init_environ_numerical_cores(this, gcutm)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP), INTENT(IN) :: gcutm
        !
        CLASS(environ_setup), INTENT(INOUT) :: this
        !
        !--------------------------------------------------------------------------------
        !
        IF (this%lfft_system) &
            CALL this%ref_fft%init(gcutm, this%system_cell, this%use_inter_corr)
        !
        IF (this%lfft_environment) &
            CALL this%env_fft%init(gcutm, this%environment_cell, this%use_inter_corr)
        !
        IF (this%l1da) CALL this%env_1da%init(pbc_dim, pbc_axis, this%environment_cell)
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
        IF (this%lfft_system) CALL this%ref_fft%update_cell(this%system_cell)
        !
        IF (this%lfft_environment) CALL this%env_fft%update_cell(this%environment_cell)
        !
        IF (this%l1da) CALL this%env_1da%update_cell(this%environment_cell)
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
    !                                   ACCESS METHODS
    !
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE set_tddfpt(this, flag)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        LOGICAL, INTENT(IN) :: flag
        !
        CLASS(environ_setup), INTENT(INOUT) :: this
        !
        !--------------------------------------------------------------------------------
        !
        this%ltddfpt = flag
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE set_tddfpt
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE set_restart(this, flag)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        LOGICAL, INTENT(IN) :: flag
        !
        CLASS(environ_setup), INTENT(INOUT) :: this
        !
        !--------------------------------------------------------------------------------
        !
        this%restart = flag
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE set_restart
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    REAL(DP) FUNCTION get_threshold(this)
        !--------------------------------------------------------------------------------
        !
        CLASS(environ_setup), INTENT(INOUT) :: this
        !
        !--------------------------------------------------------------------------------
        !
        get_threshold = this%threshold
        !
        !--------------------------------------------------------------------------------
    END FUNCTION get_threshold
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    INTEGER FUNCTION get_nskip(this)
        !--------------------------------------------------------------------------------
        !
        CLASS(environ_setup), INTENT(INOUT) :: this
        !
        !--------------------------------------------------------------------------------
        !
        get_nskip = this%nskip
        !
        !--------------------------------------------------------------------------------
    END FUNCTION get_nskip
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    LOGICAL FUNCTION is_tddfpt(this)
        !--------------------------------------------------------------------------------
        !
        CLASS(environ_setup), INTENT(INOUT) :: this
        !
        !--------------------------------------------------------------------------------
        !
        is_tddfpt = this%ltddfpt
        !
        !--------------------------------------------------------------------------------
    END FUNCTION is_tddfpt
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    LOGICAL FUNCTION is_restart(this)
        !--------------------------------------------------------------------------------
        !
        CLASS(environ_setup), INTENT(INOUT) :: this
        !
        !--------------------------------------------------------------------------------
        !
        is_restart = this%restart
        !
        !--------------------------------------------------------------------------------
    END FUNCTION is_restart
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
        CALL this%set_execution_flags()
        !
        CALL this%set_simulation_flags()
        !
        CALL this%set_environment_flags()
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
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE set_execution_flags
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE set_simulation_flags(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_setup), INTENT(INOUT) :: this
        !
        !--------------------------------------------------------------------------------
        !
        this%ldoublecell = SUM(env_nrep) > 0
        !
        !--------------------------------------------------------------------------------
        ! Correction flags
        !
        SELECT CASE (TRIM(ADJUSTL(pbc_correction)))
            !
        CASE ('none')
            !
        CASE ('parabolic')
            this%lperiodic = .TRUE.
            !
        CASE ('gcs') ! gouy-chapman-stern
            this%lperiodic = .TRUE.
            this%lelectrolyte = .TRUE.
            !
        CASE ('ms') ! mott-schottky
            this%lperiodic = .TRUE.
            this%lsemiconductor = .TRUE.
            !
        END SELECT
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE set_simulation_flags
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
        INTEGER :: i
        REAL(DP) :: factor
        !
        !--------------------------------------------------------------------------------
        !
        factor = 1.D-3 / RYDBERG_SI * BOHR_RADIUS_SI**2
        this%surface_tension = env_surface_tension * factor
        this%lsurface = this%surface_tension > 0.D0
        !
        factor = 1.D9 / RYDBERG_SI * BOHR_RADIUS_SI**3
        this%pressure = env_pressure * factor
        this%lvolume = this%pressure /= 0.D0
        !
        this%confine = env_confine
        this%lconfine = this%confine /= 0.D0
        !
        this%lexternals = env_external_charges > 0
        this%lelectrolyte = env_electrolyte_ntyp > 0
        !
        !--------------------------------------------------------------------------------
        ! Dielectric flags
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
    END SUBROUTINE set_environment_flags
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
        CLASS(environ_setup), TARGET, INTENT(INOUT) :: this
        !
        CLASS(environ_core), POINTER :: &
            local_outer_core, local_inner_core, local_pbc_core, local_deriv_core
        !
        CHARACTER(LEN=80) :: local_label
        !
        CHARACTER(LEN=80) :: sub_name = 'init_environ_core_containers'
        !
        !--------------------------------------------------------------------------------
        ! Calling program reference core
        !
        this%lfft_system = .TRUE.
        local_label = 'system'
        !
        CALL this%reference_container%init(local_label, &
                                           elect_core=this%ref_fft, &
                                           inter_corr=this%use_inter_corr)
        !
        !--------------------------------------------------------------------------------
        ! Environment core containers
        !
        local_label = 'environment'
        !
        CALL this%outer_container%init(local_label, inter_corr=this%use_inter_corr)
        !
        IF (this%need_inner) THEN
            local_label = local_label//'_inner'
            !
            CALL this%inner_container%init(local_label, inter_corr=this%use_inter_corr)
            !
        END IF
        !
        !--------------------------------------------------------------------------------
        ! Derivatives core
        !
        IF (this%lboundary) THEN
            !
            SELECT CASE (deriv_core)
                !
            CASE ('fft')
                this%lfft_environment = .TRUE.
                local_deriv_core => this%env_fft
                !
            END SELECT
            !
            CALL this%outer_container%set_derivatives(local_deriv_core, deriv_method)
            !
        END IF
        !
        !--------------------------------------------------------------------------------
        ! Electrostatic cores
        !
        this%need_inner = inner_solver /= 'none'
        !
        IF (this%lelectrostatic) THEN
            !
            !----------------------------------------------------------------------------
            ! Outer core
            !
            SELECT CASE (core)
                !
            CASE ('fft')
                this%lfft_environment = .TRUE.
                local_outer_core => this%env_fft
                !
            END SELECT
            !
            CALL this%outer_container%set_electrostatics(local_outer_core)
            !
            !----------------------------------------------------------------------------
            ! Inner core
            !
            IF (this%need_inner) THEN
                !
                SELECT CASE (inner_core)
                    !
                CASE ('fft')
                    local_inner_core => this%env_fft
                    !
                END SELECT
                !
                CALL this%inner_container%set_electrostatics(local_inner_core)
                !
            END IF
            !
        END IF
        !
        !--------------------------------------------------------------------------------
        ! Corrections core
        !
        IF (this%lperiodic) THEN
            !
            SELECT CASE (pbc_core)
                !
            CASE ('1da')
                this%l1da = .TRUE.
                local_pbc_core => this%env_1da
                !
            END SELECT
            !
            CALL this%outer_container%set_corrections(local_pbc_core, pbc_correction)
            !
            IF (this%need_inner) &
                CALL this%inner_container%set_corrections(local_pbc_core, &
                                                          pbc_correction)
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
        CALL this%reference_direct%init_direct(this%reference_container)
        !
        !--------------------------------------------------------------------------------
        ! Outer solver
        !
        SELECT CASE (solver)
            !
        CASE ('direct')
            CALL this%direct%init_direct(this%outer_container)
            !
            local_outer_solver => this%direct
            !
        CASE ('cg', 'sd')
            !
            IF (TRIM(ADJUSTL(solver)) == 'cg') lconjugate = .TRUE.
            !
            CALL this%gradient%init(lconjugate, step_type, step, preconditioner, &
                                    screening_type, screening, this%outer_container, &
                                    maxstep, tol, auxiliary)
            !
            local_outer_solver => this%gradient
            !
        CASE ('fixed-point')
            !
            CALL this%fixedpoint%init(mix_type, mix, ndiis, this%outer_container, &
                                      maxstep, tol, auxiliary)
            !
            local_outer_solver => this%fixedpoint
            !
        CASE ('newton')
            CALL this%newton%init(this%outer_container, maxstep, tol, auxiliary)
            !
            local_outer_solver => this%newton
            !
        END SELECT
        !
        !--------------------------------------------------------------------------------
        ! Inner solver
        !
        IF (this%need_inner) THEN
            lconjugate = .FALSE.
            !
            SELECT CASE (solver)
                !
            CASE ('fixed-point')
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
                            screening_type, screening, this%inner_container, &
                            inner_maxstep, inner_tol, auxiliary)
                        !
                        local_inner_solver => this%inner_gradient
                        !
                    CASE ('fixed-point')
                        local_auxiliary = 'full'
                        !
                        CALL this%inner_fixedpoint%init( &
                            mix_type, inner_mix, ndiis, this%inner_container, &
                            inner_maxstep, inner_tol, local_auxiliary)
                        !
                        local_inner_solver => this%inner_fixedpoint
                        !
                    END SELECT
                    !
                ELSE
                    !
                    CALL io%error(sub_name, &
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
                    screening, this%inner_container, inner_maxstep, inner_tol, auxiliary)
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
        local_problem = 'poisson'
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
        IF (this%need_inner) THEN
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
        IF (this%need_inner) &
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
        CHARACTER(LEN=80) :: local_label
        !
        !--------------------------------------------------------------------------------
        !
        IF (.NOT. io%lnode) RETURN
        !
        !--------------------------------------------------------------------------------
        ! Banner
        !
        WRITE (io%unit, *)
        WRITE (io%unit, 1000)
        WRITE (io%unit, 1001) bibliography(1)
        !
        !--------------------------------------------------------------------------------
        ! Environ Summary
        !
        WRITE (io%unit, 1002) environ_thr
        !
        IF (this%lsolvent) THEN
            !
            IF (stype == 0) THEN
                WRITE (io%unit, 1003) 'Fatteber-Gygi'
                WRITE (io%unit, 1004) (rhomax + rhomin) * 0.5_DP, tbeta
            ELSE IF (stype == 1) THEN
                WRITE (io%unit, 1003) 'SCCS'
                WRITE (io%unit, 1005) rhomax, rhomin
            END IF
            !
            IF (solvent_radius > 0.D0) WRITE (io%unit, 1006)
            !
            IF (field_awareness > 0.D0) THEN
                WRITE (io%unit, 1007)
                WRITE (io%unit, 1008) field_awareness, charge_asymmetry
                WRITE (io%unit, 1009) field_min, field_max
            END IF
            !
        END IF
        !
        IF (env_static_permittivity > 1.D0) THEN
            WRITE (io%unit, 1010) env_static_permittivity
            !
            IF (this%ltddfpt) WRITE (io%unit, 1011) env_optical_permittivity
            !
            WRITE (io%unit, 1012) TRIM(solvent_mode)
        END IF
        !
        IF (this%surface_tension > 0.D0) &
            WRITE (io%unit, 1013) env_surface_tension, this%surface_tension
        !
        IF (this%pressure /= 0.D0) &
            WRITE (io%unit, 1014) env_pressure, this%pressure
        !
        !--------------------------------------------------------------------------------
        ! Electrostatic Summary
        !
        IF (this%lelectrostatic) THEN
            WRITE (io%unit, 1015)
            WRITE (io%unit, 1016) problem, solver
            WRITE (io%unit, 1017) auxiliary
            WRITE (io%unit, 1018) core
            !
            local_label = '1d-analytic'
            !
            IF (this%lperiodic) WRITE (io%unit, 1019) ADJUSTL(local_label)
            !
        END IF
        !
        WRITE (io%unit, 1020)
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
1020    FORMAT(/)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE environ_setup_summary
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE print_environ_potential_warning(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_setup), INTENT(INOUT) :: this
        !
        !--------------------------------------------------------------------------------
        !
        IF (.NOT. io%lnode) RETURN
        !
        IF (this%lperiodic) WRITE (io%unit, 1100)
        !
1100    FORMAT(/, &
                5(' '), 'WARNING: you are using the parabolic pbc correction;', /, &
                5(' '), '         the potential shift above must be added to ', /, &
                5(' '), '         band and Fermi energies.')
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE print_environ_potential_warning
    !------------------------------------------------------------------------------------
    !>
    !! Write out the time informations of the Environ dependent calculations.
    !! Called by print_clock_pw.f90
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE print_environ_clocks(this, passed_unit)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_setup), INTENT(INOUT) :: this
        INTEGER, INTENT(IN), OPTIONAL :: passed_unit
        !
        INTEGER :: actual_unit
        !
        !--------------------------------------------------------------------------------
        !
        IF (.NOT. io%lnode) RETURN
        !
        IF (PRESENT(passed_unit)) THEN
            actual_unit = passed_unit
        ELSE
            actual_unit = io%unit
        END IF
        !
        WRITE (actual_unit, *)
        WRITE (actual_unit, '(5X,"Environ routines")')
        !
        !--------------------------------------------------------------------------------
        ! Dielectric subroutines
        !
        IF (this%lelectrostatic) THEN
            !
            CALL env_print_clock('calc_eelect')
            !
            CALL env_print_clock('calc_velect')
            !
            CALL env_print_clock('calc_vgcs')
            !
            CALL env_print_clock('dielectric')
            !
            CALL env_print_clock('electrolyte')
            !
            CALL env_print_clock('calc_felect')
            !
        END IF
        !
        !--------------------------------------------------------------------------------
        !
        IF (this%lsemiconductor) CALL env_print_clock('calc_vms')
        !
        !--------------------------------------------------------------------------------
        ! TDDFT
        !
        IF (this%ltddfpt) CALL env_print_clock('calc_vsolvent_tddfpt')
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE print_environ_clocks
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE class_setup
!----------------------------------------------------------------------------------------
