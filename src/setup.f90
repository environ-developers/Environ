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
    USE env_mp, ONLY: env_mp_start
    !
    USE class_io, ONLY: io
    !
    USE environ_param, ONLY: DP, tpi2, BOHR_RADIUS_SI, RYDBERG_SI
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
    USE class_core_fft_lowpass
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
        LOGICAL :: has_numerical_setup = .FALSE.
        !
        !--------------------------------------------------------------------------------
        ! Main flags
        !
        LOGICAL :: restart = .FALSE.
        REAL(DP) :: threshold = 0.0_DP
        INTEGER :: nskip = 0
        INTEGER :: niter_scf = 0
        INTEGER :: niter_ionic = 0
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
        LOGICAL :: lmsgcs = .FALSE.
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
        CLASS(core_fft), ALLOCATABLE :: env_fft, ref_fft
        TYPE(core_1da) :: env_1da
        !
        !--------------------------------------------------------------------------------
    CONTAINS
        !--------------------------------------------------------------------------------
        !
        PROCEDURE, PRIVATE :: create => create_environ_setup
        PROCEDURE :: init => init_environ_setup
        !
        PROCEDURE :: init_cell => environ_init_cell
        PROCEDURE :: update_cell => environ_update_cell
        PROCEDURE :: end_cell_update => environ_end_cell_update
        !
        PROCEDURE :: init_numerical => init_environ_numerical_base
        PROCEDURE :: update_cores => update_environ_numerical_cores
        !
        PROCEDURE :: update_mapping => environ_update_mapping
        !
        PROCEDURE :: set_tddfpt
        PROCEDURE :: set_restart
        PROCEDURE :: get_threshold
        PROCEDURE :: get_nskip
        PROCEDURE :: get_nnt
        PROCEDURE :: get_nnr
        PROCEDURE :: get_ir_end
        PROCEDURE :: get_nr
        PROCEDURE :: get_coords
        PROCEDURE :: is_tddfpt
        PROCEDURE :: is_msgcs
        PROCEDURE :: is_restart
        PROCEDURE :: has_solvent
        PROCEDURE :: has_electrostatics
        !
        PROCEDURE, PRIVATE :: set_flags => set_environ_flags
        !
        PROCEDURE, PRIVATE :: set_core_containers => set_environ_core_containers
        PROCEDURE, PRIVATE :: set_electrostatics => set_environ_electrostatic
        !
        PROCEDURE, PRIVATE :: set_execution_flags
        PROCEDURE, PRIVATE :: set_simulation_flags
        PROCEDURE, PRIVATE :: set_environment_flags
        PROCEDURE, PRIVATE :: set_derived_flags
        PROCEDURE, PRIVATE :: set_numerical_flags
        PROCEDURE, PRIVATE :: set_electrostatic_flags
        !
        PROCEDURE :: print_summary => environ_setup_summary
        PROCEDURE :: print_potential_warning => print_environ_potential_warning
        PROCEDURE :: print_clocks => print_environ_clocks
        !
        !--------------------------------------------------------------------------------
    END TYPE environ_setup
    !------------------------------------------------------------------------------------
    !
    CHARACTER(LEN=256) :: bibliography(7)
    !
    DATA bibliography/ &
        "O. Andreussi, I. Dabo and N. Marzari, J. Chem. Phys. 136, 064102 (2012)", &
        "I. Timrov, O. Andreussi, A. Biancardi, N. Marzari, and S. Baroni, J. Chem. Phys. 142, 034111 (2015)", &
        "O. Andreussi, N.G. Hoermann, F. Nattino, G. Fisicaro, S. Goedecker, and N. Marzari, &
        &J. Chem. Theory Comput. 15, 1996 (2019)", &
        "F. Nattino, M. Truscott, N. Marzari, and O. Andreussi, J. Chem. Phys. 150, 041722 (2019)", &
        "M. Truscott, O. Andreussi, J. Phys. Chem. B, 123, 16, 3513–3524 (2019)", &
        "Q. Campbell and I. Dabo, Phys. Rev. B 95, 205308 (2017)", &
        "Q. Campbell, D. Fisher and I. Dabo, Phys. Rev. Mat. 3, 015404 (2019)"/
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
    SUBROUTINE create_environ_setup(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_setup), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: routine = 'create_environ_setup'
        !
        !--------------------------------------------------------------------------------
        !
        IF (ASSOCIATED(this%environment_cell)) CALL io%create_error(routine)
        !
        !--------------------------------------------------------------------------------
        !
        this%restart = .FALSE.
        this%threshold = 0.0_DP
        this%nskip = 0
        this%niter_scf = 0
        this%niter_ionic = 0
        this%nrep = 1
        this%static_permittivity = 0.0_DP
        this%optical_permittivity = 0.0_DP
        this%surface_tension = 0.0_DP
        this%pressure = 0.0_DP
        this%confine = 0.0_DP
        this%lstatic = .FALSE.
        this%loptical = .FALSE.
        this%lsurface = .FALSE.
        this%lvolume = .FALSE.
        this%lconfine = .FALSE.
        this%lexternals = .FALSE.
        this%lregions = .FALSE.
        this%lelectrolyte = .FALSE.
        this%lsemiconductor = .FALSE.
        this%lperiodic = .FALSE.
        this%ldoublecell = .FALSE.
        this%ltddfpt = .FALSE.
        this%laddcharges = .FALSE.
        this%ldielectric = .FALSE.
        this%lsolvent = .FALSE.
        this%lelectrostatic = .FALSE.
        this%lsoftsolvent = .FALSE.
        this%lsoftelectrolyte = .FALSE.
        this%lsoftcavity = .FALSE.
        this%lrigidsolvent = .FALSE.
        this%lrigidelectrolyte = .FALSE.
        this%lrigidcavity = .FALSE.
        this%lcoredensity = .FALSE.
        this%lsmearedions = .FALSE.
        this%lboundary = .FALSE.
        this%lgradient = .FALSE.
        this%l1da = .FALSE.
        this%lfft_system = .FALSE.
        this%lfft_environment = .FALSE.
        this%need_inner = .FALSE.
        this%need_gradient = .FALSE.
        this%need_factsqrt = .FALSE.
        this%need_auxiliary = .FALSE.
        this%has_numerical_setup = .FALSE.
        !
        NULLIFY (this%environment_cell)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE create_environ_setup
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_environ_setup(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_setup), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: routine = 'init_environ_setup'
        !
        !--------------------------------------------------------------------------------
        !
        CALL this%create()
        !
        CALL this%set_flags()
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
    SUBROUTINE environ_init_cell(this, comm_in, at, gcutm, nr)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: comm_in
        REAL(DP), INTENT(IN) :: at(3, 3)
        REAL(DP), OPTIONAL, INTENT(IN) :: gcutm
        INTEGER, OPTIONAL, INTENT(IN) :: nr(3)
        !
        CLASS(environ_setup), TARGET, INTENT(INOUT) :: this
        !
        INTEGER :: i
        !
        INTEGER :: nproc, mpime
        !
        INTEGER :: environment_nr(3)
        REAL(DP) :: environment_at(3, 3)
        !
        REAL(DP) :: local_gcutm, at2
        !
        CHARACTER(LEN=80) :: routine = 'environ_init_cell'
        !
        !--------------------------------------------------------------------------------
        !
        IF (at(1, 1) < 1.D0) CALL io%warning("strange lattice parameter", 1002)
        !
        !--------------------------------------------------------------------------------
        ! Set G-vector cutoff value
        !
        IF (PRESENT(gcutm)) THEN ! passed from calling program
            local_gcutm = gcutm
        ELSE IF (env_ecut /= 0.D0) THEN ! derived from user-defined energy cutoff
            local_gcutm = env_ecut / tpi2
        ELSE IF (PRESENT(nr)) THEN ! overestimated from calling program FFT-grid
            at2 = SUM(at(:, 1)**2)
            local_gcutm = CEILING((nr(1) - 3)**2 * 0.25 / at2 + 0.5 / SQRT(at2) * nr(1))
        ELSE
            CALL io%error(routine, "Missing FFT-grid information", 1003)
        END IF
        !
        !--------------------------------------------------------------------------------
        ! Initializing necessary mp buffers
        !
        CALL env_mp_start(nproc, mpime, io%comm)
        !
        !--------------------------------------------------------------------------------
        !
        CALL this%system_cell%init(comm_in, at, local_gcutm, nr, 'system')
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
            DO i = 1, 3
                environment_at(:, i) = at(:, i) * (2.D0 * env_nrep(i) + 1.D0)
            END DO
            !
            environment_nr = this%system_cell%nr * (2 * env_nrep + 1)
            !
            CALL this%environment_cell%init(comm_in, environment_at, local_gcutm, &
                                            environment_nr, 'environment')
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
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_environ_numerical_base(this, use_internal_pbc_corr)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        LOGICAL, OPTIONAL, INTENT(IN) :: use_internal_pbc_corr
        !
        CLASS(environ_setup), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: routine = 'init_environ_numerical_base'
        !
        !--------------------------------------------------------------------------------
        ! Initialize cores
        !
        IF (deriv_lowpass_p1 .gt. 0.D0 .or. deriv_lowpass_p2 .gt. 0.D0) THEN
            IF (deriv_lowpass_p1 .le. 0.D0 .or. deriv_lowpass_p2 .le. 0.D0) &
                CALL io%error(routine, "Both lowpass parameters need to be positive", 1)
            ALLOCATE(core_fft_lowpass :: this%ref_fft)
            ALLOCATE(core_fft_lowpass :: this%env_fft)
        ELSE
            ALLOCATE(core_fft :: this%ref_fft)
            ALLOCATE(core_fft :: this%env_fft)
        ENDIF
        !
        IF (this%lfft_system) &
            CALL this%ref_fft%init(this%system_cell, use_internal_pbc_corr)
        !
        IF (this%lfft_environment) &
            CALL this%env_fft%init(this%environment_cell, use_internal_pbc_corr)
        !
        IF (this%l1da) CALL this%env_1da%init(pbc_dim, pbc_axis, this%environment_cell)
        !
        !--------------------------------------------------------------------------------
        ! Initialize core containers and solvers
        !
        CALL this%set_core_containers(use_internal_pbc_corr)
        !
        IF (this%lelectrostatic) CALL this%set_electrostatics()
        !
        this%has_numerical_setup = .TRUE.
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_environ_numerical_base
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
        INTEGER :: i
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
            DO i = 1, 3
                !
                environment_at(:, i) = at(:, i) * &
                                       (2.D0 * this%mapping%nrep(i) + 1.D0)
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
        CHARACTER(LEN=80) :: routine = 'environ_end_cell_update'
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
        CHARACTER(LEN=80) :: routine = 'update_environ_numerical_cores'
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
        CHARACTER(LEN=80) :: routine = 'environ_update_mapping'
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
        IMPLICIT NONE
        !
        CLASS(environ_setup), INTENT(IN) :: this
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
        IMPLICIT NONE
        !
        CLASS(environ_setup), INTENT(IN) :: this
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
    INTEGER FUNCTION get_nnt(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_setup), INTENT(IN) :: this
        !
        CHARACTER(LEN=80) :: routine = 'get_nnt'
        !
        !--------------------------------------------------------------------------------
        !
        get_nnt = this%system_cell%nnt
        !
        !--------------------------------------------------------------------------------
    END FUNCTION get_nnt
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    INTEGER FUNCTION get_nnr(this)
    !--------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    CLASS(environ_setup), INTENT(IN) :: this
    !
    CHARACTER(LEN=80) :: routine = 'get_nnr'
    !
    !--------------------------------------------------------------------------------
    !
    get_nnr = this%system_cell%nnr
    !
    !--------------------------------------------------------------------------------
    END FUNCTION get_nnr
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    INTEGER FUNCTION get_ir_end(this)
    !--------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    CLASS(environ_setup), INTENT(IN) :: this
    !
    CHARACTER(LEN=80) :: routine = 'get_ir_end'
    !
    !--------------------------------------------------------------------------------
    !
    get_ir_end = this%system_cell%ir_end
    !
    !--------------------------------------------------------------------------------
    END FUNCTION get_ir_end
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    INTEGER FUNCTION get_nr(this, i)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_setup), INTENT(IN) :: this
        INTEGER, INTENT(IN) :: i
        !
        CHARACTER(LEN=80) :: routine = 'get_nr'
        !
        !--------------------------------------------------------------------------------
        !
        get_nr = this%system_cell%nr(i)
        !
        !--------------------------------------------------------------------------------
    END FUNCTION get_nr
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    FUNCTION get_coords(this, nnr) RESULT(coords)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_setup), INTENT(IN) :: this
        INTEGER, INTENT(IN) :: nnr
        !
        REAL(DP) :: coords(3, nnr)
        !
        !--------------------------------------------------------------------------------
        !
        coords = this%system_cell%coords
        !
        !--------------------------------------------------------------------------------
    END FUNCTION get_coords
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    LOGICAL FUNCTION is_tddfpt(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_setup), INTENT(IN) :: this
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
    LOGICAL FUNCTION is_msgcs(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_setup), INTENT(IN) :: this
        !
        !--------------------------------------------------------------------------------
        !
        is_msgcs = this%lmsgcs
        !
        !--------------------------------------------------------------------------------
    END FUNCTION is_msgcs
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    LOGICAL FUNCTION is_restart(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_setup), INTENT(IN) :: this
        !
        !--------------------------------------------------------------------------------
        !
        is_restart = this%restart
        !
        !--------------------------------------------------------------------------------
    END FUNCTION is_restart
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    LOGICAL FUNCTION has_solvent(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_setup), INTENT(IN) :: this
        !
        !--------------------------------------------------------------------------------
        !
        has_solvent = this%lsolvent
        !
        !--------------------------------------------------------------------------------
    END FUNCTION has_solvent
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    LOGICAL FUNCTION has_electrostatics(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_setup), INTENT(IN) :: this
        !
        !--------------------------------------------------------------------------------
        !
        has_electrostatics = this%lelectrostatic
        !
        !--------------------------------------------------------------------------------
    END FUNCTION has_electrostatics
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
        CHARACTER(LEN=80) :: routine = 'set_environ_flags'
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
        CALL this%set_numerical_flags()
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE set_environ_flags
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
        CHARACTER(LEN=80) :: routine = 'set_simulation_flags'
        !
        !--------------------------------------------------------------------------------
        !
        this%ldoublecell = SUM(env_nrep) > 0
        !
        !--------------------------------------------------------------------------------
        ! Correction flags
        !
        SELECT CASE (pbc_correction)
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
        CASE ('ms-gcs') ! mott-schottky + gouy-chapman-stern
            this%lperiodic = .TRUE.
            this%lsemiconductor = .TRUE.
            this%lmsgcs = .TRUE.
            !
        CASE DEFAULT
            CALL io%error(routine, "Unexpected correction type", 1)
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
        this%lsurface = this%surface_tension /= 0.D0 .OR. lsurface
        !
        factor = 1.D9 / RYDBERG_SI * BOHR_RADIUS_SI**3
        this%pressure = env_pressure * factor
        this%lvolume = this%pressure /= 0.D0 .OR. lvolume
        !
        this%confine = env_confine
        this%lconfine = this%confine /= 0.D0
        !
        this%lexternals = env_external_charges > 0
        this%lelectrolyte = this%lelectrolyte .OR. env_electrolyte_ntyp > 0
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
        this%lelectrostatic = (this%ldielectric .OR. this%lelectrolyte .OR. &
                               this%lexternals .OR. this%lperiodic .OR. field_aware) &
                              .AND. (.NOT. no_electrostatics)
        !
        this%lsoftsolvent = this%lsolvent .AND. (solvent_mode == 'electronic' .OR. &
                                                 solvent_mode == 'full' .OR. &
                                                 field_aware)
        !
        this%lsoftelectrolyte = this%lelectrolyte .AND. &
                                (electrolyte_mode == 'electronic' .OR. &
                                 electrolyte_mode == 'full' .OR. &
                                 field_aware)
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
        this%lgradient = this%ldielectric .OR. field_aware .OR. lsurface
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE set_derived_flags
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE set_numerical_flags(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_setup), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: routine = 'set_numerical_flags'
        !
        !--------------------------------------------------------------------------------
        !
        IF (this%lelectrostatic) this%lfft_system = .TRUE.
        !
        IF (this%lboundary .OR. this%lelectrostatic) THEN
            IF (deriv_core == 'fft') this%lfft_environment = .TRUE.
        END IF
        !
        this%need_inner = inner_solver /= 'none'
        !
        IF (this%lperiodic .AND. pbc_core == '1da') this%l1da = .TRUE.
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE set_numerical_flags
    !------------------------------------------------------------------------------------
    !>
    !! Initialize derivative and electrostatic core containers
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE set_environ_core_containers(this, use_internal_pbc_corr)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        LOGICAL, OPTIONAL, INTENT(IN) :: use_internal_pbc_corr
        !
        CLASS(environ_setup), TARGET, INTENT(INOUT) :: this
        !
        CLASS(environ_core), POINTER :: &
            local_outer_core, local_inner_core, local_pbc_core, local_deriv_core
        !
        CHARACTER(LEN=80) :: routine = 'init_environ_core_containers'
        !
        !--------------------------------------------------------------------------------
        ! Environment core containers
        !
        CALL this%outer_container%init('environment', inter_corr=use_internal_pbc_corr)
        !
        !--------------------------------------------------------------------------------
        ! Derivatives core
        !
        IF (this%lboundary) THEN
            !
            SELECT CASE (deriv_core)
                !
            CASE ('fft')
                local_deriv_core => this%env_fft
                !
            CASE DEFAULT
                CALL io%error(routine, "Unexpected derivatives core", 1)
                !
            END SELECT
            !
            CALL this%outer_container%set_derivatives(local_deriv_core)
            !
        END IF
        !
        !--------------------------------------------------------------------------------
        ! Electrostatic cores
        !
        IF (this%lelectrostatic) THEN
            !
            !----------------------------------------------------------------------------
            ! Calling program reference core container
            !
            CALL this%reference_container%init('system', &
                                               elect_core=this%ref_fft, &
                                               inter_corr=use_internal_pbc_corr)
            !
            !----------------------------------------------------------------------------
            ! Inner core container
            !
            IF (this%need_inner) &
                CALL this%inner_container%init('inner', inter_corr=use_internal_pbc_corr)
            !
            !----------------------------------------------------------------------------
            ! Outer core
            !
            SELECT CASE (core)
                !
            CASE ('fft')
                local_outer_core => this%env_fft
                !
            CASE DEFAULT
                CALL io%error(routine, "Unexpected outer core", 1)
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
                CASE DEFAULT
                    CALL io%error(routine, "Unexpected inner core", 1)
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
                local_pbc_core => this%env_1da
                !
            CASE DEFAULT
                CALL io%error(routine, "Unexpected corrections core", 1)
                !
            END SELECT
            !
            CALL this%outer_container%set_corrections(local_pbc_core)
            !
            IF (this%need_inner) &
                CALL this%inner_container%set_corrections(local_pbc_core)
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
        CHARACTER(LEN=80) :: local_auxiliary, local_problem
        !
        CHARACTER(LEN=80) :: routine = 'init_environ_electrostatic'
        !
        !--------------------------------------------------------------------------------
        ! Calling program reference solver
        !
        CALL this%reference_direct%init(this%reference_container)
        !
        !--------------------------------------------------------------------------------
        ! Outer solver
        !
        CALL this%direct%init(this%outer_container, pbc_correction)
        !
        SELECT CASE (solver)
            !
        CASE ('none')
            !
        CASE ('direct')
            local_outer_solver => this%direct
            !
        CASE ('cg', 'sd')
            !
            IF (solver == 'cg') lconjugate = .TRUE.
            !
            CALL this%gradient%init(lconjugate, step_type, step, preconditioner, &
                                    screening_type, screening, this%outer_container, &
                                    this%direct, maxstep, tol, auxiliary)
            !
            local_outer_solver => this%gradient
            !
        CASE ('fixed-point')
            !
            CALL this%fixedpoint%init(mix_type, mix, ndiis, this%outer_container, &
                                      this%direct, maxstep, tol, auxiliary)
            !
            local_outer_solver => this%fixedpoint
            !
        CASE ('newton')
            CALL this%newton%init(this%outer_container, this%direct, maxstep, tol, &
                                  auxiliary)
            !
            local_outer_solver => this%newton
            !
        CASE DEFAULT
            CALL io%error(routine, "Unexpected outer solver", 1)
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
                    !
                    SELECT CASE (inner_solver)
                        !
                    CASE ('cg', 'sd')
                        !
                        IF (inner_solver == 'cg') lconjugate = .TRUE.
                        !
                        CALL this%inner_gradient%init( &
                            lconjugate, step_type, step, preconditioner, &
                            screening_type, screening, this%inner_container, &
                            this%direct, inner_maxstep, inner_tol, auxiliary)
                        !
                        local_inner_solver => this%inner_gradient
                        !
                    CASE ('fixed-point')
                        local_auxiliary = 'full'
                        !
                        CALL this%inner_fixedpoint%init( &
                            mix_type, inner_mix, ndiis, this%inner_container, &
                            this%direct, inner_maxstep, inner_tol, local_auxiliary)
                        !
                        local_inner_solver => this%inner_fixedpoint
                        !
                    END SELECT
                    !
                ELSE
                    !
                    CALL io%error(routine, &
                                  'Unexpected value for auxiliary charge in nested solver', 1)
                    !
                END IF
                !
            CASE ('newton')
                lconjugate = .TRUE.
                !
                CALL this%inner_gradient%init( &
                    lconjugate, step_type, step, preconditioner, screening_type, &
                    screening, this%inner_container, this%direct, inner_maxstep, &
                    inner_tol, auxiliary)
                !
                local_inner_solver => this%inner_gradient
                !
            CASE DEFAULT
                CALL io%error(routine, "Unexpected inner solver", 1)
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
        CALL this%set_electrostatic_flags(this%reference)
        !
        CALL this%set_electrostatic_flags(this%outer)
        !
        IF (this%need_inner) CALL this%set_electrostatic_flags(this%inner)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE set_environ_electrostatic
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE set_electrostatic_flags(this, solver_setup)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(electrostatic_setup), INTENT(IN) :: solver_setup
        !
        CLASS(environ_setup), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: routine = 'set_electrostatic_flags'
        !
        !--------------------------------------------------------------------------------
        !
        SELECT CASE (solver_setup%problem)
            !
        CASE ('none', 'poisson')
            !
        CASE ('generalized', 'linpb', 'linmodpb', 'pb', 'modpb')
            !
            SELECT TYPE (solver => solver_setup%solver)
                !
            TYPE IS (solver_gradient)
                !
                SELECT CASE (preconditioner)
                    !
                CASE ('sqrt')
                    this%need_factsqrt = .TRUE.
                    !
                CASE ('left', 'none')
                    this%need_gradient = .TRUE.
                    !
                CASE DEFAULT
                    CALL io%error(routine, "Unexpected 'preconditioner'", 1)
                    !
                END SELECT
                !
            END SELECT
            !
            SELECT TYPE (solver => solver_setup%solver)
                !
            CLASS IS (solver_iterative)
                !
                IF (solver%auxiliary /= 'none') this%need_auxiliary = .TRUE.
                !
            CLASS DEFAULT
                CALL io%error(routine, "Unexpected solver", 1)
                !
            END SELECT
            !
        CASE DEFAULT
            CALL io%error(routine, "Unexpected 'problem'", 1)
            !
        END SELECT
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE set_electrostatic_flags
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
    SUBROUTINE environ_setup_summary(this, stdout)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_setup), INTENT(IN) :: this
        INTEGER, OPTIONAL, INTENT(IN) :: stdout
        !
        !--------------------------------------------------------------------------------
        !
        IF (.NOT. io%lnode) RETURN
        !
        IF (PRESENT(stdout)) CALL io%update_unit(stdout)
        !
        !--------------------------------------------------------------------------------
        ! Banner
        !
        CALL io%divider(.TRUE.)
        !
        WRITE (io%unit, 1000)
        !
        CALL io%divider(.TRUE.)
        !
        !--------------------------------------------------------------------------------
        ! Citations
        !
        WRITE (io%unit, 1001)
        !
        WRITE (io%unit, 1002) TRIM(bibliography(1))
        !
        IF (this%ltddfpt) WRITE (io%unit, 1002) TRIM(bibliography(2))
        !
        IF (solvent_radius > 0.D0) WRITE (io%unit, 1002) TRIM(bibliography(3))
        !
        IF (this%lelectrolyte) WRITE (io%unit, 1002) TRIM(bibliography(4))
        !
        IF (field_aware) WRITE (io%unit, 1002) TRIM(bibliography(5))
        !
        IF (this%lsemiconductor) WRITE (io%unit, 1002) TRIM(bibliography(6))
        !
        IF (this%lmsgcs) WRITE (io%unit, 1002) TRIM(bibliography(7))
        !
        WRITE (io%unit, 1003)
        !
        CALL io%divider(.TRUE.)
        !
        !--------------------------------------------------------------------------------
        ! General parameters
        !
        WRITE (io%unit, 1004)
        !
        WRITE (io%unit, 1005) environ_thr
        !
        IF (env_static_permittivity > 1.D0) THEN
            WRITE (io%unit, 1006) env_static_permittivity
            !
            IF (this%ltddfpt) WRITE (io%unit, 1007) env_optical_permittivity
            !
        END IF
        !
        IF (this%surface_tension /= 0.D0) &
            WRITE (io%unit, 1008) env_surface_tension, this%surface_tension
        !
        IF (this%pressure /= 0.D0) &
            WRITE (io%unit, 1009) env_pressure, this%pressure
        !
        !--------------------------------------------------------------------------------
        ! Boundary parameters
        !
        IF (this%lsolvent) THEN
            WRITE (io%unit, 1010)
            WRITE (io%unit, 1011) TRIM(solvent_mode)
            WRITE (io%unit, 1012) TRIM(deriv_method)
            !
            IF (this%lelectrolyte) WRITE (io%unit, 1013) TRIM(electrolyte_deriv_method)
            !
            WRITE (io%unit, 1014) TRIM(deriv_core)
            !
            IF (solvent_mode == 'ionic' .OR. solvent_mode == 'full') THEN
                WRITE (io%unit, 1015) 'SCCS'
                WRITE (io%unit, 1016) rhomax, rhomin
            END IF
            !
            IF (solvent_mode == 'ionic') THEN
                WRITE (io%unit, 1015) 'SSCS'
                WRITE (io%unit, 1017) TRIM(radius_mode), softness, alpha
            END IF
            !
            IF (solvent_radius > 0.D0) WRITE (io%unit, 1018)
            !
            IF (field_aware) THEN
                WRITE (io%unit, 1019)
                WRITE (io%unit, 1020) field_factor, field_asymmetry, field_min, field_max
            END IF
            !
        END IF
        !
        !--------------------------------------------------------------------------------
        ! Electrostatic Summary
        !
        IF (this%lelectrostatic) THEN
            WRITE (io%unit, 1021)
            !
            WRITE (io%unit, 1022) &
                TRIM(problem), TRIM(solver), TRIM(auxiliary), TRIM(core)
            !
            IF (inner_solver /= 'none') THEN
                WRITE (io%unit, 1023)
                WRITE (io%unit, 1024) TRIM(inner_solver), TRIM(inner_core)
            END IF
            !
            IF (this%lperiodic) &
                WRITE (io%unit, 1025) TRIM(pbc_correction), TRIM(pbc_core)
            !
        END IF
        !
        CALL io%divider(.TRUE.)
        !
        !--------------------------------------------------------------------------------
        !
1000    FORMAT(34X, "Environ Setup Summary")
        !
1001    FORMAT(5X, "Please cite",/)
        !
1002    FORMAT(10X, A)
        !
1003    FORMAT(/, 5X, "in publications or presentations arising from this work.")
        !
1004    FORMAT(5X, "Parameters", /, 5X, 10('='),/)
        !
1005    FORMAT(5X, "compensation onset threshold      = ", E24.4)
        !
1006    FORMAT(5X, "static permittivity               = ", F24.2)
        !
1007    FORMAT(5X, "optical permittivity              = ", F24.4)
        !
1008    FORMAT(5X, "surface tension in input (dyn/cm) = ", F24.2, /, &
               5X, "surface tension in internal units = ", E24.4)
        !
1009    FORMAT(5X, "external pressure in input (GPa)  = ", F24.2, /, &
               5X, "external pressure in inter. units = ", E24.4)
        !
1010    FORMAT(/, 5X, "Solvent Boundary", /, 5X, 16('='),/)
        !
1011    FORMAT(5X, "solvent mode                      = ", A24)
        !
1012    FORMAT(5X, "derivatives method                = ", A24)
1013    FORMAT(5X, "electrolyte derivatives method    = ", A24)
1014    FORMAT(5X, "numerical core for derivatives    = ", A24)
        !
1015    FORMAT(5X, "switching function adopted        = ", A24)
        !
1016    FORMAT(5X, "density limit for vacuum region   = ", E24.4, /, &
               5X, "density limit for bulk solvent    = ", E24.4)
        !
1017    FORMAT(5X, "soft-sphere radius mode           = ", A24, /, &
               5X, "soft-sphere softness              = ", F24.2, /, &
               5X, "alpha                             = ", F24.2)
        !
1018    FORMAT(5X, "interface is solvent aware")
        !
1019    FORMAT(5X, "interface is field aware")
        !
1020    FORMAT(5X, "field aware factor                = ", F24.2, /, &
               5X, "asymmetry of field-awareness      = ", F24.2, /, &
               5X, "field limit for no correction     = ", F24.2, /, &
               5X, "field limit for full correction   = ", F24.2)
        !
1021    FORMAT(/, 5X, "Electrostatic Setup", /, 5X, 19('='),/)
        !
1022    FORMAT(5X, "electrostatic problem to solve    = ", A24, /, &
               5X, "numerical solver adopted          = ", A24, /, &
               5X, "type of auxiliary density adopted = ", A24, /, &
               5X, "numerical core for poisson        = ", A24)
        !
1023    FORMAT(5X, "adopting a nested solver scheme")
        !
1024    FORMAT(5X, "inner solver                      = ", A24, /, &
               5X, "inner core                        = ", A24)
        !
1025    FORMAT(5X, "type of pbc corrections           = ", A24, /, &
               5X, "numerical core for corrections    = ", A24)
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
        CLASS(environ_setup), INTENT(IN) :: this
        !
        !--------------------------------------------------------------------------------
        !
        IF (.NOT. io%lnode) RETURN
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
        CLASS(environ_setup), INTENT(IN) :: this
        INTEGER, OPTIONAL, INTENT(IN) :: passed_unit
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
