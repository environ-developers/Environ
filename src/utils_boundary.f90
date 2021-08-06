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
MODULE class_boundary
    !------------------------------------------------------------------------------------
    !
    USE env_base_io, ONLY: ionode, environ_unit, program_unit, verbose, depth
    !
    USE environ_param, ONLY: DP, e2
    !
    USE class_cell
    USE class_density
    USE class_functions
    USE class_gradient
    USE class_hessian
    !
    USE class_core_container_derivatives
    USE class_core_fd
    USE class_core_fft
    !
    USE class_electrons
    USE class_ions
    USE class_system
    !
    USE generate_boundary
    !
    ! USE environ_debugging
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
    TYPE, PUBLIC :: environ_boundary
        !--------------------------------------------------------------------------------
        !
        CHARACTER(LEN=80) :: label ! boundary label
        CHARACTER(LEN=80) :: mode ! choice of the interface
        INTEGER :: update_status = 0
        LOGICAL :: initialized = .FALSE.
        !
        !--------------------------------------------------------------------------------
        ! Parameters for the electrons-dependent interface
        !
        LOGICAL :: need_electrons
        TYPE(environ_electrons), POINTER :: electrons
        !
        !--------------------------------------------------------------------------------
        ! Parameters for the ions-dependent interface
        !
        LOGICAL :: need_ions
        TYPE(environ_ions), POINTER :: ions
        !
        !--------------------------------------------------------------------------------
        ! Parameters for the system-dependent interface
        !
        LOGICAL :: need_system
        TYPE(environ_system), POINTER :: system
        !
        !--------------------------------------------------------------------------------
        !
        TYPE(environ_density) :: scaled ! scaled switching function of interface
        ! varying from 1 (QM region) to 0 (environment region)
        !
        INTEGER :: deriv = 0
        TYPE(environ_gradient) :: gradient
        TYPE(environ_density) :: laplacian
        TYPE(environ_density) :: dsurface
        TYPE(environ_hessian) :: hessian
        !
        TYPE(container_derivatives), POINTER :: derivatives
        !
        !--------------------------------------------------------------------------------
        ! Global properties of the boundary
        !
        REAL(DP) :: volume
        REAL(DP) :: surface
        !
        !--------------------------------------------------------------------------------
        ! Components needed for boundary of density
        !
        INTEGER :: b_type
        REAL(DP) :: rhomax, rhomin, fact
        REAL(DP) :: rhozero, deltarho, tbeta
        REAL(DP) :: const
        TYPE(environ_density) :: density
        !
        TYPE(environ_density) :: dscaled
        TYPE(environ_density) :: d2scaled
        !
        !--------------------------------------------------------------------------------
        ! Components needed for boundary of functions
        !
        REAL(DP) :: alpha ! solvent-dependent scaling factor
        REAL(DP) :: softness ! sharpness of the interface
        TYPE(environ_function), ALLOCATABLE :: soft_spheres(:)
        !
        !--------------------------------------------------------------------------------
        !
        TYPE(environ_function) :: simple ! components needed for boundary of system
        !
        !--------------------------------------------------------------------------------
        ! Components needed for solvent-aware boundary
        !
        LOGICAL :: solvent_aware
        TYPE(environ_function) :: solvent_probe
        REAL(DP) :: filling_threshold, filling_spread
        !
        TYPE(environ_density) :: local
        TYPE(environ_density) :: probe
        TYPE(environ_density) :: filling
        TYPE(environ_density) :: dfilling
        !
        !--------------------------------------------------------------------------------
        ! Components needed for field-aware boundary
        !
        LOGICAL :: field_aware
        REAL(DP) :: field_factor, charge_asymmetry, field_max, field_min
        !
        TYPE(environ_density) :: normal_field
        REAL(DP), ALLOCATABLE :: ion_field(:)
        TYPE(environ_function), ALLOCATABLE :: local_spheres(:)
        TYPE(environ_density), ALLOCATABLE :: dion_field_drho(:)
        REAL(DP), ALLOCATABLE :: partial_of_ion_field(:, :, :)
        !
        !--------------------------------------------------------------------------------
    CONTAINS
        !--------------------------------------------------------------------------------
        !
        PROCEDURE :: create => create_environ_boundary
        PROCEDURE :: init_first => init_environ_boundary_first
        PROCEDURE :: init_second => init_environ_boundary_second
        PROCEDURE :: copy => copy_environ_boundary
        PROCEDURE :: update => update_environ_boundary
        PROCEDURE :: destroy => destroy_environ_boundary
        !
        PROCEDURE :: set_soft_spheres
        !
        PROCEDURE :: vconfine => calc_vconfine
        PROCEDURE :: evolume => calc_evolume
        PROCEDURE :: esurface => calc_esurface
        !
        PROCEDURE, NOPASS :: deconfine_dboundary => calc_deconfine_dboundary
        PROCEDURE, NOPASS :: devolume_dboundary => calc_devolume_dboundary
        PROCEDURE :: desurface_dboundary => calc_desurface_dboundary
        PROCEDURE :: dboundary_dions => calc_dboundary_dions
        PROCEDURE :: sa_de_dboundary => calc_solvent_aware_de_dboundary
        !
        PROCEDURE :: of_density => boundary_of_density
        PROCEDURE :: of_functions => boundary_of_functions
        PROCEDURE :: of_system => boundary_of_system
        !
        PROCEDURE :: convolution_deriv => compute_convolution_deriv
        PROCEDURE :: solvent_aware_boundary
        PROCEDURE :: calc_dsurface ! #TODO do we need this?
        PROCEDURE :: invert => invert_boundary
        !
        PROCEDURE :: printout => print_environ_boundary
        !
        !--------------------------------------------------------------------------------
    END TYPE environ_boundary
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
    SUBROUTINE create_environ_boundary(this, local_label)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CHARACTER(LEN=80), INTENT(IN) :: local_label
        !
        CLASS(environ_boundary), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: label = ' '
        !
        CHARACTER(LEN=80) :: sub_name = 'create_environ_boundary'
        !
        !--------------------------------------------------------------------------------
        !
        this%update_status = 0
        this%label = local_label
        !
        label = 'boundary_'//TRIM(ADJUSTL(local_label))
        !
        CALL this%scaled%create(label)
        !
        this%volume = 0.D0
        this%surface = 0.D0
        !
        this%need_electrons = .FALSE.
        NULLIFY (this%electrons)
        this%need_ions = .FALSE.
        NULLIFY (this%ions)
        this%need_system = .FALSE.
        NULLIFY (this%system)
        !
        !--------------------------------------------------------------------------------
        ! Optional components
        !
        this%deriv = 0
        label = 'gradboundary_'//TRIM(ADJUSTL(local_label))
        !
        CALL this%gradient%create(label)
        !
        label = 'laplboundary_'//TRIM(ADJUSTL(local_label))
        !
        CALL this%laplacian%create(label)
        !
        label = 'dsurface_'//TRIM(ADJUSTL(local_label))
        !
        CALL this%dsurface%create(label)
        !
        label = 'hessboundary_'//TRIM(ADJUSTL(local_label))
        !
        CALL this%hessian%create(label)
        !
        !--------------------------------------------------------------------------------
        ! Components required for boundary of density
        !
        label = 'boundary_density_'//TRIM(ADJUSTL(local_label))
        !
        CALL this%density%create(label)
        !
        label = 'dboundary_'//TRIM(ADJUSTL(local_label))
        !
        CALL this%dscaled%create(label)
        !
        label = 'd2boundary_'//TRIM(ADJUSTL(local_label))
        !
        CALL this%d2scaled%create(label)
        !
        !--------------------------------------------------------------------------------
        ! Components required for boundary of functions
        !
        IF (ALLOCATED(this%soft_spheres)) &
            CALL env_errore(sub_name, 'Trying to create an already allocated object', 1)
        !
        !--------------------------------------------------------------------------------
        ! Components required for solvent-aware interface
        !
        this%solvent_aware = .FALSE.
        label = 'local_'//TRIM(ADJUSTL(local_label))
        !
        CALL this%local%create(label)
        !
        label = 'probe_'//TRIM(ADJUSTL(local_label))
        !
        CALL this%probe%create(label)
        !
        label = 'filling_'//TRIM(ADJUSTL(local_label))
        !
        CALL this%filling%create(label)
        !
        label = 'dfilling_'//TRIM(ADJUSTL(local_label))
        !
        CALL this%dfilling%create(label)
        !
        !--------------------------------------------------------------------------------
        ! Components required for field-aware interface
        !
        this%field_aware = .FALSE.
        label = 'normal_field_'//TRIM(ADJUSTL(local_label))
        !
        CALL this%normal_field%create(label)
        !
        IF (ALLOCATED(this%ion_field)) &
            CALL env_errore(sub_name, 'Trying to create an already allocated object', 1)
        !
        IF (ALLOCATED(this%local_spheres)) &
            CALL env_errore(sub_name, 'Trying to create an already allocated object', 1)
        !
        IF (ALLOCATED(this%dion_field_drho)) &
            CALL env_errore(sub_name, 'Trying to create an already allocated object', 1)
        !
        IF (ALLOCATED(this%partial_of_ion_field)) &
            CALL env_errore(sub_name, 'Trying to create an already allocated object', 1)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE create_environ_boundary
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_environ_boundary_first(this, need_gradient, need_laplacian, &
                                           need_hessian, mode, stype, rhomax, rhomin, &
                                           tbeta, const, alpha, softness, &
                                           system_distance, system_spread, &
                                           solvent_radius, radial_scale, &
                                           radial_spread, filling_threshold, &
                                           filling_spread, field_factor, &
                                           charge_asymmetry, field_max, field_min, &
                                           electrons, ions, system, derivatives)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CHARACTER(LEN=80), INTENT(IN) :: mode
        INTEGER, INTENT(IN) :: stype
        REAL(DP), INTENT(IN) :: rhomax, rhomin, tbeta, const
        LOGICAL, INTENT(IN) :: need_gradient, need_laplacian, need_hessian
        !
        REAL(DP), INTENT(IN) :: alpha, softness, system_distance, system_spread, &
                                solvent_radius, radial_scale, radial_spread, &
                                filling_threshold, filling_spread, field_factor, &
                                charge_asymmetry, field_max, field_min
        !
        TYPE(environ_electrons), TARGET, INTENT(IN) :: electrons
        TYPE(environ_ions), TARGET, INTENT(IN) :: ions
        TYPE(environ_system), TARGET, INTENT(IN) :: system
        TYPE(container_derivatives), TARGET, INTENT(IN) :: derivatives
        !
        CLASS(environ_boundary), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: sub_name = 'init_environ_boundary_first'
        !
        !--------------------------------------------------------------------------------
        !
        IF (need_hessian) THEN
            this%deriv = 3
        ELSE IF (need_laplacian) THEN
            this%deriv = 2
        ELSE IF (need_gradient) THEN
            this%deriv = 1
        END IF
        !
        this%mode = mode
        !
        this%need_electrons = (mode == 'electronic') .OR. &
                              (mode == 'full') .OR. &
                              (mode == 'fa-ionic') .OR. &
                              (mode == 'fa-electronic')
        !
        IF (this%need_electrons) this%electrons => electrons
        !
        this%need_ions = (mode == 'ionic') .OR. &
                         (mode == 'full') .OR. &
                         (mode == 'fa-ionic') .OR. &
                         (mode == 'fa-electronic')
        !
        IF (this%need_ions) this%ions => ions
        !
        this%need_system = (mode == 'system')
        !
        IF (this%need_system) this%system => system
        !
        this%b_type = stype
        this%rhomax = rhomax
        this%rhomin = rhomin
        this%fact = LOG(rhomax / rhomin)
        this%rhozero = (rhomax + rhomin) * 0.5_DP
        this%tbeta = tbeta
        this%deltarho = rhomax - rhomin
        !
        IF (const == 1.D0 .AND. this%need_electrons .AND. stype == 2) &
            CALL env_errore(sub_name, &
                            'stype=2 boundary requires dielectric constant > 1', 1)
        !
        this%const = const
        !
        this%alpha = alpha
        this%softness = softness
        !
        IF (this%mode == 'ionic' .OR. this%mode == 'fa-ionic') &
            ALLOCATE (this%soft_spheres(this%ions%number))
        !
        this%simple%f_type = 4
        this%simple%pos => system%pos
        this%simple%volume = 1.D0
        this%simple%dim = system%dim
        this%simple%axis = system%axis
        this%simple%width = system_distance
        this%simple%spread = system_spread
        !
        this%solvent_aware = solvent_radius > 0.D0
        !
        IF (this%solvent_aware) THEN
            this%solvent_probe%f_type = 2
            ALLOCATE (this%solvent_probe%pos(3))
            this%solvent_probe%pos = 0.D0
            this%solvent_probe%volume = 1.D0
            this%solvent_probe%dim = 0
            this%solvent_probe%axis = 1
            this%solvent_probe%spread = radial_spread
            this%solvent_probe%width = solvent_radius * radial_scale
        END IF
        !
        this%filling_threshold = filling_threshold
        this%filling_spread = filling_spread
        !
        this%derivatives => derivatives
        this%field_aware = field_factor > 0.D0
        this%field_factor = field_factor
        this%charge_asymmetry = charge_asymmetry
        this%field_max = field_max
        this%field_min = field_min
        !
        IF (this%field_aware .AND. this%mode == 'fa-ionic') THEN
            ALLOCATE (this%ion_field(this%ions%number))
            ALLOCATE (this%dion_field_drho(this%ions%number))
            !
            ALLOCATE (this%partial_of_ion_field(3, this%ions%number, &
                                                this%ions%number))
            !
            ALLOCATE (this%local_spheres(this%ions%number))
        END IF
        !
        this%initialized = .FALSE.
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_environ_boundary_first
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_environ_boundary_second(this, cell)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_cell), INTENT(IN) :: cell
        !
        CLASS(environ_boundary), INTENT(INOUT) :: this
        !
        INTEGER :: i
        !
        CHARACTER(LEN=80) :: sub_name = 'init_environ_boundary_second'
        !
        !--------------------------------------------------------------------------------
        !
        CALL this%scaled%init(cell)
        !
        IF (this%mode == 'electronic' .OR. this%mode == 'full' .OR. &
            this%mode == 'fa-electronic' .OR. this%mode == 'fa-full') THEN
            !
            CALL this%density%init(cell)
            !
            CALL this%dscaled%init(cell)
            !
            CALL this%d2scaled%init(cell)
            !
        END IF
        !
        IF (this%deriv >= 1) CALL this%gradient%init(cell)
        !
        IF (this%deriv >= 2) CALL this%laplacian%init(cell)
        !
        IF (this%deriv >= 3) CALL this%dsurface%init(cell)
        !
        IF (this%solvent_aware) THEN
            !
            CALL this%local%init(cell)
            !
            CALL this%probe%init(cell)
            !
            CALL this%filling%init(cell)
            !
            CALL this%dfilling%init(cell)
            !
            IF (this%deriv >= 3) CALL this%hessian%init(cell)
            !
        END IF
        !
        IF (this%field_aware) THEN
            !
            CALL env_errore(sub_name, 'field-aware not yet implimented ', 1)
            !
            IF (this%mode == 'fa-electronic' .OR. &
                !
                this%mode == 'fa-full') THEN
                !
                CALL this%normal_field%init(cell)
                !
            ELSE IF (this%mode == 'fa-ionic') THEN
                !
                DO i = 1, this%ions%number
                    CALL this%dion_field_drho(i)%init(cell)
                END DO
                !
            ELSE
                CALL env_errore(sub_name, 'Boundary must be field-aware', 1)
            END IF
            !
        END IF
        !
        this%initialized = .TRUE.
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_environ_boundary_second
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE copy_environ_boundary(this, copy)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_boundary), INTENT(IN) :: this
        !
        TYPE(environ_boundary), INTENT(OUT) :: copy
        !
        INTEGER :: i, n, m
        !
        !--------------------------------------------------------------------------------
        !
        copy%electrons => this%electrons
        copy%ions => this%ions
        copy%system => this%system
        copy%derivatives => this%derivatives
        !
        copy%mode = this%mode
        copy%update_status = this%update_status
        copy%need_electrons = this%need_electrons
        copy%need_ions = this%need_ions
        copy%need_system = this%need_system
        copy%deriv = this%deriv
        copy%volume = this%volume
        copy%surface = this%surface
        copy%b_type = this%b_type
        copy%rhomax = this%rhomax
        copy%rhomin = this%rhomin
        copy%fact = this%fact
        copy%rhozero = this%rhozero
        copy%deltarho = this%deltarho
        copy%tbeta = this%tbeta
        copy%const = this%const
        copy%alpha = this%alpha
        copy%softness = this%softness
        copy%solvent_aware = this%solvent_aware
        copy%filling_threshold = this%filling_threshold
        copy%filling_spread = this%filling_spread
        copy%field_aware = this%field_aware
        copy%field_factor = this%field_factor
        copy%charge_asymmetry = this%charge_asymmetry
        copy%field_max = this%field_max
        copy%field_min = this%field_min
        copy%initialized = this%initialized
        !
        IF (ASSOCIATED(this%scaled%cell)) CALL this%scaled%copy(copy%scaled)
        !
        IF (ASSOCIATED(this%gradient%cell)) CALL this%gradient%copy(copy%gradient)
        !
        IF (ASSOCIATED(this%laplacian%cell)) CALL this%laplacian%copy(copy%laplacian)
        !
        IF (ASSOCIATED(this%dsurface%cell)) CALL this%dsurface%copy(copy%dsurface)
        !
        IF (ASSOCIATED(this%hessian%cell)) CALL this%hessian%copy(copy%hessian)
        !
        IF (ASSOCIATED(this%density%cell)) CALL this%density%copy(copy%density)
        !
        IF (ASSOCIATED(this%dscaled%cell)) CALL this%dscaled%copy(copy%dscaled)
        !
        IF (ASSOCIATED(this%d2scaled%cell)) CALL this%d2scaled%copy(copy%d2scaled)
        !
        CALL this%simple%copy(copy%simple)
        !
        CALL this%solvent_probe%copy(copy%solvent_probe)
        !
        IF (ASSOCIATED(this%local%cell)) CALL this%local%copy(copy%local)
        !
        IF (ASSOCIATED(this%probe%cell)) CALL this%probe%copy(copy%probe)
        !
        IF (ASSOCIATED(this%filling%cell)) CALL this%filling%copy(copy%filling)
        !
        IF (ASSOCIATED(this%dfilling%cell)) CALL this%dfilling%copy(copy%dfilling)
        !
        IF (ASSOCIATED(this%normal_field%cell)) &
            CALL this%normal_field%copy(copy%normal_field)
        !
        IF (ALLOCATED(this%soft_spheres)) THEN
            n = SIZE(this%soft_spheres)
            !
            IF (ALLOCATED(copy%soft_spheres)) THEN
                m = SIZE(copy%soft_spheres)
                !
                CALL destroy_environ_functions(copy%soft_spheres, m)
                !
            END IF
            !
            CALL copy_environ_functions(this%soft_spheres, n, copy%soft_spheres)
            !
        ELSE
            IF (ALLOCATED(copy%soft_spheres)) DEALLOCATE (copy%soft_spheres)
        END IF
        !
        IF (ALLOCATED(this%ion_field)) THEN
            n = SIZE(this%ion_field)
            !
            IF (ALLOCATED(copy%ion_field)) DEALLOCATE (copy%ion_field)
            !
            IF (ALLOCATED(copy%partial_of_ion_field)) &
                DEALLOCATE (copy%partial_of_ion_field)
            !
            ALLOCATE (copy%ion_field(n))
            ALLOCATE (copy%partial_of_ion_field(3, n, n))
            !
            IF (ALLOCATED(copy%dion_field_drho)) THEN
                m = SIZE(copy%dion_field_drho)
                !
                DO i = 1, m
                    CALL copy%dion_field_drho(i)%destroy()
                END DO
                !
                DEALLOCATE (copy%dion_field_drho)
            END IF
            !
            ALLOCATE (copy%dion_field_drho(n))
            !
            DO i = 1, n
                CALL this%dion_field_drho(i)%copy(copy%dion_field_drho(i))
            END DO
            !
        ELSE
            !
            IF (ALLOCATED(copy%ion_field)) DEALLOCATE (copy%ion_field)
            !
            IF (ALLOCATED(copy%partial_of_ion_field)) &
                DEALLOCATE (copy%partial_of_ion_field)
            !
            IF (ALLOCATED(copy%dion_field_drho)) DEALLOCATE (copy%dion_field_drho)
            !
        END IF
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE copy_environ_boundary
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE update_environ_boundary(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_boundary), INTENT(INOUT) :: this
        !
        LOGICAL :: update_anything
        !
        INTEGER :: i
        TYPE(environ_cell), POINTER :: cell
        CHARACTER(LEN=80) :: label
        !
        CHARACTER(LEN=80) :: sub_name = 'update_environ_boundary'
        !
        !--------------------------------------------------------------------------------
        !
        cell => this%scaled%cell
        !
        update_anything = .FALSE.
        !
        IF (this%need_ions) update_anything = this%ions%lupdate
        !
        IF (this%need_electrons) &
            update_anything = update_anything .OR. this%electrons%lupdate
        !
        IF (this%need_system) &
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
        SELECT CASE (this%mode)
            !
        CASE ('full')
            !
            IF (this%ions%lupdate) THEN
                !
                !------------------------------------------------------------------------
                ! Compute the ionic part
                !
                CALL density_of_functions(this%ions%core_electrons, this%ions%number, &
                                          this%ions%core, .TRUE.)
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
                    CALL env_errore(sub_name, &
                                    'Wrong update status, possibly &
                                    &missing ionic update', 1)
                !
                this%density%of_r = this%electrons%density%of_r + this%ions%core%of_r
                !
                CALL this%of_density()
                !
                this%update_status = 2 ! boundary has changed and is ready
                !
            END IF
            !
        CASE ('electronic')
            !
            IF (this%electrons%lupdate) THEN
                !
                this%density%of_r = this%electrons%density%of_r
                !
                CALL this%of_density()
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
        CASE ('ionic')
            !
            IF (this%ions%lupdate) THEN
                !
                !------------------------------------------------------------------------
                ! Only ions are needed, fully update the boundary
                !
                CALL this%of_functions()
                !
                this%update_status = 2 ! boundary has changed and is ready
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
        CASE ('system')
            !
            IF (this%system%lupdate) THEN
                !
                !------------------------------------------------------------------------
                ! Only ions are needed, fully update the boundary
                !
                CALL this%of_system()
                !
                ! TO DEBUG SOLVENT-AWARE
                ! !
                ! CALL this%invert()
                ! !
                ! CALL test_de_dboundary(this)
                !
                this%update_status = 2 ! boundary has changed and is ready
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
            CALL env_errore(sub_name, 'Unrecognized boundary mode', 1)
            !
        END SELECT
        !
        !--------------------------------------------------------------------------------
        ! Solvent-aware interface
        !
        IF (this%update_status == 2 .AND. this%solvent_aware) &
            CALL this%solvent_aware_boundary()
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE update_environ_boundary
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE destroy_environ_boundary(this, lflag)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        LOGICAL, INTENT(IN) :: lflag
        !
        CLASS(environ_boundary), INTENT(INOUT) :: this
        !
        INTEGER :: i
        !
        CHARACTER(LEN=80) :: sub_name = 'destroy_environ_boundary'
        !
        !--------------------------------------------------------------------------------
        !
        IF (this%initialized) THEN
            !
            CALL this%scaled%destroy()
            !
            IF (this%mode == 'electronic' .OR. this%mode == 'full') THEN
                !
                CALL this%density%destroy()
                !
                CALL this%dscaled%destroy()
                !
                CALL this%d2scaled%destroy()
                !
            END IF
            !
            IF (this%deriv >= 1) CALL this%gradient%destroy()
            !
            IF (this%deriv >= 2) CALL this%laplacian%destroy()
            !
            IF (this%deriv >= 3) CALL this%dsurface%destroy()
            !
            IF (this%solvent_aware) THEN
                !
                CALL this%local%destroy()
                !
                CALL this%probe%destroy()
                !
                CALL this%filling%destroy()
                !
                CALL this%dfilling%destroy()
                !
                IF (this%deriv >= 3) CALL this%hessian%destroy()
                !
            END IF
            !
            this%initialized = .FALSE.
        END IF
        !
        IF (lflag) THEN
            !
            !----------------------------------------------------------------------------
            ! These components were allocated first, destroy only if lflag = .TRUE.
            !
            IF (this%need_ions) THEN
                IF (this%mode == 'ionic' .OR. this%mode == 'fa-ionic') THEN
                    !
                    CALL destroy_environ_functions(this%soft_spheres, this%ions%number)
                    !
                    IF (this%field_aware .AND. this%mode == 'fa-ionic') THEN
                        !
                        CALL env_errore(sub_name, 'field-aware not yet implimented ', 1)
                        !
                        DEALLOCATE (this%ion_field)
                        DEALLOCATE (this%partial_of_ion_field)
                        !
                        CALL destroy_environ_functions(this%local_spheres, &
                                                       this%ions%number)
                        !
                        DEALLOCATE (this%dion_field_drho)
                    END IF
                    !
                END IF
                !
                IF (.NOT. ASSOCIATED(this%ions)) &
                    CALL env_errore(sub_name, &
                                    'Trying to destroy a non associated object', 1)
                !
                NULLIFY (this%ions)
            ELSE
                !
                IF (ASSOCIATED(this%ions)) &
                    CALL env_errore(sub_name, 'Found an unexpected associated object', 1)
                !
            END IF
            !
            IF (this%need_electrons) THEN
                IF (ASSOCIATED(this%electrons)) NULLIFY (this%electrons)
            END IF
            !
            IF (this%solvent_aware) DEALLOCATE (this%solvent_probe%pos)
            !
            IF (this%need_system) THEN
                IF (ASSOCIATED(this%system)) NULLIFY (this%system)
            END IF
            !
        END IF
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE destroy_environ_boundary
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE set_soft_spheres(this, scale)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        LOGICAL, INTENT(IN), OPTIONAL :: scale
        !
        CLASS(environ_boundary), INTENT(INOUT) :: this
        !
        INTEGER :: i
        REAL(DP) :: radius
        !
        !--------------------------------------------------------------------------------
        !
        IF (.NOT. (this%mode == 'ionic' .OR. this%mode == 'fa-ionic')) RETURN
        !
        DO i = 1, this%ions%number
            radius = this%ions%iontype(this%ions%ityp(i))%solvationrad * this%alpha
            !
            CALL this%soft_spheres(i)%init(5, 1, 0, radius, this%softness, 1.D0, &
                                           this%ions%tau(:, i))
            !
        END DO
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE set_soft_spheres
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !                                 EMBEDDING METHODS
    !
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !>
    !! Calculates the confine contribution to the potential
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_vconfine(this, confine, vconfine)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP), INTENT(IN) :: confine
        !
        CLASS(environ_boundary), INTENT(INOUT) :: this
        TYPE(environ_density), INTENT(INOUT) :: vconfine
        !
        !--------------------------------------------------------------------------------
        ! The confine potetial is defined as confine * ( 1 - s(r) )
        !
        vconfine%of_r = 0.D0
        vconfine%of_r = confine * (1.D0 - this%scaled%of_r)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_vconfine
    !------------------------------------------------------------------------------------
    !>
    !! Calculates the confine contribution to the interface potential
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_deconfine_dboundary(confine, rhoelec, de_dboundary)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP), INTENT(IN) :: confine
        CLASS(environ_density), INTENT(IN) :: rhoelec
        !
        CLASS(environ_density), INTENT(INOUT) :: de_dboundary
        !
        !--------------------------------------------------------------------------------
        !
        de_dboundary%of_r = de_dboundary%of_r - confine * rhoelec%of_r
        ! the functional derivative of the confine term is - confine * rho^elec(r)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_deconfine_dboundary
    !------------------------------------------------------------------------------------
    !>
    !! Calculates the PV contribution to the energy
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_evolume(this, pressure, evolume)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP), INTENT(IN) :: pressure
        !
        CLASS(environ_boundary), TARGET, INTENT(INOUT) :: this
        REAL(DP), INTENT(OUT) :: evolume
        !
        !--------------------------------------------------------------------------------
        !
        evolume = pressure * this%volume * e2 / 2.D0 ! computes the PV energy
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_evolume
    !------------------------------------------------------------------------------------
    !>
    !! Calculates the PV contribution to the potential
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_devolume_dboundary(pressure, de_dboundary)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP), INTENT(IN) :: pressure
        !
        TYPE(environ_density), INTENT(INOUT) :: de_dboundary
        !
        !--------------------------------------------------------------------------------
        !
        de_dboundary%of_r = de_dboundary%of_r + pressure
        ! the functional derivative of the volume term is just unity
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_devolume_dboundary
    !------------------------------------------------------------------------------------
    !>
    !! Calculates the cavitation contribution to the energy
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_esurface(this, surface_tension, esurface)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP), INTENT(IN) :: surface_tension
        !
        CLASS(environ_boundary), TARGET, INTENT(INOUT) :: this
        REAL(DP), INTENT(OUT) :: esurface
        !
        !--------------------------------------------------------------------------------
        !
        esurface = surface_tension * this%surface * e2 / 2.D0
        ! computes the cavitation energy
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_esurface
    !------------------------------------------------------------------------------------
    !>
    !! Calculates the cavitation contribution to the potential
    !
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_desurface_dboundary(this, surface_tension, de_dboundary)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP), INTENT(IN) :: surface_tension
        CLASS(environ_boundary), TARGET, INTENT(IN) :: this
        !
        TYPE(environ_density), TARGET, INTENT(INOUT) :: de_dboundary
        !
        !--------------------------------------------------------------------------------
        !
        de_dboundary%of_r = de_dboundary%of_r + surface_tension * this%dsurface%of_r
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_desurface_dboundary
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !                                  GENERAL METHODS
    !
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !>
    !! Calculates the dielectric constant as a function of the charge
    !! density, and the derivatives of the the dielectric constant
    !! with respect to the charge density. Additionally calculates the
    !! volume and surface components.
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE boundary_of_density(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_boundary), TARGET, INTENT(INOUT) :: this
        !
        INTEGER, POINTER :: ir_end, nnr, stype, deriv
        REAL(DP), POINTER :: const, rhomax, rhomin, tbeta
        REAL(DP), DIMENSION(:), POINTER :: rho, eps, deps, d2eps, lapleps, dsurface
        REAL(DP), POINTER :: gradeps(:, :)
        TYPE(environ_hessian), POINTER :: hessian
        !
        INTEGER :: ir, ipol, jpol
        !
        CHARACTER(LEN=80) :: sub_name = 'boundary_of_density'
        !
        !--------------------------------------------------------------------------------
        !
        IF (.NOT. ASSOCIATED(this%density%cell, this%scaled%cell)) &
            CALL env_errore(sub_name, 'Inconsistent domains', 1)
        !
        ir_end => this%density%cell%ir_end
        nnr => this%density%cell%nnr
        rho => this%density%of_r
        !
        stype => this%b_type
        eps => this%scaled%of_r
        deps => this%dscaled%of_r
        d2eps => this%d2scaled%of_r
        !
        IF (stype == 1 .OR. stype == 2) THEN
            rhomax => this%rhomax
            rhomin => this%rhomin
            tbeta => this%fact
            const => this%const
        ELSE IF (stype == 0) THEN
            rhomax => this%rhozero
            rhomin => this%deltarho
            tbeta => this%tbeta
            const => this%const
        END IF
        !
        DO ir = 1, ir_end
            eps(ir) = boundfunct(rho(ir), rhomax, rhomin, tbeta, const, stype)
            deps(ir) = dboundfunct(rho(ir), rhomax, rhomin, tbeta, const, stype)
            d2eps(ir) = d2boundfunct(rho(ir), rhomax, rhomin, tbeta, const, stype)
        END DO
        !
        !--------------------------------------------------------------------------------
        ! Compute boundary derivatives, if needed
        !
        deriv => this%deriv
        !
        IF (deriv >= 1) gradeps => this%gradient%of_r
        !
        IF (deriv >= 2) lapleps => this%laplacian%of_r
        !
        IF (deriv >= 3) THEN
            dsurface => this%dsurface%of_r
            !
            IF (this%solvent_aware) THEN
                hessian => this%hessian
            ELSE
                ALLOCATE (hessian)
                !
                CALL hessian%init(this%density%cell)
                !
            END IF
            !
        END IF
        !
        SELECT CASE (this%derivatives%type_)
            !
        CASE ('fft')
            !
            IF (deriv == 1 .OR. deriv == 2) &
                CALL this%derivatives%gradient(this%scaled, this%gradient)
            !
            IF (deriv == 2) CALL this%derivatives%laplacian(this%scaled, this%laplacian)
            !
            IF (deriv == 3) &
                CALL this%calc_dsurface(this%scaled, this%gradient, this%laplacian, &
                                        hessian, this%dsurface)
            !
        CASE ('chain', 'fd')
            !
            IF (deriv == 1 .OR. deriv == 2) &
                CALL this%derivatives%gradient(this%density, this%gradient)
            !
            IF (deriv == 2) CALL this%derivatives%laplacian(this%density, this%laplacian)
            !
            IF (deriv == 3) THEN
                !
                CALL this%calc_dsurface(this%density, this%gradient, this%laplacian, &
                                        hessian, this%dsurface)
                !
                IF (this%solvent_aware) THEN
                    !
                    DO ipol = 1, 3
                        !
                        DO jpol = 1, 3
                            !
                            hessian%of_r(ipol, jpol, :) = &
                                hessian%of_r(ipol, jpol, :) * deps(:) + &
                                gradeps(ipol, :) * gradeps(jpol, :) * d2eps(:)
                            !
                        END DO
                        !
                    END DO
                    !
                END IF
                !
            END IF
            !
            IF (deriv > 1) &
                lapleps(:) = lapleps(:) * deps(:) + &
                             (gradeps(1, :)**2 + gradeps(2, :)**2 + &
                              gradeps(3, :)**2) * d2eps(:)
            !
            IF (deriv >= 1) THEN
                !
                IF (this%derivatives%type_ == 'chain') THEN
                    !
                    DO ipol = 1, 3
                        gradeps(ipol, :) = gradeps(ipol, :) * deps(:)
                    END DO
                    !
                    ! ELSE IF (this%derivatives%type_ == 'fd') THEN
                    !     CALL fd%gradient(this%scaled, this%gradient)
                END IF
                !
            END IF
            !
        END SELECT
        !
        !--------------------------------------------------------------------------------
        ! Final updates
        !
        this%volume = this%scaled%integrate()
        !
        IF (deriv >= 1) THEN
            !
            CALL this%gradient%update_modulus()
            !
            this%surface = this%gradient%modulus%integrate()
        END IF
        !
        IF (deriv >= 3 .AND. .NOT. this%solvent_aware) CALL hessian%destroy()
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE boundary_of_density
    !------------------------------------------------------------------------------------
    !>
    !! @brief Updates boundary object using function objects
    !!
    !! Calculates the dielectric constant as a function of the charge
    !! density, and derivatives of the dielectric constant with respect
    !! to the charge density. Also updates the volume and surface
    !! components. This function is implemented for the soft-spheres
    !! interface model. It expects a series of environ_functions of
    !! dimension equal to nsoft_spheres.
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE boundary_of_functions(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_boundary), TARGET, INTENT(INOUT) :: this
        !
        INTEGER, POINTER :: nnr, ir_end, deriv
        TYPE(environ_cell), POINTER :: cell
        !
        INTEGER :: i
        !
        INTEGER, POINTER :: nsoft_spheres
        !
        TYPE(environ_density), ALLOCATABLE :: local(:)
        TYPE(environ_gradient), ALLOCATABLE :: gradlocal(:)
        TYPE(environ_density), ALLOCATABLE :: lapllocal(:)
        TYPE(environ_hessian), ALLOCATABLE :: hesslocal(:)
        TYPE(environ_hessian), POINTER :: hessian
        !
        CHARACTER(LEN=80) :: label
        !
        CHARACTER(LEN=80) :: sub_name = 'boundary_of_functions'
        !
        !--------------------------------------------------------------------------------
        !
        cell => this%scaled%cell
        nnr => cell%nnr
        ir_end => cell%ir_end
        nsoft_spheres => this%ions%number
        !
        ALLOCATE (local(nsoft_spheres))
        !
        !--------------------------------------------------------------------------------
        ! Compute soft spheres and generate boundary
        !
        this%scaled%of_r = 1.D0
        !
        DO i = 1, nsoft_spheres
            !
            CALL local(i)%init(cell)
            !
            CALL this%soft_spheres(i)%density(local(i), .FALSE.)
            !
            this%scaled%of_r = this%scaled%of_r * local(i)%of_r
        END DO
        !
        !--------------------------------------------------------------------------------
        ! Generate boundary derivatives, if needed
        !
        deriv => this%deriv
        !
        IF (deriv == 3) THEN
            !
            IF (this%solvent_aware) THEN
                hessian => this%hessian
                hessian%of_r = 0.D0
            ELSE
                ALLOCATE (hessian)
                !
                CALL hessian%init(cell)
                !
            END IF
            !
        END IF
        !
        SELECT CASE (this%derivatives%type_)
            !
        CASE ('fft')
            !
            IF (deriv == 1 .OR. deriv == 2) &
                CALL this%derivatives%gradient(this%scaled, this%gradient)
            !
            IF (deriv == 2) CALL this%derivatives%laplacian(this%scaled, this%laplacian)
            !
            IF (deriv == 3) &
                CALL this%calc_dsurface(this%scaled, this%gradient, this%laplacian, &
                                        hessian, this%dsurface)
            !
        CASE ('highmem')
            !
            IF (deriv >= 1) ALLOCATE (gradlocal(nsoft_spheres))
            !
            IF (deriv == 2) ALLOCATE (lapllocal(nsoft_spheres))
            !
            IF (deriv == 3) ALLOCATE (hesslocal(nsoft_spheres))
            !
            !----------------------------------------------------------------------------
            ! Compute and temporarily store soft spheres derivatives
            !
            DO i = 1, nsoft_spheres
                !
                IF (deriv >= 1) CALL gradlocal(i)%init(cell)
                !
                IF (deriv == 2) CALL lapllocal(i)%init(cell)
                !
                IF (deriv == 3) CALL hesslocal(i)%init(cell)
                !
                IF (deriv >= 1) CALL this%soft_spheres(i)%gradient(gradlocal(i), .FALSE.)
                !
                IF (deriv == 2) CALL this%soft_spheres(i)%laplacian(lapllocal(i), .FALSE.)
                !
                IF (deriv == 3) CALL this%soft_spheres(i)%hessian(hesslocal(i), .FALSE.)
                !
            END DO
            !
            IF (deriv == 1 .OR. deriv == 2) &
                CALL calc_gradient_of_boundary_highmem(nsoft_spheres, local, &
                                                       gradlocal, this%gradient)
            !
            IF (deriv == 2) &
                CALL calc_laplacian_of_boundary_highmem(nsoft_spheres, local, &
                                                        gradlocal, lapllocal, &
                                                        this%laplacian)
            !
            IF (deriv == 3) &
                CALL calc_dsurface_of_boundary_highmem(nsoft_spheres, local, &
                                                       gradlocal, hesslocal, &
                                                       this%gradient, &
                                                       this%laplacian, hessian, &
                                                       this%dsurface)
            !
            DO i = 1, nsoft_spheres
                !
                IF (deriv >= 1) CALL gradlocal(i)%destroy()
                !
                IF (deriv == 2) CALL lapllocal(i)%destroy()
                !
                IF (deriv == 3) CALL hesslocal(i)%destroy()
                !
            END DO
            !
            IF (deriv >= 1) DEALLOCATE (gradlocal)
            !
            IF (deriv == 2) DEALLOCATE (lapllocal)
            !
            IF (deriv == 3) DEALLOCATE (hesslocal)
            !
        CASE ('lowmem')
            !
            IF (deriv >= 1) ALLOCATE (gradlocal(nsoft_spheres))
            !
            IF (deriv == 2) ALLOCATE (lapllocal(nsoft_spheres))
            !
            IF (deriv == 3) ALLOCATE (hesslocal(nsoft_spheres))
            !
            !----------------------------------------------------------------------------
            ! Compute and temporarily store soft spheres derivatives
            !
            DO i = 1, nsoft_spheres
                !
                IF (deriv >= 1) CALL gradlocal(i)%init(cell)
                !
                IF (deriv == 2) CALL lapllocal(i)%init(cell)
                !
                IF (deriv == 3) CALL hesslocal(i)%init(cell)
                !
                IF (deriv >= 1) CALL this%soft_spheres(i)%gradient(gradlocal(i), .FALSE.)
                !
                IF (deriv == 2) CALL this%soft_spheres(i)%laplacian(lapllocal(i), .FALSE.)
                !
                IF (deriv == 3) CALL this%soft_spheres(i)%hessian(hesslocal(i), .FALSE.)
                !
            END DO
            !
            IF (deriv >= 1) &
                CALL calc_gradient_of_boundary_lowmem(nsoft_spheres, local, &
                                                      gradlocal, this%scaled, &
                                                      this%gradient)
            !
            IF (deriv == 2) &
                CALL calc_laplacian_of_boundary_lowmem(nsoft_spheres, local, &
                                                       gradlocal, lapllocal, &
                                                       this%scaled, &
                                                       this%gradient, &
                                                       this%laplacian)
            !
            IF (deriv == 3) &
                CALL calc_dsurface_of_boundary_lowmem(nsoft_spheres, local, &
                                                      gradlocal, hesslocal, &
                                                      this%gradient, &
                                                      this%laplacian, hessian, &
                                                      this%scaled, &
                                                      this%dsurface)
            !
            DO i = 1, nsoft_spheres
                !
                IF (deriv >= 1) CALL gradlocal(i)%destroy()
                !
                IF (deriv == 2) CALL lapllocal(i)%destroy()
                !
                IF (deriv == 3) CALL hesslocal(i)%destroy()
                !
            END DO
            !
            IF (deriv >= 1) DEALLOCATE (gradlocal)
            !
            IF (deriv == 2) DEALLOCATE (lapllocal)
            !
            IF (deriv == 3) DEALLOCATE (hesslocal)
            !
        END SELECT
        !
        !--------------------------------------------------------------------------------
        ! Final updates
        !
        this%scaled%of_r = 1.D0 - this%scaled%of_r
        this%volume = this%scaled%integrate()
        !
        IF (deriv >= 1) THEN
            this%gradient%of_r = -this%gradient%of_r
            !
            CALL this%gradient%update_modulus()
            !
            this%surface = this%gradient%modulus%integrate()
            !
            IF (deriv >= 2) this%laplacian%of_r = -this%laplacian%of_r
            !
            IF (deriv == 3) THEN
                this%dsurface%of_r = -this%dsurface%of_r
                !
                IF (this%solvent_aware) THEN
                    this%hessian%of_r = -this%hessian%of_r
                ELSE
                    !
                    CALL hessian%destroy()
                    !
                    DEALLOCATE (hessian)
                END IF
                !
            END IF
            !
        END IF
        !
        DO i = 1, nsoft_spheres
            CALL local(i)%destroy()
        END DO
        !
        DEALLOCATE (local)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE boundary_of_functions
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
    SUBROUTINE boundary_of_system(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_boundary), TARGET, INTENT(INOUT) :: this
        !
        INTEGER, POINTER :: nnr, ir_end, deriv
        TYPE(environ_cell), POINTER :: cell
        !
        TYPE(environ_hessian), POINTER :: hesslocal
        !
        INTEGER :: i
        !
        CHARACTER(LEN=80) :: sub_name = 'boundary_of_system'
        !
        !--------------------------------------------------------------------------------
        !
        cell => this%scaled%cell
        nnr => cell%nnr
        ir_end => cell%ir_end
        !
        CALL this%simple%density(this%scaled, .TRUE.)
        ! compute soft spheres and generate boundary
        !
        !--------------------------------------------------------------------------------
        ! Generate boundary derivatives, if needed
        !
        deriv => this%deriv
        !
        IF (deriv >= 3) THEN
            !
            IF (this%solvent_aware) THEN
                hesslocal => this%hessian
            ELSE
                ALLOCATE (hesslocal)
                !
                CALL hesslocal%init(cell)
                !
            END IF
            !
        END IF
        !
        SELECT CASE (this%derivatives%type_)
            !
        CASE ('fft')
            !
            IF (deriv == 1 .OR. deriv == 2) &
                CALL this%derivatives%gradient(this%scaled, this%gradient)
            !
            IF (deriv == 2) CALL this%derivatives%laplacian(this%scaled, this%laplacian)
            !
            IF (deriv == 3) &
                CALL this%calc_dsurface(this%scaled, this%gradient, this%laplacian, &
                                        hesslocal, this%dsurface)
            !
        CASE ('chain')
            !
            IF (deriv >= 1) CALL this%simple%gradient(this%gradient, .TRUE.)
            !
            IF (deriv >= 2) CALL this%simple%laplacian(this%laplacian, .TRUE.)
            !
            IF (deriv >= 3) THEN
                !
                CALL this%simple%hessian(hesslocal, .TRUE.)
                !
                CALL calc_dsurface_no_pre(nnr, ir_end, this%gradient%of_r, &
                                          hesslocal%of_r, this%dsurface%of_r)
                !
            END IF
            !
        END SELECT
        !
        IF (deriv >= 3) THEN
            !
            IF (.NOT. this%solvent_aware) THEN
                !
                CALL hesslocal%destroy()
                !
                DEALLOCATE (hesslocal)
            END IF
            !
        END IF
        !
        this%volume = this%scaled%integrate()
        !
        IF (deriv >= 1) THEN
            !
            CALL this%gradient%update_modulus()
            !
            this%surface = this%gradient%modulus%integrate()
        END IF
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE boundary_of_system
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_dboundary_dions(this, index, partial)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: index
        CLASS(environ_boundary), INTENT(IN), TARGET :: this
        !
        TYPE(environ_gradient), INTENT(INOUT) :: partial
        !
        REAL(DP), PARAMETER :: tolspuriousforce = 1.D-5
        !
        INTEGER, POINTER :: number
        TYPE(environ_cell), POINTER :: cell
        !
        INTEGER :: i, ipol
        REAL(DP) :: spurious_force
        TYPE(environ_density) :: local
        !
        CHARACTER(LEN=80) :: sub_name = 'calc_dboundary_dions'
        !
        !--------------------------------------------------------------------------------
        !
        IF (this%mode == 'electronic') RETURN
        ! exit if boundary is only defined on electronic density
        !
        cell => partial%cell
        !
        IF (this%need_ions) THEN
            number => this%ions%number
        ELSE IF (this%need_system) THEN
            number => this%system%ions%number
        ELSE
            CALL env_errore(sub_name, 'Missing details of ions', 1)
        END IF
        !
        IF (index > number) &
            CALL env_errore(sub_name, 'Index greater than number of ions', 1)
        !
        IF (index <= 0) &
            CALL env_errore(sub_name, 'Index of ion is zero or lower', 1)
        !
        IF (this%mode == 'ionic' .AND. &
            .NOT. ALLOCATED(this%soft_spheres)) &
            CALL env_errore(sub_name, 'Missing details of ionic boundary', 1)
        !
        IF (this%mode == 'full' .AND. .NOT. ALLOCATED(this%ions%core_electrons)) &
            CALL env_errore(sub_name, 'Missing details of core electrons', 1)
        !
        IF (this%mode == 'full' .AND. .NOT. ASSOCIATED(this%dscaled%cell, cell)) &
            CALL env_errore(sub_name, 'Mismatch or unassociated boundary derivative', 1)
        !
        IF (this%mode == 'ionic' .OR. this%mode == 'fa-ionic') THEN
            !
            CALL this%soft_spheres(index)%gradient(partial, .TRUE.)
            !
            CALL local%init(cell)
            !
            DO i = 1, number
                !
                IF (i == index) CYCLE
                !
                CALL this%soft_spheres(i)%density(local, .TRUE.)
                !
                DO ipol = 1, 3
                    partial%of_r(ipol, :) = partial%of_r(ipol, :) * local%of_r(:)
                END DO
                !
            END DO
            !
            CALL local%destroy()
            !
        ELSE IF (this%mode == 'full') THEN
            !
            CALL this%ions%core_electrons(index)%gradient(partial, .TRUE.)
            !
            DO ipol = 1, 3
                partial%of_r(ipol, :) = -partial%of_r(ipol, :) * this%dscaled%of_r(:)
            END DO
            !
            CALL partial%update_modulus()
            !
            spurious_force = partial%modulus%integrate()
            !
            IF (spurious_force > tolspuriousforce .AND. ionode) &
                WRITE (program_unit, 1000) index, spurious_force
            !
1000        FORMAT(' WARNING: Unphysical forces due to core electrons are non-negligible ', /, &
                   ' atom type ', I3, ' is subject to a spurious force of ', F12.6)
            !
        ELSE IF (this%mode == 'system') THEN
            !
            ! PROBABLY THERE IS A UNIFORM CONTRIBUTION TO THE FORCES
            ! WHICH SHOULD ONLY AFFECT THE COM OF THE SYSTEM, POSSIBLY NEED TO ADD
            ! A CHECK ON ATOMS THAT BELONG TO THE SYSTEM
            !
            partial%of_r = 0.D0
        END IF
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_dboundary_dions
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE invert_boundary(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_boundary), INTENT(INOUT) :: this
        !
        !--------------------------------------------------------------------------------
        !
        this%scaled%of_r = 1.D0 - this%scaled%of_r
        !
        this%volume = this%scaled%integrate()
        !
        IF (this%deriv >= 1) this%gradient%of_r = -this%gradient%of_r
        !
        IF (this%deriv >= 2) this%laplacian%of_r = -this%laplacian%of_r
        !
        IF (this%deriv >= 3) THEN
            this%dsurface%of_r = -this%dsurface%of_r
            !
            IF (this%solvent_aware) this%hessian%of_r = -this%hessian%of_r
            !
        END IF
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE invert_boundary
    !------------------------------------------------------------------------------------
    !>
    !! Fill voids of the continuum interface that are too small
    !! to fit a solvent molecule
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE solvent_aware_boundary(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_boundary), INTENT(INOUT), TARGET :: this
        !
        INTEGER, POINTER :: nnr, ir_end, deriv
        REAL(DP), POINTER :: thr, spr
        TYPE(environ_cell), POINTER :: cell
        !
        INTEGER :: ir, ipol, jpol
        TYPE(environ_density) :: filled_fraction
        TYPE(environ_density) :: d2filling
        !
        TYPE(environ_density) :: local
        TYPE(environ_gradient) :: gradlocal
        TYPE(environ_density) :: lapllocal
        TYPE(environ_hessian) :: hesslocal
        !
        CHARACTER(LEN=80) :: label
        !
        CHARACTER(LEN=80) :: sub_name = 'solvent_aware_boundary'
        !
        !--------------------------------------------------------------------------------
        !
        cell => this%scaled%cell
        nnr => this%scaled%cell%nnr
        ir_end => this%scaled%cell%ir_end
        deriv => this%deriv
        !
        thr => this%filling_threshold
        spr => this%filling_spread
        !
        CALL filled_fraction%init(cell)
        !
        IF (deriv >= 2 .AND. this%derivatives%type_ /= 'fft') &
            CALL d2filling%init(cell)
        !
        !--------------------------------------------------------------------------------
        ! Step 0: save local interface function for later use
        !
        this%local%of_r = this%scaled%of_r
        !
        !--------------------------------------------------------------------------------
        ! Step 1: compute the convolution function, this may be made moved out of here
        !
        CALL this%solvent_probe%density(this%probe, .TRUE.)
        !
        this%probe%of_r = this%probe%of_r / this%probe%integrate()
        !
        !--------------------------------------------------------------------------------
        ! Step 2: compute filled fraction, i.e. convolution of local boundary with probe
        !
        CALL this%derivatives%convolution(this%local, this%probe, filled_fraction)
        !
        !--------------------------------------------------------------------------------
        ! Step 3: compute the filling function and its derivative
        !
        this%filling%of_r = 0.D0
        this%dfilling%of_r = 0.D0
        !
        DO ir = 1, ir_end
            this%filling%of_r(ir) = 1.D0 - sfunct2(filled_fraction%of_r(ir), thr, spr)
            this%dfilling%of_r(ir) = -dsfunct2(filled_fraction%of_r(ir), thr, spr)
            !
            IF (deriv >= 2 .AND. this%derivatives%type_ /= 'fft') &
                d2filling%of_r(ir) = -d2sfunct2(filled_fraction%of_r(ir), thr, spr)
            !
        END DO
        !
        !--------------------------------------------------------------------------------
        ! Step 4: compute solvent-aware interface
        !
        this%scaled%of_r = this%local%of_r + (1.D0 - this%local%of_r) * this%filling%of_r
        !
        !--------------------------------------------------------------------------------
        ! Step 5: compute boundary derivatives, if needed
        !
        SELECT CASE (this%derivatives%type_)
            !
        CASE ('fft')
            !
            IF (deriv == 1 .OR. deriv == 2) &
                CALL this%derivatives%gradient(this%scaled, this%gradient)
            !
            IF (deriv == 2) CALL this%derivatives%laplacian(this%scaled, this%laplacian)

            IF (deriv == 3) &
                CALL this%calc_dsurface(this%scaled, this%gradient, this%laplacian, &
                                        this%hessian, this%dsurface)
            !
        CASE ('chain', 'highmem')
            !
            !----------------------------------------------------------------------------
            ! Allocate local fields for derivatives of convolution
            !
            IF (deriv >= 1) CALL gradlocal%init(cell)
            !
            IF (deriv >= 2) CALL lapllocal%init(cell)
            !
            IF (deriv >= 3) CALL hesslocal%init(cell)
            !
            !----------------------------------------------------------------------------
            ! Compute derivative of convolution with probe
            !
            IF (deriv > 1) &
                CALL this%convolution_deriv(deriv, gradlocal, lapllocal, hesslocal)
            !
            !----------------------------------------------------------------------------
            ! Update derivatives of interface function in reverse order
            !
            IF (deriv >= 3) THEN
                !
                DO ipol = 1, 3
                    !
                    DO jpol = 1, 3
                        !
                        this%hessian%of_r(ipol, jpol, :) = &
                            this%hessian%of_r(ipol, jpol, :) * &
                            (1.D0 - this%filling%of_r) - this%dfilling%of_r * &
                            (this%gradient%of_r(ipol, :) * &
                             gradlocal%of_r(jpol, :) + &
                             this%gradient%of_r(jpol, :) * &
                             gradlocal%of_r(ipol, :)) + &
                            (1.D0 - this%local%of_r) * &
                            (d2filling%of_r * gradlocal%of_r(ipol, :) * &
                             gradlocal%of_r(jpol, :) + &
                             this%dfilling%of_r * hesslocal%of_r(ipol, jpol, :))
                        !
                    END DO
                    !
                END DO
                !
                CALL hesslocal%destroy()
                !
            END IF
            !
            IF (deriv >= 2) THEN
                !
                CALL local%init(cell)
                !
                CALL this%gradient%scalar_product(gradlocal, local)
                !
                this%laplacian%of_r = &
                    this%laplacian%of_r * (1.D0 - this%filling%of_r) - &
                    2.D0 * local%of_r * this%dfilling%of_r + &
                    (1.D0 - this%local%of_r) * &
                    (d2filling%of_r * gradlocal%modulus%of_r**2 + &
                     this%dfilling%of_r * lapllocal%of_r)
                !
                CALL local%destroy()
                !
                CALL lapllocal%destroy()
                !
                CALL d2filling%destroy()
                !
            END IF
            !
            IF (deriv >= 1) THEN
                !
                DO ipol = 1, 3
                    !
                    this%gradient%of_r(ipol, :) = &
                        this%gradient%of_r(ipol, :) * &
                        (1.D0 - this%filling%of_r(:)) + &
                        gradlocal%of_r(ipol, :) * &
                        (1.D0 - this%local%of_r(:)) * &
                        this%dfilling%of_r(:)
                    !
                END DO
                !
                CALL gradlocal%destroy()
                !
            END IF
            !
            !----------------------------------------------------------------------------
            ! Recompute dsurface, if needed
            !
            IF (deriv >= 3) THEN
                !
                CALL calc_dsurface_no_pre(nnr, ir_end, this%gradient%of_r, &
                                          this%hessian%of_r, this%dsurface%of_r)
                !
            END IF
            !
        END SELECT
        !
        !--------------------------------------------------------------------------------
        ! Final updates
        !
        this%volume = this%scaled%integrate()
        !
        IF (deriv >= 1) THEN
            !
            CALL this%gradient%update_modulus()
            !
            this%surface = this%gradient%modulus%integrate()
        END IF
        !
        CALL filled_fraction%destroy()
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE solvent_aware_boundary
    !------------------------------------------------------------------------------------
    !>
    !! @brief Compute the functional derivative of the energy w.r.t the boundary
    !!
    !! @param[out]  de_dboundary  the computed derivative
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_solvent_aware_de_dboundary(this, de_dboundary)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_boundary), INTENT(IN), TARGET :: this
        !
        TYPE(environ_density), INTENT(INOUT) :: de_dboundary
        !
        INTEGER, POINTER :: nnr, ir_end
        REAL(DP), POINTER :: thr, spr
        TYPE(environ_cell), POINTER :: cell
        !
        TYPE(environ_density) :: local
        !
        CHARACTER(LEN=80) :: sub_name = 'calc_solvent_aware_de_dboundary'
        !
        !--------------------------------------------------------------------------------
        !
        cell => this%scaled%cell
        nnr => this%scaled%cell%nnr
        !
        CALL local%init(cell)
        !
        !--------------------------------------------------------------------------------
        ! Step 1: compute (1-s)*de_dboudary*dfilling
        !
        local%of_r = (1.D0 - this%local%of_r) * de_dboundary%of_r * this%dfilling%of_r
        !
        !--------------------------------------------------------------------------------
        ! Step 2: compute convolution with the probe function
        !
        CALL this%derivatives%convolution(this%probe, local, local)
        !
        !--------------------------------------------------------------------------------
        ! Step 3: update the functional derivative of the energy wrt boundary
        !
        de_dboundary%of_r = de_dboundary%of_r * (1.D0 - this%filling%of_r) + local%of_r
        !
        CALL local%destroy()
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_solvent_aware_de_dboundary
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE compute_convolution_deriv(this, deriv, grad, lapl, hess)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: deriv
        CLASS(environ_boundary), INTENT(IN) :: this
        !
        TYPE(environ_gradient), INTENT(INOUT) :: grad
        TYPE(environ_density), INTENT(INOUT) :: lapl
        TYPE(environ_hessian), INTENT(INOUT) :: hess
        !
        CHARACTER(LEN=80) :: sub_name = 'compute_convolution_deriv'
        !
        !--------------------------------------------------------------------------------
        !
        IF (deriv <= 0) RETURN
        !
        IF (deriv >= 1) THEN
            !
            CALL this%derivatives%convolution(this%probe, this%gradient, grad)
            !
            CALL grad%update_modulus()
            !
        END IF
        !
        IF (deriv >= 2) &
            CALL this%derivatives%convolution(this%probe, this%laplacian, lapl)
        !
        IF (deriv >= 3) CALL this%derivatives%convolution(this%probe, this%hessian, hess)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE compute_convolution_deriv
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_dsurface(this, dens, grad, lapl, hess, dsurface)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_density), INTENT(IN) :: dens
        !
        CLASS(environ_boundary), INTENT(INOUT) :: this
        TYPE(environ_gradient), INTENT(INOUT) :: grad
        TYPE(environ_density), INTENT(INOUT) :: lapl
        TYPE(environ_hessian), INTENT(INOUT) :: hess
        TYPE(environ_density), INTENT(INOUT) :: dsurface
        !
        CHARACTER(LEN=80) :: sub_name = 'calc_dsurface'
        !
        !--------------------------------------------------------------------------------
        !
        CALL this%derivatives%hessian(dens, grad, hess)
        !
        lapl%of_r(:) = hess%of_r(1, 1, :) + hess%of_r(2, 2, :) + hess%of_r(3, 3, :)
        !
        CALL calc_dsurface_no_pre(dens%cell%nnr, dens%cell%ir_end, grad%of_r, &
                                  hess%of_r, dsurface%of_r)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_dsurface
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !                                   OUTPUT METHODS
    !
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE print_environ_boundary(this, local_verbose, local_depth)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_boundary), INTENT(IN) :: this
        INTEGER, INTENT(IN), OPTIONAL :: local_verbose
        INTEGER, INTENT(IN), OPTIONAL :: local_depth
        !
        INTEGER :: verbosity, passed_verbosity, passed_depth
        !
        CHARACTER(LEN=80) :: sub_name = 'print_environ_boundary'
        !
        !--------------------------------------------------------------------------------
        !
        IF (verbose == 0) RETURN
        !
        IF (PRESENT(local_verbose)) THEN
            verbosity = verbose + local_verbose
        ELSE
            verbosity = verbose
        END IF
        !
        IF (verbosity == 0) RETURN
        !
        IF (PRESENT(local_depth)) THEN
            passed_verbosity = verbosity - verbose - local_depth
            passed_depth = local_depth
        ELSE
            passed_verbosity = verbosity - verbose - depth
            passed_depth = depth
        END IF
        !
        IF (verbosity >= 1) THEN
            !
            IF (ionode) THEN
                !
                IF (verbosity >= verbose) THEN ! header
                    WRITE (environ_unit, 1100)
                ELSE
                    !
                    CALL env_block_divider(verbosity)
                    !
                    WRITE (environ_unit, 1101)
                END IF
                !
                WRITE (environ_unit, 1102) this%label, this%mode
                !
            END IF
            !
            IF (this%need_electrons) THEN
                !
                IF (ionode) THEN
                    !
                    WRITE (environ_unit, 1103) this%b_type
                    !
                    SELECT CASE (this%b_type)
                        !
                    CASE (0)
                        WRITE (environ_unit, 1104) this%rhozero, this%tbeta
                        !
                    CASE (1)
                        WRITE (environ_unit, 1105) this%rhomax, this%rhomin
                        !
                        IF (verbosity >= 3) WRITE (environ_unit, 1106) this%fact
                        !
                    CASE (2)
                        WRITE (environ_unit, 1107) this%rhomax, this%rhomin
                        !
                    END SELECT
                    !
                END IF
                !
                IF (verbosity >= 4) THEN
                    !
                    CALL this%density%printout(passed_verbosity, passed_depth)
                    !
                    IF (ionode .AND. this%need_ions) WRITE (environ_unit, 1108)
                    !
                END IF
                !
                IF (verbosity >= 5) THEN
                    !
                    CALL this%dscaled%printout(passed_verbosity, passed_depth)
                    !
                    CALL this%d2scaled%printout(passed_verbosity, passed_depth)
                    !
                END IF
                !
            ELSE IF (this%need_ions) THEN
                !
                IF (ionode) WRITE (environ_unit, 1109) this%alpha, this%softness
                !
                IF (verbosity >= 3) &
                    CALL print_environ_functions(this%soft_spheres, this%ions%number, &
                                                 passed_verbosity, passed_depth)
                !
            ELSE IF (ionode .AND. this%need_system) THEN
                !
                WRITE (environ_unit, 1110) &
                    this%simple%pos, this%simple%width, &
                    this%simple%spread, this%simple%dim, &
                    this%simple%axis
                !
            END IF
            !
            IF (ionode) THEN
                WRITE (environ_unit, 1111) this%volume
                !
                IF (this%deriv >= 1) WRITE (environ_unit, 1112) this%surface
                !
            END IF
            !
            IF (verbosity >= 4) CALL this%scaled%printout(passed_verbosity, passed_depth)
            !
            IF (this%solvent_aware) THEN
                !
                IF (ionode) &
                    WRITE (environ_unit, 1113) &
                    this%filling_threshold, this%filling_spread, &
                    this%solvent_probe%width, this%solvent_probe%spread
                !
                IF (verbosity >= 4) THEN
                    !
                    CALL this%local%printout(passed_verbosity, passed_depth)
                    !
                    CALL this%filling%printout(passed_verbosity, passed_depth)
                    !
                END IF
                !
                IF (verbosity >= 5) THEN
                    !
                    CALL this%dfilling%printout(passed_verbosity, passed_depth)
                    !
                    CALL this%probe%printout(passed_verbosity, passed_depth)
                    !
                END IF
                !
            END IF
            !
            IF (verbosity >= 5) THEN
                !
                IF (this%deriv >= 1) &
                    CALL this%gradient%printout(passed_verbosity, passed_depth)
                !
                IF (this%deriv >= 2) &
                    CALL this%laplacian%printout(passed_verbosity, passed_depth)
                !
                IF (this%deriv == 3) &
                    CALL this%dsurface%printout(passed_verbosity, passed_depth)
                !
            END IF
            !
            IF (verbosity < verbose) CALL env_block_divider(verbosity)
            !
        END IF
        !
        FLUSH (environ_unit)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
        !
1100    FORMAT(/, 4('%'), ' BOUNDARY ', 66('%'))
1101    FORMAT(/, ' BOUNDARY', /, '========')
        !
1102    FORMAT(/, ' boundary label             = ', A20, /, &
                ' boundary mode              = ', A20)
        !
1103    FORMAT(/, ' boundary is built as a type-', I1, ' function of a smooth density')
        !
1104    FORMAT(/, ' using the Fattebert-Gygi function:', /, &
                ' rhozero                    = ', F14.7, /, &
                ' 2*beta                     = ', F14.7)
        !
1105    FORMAT(/, ' using the optimal SCCS function:', /, &
                ' rhomax                     = ', F14.7, /, &
                ' rhomin                     = ', F14.7)
        !
1106    FORMAT(' log(rhomax/rhomin)         = ', F14.7)
        !
1107    FORMAT(/, ' using the modified SCCS function:', /, &
                ' rhomax                     = ', F14.7, /, &
                ' rhomin                     = ', F14.7)
        !
1108    FORMAT(/, ' adding fictitious core-electrons')
        !
1109    FORMAT(/, ' boundary is built from soft-spheres centered on ionic positions:', /, &
                ' solvent-dependent scaling  = ', F14.7, /, &
                ' softness parameter         = ', F14.7)
        !
1110    FORMAT(/, ' boundary is built as an analytic function centered on system position:', /, &
                ' center of the boundary     = ', 3F14.7, /, &
                ' distance from the center   = ', F14.7, /, &
                ' spread of the interface    = ', F14.7, /, &
                ' dimensionality             = ', I2, /, &
                ' axis                       = ', I2)
        !
1111    FORMAT(/, ' volume of the QM region    = ', F14.7)
        !
1112    FORMAT(/, ' surface of the QM region   = ', F14.7)
        !
1113    FORMAT(/, ' using solvent-aware boundary:', /, &
                ' filling threshold          = ', F14.7, /, &
                ' filling spread             = ', F14.7, /, &
                ' solvent radius x rad scale = ', F14.7, /, &
                ' spread of solvent probe    = ', F14.7)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE print_environ_boundary
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE class_boundary
!----------------------------------------------------------------------------------------
