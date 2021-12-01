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
    USE class_io, ONLY: io
    !
    USE environ_param, ONLY: DP, e2, sqrtpi, tpi
    !
    USE class_cell
    USE class_density
    USE class_function
    USE class_function_erfc
    USE class_functions
    USE class_gradient
    USE class_hessian
    !
    USE class_core_container
    !
    USE class_core_fft
    !
    USE class_electrons
    USE class_ions
    USE class_system
    !
    USE tools_math, ONLY: environ_erfc
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
        !
        !--------------------------------------------------------------------------------
        ! Parameters for the electrons-dependent interface
        !
        LOGICAL :: need_electrons = .FALSE.
        TYPE(environ_electrons), POINTER :: electrons => NULL()
        !
        !--------------------------------------------------------------------------------
        ! Parameters for the ions-dependent interface
        !
        LOGICAL :: need_ions = .FALSE.
        TYPE(environ_ions), POINTER :: ions => NULL()
        !
        !--------------------------------------------------------------------------------
        ! Parameters for the system-dependent interface
        !
        LOGICAL :: need_system = .FALSE.
        TYPE(environ_system), POINTER :: system => NULL()
        !
        !--------------------------------------------------------------------------------
        !
        TYPE(environ_density) :: scaled ! scaled switching function of interface
        ! varying from 1 (QM region) to 0 (environment region)
        !
        INTEGER :: deriv = 0
        TYPE(environ_gradient) :: gradient
        TYPE(environ_density) :: laplacian
        TYPE(environ_hessian) :: hessian
        TYPE(environ_density) :: dsurface
        !
        TYPE(core_container), POINTER :: cores => NULL()
        !
        CHARACTER(LEN=80) :: derivatives_method
        !
        !--------------------------------------------------------------------------------
        ! Global properties of the boundary
        !
        REAL(DP) :: volume = 0.D0
        REAL(DP) :: surface = 0.D0
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
        CLASS(environ_function), ALLOCATABLE :: soft_spheres(:)
        !
        !--------------------------------------------------------------------------------
        !
        TYPE(environ_function_erfc) :: simple ! components needed for boundary of system
        !
        !--------------------------------------------------------------------------------
        ! Components needed for solvent-aware boundary
        !
        LOGICAL :: solvent_aware = .FALSE.
        TYPE(environ_function_erfc) :: solvent_probe
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
        LOGICAL :: field_aware = .FALSE.
        REAL(DP) :: field_factor, charge_asymmetry, field_max, field_min

        TYPE(environ_density) :: normal_field
        REAL(DP), ALLOCATABLE :: ion_field(:)
        CLASS(environ_function), ALLOCATABLE :: local_spheres(:)
        TYPE(environ_density), ALLOCATABLE :: dion_field_drho(:)
        REAL(DP), ALLOCATABLE :: partial_of_ion_field(:, :, :)
        !
        !--------------------------------------------------------------------------------
    CONTAINS
        !--------------------------------------------------------------------------------
        !
        PROCEDURE, PRIVATE :: create => create_environ_boundary
        PROCEDURE :: init => init_environ_boundary
        PROCEDURE :: copy => copy_environ_boundary
        PROCEDURE :: update => update_environ_boundary
        PROCEDURE :: destroy => destroy_environ_boundary
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
        PROCEDURE :: boundary_of_density
        PROCEDURE :: boundary_of_functions
        PROCEDURE :: boundary_of_system
        !
        PROCEDURE :: convolution_deriv => compute_convolution_deriv
        PROCEDURE :: solvent_aware_boundary
        PROCEDURE :: calc_dsurface ! #TODO do we need this?
        PROCEDURE :: invert => invert_boundary
        !
        PROCEDURE, PRIVATE :: set_soft_spheres
        PROCEDURE, PRIVATE :: update_soft_spheres
        !
        PROCEDURE :: printout => print_environ_boundary
        !
        !--------------------------------------------------------------------------------
    END TYPE environ_boundary
    !------------------------------------------------------------------------------------
    !
    INTEGER :: bound_tol = 1.D-60
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
        CLASS(environ_boundary), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: sub_name = 'create_environ_boundary'
        !
        !--------------------------------------------------------------------------------
        !
        IF (ASSOCIATED(this%electrons)) CALL io%create_error(sub_name)
        !
        IF (ASSOCIATED(this%ions)) CALL io%create_error(sub_name)
        !
        IF (ASSOCIATED(this%system)) CALL io%create_error(sub_name)
        !
        IF (ASSOCIATED(this%cores)) CALL io%create_error(sub_name)
        !
        IF (ALLOCATED(this%soft_spheres)) CALL io%create_error(sub_name)
        !
        IF (ALLOCATED(this%ion_field)) CALL io%create_error(sub_name)
        !
        IF (ALLOCATED(this%local_spheres)) CALL io%create_error(sub_name)
        !
        IF (ALLOCATED(this%dion_field_drho)) CALL io%create_error(sub_name)
        !
        IF (ALLOCATED(this%partial_of_ion_field)) CALL io%create_error(sub_name)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE create_environ_boundary
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_environ_boundary(this, need_gradient, need_laplacian, need_hessian, &
                                     mode, stype, rhomax, rhomin, tbeta, const, alpha, &
                                     softness, system_distance, system_spread, &
                                     solvent_radius, radial_scale, radial_spread, &
                                     filling_threshold, filling_spread, field_factor, &
                                     charge_asymmetry, field_max, field_min, electrons, &
                                     ions, system, cores, deriv_method, cell, label)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: stype
        CHARACTER(LEN=80), INTENT(IN) :: mode, deriv_method
        LOGICAL, INTENT(IN) :: need_gradient, need_laplacian, need_hessian
        !
        REAL(DP), INTENT(IN) :: rhomax, rhomin, tbeta, const, alpha, softness, &
                                system_distance, system_spread, solvent_radius, &
                                radial_scale, radial_spread, filling_threshold, &
                                filling_spread, field_factor, charge_asymmetry, &
                                field_max, field_min
        !
        TYPE(environ_electrons), TARGET, INTENT(IN) :: electrons
        TYPE(environ_ions), TARGET, INTENT(IN) :: ions
        TYPE(environ_system), TARGET, INTENT(IN) :: system
        TYPE(environ_cell), INTENT(IN) :: cell
        TYPE(core_container), TARGET, INTENT(IN) :: cores
        CHARACTER(LEN=80), INTENT(IN), OPTIONAL :: label
        !
        CLASS(environ_boundary), INTENT(INOUT) :: this
        !
        INTEGER :: i
        !
        CHARACTER(LEN=80) :: local_label
        !
        CHARACTER(LEN=80) :: sub_name = 'init_environ_boundary'
        !
        !--------------------------------------------------------------------------------
        !
        CALL this%create()
        !
        IF (need_hessian) THEN
            this%deriv = 3
        ELSE IF (need_laplacian) THEN
            this%deriv = 2
        ELSE IF (need_gradient) THEN
            this%deriv = 1
        END IF
        !
        !--------------------------------------------------------------------------------
        !
        this%mode = mode
        !
        this%need_electrons = mode == 'electronic' .OR. mode == 'full' .OR. &
                              mode == 'fa-ionic' .OR. mode == 'fa-electronic'
        !
        IF (this%need_electrons) this%electrons => electrons
        !
        this%need_ions = mode == 'ionic' .OR. mode == 'full' .OR. &
                         mode == 'fa-ionic' .OR. mode == 'fa-electronic'
        !
        IF (this%need_ions) this%ions => ions
        !
        this%need_system = mode == 'system'
        !
        IF (this%need_system) this%system => system
        !
        !--------------------------------------------------------------------------------
        !
        this%label = label
        this%b_type = stype
        this%rhomax = rhomax
        this%rhomin = rhomin
        this%fact = LOG(rhomax / rhomin)
        this%rhozero = (rhomax + rhomin) * 0.5_DP
        this%tbeta = tbeta
        this%deltarho = rhomax - rhomin
        !
        IF (const == 1.D0 .AND. this%need_electrons .AND. stype == 2) &
            CALL io%error(sub_name, &
                          'stype=2 boundary requires dielectric constant > 1', 1)
        !
        this%const = const
        this%alpha = alpha
        this%softness = softness
        !
        !--------------------------------------------------------------------------------
        !
        IF (this%need_system) &
            CALL this%simple%init(3, system%axis, system%dim, system_distance, &
                                  system_spread, 1.D0, system%pos)
        !
        !--------------------------------------------------------------------------------
        ! Derivatives
        !
        this%cores => cores
        this%derivatives_method = deriv_method
        !
        !--------------------------------------------------------------------------------
        ! Solvent aware
        !
        this%solvent_aware = solvent_radius > 0.D0
        !
        IF (this%solvent_aware) &
            CALL this%solvent_probe%init(2, 1, 0, solvent_radius * radial_scale, &
                                         radial_spread, 1.D0)
        !
        this%filling_threshold = filling_threshold
        this%filling_spread = filling_spread
        !
        !--------------------------------------------------------------------------------
        ! Field aware
        !
        this%field_aware = field_factor > 0.D0
        this%field_factor = field_factor
        this%charge_asymmetry = charge_asymmetry
        this%field_max = field_max
        this%field_min = field_min
        !
        IF (this%field_aware .AND. this%mode == 'fa-ionic') THEN
            ALLOCATE (this%ion_field(this%ions%number))
            ALLOCATE (this%dion_field_drho(this%ions%number))
            ALLOCATE (this%partial_of_ion_field(3, this%ions%number, this%ions%number))
            ALLOCATE (environ_function_erfc :: this%local_spheres(this%ions%number))
        END IF
        !
        !--------------------------------------------------------------------------------
        ! Soft spheres
        !
        IF (this%mode == 'ionic' .OR. this%mode == 'fa-ionic') &
            CALL this%set_soft_spheres()
        !
        !--------------------------------------------------------------------------------
        ! Densities
        !
        local_label = 'boundary_'//TRIM(ADJUSTL(label))
        !
        CALL this%scaled%init(cell, local_label)
        !
        IF (this%mode == 'electronic' .OR. this%mode == 'full' .OR. &
            this%mode == 'fa-electronic' .OR. this%mode == 'fa-full') THEN
            !
            local_label = 'boundary_density_'//TRIM(ADJUSTL(label))
            !
            CALL this%density%init(cell, local_label)
            !
            local_label = 'dboundary_'//TRIM(ADJUSTL(label))
            !
            CALL this%dscaled%init(cell, local_label)
            !
            local_label = 'd2boundary_'//TRIM(ADJUSTL(label))
            !
            CALL this%d2scaled%init(cell, local_label)
            !
        END IF
        !
        IF (this%deriv >= 1) THEN
            local_label = 'gradboundary_'//TRIM(ADJUSTL(label))
            !
            CALL this%gradient%init(cell, local_label)
            !
        END IF
        !
        IF (this%deriv >= 2) THEN
            local_label = 'laplboundary_'//TRIM(ADJUSTL(label))
            !
            CALL this%laplacian%init(cell, local_label)
            !
        END IF
        !
        IF (this%deriv >= 3) THEN
            local_label = 'dsurface_'//TRIM(ADJUSTL(label))
            !
            CALL this%dsurface%init(cell, local_label)
            !
        END IF
        !
        IF (this%solvent_aware) THEN
            local_label = 'local_'//TRIM(ADJUSTL(label))
            !
            CALL this%local%init(cell, local_label)
            !
            local_label = 'probe_'//TRIM(ADJUSTL(label))
            !
            CALL this%probe%init(cell, local_label)
            !
            local_label = 'filling_'//TRIM(ADJUSTL(label))
            !
            CALL this%filling%init(cell, local_label)
            !
            local_label = 'dfilling_'//TRIM(ADJUSTL(label))
            !
            CALL this%dfilling%init(cell, local_label)
            !
            IF (this%deriv >= 3) THEN
                local_label = 'hessboundary_'//TRIM(ADJUSTL(label))
                !
                CALL this%hessian%init(cell, local_label)
                !
            END IF
            !
        END IF
        !
        IF (this%field_aware) THEN
            !
            CALL io%error(sub_name, 'field-aware not yet implimented', 1)
            !
            IF (this%mode == 'fa-electronic' .OR. &
                this%mode == 'fa-full') THEN
                !
                local_label = 'normal_field_'//TRIM(ADJUSTL(label))
                !
                CALL this%normal_field%init(cell, local_label)
                !
            ELSE IF (this%mode == 'fa-ionic') THEN
                !
                DO i = 1, this%ions%number
                    CALL this%dion_field_drho(i)%init(cell)
                END DO
                !
            ELSE
                CALL io%error(sub_name, 'Boundary must be field-aware', 1)
            END IF
            !
        END IF
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_environ_boundary
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
        copy%cores => this%cores
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
        !
        copy%field_aware = this%field_aware
        copy%field_factor = this%field_factor
        copy%charge_asymmetry = this%charge_asymmetry
        copy%field_max = this%field_max
        copy%field_min = this%field_min
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
        IF (this%need_system) CALL this%simple%copy(copy%simple)
        !
        IF (this%solvent_aware) CALL this%solvent_probe%copy(copy%solvent_probe)
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
                    CALL io%error(sub_name, &
                                  'Wrong update status, possibly &
                                  &missing ionic update', 1)
                !
                this%density%of_r = this%electrons%density%of_r + this%ions%core%of_r
                !
                CALL this%boundary_of_density()
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
                CALL this%boundary_of_density()
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
                CALL this%update_soft_spheres()
                !
                CALL this%boundary_of_functions()
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
                CALL this%boundary_of_system()
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
            CALL io%error(sub_name, 'Unrecognized boundary mode', 1)
            !
        END SELECT
        !
        !--------------------------------------------------------------------------------
        ! Solvent-aware interface
        !
        IF (this%update_status == 2 .AND. this%solvent_aware) &
            CALL this%solvent_aware_boundary()
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
        CLASS(environ_boundary), INTENT(INOUT) :: this
        !
        INTEGER :: i
        !
        CHARACTER(LEN=80) :: sub_name = 'destroy_environ_boundary'
        !
        !--------------------------------------------------------------------------------
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
        IF (this%need_ions) THEN
            IF (this%mode == 'ionic' .OR. this%mode == 'fa-ionic') THEN
                !
                CALL destroy_environ_functions(this%soft_spheres, this%ions%number)
                !
                IF (this%field_aware .AND. this%mode == 'fa-ionic') THEN
                    !
                    CALL io%error(sub_name, 'field-aware not yet implimented ', 1)
                    !
                    DEALLOCATE (this%ion_field)
                    DEALLOCATE (this%partial_of_ion_field)
                    !
                    CALL destroy_environ_functions(this%local_spheres, this%ions%number)
                    !
                    DEALLOCATE (this%dion_field_drho)
                END IF
                !
            END IF
            !
            IF (.NOT. ASSOCIATED(this%ions)) CALL io%destroy_error(sub_name)
            !
            NULLIFY (this%ions)
        ELSE
            !
            IF (ASSOCIATED(this%ions)) &
                CALL io%error(sub_name, 'Found an unexpected associated object', 1)
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
            !
            CALL this%simple%destroy()
            !
            IF (ASSOCIATED(this%system)) NULLIFY (this%system)
            !
        END IF
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE destroy_environ_boundary
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
        de_dboundary%of_r = de_dboundary%of_r - confine * rhoelec%of_r * e2 / 2.D0
        ! the functional derivative of the confine term is - confine * rho^elec(r)
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
        de_dboundary%of_r = de_dboundary%of_r + pressure * e2 / 2.D0
        ! the functional derivative of the volume term is just unity
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
        de_dboundary%of_r = de_dboundary%of_r + &
                            surface_tension * this%dsurface%of_r * e2 / 2.D0
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
    SUBROUTINE boundary_of_density(this, density)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_density), TARGET, INTENT(IN), OPTIONAL :: density
        !
        CLASS(environ_boundary), TARGET, INTENT(INOUT) :: this
        !
        INTEGER, POINTER :: ir_end, stype, deriv
        REAL(DP), POINTER :: const, rhomax, rhomin, tbeta
        REAL(DP), DIMENSION(:), POINTER :: rho, eps, deps, d2eps, lapleps, dsurface
        REAL(DP), POINTER :: gradeps(:, :)
        !
        TYPE(environ_hessian), POINTER :: hessian
        TYPE(environ_density), POINTER :: local_density
        !
        INTEGER :: ir, ipol, jpol
        !
        CHARACTER(LEN=80) :: sub_name = 'boundary_of_density'
        !
        !--------------------------------------------------------------------------------
        !
        IF (PRESENT(density)) THEN
            local_density => density
        ELSE
            local_density => this%density
        END IF
        !
        IF (.NOT. ASSOCIATED(local_density%cell, this%scaled%cell)) &
            CALL io%error(sub_name, 'Inconsistent domains', 1)
        !
        ir_end => local_density%cell%ir_end
        rho => local_density%of_r
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
                CALL hessian%init(local_density%cell)
                !
            END IF
            !
        END IF
        !
        ASSOCIATE (derivatives => this%cores%derivatives)
            !
            SELECT CASE (this%derivatives_method)
                !
            CASE ('fft')
                !
                IF (deriv == 1 .OR. deriv == 2) &
                    CALL derivatives%gradient(this%scaled, this%gradient)
                !
                IF (deriv == 2) CALL derivatives%laplacian(this%scaled, this%laplacian)
                !
                IF (deriv == 3) &
                    CALL this%calc_dsurface(this%scaled, this%gradient, &
                                            this%laplacian, hessian, this%dsurface)
                !
            CASE ('chain')
                !
                IF (deriv == 1 .OR. deriv == 2) &
                    CALL derivatives%gradient(local_density, this%gradient)
                !
                IF (deriv == 2) &
                    CALL derivatives%laplacian(local_density, this%laplacian)
                !
                IF (deriv == 3) THEN
                    !
                    CALL this%calc_dsurface(local_density, this%gradient, &
                                            this%laplacian, hessian, this%dsurface)
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
                    DO ipol = 1, 3
                        gradeps(ipol, :) = gradeps(ipol, :) * deps(:)
                    END DO
                    !
                END IF
                !
            CASE DEFAULT
                CALL io%error(sub_name, "Unexpected derivatives method", 1)
                !
            END SELECT
            !
        END ASSOCIATE
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
        INTEGER, POINTER :: ir_end, deriv
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
        SELECT CASE (this%derivatives_method)
            !
        CASE ('fft')
            !
            IF (deriv == 1 .OR. deriv == 2) &
                CALL this%cores%derivatives%gradient(this%scaled, this%gradient)
            !
            IF (deriv == 2) &
                CALL this%cores%derivatives%laplacian(this%scaled, this%laplacian)
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
        CASE DEFAULT
            CALL io%error(sub_name, "Unexpected derivatives method", 1)
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
        INTEGER, POINTER :: ir_end, deriv
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
        SELECT CASE (this%derivatives_method)
            !
        CASE ('fft')
            !
            IF (deriv == 1 .OR. deriv == 2) &
                CALL this%cores%derivatives%gradient(this%scaled, this%gradient)
            !
            IF (deriv == 2) &
                CALL this%cores%derivatives%laplacian(this%scaled, this%laplacian)
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
                CALL calc_dsurface_no_pre(cell%nnr, ir_end, this%gradient%of_r, &
                                          hesslocal%of_r, this%dsurface%of_r)
                !
            END IF
            !
        CASE DEFAULT
            CALL io%error(sub_name, "Unexpected derivatives method", 1)
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
            CALL io%error(sub_name, 'Missing details of ions', 1)
        END IF
        !
        IF (index > number) &
            CALL io%error(sub_name, 'Index greater than number of ions', 1)
        !
        IF (index <= 0) &
            CALL io%error(sub_name, 'Index of ion is zero or lower', 1)
        !
        IF (this%mode == 'ionic' .AND. &
            .NOT. ALLOCATED(this%soft_spheres)) &
            CALL io%error(sub_name, 'Missing details of ionic boundary', 1)
        !
        IF (this%mode == 'full') THEN
            !
            IF (.NOT. ALLOCATED(this%ions%core_electrons)) &
                CALL io%error(sub_name, 'Missing details of core electrons', 1)
            !
            IF (.NOT. ASSOCIATED(this%dscaled%cell, cell)) &
                CALL io%error(sub_name, 'Mismatch or unassociated boundary derivative', 1)
            !
        END IF
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
            IF (spurious_force > tolspuriousforce .AND. io%lnode) &
                WRITE (io%unit, 1000) index, spurious_force
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
        !--------------------------------------------------------------------------------
        !
1000    FORMAT(' WARNING: Unphysical forces due to core electrons are non-negligible ', /, &
               ' atom type ', I3, ' is subject to a spurious force of ', F12.6)
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
        ASSOCIATE (cell => this%scaled%cell, &
                   derivatives => this%cores%derivatives, &
                   derivatives_method => this%derivatives_method, &
                   ir_end => this%scaled%cell%ir_end, &
                   deriv => this%deriv, &
                   thr => this%filling_threshold, &
                   spr => this%filling_spread)
            !
            !----------------------------------------------------------------------------
            !
            CALL filled_fraction%init(cell)
            !
            IF (deriv >= 2 .AND. derivatives_method /= 'fft') CALL d2filling%init(cell)
            !
            !----------------------------------------------------------------------------
            ! Step 0: save local interface function for later use
            !
            this%local%of_r = this%scaled%of_r
            !
            !----------------------------------------------------------------------------
            ! Step 1: compute the convolution function,
            !         this may be made moved out of here
            !
            CALL this%solvent_probe%density(this%probe, .TRUE.)
            !
            this%probe%of_r = this%probe%of_r / this%probe%integrate()
            !
            !----------------------------------------------------------------------------
            ! Step 2: compute filled fraction,
            !         i.e. convolution of local boundary with probe
            !
            CALL derivatives%convolution(this%local, this%probe, filled_fraction)
            !
            !----------------------------------------------------------------------------
            ! Step 3: compute the filling function and its derivative
            !
            this%filling%of_r = 0.D0
            this%dfilling%of_r = 0.D0
            !
            DO ir = 1, ir_end
                !
                this%filling%of_r(ir) = 1.D0 - &
                                        sfunct2(filled_fraction%of_r(ir), thr, spr)
                !
                this%dfilling%of_r(ir) = -dsfunct2(filled_fraction%of_r(ir), thr, spr)
                !
                IF (deriv >= 2 .AND. derivatives_method /= 'fft') &
                    d2filling%of_r(ir) = -d2sfunct2(filled_fraction%of_r(ir), thr, spr)
                !
            END DO
            !
            !----------------------------------------------------------------------------
            ! Step 4: compute solvent-aware interface
            !
            this%scaled%of_r = this%local%of_r + &
                               (1.D0 - this%local%of_r) * this%filling%of_r
            !
            !----------------------------------------------------------------------------
            ! Step 5: compute boundary derivatives, if needed
            !
            SELECT CASE (derivatives_method)
                !
            CASE ('fft')
                !
                IF (deriv == 1 .OR. deriv == 2) &
                    CALL derivatives%gradient(this%scaled, this%gradient)
                !
                IF (deriv == 2) CALL derivatives%laplacian(this%scaled, this%laplacian)

                IF (deriv == 3) &
                    CALL this%calc_dsurface(this%scaled, this%gradient, this%laplacian, &
                                            this%hessian, this%dsurface)
                !
            CASE ('chain', 'highmem', 'lowmem')
                !
                !------------------------------------------------------------------------
                ! Allocate local fields for derivatives of convolution
                !
                IF (deriv >= 1) CALL gradlocal%init(cell)
                !
                IF (deriv >= 2) CALL lapllocal%init(cell)
                !
                IF (deriv >= 3) CALL hesslocal%init(cell)
                !
                !------------------------------------------------------------------------
                ! Compute derivative of convolution with probe
                !
                IF (deriv > 1) &
                    CALL this%convolution_deriv(deriv, gradlocal, lapllocal, hesslocal)
                !
                !------------------------------------------------------------------------
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
                !------------------------------------------------------------------------
                ! Recompute dsurface, if needed
                !
                IF (deriv >= 3) THEN
                    !
                    CALL calc_dsurface_no_pre(this%scaled%cell%nnr, ir_end, &
                                              this%gradient%of_r, this%hessian%of_r, &
                                              this%dsurface%of_r)
                    !
                END IF
                !
            CASE DEFAULT
                CALL io%error(sub_name, "Unexpected derivatives method", 1)
                !
            END SELECT
            !
            !----------------------------------------------------------------------------
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
        END ASSOCIATE
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
        TYPE(environ_density) :: local
        !
        CHARACTER(LEN=80) :: sub_name = 'calc_solvent_aware_de_dboundary'
        !
        !--------------------------------------------------------------------------------
        !
        CALL local%init(this%scaled%cell)
        !
        !--------------------------------------------------------------------------------
        ! Step 1: compute (1-s)*de_dboudary*dfilling
        !
        local%of_r = (1.D0 - this%local%of_r) * de_dboundary%of_r * this%dfilling%of_r
        !
        !--------------------------------------------------------------------------------
        ! Step 2: compute convolution with the probe function
        !
        CALL this%cores%derivatives%convolution(this%probe, local, local)
        !
        !--------------------------------------------------------------------------------
        ! Step 3: update the functional derivative of the energy wrt boundary
        !
        de_dboundary%of_r = de_dboundary%of_r * (1.D0 - this%filling%of_r) + local%of_r
        !
        CALL local%destroy()
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
        ASSOCIATE (derivatives => this%cores%derivatives)
            !
            IF (deriv >= 1) THEN
                !
                CALL derivatives%convolution(this%probe, this%gradient, grad)
                !
                CALL grad%update_modulus()
                !
            END IF
            !
            IF (deriv >= 2) &
                CALL derivatives%convolution(this%probe, this%laplacian, lapl)
            !
            IF (deriv >= 3) CALL derivatives%convolution(this%probe, this%hessian, hess)
            !
        END ASSOCIATE
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
        CALL this%cores%derivatives%hessian(dens, grad, hess)
        !
        lapl%of_r(:) = hess%of_r(1, 1, :) + hess%of_r(2, 2, :) + hess%of_r(3, 3, :)
        !
        CALL calc_dsurface_no_pre(dens%cell%nnr, dens%cell%ir_end, grad%of_r, &
                                  hess%of_r, dsurface%of_r)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_dsurface
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
    SUBROUTINE set_soft_spheres(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_boundary), INTENT(INOUT) :: this
        !
        INTEGER :: i
        !
        INTEGER, DIMENSION(this%ions%number) :: axes, dims
        REAL(DP), DIMENSION(this%ions%number) :: spreads, volumes
        !
        REAL(DP), ALLOCATABLE :: radii(:)
        !
        TYPE(environ_function_erfc) :: fsrc
        !
        CHARACTER(LEN=20) :: local_item = 'solvationrad'
        !
        !--------------------------------------------------------------------------------
        !
        axes = 1
        dims = 0
        spreads = this%softness
        volumes = 1.D0
        !
        CALL this%ions%get_iontype_array(radii, local_item)
        !
        radii = radii * this%alpha
        !
        CALL init_environ_functions(this%soft_spheres, fsrc, this%ions%number, 4, &
                                    axes, dims, radii, spreads, volumes, this%ions%tau)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE set_soft_spheres
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE update_soft_spheres(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_boundary), INTENT(INOUT) :: this
        !
        INTEGER :: i
        !
        !--------------------------------------------------------------------------------
        !
        DO i = 1, this%ions%number
            this%soft_spheres(i)%pos = this%ions%tau(:, i)
        END DO
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE update_soft_spheres
    !------------------------------------------------------------------------------------
    !>
    !! Switching function 0: goes from 1 to 0 when passing through the
    !! threshold
    !!
    !! \f[
    !!    1 + \frac{1 - (x/x_t)^k}{1 + (x/x_t)^k}
    !! \f]
    !! where \f$x_t\f$ is the threshold
    !!
    !------------------------------------------------------------------------------------
    FUNCTION sfunct0(x, xthr, fact)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP) :: sfunct0
        REAL(DP) :: x, xthr, fact
        !
        REAL(DP) :: arg
        !
        !--------------------------------------------------------------------------------
        !
        arg = (ABS(x) / xthr)**fact
        sfunct0 = 0.5D0 * (1.D0 + (1.D0 - arg) / (1.D0 + arg))
        !
        !--------------------------------------------------------------------------------
    END FUNCTION sfunct0
    !------------------------------------------------------------------------------------
    !>
    !! Derivative of switching function 0
    !!
    !------------------------------------------------------------------------------------
    FUNCTION dsfunct0(x, xthr, fact)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP) :: dsfunct0
        REAL(DP) :: x, xthr, fact
        !
        REAL(DP) :: arg
        !
        !--------------------------------------------------------------------------------
        !
        arg = (ABS(x) / xthr)**fact
        dsfunct0 = -fact * ABS(x)**(fact - 1.D0) / xthr**fact / (1.D0 + arg)**2
        !
        !--------------------------------------------------------------------------------
    END FUNCTION dsfunct0
    !------------------------------------------------------------------------------------
    !>
    !! Switching function 1 that goes from 1 to 0 when passing from
    !! xmin to xmax.
    !!
    !! NOTE: fact should be equal to LOG(xmax/xmin) but is
    !! passed in input to save time
    !!
    !! \f[
    !!    x - \sin(x)
    !! \f]
    !!
    !------------------------------------------------------------------------------------
    FUNCTION sfunct1(x, xmax, xmin, fact)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP) :: sfunct1
        REAL(DP) :: x, xmax, xmin, fact
        !
        REAL(DP) :: arg
        !
        !--------------------------------------------------------------------------------
        !
        IF (x <= xmin) THEN
            sfunct1 = 1.D0
        ELSE IF (x < xmax) THEN
            arg = tpi * LOG(xmax / ABS(x)) / fact
            sfunct1 = (arg - SIN(arg)) / tpi
        ELSE
            sfunct1 = 0.D0
        END IF
        !
        !--------------------------------------------------------------------------------
    END FUNCTION sfunct1
    !------------------------------------------------------------------------------------
    !>
    !! @brief Derivative of switching function 1
    !!
    !! NOTE: fact should be equal to LOG(xmax/xmin) but is passed in
    !! input to save time.
    !!
    !------------------------------------------------------------------------------------
    FUNCTION dsfunct1(x, xmax, xmin, fact)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP) :: dsfunct1
        REAL(DP) :: x, xmax, xmin, fact
        !
        REAL(DP) :: arg
        !
        !--------------------------------------------------------------------------------
        !
        IF (x <= xmin) THEN
            dsfunct1 = 0.D0
        ELSE IF (x < xmax) THEN
            arg = tpi * LOG(xmax / ABS(x)) / fact
            dsfunct1 = (COS(arg) - 1.D0) / ABS(x) / fact ! #TODO in fact should not use ABS(x)
        ELSE
            dsfunct1 = 0.D0
        END IF
        !
        !--------------------------------------------------------------------------------
    END FUNCTION dsfunct1
    !------------------------------------------------------------------------------------
    !>
    !! @brief Second derivative of switching function 1
    !!
    !! Note: fact should be equal to LOG(xmax/xmin) but is passed in
    !! input to save time
    !!
    !------------------------------------------------------------------------------------
    FUNCTION d2sfunct1(x, xmax, xmin, fact)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP) :: d2sfunct1
        REAL(DP) :: x, xmax, xmin, fact
        !
        REAL(DP) :: arg
        !
        !--------------------------------------------------------------------------------
        !
        IF (x <= xmin) THEN
            d2sfunct1 = 0.D0
        ELSE IF (x < xmax) THEN
            arg = tpi * LOG(xmax / ABS(x)) / fact
            d2sfunct1 = (tpi * SIN(arg) + fact * (1.D0 - COS(arg))) / (x * fact)**2
        ELSE
            d2sfunct1 = 0.D0
        END IF
        !
        !--------------------------------------------------------------------------------
    END FUNCTION d2sfunct1
    !------------------------------------------------------------------------------------
    !>
    !! Switching function 2, erfc() that goes from 1 to 0 when passing
    !! through xthr.
    !!
    !------------------------------------------------------------------------------------
    FUNCTION sfunct2(x, xthr, spread)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP) :: sfunct2
        REAL(DP) :: x, xthr, spread
        !
        REAL(DP) :: arg
        !
        !--------------------------------------------------------------------------------
        !
        arg = (x - xthr) / spread
        sfunct2 = 0.5D0 * environ_erfc(arg)
        !
        !--------------------------------------------------------------------------------
    END FUNCTION sfunct2
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    FUNCTION dsfunct2(x, xthr, spread)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP) :: dsfunct2
        REAL(DP) :: x, xthr, spread
        !
        REAL(DP) :: arg
        !
        !--------------------------------------------------------------------------------
        !
        arg = (x - xthr) / spread
        !
        IF (ABS(arg) > 6.D0) THEN ! 6.D0 is the threshold of environ_erfc(x)
            dsfunct2 = 0.D0
        ELSE
            dsfunct2 = -EXP(-arg**2) / sqrtpi / spread
        END IF
        !
        !--------------------------------------------------------------------------------
    END FUNCTION dsfunct2
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    FUNCTION d2sfunct2(x, xthr, spread)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP) :: d2sfunct2
        REAL(DP) :: x, xthr, spread
        !
        REAL(DP) :: arg
        !
        !--------------------------------------------------------------------------------
        !
        arg = (x - xthr) / spread
        IF (ABS(arg) > 6.D0) THEN
            d2sfunct2 = 0.D0
        ELSE
            d2sfunct2 = EXP(-arg**2) / sqrtpi / spread**2 * 2.D0 * arg
        END IF
        !
        !--------------------------------------------------------------------------------
    END FUNCTION d2sfunct2
    !------------------------------------------------------------------------------------
    !>
    !! Calculates the density-dependent dielectric constant
    !!
    !! ifunct = 0 => original Fattebert and Gygi function
    !!
    !------------------------------------------------------------------------------------
    FUNCTION boundfunct(rho, rhomax, rhomin, tbeta, const, ifunct)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP) :: boundfunct
        REAL(DP) :: rho
        REAL(DP) :: rhomax
        REAL(DP) :: rhomin
        REAL(DP) :: tbeta
        REAL(DP) :: const
        !
        INTEGER :: ifunct
        !
        REAL(DP) :: arg
        !
        CHARACTER(LEN=80) :: fun_name = 'boundfunct'
        !
        !--------------------------------------------------------------------------------
        !
        SELECT CASE (ifunct)
            !
        CASE (0)
            boundfunct = 1.D0 - sfunct0(rho, rhomax, tbeta)
            !
        CASE (1)
            boundfunct = 1.D0 - sfunct1(rho, rhomax, rhomin, tbeta)
            !
        CASE (2)
            !
            boundfunct = &
                (const - EXP(LOG(const) * sfunct1(rho, rhomax, rhomin, tbeta))) / &
                (const - 1.D0)
            !
        CASE DEFAULT
            CALL io%error(fun_name, 'Unknown boundary type', 1)
            !
        END SELECT
        !
        !--------------------------------------------------------------------------------
    END FUNCTION boundfunct
    !------------------------------------------------------------------------------------
    !>
    !! @brief Calculates the derivative of the density-dependent dielectric
    !! constant
    !!
    !! ifunct = 0 => original Fattebert and Gygi function
    !!
    !! @param[in]    rho      electrostatic density
    !! @param[in]    rhomax   maximum density cutoff
    !! @param[in]    rhomin   minimum density cutoff
    !! @param[in]    tbeta
    !! @param[in]    const
    !! @param[in]    ifunct
    !! @return       the second derivative of the boundary function
    !!
    !------------------------------------------------------------------------------------
    FUNCTION dboundfunct(rho, rhomax, rhomin, tbeta, const, ifunct)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP) :: dboundfunct
        REAL(DP) :: rho
        REAL(DP) :: rhomax
        REAL(DP) :: rhomin
        REAL(DP) :: tbeta
        REAL(DP) :: const
        !
        INTEGER :: ifunct
        !
        REAL(DP) :: arg
        !
        CHARACTER(LEN=80) :: fun_name = 'dboundfunct'
        !
        !--------------------------------------------------------------------------------
        !
        SELECT CASE (ifunct)
            !
        CASE (0)
            dboundfunct = -dsfunct0(rho, rhomax, tbeta)
            !
        CASE (1)
            dboundfunct = -dsfunct1(rho, rhomax, rhomin, tbeta)
            !
        CASE (2)
            !
            dboundfunct = -EXP(LOG(const) * sfunct1(rho, rhomax, rhomin, tbeta)) / &
                          (const - 1.D0) * LOG(const) * &
                          dsfunct1(rho, rhomax, rhomin, tbeta)
            !
        CASE DEFAULT
            CALL io%error(fun_name, 'Unknown boundary type', 1)
            !
        END SELECT
        !
        !--------------------------------------------------------------------------------
    END FUNCTION dboundfunct
    !------------------------------------------------------------------------------------
    !>
    !! @brief Calculates the second derivative of the density-dependent
    !! dielectric constant
    !!
    !! ifunct = 0 => original Fattebery and Gygi function
    !!
    !! @param[in]    rho      electrostatic density
    !! @param[in]    rhomax   maximum density cutoff
    !! @param[in]    rhomin   minimum density cutoff
    !! @param[in]    tbeta
    !! @param[in]    const
    !! @param[in]    ifunct
    !! @return       the second derivative of the boundary function
    !!
    !------------------------------------------------------------------------------------
    FUNCTION d2boundfunct(rho, rhomax, rhomin, tbeta, const, ifunct)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP) :: d2boundfunct
        REAL(DP) :: rho
        REAL(DP) :: rhomax
        REAL(DP) :: rhomin
        REAL(DP) :: tbeta
        REAL(DP) :: const
        !
        INTEGER :: ifunct
        !
        REAL(DP) :: arg, arg2
        !
        CHARACTER(LEN=80) :: fun_name = 'd2boundfunct'
        !
        !--------------------------------------------------------------------------------
        !
        SELECT CASE (ifunct)
            !
        CASE (0)
            CALL io%error(fun_name, 'Option not yet implemented', 1)
            !
        CASE (1)
            d2boundfunct = -d2sfunct1(rho, rhomax, rhomin, tbeta)
            !
        CASE (2)
            !
            d2boundfunct = -EXP(LOG(const) * sfunct1(rho, rhomax, rhomin, tbeta)) / &
                           (const - 1.D0) * LOG(const) * &
                           (LOG(const) * dsfunct1(rho, rhomax, rhomin, tbeta)**2 + &
                            d2sfunct1(rho, rhomax, rhomin, tbeta))
            !
        CASE DEFAULT
            CALL io%error(fun_name, 'Unknown boundary type', 1)
            !
        END SELECT
        !
        !--------------------------------------------------------------------------------
    END FUNCTION d2boundfunct
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_partial_of_boundary(n, i, local, gradlocal, partial)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: n, i
        TYPE(environ_density), INTENT(IN) :: local(n)
        TYPE(environ_gradient), INTENT(IN) :: gradlocal(n)
        !
        TYPE(environ_gradient), INTENT(INOUT) :: partial
        !
        INTEGER :: j, ipol
        !
        CHARACTER(LEN=80) :: sub_name = 'calc_partial_of_boundary'
        !
        !--------------------------------------------------------------------------------
        !
        IF (i > n) CALL io%error(sub_name, 'Index out of bound', 1)
        !
        DO ipol = 1, 3
            partial%of_r(ipol, :) = gradlocal(i)%of_r(ipol, :)
            !
            DO j = 1, n
                !
                IF (j == i) CYCLE
                !
                partial%of_r(ipol, :) = partial%of_r(ipol, :) * local(j)%of_r(:)
            END DO
            !
        END DO
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_partial_of_boundary
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_gradient_of_boundary_highmem(n, local, gradlocal, gradient)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: n
        TYPE(environ_density), INTENT(IN) :: local(n)
        TYPE(environ_gradient), INTENT(IN) :: gradlocal(n)
        !
        TYPE(environ_gradient), INTENT(INOUT) :: gradient
        !
        INTEGER :: i, j, ipol
        TYPE(environ_cell), POINTER :: cell
        TYPE(environ_gradient) :: partial
        !
        !--------------------------------------------------------------------------------
        !
        cell => gradient%cell
        !
        CALL partial%init(cell)
        !
        gradient%of_r = 0.D0
        !
        DO i = 1, n
            !
            CALL calc_partial_of_boundary(n, i, local, gradlocal, partial)
            !
            gradient%of_r = gradient%of_r + partial%of_r
        END DO
        !
        CALL partial%destroy()
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_gradient_of_boundary_highmem
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_laplacian_of_boundary_highmem(n, local, gradlocal, lapllocal, &
                                                  laplacian)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: n
        TYPE(environ_density), INTENT(IN) :: local(n)
        TYPE(environ_gradient), INTENT(IN) :: gradlocal(n)
        TYPE(environ_density), INTENT(IN) :: lapllocal(n)
        !
        TYPE(environ_density), INTENT(INOUT) :: laplacian
        !
        INTEGER :: i, j, k, ipol
        TYPE(environ_cell), POINTER :: cell
        TYPE(environ_density) :: tmp
        !
        !--------------------------------------------------------------------------------
        !
        cell => laplacian%cell
        !
        CALL tmp%init(cell)
        !
        laplacian%of_r = 0.D0
        !
        DO i = 1, n
            !
            DO j = 1, n
                !
                IF (j == i) THEN
                    tmp%of_r = lapllocal(i)%of_r
                ELSE
                    CALL gradlocal(i)%scalar_product(gradlocal(j), tmp)
                END IF
                !
                DO k = 1, n
                    !
                    IF (k == j .OR. k == i) CYCLE
                    !
                    tmp%of_r = tmp%of_r * local(k)%of_r
                END DO
                !
                laplacian%of_r = laplacian%of_r + tmp%of_r
            END DO
            !
        END DO
        !
        CALL tmp%destroy()
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_laplacian_of_boundary_highmem
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_dsurface_of_boundary_highmem(n, local, gradlocal, hesslocal, &
                                                 gradient, laplacian, hessian, dsurface)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: n
        TYPE(environ_density), INTENT(IN) :: local(n)
        TYPE(environ_gradient), INTENT(IN) :: gradlocal(n)
        TYPE(environ_hessian), INTENT(IN) :: hesslocal(n)
        !
        TYPE(environ_gradient), INTENT(INOUT) :: gradient
        TYPE(environ_density), INTENT(INOUT) :: laplacian, dsurface
        TYPE(environ_hessian), INTENT(INOUT) :: hessian
        !
        INTEGER :: i, j, k, ipol, jpol
        TYPE(environ_cell), POINTER :: cell
        TYPE(environ_density) :: dens
        TYPE(environ_gradient) :: partial
        !
        !--------------------------------------------------------------------------------
        !
        cell => laplacian%cell
        !
        CALL dens%init(cell)
        !
        CALL partial%init(cell)
        !
        gradient%of_r = 0.D0
        !
        DO i = 1, n
            !
            CALL calc_partial_of_boundary(n, i, local, gradlocal, partial)
            !
            gradient%of_r = gradient%of_r + partial%of_r
            !
            DO j = 1, n
                !
                DO ipol = 1, 3
                    !
                    DO jpol = 1, 3
                        !
                        IF (j == i) THEN
                            dens%of_r(:) = hesslocal(i)%of_r(ipol, jpol, :)
                        ELSE
                            !
                            dens%of_r(:) = gradlocal(i)%of_r(ipol, :) * &
                                           gradlocal(j)%of_r(jpol, :)
                            !
                        END IF
                        !
                        DO k = 1, n
                            !
                            IF (k == j .OR. k == i) CYCLE
                            !
                            dens%of_r = dens%of_r * local(k)%of_r
                        END DO
                        !
                        hessian%of_r(ipol, jpol, :) = hessian%of_r(ipol, jpol, :) + &
                                                      dens%of_r(:)
                        !
                    END DO
                    !
                END DO
                !
            END DO
            !
        END DO
        !
        !--------------------------------------------------------------------------------
        ! Final operations
        !
        laplacian%of_r = hessian%of_r(1, 1, :) + hessian%of_r(2, 2, :) + &
                         hessian%of_r(3, 3, :)
        !
        CALL calc_dsurface_no_pre(cell%nnr, cell%ir_end, gradient%of_r, hessian%of_r, &
                                  dsurface%of_r)
        !
        CALL dens%destroy()
        !
        CALL partial%destroy()
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_dsurface_of_boundary_highmem
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_gradient_of_boundary_lowmem(n, local, gradlocal, scaled, gradient)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: n
        TYPE(environ_density), INTENT(IN) :: scaled ! soft sphere interface function
        TYPE(environ_density), INTENT(IN) :: local(n)
        TYPE(environ_gradient), INTENT(IN) :: gradlocal(n)
        !
        TYPE(environ_gradient), INTENT(INOUT) :: gradient
        !
        INTEGER :: i, j, ipol
        TYPE(environ_cell), POINTER :: cell
        !
        !--------------------------------------------------------------------------------
        !
        cell => gradient%cell
        !
        gradient%of_r = 0.D0
        !
        !--------------------------------------------------------------------------------
        ! Temporary quotient
        !
        DO i = 1, n
            !
            DO j = 1, cell%nnr
                !
                IF (ABS(local(i)%of_r(j)) <= bound_tol) CYCLE
                !
                DO ipol = 1, 3
                    gradient%of_r(ipol, j) = gradient%of_r(ipol, j) + &
                                             (gradlocal(i)%of_r(ipol, j) / &
                                              local(i)%of_r(j) * scaled%of_r(j))
                    !
                END DO
                !
            END DO
            !
        END DO
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_gradient_of_boundary_lowmem
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_laplacian_of_boundary_lowmem(n, local, gradlocal, lapllocal, &
                                                 scaled, gradient, laplacian)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: n
        TYPE(environ_density), INTENT(IN) :: scaled ! soft sphere interface function
        TYPE(environ_density), INTENT(IN) :: local(n)
        TYPE(environ_gradient), INTENT(IN) :: gradlocal(n)
        TYPE(environ_density), INTENT(IN) :: lapllocal(n)
        TYPE(environ_gradient), INTENT(IN) :: gradient
        !
        TYPE(environ_density), INTENT(INOUT) :: laplacian
        !
        INTEGER :: i, j, k, ipol
        TYPE(environ_cell), POINTER :: cell
        !
        !--------------------------------------------------------------------------------
        !
        cell => laplacian%cell
        !
        DO i = 1, n
            !
            DO j = 1, cell%nnr
                !
                IF (ABS(local(i)%of_r(j)) <= bound_tol) CYCLE
                !
                laplacian%of_r(j) = laplacian%of_r(j) + &
                                    (lapllocal(i)%of_r(j) / &
                                     local(i)%of_r(j) * scaled%of_r(j))
                !
                DO ipol = 1, 3
                    !
                    laplacian%of_r(j) = laplacian%of_r(j) - &
                                        ((gradlocal(i)%of_r(ipol, j)**2 / &
                                          local(i)%of_r(j)**2) * scaled%of_r(j))
                    !
                    laplacian%of_r(j) = laplacian%of_r(j) + &
                                        (gradient%of_r(ipol, j) * &
                                         gradlocal(i)%of_r(ipol, j) / local(i)%of_r(j))
                    !
                END DO
                !
            END DO
            !
        END DO
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_laplacian_of_boundary_lowmem
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_dsurface_of_boundary_lowmem(n, local, gradlocal, hesslocal, &
                                                gradient, laplacian, hessian, &
                                                scaled, dsurface)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: n
        TYPE(environ_density), INTENT(IN) :: scaled
        TYPE(environ_density), INTENT(IN) :: local(n)
        TYPE(environ_gradient), INTENT(IN) :: gradlocal(n)
        TYPE(environ_hessian), INTENT(IN) :: hesslocal(n)
        TYPE(environ_gradient), INTENT(IN) :: gradient
        !
        TYPE(environ_density), INTENT(INOUT) :: laplacian
        TYPE(environ_density), INTENT(INOUT) :: dsurface
        TYPE(environ_hessian), INTENT(INOUT) :: hessian
        !
        INTEGER :: i, j, k, ipol, jpol
        TYPE(environ_cell), POINTER :: cell
        !
        !--------------------------------------------------------------------------------
        !
        cell => laplacian%cell
        !
        DO i = 1, n
            !
            DO j = 1, cell%nnr
                !
                IF (ABS(local(i)%of_r(j)) <= bound_tol) CYCLE
                !
                DO ipol = 1, 3
                    !
                    DO jpol = 1, 3
                        !
                        hessian%of_r(ipol, jpol, j) = &
                            hessian%of_r(ipol, jpol, j) + &
                            (hesslocal(i)%of_r(ipol, jpol, j) / &
                             local(i)%of_r(j) * scaled%of_r(j))
                        !
                        hessian%of_r(ipol, jpol, j) = &
                            hessian%of_r(ipol, jpol, j) - &
                            ((gradlocal(i)%of_r(ipol, j) * gradlocal(i)%of_r(jpol, j) / &
                              local(i)%of_r(j)**2) * scaled%of_r(j))
                        !
                        hessian%of_r(ipol, jpol, j) = &
                            hessian%of_r(ipol, jpol, j) + &
                            (gradient%of_r(ipol, j) * gradlocal(i)%of_r(jpol, j) / &
                             local(i)%of_r(j))
                        !
                    END DO
                    !
                END DO
                !
            END DO
            !
        END DO
        !
        laplacian%of_r = hessian%of_r(1, 1, :) + hessian%of_r(2, 2, :) + &
                         hessian%of_r(3, 3, :)
        !
        CALL calc_dsurface_no_pre(cell%nnr, cell%ir_end, gradient%of_r, hessian%of_r, &
                                  dsurface%of_r)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_dsurface_of_boundary_lowmem
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_dsurface_no_pre(n, iend, grad, hess, dsurface)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: n, iend
        REAL(DP), INTENT(IN) :: grad(3, n)
        REAL(DP), INTENT(IN) :: hess(3, 3, n)
        !
        REAL(DP), INTENT(OUT) :: dsurface(n)
        !
        REAL(DP), PARAMETER :: toldsurface = 1.D-50
        !
        INTEGER :: ipol, jpol, i
        REAL(DP) :: gmod
        !
        !--------------------------------------------------------------------------------
        !
        DO i = 1, iend
            dsurface(i) = 0.D0
            gmod = SUM(grad(:, i)**2)
            !
            IF (gmod < toldsurface) CYCLE
            !
            DO ipol = 1, 3
                !
                DO jpol = 1, 3
                    !
                    IF (ipol == jpol) CYCLE
                    !
                    dsurface(i) = dsurface(i) + &
                                  grad(ipol, i) * grad(jpol, i) * hess(ipol, jpol, i) - &
                                  grad(ipol, i) * grad(ipol, i) * hess(jpol, jpol, i)
                    !
                END DO
                !
            END DO
            !
            dsurface(i) = dsurface(i) / gmod / SQRT(gmod)
        END DO
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_dsurface_no_pre
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
        CLASS(environ_boundary), INTENT(IN) :: this
        INTEGER, INTENT(IN), OPTIONAL :: verbose, debug_verbose, unit
        !
        INTEGER :: base_verbose, local_verbose, passed_verbose, local_unit
        !
        CHARACTER(LEN=80) :: sub_name = 'print_environ_boundary'
        !
        !--------------------------------------------------------------------------------
        !
        IF (PRESENT(debug_verbose)) THEN
            base_verbose = debug_verbose
            !
            IF (PRESENT(verbose)) THEN
                local_verbose = verbose
            ELSE
                local_verbose = debug_verbose
            END IF
            !
            passed_verbose = verbose - 1
            !
        ELSE IF (io%verbosity > 0) THEN
            base_verbose = io%verbosity
            !
            IF (PRESENT(verbose)) THEN
                local_verbose = base_verbose + verbose
            ELSE
                local_verbose = base_verbose
            END IF
            !
            passed_verbose = local_verbose - base_verbose - 1
            !
        ELSE
            RETURN
        END IF
        !
        IF (PRESENT(unit)) THEN
            local_unit = unit
        ELSE
            local_unit = io%debug_unit
        END IF
        !
        IF (local_verbose >= 1) THEN
            !
            IF (io%lnode) THEN
                WRITE (local_unit, 1100)
                WRITE (local_unit, 1101) this%label, this%mode
            END IF
            !
            IF (this%need_electrons) THEN
                !
                IF (io%lnode) THEN
                    WRITE (local_unit, 1102) this%b_type
                    !
                    SELECT CASE (this%b_type)
                        !
                    CASE (0)
                        WRITE (local_unit, 1103) this%rhozero, this%tbeta
                        !
                    CASE (1)
                        WRITE (local_unit, 1104) this%rhomax, this%rhomin
                        !
                        IF (local_verbose >= 3) WRITE (local_unit, 1105) this%fact
                        !
                    CASE (2)
                        WRITE (local_unit, 1106) this%rhomax, this%rhomin
                        !
                    CASE DEFAULT
                        CALL io%error(sub_name, 'Unexpected boundary type', 1)
                        !
                    END SELECT
                    !
                END IF
                !
                IF (local_verbose >= 4) THEN
                    !
                    CALL this%density%printout(passed_verbose, debug_verbose, local_unit)
                    !
                    IF (io%lnode .AND. this%need_ions) WRITE (local_unit, 1107)
                    !
                END IF
                !
                IF (local_verbose >= 5) THEN
                    !
                    CALL this%dscaled%printout(passed_verbose, debug_verbose, local_unit)
                    !
                    CALL this%d2scaled%printout(passed_verbose, debug_verbose, &
                                                local_unit)
                    !
                END IF
                !
            ELSE IF (this%need_ions) THEN
                !
                IF (io%lnode) WRITE (local_unit, 1108) this%alpha, this%softness
                !
                IF (local_verbose >= 3) &
                    CALL print_environ_functions(this%soft_spheres, this%ions%number, &
                                                 passed_verbose, debug_verbose, &
                                                 local_unit)
                !
            ELSE IF (io%lnode .AND. this%need_system) THEN
                !
                WRITE (local_unit, 1109) &
                    this%simple%pos, this%simple%width, &
                    this%simple%spread, this%simple%dim, &
                    this%simple%axis
                !
            END IF
            !
            IF (io%lnode) THEN
                WRITE (local_unit, 1110) this%volume
                !
                IF (this%deriv >= 1) WRITE (local_unit, 1111) this%surface
                !
            END IF
            !
            IF (local_verbose >= 4) &
                CALL this%scaled%printout(passed_verbose, debug_verbose, local_unit)
            !
            IF (this%solvent_aware) THEN
                !
                IF (io%lnode) &
                    WRITE (local_unit, 1112) &
                    this%filling_threshold, this%filling_spread, &
                    this%solvent_probe%width, this%solvent_probe%spread
                !
                IF (local_verbose >= 4) THEN
                    !
                    CALL this%local%printout(passed_verbose, debug_verbose, local_unit)
                    !
                    CALL this%filling%printout(passed_verbose, debug_verbose, local_unit)
                    !
                END IF
                !
                IF (local_verbose >= 5) THEN
                    !
                    CALL this%dfilling%printout(passed_verbose, debug_verbose, &
                                                local_unit)
                    !
                    CALL this%probe%printout(passed_verbose, debug_verbose, local_unit)
                    !
                END IF
                !
            END IF
            !
            IF (local_verbose >= 5) THEN
                !
                IF (this%deriv >= 1) &
                    CALL this%gradient%printout(passed_verbose, debug_verbose, &
                                                local_unit)
                !
                IF (this%deriv >= 2) &
                    CALL this%laplacian%printout(passed_verbose, debug_verbose, &
                                                 local_unit)
                !
                IF (this%deriv == 3) &
                    CALL this%dsurface%printout(passed_verbose, debug_verbose, &
                                                local_unit)
                !
            END IF
            !
        END IF
        !
        FLUSH (local_unit)
        !
        !--------------------------------------------------------------------------------
        !
1100    FORMAT(/, 4('%'), ' BOUNDARY ', 66('%'))
        !
1101    FORMAT(/, ' boundary label             = ', A20, /, &
                ' boundary mode              = ', A20)
        !
1102    FORMAT(/, ' boundary is built as a type-', I1, ' function of a smooth density')
        !
1103    FORMAT(/, ' using the Fattebert-Gygi function:', /, &
                ' rhozero                    = ', F14.7, /, &
                ' 2*beta                     = ', F14.7)
        !
1104    FORMAT(/, ' using the optimal SCCS function:', /, &
                ' rhomax                     = ', F14.7, /, &
                ' rhomin                     = ', F14.7)
        !
1105    FORMAT(' log(rhomax/rhomin)         = ', F14.7)
        !
1106    FORMAT(/, ' using the modified SCCS function:', /, &
                ' rhomax                     = ', F14.7, /, &
                ' rhomin                     = ', F14.7)
        !
1107    FORMAT(/, ' adding fictitious core-electrons')
        !
1108    FORMAT(/, ' boundary is built from soft-spheres centered on ionic positions:', /, &
                ' solvent-dependent scaling  = ', F14.7, /, &
                ' softness parameter         = ', F14.7)
        !
1109    FORMAT(/, ' boundary is built as an analytic function centered on system position:', /, &
                ' center of the boundary     = ', 3F14.7, /, &
                ' distance from the center   = ', F14.7, /, &
                ' spread of the interface    = ', F14.7, /, &
                ' dimensionality             = ', I14, /, &
                ' axis                       = ', I14)
        !
1110    FORMAT(/, ' volume of the QM region    = ', F14.7)
        !
1111    FORMAT(/, ' surface of the QM region   = ', F14.7)
        !
1112    FORMAT(/, ' using solvent-aware boundary:', /, &
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
