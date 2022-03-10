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
    USE environ_param, ONLY: DP, e2, tpi
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
    USE boundary_tools
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
        TYPE(environ_density), ALLOCATABLE :: denloc(:)
        TYPE(environ_gradient), ALLOCATABLE :: gradloc(:)
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
        TYPE(environ_functions) :: soft_spheres
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
        REAL(DP) :: field_factor, field_asymmetry, field_max, field_min
        !
        TYPE(environ_functions) :: unscaled_spheres
        !
        REAL(DP), ALLOCATABLE :: ion_field(:)
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
        PROCEDURE :: fa_de_drho => calc_field_aware_de_drho
        PROCEDURE :: fa_dboundary_dions => calc_field_aware_dboundary_dions
        PROCEDURE :: ion_field_partial => calc_ion_field_partial
        !
        PROCEDURE :: boundary_of_density
        PROCEDURE :: boundary_of_functions
        PROCEDURE :: boundary_of_system
        !
        PROCEDURE :: convolution => compute_convolution_deriv
        PROCEDURE :: solvent_aware_boundary
        PROCEDURE :: calc_dsurface ! #TODO do we need this?
        PROCEDURE :: invert => invert_boundary
        !
        PROCEDURE, PRIVATE :: set_soft_spheres
        PROCEDURE, PRIVATE :: update_soft_spheres
        PROCEDURE, PRIVATE :: calc_ion_field
        PROCEDURE, PRIVATE :: calc_dion_field_drho
        PROCEDURE, PRIVATE :: scaling_of_field
        PROCEDURE, PRIVATE :: dscaling_of_field
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
        IF (ALLOCATED(this%ion_field)) CALL io%create_error(sub_name)
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
                                     filling_threshold, filling_spread, field_aware, &
                                     field_factor, field_asymmetry, field_max, &
                                     field_min, electrons, ions, system, cores, &
                                     deriv_method, cell, label)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: stype
        CHARACTER(LEN=*), INTENT(IN) :: mode, deriv_method
        LOGICAL, INTENT(IN) :: need_gradient, need_laplacian, need_hessian, field_aware
        !
        REAL(DP), INTENT(IN) :: rhomax, rhomin, tbeta, const, alpha, softness, &
                                system_distance, system_spread, solvent_radius, &
                                radial_scale, radial_spread, filling_threshold, &
                                filling_spread, field_factor, field_asymmetry, &
                                field_max, field_min
        !
        TYPE(environ_electrons), TARGET, INTENT(IN) :: electrons
        TYPE(environ_ions), TARGET, INTENT(IN) :: ions
        TYPE(environ_system), TARGET, INTENT(IN) :: system
        TYPE(environ_cell), INTENT(IN) :: cell
        TYPE(core_container), TARGET, INTENT(IN) :: cores
        CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: label
        !
        CLASS(environ_boundary), INTENT(INOUT) :: this
        !
        INTEGER :: i
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
        this%need_electrons = mode == 'electronic' .OR. mode == 'full' .OR. field_aware
        !
        IF (this%need_electrons) this%electrons => electrons
        !
        this%need_ions = mode == 'ionic' .OR. mode == 'full' .OR. field_aware
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
            CALL io%error(sub_name, "stype=2 boundary requires dielectric constant > 1", 1)
        !
        this%const = const
        this%alpha = alpha
        this%softness = softness
        !
        !--------------------------------------------------------------------------------
        !
        IF (this%need_system) &
            CALL this%simple%init(3, system%axis, system%dim, system_distance, &
                                  system_spread, 1.D0, system%com)
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
        this%field_aware = field_aware
        this%field_factor = field_factor
        this%field_asymmetry = field_asymmetry
        this%field_max = field_max
        this%field_min = field_min
        !
        IF (this%field_aware .AND. this%mode == 'ionic') THEN
            !
            ASSOCIATE (n => this%ions%number)
                !
                ALLOCATE (this%ion_field(n))
                ALLOCATE (this%dion_field_drho(n))
                ALLOCATE (this%partial_of_ion_field(3, n, n))
                ALLOCATE (environ_function_erfc :: this%unscaled_spheres%array(n))
                !
            END ASSOCIATE
            !
        END IF
        !
        !--------------------------------------------------------------------------------
        ! Soft spheres
        !
        IF (this%mode == 'ionic') THEN
            !
            CALL this%set_soft_spheres()
            !
            ALLOCATE (this%denloc(this%ions%number))
            ALLOCATE (this%gradloc(this%ions%number))
            !
            DO i=1,this%ions%number
                CALL this%denloc(i)%init(cell)
                CALL this%gradloc(i)%init(cell)
            END DO
            !
        END IF
        !
        !--------------------------------------------------------------------------------
        ! Densities
        !
        CALL this%scaled%init(cell, 'boundary_'//label)
        !
        IF (this%mode == 'electronic' .OR. this%mode == 'full') THEN
            !
            CALL this%density%init(cell, 'boundary_density_'//label)
            !
            CALL this%dscaled%init(cell, 'dboundary_'//label)
            !
            CALL this%d2scaled%init(cell, 'd2boundary_'//label)
            !
        END IF
        !
        IF (this%deriv >= 1) CALL this%gradient%init(cell, 'gradboundary_'//label)
        !
        IF (this%deriv >= 2) CALL this%laplacian%init(cell, 'laplboundary_'//label)
        !
        IF (this%deriv >= 3) CALL this%dsurface%init(cell, 'dsurface_'//label)
        !
        IF (this%solvent_aware) THEN
            !
            CALL this%local%init(cell, 'local_'//label)
            !
            CALL this%probe%init(cell, 'probe_'//label)
            !
            CALL this%filling%init(cell, 'filling_'//label)
            !
            CALL this%dfilling%init(cell, 'dfilling_'//label)
            !
            IF (this%deriv >= 3) CALL this%hessian%init(cell, 'hessboundary_'//label)
            !
        END IF
        !
        IF (this%field_aware) THEN
            !
            IF (this%mode == 'ionic') THEN
                !
                DO i = 1, this%ions%number
                    CALL this%dion_field_drho(i)%init(cell)
                END DO
                !
            ELSE
                CALL io%error(sub_name, "field-aware not implemented for specified mode", 1)
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
        copy%field_asymmetry = this%field_asymmetry
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
        IF (this%soft_spheres%number /= 0) THEN
            !
            IF (copy%soft_spheres%number /= 0) CALL copy%soft_spheres%destroy()
            !
            CALL this%soft_spheres%copy(copy%soft_spheres)
            !
        ELSE
            IF (copy%soft_spheres%number /= 0) DEALLOCATE (copy%soft_spheres%array)
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
        CHARACTER(LEN=80) :: sub_name = 'update_environ_boundary'
        !
        !--------------------------------------------------------------------------------
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
                CALL this%ions%core_electrons%density(this%ions%core, .TRUE.)
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
                                  "Wrong update status, possibly missing ionic update", 1)
                !
                this%density%of_r = this%electrons%density%of_r + this%ions%core%of_r
                !
                CALL this%boundary_of_density()
                !
                this%update_status = 2 ! boundary has changed and is ready
            END IF
            !
        CASE ('electronic')
            !
            IF (this%electrons%lupdate) THEN
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
            IF (this%field_aware) THEN
                !
                IF (this%ions%lupdate) THEN
                    !
                    CALL this%calc_dion_field_drho()
                    !
                    this%update_status = 1
                ELSE IF (this%electrons%lupdate) THEN
                    !
                    CALL this%calc_ion_field()
                    !
                    CALL this%update_soft_spheres(this%field_aware)
                    !
                    CALL this%boundary_of_functions()
                    !
                    this%update_status = 2
                END IF
                !
            ELSE IF (this%ions%lupdate) THEN
                !
                !------------------------------------------------------------------------
                ! Only ions are needed, fully update the boundary
                !
                CALL this%soft_spheres%update(this%ions%number, this%ions%tau)
                !
                CALL this%boundary_of_functions()
                !
                this%update_status = 2 ! boundary has changed and is ready
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
            CALL io%error(sub_name, "Unrecognized boundary mode", 1)
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
        CHARACTER(LEN=80) :: sub_name = 'destroy_environ_boundary'
        !
        INTEGER :: i
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
            !
            IF (this%mode == 'ionic') THEN
                !
                IF (ALLOCATED(this%denloc)) THEN
                    !
                    DO i=1,this%soft_spheres%number
                        !
                        CALL this%denloc(i)%destroy()
                        CALL this%gradloc(i)%destroy()
                        !
                    END DO
                    !
                    DEALLOCATE(this%denloc)
                    DEALLOCATE(this%gradloc)
                    !
                END IF
                !
                CALL this%soft_spheres%destroy()
                !
                IF (this%field_aware) THEN
                    !
                    CALL this%unscaled_spheres%destroy()
                    !
                    DEALLOCATE (this%ion_field)
                    DEALLOCATE (this%dion_field_drho)
                    DEALLOCATE (this%partial_of_ion_field)
                END IF
                !
            END IF
            !
            NULLIFY (this%ions)
        END IF
        !
        IF (this%need_electrons) NULLIFY (this%electrons)
        !
        IF (this%solvent_aware) DEALLOCATE (this%solvent_probe%pos)
        !
        IF (this%need_system) THEN
            !
            CALL this%simple%destroy()
            !
            NULLIFY (this%system)
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
        CLASS(environ_boundary), INTENT(IN) :: this
        REAL(DP), INTENT(IN) :: confine
        !
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
        CLASS(environ_boundary), INTENT(IN) :: this
        REAL(DP), INTENT(IN) :: pressure
        !
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
        CLASS(environ_boundary), INTENT(IN) :: this
        REAL(DP), INTENT(IN) :: surface_tension
        !
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
        CLASS(environ_boundary), INTENT(IN) :: this
        REAL(DP), INTENT(IN) :: surface_tension
        !
        TYPE(environ_density), INTENT(INOUT) :: de_dboundary
        !
        !--------------------------------------------------------------------------------
        !
        de_dboundary%of_r = de_dboundary%of_r + &
                            surface_tension * this%dsurface%of_r * e2 / 2.D0
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_desurface_dboundary
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_dboundary_dions(this, index, partial)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_boundary), TARGET, INTENT(IN) :: this
        INTEGER, INTENT(IN) :: index
        !
        TYPE(environ_gradient), INTENT(INOUT) :: partial
        !
        REAL(DP), PARAMETER :: tolspuriousforce = 1.D-5
        !
        INTEGER, POINTER :: number
        !
        INTEGER :: i, j
        REAL(DP) :: spurious_force
        TYPE(environ_density) :: denlocal
        !
        CHARACTER(LEN=80) :: sub_name = 'calc_dboundary_dions'
        !
        !--------------------------------------------------------------------------------
        !
        IF (this%mode == 'electronic') RETURN
        ! exit if boundary is only defined on electronic density
        !
        IF (this%need_ions) THEN
            number => this%ions%number
        ELSE IF (this%need_system) THEN
            number => this%system%ions%number
        ELSE
            CALL io%error(sub_name, "Missing details of ions", 1)
        END IF
        !
        IF (index > number) &
            CALL io%error(sub_name, "Index greater than number of ions", 1)
        !
        IF (index <= 0) &
            CALL io%error(sub_name, "Index of ion is zero or lower", 1)
        !
        IF (this%mode == 'ionic' .AND. this%soft_spheres%number == 0) &
            CALL io%error(sub_name, "Missing details of ionic boundary", 1)
        !
        IF (this%mode == 'full') THEN
            !
            IF (this%ions%core_electrons%number == 0) &
                CALL io%error(sub_name, "Missing details of core electrons", 1)
            !
            IF (.NOT. ASSOCIATED(this%dscaled%cell, partial%cell)) &
                CALL io%error(sub_name, "Mismatch or unassociated boundary derivative", 1)
            !
        END IF
        !
        IF (this%mode == 'ionic' .OR. this%mode == 'fa-ionic') THEN
            !
            IF (this%derivatives_method == 'fft') THEN
                CALL this%soft_spheres%array(index)%gradient(partial, .TRUE.)
            ELSE
                partial%of_r = this%gradloc(index)%of_r
            END IF
            !
            DO i = 1, number
                !
                IF (i == index) CYCLE
                !
                DO j = 1, 3
                    partial%of_r(j, :) = partial%of_r(j, :) * this%denloc(i)%of_r
                END DO
                !
            END DO
            !
        ELSE IF (this%mode == 'full') THEN
            !
            CALL this%ions%core_electrons%array(index)%gradient(partial, .TRUE.)
            !
            DO j = 1, 3
                partial%of_r(j, :) = -partial%of_r(j, :) * this%dscaled%of_r
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
1000    FORMAT(" WARNING: Unphysical forces due to core electrons are non-negligible ", /, &
               " atom type ", I3, " is subject to a spurious force of ", F12.6)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_dboundary_dions
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
        CLASS(environ_boundary), INTENT(IN) :: this
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
        ! Step 3: update the functional derivative of the energy w.r.t boundary
        !
        de_dboundary%of_r = de_dboundary%of_r * (1.D0 - this%filling%of_r) + local%of_r
        !
        CALL local%destroy()
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_solvent_aware_de_dboundary
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
        TYPE(environ_density), OPTIONAL, TARGET, INTENT(IN) :: density
        !
        CLASS(environ_boundary), TARGET, INTENT(INOUT) :: this
        !
        INTEGER, POINTER :: stype
        REAL(DP), POINTER :: rhomax, rhomin, tbeta, eps
        !
        TYPE(environ_density), POINTER :: denloc
        TYPE(environ_hessian), POINTER :: hessloc
        !
        INTEGER :: i, j
        !
        CHARACTER(LEN=80) :: sub_name = 'boundary_of_density'
        !
        !--------------------------------------------------------------------------------
        !
        IF (PRESENT(density)) THEN
            denloc => density
        ELSE
            denloc => this%density
        END IF
        !
        IF (.NOT. ASSOCIATED(denloc%cell, this%scaled%cell)) &
            CALL io%error(sub_name, "Inconsistent domains", 1)
        !
        stype => this%b_type
        !
        IF (stype == 1 .OR. stype == 2) THEN
            rhomax => this%rhomax
            rhomin => this%rhomin
            tbeta => this%fact
            eps => this%const
        ELSE IF (stype == 0) THEN
            rhomax => this%rhozero
            rhomin => this%deltarho
            tbeta => this%tbeta
            eps => this%const
        END IF
        !
        !--------------------------------------------------------------------------------
        !
        ASSOCIATE (cell => denloc%cell, &
                   deriv => this%deriv, &
                   derivatives => this%cores%derivatives, &
                   rho => denloc%of_r, &
                   scal => this%scaled, &
                   dscal => this%dscaled, &
                   d2scal => this%d2scaled, &
                   grad => this%gradient, &
                   lapl => this%laplacian, &
                   dsurf => this%dsurface)
            !
            !----------------------------------------------------------------------------
            !
            DO i = 1, cell%ir_end
                scal%of_r(i) = boundfunct(rho(i), rhomax, rhomin, tbeta, eps, stype)
                dscal%of_r(i) = dboundfunct(rho(i), rhomax, rhomin, tbeta, eps, stype)
                d2scal%of_r(i) = d2boundfunct(rho(i), rhomax, rhomin, tbeta, eps, stype)
            END DO
            !
            !----------------------------------------------------------------------------
            ! Compute boundary derivatives, if needed
            !
            IF (deriv >= 3) THEN
                !
                IF (this%solvent_aware) THEN
                    hessloc => this%hessian
                ELSE
                    ALLOCATE (hessloc)
                    !
                    CALL hessloc%init(cell)
                    !
                END IF
                !
            END IF
            !
            SELECT CASE (this%derivatives_method)
                !
            CASE ('fft')
                !
                IF (deriv == 1 .OR. deriv == 2) CALL derivatives%gradient(scal, grad)
                !
                IF (deriv == 2) CALL derivatives%laplacian(scal, lapl)
                !
                IF (deriv == 3) CALL this%calc_dsurface(scal, grad, lapl, hessloc, dsurf)
                !
            CASE ('chain')
                !
                IF (deriv == 1 .OR. deriv == 2) CALL derivatives%gradient(denloc, grad)
                !
                IF (deriv == 2) CALL derivatives%laplacian(denloc, lapl)
                !
                IF (deriv == 3) THEN
                    !
                    CALL this%calc_dsurface(denloc, grad, lapl, hessloc, dsurf)
                    !
                    IF (this%solvent_aware) THEN
                        !
                        DO i = 1, 3
                            !
                            DO j = 1, 3
                                !
                                hessloc%of_r(i, j, :) = &
                                    hessloc%of_r(i, j, :) * dscal%of_r + &
                                    grad%of_r(i, :) * grad%of_r(j, :) * d2scal%of_r
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
                    lapl%of_r = lapl%of_r * dscal%of_r + &
                                (grad%of_r(1, :)**2 + &
                                 grad%of_r(2, :)**2 + &
                                 grad%of_r(3, :)**2) * d2scal%of_r
                !
                IF (deriv >= 1) THEN
                    !
                    DO i = 1, 3
                        grad%of_r(i, :) = grad%of_r(i, :) * dscal%of_r
                    END DO
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
            this%volume = scal%integrate()
            !
            IF (deriv >= 1) THEN
                !
                CALL grad%update_modulus()
                !
                this%surface = grad%modulus%integrate()
            END IF
            !
            IF (deriv >= 3 .AND. .NOT. this%solvent_aware) CALL hessloc%destroy()
            !
        END ASSOCIATE
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
        INTEGER :: i
        !
        TYPE(environ_density), ALLOCATABLE :: laplloc(:)
        TYPE(environ_hessian), ALLOCATABLE :: hessloc(:)
        TYPE(environ_hessian), POINTER :: hess
        !
        CHARACTER(LEN=80) :: sub_name = 'boundary_of_functions'
        !
        !--------------------------------------------------------------------------------
        !
        ASSOCIATE (cell => this%scaled%cell, &
                   nss => this%soft_spheres%number, &
                   soft_spheres => this%soft_spheres%array, &
                   deriv => this%deriv, &
                   derivatives => this%cores%derivatives, &
                   scal => this%scaled, &
                   grad => this%gradient, &
                   lapl => this%laplacian, &
                   dsurf => this%dsurface)
            !
            !----------------------------------------------------------------------------
            ! Compute soft spheres and generate boundary
            !
            scal%of_r = 1.D0
            !
            DO i = 1, nss
                !
                CALL soft_spheres(i)%density(this%denloc(i), .TRUE.)
                !
                scal%of_r = scal%of_r * this%denloc(i)%of_r
            END DO
            !
            !----------------------------------------------------------------------------
            ! Generate boundary derivatives, if needed
            !
            IF (deriv == 3) THEN
                !
                IF (this%solvent_aware) THEN
                    hess => this%hessian
                    hess%of_r = 0.D0
                ELSE
                    ALLOCATE (hess)
                    !
                    CALL hess%init(cell)
                    !
                END IF
                !
            END IF
            !
            SELECT CASE (this%derivatives_method)
                !
            CASE ('fft')
                !
                IF (deriv == 1 .OR. deriv == 2) CALL derivatives%gradient(scal, grad)
                !
                IF (deriv == 2) CALL derivatives%laplacian(scal, lapl)
                !
                IF (deriv == 3) CALL this%calc_dsurface(scal, grad, lapl, hess, dsurf)
                !
            CASE ('highmem')
                !
                IF (deriv == 2) ALLOCATE (laplloc(nss))
                !
                IF (deriv == 3) ALLOCATE (hessloc(nss))
                !
                !------------------------------------------------------------------------
                ! Compute and temporarily store soft spheres derivatives
                !
                DO i = 1, nss
                    !
                    IF (deriv == 2) CALL laplloc(i)%init(cell)
                    !
                    IF (deriv == 3) CALL hessloc(i)%init(cell)
                    !
                    IF (deriv >= 1) CALL soft_spheres(i)%gradient(this%gradloc(i), .TRUE.)
                    !
                    IF (deriv == 2) CALL soft_spheres(i)%laplacian(laplloc(i), .FALSE.)
                    !
                    IF (deriv == 3) CALL soft_spheres(i)%hessian(hessloc(i), .FALSE.)
                    !
                END DO
                !
                IF (deriv == 1 .OR. deriv == 2) &
                    CALL gradient_of_boundary(nss, this%denloc, this%gradloc, grad)
                !
                IF (deriv == 2) &
                    CALL laplacian_of_boundary(nss, this%denloc, this%gradloc, laplloc, lapl)
                !
                IF (deriv == 3) &
                    CALL dsurface_of_boundary(nss, this%denloc, this%gradloc, hessloc, grad, &
                                              lapl, hess, dsurf)
                !
                DO i = 1, nss
                    !
                    IF (deriv == 2) CALL laplloc(i)%destroy()
                    !
                    IF (deriv == 3) CALL hessloc(i)%destroy()
                    !
                END DO
                !
                IF (deriv == 2) DEALLOCATE (laplloc)
                !
                IF (deriv == 3) DEALLOCATE (hessloc)
                !
            CASE ('lowmem')
                !
                IF (deriv == 2) ALLOCATE (laplloc(nss))
                !
                IF (deriv == 3) ALLOCATE (hessloc(nss))
                !
                !------------------------------------------------------------------------
                ! Compute and store soft spheres derivatives
                !
                DO i = 1, nss
                    !
                    IF (deriv == 2) CALL laplloc(i)%init(cell)
                    !
                    IF (deriv == 3) CALL hessloc(i)%init(cell)
                    !
                    IF (deriv >= 1) CALL soft_spheres(i)%gradient(this%gradloc(i), .TRUE.)
                    !
                    IF (deriv == 2) CALL soft_spheres(i)%laplacian(laplloc(i), .FALSE.)
                    !
                    IF (deriv == 3) CALL soft_spheres(i)%hessian(hessloc(i), .FALSE.)
                    !
                END DO
                !
                IF (deriv >= 1) &
                    CALL gradient_of_boundary(nss, this%denloc, this%gradloc, scal, grad)
                !
                IF (deriv == 2) &
                    CALL laplacian_of_boundary(nss, this%denloc, this%gradloc, laplloc, scal, &
                                               grad, lapl)
                !
                IF (deriv == 3) &
                    CALL dsurface_of_boundary(nss, this%denloc, this%gradloc, hessloc, grad, &
                                              lapl, hess, scal, dsurf)
                !
                DO i = 1, nss
                    !
                    IF (deriv == 2) CALL laplloc(i)%destroy()
                    !
                    IF (deriv == 3) CALL hessloc(i)%destroy()
                    !
                END DO
                !
                IF (deriv == 2) DEALLOCATE (laplloc)
                !
                IF (deriv == 3) DEALLOCATE (hessloc)
                !
            CASE DEFAULT
                CALL io%error(sub_name, "Unexpected derivatives method", 1)
                !
            END SELECT
            !
            !----------------------------------------------------------------------------
            ! Final updates
            !
            scal%of_r = 1.D0 - scal%of_r
            this%volume = scal%integrate()
            !
            IF (deriv >= 1) THEN
                grad%of_r = -grad%of_r
                !
                CALL grad%update_modulus()
                !
                this%surface = grad%modulus%integrate()
                !
                IF (deriv >= 2) lapl%of_r = -lapl%of_r
                !
                IF (deriv == 3) THEN
                    dsurf%of_r = -dsurf%of_r
                    !
                    IF (this%solvent_aware) THEN
                        hess%of_r = -hess%of_r
                    ELSE
                        !
                        CALL hess%destroy()
                        !
                        DEALLOCATE (hess)
                    END IF
                    !
                END IF
                !
            END IF
            !
        END ASSOCIATE
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
        TYPE(environ_hessian), POINTER :: hessloc
        !
        CHARACTER(LEN=80) :: sub_name = 'boundary_of_system'
        !
        !--------------------------------------------------------------------------------
        !
        ASSOCIATE (cell => this%scaled%cell, &
                   deriv => this%deriv, &
                   derivatives => this%cores%derivatives, &
                   simple => this%simple, &
                   scal => this%scaled, &
                   grad => this%gradient, &
                   lapl => this%laplacian, &
                   dsurf => this%dsurface)
            !
            !----------------------------------------------------------------------------
            !
            CALL simple%density(scal, .TRUE.)
            ! compute soft spheres and generate boundary
            !
            !----------------------------------------------------------------------------
            ! Generate boundary derivatives, if needed
            !
            IF (deriv >= 3) THEN
                !
                IF (this%solvent_aware) THEN
                    hessloc => this%hessian
                ELSE
                    ALLOCATE (hessloc)
                    !
                    CALL hessloc%init(cell)
                    !
                END IF
                !
            END IF
            !
            SELECT CASE (this%derivatives_method)
                !
            CASE ('fft')
                !
                IF (deriv == 1 .OR. deriv == 2) CALL derivatives%gradient(scal, grad)
                !
                IF (deriv == 2) CALL derivatives%laplacian(scal, lapl)
                !
                IF (deriv == 3) CALL this%calc_dsurface(scal, grad, lapl, hessloc, dsurf)
                !
            CASE ('chain')
                !
                IF (deriv >= 1) CALL simple%gradient(grad, .TRUE.)
                !
                IF (deriv >= 2) CALL simple%laplacian(lapl, .TRUE.)
                !
                IF (deriv >= 3) THEN
                    !
                    CALL simple%hessian(hessloc, .TRUE.)
                    !
                    CALL calc_dsurface_no_pre(cell, grad%of_r, hessloc%of_r, dsurf%of_r)
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
                    CALL hessloc%destroy()
                    !
                    DEALLOCATE (hessloc)
                END IF
                !
            END IF
            !
            this%volume = scal%integrate()
            !
            IF (deriv >= 1) THEN
                !
                CALL grad%update_modulus()
                !
                this%surface = grad%modulus%integrate()
            END IF
            !
        END ASSOCIATE
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE boundary_of_system
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
        CLASS(environ_boundary), INTENT(INOUT) :: this
        !
        INTEGER :: i, j
        TYPE(environ_density) :: fillfrac
        TYPE(environ_density) :: d2fill
        !
        TYPE(environ_density) :: denloc
        TYPE(environ_gradient) :: gradloc
        TYPE(environ_density) :: laplloc
        TYPE(environ_hessian) :: hessloc
        !
        CHARACTER(LEN=80) :: sub_name = 'solvent_aware_boundary'
        !
        !--------------------------------------------------------------------------------
        !
        ASSOCIATE (cell => this%scaled%cell, &
                   deriv => this%deriv, &
                   derivatives => this%cores%derivatives, &
                   derivatives_method => this%derivatives_method, &
                   thr => this%filling_threshold, &
                   spr => this%filling_spread, &
                   fill => this%filling, &
                   dfill => this%dfilling, &
                   probe => this%probe, &
                   loc => this%local, &
                   scal => this%scaled, &
                   grad => this%gradient, &
                   lapl => this%laplacian, &
                   hess => this%hessian, &
                   dsurf => this%dsurface)
            !
            !----------------------------------------------------------------------------
            !
            CALL fillfrac%init(cell)
            !
            IF (deriv >= 2 .AND. derivatives_method /= 'fft') CALL d2fill%init(cell)
            !
            !----------------------------------------------------------------------------
            ! Step 0: save local interface function for later use
            !
            loc%of_r = scal%of_r
            !
            !----------------------------------------------------------------------------
            ! Step 1: compute the convolution function,
            !         this may be made moved out of here
            !
            CALL this%solvent_probe%density(probe, .TRUE.)
            !
            probe%of_r = probe%of_r / probe%integrate()
            !
            !----------------------------------------------------------------------------
            ! Step 2: compute filled fraction,
            !         i.e. convolution of local boundary with probe
            !
            CALL derivatives%convolution(loc, probe, fillfrac)
            !
            !----------------------------------------------------------------------------
            ! Step 3: compute the filling function and its derivative
            !
            fill%of_r = 0.D0
            dfill%of_r = 0.D0
            !
            DO i = 1, cell%ir_end
                fill%of_r(i) = 1.D0 - sfunct2(fillfrac%of_r(i), thr, spr)
                dfill%of_r(i) = -dsfunct2(fillfrac%of_r(i), thr, spr)
                !
                IF (deriv >= 2 .AND. derivatives_method /= 'fft') &
                    d2fill%of_r(i) = -d2sfunct2(fillfrac%of_r(i), thr, spr)
                !
            END DO
            !
            !----------------------------------------------------------------------------
            ! Step 4: compute solvent-aware interface
            !
            scal%of_r = loc%of_r + (1.D0 - loc%of_r) * fill%of_r
            !
            !----------------------------------------------------------------------------
            ! Step 5: compute boundary derivatives, if needed
            !
            SELECT CASE (derivatives_method)
                !
            CASE ('fft')
                !
                IF (deriv == 1 .OR. deriv == 2) CALL derivatives%gradient(scal, grad)
                !
                IF (deriv == 2) CALL derivatives%laplacian(scal, lapl)

                IF (deriv == 3) CALL this%calc_dsurface(scal, grad, lapl, hess, dsurf)
                !
            CASE ('chain', 'highmem', 'lowmem')
                !
                !------------------------------------------------------------------------
                ! Allocate local fields for derivatives of convolution
                !
                IF (deriv >= 1) CALL gradloc%init(cell)
                !
                IF (deriv >= 2) CALL laplloc%init(cell)
                !
                IF (deriv >= 3) CALL hessloc%init(cell)
                !
                !------------------------------------------------------------------------
                ! Compute derivative of convolution with probe
                !
                IF (deriv > 1) CALL this%convolution(deriv, gradloc, laplloc, hessloc)
                !
                !------------------------------------------------------------------------
                ! Update derivatives of interface function in reverse order
                !
                IF (deriv >= 3) THEN
                    !
                    DO i = 1, 3
                        !
                        DO j = 1, 3
                            !
                            hess%of_r(i, j, :) = &
                                hess%of_r(i, j, :) * (1.D0 - fill%of_r) - &
                                dfill%of_r * &
                                (grad%of_r(i, :) * gradloc%of_r(j, :) + &
                                 grad%of_r(j, :) * gradloc%of_r(i, :)) + &
                                (1.D0 - loc%of_r) * &
                                (d2fill%of_r * gradloc%of_r(i, :) * gradloc%of_r(j, :) + &
                                 dfill%of_r * hessloc%of_r(i, j, :))
                            !
                        END DO
                        !
                    END DO
                    !
                    CALL hessloc%destroy()
                    !
                END IF
                !
                IF (deriv >= 2) THEN
                    !
                    CALL denloc%init(cell)
                    !
                    CALL grad%scalar_product(gradloc, denloc)
                    !
                    lapl%of_r = &
                        lapl%of_r * (1.D0 - fill%of_r) - &
                        2.D0 * denloc%of_r * dfill%of_r + &
                        (1.D0 - loc%of_r) * (d2fill%of_r * gradloc%modulus%of_r**2 + &
                                             dfill%of_r * laplloc%of_r)
                    !
                    CALL denloc%destroy()
                    !
                    CALL laplloc%destroy()
                    !
                    CALL d2fill%destroy()
                    !
                END IF
                !
                IF (deriv >= 1) THEN
                    !
                    DO i = 1, 3
                        !
                        grad%of_r(i, :) = &
                            grad%of_r(i, :) * (1.D0 - fill%of_r) + &
                            gradloc%of_r(i, :) * (1.D0 - loc%of_r) * dfill%of_r
                        !
                    END DO
                    !
                    CALL gradloc%destroy()
                    !
                END IF
                !
                !------------------------------------------------------------------------
                ! Recompute dsurface, if needed
                !
                IF (deriv >= 3) &
                    CALL calc_dsurface_no_pre(cell, grad%of_r, hess%of_r, dsurf%of_r)
                !
            CASE DEFAULT
                CALL io%error(sub_name, "Unexpected derivatives method", 1)
                !
            END SELECT
            !
            !----------------------------------------------------------------------------
            ! Final updates
            !
            this%volume = scal%integrate()
            !
            IF (deriv >= 1) THEN
                !
                CALL grad%update_modulus()
                !
                this%surface = grad%modulus%integrate()
            END IF
            !
            CALL fillfrac%destroy()
            !
        END ASSOCIATE
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE solvent_aware_boundary
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !                                   FIELD-AWARE METHODS
    !
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE update_soft_spheres(this, field_scaling)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        LOGICAL, INTENT(IN), OPTIONAL :: field_scaling
        !
        CLASS(environ_boundary), INTENT(INOUT) :: this
        !
        INTEGER :: i
        REAL(DP) :: field_scale
        !
        CHARACTER(LEN=80) :: sub_name = 'update_soft_spheres'
        !
        !--------------------------------------------------------------------------------
        !
        DO i = 1, this%ions%number
            !
            ASSOCIATE (soft_sphere => this%soft_spheres%array(i), &
                       solvationrad => this%ions%iontype(this%ions%ityp(i))%solvationrad)
                !
                !------------------------------------------------------------------------
                ! field-aware scaling of soft-sphere radii
                !
                IF (PRESENT(field_scaling)) THEN
                    !
                    IF (field_scaling) THEN
                        field_scale = this%scaling_of_field(i)
                    ELSE
                        field_scale = 1.D0
                    END IF
                ELSE
                    field_scale = 1.D0
                END IF
                !
                soft_sphere%pos = this%ions%tau(:, i)
                soft_sphere%width = solvationrad * this%alpha * field_scale
                !
            END ASSOCIATE
            !
        END DO
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE update_soft_spheres
    !------------------------------------------------------------------------------------
    !>
    !! Computes the flux due to the ions
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_ion_field(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_boundary), INTENT(INOUT) :: this
        !
        INTEGER :: i, j
        !
        TYPE(environ_density), ALLOCATABLE :: local(:)
        !
        TYPE(environ_density) :: aux, prod
        TYPE(environ_gradient) :: auxg, field
        !
        CHARACTER(LEN=80) :: sub_name = 'calc_ion_field'
        !
        !--------------------------------------------------------------------------------
        !
        ASSOCIATE (cell => this%scaled%cell, &
                   n => this%ions%number, &
                   electrostatics => this%cores%electrostatics)
            !
            !----------------------------------------------------------------------------
            !
            ALLOCATE (local(n))
            !
            !----------------------------------------------------------------------------
            !
            DO i = 1, n
                !
                CALL local(i)%init(cell)
                !
                CALL this%unscaled_spheres%array(i)%density(local(i), .FALSE.)
                !
            END DO
            !
            !----------------------------------------------------------------------------
            ! Compute field
            !
            CALL aux%init(cell)
            !
            aux%of_r = this%electrons%density%of_r + this%ions%density%of_r
            !
            CALL field%init(cell)
            !
            CALL electrostatics%grad_v_h_of_rho_r(cell%nnr, aux%of_r, field%of_r)
            !
            !----------------------------------------------------------------------------
            ! Compute ion flux
            !
            this%ion_field = 0.D0
            !
            CALL prod%init(cell)
            !
            CALL auxg%init(cell)
            !
            DO i = 1, n
                prod%of_r = 1.D0
                !
                DO j = 1, n
                    !
                    IF (i == j) CYCLE
                    !
                    prod%of_r = prod%of_r * local(j)%of_r
                END DO
                !
                !------------------------------------------------------------------------
                ! Compute field flux through soft-sphere interface
                !
                CALL this%unscaled_spheres%array(i)%gradient(auxg, .TRUE.)
                !
                CALL field%scalar_product(auxg, aux)
                !
                aux%of_r = -aux%of_r * prod%of_r
                this%ion_field(i) = aux%integrate()
            END DO
            !
            CALL auxg%destroy()
            !
            CALL prod%destroy()
            !
            CALL field%destroy()
            !
            CALL aux%destroy()
            !
        END ASSOCIATE
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_ion_field
    !------------------------------------------------------------------------------------
    !>
    !! Computes the derivative of the flux due to the ions w.r.t ionic position
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_ion_field_partial(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_boundary), INTENT(INOUT) :: this
        !
        INTEGER :: i, j, k
        !
        TYPE(environ_density) :: aux, prod
        TYPE(environ_gradient) :: auxg, field
        TYPE(environ_hessian) :: hessloc, auxh
        !
        TYPE(environ_density), ALLOCATABLE :: local(:)
        TYPE(environ_gradient), ALLOCATABLE :: gradloc(:)
        REAL(DP), ALLOCATABLE :: ion_field(:)
        !
        CHARACTER(LEN=80) :: sub_name = 'calc_ion_field_partial'
        !
        !--------------------------------------------------------------------------------
        !
        ASSOCIATE (cell => this%scaled%cell, &
                   n => this%ions%number, &
                   electrostatics => this%cores%electrostatics)
            !
            !----------------------------------------------------------------------------
            !
            ALLOCATE (local(n))
            ALLOCATE (gradloc(n))
            ALLOCATE (ion_field(n))
            !
            !----------------------------------------------------------------------------
            !
            DO i = 1, n
                !
                CALL local(i)%init(cell)
                !
                CALL gradloc(i)%init(cell)
                !
                CALL this%soft_spheres%array(i)%density(local(i), .FALSE.)
                !
                CALL this%soft_spheres%array(i)%gradient(gradloc(i), .FALSE.)
                !
            END DO
            !
            CALL hessloc%init(cell)
            !
            !----------------------------------------------------------------------------
            ! Compute field
            !
            CALL aux%init(cell)
            !
            aux%of_r = this%electrons%density%of_r + this%ions%density%of_r
            !
            CALL field%init(cell)
            !
            CALL electrostatics%grad_v_h_of_rho_r(cell%nnr, aux%of_r, field%of_r)
            !
            !----------------------------------------------------------------------------
            ! Compute field flux
            !
            ion_field = 0.D0
            this%partial_of_ion_field = 0.D0
            !
            CALL prod%init(cell)
            !
            CALL auxg%init(cell)
            !
            CALL auxh%init(cell)
            !
            DO i = 1, n
                prod%of_r = 1.D0
                !
                DO j = 1, n
                    !
                    IF (i == j) CYCLE
                    !
                    prod%of_r = prod%of_r * local(j)%of_r
                END DO
                !
                CALL field%scalar_product(gradloc(i), aux) ! here aux is the normal field
                !
                aux%of_r = -aux%of_r * prod%of_r
                ion_field(i) = aux%integrate()
                !
                DO j = 1, n
                    !
                    !--------------------------------------------------------------------
                    ! This is pretty ugly, is there a faster way to implement this?
                    !
                    CALL this%ions%smeared_ions%array(j)%density(aux, .TRUE.)
                    !
                    CALL electrostatics%hess_v_h_of_rho_r(cell%nnr, aux%of_r, &
                                                          hessloc%of_r)
                    !
                    CALL hessloc%scalar_product(gradloc(i), auxg)
                    !
                    this%partial_of_ion_field(:, i, j) = &
                        this%partial_of_ion_field(:, i, j) - &
                        auxg%scalar_product_density(prod)
                    !
                    IF (i == j) THEN
                        !
                        !----------------------------------------------------------------
                        ! Hessian of soft-sphere times the field
                        !
                        CALL this%soft_spheres%array(i)%hessian(auxh, .TRUE.)
                        !
                        CALL auxh%scalar_product(field, auxg)
                        !
                        this%partial_of_ion_field(:, i, j) = &
                            this%partial_of_ion_field(:, i, j) + &
                            auxg%scalar_product_density(prod)
                        !
                    ELSE
                        !
                        !----------------------------------------------------------------
                        ! Ion field times gradient of differential soft-sphere
                        !
                        CALL gradloc(i)%scalar_product(field, aux)
                        !
                        DO k = 1, n
                            !
                            IF (i == k) CYCLE
                            !
                            IF (j == k) CYCLE
                            !
                            aux%of_r = aux%of_r * local(k)%of_r
                        END DO
                        !
                        this%partial_of_ion_field(:, i, j) = &
                            this%partial_of_ion_field(:, i, j) + &
                            gradloc(j)%scalar_product_density(aux)
                        !
                    END IF
                    !
                END DO
                !
            END DO
            !
            CALL field%destroy()
            !
            CALL prod%destroy()
            !
            CALL aux%destroy()
            !
            CALL auxg%destroy()
            !
            CALL auxh%destroy()
            !
            CALL hessloc%destroy()
            !
            DO i = 1, n
                !
                CALL local(i)%destroy()
                !
                CALL gradloc(i)%destroy()
                !
            END DO
            !
        END ASSOCIATE
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_ion_field_partial
    !------------------------------------------------------------------------------------
    !>
    !! Computes the functional derivative of the flux due to the ions w.r.t the
    !! electronic density
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_dion_field_drho(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_boundary), INTENT(INOUT) :: this
        !
        INTEGER :: i, j, k
        !
        TYPE(environ_density) :: prod
        TYPE(environ_gradient) :: auxg
        !
        TYPE(environ_density), ALLOCATABLE :: local(:)
        !
        CHARACTER(LEN=80) :: sub_name = 'calc_dion_field_drho'
        !
        !--------------------------------------------------------------------------------
        !
        ASSOCIATE (cell => this%scaled%cell, &
                   n => this%ions%number, &
                   electrostatics => this%cores%electrostatics)
            !
            !----------------------------------------------------------------------------
            !
            ALLOCATE (local(n))
            !
            !----------------------------------------------------------------------------
            !
            DO i = 1, n
                !
                CALL local(i)%init(cell)
                !
                CALL this%unscaled_spheres%array(i)%density(local(i), .FALSE.)
                !
            END DO
            !
            !----------------------------------------------------------------------------
            ! Compute field flux
            !
            CALL prod%init(cell)
            !
            CALL auxg%init(cell)
            !
            DO i = 1, n
                !
                !------------------------------------------------------------------------
                ! Compute product of other soft-spheres
                !
                prod%of_r = 1.D0
                !
                DO j = 1, n
                    !
                    IF (i == j) CYCLE
                    !
                    prod%of_r = prod%of_r * local(j)%of_r
                END DO
                !
                !------------------------------------------------------------------------
                ! Compute functional derivative of field w.r.t electric density
                !
                CALL this%unscaled_spheres%array(i)%gradient(auxg, .TRUE.)
                !
                DO k = 1, 3
                    auxg%of_r(k, :) = auxg%of_r(k, :) * prod%of_r
                END DO
                !
                CALL electrostatics%field_of_grad_rho(cell%nnr, auxg%of_r, &
                                                      this%dion_field_drho(i)%of_r)
                !
            END DO
            !
            CALL auxg%destroy()
            !
            DO i = 1, n
                CALL local(i)%destroy()
            END DO
            !
        END ASSOCIATE
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_dion_field_drho
    !------------------------------------------------------------------------------------
    !>
    !! Computes the functional derivative of the energy w.r.t the electronic density
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_field_aware_de_drho(this, de_dboundary, de_drho)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_density), INTENT(IN) :: de_dboundary
        !
        CLASS(environ_boundary), INTENT(INOUT) :: this
        TYPE(environ_density), INTENT(INOUT) :: de_drho
        !
        INTEGER :: i, j
        REAL(DP) :: df
        !
        TYPE(environ_density) :: aux
        !
        TYPE(environ_density), ALLOCATABLE :: local(:)
        !
        REAL(DP), POINTER :: solvationrad
        !
        CHARACTER(LEN=80) :: sub_name = 'calc_field_aware_de_drho'
        !
        !--------------------------------------------------------------------------------
        !
        ASSOCIATE (cell => this%scaled%cell, &
                   n => this%ions%number)
            !
            !----------------------------------------------------------------------------
            !
            ALLOCATE (local(n))
            !
            !----------------------------------------------------------------------------
            !
            IF (this%mode == 'ionic') THEN
                !
                DO i = 1, n
                    !
                    CALL local(i)%init(cell)
                    !
                    CALL this%soft_spheres%array(i)%density(local(i), .FALSE.)
                    !
                END DO
                !
                CALL aux%init(cell)
                !
                DO i = 1, n
                    solvationrad => this%ions%iontype(this%ions%ityp(i))%solvationrad
                    !
                    CALL this%soft_spheres%array(i)%derivative(aux, .TRUE.)
                    !
                    DO j = 1, n
                        !
                        IF (i == j) CYCLE
                        !
                        aux%of_r = aux%of_r * local(j)%of_r
                    END DO
                    !
                    df = this%dscaling_of_field(i) * solvationrad * this%alpha * &
                         aux%scalar_product(de_dboundary)
                    !
                    de_drho%of_r = de_drho%of_r + this%dion_field_drho(i)%of_r * df
                END DO
                !
                CALL aux%destroy()
                !
                DO i = 1, n
                    CALL local(i)%destroy()
                END DO
                !
            ELSE
                CALL io%error(sub_name, "boundary mode not implemented", 1)
            END IF
            !
        END ASSOCIATE
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_field_aware_de_drho
    !------------------------------------------------------------------------------------
    !>
    !! Computes the functional derivative of the boundary w.r.t the ionic positions
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_field_aware_dboundary_dions(this, index, partial)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: index
        !
        CLASS(environ_boundary), INTENT(INOUT) :: this
        TYPE(environ_gradient), INTENT(INOUT) :: partial
        !
        INTEGER :: i, j, k
        REAL(DP) :: df
        !
        TYPE(environ_density) :: aux
        TYPE(environ_gradient) :: auxg
        !
        TYPE(environ_density), ALLOCATABLE :: local(:)
        !
        REAL(DP), POINTER :: solvationrad
        !
        CHARACTER(LEN=80) :: sub_name = 'calc_field_aware_dboundary_dions'
        !
        !--------------------------------------------------------------------------------
        !
        ASSOCIATE (cell => this%scaled%cell, &
                   n => this%ions%number)
            !
            !----------------------------------------------------------------------------
            !
            ALLOCATE (local(n))
            !
            !----------------------------------------------------------------------------
            !
            IF (this%mode == 'ionic') THEN
                !
                DO i = 1, n
                    !
                    CALL local(i)%init(cell)
                    !
                    CALL this%soft_spheres%array(i)%density(local(i), .FALSE.)
                    !
                END DO
                !
                CALL aux%init(cell)
                !
                CALL auxg%init(cell)
                !
                DO i = 1, n
                    solvationrad => this%ions%iontype(this%ions%ityp(i))%solvationrad
                    !
                    CALL this%soft_spheres%array(i)%derivative(aux, .TRUE.)
                    !
                    DO j = 1, n
                        !
                        IF (i == j) CYCLE
                        aux%of_r = aux%of_r + local(j)%of_r
                        !
                    END DO
                    !
                    df = this%dscaling_of_field(i) * solvationrad * this%alpha
                    aux%of_r = aux%of_r * df
                    !
                    DO k = 1, 3
                        !
                        auxg%of_r(k, :) = &
                            auxg%of_r(k, :) + &
                            aux%of_r * this%partial_of_ion_field(k, i, index)
                        !
                    END DO
                    !
                END DO
                !
                partial%of_r = partial%of_r * auxg%of_r
                !
                CALL aux%destroy()
                !
                CALL auxg%destroy()
                !
                DO i = 1, n
                    CALL local(i)%destroy()
                END DO
                !
            ELSE
                CALL io%error(sub_name, "boundary mode not implemented", 1)
            END IF
            !
        END ASSOCIATE
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_field_aware_dboundary_dions
    !------------------------------------------------------------------------------------
    !>
    !! Returns field-aware scaling function with given ion_field and field aware
    !! boundary parameters
    !!
    !------------------------------------------------------------------------------------
    FUNCTION scaling_of_field(this, i) RESULT(scaling)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_boundary), INTENT(IN) :: this
        INTEGER, INTENT(IN) :: i
        !
        REAL(DP) :: scaling, multiplier, arg, diff
        !
        CHARACTER(LEN=80) :: fun_name = 'scaling_of_field'
        !
        !--------------------------------------------------------------------------------
        !
        multiplier = (this%field_asymmetry - SIGN(1.D0, this%ion_field(i)))**2 * &
                     this%field_factor
        !
        IF (ABS(this%ion_field(i)) < this%field_min) THEN
            scaling = 0.D0
        ELSE IF (ABS(this%ion_field(i)) > this%field_max) THEN
            scaling = 1.D0
        ELSE
            diff = this%field_max - this%field_min
            arg = tpi * (ABS(this%ion_field(i)) - this%field_min) / diff
            scaling = (arg - SIN(arg)) / tpi
        END IF
        !
        scaling = 1.D0 - scaling * multiplier
        !
        !--------------------------------------------------------------------------------
    END FUNCTION scaling_of_field
    !------------------------------------------------------------------------------------
    !>
    !! Returns field-aware scaling function with given ion_field and field aware
    !! boundary parameters
    !!
    !------------------------------------------------------------------------------------
    FUNCTION dscaling_of_field(this, i) RESULT(dscaling)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_boundary), INTENT(IN) :: this
        INTEGER, INTENT(IN) :: i
        !
        REAL(DP) :: dscaling, multiplier, arg, diff
        !
        !--------------------------------------------------------------------------------
        !
        multiplier = (this%field_asymmetry - SIGN(1.D0, this%ion_field(i)))**2 * &
                     this%field_factor
        !
        IF (ABS(this%ion_field(i)) < this%field_min) THEN
            dscaling = 0.D0
        ELSE IF (ABS(this%ion_field(i)) > this%field_max) THEN
            dscaling = 0.D0
        ELSE
            diff = this%field_max - this%field_min
            arg = tpi * (ABS(this%ion_field(i)) - this%field_min) / diff
            dscaling = (1.D0 - COS(arg)) / diff
        END IF
        !
        dscaling = -dscaling * multiplier * SIGN(1.D0, this%ion_field(i))
        !
        !--------------------------------------------------------------------------------
    END FUNCTION dscaling_of_field
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
        INTEGER, DIMENSION(this%ions%number) :: axes, dims
        REAL(DP), DIMENSION(this%ions%number) :: spreads, volumes
        !
        REAL(DP), ALLOCATABLE :: radii(:)
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
        CALL this%soft_spheres%init(this%ions%number, 4, axes, dims, radii, spreads, &
                                    volumes, this%ions%tau)
        !
        IF (this%field_aware) CALL this%soft_spheres%copy(this%unscaled_spheres)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE set_soft_spheres
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE compute_convolution_deriv(this, deriv, grad, lapl, hess)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_boundary), INTENT(IN) :: this
        INTEGER, INTENT(IN) :: deriv
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
        ASSOCIATE (derivatives => this%cores%derivatives, &
                   probe => this%probe)
            !
            IF (deriv >= 1) THEN
                !
                CALL derivatives%convolution(probe, this%gradient, grad)
                !
                CALL grad%update_modulus()
                !
            END IF
            !
            IF (deriv >= 2) CALL derivatives%convolution(probe, this%laplacian, lapl)
            !
            IF (deriv >= 3) CALL derivatives%convolution(probe, this%hessian, hess)
            !
        END ASSOCIATE
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE compute_convolution_deriv
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_dsurface(this, dens, grad, lapl, hess, dsurf)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_boundary), INTENT(IN) :: this
        TYPE(environ_density), INTENT(IN) :: dens
        !
        TYPE(environ_gradient), INTENT(INOUT) :: grad
        TYPE(environ_density), INTENT(INOUT) :: lapl
        TYPE(environ_hessian), INTENT(INOUT) :: hess
        TYPE(environ_density), INTENT(INOUT) :: dsurf
        !
        CHARACTER(LEN=80) :: sub_name = 'calc_dsurface'
        !
        !--------------------------------------------------------------------------------
        !
        CALL this%cores%derivatives%hessian(dens, grad, hess)
        !
        lapl%of_r = hess%of_r(1, 1, :) + hess%of_r(2, 2, :) + hess%of_r(3, 3, :)
        !
        CALL calc_dsurface_no_pre(dens%cell, grad%of_r, hess%of_r, dsurf%of_r)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_dsurface
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
        INTEGER, OPTIONAL, INTENT(IN) :: verbose, debug_verbose, unit
        !
        INTEGER :: base_verbose, local_verbose, passed_verbose, local_unit, i
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
                        CALL io%error(sub_name, "Unexpected boundary type", 1)
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
                    CALL this%soft_spheres%printout(passed_verbose, debug_verbose, &
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
            IF (this%field_aware) THEN
                !
                IF (io%lnode .AND. local_verbose >= 1) THEN
                    !
                    WRITE (local_unit, 1113)
                    !
                    DO i = 1, this%ions%number
                        !
                        WRITE (local_unit, 1114) i, &
                            this%ions%iontype(this%ions%ityp(i))%label, &
                            this%ions%iontype(this%ions%ityp(i))%solvationrad, &
                            this%ion_field(i), this%scaling_of_field(i)
                        !
                    END DO
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
1100    FORMAT(/, 4('%'), " BOUNDARY ", 66('%'))
        !
1101    FORMAT(/, " boundary label             = ", A20, /, &
                " boundary mode              = ", A20)
        !
1102    FORMAT(/, " boundary is built as a type-", I1, " function of a smooth density")
        !
1103    FORMAT(/, " using the Fattebert-Gygi function:", /, &
                " rhozero                    = ", F14.7, /, &
                " 2*beta                     = ", F14.7)
        !
1104    FORMAT(/, " using the optimal SCCS function:", /, &
                " rhomax                     = ", F14.7, /, &
                " rhomin                     = ", F14.7)
        !
1105    FORMAT(" log(rhomax/rhomin)         = ", F14.7)
        !
1106    FORMAT(/, " using the modified SCCS function:", /, &
                " rhomax                     = ", F14.7, /, &
                " rhomin                     = ", F14.7)
        !
1107    FORMAT(/, " adding fictitious core-electrons")
        !
1108    FORMAT(/, " boundary is built from soft-spheres centered on ionic positions:", /, &
                " solvent-dependent scaling  = ", F14.7, /, &
                " softness parameter         = ", F14.7)
        !
1109    FORMAT(/, " boundary is built as an analytic function centered on system position:", /, &
                " center of the boundary     = ", 3F14.7, /, &
                " distance from the center   = ", F14.7, /, &
                " spread of the interface    = ", F14.7, /, &
                " dimensionality             = ", I14, /, &
                " axis                       = ", I14)
        !
1110    FORMAT(/, " volume of the QM region    = ", F14.7)
        !
1111    FORMAT(/, " surface of the QM region   = ", F14.7)
        !
1112    FORMAT(/, " using solvent-aware boundary:", /, &
                " filling threshold          = ", F14.7, /, &
                " filling spread             = ", F14.7, /, &
                " solvent radius x rad scale = ", F14.7, /, &
                " spread of solvent probe    = ", F14.7)
        !
1113    FORMAT(/, "                solvation                scaling of", /, &
                "   i | label |     radius | field flux |      field", /, &
                1X, 50('-'))
        !
1114    FORMAT(1X, I3, " | ", A5, 3(" | ", F10.4))
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE print_environ_boundary
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE class_boundary
!----------------------------------------------------------------------------------------
