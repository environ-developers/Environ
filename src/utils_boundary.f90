! Copyright (C) 2018 ENVIRON (www.quantum-environment.org)
!
!    This file is part of Environ version 1.1
!
!    Environ 1.1 is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 2 of the License, or
!    (at your option) any later version.
!
!    Environ 1.1 is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more detail, either the file
!    `License' in the root directory of the present distribution, or
!    online at <http://www.gnu.org/licenses/>.
!
!----------------------------------------------------------------------------------------
!
! Authors: Oliviero Andreussi (Department of Physics, UNT)
!          Francesco Nattino  (THEOS and NCCR-MARVEL, EPFL)
!          Nicola Marzari     (THEOS and NCCR-MARVEL, EPFL)
!
!----------------------------------------------------------------------------------------
!>
!! Module containing the main routines to handle environ_boundary
!! derived data types.
!!
!! Environ_boundary contains all the specifications and the details of
!! the smooth interface between the QM and the continuum regions of the
!! simulation cell. The main interface function is stored in the %scaled
!! component, the type also stores boundary real-space derivatives (gradient,
!! laplacian, dsurface, hessian) and other quantities needed by Environ
!! modules.
!!
!----------------------------------------------------------------------------------------
MODULE utils_boundary
    !------------------------------------------------------------------------------------
    !
    USE env_base_io, ONLY: ionode, environ_unit, verbose, depth
    !
    USE environ_param, ONLY: DP, e2
    !
    USE types_representation, ONLY: environ_functions, environ_gradient, environ_density
    USE types_core, ONLY: fft_core, fd_core, core_container
    USE types_cell, ONLY: environ_cell
    !
    USE types_physical, ONLY: environ_boundary, environ_electrons, environ_ions, &
                              environ_system
    !
    USE base_environ, ONLY: niter
    !
    USE utils_functions, ONLY: copy_environ_functions, destroy_environ_functions
    !
    USE utils_density, ONLY: create_environ_density, init_environ_density, &
                             copy_environ_density, destroy_environ_density
    !
    USE utils_gradient, ONLY: create_environ_gradient, init_environ_gradient, &
                              copy_environ_gradient, destroy_environ_gradient, &
                              update_gradient_modulus
    !
    USE utils_hessian, ONLY: create_environ_hessian, init_environ_hessian, &
                             copy_environ_hessian, destroy_environ_hessian
    !
    USE tools_functions, ONLY: density_of_functions, gradient_of_functions
    USE tools_math, ONLY: integrate_environ_density, scalar_product_environ_density
    !
    ! USE tools_field_aware, ONLY: compute_normal_field, field_aware_density, &
    !                                 compute_dion_field_drho, compute_ion_field ! #TODO field_aware
    !
    USE generate_boundary, ONLY: boundary_of_density, boundary_of_functions, &
                                 boundary_of_system, solvent_aware_boundary, &
                                 invert_boundary ! #TODO is this for DEBUGGING?
    !
    ! USE environ_debugging
    !
    !------------------------------------------------------------------------------------
    !
    PRIVATE
    !
    PUBLIC :: create_environ_boundary, init_environ_boundary_first, &
              init_environ_boundary_second, copy_environ_boundary, &
              update_environ_boundary, destroy_environ_boundary, set_soft_spheres
    !
    !------------------------------------------------------------------------------------
CONTAINS
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE create_environ_boundary(boundary, local_label)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CHARACTER(LEN=80), INTENT(IN) :: local_label
        !
        TYPE(environ_boundary), INTENT(INOUT) :: boundary
        !
        CHARACTER(LEN=80) :: label = ' '
        !
        CHARACTER(LEN=80) :: sub_name = 'create_environ_boundary'
        !
        !--------------------------------------------------------------------------------
        !
        boundary%update_status = 0
        boundary%label = local_label
        !
        label = 'boundary_'//TRIM(ADJUSTL(local_label))
        !
        CALL create_environ_density(boundary%scaled, label)
        !
        boundary%volume = 0.D0
        boundary%surface = 0.D0
        !
        boundary%need_electrons = .FALSE.
        NULLIFY (boundary%electrons)
        boundary%need_ions = .FALSE.
        NULLIFY (boundary%ions)
        boundary%need_system = .FALSE.
        NULLIFY (boundary%system)
        !
        !--------------------------------------------------------------------------------
        ! Optional components
        !
        boundary%deriv = 0
        label = 'gradboundary_'//TRIM(ADJUSTL(local_label))
        !
        CALL create_environ_gradient(boundary%gradient, label)
        !
        label = 'laplboundary_'//TRIM(ADJUSTL(local_label))
        !
        CALL create_environ_density(boundary%laplacian, label)
        !
        label = 'dsurface_'//TRIM(ADJUSTL(local_label))
        !
        CALL create_environ_density(boundary%dsurface, label)
        !
        label = 'hessboundary_'//TRIM(ADJUSTL(local_label))
        !
        CALL create_environ_hessian(boundary%hessian, label)
        !
        !--------------------------------------------------------------------------------
        ! Components required for boundary of density
        !
        label = 'boundary_density_'//TRIM(ADJUSTL(local_label))
        !
        CALL create_environ_density(boundary%density, label)
        !
        label = 'dboundary_'//TRIM(ADJUSTL(local_label))
        !
        CALL create_environ_density(boundary%dscaled, label)
        !
        label = 'd2boundary_'//TRIM(ADJUSTL(local_label))
        !
        CALL create_environ_density(boundary%d2scaled, label)
        !
        !--------------------------------------------------------------------------------
        ! Components required for boundary of functions
        !
        IF (ALLOCATED(boundary%soft_spheres)) &
            CALL env_errore(sub_name, 'Trying to create an already allocated object', 1)
        !
        !--------------------------------------------------------------------------------
        ! Components required for solvent-aware interface
        !
        boundary%solvent_aware = .FALSE.
        label = 'local_'//TRIM(ADJUSTL(local_label))
        !
        CALL create_environ_density(boundary%local, label)
        !
        label = 'probe_'//TRIM(ADJUSTL(local_label))
        !
        CALL create_environ_density(boundary%probe, label)
        !
        label = 'filling_'//TRIM(ADJUSTL(local_label))
        !
        CALL create_environ_density(boundary%filling, label)
        !
        label = 'dfilling_'//TRIM(ADJUSTL(local_label))
        !
        CALL create_environ_density(boundary%dfilling, label)
        !
        !--------------------------------------------------------------------------------
        ! Components required for field-aware interface
        !
        boundary%field_aware = .FALSE.
        label = 'normal_field_'//TRIM(ADJUSTL(local_label))
        !
        CALL create_environ_density(boundary%normal_field, label)
        !
        IF (ALLOCATED(boundary%ion_field)) &
            CALL env_errore(sub_name, 'Trying to create an already allocated object', 1)
        !
        IF (ALLOCATED(boundary%local_spheres)) &
            CALL env_errore(sub_name, 'Trying to create an already allocated object', 1)
        !
        IF (ALLOCATED(boundary%dion_field_drho)) &
            CALL env_errore(sub_name, 'Trying to create an already allocated object', 1)
        !
        IF (ALLOCATED(boundary%partial_of_ion_field)) &
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
    SUBROUTINE init_environ_boundary_first(need_gradient, need_laplacian, &
                                           need_hessian, mode, stype, rhomax, rhomin, &
                                           tbeta, const, alpha, softness, &
                                           system_distance, system_spread, &
                                           solvent_radius, radial_scale, &
                                           radial_spread, filling_threshold, &
                                           filling_spread, field_factor, &
                                           charge_asymmetry, field_max, field_min, &
                                           electrons, ions, system, core, boundary)
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
        TYPE(core_container), TARGET, INTENT(IN) :: core
        !
        TYPE(environ_boundary), INTENT(INOUT) :: boundary
        !
        CHARACTER(LEN=80) :: sub_name = 'init_environ_boundary_first'
        !
        !--------------------------------------------------------------------------------
        !
        IF (need_hessian) THEN
            boundary%deriv = 3
        ELSE IF (need_laplacian) THEN
            boundary%deriv = 2
        ELSE IF (need_gradient) THEN
            boundary%deriv = 1
        END IF
        !
        boundary%mode = mode
        !
        boundary%need_electrons = (mode == 'electronic') .OR. &
                                  (mode == 'full') .OR. &
                                  (mode == 'fa-ionic') .OR. &
                                  (mode == 'fa-electronic')
        !
        IF (boundary%need_electrons) boundary%electrons => electrons
        !
        boundary%need_ions = (mode == 'ionic') .OR. &
                             (mode == 'full') .OR. &
                             (mode == 'fa-ionic') .OR. &
                             (mode == 'fa-electronic')
        !
        IF (boundary%need_ions) boundary%ions => ions
        !
        boundary%need_system = (mode == 'system')
        !
        IF (boundary%need_system) boundary%system => system
        !
        boundary%type_ = stype
        boundary%rhomax = rhomax
        boundary%rhomin = rhomin
        boundary%fact = LOG(rhomax / rhomin)
        boundary%rhozero = (rhomax + rhomin) * 0.5_DP
        boundary%tbeta = tbeta
        boundary%deltarho = rhomax - rhomin
        !
        IF (const == 1.D0 .AND. boundary%need_electrons .AND. stype == 2) &
            CALL env_errore(sub_name, &
                            'stype=2 boundary requires dielectric constant > 1', 1)
        !
        boundary%const = const
        !
        boundary%alpha = alpha
        boundary%softness = softness
        !
        IF (boundary%mode == 'ionic' .OR. boundary%mode == 'fa-ionic') &
            ALLOCATE (boundary%soft_spheres(boundary%ions%number))
        !
        boundary%simple%type_ = 4
        boundary%simple%pos => system%pos
        boundary%simple%volume = 1.D0
        boundary%simple%dim = system%dim
        boundary%simple%axis = system%axis
        boundary%simple%width = system_distance
        boundary%simple%spread = system_spread
        !
        boundary%solvent_aware = solvent_radius > 0.D0
        !
        IF (boundary%solvent_aware) THEN
            boundary%solvent_probe%type_ = 2
            ALLOCATE (boundary%solvent_probe%pos(3))
            boundary%solvent_probe%pos = 0.D0
            boundary%solvent_probe%volume = 1.D0
            boundary%solvent_probe%dim = 0
            boundary%solvent_probe%axis = 1
            boundary%solvent_probe%spread = radial_spread
            boundary%solvent_probe%width = solvent_radius * radial_scale
        END IF
        !
        boundary%filling_threshold = filling_threshold
        boundary%filling_spread = filling_spread
        !
        boundary%core => core
        boundary%field_aware = field_factor > 0.D0
        boundary%field_factor = field_factor
        boundary%charge_asymmetry = charge_asymmetry
        boundary%field_max = field_max
        boundary%field_min = field_min
        !
        IF (boundary%field_aware .AND. boundary%mode == 'fa-ionic') THEN
            ALLOCATE (boundary%ion_field(boundary%ions%number))
            ALLOCATE (boundary%dion_field_drho(boundary%ions%number))
            !
            ALLOCATE (boundary%partial_of_ion_field(3, boundary%ions%number, &
                                                    boundary%ions%number))
            !
            ALLOCATE (boundary%local_spheres(boundary%ions%number))
        END IF
        !
        boundary%initialized = .FALSE.
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_environ_boundary_first
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_environ_boundary_second(cell, boundary)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_cell), INTENT(IN) :: cell
        !
        TYPE(environ_boundary), INTENT(INOUT) :: boundary
        !
        INTEGER :: i
        !
        CHARACTER(LEN=80) :: sub_name = 'init_environ_boundary_second'
        !
        !--------------------------------------------------------------------------------
        !
        CALL init_environ_density(cell, boundary%scaled)
        !
        IF (boundary%mode == 'electronic' .OR. boundary%mode == 'full' .OR. &
            boundary%mode == 'fa-electronic' .OR. boundary%mode == 'fa-full') THEN
            !
            CALL init_environ_density(cell, boundary%density)
            !
            CALL init_environ_density(cell, boundary%dscaled)
            !
            CALL init_environ_density(cell, boundary%d2scaled)
            !
        END IF
        !
        IF (boundary%deriv >= 1) CALL init_environ_gradient(cell, boundary%gradient)
        !
        IF (boundary%deriv >= 2) CALL init_environ_density(cell, boundary%laplacian)
        !
        IF (boundary%deriv >= 3) CALL init_environ_density(cell, boundary%dsurface)
        !
        IF (boundary%solvent_aware) THEN
            !
            CALL init_environ_density(cell, boundary%local)
            !
            CALL init_environ_density(cell, boundary%probe)
            !
            CALL init_environ_density(cell, boundary%filling)
            !
            CALL init_environ_density(cell, boundary%dfilling)
            !
            IF (boundary%deriv >= 3) CALL init_environ_hessian(cell, boundary%hessian)
            !
        END IF
        !
        IF (boundary%field_aware) THEN
            !
            CALL env_errore(sub_name, 'field-aware not yet implimented ', 1)
            !
            IF (boundary%mode == 'fa-electronic' .OR. &
                !
                boundary%mode == 'fa-full') THEN
                !
                CALL init_environ_density(cell, boundary%normal_field)
                !
            ELSE IF (boundary%mode == 'fa-ionic') THEN
                !
                DO i = 1, boundary%ions%number
                    CALL init_environ_density(cell, boundary%dion_field_drho(i))
                END DO
                !
            ELSE
                CALL env_errore(sub_name, 'Boundary must be field-aware', 1)
            END IF
            !
        END IF
        !
        boundary%initialized = .TRUE.
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_environ_boundary_second
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE copy_environ_boundary(boriginal, bcopy)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_boundary), INTENT(IN) :: boriginal
        !
        TYPE(environ_boundary), INTENT(OUT) :: bcopy
        !
        INTEGER :: i, n, m
        !
        !--------------------------------------------------------------------------------
        !
        bcopy%electrons => boriginal%electrons
        bcopy%ions => boriginal%ions
        bcopy%system => boriginal%system
        bcopy%core => boriginal%core
        !
        bcopy%mode = boriginal%mode
        bcopy%update_status = boriginal%update_status
        bcopy%need_electrons = boriginal%need_electrons
        bcopy%need_ions = boriginal%need_ions
        bcopy%need_system = boriginal%need_system
        bcopy%deriv = boriginal%deriv
        bcopy%volume = boriginal%volume
        bcopy%surface = boriginal%surface
        bcopy%type_ = boriginal%type_
        bcopy%rhomax = boriginal%rhomax
        bcopy%rhomin = boriginal%rhomin
        bcopy%fact = boriginal%fact
        bcopy%rhozero = boriginal%rhozero
        bcopy%deltarho = boriginal%deltarho
        bcopy%tbeta = boriginal%tbeta
        bcopy%const = boriginal%const
        bcopy%alpha = boriginal%alpha
        bcopy%softness = boriginal%softness
        bcopy%solvent_aware = boriginal%solvent_aware
        bcopy%filling_threshold = boriginal%filling_threshold
        bcopy%filling_spread = boriginal%filling_spread
        bcopy%field_aware = boriginal%field_aware
        bcopy%field_factor = boriginal%field_factor
        bcopy%charge_asymmetry = boriginal%charge_asymmetry
        bcopy%field_max = boriginal%field_max
        bcopy%field_min = boriginal%field_min
        bcopy%initialized = boriginal%initialized
        !
        IF (ASSOCIATED(boriginal%scaled%cell)) &
            CALL copy_environ_density(boriginal%scaled, bcopy%scaled)
        !
        IF (ASSOCIATED(boriginal%gradient%cell)) &
            CALL copy_environ_gradient(boriginal%gradient, bcopy%gradient)
        !
        IF (ASSOCIATED(boriginal%laplacian%cell)) &
            CALL copy_environ_density(boriginal%laplacian, bcopy%laplacian)
        !
        IF (ASSOCIATED(boriginal%dsurface%cell)) &
            CALL copy_environ_density(boriginal%dsurface, bcopy%dsurface)
        !
        IF (ASSOCIATED(boriginal%hessian%cell)) &
            CALL copy_environ_hessian(boriginal%hessian, bcopy%hessian)
        !
        IF (ASSOCIATED(boriginal%density%cell)) &
            CALL copy_environ_density(boriginal%density, bcopy%density)
        !
        IF (ASSOCIATED(boriginal%dscaled%cell)) &
            CALL copy_environ_density(boriginal%dscaled, bcopy%dscaled)
        !
        IF (ASSOCIATED(boriginal%d2scaled%cell)) &
            CALL copy_environ_density(boriginal%d2scaled, bcopy%d2scaled)
        !
        CALL copy_environ_functions(boriginal%simple, bcopy%simple)
        !
        CALL copy_environ_functions(boriginal%solvent_probe, bcopy%solvent_probe)
        !
        IF (ASSOCIATED(boriginal%local%cell)) &
            CALL copy_environ_density(boriginal%local, bcopy%local)
        !
        IF (ASSOCIATED(boriginal%probe%cell)) &
            CALL copy_environ_density(boriginal%probe, bcopy%probe)
        !
        IF (ASSOCIATED(boriginal%filling%cell)) &
            CALL copy_environ_density(boriginal%filling, bcopy%filling)
        !
        IF (ASSOCIATED(boriginal%dfilling%cell)) &
            CALL copy_environ_density(boriginal%dfilling, bcopy%dfilling)
        !
        IF (ASSOCIATED(boriginal%normal_field%cell)) &
            CALL copy_environ_density(boriginal%normal_field, bcopy%normal_field)
        !
        IF (ALLOCATED(boriginal%soft_spheres)) THEN
            n = SIZE(boriginal%soft_spheres)
            !
            IF (ALLOCATED(bcopy%soft_spheres)) THEN
                m = SIZE(bcopy%soft_spheres)
                !
                CALL destroy_environ_functions(m, bcopy%soft_spheres)
                !
            END IF
            !
            ALLOCATE (bcopy%soft_spheres(n))
            !
            DO i = 1, n
                !
                CALL copy_environ_functions(boriginal%soft_spheres(i), &
                                            bcopy%soft_spheres(i))
                !
            END DO
            !
        ELSE
            IF (ALLOCATED(bcopy%soft_spheres)) DEALLOCATE (bcopy%soft_spheres)
        END IF
        !
        IF (ALLOCATED(boriginal%ion_field)) THEN
            n = SIZE(boriginal%ion_field)
            !
            IF (ALLOCATED(bcopy%ion_field)) DEALLOCATE (bcopy%ion_field)
            !
            IF (ALLOCATED(bcopy%partial_of_ion_field)) &
                DEALLOCATE (bcopy%partial_of_ion_field)
            !
            ALLOCATE (bcopy%ion_field(n))
            ALLOCATE (bcopy%partial_of_ion_field(3, n, n))
            !
            IF (ALLOCATED(bcopy%dion_field_drho)) THEN
                m = SIZE(bcopy%dion_field_drho)
                !
                DO i = 1, m
                    CALL destroy_environ_density(bcopy%dion_field_drho(i))
                END DO
                !
                DEALLOCATE (bcopy%dion_field_drho)
            END IF
            !
            ALLOCATE (bcopy%dion_field_drho(n))
            !
            DO i = 1, n
                !
                CALL copy_environ_density(boriginal%dion_field_drho(i), &
                                          bcopy%dion_field_drho(i))
                !
            END DO
            !
        ELSE
            !
            IF (ALLOCATED(bcopy%ion_field)) DEALLOCATE (bcopy%ion_field)
            !
            IF (ALLOCATED(bcopy%partial_of_ion_field)) &
                DEALLOCATE (bcopy%partial_of_ion_field)
            !
            IF (ALLOCATED(bcopy%dion_field_drho)) DEALLOCATE (bcopy%dion_field_drho)
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
    SUBROUTINE update_environ_boundary(bound)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_boundary), INTENT(INOUT) :: bound
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
        cell => bound%scaled%cell
        !
        update_anything = .FALSE.
        IF (bound%need_ions) update_anything = bound%ions%update
        !
        IF (bound%need_electrons) &
            update_anything = update_anything .OR. bound%electrons%update
        !
        IF (bound%need_system) &
            update_anything = update_anything .OR. bound%system%update
        !
        IF (.NOT. update_anything) THEN
            !
            IF (bound%update_status == 2) bound%update_status = 0
            ! nothing is under update, change update_status and exit
            !
            RETURN
            !
        END IF
        !
        SELECT CASE (bound%mode)
            !
        CASE ('full')
            !
            IF (bound%ions%update) THEN
                !
                !------------------------------------------------------------------------
                ! Compute the ionic part
                !
                CALL density_of_functions(bound%ions%number, bound%ions%core_electrons, &
                                          bound%ions%core, .TRUE.)
                !
                bound%update_status = 1 ! waiting to finish update
                !
            END IF
            !
            IF (bound%electrons%update) THEN
                !
                !------------------------------------------------------------------------
                ! Check if the ionic part has been updated
                !
                IF (bound%update_status == 0) &
                    CALL env_errore(sub_name, &
                                    'Wrong update status, possibly &
                                    &missing ionic update', 1)
                !
                bound%density%of_r = bound%electrons%density%of_r + bound%ions%core%of_r
                !
                CALL boundary_of_density(bound%density, bound)
                !
                bound%update_status = 2 ! boundary has changed and is ready
                !
            END IF
            !
        CASE ('electronic')
            !
            IF (bound%electrons%update) THEN
                !
                bound%density%of_r = bound%electrons%density%of_r
                !
                CALL boundary_of_density(bound%density, bound)
                !
                bound%update_status = 2 ! boundary has changes and is ready
                !
                ! CALL test_energy_derivatives(1, bound) ! DEBUGGING
                !
            ELSE
                !
                IF (bound%update_status == 2) bound%update_status = 0
                ! boundary has not changed
                !
                RETURN
                !
            END IF
            !
        CASE ('ionic')
            !
            IF (bound%ions%update) THEN
                !
                !------------------------------------------------------------------------
                ! Only ions are needed, fully update the boundary
                !
                CALL boundary_of_functions(bound%ions%number, bound%soft_spheres, bound)
                !
                bound%update_status = 2 ! boundary has changed and is ready
                !
            ELSE
                !
                IF (bound%update_status == 2) bound%update_status = 0
                ! boundary has not changed
                !
                RETURN
                !
            END IF
            !
            ! CASE ('fa-electronic') ! #TODO field-aware
            !     !
            !     IF (bound%ions%update) THEN
            !         bound%update_status = 1 ! waiting to finish update
            !     END IF
            !     !
            !     IF (bound%electrons%update) THEN
            !         !
            !         CALL compute_normal_field(bound%ions, bound%electrons, &
            !                                   bound%normal_field)
            !         !
            !         CALL field_aware_density(bound%electrons, bound)
            !         !
            !         CALL boundary_of_density(bound%density, bound)
            !         !
            !         ! TO DEBUG FIELD-AWARE: testing energy derivatives
            !         !
            !         ! CALL extract_boundary_data(bound)
            !         !
            !         ! CALL test_energy_derivatives(2, bound)
            !         !
            !         ! IF (niter == 1) CALL test_energy_derivatives(2, bound)
            !         !
            !         ! CALL test_normal_field_derivatives(bound)
            !         !
            !         bound%update_status = 2 ! boundary has changes and is ready
            !     END IF
            !     !
            ! CASE ('fa-ionic')
            !     !
            !     IF (bound%ions%update) THEN
            !         !
            !         CALL compute_dion_field_drho(bound%ions%number, bound%local_spheres, &
            !                                      bound%dion_field_drho, bound%core%fft &
            !                                      )
            !         !
            !         bound%update_status = 1 ! waiting to finish update
            !     END IF
            !     !
            !     IF (bound%electrons%update) THEN
            !         !
            !         CALL compute_ion_field(bound%ions%number, bound%local_spheres, &
            !                                bound%ions, bound%electrons, bound%ion_field)
            !         !
            !         CALL set_soft_spheres(bound, .TRUE.)
            !         !
            !         CALL boundary_of_functions(bound%ions%number, bound%soft_spheres, bound)
            !         !
            !         ! TO DEBUG FIELD - AWARE:testing ion_field derivatives
            !         !
            !         CALL test_ion_field_derivatives(2, bound)
            !         !
            !         CALL test_energy_derivatives(2, bound)
            !         !
            !         ! TO DEBUG FIELD - AWARE:testing energy derivatives
            !         !
            !         IF (ionode) WRITE (program_unit, '(1X,a,i14.7)') ' niter = ', niter
            !         !
            !         IF (ionode) WRITE (environ_unit, '(a,i14.7)') ' niter = ', niter
            !         !
            !         IF (niter == 32) CALL test_energy_derivatives(2, bound)
            !         !
            !         IF (niter == 32) CALL test_ion_field_derivatives(2, bound)
            !         !
            !         bound%update_status = 2 ! boundary has changes and is ready
            !         !
            !     END IF
            !
        CASE ('system')
            !
            IF (bound%system%update) THEN
                !
                !------------------------------------------------------------------------
                ! Only ions are needed, fully update the boundary
                !
                CALL boundary_of_system(bound%simple, bound)
                !
                ! TO DEBUG SOLVENT-AWARE
                ! !
                ! CALL invert_boundary(bound)
                ! !
                ! CALL test_de_dboundary(bound)
                !
                bound%update_status = 2 ! boundary has changed and is ready
                !
            ELSE
                !
                IF (bound%update_status == 2) bound%update_status = 0
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
        IF (bound%update_status == 2 .AND. bound%solvent_aware) &
            CALL solvent_aware_boundary(bound)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE update_environ_boundary
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE destroy_environ_boundary(lflag, boundary)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        LOGICAL, INTENT(IN) :: lflag
        !
        TYPE(environ_boundary), INTENT(INOUT) :: boundary
        !
        INTEGER :: i
        !
        CHARACTER(LEN=80) :: sub_name = 'destroy_environ_boundary'
        !
        !--------------------------------------------------------------------------------
        !
        IF (boundary%initialized) THEN
            !
            CALL destroy_environ_density(boundary%scaled)
            !
            IF (boundary%mode == 'electronic' .OR. boundary%mode == 'full') THEN
                !
                CALL destroy_environ_density(boundary%density)
                !
                CALL destroy_environ_density(boundary%dscaled)
                !
                CALL destroy_environ_density(boundary%d2scaled)
                !
            END IF
            !
            IF (boundary%deriv >= 1) CALL destroy_environ_gradient(boundary%gradient)
            !
            IF (boundary%deriv >= 2) CALL destroy_environ_density(boundary%laplacian)
            !
            IF (boundary%deriv >= 3) CALL destroy_environ_density(boundary%dsurface)
            !
            IF (boundary%solvent_aware) THEN
                !
                CALL destroy_environ_density(boundary%local)
                !
                CALL destroy_environ_density(boundary%probe)
                !
                CALL destroy_environ_density(boundary%filling)
                !
                CALL destroy_environ_density(boundary%dfilling)
                !
                IF (boundary%deriv >= 3) &
                    CALL destroy_environ_hessian(boundary%hessian)
                !
            END IF
            !
            ! IF (boundary%field_aware) THEN ! #TODO field-aware
            !     !
            !     IF (boundary%mode == 'fa-electronic' .OR. &
            !         boundary%mode == 'fa-full') THEN
            !         !
            !         CALL destroy_environ_density(boundary%normal_field)
            !         !
            !     ELSE IF (boundary%mode == 'fa-ionic') THEN
            !         !
            !         DO i = 1, boundary%ions%number
            !             CALL destroy_environ_density(boundary%dion_field_drho(i))
            !         END DO
            !         !
            !     END IF
            !     !
            ! END IF
            !
            boundary%initialized = .FALSE.
            !
        END IF
        !
        IF (lflag) THEN
            !
            !----------------------------------------------------------------------------
            ! These components were allocated first, destroy only if lflag = .TRUE.
            !
            IF (boundary%need_ions) THEN
                IF (boundary%mode == 'ionic' .OR. boundary%mode == 'fa-ionic') THEN
                    !
                    CALL destroy_environ_functions(boundary%ions%number, &
                                                   boundary%soft_spheres)
                    !
                    IF (boundary%field_aware .AND. boundary%mode == 'fa-ionic') THEN
                        !
                        CALL env_errore(sub_name, 'field-aware not yet implimented ', 1)
                        !
                        DEALLOCATE (boundary%ion_field)
                        DEALLOCATE (boundary%partial_of_ion_field)
                        !
                        CALL destroy_environ_functions(boundary%ions%number, &
                                                       boundary%local_spheres)
                        !
                        DEALLOCATE (boundary%dion_field_drho)
                    END IF
                    !
                END IF
                !
                IF (.NOT. ASSOCIATED(boundary%ions)) &
                    CALL env_errore(sub_name, &
                                    'Trying to destroy a non associated object', 1)
                !
                NULLIFY (boundary%ions)
            ELSE
                !
                IF (ASSOCIATED(boundary%ions)) &
                    CALL env_errore(sub_name, 'Found an unexpected associated object', 1)
                !
            END IF
            !
            IF (boundary%need_electrons) THEN
                IF (ASSOCIATED(boundary%electrons)) NULLIFY (boundary%electrons)
            END IF
            !
            IF (boundary%solvent_aware) DEALLOCATE (boundary%solvent_probe%pos)
            !
            IF (boundary%need_system) THEN
                IF (ASSOCIATED(boundary%system)) NULLIFY (boundary%system)
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
    SUBROUTINE set_soft_spheres(boundary, scale)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        LOGICAL, INTENT(IN), OPTIONAL :: scale
        !
        TYPE(environ_boundary), INTENT(INOUT) :: boundary
        !
        ! LOGICAL :: lscale1, lscale2 ! #TODO field-aware
        INTEGER :: i
        REAL(DP) :: radius !, f ! #TODO field-aware
        !
        !--------------------------------------------------------------------------------
        !
        IF (.NOT. (boundary%mode == 'ionic' .OR. boundary%mode == 'fa-ionic')) RETURN
        !
        ! #TODO field-aware
        !
        ! lscale1 = .FALSE.
        ! lscale2 = .FALSE.
        !
        ! IF (PRESENT(scale)) THEN
        !     lscale1 = scale .AND. (boundary%mode == 'fa-ionic')
        ! ELSE
        !     lscale2 = (boundary%mode == 'fa-ionic')
        ! END IF

        ! f = 1.D0
        !
        DO i = 1, boundary%ions%number ! #TODO field-aware
            !
            ! IF (lscale1) f = scaling_of_field(boundary%field_factor, &
            !                                   boundary%charge_asymmetry, &
            !                                   boundary%field_max, &
            !                                   boundary%field_min, &
            !                                   boundary%ion_field(i))
            !
            radius = boundary%ions%iontype(boundary%ions%ityp(i))%solvationrad * &
                     boundary%alpha ! * f
            !
            boundary%soft_spheres(i) = environ_functions(5, 1, 0, radius, &
                                                         boundary%softness, 1.D0, &
                                                         boundary%ions%tau(:, i))
            !
            ! IF (lscale2) boundary%local_spheres(i) = boundary%soft_spheres(i)

            ! IF (lscale1 .AND. verbose >= 1) &
            !     WRITE (environ_unit, 6100) &
            !     i, boundary%ions%iontype(boundary%ions%ityp(i))%label, &
            !     boundary%ions%iontype(boundary%ions%ityp(i))%solvationrad, &
            !     boundary%alpha, boundary%ion_field(i), f, radius
            !
! 6100        FORMAT("atom numer = ", i3, " atom label = ", a3, &
            !    " solvation radius = ", f8.4, " scaling = ", f8.4, " field flux = ", &
            !    f8.4, " scaling of field = ", f8.4, " final radius = ", f8.4)
        END DO
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE set_soft_spheres
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE utils_boundary
!----------------------------------------------------------------------------------------
