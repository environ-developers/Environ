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
MODULE class_boundary_ionic
    !------------------------------------------------------------------------------------
    !
    USE class_io, ONLY: io
    USE env_mp, ONLY: env_mp_bcast, env_mp_sum
    !
    USE environ_param, ONLY: DP, tpi
    !
    USE class_cell
    USE class_density
    USE class_function_erfc
    USE class_functions
    USE class_gradient
    USE class_hessian
    !
    USE class_boundary
    USE class_electrons
    USE class_ions
    !
    USE boundary_tools
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
    TYPE, EXTENDS(environ_boundary), PUBLIC :: environ_boundary_ionic
        !--------------------------------------------------------------------------------
        !
        TYPE(environ_ions), POINTER :: ions => NULL()
        !
        REAL(DP) :: alpha ! solvent-dependent scaling factor
        REAL(DP) :: softness ! sharpness of the interface
        !
        TYPE(environ_functions) :: soft_spheres
        !
        !--------------------------------------------------------------------------------
        ! Field aware
        !
        TYPE(environ_electrons), POINTER :: electrons => NULL()
        !
        TYPE(environ_functions), ALLOCATABLE :: unscaled_spheres
        !
        REAL(DP), ALLOCATABLE :: ion_field(:)
        TYPE(environ_density), ALLOCATABLE :: dion_field_drho(:)
        REAL(DP), ALLOCATABLE :: partial_of_ion_field(:, :, :)
        !
        !--------------------------------------------------------------------------------
    CONTAINS
        !--------------------------------------------------------------------------------
        !
        PROCEDURE :: create => create_environ_boundary
        PROCEDURE :: init => init_environ_boundary
        PROCEDURE :: update => update_environ_boundary
        PROCEDURE :: destroy => destroy_environ_boundary
        PROCEDURE :: build => boundary_of_functions
        !
        PROCEDURE :: dboundary_dions => calc_dboundary_dions
        !
        PROCEDURE :: fa_dboundary_dions => calc_field_aware_dboundary_dions
        PROCEDURE :: ion_field_partial => calc_ion_field_partial
        !
        PROCEDURE, PRIVATE :: set_soft_spheres
        !
        PROCEDURE, PRIVATE :: update_soft_spheres
        PROCEDURE, PRIVATE :: calc_ion_field
        PROCEDURE, PRIVATE :: calc_dion_field_drho
        PROCEDURE, PRIVATE :: scaling_of_field
        PROCEDURE, PRIVATE :: dscaling_of_field
        !
        PROCEDURE :: printout => print_environ_boundary
        !
        !--------------------------------------------------------------------------------
    END TYPE environ_boundary_ionic
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
        CLASS(environ_boundary_ionic), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: routine = 'create_environ_boundary'
        !
        !--------------------------------------------------------------------------------
        !
        IF (ASSOCIATED(this%ions)) CALL io%create_error(routine)
        !
        IF (ASSOCIATED(this%electrons)) CALL io%create_error(routine)
        !
        IF (ALLOCATED(this%ion_field)) CALL io%create_error(routine)
        !
        IF (ALLOCATED(this%dion_field_drho)) CALL io%create_error(routine)
        !
        IF (ALLOCATED(this%partial_of_ion_field)) CALL io%create_error(routine)
        !
        IF (ALLOCATED(this%unscaled_spheres)) CALL io%create_error(routine)
        !
        !--------------------------------------------------------------------------------
        !
        this%alpha = 0.D0
        this%softness = 0.D0
        !
        NULLIFY (this%ions)
        NULLIFY (this%electrons)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE create_environ_boundary
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_environ_boundary(this, alpha, softness, ions, electrons)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP), INTENT(IN) :: alpha
        REAL(DP), INTENT(IN) :: softness
        !
        TYPE(environ_ions), TARGET, INTENT(IN) :: ions
        TYPE(environ_electrons), OPTIONAL, TARGET, INTENT(IN) :: electrons
        !
        CLASS(environ_boundary_ionic), TARGET, INTENT(INOUT) :: this
        !
        INTEGER :: i
        !
        INTEGER, POINTER :: nat, nnr
        !
        CHARACTER(LEN=80) :: routine = 'init_environ_boundary'
        !
        !--------------------------------------------------------------------------------
        !
        CALL this%create()
        !
        this%ions => ions
        nat => this%ions%number
        !
        this%alpha = alpha
        this%softness = softness
        !
        !--------------------------------------------------------------------------------
        ! Field aware
        !
        IF (this%field_aware) THEN
            !
            IF (.NOT. PRESENT(electrons)) CALL io%error(routine, "Missing electrons", 1)
            !
            this%electrons => electrons
            !
            ALLOCATE (this%ion_field(nat))
            ALLOCATE (this%dion_field_drho(nat))
            ALLOCATE (this%partial_of_ion_field(3, nat, nat))
            !
            DO i = 1, nat
                CALL this%dion_field_drho(i)%init(this%cell)
            END DO
            !
        END IF
        !
        !--------------------------------------------------------------------------------
        !
        CALL this%set_soft_spheres()
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_environ_boundary
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE update_environ_boundary(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_boundary_ionic), INTENT(INOUT) :: this
        !
        LOGICAL :: update_anything
        !
        CHARACTER(LEN=80) :: routine = 'update_environ_boundary'
        !
        !--------------------------------------------------------------------------------
        !
        update_anything = .FALSE.
        !
        update_anything = this%ions%lupdate
        !
        IF (this%field_aware) &
            update_anything = update_anything .OR. this%electrons%lupdate
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
                CALL this%update_soft_spheres()
                !
                CALL this%build()
                !
                this%update_status = 2
            END IF
            !
        ELSE IF (this%ions%lupdate) THEN
            !
            !----------------------------------------------------------------------------
            ! Only ions are needed, fully update the boundary
            !
            CALL this%build()
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
        CALL this%update_solvent_aware() ! update solvent aware if applicable
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
        CLASS(environ_boundary_ionic), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: routine = 'destroy_environ_boundary'
        !
        !--------------------------------------------------------------------------------
        !
        CALL this%pre_destroy()
        !
        CALL this%soft_spheres%destroy()
        !
        NULLIFY (this%ions)
        !
        IF (this%field_aware) THEN
            DEALLOCATE (this%ion_field)
            DEALLOCATE (this%dion_field_drho)
            DEALLOCATE (this%partial_of_ion_field)
            !
            CALL this%unscaled_spheres%destroy()
            !
            NULLIFY (this%electrons)
        END IF
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE destroy_environ_boundary
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
    SUBROUTINE boundary_of_functions(this, density)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_density), OPTIONAL, TARGET, INTENT(IN) :: density
        !
        CLASS(environ_boundary_ionic), TARGET, INTENT(INOUT) :: this
        !
        INTEGER :: i, j, count
        !
        TYPE(environ_density), ALLOCATABLE :: local(:)
        TYPE(environ_gradient), ALLOCATABLE :: gradloc(:)
        TYPE(environ_density), ALLOCATABLE :: laplloc(:)
        TYPE(environ_hessian), ALLOCATABLE :: hessloc(:)
        TYPE(environ_hessian), POINTER :: hess
        !
        REAL(DP), ALLOCATABLE :: r(:, :, :), dist(:, :)
        !
        CHARACTER(LEN=80) :: routine = 'boundary_of_functions'
        !
        !--------------------------------------------------------------------------------
        !
        ASSOCIATE (cell => this%scaled%cell, &
                   nss => this%soft_spheres%number, &
                   soft_spheres => this%soft_spheres%array, &
                   derivatives => this%cores%derivatives, &
                   ng => this%need_gradient, &
                   nl => this%need_laplacian, &
                   nh => this%need_hessian, &
                   scal => this%scaled, &
                   grad => this%gradient, &
                   lapl => this%laplacian, &
                   dsurf => this%dsurface, &
                   nat => this%ions%number)
            !
            !----------------------------------------------------------------------------
            ! Generate boundary derivatives, if needed
            !
            IF (nh) THEN
                !
                IF (this%solvent_aware) THEN
                    hess => this%hessian
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
                ALLOCATE (local(1))
                !
                !------------------------------------------------------------------------
                ! Compute soft spheres and generate boundary
                !
                scal%of_r = 1.D0
                !
                CALL local(1)%init(cell)
                !
                DO i = 1, nss
                    !
                    CALL soft_spheres(i)%density(local(1), .TRUE.)
                    !
                    scal%of_r = scal%of_r * local(1)%of_r
                END DO
                !
                CALL this%compute_boundary_derivatives_fft(scal, hess)
                !
                CALL local(1)%destroy()
                !
                DEALLOCATE(local)
                !
            CASE ('highmem')
                !
                ALLOCATE (local(nss))
                !
                IF (ng) ALLOCATE(gradloc(nss))
                !
                IF (nl .AND. .NOT. nh) ALLOCATE (laplloc(nss))
                !
                IF (nh) ALLOCATE (hessloc(nss))
                !
                !------------------------------------------------------------------------
                ! Compute soft spheres and generate boundary
                !
                scal%of_r = 1.D0
                !
                DO i = 1, nss
                    !
                    CALL local(i)%init(cell)
                    !
                    CALL soft_spheres(i)%density(local(i), .FALSE.)
                    !
                    scal%of_r = scal%of_r * local(i)%of_r
                END DO
                !
                !------------------------------------------------------------------------
                ! Compute and temporarily store soft spheres derivatives
                !
                DO i = 1, nss
                    !
                    IF (ng) CALL gradloc(i)%init(cell)
                    !
                    IF (nl .AND. .NOT. nh) CALL laplloc(i)%init(cell)
                    !
                    IF (nh) CALL hessloc(i)%init(cell)
                    !
                    IF (ng) &
                        CALL soft_spheres(i)%gradient(gradloc(i), .TRUE.)
                    !
                    IF (nl .AND. .NOT. nh) &
                        CALL soft_spheres(i)%laplacian(laplloc(i), .FALSE.)
                    !
                    IF (nh) &
                        CALL soft_spheres(i)%hessian(hessloc(i), .FALSE.)
                    !
                END DO
                !
                IF (ng .AND. .NOT. nh) &
                    CALL gradient_of_boundary(nss, local, gradloc, grad)
                !
                IF (nl .AND. .NOT. nh) &
                    CALL laplacian_of_boundary(nss, local, gradloc, laplloc, lapl)
                !
                IF (nh) THEN
                    !
                    CALL dsurface_of_boundary(nss, local, gradloc, hessloc, grad, hess, dsurf)
                    !
                    IF (nl) lapl%of_r = hess%trace()
                    !
                END IF
                !
                DO i = 1, nss
                    !
                    CALL local(i)%destroy()
                    !
                    IF (ng) CALL gradloc(i)%destroy()
                    !
                    IF (nl .AND. .NOT. nh) CALL laplloc(i)%destroy()
                    !
                    IF (nh) CALL hessloc(i)%destroy()
                    !
                END DO
                !
                DEALLOCATE(local)
                !
                IF (ng) DEALLOCATE(gradloc)
                !
                IF (nl .AND. .NOT. nh) DEALLOCATE (laplloc)
                !
                IF (nh) DEALLOCATE (hessloc)
                !
            CASE ('lowmem')
                !
                ALLOCATE (local(1))
                !
                IF (ng) ALLOCATE(gradloc(1))
                !
                IF (nl .AND. .NOT. nh) ALLOCATE (laplloc(1))
                !
                IF (nh) ALLOCATE (hessloc(1))
                !
                !------------------------------------------------------------------------
                ! Compute soft spheres and generate boundary
                !
                scal%of_r = 1.D0
                !
                CALL local(1)%init(cell)
                !
                DO i = 1, nss
                    !
                    CALL soft_spheres(i)%density(local(1), .TRUE.)
                    !
                    scal%of_r = scal%of_r * local(1)%of_r
                END DO
                !
                !------------------------------------------------------------------------
                ! Compute soft spheres derivatives
                !
                IF (ng) THEN
                    !
                    CALL gradloc(1)%init(cell)
                    !
                    grad%of_r = 0.D0
                    !
                    DO i = 1, nss
                        !
                        !----------------------------------------------------------------
                        ! Re-compute local
                        !
                        CALL soft_spheres(i)%density(local(1), .TRUE.)
                        !
                        CALL soft_spheres(i)%gradient(gradloc(1), .TRUE.)
                        !
                        CALL gradient_of_boundary(local(1), gradloc(1), scal, grad)
                        !
                    END DO
                    !
                ENDIF
                !
                IF (nl .or. nh) THEN
                    !
                    IF (nl .AND. .NOT. nh) CALL laplloc(1)%init(cell)
                    !
                    IF (nh) CALL hessloc(1)%init(cell)
                    !
                    DO i = 1, nss
                        !
                        !----------------------------------------------------------------
                        ! Re-compute local and gradloc
                        !
                        CALL soft_spheres(i)%density(local(1), .TRUE.)
                        !
                        CALL soft_spheres(i)%gradient(gradloc(1), .TRUE.)
                        !
                        IF (nl .AND. .NOT. nh) THEN
                            !
                            CALL soft_spheres(i)%laplacian(laplloc(1), .TRUE.)
                            !
                            CALL laplacian_of_boundary(local(1), gradloc(1), &
                                                       laplloc(1), scal, grad, lapl)
                            !
                        ENDIF
                        !
                        IF (nh) THEN
                            !
                            CALL soft_spheres(i)%hessian(hessloc(1), .TRUE.)
                            !
                            CALL dsurface_of_boundary(local(1), gradloc(1), hessloc(1), &
                                                      grad, hess, scal, dsurf)
                            !
                            IF (nl) lapl%of_r = hess%trace()
                            !
                        ENDIF
                        !
                    ENDDO
                    !
                ENDIF
                !
                CALL local(1)%destroy()
                !
                IF (ng) CALL gradloc(1)%destroy()
                !
                IF (nl .AND. .NOT. nh) CALL laplloc(1)%destroy()
                !
                IF (nh) CALL hessloc(1)%destroy()
                !
                DEALLOCATE(local)
                !
                IF (ng) DEALLOCATE(gradloc)
                !
                IF (nl .AND. .NOT. nh) DEALLOCATE (laplloc)
                !
                IF (nh) DEALLOCATE (hessloc)
                !
            CASE DEFAULT
                CALL io%error(routine, "Unexpected derivatives method", 1)
                !
            END SELECT
            !
            !----------------------------------------------------------------------------
            ! Final updates
            !
            scal%of_r = 1.D0 - scal%of_r
            this%volume = scal%integrate()
            !
            IF (ng) THEN
                grad%of_r = -grad%of_r
                !
                CALL grad%update_modulus()
                !
                this%surface = grad%modulus%integrate()
                !
                IF (nl) lapl%of_r = -lapl%of_r
                !
                IF (nh) THEN
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
    !------------------------------------------------------------------------------------
    !
    !                                  GENERAL METHODS
    !
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_dboundary_dions(this, index, partial)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_boundary_ionic), TARGET, INTENT(IN) :: this
        INTEGER, INTENT(IN) :: index
        !
        TYPE(environ_gradient), INTENT(INOUT) :: partial
        !
        INTEGER, POINTER :: nat
        !
        INTEGER :: i, j
        TYPE(environ_density) :: denlocal
        !
        CHARACTER(LEN=80) :: routine = 'calc_dboundary_dions'
        !
        !--------------------------------------------------------------------------------
        !
        nat => this%ions%number
        !
        IF (index > nat) CALL io%error(routine, "Index greater than number of ions", 1)
        !
        IF (index <= 0) CALL io%error(routine, "Index of ion is zero or lower", 1)
        !
        IF (this%soft_spheres%number == 0) &
            CALL io%error(routine, "Missing details of ionic boundary", 1)
        !
        CALL this%soft_spheres%array(index)%gradient(partial, .TRUE.)
        !
        CALL denlocal%init(partial%cell)
        !
        DO i = 1, nat
            !
            IF (i == index) CYCLE
            !
            CALL this%soft_spheres%array(i)%density(denlocal, .TRUE.)
            !
            DO j = 1, 3
                partial%of_r(j, :) = partial%of_r(j, :) * denlocal%of_r
            END DO
            !
        END DO
        !
        CALL denlocal%destroy()
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_dboundary_dions
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
        CLASS(environ_boundary_ionic), INTENT(INOUT) :: this
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
        CHARACTER(LEN=80) :: routine = 'calc_field_aware_dboundary_dions'
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
        END ASSOCIATE
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_field_aware_dboundary_dions
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
        CLASS(environ_boundary_ionic), INTENT(INOUT) :: this
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
        CHARACTER(LEN=80) :: routine = 'calc_ion_field_partial'
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
        CLASS(environ_boundary_ionic), INTENT(INOUT) :: this
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
        IF (this%field_aware) ALLOCATE (this%unscaled_spheres, source=this%soft_spheres)
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
        CLASS(environ_boundary_ionic), INTENT(INOUT) :: this
        !
        INTEGER :: i
        REAL(DP) :: field_scale
        !
        CHARACTER(LEN=80) :: routine = 'update_soft_spheres'
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
                IF (this%field_aware) THEN
                    field_scale = this%scaling_of_field(i)
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
        CLASS(environ_boundary_ionic), INTENT(INOUT) :: this
        !
        INTEGER :: i, j
        !
        TYPE(environ_density), ALLOCATABLE :: local(:)
        !
        TYPE(environ_density) :: aux, prod
        TYPE(environ_gradient) :: auxg, field
        !
        CHARACTER(LEN=80) :: routine = 'calc_ion_field'
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
    !! Computes the functional derivative of the flux due to the ions w.r.t the
    !! electronic density
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_dion_field_drho(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_boundary_ionic), INTENT(INOUT) :: this
        !
        INTEGER :: i, j, k
        !
        TYPE(environ_density) :: prod
        TYPE(environ_gradient) :: auxg
        !
        TYPE(environ_density), ALLOCATABLE :: local(:)
        !
        CHARACTER(LEN=80) :: routine = 'calc_dion_field_drho'
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
    !! Returns field-aware scaling function with given ion_field and field aware
    !! boundary parameters
    !!
    !------------------------------------------------------------------------------------
    FUNCTION scaling_of_field(this, i) RESULT(scaling)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_boundary_ionic), INTENT(IN) :: this
        INTEGER, INTENT(IN) :: i
        !
        REAL(DP) :: scaling, multiplier, arg, diff
        !
        CHARACTER(LEN=80) :: routine = 'scaling_of_field'
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
        CLASS(environ_boundary_ionic), INTENT(IN) :: this
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
        CLASS(environ_boundary_ionic), INTENT(IN) :: this
        INTEGER, OPTIONAL, INTENT(IN) :: verbose, debug_verbose, unit
        !
        INTEGER :: base_verbose, local_verbose, passed_verbose, local_unit, i
        !
        CHARACTER(LEN=80) :: routine = 'print_environ_boundary'
        !
        !--------------------------------------------------------------------------------
        !
        IF (.NOT. PRESENT(debug_verbose) .AND. io%verbosity <= 0) RETURN
        !
        CALL this%pre_printout(verbose, debug_verbose, unit)
        !
        CALL this%print_setup(base_verbose, local_verbose, passed_verbose, local_unit, &
                              verbose, debug_verbose, unit)
        !
        !--------------------------------------------------------------------------------
        !
        IF (local_verbose >= 1) THEN
            !
            IF (io%lnode) WRITE (local_unit, 1100) this%alpha, this%softness
            !
            IF (local_verbose >= 3) &
                CALL this%soft_spheres%printout(passed_verbose, debug_verbose, local_unit)
            !
            IF (this%field_aware) THEN
                !
                IF (io%lnode .AND. local_verbose >= 1) THEN
                    WRITE (local_unit, 1101)
                    !
                    DO i = 1, this%ions%number
                        !
                        WRITE (local_unit, 1102) i, &
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
        END IF
        !
        FLUSH (local_unit)
        !
        !--------------------------------------------------------------------------------
        !
1100    FORMAT(/, " boundary is built from soft-spheres centered on ionic positions:", /, &
                " solvent-dependent scaling  = ", F14.7, /, &
                " softness parameter         = ", F14.7)
        !
1101    FORMAT(/, "                solvation                scaling of", /, &
                "   i | label |     radius | field flux |      field", /, &
                1X, 50('-'))
        !
1102    FORMAT(1X, I3, " | ", A5, 3(" | ", F10.4))
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE print_environ_boundary
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE class_boundary_ionic
!----------------------------------------------------------------------------------------
