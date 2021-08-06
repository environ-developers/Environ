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
!          Ismaila Dabo       (DMSE, Penn State)
!          Nicola Marzari     (THEOS and NCCR-MARVEL, EPFL)
!          Edan Bainglass     (Department of Physics, UNT)
!
!----------------------------------------------------------------------------------------
!>
!! Module containing the main routines to handle environ_dielectric
!! derived data types.
!!
!! Environ_dielectric is the type to store the details of the dielectric
!! embedding. It contains the specifics of externally-defined dielectric
!! regions and it links the boundary details. Starting from these quantities,
!! It builds the dielectric function in space (stored in %epsilon component)
!! and the factors derived from it that are required by the generalized
!! Poisson solver (gradient of the logarithm, sqrt factor needed by
!! preconditioned conjugate gradient, etc.).
!!
!----------------------------------------------------------------------------------------
MODULE class_dielectric
    !------------------------------------------------------------------------------------
    !
    USE env_base_io, ONLY: ionode, environ_unit, verbose, depth
    !
    USE environ_param, ONLY: DP, e2, fpi
    !
    USE class_cell
    USE class_density
    USE class_functions
    USE class_gradient
    !
    USE class_boundary
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
    TYPE, PUBLIC :: environ_dielectric
        !--------------------------------------------------------------------------------
        !
        LOGICAL :: lupdate = .FALSE.
        LOGICAL :: initialized = .FALSE.
        !
        !--------------------------------------------------------------------------------
        ! Basic properties of the dielectric space from input
        !
        INTEGER :: nregions
        TYPE(environ_function), ALLOCATABLE :: regions(:)
        !
        REAL(DP) :: constant
        TYPE(environ_density) :: background
        TYPE(environ_gradient) :: gradbackground
        TYPE(environ_density) :: laplbackground
        !
        !--------------------------------------------------------------------------------
        !
        TYPE(environ_boundary), POINTER :: boundary
        ! boundary is the pointer to the object controlling the interface
        ! between the QM and the continuum region
        !
        !--------------------------------------------------------------------------------
        ! The dielectric function over space is built from the boundary of the
        ! continuum environment and the basic dielectric properties of space
        !
        TYPE(environ_density) :: epsilon
        !
        TYPE(environ_density) :: depsilon
        ! this is needed in the extra term of kohn-sham/forces
        !
        !--------------------------------------------------------------------------------
        ! Quantities related to the dielectric permittivity and
        ! they may be needed by the different solvers
        !
        LOGICAL :: need_gradient = .FALSE.
        TYPE(environ_gradient) :: gradient
        !
        LOGICAL :: need_factsqrt = .FALSE.
        TYPE(environ_density) :: factsqrt
        !
        LOGICAL :: need_gradlog = .FALSE.
        TYPE(environ_gradient) :: gradlog
        !
        !--------------------------------------------------------------------------------
        ! Dielectric polarization charges and individual components
        !
        TYPE(environ_density) :: density
        LOGICAL :: need_auxiliary = .FALSE.
        TYPE(environ_density) :: iterative
        REAL(DP) :: charge = 0.0_DP
        !
        !--------------------------------------------------------------------------------
    CONTAINS
        !--------------------------------------------------------------------------------
        !
        PROCEDURE :: create => create_environ_dielectric
        PROCEDURE :: init_first => init_environ_dielectric_first
        PROCEDURE :: init_second => init_environ_dielectric_second
        PROCEDURE :: update => update_environ_dielectric
        PROCEDURE :: destroy => destroy_environ_dielectric
        !
        PROCEDURE :: set_regions => set_dielectric_regions
        !
        PROCEDURE :: of_boundary => dielectric_of_boundary
        PROCEDURE :: of_potential => dielectric_of_potential
        PROCEDURE :: de_dboundary => calc_dedielectric_dboundary
        PROCEDURE :: dv_dboundary => calc_dvdielectric_dboundary
        !
        PROCEDURE, PRIVATE :: update_background => update_dielectric_background
        !
        PROCEDURE :: printout => print_environ_dielectric
        !
        !--------------------------------------------------------------------------------
    END TYPE environ_dielectric
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
    SUBROUTINE create_environ_dielectric(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_dielectric), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: label = ' '
        !
        CHARACTER(LEN=80) :: sub_name = 'create_environ_dielectric'
        !
        this%constant = 1.0_DP
        !
        !--------------------------------------------------------------------------------
        !
        IF (ALLOCATED(this%regions)) &
            CALL env_errore(sub_name, 'Trying to create an already allocated object', 1)
        !
        label = 'background'
        !
        CALL this%background%create(label)
        !
        label = 'gradbackground'
        !
        CALL this%gradbackground%create(label)
        !
        label = 'laplbackground'
        !
        CALL this%laplbackground%create(label)
        !
        label = 'epsilon'
        !
        CALL this%epsilon%create(label)
        !
        label = 'depsilon'
        !
        CALL this%depsilon%create(label)
        !
        NULLIFY (this%boundary)
        !
        label = 'epsilon_gradlog'
        !
        CALL this%gradlog%create(label)
        !
        this%need_gradient = .FALSE.
        label = 'epsilon_gradient'
        !
        CALL this%gradient%create(label)
        !
        this%need_factsqrt = .FALSE.
        label = 'epsilon_factsqrt'
        !
        CALL this%factsqrt%create(label)
        !
        label = 'polarization_density'
        !
        CALL this%density%create(label)
        !
        this%need_auxiliary = .FALSE.
        label = 'iterative'
        !
        CALL this%iterative%create(label)
        !
        this%charge = 0.D0
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE create_environ_dielectric
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_environ_dielectric_first(this, constant, boundary, need_gradient, &
                                             need_factsqrt, need_auxiliary)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        LOGICAL, INTENT(IN) :: need_gradient, need_factsqrt, need_auxiliary
        TYPE(environ_boundary), TARGET, INTENT(IN) :: boundary
        !
        CLASS(environ_dielectric), INTENT(INOUT) :: this
        !
        REAL(DP) :: constant
        !
        !--------------------------------------------------------------------------------
        !
        this%constant = constant
        !
        this%boundary => boundary
        !
        this%need_gradient = need_gradient
        this%need_factsqrt = need_factsqrt
        !
        this%need_auxiliary = need_auxiliary
        !
        this%initialized = .FALSE.
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_environ_dielectric_first
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_environ_dielectric_second(this, cell)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_dielectric), INTENT(INOUT) :: this
        TYPE(environ_cell), INTENT(INOUT) :: cell
        !
        INTEGER :: i
        !
        CHARACTER(LEN=80) :: sub_name = 'init_environ_dielectric_second'
        !
        !--------------------------------------------------------------------------------
        !
        IF (this%nregions > 0) THEN
            !
            DO i = 1, this%nregions
                this%regions(i)%pos = this%regions(i)%pos
            END DO
            !
        END IF
        !
        CALL this%background%init(cell)
        !
        this%background%of_r(:) = this%constant
        !
        IF (this%nregions > 0) THEN
            !
            CALL this%gradbackground%init(cell)
            !
            IF (this%need_factsqrt) CALL this%laplbackground%init(cell)
            !
        END IF
        !
        CALL this%epsilon%init(cell)
        !
        CALL this%depsilon%init(cell)
        !
        CALL this%gradlog%init(cell)
        !
        IF (this%need_gradient) CALL this%gradient%init(cell)
        !
        IF (this%need_factsqrt) CALL this%factsqrt%init(cell)
        !
        CALL this%density%init(cell)
        !
        IF (this%need_auxiliary) CALL this%iterative%init(cell)
        !
        this%initialized = .TRUE.
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_environ_dielectric_second
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE update_environ_dielectric(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_dielectric), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: sub_name = 'update_environ_dielectric'
        !
        !--------------------------------------------------------------------------------
        !
        CALL env_start_clock(sub_name)
        !
        IF (this%epsilon%cell%lupdate) THEN
            !
            !----------------------------------------------------------------------------
            ! Cells has changed, may need to update the background
            !
            IF (this%nregions > 0) THEN
                !
                CALL this%update_background()
                ! recompute background dielectric and its derivative
                !
                this%lupdate = .TRUE.
                ! background has changed, need to update the dielectric when ready
                !
            END IF
            !
        END IF
        !
        !--------------------------------------------------------------------------------
        ! Check if the boundary is under update (status = 1) or has been
        ! fully updated (status = 2)
        !
        IF (this%boundary%update_status > 0) this%lupdate = .TRUE.
        !
        IF (this%lupdate) THEN
            !
            !----------------------------------------------------------------------------
            ! Update the dielectric in space and its derivatives if
            ! the boundary is ready
            !
            IF (this%boundary%update_status == 2) THEN
                !
                CALL this%of_boundary()
                !
                this%lupdate = .FALSE.
            END IF
            !
        END IF
        !
        CALL env_stop_clock(sub_name)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE update_environ_dielectric
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE destroy_environ_dielectric(this, lflag)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        LOGICAL, INTENT(IN) :: lflag
        !
        CLASS(environ_dielectric), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: sub_name = 'destroy_environ_dielectric'
        !
        !--------------------------------------------------------------------------------
        !
        IF (lflag) THEN
            !
            IF (this%nregions > 0) THEN
                CALL destroy_environ_functions(this%regions, this%nregions)
            ELSE
                !
                IF (ALLOCATED(this%regions)) &
                    CALL env_errore(sub_name, 'Found unexpected allocated object', 1)
                !
            END IF
            !
            IF (.NOT. ASSOCIATED(this%boundary)) &
                CALL env_errore(sub_name, 'Trying to destroy a non associated object', 1)
            !
            NULLIFY (this%boundary)
        END IF
        !
        IF (this%initialized) THEN
            !
            CALL this%background%destroy()
            !
            IF (this%nregions > 0) THEN
                !
                CALL this%gradbackground%destroy()
                !
                IF (this%need_factsqrt) CALL this%laplbackground%destroy()
                !
            END IF
            !
            CALL this%epsilon%destroy()
            !
            CALL this%depsilon%destroy()
            !
            CALL this%gradlog%destroy()
            !
            IF (this%need_gradient) CALL this%gradient%destroy()
            !
            IF (this%need_factsqrt) CALL this%factsqrt%destroy()
            !
            CALL this%density%destroy()
            !
            IF (this%need_auxiliary) CALL this%iterative%destroy()
            !
            this%initialized = .FALSE.
            !
        END IF
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE destroy_environ_dielectric
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
    SUBROUTINE set_dielectric_regions(this, nregions, dims, axes, pos, widths, &
                                      spreads, eps)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: nregions
        INTEGER, DIMENSION(nregions), INTENT(IN) :: dims, axes
        REAL(DP), DIMENSION(nregions), INTENT(IN) :: widths, spreads, eps
        REAL(DP), INTENT(IN) :: pos(3, nregions)
        !
        CLASS(environ_dielectric), INTENT(INOUT) :: this
        !
        !--------------------------------------------------------------------------------
        !
        this%nregions = nregions
        !
        IF (this%nregions > 0) THEN
            !
            CALL init_environ_functions(this%regions, nregions, 4, &
                                        axes, dims, widths, spreads, eps, pos)
            !
        END IF
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE set_dielectric_regions
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE dielectric_of_boundary(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_dielectric), TARGET, INTENT(INOUT) :: this
        !
        INTEGER :: ipol
        INTEGER, POINTER :: nnr
        REAL(DP), DIMENSION(:), POINTER :: factsqrteps, eps, deps, const
        REAL(DP), DIMENSION(:), POINTER :: scaled, gradscaledmod, laplscaled
        REAL(DP), DIMENSION(:, :), POINTER :: gradeps, gradlogeps, gradscaled
        !
        REAL(DP), DIMENSION(:), POINTER :: laplback, gradbackmod
        REAL(DP), POINTER :: gradback(:, :)
        !
        REAL(DP), ALLOCATABLE :: dlogeps(:)
        REAL(DP), ALLOCATABLE :: d2eps(:)
        !
        REAL(DP), ALLOCATABLE :: deps_dback(:)
        REAL(DP), ALLOCATABLE :: dlogeps_dback(:)
        REAL(DP), ALLOCATABLE :: d2eps_dback2(:)
        REAL(DP), ALLOCATABLE :: d2eps_dbackdbound(:)
        REAL(DP), ALLOCATABLE :: gradepsmod2(:)
        !
        CHARACTER(LEN=80) :: sub_name = 'dielectric_of_boundary'
        !
        !--------------------------------------------------------------------------------
        !
        nnr => this%epsilon%cell%nnr
        eps => this%epsilon%of_r
        deps => this%depsilon%of_r
        const => this%background%of_r
        scaled => this%boundary%scaled%of_r
        !
        IF (.NOT. ALLOCATED(this%boundary%gradient%of_r)) &
            CALL env_errore(sub_name, 'Missing required gradient of boundary', 1)
        !
        gradscaled => this%boundary%gradient%of_r
        gradlogeps => this%gradlog%of_r
        ALLOCATE (dlogeps(nnr))
        !
        IF (this%need_gradient) THEN
            !
            IF (.NOT. ALLOCATED(this%boundary%gradient%of_r)) &
                CALL env_errore(sub_name, 'Missing required gradient of boundary', 1)
            !
            gradscaled => this%boundary%gradient%of_r
            gradeps => this%gradient%of_r
        END IF
        !
        IF (this%need_factsqrt) THEN
            !
            IF (.NOT. ALLOCATED(this%boundary%gradient%of_r)) &
                CALL env_errore(sub_name, 'Missing required gradient of boundary', 1)
            !
            gradscaledmod => this%boundary%gradient%modulus%of_r
            !
            IF (.NOT. ALLOCATED(this%boundary%laplacian%of_r)) &
                CALL env_errore(sub_name, 'Missing required laplacian of boundary', 1)
            !
            laplscaled => this%boundary%laplacian%of_r
            factsqrteps => this%factsqrt%of_r
            ALLOCATE (d2eps(nnr))
        END IF
        !
        IF (this%nregions > 0) THEN
            gradback => this%gradbackground%of_r
            ALLOCATE (deps_dback(nnr))
            ALLOCATE (dlogeps_dback(nnr))
            !
            IF (this%need_factsqrt) THEN
                laplback => this%laplbackground%of_r
                gradbackmod => this%gradbackground%modulus%of_r
                ALLOCATE (gradepsmod2(nnr))
                ALLOCATE (d2eps_dbackdbound(nnr))
                ALLOCATE (d2eps_dback2(nnr))
            END IF
            !
        END IF
        !
        !--------------------------------------------------------------------------------
        ! Compute epsilon(r) and its derivative wrt boundary
        !
        SELECT CASE (this%boundary%b_type)
            !
        CASE (0, 2)
            eps = 1.D0 + (const - 1.D0) * (1.D0 - scaled)
            deps = (1.D0 - const)
            dlogeps = deps / eps
            !
            IF (this%need_factsqrt) d2eps = 0.D0
            !
            IF (this%nregions > 0) THEN
                deps_dback = 1.D0 - scaled
                dlogeps_dback = deps_dback / eps
                !
                IF (this%need_factsqrt) THEN
                    d2eps_dback2 = 0.D0
                    d2eps_dbackdbound = -1.D0
                END IF
                !
            END IF
            !
        CASE (1)
            eps = EXP(LOG(const) * (1.D0 - scaled))
            deps = -eps * LOG(const)
            dlogeps = -LOG(const)
            !
            IF (this%need_factsqrt) d2eps = eps * LOG(const)**2
            !
            IF (this%nregions > 0) THEN
                deps_dback = eps * (1.D0 - scaled) / const
                dlogeps_dback = (1.D0 - scaled) / const
                !
                IF (this%need_factsqrt) THEN
                    d2eps_dback2 = -deps_dback * scaled / const
                    !
                    d2eps_dbackdbound = eps / const * &
                                        (1.D0 - (1.D0 - scaled) * LOG(const))
                    !
                END IF
                !
            END IF
            !
        CASE DEFAULT
            CALL env_errore(sub_name, 'Unkown boundary type', 1)
            !
        END SELECT
        !
        !--------------------------------------------------------------------------------
        ! If needed, compute derived quantites
        !
        DO ipol = 1, 3
            gradlogeps(ipol, :) = dlogeps(:) * gradscaled(ipol, :)
            !
            IF (this%nregions > 0) &
                gradlogeps(ipol, :) = gradlogeps(ipol, :) + &
                                      dlogeps_dback(:) * gradback(ipol, :)
            !
        END DO
        !
        CALL this%gradlog%update_modulus()
        !
        DEALLOCATE (dlogeps)
        !
        IF (this%need_gradient) THEN
            !
            DO ipol = 1, 3
                gradeps(ipol, :) = deps(:) * gradscaled(ipol, :)
                !
                IF (this%nregions > 0) &
                    gradeps(ipol, :) = gradeps(ipol, :) + &
                                       deps_dback(:) * gradback(ipol, :)
                !
            END DO
            !
            CALL this%gradient%update_modulus()
            !
        END IF
        !
        IF (this%need_factsqrt) THEN
            !
            IF (this%nregions <= 0) THEN
                !
                factsqrteps = (d2eps - 0.5D0 * deps**2 / eps) * gradscaledmod**2 + &
                              deps * laplscaled
                !
            ELSE
                !
                CALL this%boundary%gradient%scalar_product(this%gradbackground, &
                                                           this%factsqrt)
                !
                IF (this%need_gradient) THEN
                    gradepsmod2 = this%gradient%modulus%of_r**2
                ELSE
                    gradepsmod2 = 0.D0
                    !
                    DO ipol = 1, 3
                        !
                        gradepsmod2(:) = gradepsmod2(:) + &
                                         (deps(:) * gradscaled(ipol, :) + &
                                          deps_dback(:) * gradback(ipol, :))**2
                        !
                    END DO
                    !
                END IF
                !
                factsqrteps = 2.D0 * d2eps_dbackdbound * factsqrteps + &
                              d2eps_dback2 * gradbackmod**2 + deps_dback * laplback + &
                              d2eps * gradscaledmod**2 + deps * laplscaled - &
                              0.5D0 * gradepsmod2 / eps
                !
            END IF
            !
            factsqrteps = factsqrteps * 0.5D0 / e2 / fpi
            DEALLOCATE (d2eps)
        END IF
        !
        IF (this%nregions > 0) THEN
            DEALLOCATE (deps_dback)
            DEALLOCATE (dlogeps_dback)
            !
            IF (this%need_factsqrt) THEN
                DEALLOCATE (d2eps_dback2)
                DEALLOCATE (d2eps_dbackdbound)
            END IF
            !
        END IF
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE dielectric_of_boundary
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE dielectric_of_potential(this, charges, potential)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_density), INTENT(IN) :: charges, potential
        !
        CLASS(environ_dielectric), INTENT(INOUT) :: this
        !
        TYPE(environ_cell), POINTER :: cell
        !
        TYPE(environ_gradient) :: gradient
        !
        CHARACTER(LEN=80) :: sub_name = 'dielectric_of_potential'
        !
        !--------------------------------------------------------------------------------
        !
        IF (.NOT. ASSOCIATED(potential%cell, charges%cell)) &
            CALL env_errore(sub_name, &
                            'Mismatch in domains of potential and charges', 1)
        !
        IF (.NOT. ASSOCIATED(potential%cell, this%density%cell)) &
            CALL env_errore(sub_name, &
                            'Mismatch in domains of potential and dielectric', 1)
        !
        cell => charges%cell
        !
        CALL gradient%init(cell)
        !
        CALL this%boundary%derivatives%gradient(potential, gradient)
        !
        CALL this%gradlog%scalar_product(gradient, this%density)
        !
        this%density%of_r = this%density%of_r / fpi / e2 + charges%of_r * &
                            (1.D0 - this%epsilon%of_r) / &
                            this%epsilon%of_r
        !
        CALL gradient%destroy()
        !
        this%charge = this%density%integrate()
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE dielectric_of_potential
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_dedielectric_dboundary(this, velectrostatic, de_dboundary)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_dielectric), INTENT(IN) :: this
        TYPE(environ_density), INTENT(IN) :: velectrostatic
        !
        TYPE(environ_density), INTENT(INOUT) :: de_dboundary
        !
        TYPE(environ_cell), POINTER :: cell
        TYPE(environ_gradient) :: gradient
        !
        !--------------------------------------------------------------------------------
        !
        cell => de_dboundary%cell
        !
        CALL gradient%init(cell)
        !
        CALL this%boundary%derivatives%gradient(velectrostatic, gradient)
        !
        CALL gradient%update_modulus()
        !
        de_dboundary%of_r = de_dboundary%of_r - &
                            gradient%modulus%of_r**2 * this%depsilon%of_r * &
                            0.5D0 / fpi / e2
        !
        CALL gradient%destroy()
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_dedielectric_dboundary
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_dvdielectric_dboundary(this, velectrostatic, dvelectrostatic, &
                                           dv_dboundary)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_dielectric), INTENT(IN) :: this
        TYPE(environ_density), INTENT(IN) :: velectrostatic
        TYPE(environ_density), INTENT(IN) :: dvelectrostatic
        !
        TYPE(environ_density), INTENT(INOUT) :: dv_dboundary
        !
        TYPE(environ_cell), POINTER :: cell
        TYPE(environ_density) :: aux
        TYPE(environ_gradient) :: gradient, dgradient
        !
        !--------------------------------------------------------------------------------
        !
        cell => dv_dboundary%cell
        !
        CALL gradient%init(cell)
        !
        CALL this%boundary%derivatives%gradient(velectrostatic, gradient)
        !
        CALL dgradient%init(cell)
        !
        CALL this%boundary%derivatives%gradient(dvelectrostatic, dgradient)
        !
        CALL aux%init(cell)
        !
        CALL gradient%scalar_product(dgradient, aux)
        !
        CALL gradient%destroy()
        !
        CALL dgradient%destroy()
        !
        dv_dboundary%of_r = dv_dboundary%of_r - &
                            aux%of_r * this%depsilon%of_r / (fpi * e2)
        !
        CALL aux%destroy()
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_dvdielectric_dboundary
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
    SUBROUTINE update_dielectric_background(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_dielectric), INTENT(INOUT) :: this
        !
        INTEGER :: i, ipol
        TYPE(environ_density) :: local
        TYPE(environ_gradient) :: gradlocal
        TYPE(environ_density) :: lapllocal
        TYPE(environ_cell), POINTER :: cell
        !
        !--------------------------------------------------------------------------------
        !
        this%background%of_r = this%constant
        !
        IF (this%nregions <= 0) RETURN
        !
        this%gradbackground%of_r = 0.D0
        !
        IF (this%need_factsqrt) this%laplbackground%of_r = 0.D0
        !
        cell => this%background%cell
        !
        CALL local%init(cell)
        !
        CALL gradlocal%init(cell)
        !
        IF (this%need_factsqrt) CALL lapllocal%init(cell)
        !
        DO i = 1, this%nregions
            !
            CALL this%regions(i)%density(local, .TRUE.)
            !
            CALL this%regions(i)%gradient(gradlocal, .TRUE.)
            !
            !----------------------------------------------------------------------------
            ! Update background and derivatives in reverse order
            !
            IF (this%need_factsqrt) THEN
                !
                CALL this%gradbackground%scalar_product(gradlocal, lapllocal)
                !
                this%laplbackground%of_r = &
                    this%laplbackground%of_r(:) * &
                    (1.D0 - local%of_r(:) / this%regions(i)%volume) - &
                    2.D0 * lapllocal%of_r / this%regions(i)%volume
                !
                CALL this%regions(i)%laplacian(lapllocal, .TRUE.)
                !
                this%laplbackground%of_r(:) = &
                    this%laplbackground%of_r(:) + &
                    lapllocal%of_r(:) * (1.D0 - this%background%of_r(:) / &
                                         this%regions(i)%volume)
                !
            END IF
            !
            DO ipol = 1, 3
                !
                this%gradbackground%of_r(ipol, :) = &
                    this%gradbackground%of_r(ipol, :) * &
                    (1.D0 - local%of_r(:) / this%regions(i)%volume) + &
                    gradlocal%of_r(ipol, :) * (1.D0 - this%background%of_r(:) / &
                                               this%regions(i)%volume)
                !
            END DO
            !
            this%background%of_r(:) = &
                this%background%of_r(:) + &
                local%of_r(:) * (1.D0 - this%background%of_r(:) / &
                                 this%regions(i)%volume)
            !
        END DO
        !
        CALL local%destroy()
        !
        CALL gradlocal%destroy()
        !
        IF (this%need_factsqrt) CALL lapllocal%destroy()
        !
        IF (verbose >= 3) CALL this%background%printout()
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE update_dielectric_background
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
    SUBROUTINE print_environ_dielectric(this, local_verbose, local_depth)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_dielectric), INTENT(IN) :: this
        INTEGER, INTENT(IN), OPTIONAL :: local_verbose
        INTEGER, INTENT(IN), OPTIONAL :: local_depth
        !
        INTEGER :: verbosity, passed_verbosity, passed_depth
        !
        CHARACTER(LEN=80) :: sub_name = 'print_environ_dielectric'
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
                    WRITE (environ_unit, 1000)
                ELSE
                    !
                    CALL env_block_divider(verbosity)
                    !
                    WRITE (environ_unit, 1001)
                END IF
                !
            END IF
            !
            IF (this%nregions == 0) THEN
                IF (ionode) WRITE (environ_unit, 1002) this%constant
            ELSE
                IF (ionode) WRITE (environ_unit, 1003) this%constant, this%nregions
                !
                CALL print_environ_functions(this%regions, this%nregions, &
                                             passed_verbosity, passed_depth)
                !
                IF (verbosity >= 4) &
                    CALL this%background%printout(passed_verbosity, passed_depth)
                !
            END IF
            !
            IF (verbosity >= 3) THEN
                !
                CALL this%density%printout(passed_verbosity, passed_depth)
                !
                CALL this%epsilon%printout(passed_verbosity, passed_depth)
                !
            END IF
            !
            IF (verbosity >= 5) &
                CALL this%depsilon%printout(passed_verbosity, passed_depth)
            !
            IF (ionode) THEN
                WRITE (environ_unit, 1004) this%need_gradient, this%need_factsqrt
                WRITE (environ_unit, 1005) this%charge
            END IF
            !
            IF (verbosity >= 5) THEN
                !
                CALL this%gradlog%printout(passed_verbosity, passed_depth)
                !
                IF (this%need_gradient) &
                    CALL this%gradient%printout(passed_verbosity, passed_depth)
                !
                IF (this%need_factsqrt) &
                    CALL this%factsqrt%printout(passed_verbosity, passed_depth)
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
1000    FORMAT(/, 4('%'), ' DIELECTRIC ', 65('%'))
1001    FORMAT(/, ' DIELECTRIC', /, ' ==========')
        !
1002    FORMAT(/, ' dielectric build on homogeneous background:', /, &
                ' environment bulk permitt.  = ', F14.7)
        !
1003    FORMAT(/, ' dielectric build in the presence of dielectric regions:', /, &
                ' environment bulk permitt.  = ', F14.7, /, &
                ' number of dielec. regions  = ', I4)
        !
1004    FORMAT(/, ' dielectric flags:', /, &
                ' need gradient              = ', L2, /, &
                ' need factor depend. sqrt   = ', L2)
        !
1005    FORMAT(/, ' total dielectric charge    = ', F14.7)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE print_environ_dielectric
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE class_dielectric
!----------------------------------------------------------------------------------------
