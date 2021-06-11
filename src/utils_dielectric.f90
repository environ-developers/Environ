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
!          Ismaila Dabo       (DMSE, Penn State)
!          Nicola Marzari     (THEOS and NCCR-MARVEL, EPFL)
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
MODULE utils_dielectric
    !------------------------------------------------------------------------------------
    !
    USE environ_param, ONLY: DP, e2, fpi
    !
    USE types_physical, ONLY: environ_dielectric, environ_boundary
    USE types_representation, ONLY: environ_density, environ_gradient
    USE types_cell, ONLY: environ_cell
    !
    USE utils_functions, ONLY: destroy_environ_functions
    !
    USE utils_density, ONLY: create_environ_density, init_environ_density, &
                             destroy_environ_density
    !
    USE utils_gradient, ONLY: create_environ_gradient, init_environ_gradient, &
                              update_gradient_modulus, destroy_environ_gradient
    !
    USE tools_math, ONLY: scalar_product_environ_gradient
    !
    USE tools_functions, ONLY: density_of_functions, gradient_of_functions, &
                               laplacian_of_functions
    !
    USE environ_output, ONLY: verbose, print_environ_density
    !
    !------------------------------------------------------------------------------------
    !
    PRIVATE
    !
    PUBLIC :: create_environ_dielectric, init_environ_dielectric_first, &
              init_environ_dielectric_second, update_environ_dielectric, &
              destroy_environ_dielectric, set_dielectric_regions
    !
    !------------------------------------------------------------------------------------
CONTAINS
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE create_environ_dielectric(dielectric)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_dielectric), INTENT(INOUT) :: dielectric
        !
        CHARACTER(LEN=80) :: label = ' '
        !
        CHARACTER(LEN=80) :: sub_name = 'create_environ_dielectric'
        !
        dielectric%constant = 1.0_DP
        !
        !--------------------------------------------------------------------------------
        !
        IF (ALLOCATED(dielectric%regions)) &
            CALL env_errore(sub_name, 'Trying to create an already allocated object', 1)
        !
        label = 'background'
        !
        CALL create_environ_density(dielectric%background, label)
        !
        label = 'gradbackground'
        !
        CALL create_environ_gradient(dielectric%gradbackground, label)
        !
        label = 'laplbackground'
        !
        CALL create_environ_density(dielectric%laplbackground, label)
        !
        label = 'epsilon'
        !
        CALL create_environ_density(dielectric%epsilon, label)
        !
        label = 'depsilon'
        !
        CALL create_environ_density(dielectric%depsilon, label)
        !
        NULLIFY (dielectric%boundary)
        !
        label = 'epsilon_gradlog'
        !
        CALL create_environ_gradient(dielectric%gradlog, label)
        !
        dielectric%need_gradient = .FALSE.
        label = 'epsilon_gradient'
        !
        CALL create_environ_gradient(dielectric%gradient, label)
        !
        dielectric%need_factsqrt = .FALSE.
        label = 'epsilon_factsqrt'
        !
        CALL create_environ_density(dielectric%factsqrt, label)
        !
        label = 'polarization_density'
        !
        CALL create_environ_density(dielectric%density, label)
        !
        dielectric%need_auxiliary = .FALSE.
        label = 'iterative'
        !
        CALL create_environ_density(dielectric%iterative, label)
        !
        dielectric%charge = 0.D0
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE create_environ_dielectric
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_environ_dielectric_first(constant, boundary, need_gradient, &
                                             need_factsqrt, need_auxiliary, dielectric)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        LOGICAL, INTENT(IN) :: need_gradient, need_factsqrt, need_auxiliary
        TYPE(environ_boundary), TARGET, INTENT(IN) :: boundary
        !
        TYPE(environ_dielectric), INTENT(INOUT) :: dielectric
        !
        REAL(DP) :: constant
        !
        !--------------------------------------------------------------------------------
        !
        dielectric%constant = constant
        !
        dielectric%boundary => boundary
        !
        dielectric%need_gradient = need_gradient
        dielectric%need_factsqrt = need_factsqrt
        !
        dielectric%need_auxiliary = need_auxiliary
        !
        dielectric%initialized = .FALSE.
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_environ_dielectric_first
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_environ_dielectric_second(cell, dielectric)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_cell), INTENT(INOUT) :: cell
        TYPE(environ_dielectric), INTENT(INOUT) :: dielectric
        !
        INTEGER :: i
        !
        CHARACTER(LEN=80) :: sub_name = 'init_environ_dielectric_second'
        !
        !--------------------------------------------------------------------------------
        !
        IF (dielectric%nregions > 0) THEN
            !
            DO i = 1, dielectric%nregions
                dielectric%regions(i)%pos = dielectric%regions(i)%pos / cell%alat
            END DO
            !
        END IF
        !
        CALL init_environ_density(cell, dielectric%background)
        !
        dielectric%background%of_r(:) = dielectric%constant
        !
        IF (dielectric%nregions > 0) THEN
            !
            CALL init_environ_gradient(cell, dielectric%gradbackground)
            !
            IF (dielectric%need_factsqrt) &
                CALL init_environ_density(cell, dielectric%laplbackground)
            !
        END IF
        !
        CALL init_environ_density(cell, dielectric%epsilon)
        !
        CALL init_environ_density(cell, dielectric%depsilon)
        !
        CALL init_environ_gradient(cell, dielectric%gradlog)
        !
        IF (dielectric%need_gradient) &
            CALL init_environ_gradient(cell, dielectric%gradient)
        !
        IF (dielectric%need_factsqrt) &
            CALL init_environ_density(cell, dielectric%factsqrt)
        !
        CALL init_environ_density(cell, dielectric%density)
        !
        IF (dielectric%need_auxiliary) &
            CALL init_environ_density(cell, dielectric%iterative)
        !
        dielectric%initialized = .TRUE.
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_environ_dielectric_second
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE update_environ_dielectric(dielectric)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_dielectric), INTENT(INOUT) :: dielectric
        !
        CHARACTER(LEN=80) :: sub_name = 'update_environ_dielectric'
        !
        !--------------------------------------------------------------------------------
        !
        CALL env_start_clock(sub_name)
        !
        IF (dielectric%epsilon%cell%update) THEN
            !
            !----------------------------------------------------------------------------
            ! Cells has changed, may need to update the background
            !
            IF (dielectric%nregions > 0) THEN
                !
                CALL update_dielectric_background(dielectric)
                ! recompute background dielectric and its derivative
                !
                dielectric%update = .TRUE.
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
        IF (dielectric%boundary%update_status > 0) dielectric%update = .TRUE.
        !
        IF (dielectric%update) THEN
            !
            !----------------------------------------------------------------------------
            ! Update the dielectric in space and its derivatives if
            ! the boundary is ready
            !
            IF (dielectric%boundary%update_status == 2) THEN
                !
                CALL dielectric_of_boundary(dielectric)
                !
                dielectric%update = .FALSE.
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
    SUBROUTINE destroy_environ_dielectric(lflag, dielectric)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        LOGICAL, INTENT(IN) :: lflag
        !
        TYPE(environ_dielectric), INTENT(INOUT) :: dielectric
        !
        CHARACTER(LEN=80) :: sub_name = 'destroy_environ_dielectric'
        !
        !--------------------------------------------------------------------------------
        !
        IF (lflag) THEN
            !
            IF (dielectric%nregions > 0) THEN
                CALL destroy_environ_functions(dielectric%nregions, dielectric%regions)
            ELSE
                !
                IF (ALLOCATED(dielectric%regions)) &
                    CALL env_errore(sub_name, 'Found unexpected allocated object', 1)
                !
            END IF
            !
            IF (.NOT. ASSOCIATED(dielectric%boundary)) &
                CALL env_errore(sub_name, 'Trying to destroy a non associated object', 1)
            !
            NULLIFY (dielectric%boundary)
        END IF
        !
        IF (dielectric%initialized) THEN
            !
            CALL destroy_environ_density(dielectric%background)
            !
            IF (dielectric%nregions > 0) THEN
                !
                CALL destroy_environ_gradient(dielectric%gradbackground)
                !
                IF (dielectric%need_factsqrt) &
                    CALL destroy_environ_density(dielectric%laplbackground)
                !
            END IF
            !
            CALL destroy_environ_density(dielectric%epsilon)
            !
            CALL destroy_environ_density(dielectric%depsilon)
            !
            CALL destroy_environ_gradient(dielectric%gradlog)
            !
            IF (dielectric%need_gradient) &
                CALL destroy_environ_gradient(dielectric%gradient)
            !
            IF (dielectric%need_factsqrt) &
                CALL destroy_environ_density(dielectric%factsqrt)
            !
            CALL destroy_environ_density(dielectric%density)
            !
            IF (dielectric%need_auxiliary) &
                CALL destroy_environ_density(dielectric%iterative)
            !
            dielectric%initialized = .FALSE.
            !
        END IF
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE destroy_environ_dielectric
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE set_dielectric_regions(nregions, epsregion_dim, epsregion_axis, &
                                      epsregion_pos, epsregion_width, &
                                      epsregion_spread, epsregion_eps, dielectric)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: nregions
        INTEGER, DIMENSION(nregions), INTENT(IN) :: epsregion_dim, epsregion_axis
        !
        REAL(DP), DIMENSION(nregions), INTENT(IN) :: epsregion_width, &
                                                     epsregion_spread, epsregion_eps
        !
        REAL(DP), INTENT(IN) :: epsregion_pos(3, nregions)
        !
        TYPE(environ_dielectric), INTENT(INOUT) :: dielectric
        !
        INTEGER :: i
        !
        !--------------------------------------------------------------------------------
        !
        dielectric%nregions = nregions
        !
        IF (dielectric%nregions > 0) THEN
            ALLOCATE (dielectric%regions(nregions))
            !
            DO i = 1, dielectric%nregions
                dielectric%regions(i)%type_ = 4
                dielectric%regions(i)%dim = epsregion_dim(i)
                dielectric%regions(i)%axis = epsregion_axis(i)
                dielectric%regions(i)%spread = epsregion_spread(i)
                dielectric%regions(i)%width = epsregion_width(i)
                dielectric%regions(i)%volume = epsregion_eps(i)
                ALLOCATE (dielectric%regions(i)%pos(3))
                dielectric%regions(i)%pos = epsregion_pos(:, i)
            END DO
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
    SUBROUTINE update_dielectric_background(dielectric)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_dielectric), INTENT(INOUT) :: dielectric
        !
        INTEGER :: i, ipol
        TYPE(environ_density) :: local
        TYPE(environ_gradient) :: gradlocal
        TYPE(environ_density) :: lapllocal
        TYPE(environ_cell), POINTER :: cell
        !
        !--------------------------------------------------------------------------------
        !
        dielectric%background%of_r = dielectric%constant
        !
        IF (dielectric%nregions <= 0) RETURN
        !
        dielectric%gradbackground%of_r = 0.D0
        !
        IF (dielectric%need_factsqrt) dielectric%laplbackground%of_r = 0.D0
        !
        cell => dielectric%background%cell
        !
        CALL init_environ_density(cell, local)
        !
        CALL init_environ_gradient(cell, gradlocal)
        !
        IF (dielectric%need_factsqrt) CALL init_environ_density(cell, lapllocal)
        !
        DO i = 1, dielectric%nregions
            !
            CALL density_of_functions(dielectric%regions(i), local, .TRUE.)
            !
            CALL gradient_of_functions(dielectric%regions(i), gradlocal, .TRUE.)
            !
            !----------------------------------------------------------------------------
            ! Update background and derivatives in reverse order
            !
            IF (dielectric%need_factsqrt) THEN
                !
                CALL scalar_product_environ_gradient(dielectric%gradbackground, &
                                                     gradlocal, lapllocal)
                !
                dielectric%laplbackground%of_r = &
                    dielectric%laplbackground%of_r(:) * &
                    (1.D0 - local%of_r(:) / dielectric%regions(i)%volume) - &
                    2.D0 * lapllocal%of_r / dielectric%regions(i)%volume
                !
                CALL laplacian_of_functions(dielectric%regions(i), lapllocal, .TRUE.)
                !
                dielectric%laplbackground%of_r(:) = &
                    dielectric%laplbackground%of_r(:) + &
                    lapllocal%of_r(:) * (1.D0 - dielectric%background%of_r(:) / &
                                         dielectric%regions(i)%volume)
                !
            END IF
            !
            DO ipol = 1, 3
                !
                dielectric%gradbackground%of_r(ipol, :) = &
                    dielectric%gradbackground%of_r(ipol, :) * &
                    (1.D0 - local%of_r(:) / dielectric%regions(i)%volume) + &
                    gradlocal%of_r(ipol, :) * (1.D0 - dielectric%background%of_r(:) / &
                                               dielectric%regions(i)%volume)
                !
            END DO
            !
            dielectric%background%of_r(:) = &
                dielectric%background%of_r(:) + &
                local%of_r(:) * (1.D0 - dielectric%background%of_r(:) / &
                                 dielectric%regions(i)%volume)
            !
        END DO
        !
        CALL destroy_environ_density(local)
        !
        CALL destroy_environ_gradient(gradlocal)
        !
        IF (dielectric%need_factsqrt) CALL destroy_environ_density(lapllocal)
        !
        IF (verbose >= 3) CALL print_environ_density(dielectric%background)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE update_dielectric_background
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE dielectric_of_boundary(dielectric)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_dielectric), TARGET, INTENT(INOUT) :: dielectric
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
        nnr => dielectric%epsilon%cell%nnr
        eps => dielectric%epsilon%of_r
        deps => dielectric%depsilon%of_r
        const => dielectric%background%of_r
        scaled => dielectric%boundary%scaled%of_r
        !
        IF (.NOT. ALLOCATED(dielectric%boundary%gradient%of_r)) &
            CALL env_errore(sub_name, 'Missing required gradient of boundary', 1)
        !
        gradscaled => dielectric%boundary%gradient%of_r
        gradlogeps => dielectric%gradlog%of_r
        ALLOCATE (dlogeps(nnr))
        !
        IF (dielectric%need_gradient) THEN
            !
            IF (.NOT. ALLOCATED(dielectric%boundary%gradient%of_r)) &
                CALL env_errore(sub_name, 'Missing required gradient of boundary', 1)
            !
            gradscaled => dielectric%boundary%gradient%of_r
            gradeps => dielectric%gradient%of_r
        END IF
        !
        IF (dielectric%need_factsqrt) THEN
            !
            IF (.NOT. ALLOCATED(dielectric%boundary%gradient%of_r)) &
                CALL env_errore(sub_name, 'Missing required gradient of boundary', 1)
            !
            gradscaledmod => dielectric%boundary%gradient%modulus%of_r
            !
            IF (.NOT. ALLOCATED(dielectric%boundary%laplacian%of_r)) &
                CALL env_errore(sub_name, 'Missing required laplacian of boundary', 1)
            !
            laplscaled => dielectric%boundary%laplacian%of_r
            factsqrteps => dielectric%factsqrt%of_r
            ALLOCATE (d2eps(nnr))
        END IF
        !
        IF (dielectric%nregions > 0) THEN
            gradback => dielectric%gradbackground%of_r
            ALLOCATE (deps_dback(nnr))
            ALLOCATE (dlogeps_dback(nnr))
            !
            IF (dielectric%need_factsqrt) THEN
                laplback => dielectric%laplbackground%of_r
                gradbackmod => dielectric%gradbackground%modulus%of_r
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
        SELECT CASE (dielectric%boundary%type_)
            !
        CASE (0, 2)
            eps = 1.D0 + (const - 1.D0) * (1.D0 - scaled)
            deps = (1.D0 - const)
            dlogeps = deps / eps
            !
            IF (dielectric%need_factsqrt) d2eps = 0.D0
            !
            IF (dielectric%nregions > 0) THEN
                deps_dback = 1.D0 - scaled
                dlogeps_dback = deps_dback / eps
                !
                IF (dielectric%need_factsqrt) THEN
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
            IF (dielectric%need_factsqrt) d2eps = eps * LOG(const)**2
            !
            IF (dielectric%nregions > 0) THEN
                deps_dback = eps * (1.D0 - scaled) / const
                dlogeps_dback = (1.D0 - scaled) / const
                !
                IF (dielectric%need_factsqrt) THEN
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
            IF (dielectric%nregions > 0) &
                gradlogeps(ipol, :) = gradlogeps(ipol, :) + &
                                      dlogeps_dback(:) * gradback(ipol, :)
            !
        END DO
        !
        CALL update_gradient_modulus(dielectric%gradlog)
        !
        DEALLOCATE (dlogeps)
        !
        IF (dielectric%need_gradient) THEN
            !
            DO ipol = 1, 3
                gradeps(ipol, :) = deps(:) * gradscaled(ipol, :)
                !
                IF (dielectric%nregions > 0) &
                    gradeps(ipol, :) = gradeps(ipol, :) + &
                                       deps_dback(:) * gradback(ipol, :)
                !
            END DO
            !
            CALL update_gradient_modulus(dielectric%gradient)
            !
        END IF
        !
        IF (dielectric%need_factsqrt) THEN
            !
            IF (dielectric%nregions <= 0) THEN
                !
                factsqrteps = (d2eps - 0.5D0 * deps**2 / eps) * gradscaledmod**2 + &
                              deps * laplscaled
                !
            ELSE
                !
                CALL scalar_product_environ_gradient(dielectric%boundary%gradient, &
                                                     dielectric%gradbackground, &
                                                     dielectric%factsqrt)
                !
                IF (dielectric%need_gradient) THEN
                    gradepsmod2 = dielectric%gradient%modulus%of_r**2
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
        IF (dielectric%nregions > 0) THEN
            DEALLOCATE (deps_dback)
            DEALLOCATE (dlogeps_dback)
            !
            IF (dielectric%need_factsqrt) THEN
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
    !
    !------------------------------------------------------------------------------------
END MODULE utils_dielectric
!----------------------------------------------------------------------------------------
