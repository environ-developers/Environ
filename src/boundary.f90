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
    USE environ_param, ONLY: DP, e2
    !
    USE class_cell
    USE class_density
    USE class_function_erfc
    USE class_gradient
    USE class_hessian
    !
    USE class_core_container
    !
    USE class_electrons
    USE class_ions
    USE class_system
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
    TYPE, ABSTRACT, PUBLIC :: environ_boundary
        !--------------------------------------------------------------------------------
        !
        CHARACTER(LEN=80) :: label ! boundary label
        CHARACTER(LEN=80) :: mode ! choice of the interface
        INTEGER :: update_status = 0
        !
        !--------------------------------------------------------------------------------
        !
        TYPE(environ_cell), POINTER :: cell => NULL()
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
        !
        !--------------------------------------------------------------------------------
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
        ! Components needed for solvent-aware boundary
        !
        LOGICAL :: solvent_aware = .FALSE.
        TYPE(environ_function_erfc) :: solvent_probe
        REAL(DP) :: filling_threshold, filling_spread
        !
        TYPE(environ_hessian) :: hessian
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
        !
        REAL(DP) :: field_factor
        REAL(DP) :: field_asymmetry
        REAL(DP) :: field_max
        REAL(DP) :: field_min
        !
        !--------------------------------------------------------------------------------
    CONTAINS
        !--------------------------------------------------------------------------------
        ! Admin
        !
        PROCEDURE, PRIVATE :: pre_create => pre_create_environ_boundary
        PROCEDURE :: pre_init => pre_init_environ_boundary
        PROCEDURE :: init_solvent_aware
        PROCEDURE :: update_solvent_aware
        PROCEDURE :: pre_destroy => pre_destroy_environ_boundary
        !
        PROCEDURE(create_boundary), DEFERRED :: create
        PROCEDURE(update_boundary), DEFERRED :: update
        PROCEDURE(destroy_boundary), DEFERRED :: destroy
        PROCEDURE(build_boundary), DEFERRED :: build
        !
        !--------------------------------------------------------------------------------
        ! Calculators
        !
        PROCEDURE :: vconfine => calc_vconfine
        PROCEDURE :: evolume => calc_evolume
        PROCEDURE :: esurface => calc_esurface
        !
        PROCEDURE, NOPASS :: deconfine_dboundary => calc_deconfine_dboundary
        PROCEDURE, NOPASS :: devolume_dboundary => calc_devolume_dboundary
        PROCEDURE :: desurface_dboundary => calc_desurface_dboundary
        PROCEDURE :: sa_de_dboundary => calc_solvent_aware_de_dboundary
        PROCEDURE :: dboundary_dions => calc_dboundary_dions
        PROCEDURE :: fa_dboundary_dions => calc_field_aware_dboundary_dions
        PROCEDURE :: ion_field_partial => calc_ion_field_partial
        !
        !--------------------------------------------------------------------------------
        ! Solvent aware
        !
        PROCEDURE :: solvent_aware_boundary
        !
        !--------------------------------------------------------------------------------
        ! Private helpers
        !
        PROCEDURE :: convolution => compute_convolution_deriv
        PROCEDURE :: calc_dsurface ! #TODO do we need this?
        PROCEDURE :: invert => invert_boundary
        !
        !--------------------------------------------------------------------------------
        ! Output
        !
        PROCEDURE :: pre_printout => pre_print_environ_boundary
        PROCEDURE, NOPASS :: print_setup
        !
        !--------------------------------------------------------------------------------
    END TYPE environ_boundary
    !------------------------------------------------------------------------------------
    !
    ABSTRACT INTERFACE
        SUBROUTINE create_boundary(this)
            IMPORT environ_boundary
            CLASS(environ_boundary), INTENT(INOUT) :: this
        END SUBROUTINE
        SUBROUTINE update_boundary(this)
            IMPORT environ_boundary
            CLASS(environ_boundary), INTENT(INOUT) :: this
        END SUBROUTINE
        SUBROUTINE destroy_boundary(this)
            IMPORT environ_boundary
            CLASS(environ_boundary), INTENT(INOUT) :: this
        END SUBROUTINE
        SUBROUTINE build_boundary(this, density)
            IMPORT environ_boundary, environ_density
            TYPE(environ_density), OPTIONAL, TARGET, INTENT(IN) :: density
            CLASS(environ_boundary), TARGET, INTENT(INOUT) :: this
        END SUBROUTINE
    END INTERFACE
    !
    !------------------------------------------------------------------------------------
    !
    REAL(DP), PARAMETER, PUBLIC :: tolspuriousforce = 1.D-5
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
    SUBROUTINE pre_create_environ_boundary(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_boundary), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: sub_name = 'pre_create_environ_boundary'
        !
        !--------------------------------------------------------------------------------
        !
        IF (ASSOCIATED(this%cell)) CALL io%create_error(sub_name)
        !
        IF (ASSOCIATED(this%cores)) CALL io%create_error(sub_name)
        !
        !--------------------------------------------------------------------------------
        !
        this%label = ''
        this%mode = ''
        this%update_status = 0
        this%deriv = 0
        this%derivatives_method = ''
        this%volume = 0.D0
        this%surface = 0.D0
        this%solvent_aware = .FALSE.
        this%field_aware = .FALSE.
        this%field_factor = 0.D0
        this%field_asymmetry = 0.D0
        this%field_max = 0.D0
        this%field_min = 0.D0
        !
        NULLIFY (this%cell)
        !
        NULLIFY (this%cores)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE pre_create_environ_boundary
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE pre_init_environ_boundary(this, mode, need_gradient, need_laplacian, &
                                         need_hessian, field_aware, field_factor, &
                                         field_asymmetry, field_max, field_min, &
                                         cores, deriv_method, cell, label)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CHARACTER(LEN=*), INTENT(IN) :: mode
        !
        LOGICAL, INTENT(IN) :: need_gradient
        LOGICAL, INTENT(IN) :: need_laplacian
        LOGICAL, INTENT(IN) :: need_hessian
        !
        LOGICAL, INTENT(IN) :: field_aware
        REAL(DP), INTENT(IN) :: field_factor
        REAL(DP), INTENT(IN) :: field_asymmetry
        REAL(DP), INTENT(IN) :: field_max
        REAL(DP), INTENT(IN) :: field_min
        !
        CHARACTER(LEN=*), INTENT(IN) :: deriv_method
        !
        TYPE(core_container), TARGET, INTENT(IN) :: cores
        TYPE(environ_cell), TARGET, INTENT(IN) :: cell
        CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: label
        !
        CLASS(environ_boundary), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: sub_name = 'pre_init_environ_boundary'
        !
        !--------------------------------------------------------------------------------
        !
        CALL this%pre_create()
        !
        this%mode = mode
        !
        IF (PRESENT(label)) this%label = label
        !
        IF (need_hessian) THEN
            this%deriv = 3
        ELSE IF (need_laplacian) THEN
            this%deriv = 2
        ELSE IF (need_gradient) THEN
            this%deriv = 1
        END IF
        !
        this%cores => cores
        this%derivatives_method = deriv_method
        !
        this%cell => cell
        !
        CALL this%scaled%init(cell, 'boundary_'//this%label)
        !
        IF (this%deriv >= 1) CALL this%gradient%init(cell, 'gradboundary_'//this%label)
        !
        IF (this%deriv >= 2) CALL this%laplacian%init(cell, 'laplboundary_'//this%label)
        !
        IF (this%deriv >= 3) CALL this%dsurface%init(cell, 'dsurface_'//this%label)
        !
        !--------------------------------------------------------------------------------
        ! Field awareness
        !
        IF (field_aware) THEN
            this%field_aware = field_aware
            this%field_factor = field_factor
            this%field_asymmetry = field_asymmetry
            this%field_max = field_max
            this%field_min = field_min
        END IF
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE pre_init_environ_boundary
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_solvent_aware(this, solvent_radius, radial_scale, radial_spread, &
                                  filling_threshold, filling_spread)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP), INTENT(IN) :: solvent_radius
        REAL(DP), INTENT(IN) :: radial_scale
        REAL(DP), INTENT(IN) :: radial_spread
        REAL(DP), INTENT(IN) :: filling_threshold
        REAL(DP), INTENT(IN) :: filling_spread
        !
        CLASS(environ_boundary), TARGET, INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: sub_name = 'init_solvent_aware'
        !
        !--------------------------------------------------------------------------------
        !
        this%solvent_aware = .TRUE.
        !
        CALL this%solvent_probe%init(2, 1, 0, solvent_radius * radial_scale, &
                                     radial_spread, 1.D0)
        !
        this%filling_threshold = filling_threshold
        this%filling_spread = filling_spread
        !
        CALL this%local%init(this%cell, 'local_'//this%label)
        !
        CALL this%probe%init(this%cell, 'probe_'//this%label)
        !
        CALL this%filling%init(this%cell, 'filling_'//this%label)
        !
        CALL this%dfilling%init(this%cell, 'dfilling_'//this%label)
        !
        IF (this%deriv >= 3) &
            CALL this%hessian%init(this%cell, 'hessboundary_'//this%label)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_solvent_aware
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE update_solvent_aware(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_boundary), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: sub_name = 'update_solvent_aware'
        !
        !--------------------------------------------------------------------------------
        ! Solvent-aware interface
        !
        IF (this%update_status == 2 .AND. this%solvent_aware) &
            CALL this%solvent_aware_boundary()
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE update_solvent_aware
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE pre_destroy_environ_boundary(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_boundary), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: sub_name = 'pre_destroy_environ_boundary'
        !
        !--------------------------------------------------------------------------------
        !
        CALL this%scaled%destroy()
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
            DEALLOCATE (this%solvent_probe%pos)
        END IF
        !
        NULLIFY (this%cores)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE pre_destroy_environ_boundary
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
        CHARACTER(LEN=80) :: sub_name = 'calc_dboundary_dions'
        !
        !--------------------------------------------------------------------------------
        !
        CALL io%error(sub_name, "Not implemented", 1)
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
        CLASS(environ_boundary), INTENT(INOUT) :: this
        TYPE(environ_gradient), INTENT(INOUT) :: partial
        !
        CHARACTER(LEN=80) :: sub_name = 'calc_field_aware_dboundary_dions'
        !
        !--------------------------------------------------------------------------------
        !
        CALL io%error(sub_name, "Not implemented", 1)
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
        CLASS(environ_boundary), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: sub_name = 'calc_ion_field_partial'
        !
        !--------------------------------------------------------------------------------
        !
        CALL io%error(sub_name, "Not implemented", 1)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_ion_field_partial
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !                                  GENERAL METHODS
    !
    !------------------------------------------------------------------------------------
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
    !                               PRIVATE HELPER METHODS
    !
    !------------------------------------------------------------------------------------
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
    SUBROUTINE pre_print_environ_boundary(this, verbose, debug_verbose, unit)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_boundary), INTENT(IN) :: this
        INTEGER, OPTIONAL, INTENT(IN) :: verbose, debug_verbose, unit
        !
        INTEGER :: base_verbose, local_verbose, passed_verbose, local_unit, i
        !
        CHARACTER(LEN=80) :: sub_name = 'pre_print_environ_boundary'
        !
        !--------------------------------------------------------------------------------
        !
        IF (.NOT. PRESENT(debug_verbose) .AND. io%verbosity <= 0) RETURN
        !
        CALL this%print_setup(base_verbose, local_verbose, passed_verbose, local_unit, &
                              verbose, debug_verbose, unit)
        !
        !--------------------------------------------------------------------------------
        !
        IF (local_verbose >= 1) THEN
            !
            IF (io%lnode) THEN
                WRITE (local_unit, 1100)
                WRITE (local_unit, 1101) this%label, this%mode
            END IF
            !
            IF (io%lnode) THEN
                WRITE (local_unit, 1102) this%volume
                !
                IF (this%deriv >= 1) WRITE (local_unit, 1103) this%surface
                !
            END IF
            !
            IF (local_verbose >= 4) &
                CALL this%scaled%printout(passed_verbose, debug_verbose, local_unit)
            !
            IF (this%solvent_aware) THEN
                !
                IF (io%lnode) &
                    WRITE (local_unit, 1104) &
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
1100    FORMAT(/, 4('%'), " BOUNDARY ", 66('%'))
        !
1101    FORMAT(/, " boundary label             = ", A20, /, &
                " boundary mode              = ", A20)
        !
1102    FORMAT(/, " volume of the QM region    = ", F14.7)
        !
1103    FORMAT(/, " surface of the QM region   = ", F14.7)
        !
1104    FORMAT(/, " using solvent-aware boundary:", /, &
                " filling threshold          = ", F14.7, /, &
                " filling spread             = ", F14.7, /, &
                " solvent radius x rad scale = ", F14.7, /, &
                " spread of solvent probe    = ", F14.7)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE pre_print_environ_boundary
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE print_setup(base_verbose, local_verbose, passed_verbose, local_unit, &
                           verbose, debug_verbose, unit)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, OPTIONAL, INTENT(IN) :: verbose, debug_verbose, unit
        !
        INTEGER, INTENT(OUT) :: base_verbose, local_verbose, passed_verbose, local_unit
        !
        CHARACTER(LEN=80) :: sub_name = 'print_setup'
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
        END IF
        !
        IF (PRESENT(unit)) THEN
            local_unit = unit
        ELSE
            local_unit = io%debug_unit
        END IF
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE print_setup
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE class_boundary
!----------------------------------------------------------------------------------------
