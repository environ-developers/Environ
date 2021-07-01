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
!    GNU General Public License for more detail, either the file `License'
!    in the root directory of the present distribution, or
!    online at <http://www.gnu.org/licenses/>.
!
!----------------------------------------------------------------------------------------
!
! Authors: Oliviero Andreussi (Department of Physics, UNT)
!          Ismaila Dabo       (DMSE, Penn State)
!          Nicola Marzari     (THEOS and NCCR-MARVEL, EPFL)
!
!----------------------------------------------------------------------------------------
!>
!! This module contains all the procedures to generate the boundary
!! function and its derivative, either as a functional of the electronic
!! density (self-consistent boundary), or as a function of the ionic
!! positions (soft-sphere boundary), or a simple geometric surface
!! centered on the system position (system boundary)
!!
!----------------------------------------------------------------------------------------
MODULE generate_boundary
    !------------------------------------------------------------------------------------
    !
    USE env_base_io, ONLY: ionode, program_unit
    !
    USE environ_param, ONLY: DP, sqrtpi, tpi
    !
    USE types_representation
    USE types_physical, ONLY: environ_boundary
    USE types_core, ONLY: fft_core, fd_core
    USE types_cell, ONLY: environ_cell
    !
    USE utils_hessian, ONLY: init_environ_hessian, destroy_environ_hessian
    USE utils_density, ONLY: init_environ_density, destroy_environ_density
    !
    USE utils_gradient, ONLY: init_environ_gradient, update_gradient_modulus, &
                              destroy_environ_gradient
    !
    USE tools_functions
    USE tools_fft, ONLY: gradient_fft, laplacian_fft, hessian_fft, convolution_fft
    USE tools_fd, ONLY: gradient_fd
    !
    USE tools_math, ONLY: scalar_product_environ_gradient, integrate_environ_density, &
                          environ_erfc
    !
    !------------------------------------------------------------------------------------
    !
    PRIVATE
    !
    PUBLIC :: boundary_of_density, boundary_of_functions, boundary_of_system, &
              calc_dboundary_dions, invert_boundary, solvent_aware_boundary, &
              solvent_aware_de_dboundary
    !
    !------------------------------------------------------------------------------------
CONTAINS
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
        RETURN
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
        RETURN
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
        RETURN
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
        RETURN
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
        RETURN
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
        RETURN
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
        RETURN
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
        RETURN
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
            CALL env_errore(fun_name, 'Unknown boundary type', 1)
            !
        END SELECT
        !
        RETURN
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
            CALL env_errore(fun_name, 'Unknown boundary type', 1)
            !
        END SELECT
        !
        RETURN
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
            CALL env_errore(fun_name, 'Option not yet implemented', 1)
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
            CALL env_errore(fun_name, 'Unknown boundary type', 1)
            !
        END SELECT
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END FUNCTION d2boundfunct
    !------------------------------------------------------------------------------------
    !>
    !! Calculates the dielectric constant as a function of the charge
    !! density, and the derivatives of the the dielectric constant
    !! with respect to the charge density. Additionally calculates the
    !! volume and surface components.
    !!
    !! @param[in]     density              electrostatic density to use
    !! @param[inout]  boundary             boundary to update
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE boundary_of_density(density, boundary)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_density), TARGET, INTENT(IN) :: density
        !
        TYPE(environ_boundary), TARGET, INTENT(INOUT) :: boundary
        !
        INTEGER, POINTER :: ir_end, nnr, stype, deriv
        REAL(DP), POINTER :: const, rhomax, rhomin, tbeta
        REAL(DP), DIMENSION(:), POINTER :: rho, eps, deps, d2eps, lapleps, dsurface
        REAL(DP), POINTER :: gradeps(:, :)
        TYPE(environ_hessian), POINTER :: hessian
        !
        TYPE(fft_core), POINTER :: fft
        TYPE(fd_core), POINTER :: fd
        !
        INTEGER :: ir, ipol, jpol
        !
        CHARACTER(LEN=80) :: sub_name = 'boundary_of_density'
        !
        !--------------------------------------------------------------------------------
        !
        IF (.NOT. ASSOCIATED(density%cell, boundary%scaled%cell)) &
            CALL env_errore(sub_name, 'Inconsistent domains', 1)
        !
        ir_end => density%cell%ir_end
        nnr => density%cell%nnr
        rho => density%of_r
        !
        stype => boundary%type_
        eps => boundary%scaled%of_r
        deps => boundary%dscaled%of_r
        d2eps => boundary%d2scaled%of_r
        !
        fft => boundary%core%fft
        fd => boundary%core%fd
        !
        IF (stype == 1 .OR. stype == 2) THEN
            rhomax => boundary%rhomax
            rhomin => boundary%rhomin
            tbeta => boundary%fact
            const => boundary%const
        ELSE IF (stype == 0) THEN
            rhomax => boundary%rhozero
            rhomin => boundary%deltarho
            tbeta => boundary%tbeta
            const => boundary%const
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
        deriv => boundary%deriv
        !
        IF (deriv >= 1) gradeps => boundary%gradient%of_r
        !
        IF (deriv >= 2) lapleps => boundary%laplacian%of_r
        !
        IF (deriv >= 3) THEN
            dsurface => boundary%dsurface%of_r
            !
            IF (boundary%solvent_aware) THEN
                hessian => boundary%hessian
            ELSE
                ALLOCATE (hessian)
                !
                CALL init_environ_hessian(density%cell, hessian)
                !
            END IF
            !
        END IF
        !
        SELECT CASE (boundary%core%type_)
            !
        CASE ('fft')
            !
            IF (deriv == 1 .OR. deriv == 2) &
                CALL gradient_fft(fft, boundary%scaled, boundary%gradient)
            !
            IF (deriv == 2) &
                CALL laplacian_fft(fft, boundary%scaled, boundary%laplacian)
            !
            IF (deriv == 3) &
                CALL dsurface_fft(fft, boundary%scaled, boundary%gradient, &
                                  boundary%laplacian, hessian, boundary%dsurface)
            !
        CASE ('chain', 'fd')
            !
            IF (deriv == 1 .OR. deriv == 2) &
                CALL gradient_fft(fft, density, boundary%gradient)
            !
            IF (deriv == 2) CALL laplacian_fft(fft, density, boundary%laplacian)
            !
            IF (deriv == 3) THEN
                !
                CALL dsurface_fft(fft, density, boundary%gradient, &
                                  boundary%laplacian, hessian, boundary%dsurface)
                !
                IF (boundary%solvent_aware) THEN
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
                IF (boundary%core%type_ == 'chain') THEN
                    !
                    DO ipol = 1, 3
                        gradeps(ipol, :) = gradeps(ipol, :) * deps(:)
                    END DO
                    !
                ELSE IF (boundary%core%type_ == 'fd') THEN
                    CALL gradient_fd(fd, boundary%scaled, boundary%gradient)
                END IF
                !
            END IF
            !
        END SELECT
        !
        !--------------------------------------------------------------------------------
        ! Final updates
        !
        boundary%volume = integrate_environ_density(boundary%scaled)
        !
        IF (deriv >= 1) THEN
            !
            CALL update_gradient_modulus(boundary%gradient)
            !
            boundary%surface = integrate_environ_density(boundary%gradient%modulus)
        END IF
        !
        IF (deriv >= 3 .AND. .NOT. boundary%solvent_aware) &
            CALL destroy_environ_hessian(hessian)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE boundary_of_density
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE dsurface_fft(fft, x, grad, lapl, hess, dsurface)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(fft_core), INTENT(IN) :: fft
        TYPE(environ_density), INTENT(IN) :: x
        !
        TYPE(environ_gradient), INTENT(INOUT) :: grad
        TYPE(environ_density), INTENT(INOUT) :: lapl
        TYPE(environ_hessian), INTENT(INOUT) :: hess
        TYPE(environ_density), INTENT(INOUT) :: dsurface
        !
        !--------------------------------------------------------------------------------
        !
        CALL hessian_fft(fft, x, grad, hess)
        !
        lapl%of_r(:) = hess%of_r(1, 1, :) + hess%of_r(2, 2, :) + hess%of_r(3, 3, :)
        !
        CALL calc_dsurface(x%cell%nnr, x%cell%ir_end, grad%of_r, hess%of_r, &
                           dsurface%of_r)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE dsurface_fft
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_dsurface(n, iend, grad, hess, dsurface)
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
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_dsurface
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
    !! @param[in]      nsoft_spheres     number of soft-spheres
    !! @param[in]      soft_spheres      functions for the soft-spheres
    !! @param[inout]   boundary          boundary to update
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE boundary_of_functions(nsoft_spheres, soft_spheres, boundary)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: nsoft_spheres
        TYPE(environ_functions), INTENT(IN) :: soft_spheres(nsoft_spheres)
        !
        TYPE(environ_boundary), TARGET, INTENT(INOUT) :: boundary
        !
        INTEGER, POINTER :: nnr, ir_end, deriv
        TYPE(environ_cell), POINTER :: cell
        !
        INTEGER :: i
        !
        TYPE(environ_density), ALLOCATABLE :: local(:)
        TYPE(environ_gradient), ALLOCATABLE :: gradlocal(:)
        TYPE(environ_density), ALLOCATABLE :: lapllocal(:)
        TYPE(environ_hessian), ALLOCATABLE :: hesslocal(:)
        TYPE(environ_hessian), POINTER :: hessian
        !
        TYPE(fft_core), POINTER :: fft
        !
        CHARACTER(LEN=80) :: label
        !
        CHARACTER(LEN=80) :: sub_name = 'boundary_of_functions'
        !
        !--------------------------------------------------------------------------------
        !
        cell => boundary%scaled%cell
        nnr => cell%nnr
        ir_end => cell%ir_end
        fft => boundary%core%fft
        !
        ALLOCATE (local(nsoft_spheres))
        !
        !--------------------------------------------------------------------------------
        ! Compute soft spheres and generate boundary
        !
        boundary%scaled%of_r = 1.D0
        !
        DO i = 1, nsoft_spheres
            !
            CALL init_environ_density(cell, local(i))
            !
            CALL density_of_functions(soft_spheres(i), local(i), .FALSE.)
            !
            boundary%scaled%of_r = boundary%scaled%of_r * local(i)%of_r
        END DO
        !
        !--------------------------------------------------------------------------------
        ! Generate boundary derivatives, if needed
        !
        deriv => boundary%deriv
        !
        IF (deriv == 3) THEN
            !
            IF (boundary%solvent_aware) THEN
                hessian => boundary%hessian
                hessian%of_r = 0.D0
            ELSE
                ALLOCATE (hessian)
                !
                CALL init_environ_hessian(cell, hessian)
                !
            END IF
            !
        END IF
        !
        SELECT CASE (boundary%core%type_)
            !
        CASE ('fft')
            !
            IF (deriv == 1 .OR. deriv == 2) &
                CALL gradient_fft(fft, boundary%scaled, boundary%gradient)
            !
            IF (deriv == 2) &
                CALL laplacian_fft(fft, boundary%scaled, boundary%laplacian)
            !
            IF (deriv == 3) &
                CALL dsurface_fft(fft, boundary%scaled, boundary%gradient, &
                                  boundary%laplacian, hessian, boundary%dsurface)
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
                IF (deriv >= 1) CALL init_environ_gradient(cell, gradlocal(i))
                !
                IF (deriv == 2) CALL init_environ_density(cell, lapllocal(i))
                !
                IF (deriv == 3) CALL init_environ_hessian(cell, hesslocal(i))
                !
                IF (deriv >= 1) &
                    CALL gradient_of_functions(soft_spheres(i), gradlocal(i), .FALSE.)
                !
                IF (deriv == 2) &
                    CALL laplacian_of_functions(soft_spheres(i), lapllocal(i), .FALSE.)
                !
                IF (deriv == 3) &
                    CALL hessian_of_functions(soft_spheres(i), hesslocal(i), .FALSE.)
                !
            END DO
            !
            IF (deriv == 1 .OR. deriv == 2) &
                CALL calc_gradient_of_boundary_highmem(nsoft_spheres, local, &
                                                       gradlocal, boundary%gradient)
            !
            IF (deriv == 2) &
                CALL calc_laplacian_of_boundary_highmem(nsoft_spheres, local, &
                                                        gradlocal, lapllocal, &
                                                        boundary%laplacian)
            !
            IF (deriv == 3) &
                CALL calc_dsurface_of_boundary_highmem(nsoft_spheres, local, &
                                                       gradlocal, hesslocal, &
                                                       boundary%gradient, &
                                                       boundary%laplacian, hessian, &
                                                       boundary%dsurface)
            !
            DO i = 1, nsoft_spheres
                !
                IF (deriv >= 1) CALL destroy_environ_gradient(gradlocal(i))
                !
                IF (deriv == 2) CALL destroy_environ_density(lapllocal(i))
                !
                IF (deriv == 3) CALL destroy_environ_hessian(hesslocal(i))
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
                IF (deriv >= 1) CALL init_environ_gradient(cell, gradlocal(i))
                !
                IF (deriv == 2) CALL init_environ_density(cell, lapllocal(i))
                !
                IF (deriv == 3) CALL init_environ_hessian(cell, hesslocal(i))
                !
                IF (deriv >= 1) &
                    CALL gradient_of_functions(soft_spheres(i), gradlocal(i), .FALSE.)
                !
                IF (deriv == 2) &
                    CALL laplacian_of_functions(soft_spheres(i), lapllocal(i), .FALSE.)
                !
                IF (deriv == 3) &
                    CALL hessian_of_functions(soft_spheres(i), hesslocal(i), .FALSE.)
                !
            END DO
            !
            IF (deriv >= 1) &
                CALL calc_gradient_of_boundary_lowmem(nsoft_spheres, local, &
                                                      gradlocal, boundary%scaled, &
                                                      boundary%gradient)
            !
            IF (deriv == 2) &
                CALL calc_laplacian_of_boundary_lowmem(nsoft_spheres, local, &
                                                       gradlocal, lapllocal, &
                                                       boundary%scaled, &
                                                       boundary%gradient, &
                                                       boundary%laplacian)
            !
            IF (deriv == 3) &
                CALL calc_dsurface_of_boundary_lowmem(nsoft_spheres, local, &
                                                      gradlocal, hesslocal, &
                                                      boundary%gradient, &
                                                      boundary%laplacian, hessian, &
                                                      boundary%scaled, &
                                                      boundary%dsurface)
            !
            DO i = 1, nsoft_spheres
                !
                IF (deriv >= 1) CALL destroy_environ_gradient(gradlocal(i))
                !
                IF (deriv == 2) CALL destroy_environ_density(lapllocal(i))
                !
                IF (deriv == 3) CALL destroy_environ_hessian(hesslocal(i))
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
        boundary%scaled%of_r = 1.D0 - boundary%scaled%of_r
        boundary%volume = integrate_environ_density(boundary%scaled)
        !
        IF (deriv >= 1) THEN
            boundary%gradient%of_r = -boundary%gradient%of_r
            !
            CALL update_gradient_modulus(boundary%gradient)
            !
            boundary%surface = integrate_environ_density(boundary%gradient%modulus)
            !
            IF (deriv >= 2) boundary%laplacian%of_r = -boundary%laplacian%of_r
            !
            IF (deriv == 3) THEN
                boundary%dsurface%of_r = -boundary%dsurface%of_r
                !
                IF (boundary%solvent_aware) THEN
                    boundary%hessian%of_r = -boundary%hessian%of_r
                ELSE
                    !
                    CALL destroy_environ_hessian(hessian)
                    !
                    DEALLOCATE (hessian)
                END IF
                !
            END IF
            !
        END IF
        !
        DO i = 1, nsoft_spheres
            CALL destroy_environ_density(local(i))
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
        IF (i > n) CALL env_errore(sub_name, 'Index out of bound', 1)
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
        RETURN
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
        CALL init_environ_gradient(cell, partial)
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
        CALL destroy_environ_gradient(partial)
        !
        RETURN
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
        CALL init_environ_density(cell, tmp)
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
                    !
                    CALL scalar_product_environ_gradient(gradlocal(i), gradlocal(j), &
                                                         tmp)
                    !
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
        CALL destroy_environ_density(tmp)
        !
        RETURN
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
        CALL init_environ_density(cell, dens)
        !
        CALL init_environ_gradient(cell, partial)
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
        CALL calc_dsurface(cell%nnr, cell%ir_end, gradient%of_r, hessian%of_r, &
                           dsurface%of_r)
        !
        CALL destroy_environ_density(dens)
        !
        CALL destroy_environ_gradient(partial)
        !
        RETURN
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
        INTEGER :: i, j, ipol, tol
        TYPE(environ_cell), POINTER :: cell
        !
        !--------------------------------------------------------------------------------
        !
        cell => gradient%cell
        !
        gradient%of_r = 0.D0
        tol = 1.D-60
        !
        !--------------------------------------------------------------------------------
        ! Temporary quotient
        !
        DO i = 1, n
            !
            DO j = 1, cell%nnr
                !
                IF (ABS(local(i)%of_r(j)) <= tol) CYCLE
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
        RETURN
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
        INTEGER :: i, j, k, ipol, tol
        TYPE(environ_cell), POINTER :: cell
        !
        !--------------------------------------------------------------------------------
        !
        cell => laplacian%cell
        tol = 1.D-60
        !
        DO i = 1, n
            !
            DO j = 1, cell%nnr
                !
                IF (ABS(local(i)%of_r(j)) <= tol) CYCLE
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
        RETURN
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
        INTEGER :: i, j, k, ipol, jpol, tol
        TYPE(environ_cell), POINTER :: cell
        !
        !--------------------------------------------------------------------------------
        !
        cell => laplacian%cell
        tol = 1.D-60
        !
        DO i = 1, n
            !
            DO j = 1, cell%nnr
                !
                IF (ABS(local(i)%of_r(j)) <= tol) CYCLE
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
        CALL calc_dsurface(cell%nnr, cell%ir_end, gradient%of_r, hessian%of_r, &
                           dsurface%of_r)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_dsurface_of_boundary_lowmem
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_dboundary_dions(index, boundary, partial)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: index
        TYPE(environ_boundary), INTENT(IN), TARGET :: boundary
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
        IF (boundary%mode == 'electronic') RETURN
        ! exit if boundary is only defined on electronic density
        !
        cell => partial%cell
        !
        IF (boundary%need_ions) THEN
            number => boundary%ions%number
        ELSE IF (boundary%need_system) THEN
            number => boundary%system%ions%number
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
        IF (boundary%mode == 'ionic' .AND. &
            .NOT. ALLOCATED(boundary%soft_spheres)) &
            CALL env_errore(sub_name, 'Missing details of ionic boundary', 1)
        !
        IF (boundary%mode == 'full' .AND. &
            .NOT. ALLOCATED(boundary%ions%core_electrons)) &
            CALL env_errore(sub_name, 'Missing details of core electrons', 1)
        !
        IF (boundary%mode == 'full' .AND. &
            .NOT. ASSOCIATED(boundary%dscaled%cell, cell)) &
            CALL env_errore(sub_name, 'Mismatch or unassociated boundary derivative', 1)
        !
        IF (boundary%mode == 'ionic' .OR. boundary%mode == 'fa-ionic') THEN
            !
            CALL gradient_of_functions(boundary%soft_spheres(index), partial, &
                                       .TRUE.)
            !
            CALL init_environ_density(cell, local)
            !
            DO i = 1, number
                !
                IF (i == index) CYCLE
                !
                CALL density_of_functions(boundary%soft_spheres(i), local, .TRUE.)
                !
                DO ipol = 1, 3
                    partial%of_r(ipol, :) = partial%of_r(ipol, :) * local%of_r(:)
                END DO
                !
            END DO
            !
            CALL destroy_environ_density(local)
            !
        ELSE IF (boundary%mode == 'full') THEN
            !
            CALL gradient_of_functions(boundary%ions%core_electrons(index), partial, &
                                       .TRUE.)
            !
            DO ipol = 1, 3
                !
                partial%of_r(ipol, :) = &
                    -partial%of_r(ipol, :) * boundary%dscaled%of_r(:)
                !
            END DO
            !
            CALL update_gradient_modulus(partial)
            !
            spurious_force = integrate_environ_density(partial%modulus)
            !
            IF (spurious_force > tolspuriousforce .AND. ionode) &
                WRITE (program_unit, 4001) index, spurious_force ! #TODO use env_warning?
            !
4001        FORMAT(1X, 'WARNING: Unphysical forces due to core electrons are non-negligible ' &
                   /, 1X, 'atom type ', I3, ' is subject to a spurious force of ', F12.6, ' ')
            !
        ELSE IF (boundary%mode == 'system') THEN
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
    !! Updates the boundary using a function
    !!
    !! Calculates the dielectric constant as a function of the charge
    !! density, and the derivatives of the dielectric constant with
    !! respect to the charge density. Also updates the volume and surface
    !! components. Expects an explicity defined system density function.
    !!
    !! @param[in]    simple         the input function
    !! @param[inout] boundary       the boundary to be updated
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE boundary_of_system(simple, boundary)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_functions), INTENT(IN) :: simple
        !
        TYPE(environ_boundary), TARGET, INTENT(INOUT) :: boundary
        !
        INTEGER, POINTER :: nnr, ir_end, deriv
        TYPE(environ_cell), POINTER :: cell
        TYPE(fft_core), POINTER :: fft
        !
        TYPE(environ_hessian), POINTER :: hesslocal
        !
        INTEGER :: i
        !
        CHARACTER(LEN=80) :: sub_name = 'boundary_of_system'
        !
        !--------------------------------------------------------------------------------
        !
        cell => boundary%scaled%cell
        nnr => cell%nnr
        ir_end => cell%ir_end
        !
        fft => boundary%core%fft
        !
        CALL density_of_functions(simple, boundary%scaled, .TRUE.)
        ! compute soft spheres and generate boundary
        !
        !--------------------------------------------------------------------------------
        ! Generate boundary derivatives, if needed
        !
        deriv => boundary%deriv
        !
        IF (deriv >= 3) THEN
            !
            IF (boundary%solvent_aware) THEN
                hesslocal => boundary%hessian
            ELSE
                ALLOCATE (hesslocal)
                !
                CALL init_environ_hessian(cell, hesslocal)
                !
            END IF
            !
        END IF
        !
        SELECT CASE (boundary%core%type_)
            !
        CASE ('fft')
            !
            IF (deriv == 1 .OR. deriv == 2) &
                CALL gradient_fft(fft, boundary%scaled, boundary%gradient)
            !
            IF (deriv == 2) CALL laplacian_fft(fft, boundary%scaled, boundary%laplacian)
            !
            IF (deriv == 3) &
                CALL dsurface_fft(fft, boundary%scaled, boundary%gradient, &
                                  boundary%laplacian, hesslocal, boundary%dsurface)
            !
        CASE ('chain')
            !
            IF (deriv >= 1) CALL gradient_of_functions(simple, boundary%gradient, .TRUE.)
            !
            IF (deriv >= 2) &
                CALL laplacian_of_functions(simple, boundary%laplacian, .TRUE.)
            !
            IF (deriv >= 3) THEN
                !
                CALL hessian_of_functions(simple, hesslocal, .TRUE.)
                !
                CALL calc_dsurface(nnr, ir_end, boundary%gradient%of_r, &
                                   hesslocal%of_r, boundary%dsurface%of_r)
                !
            END IF
            !
        END SELECT
        !
        IF (deriv >= 3) THEN
            !
            IF (.NOT. boundary%solvent_aware) THEN
                !
                CALL destroy_environ_hessian(hesslocal)
                !
                DEALLOCATE (hesslocal)
            END IF
            !
        END IF
        !
        boundary%volume = integrate_environ_density(boundary%scaled)
        !
        IF (deriv >= 1) THEN
            !
            CALL update_gradient_modulus(boundary%gradient)
            !
            boundary%surface = integrate_environ_density(boundary%gradient%modulus)
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
    SUBROUTINE invert_boundary(boundary)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_boundary), INTENT(INOUT) :: boundary
        !
        !--------------------------------------------------------------------------------
        !
        boundary%scaled%of_r = 1.D0 - boundary%scaled%of_r
        !
        boundary%volume = integrate_environ_density(boundary%scaled)
        !
        IF (boundary%deriv >= 1) boundary%gradient%of_r = -boundary%gradient%of_r
        !
        IF (boundary%deriv >= 2) boundary%laplacian%of_r = -boundary%laplacian%of_r
        !
        IF (boundary%deriv >= 3) THEN
            boundary%dsurface%of_r = -boundary%dsurface%of_r
            !
            IF (boundary%solvent_aware) boundary%hessian%of_r = -boundary%hessian%of_r
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
    SUBROUTINE solvent_aware_boundary(boundary)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_boundary), INTENT(INOUT), TARGET :: boundary
        !
        INTEGER, POINTER :: nnr, ir_end, deriv
        REAL(DP), POINTER :: thr, spr
        TYPE(environ_cell), POINTER :: cell
        TYPE(fft_core), POINTER :: fft
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
        !--------------------------------------------------------------------------------
        !
        cell => boundary%scaled%cell
        nnr => boundary%scaled%cell%nnr
        ir_end => boundary%scaled%cell%ir_end
        deriv => boundary%deriv
        !
        thr => boundary%filling_threshold
        spr => boundary%filling_spread
        !
        fft => boundary%core%fft
        !
        CALL init_environ_density(cell, filled_fraction)
        !
        IF (deriv >= 2 .AND. boundary%core%type_ /= 'fft') &
            CALL init_environ_density(cell, d2filling)
        !
        !--------------------------------------------------------------------------------
        ! Step 0: save local interface function for later use
        !
        boundary%local%of_r = boundary%scaled%of_r
        !
        !--------------------------------------------------------------------------------
        ! Step 1: compute the convolution function, this may be made moved out of here
        !
        CALL density_of_functions(boundary%solvent_probe, boundary%probe, .TRUE.)
        !
        boundary%probe%of_r = boundary%probe%of_r / &
                              integrate_environ_density(boundary%probe)
        !
        !--------------------------------------------------------------------------------
        ! Step 2: compute filled fraction, i.e. convolution of local boundary with probe
        !
        CALL convolution_fft(fft, boundary%local, boundary%probe, filled_fraction)
        !
        !--------------------------------------------------------------------------------
        ! Step 3: compute the filling function and its derivative
        !
        boundary%filling%of_r = 0.D0
        boundary%dfilling%of_r = 0.D0
        !
        DO ir = 1, ir_end
            !
            boundary%filling%of_r(ir) = 1.D0 - &
                                        sfunct2(filled_fraction%of_r(ir), thr, spr)
            !
            boundary%dfilling%of_r(ir) = -dsfunct2(filled_fraction%of_r(ir), thr, spr)
            !
            IF (deriv >= 2 .AND. boundary%core%type_ /= 'fft') &
                d2filling%of_r(ir) = -d2sfunct2(filled_fraction%of_r(ir), thr, spr)
            !
        END DO
        !
        !--------------------------------------------------------------------------------
        ! Step 4: compute solvent-aware interface
        !
        boundary%scaled%of_r = boundary%local%of_r + &
                               (1.D0 - boundary%local%of_r) * boundary%filling%of_r
        !
        !--------------------------------------------------------------------------------
        ! Step 5: compute boundary derivatives, if needed
        !
        SELECT CASE (boundary%core%type_)
            !
        CASE ('fft')
            !
            IF (deriv == 1 .OR. deriv == 2) &
                CALL gradient_fft(fft, boundary%scaled, boundary%gradient)
            !
            IF (deriv == 2) &
                CALL laplacian_fft(fft, boundary%scaled, boundary%laplacian)

            IF (deriv == 3) &
                CALL dsurface_fft(fft, boundary%scaled, boundary%gradient, &
                                  boundary%laplacian, boundary%hessian, &
                                  boundary%dsurface)
            !
        CASE ('chain', 'highmem')
            !
            !----------------------------------------------------------------------------
            ! Allocate local fields for derivatives of convolution
            !
            IF (deriv >= 1) CALL init_environ_gradient(cell, gradlocal)
            !
            IF (deriv >= 2) CALL init_environ_density(cell, lapllocal)
            !
            IF (deriv >= 3) CALL init_environ_hessian(cell, hesslocal)
            !
            !----------------------------------------------------------------------------
            ! Compute derivative of convolution with probe
            !
            IF (deriv > 1) &
                CALL compute_convolution_deriv(deriv, boundary, gradlocal, lapllocal, &
                                               hesslocal)
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
                        boundary%hessian%of_r(ipol, jpol, :) = &
                            boundary%hessian%of_r(ipol, jpol, :) * &
                            (1.D0 - boundary%filling%of_r) - boundary%dfilling%of_r * &
                            (boundary%gradient%of_r(ipol, :) * &
                             gradlocal%of_r(jpol, :) + &
                             boundary%gradient%of_r(jpol, :) * &
                             gradlocal%of_r(ipol, :)) + &
                            (1.D0 - boundary%local%of_r) * &
                            (d2filling%of_r * gradlocal%of_r(ipol, :) * &
                             gradlocal%of_r(jpol, :) + &
                             boundary%dfilling%of_r * hesslocal%of_r(ipol, jpol, :))
                        !
                    END DO
                    !
                END DO
                !
                CALL destroy_environ_hessian(hesslocal)
                !
            END IF
            !
            IF (deriv >= 2) THEN
                !
                CALL init_environ_density(cell, local)
                !
                CALL scalar_product_environ_gradient(boundary%gradient, gradlocal, &
                                                     local)
                !
                boundary%laplacian%of_r = &
                    boundary%laplacian%of_r * (1.D0 - boundary%filling%of_r) - &
                    2.D0 * local%of_r * boundary%dfilling%of_r + &
                    (1.D0 - boundary%local%of_r) * &
                    (d2filling%of_r * gradlocal%modulus%of_r**2 + &
                     boundary%dfilling%of_r * lapllocal%of_r)
                !
                CALL destroy_environ_density(local)
                !
                CALL destroy_environ_density(lapllocal)
                !
                CALL destroy_environ_density(d2filling)
                !
            END IF
            !
            IF (deriv >= 1) THEN
                !
                DO ipol = 1, 3
                    !
                    boundary%gradient%of_r(ipol, :) = &
                        boundary%gradient%of_r(ipol, :) * &
                        (1.D0 - boundary%filling%of_r(:)) + &
                        gradlocal%of_r(ipol, :) * &
                        (1.D0 - boundary%local%of_r(:)) * &
                        boundary%dfilling%of_r(:)
                    !
                END DO
                !
                CALL destroy_environ_gradient(gradlocal)
                !
            END IF
            !
            !----------------------------------------------------------------------------
            ! Recompute dsurface, if needed
            !
            IF (deriv >= 3) THEN
                !
                CALL calc_dsurface(nnr, ir_end, boundary%gradient%of_r, &
                                   boundary%hessian%of_r, boundary%dsurface%of_r)
                !
            END IF
            !
        END SELECT
        !
        !--------------------------------------------------------------------------------
        ! Final updates
        !
        boundary%volume = integrate_environ_density(boundary%scaled)
        !
        IF (deriv >= 1) THEN
            !
            CALL update_gradient_modulus(boundary%gradient)
            !
            boundary%surface = integrate_environ_density(boundary%gradient%modulus)
        END IF
        !
        CALL destroy_environ_density(filled_fraction)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE solvent_aware_boundary
    !------------------------------------------------------------------------------------
    !>
    !! @brief Compute the functional derivative of the energy wrt the boundary
    !!
    !! @param[in]   boundary      the boundary that contains all necessary information
    !! @param[out]  de_dboundary  the computed derivative
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE solvent_aware_de_dboundary(boundary, de_dboundary)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_boundary), INTENT(IN), TARGET :: boundary
        !
        TYPE(environ_density), INTENT(INOUT) :: de_dboundary
        !
        INTEGER, POINTER :: nnr, ir_end
        REAL(DP), POINTER :: thr, spr
        TYPE(environ_cell), POINTER :: cell
        !
        TYPE(environ_density) :: local
        !
        !--------------------------------------------------------------------------------
        !
        cell => boundary%scaled%cell
        nnr => boundary%scaled%cell%nnr
        !
        CALL init_environ_density(cell, local)
        !
        !--------------------------------------------------------------------------------
        ! Step 1: compute (1-s)*de_dboudary*dfilling
        !
        local%of_r = (1.D0 - boundary%local%of_r) * de_dboundary%of_r * &
                     boundary%dfilling%of_r
        !
        !--------------------------------------------------------------------------------
        ! Step 2: compute convolution with the probe function
        !
        CALL convolution_fft(boundary%core%fft, boundary%probe, local, local)
        !
        !--------------------------------------------------------------------------------
        ! Step 3: update the functional derivative of the energy wrt boundary
        !
        de_dboundary%of_r = de_dboundary%of_r * (1.D0 - boundary%filling%of_r) + &
                            local%of_r
        !
        CALL destroy_environ_density(local)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE solvent_aware_de_dboundary
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE compute_convolution_deriv(deriv, bound, grad, lapl, hess)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: deriv
        TYPE(environ_boundary), INTENT(IN) :: bound
        !
        TYPE(environ_gradient), INTENT(INOUT) :: grad
        TYPE(environ_density), INTENT(INOUT) :: lapl
        TYPE(environ_hessian), INTENT(INOUT) :: hess
        !
        !--------------------------------------------------------------------------------
        !
        IF (deriv <= 0) RETURN
        !
        IF (deriv >= 1) THEN
            !
            CALL convolution_fft(bound%core%fft, bound%probe, bound%gradient, grad)
            !
            CALL update_gradient_modulus(grad)
            !
        END IF
        !
        IF (deriv >= 2) &
            CALL convolution_fft(bound%core%fft, bound%probe, bound%laplacian, lapl)
        !
        IF (deriv >= 3) &
            CALL convolution_fft(bound%core%fft, bound%probe, bound%hessian, hess)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE compute_convolution_deriv
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE generate_boundary
!----------------------------------------------------------------------------------------
