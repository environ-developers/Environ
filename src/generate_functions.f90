!----------------------------------------------------------------------------------------
!
! Copyright (C) 2018-2021 ENVIRON (www.quantum-environ.org)
! Copyright (C) 2006-2010 Quantum ESPRESSO (www.quantum-espresso.org)
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
!          Edan Bainglass     (Department of Physics, UNT)
!
!----------------------------------------------------------------------------------------
!>
!! Module to generate functions on the real space dense grid
!!
!----------------------------------------------------------------------------------------
MODULE generate_functions
    !------------------------------------------------------------------------------------
    !
    USE env_mp, ONLY: env_mp_sum
    !
    USE environ_param, ONLY: DP, sqrtpi, pi, fpi
    !
    USE class_cell
    USE class_density
    USE class_gradient
    USE class_hessian
    !
    USE tools_math, ONLY: environ_erfc, environ_erf
    !
    !------------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    PRIVATE
    !
    PUBLIC :: generate_gaussian, generate_gradgaussian, generate_exponential, &
              generate_gradexponential, generate_erfc, generate_graderfc, &
              generate_laplerfc, generate_hesserfc, generate_axis, generate_distance, &
              erfcvolume
    !
    !------------------------------------------------------------------------------------
    ! NOTE: the spread of H core electrons (currently set to tol) must remain small, 
    !       but no smaller than tol. 
    !
    REAL(DP), PARAMETER :: tol = 1.D-10
    REAL(DP), PARAMETER :: exp_tol = 4.D1
    !
    !------------------------------------------------------------------------------------
CONTAINS
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE generate_gaussian(dim, axis, charge, spread, pos, density)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: dim, axis
        REAL(DP), INTENT(IN) :: charge, spread
        REAL(DP), INTENT(IN) :: pos(3)
        !
        CLASS(environ_density), INTENT(INOUT) :: density
        !
        LOGICAL :: physical
        INTEGER :: ir
        REAL(DP) :: scale, spr2, length
        REAL(DP) :: r(3), r2
        REAL(DP), ALLOCATABLE :: local(:)
        !
        CLASS(environ_cell), POINTER :: cell
        !
        CHARACTER(LEN=80) :: sub_name = 'generate_gaussian'
        !
        !--------------------------------------------------------------------------------
        !
        cell => density%cell
        !
        IF (ABS(charge) < tol) RETURN
        !
        IF (ABS(spread) < tol) &
            CALL env_errore(sub_name, 'Wrong spread for Gaussian function', 1)
        !
        IF (axis < 1 .OR. axis > 3) CALL env_errore(sub_name, 'Wrong value of axis', 1)
        !
        SELECT CASE (dim)
            !
        CASE (0)
            scale = charge / (sqrtpi * spread)**3
            !
        CASE (1)
            length = ABS(cell%at(axis, axis))
            scale = charge / length / (sqrtpi * spread)**2
            !
        CASE (2)
            length = ABS(cell%at(axis, axis))
            scale = charge * length / cell%omega / (sqrtpi * spread)
            !
        CASE default
            CALL env_errore(sub_name, 'Wrong value of dim', 1)
            !
        END SELECT
        !
        spr2 = spread**2
        !
        ALLOCATE (local(cell%nnr))
        local = 0.D0
        !
        DO ir = 1, cell%ir_end
            !
            CALL cell%get_min_distance(ir, dim, axis, pos, r, r2, physical)
            ! compute minimum distance using minimum image convention
            !
            IF (.NOT. physical) CYCLE
            !
            r2 = r2 / spr2
            !
            IF (r2 <= exp_tol) local(ir) = EXP(-r2) ! compute Gaussian function
            !
        END DO
        !
        density%of_r = density%of_r + scale * local
        DEALLOCATE (local)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE generate_gaussian
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE generate_gradgaussian(dim, axis, charge, spread, pos, gradient)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: dim, axis
        REAL(DP), INTENT(IN) :: charge, spread
        REAL(DP), INTENT(IN) :: pos(3)
        !
        CLASS(environ_gradient), INTENT(INOUT) :: gradient
        !
        LOGICAL :: physical
        INTEGER :: ir
        REAL(DP) :: scale, spr2, length
        REAL(DP) :: r(3), r2
        REAL(DP), ALLOCATABLE :: gradlocal(:, :)
        !
        CLASS(environ_cell), POINTER :: cell
        !
        CHARACTER(LEN=80) :: sub_name = 'generate_gradgaussian'
        !
        !--------------------------------------------------------------------------------
        !
        cell => gradient%cell ! sanity checks and initial setup
        !
        IF (ABS(charge) < tol) RETURN
        !
        IF (ABS(spread) < tol) &
            CALL env_errore(sub_name, 'Wrong spread for Gaussian function', 1)
        !
        IF (axis < 1 .OR. axis > 3) CALL env_errore(sub_name, 'Wrong value of axis', 1)
        !
        SELECT CASE (dim)
            !
        CASE (0)
            scale = charge / (sqrtpi * spread)**3
            !
        CASE (1)
            length = ABS(cell%at(axis, axis))
            scale = charge / length / (sqrtpi * spread)**2
            !
        CASE (2)
            length = ABS(cell%at(axis, axis))
            scale = charge * length / cell%omega / (sqrtpi * spread)
            !
        CASE default
            CALL env_errore(sub_name, 'Wrong value of dim', 1)
            !
        END SELECT
        !
        scale = scale * 2.D0 / spread**2
        !
        spr2 = spread**2
        !
        ALLOCATE (gradlocal(3, cell%nnr))
        gradlocal = 0.D0
        !
        DO ir = 1, cell%ir_end
            !
            CALL cell%get_min_distance(ir, dim, axis, pos, r, r2, physical)
            ! compute minimum distance using minimum image convention
            !
            IF (.NOT. physical) CYCLE
            !
            r2 = r2 / spr2
            !
            IF (r2 <= exp_tol) gradlocal(:, ir) = EXP(-r2) * r
            ! compute gradient of Gaussian function
            !
        END DO
        !
        gradient%of_r = gradient%of_r + scale * gradlocal
        DEALLOCATE (gradlocal)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE generate_gradgaussian
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE generate_exponential(dim, axis, width, spread, pos, density)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: dim, axis
        REAL(DP), INTENT(IN) :: width, spread
        REAL(DP), INTENT(IN) :: pos(3)
        !
        CLASS(environ_density), INTENT(INOUT) :: density
        !
        LOGICAL :: physical
        INTEGER :: ir
        REAL(DP) :: r2, dist, arg
        REAL(DP) :: r(3)
        REAL(DP), ALLOCATABLE :: local(:)
        !
        CLASS(environ_cell), POINTER :: cell
        !
        CHARACTER(LEN=80) :: sub_name = 'generate_exponential'
        !
        !--------------------------------------------------------------------------------
        !
        cell => density%cell
        !
        IF (axis < 1 .OR. axis > 3) CALL env_errore(sub_name, 'Wrong value of axis', 1)
        !
        ALLOCATE (local(cell%nnr))
        local = 0.D0
        !
        DO ir = 1, cell%ir_end
            !
            CALL cell%get_min_distance(ir, dim, axis, pos, r, r2, physical)
            ! compute minimum distance using minimum image convention
            !
            IF (.NOT. physical) CYCLE
            !
            dist = SQRT(r2)
            arg = (dist - width) / spread
            !
            IF (ABS(arg) <= exp_tol) local(ir) = EXP(-arg)
            ! compute exponentially decaying function
            !
        END DO
        !
        density%of_r = density%of_r + local
        DEALLOCATE (local)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE generate_exponential
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE generate_gradexponential(dim, axis, width, spread, pos, gradient)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: dim, axis
        REAL(DP), INTENT(IN) :: width, spread
        REAL(DP), INTENT(IN) :: pos(3)
        !
        CLASS(environ_gradient), INTENT(INOUT) :: gradient
        !
        LOGICAL :: physical
        INTEGER :: ir
        REAL(DP) :: r2, dist, arg
        REAL(DP) :: r(3)
        REAL(DP), ALLOCATABLE :: gradlocal(:, :)
        !
        CLASS(environ_cell), POINTER :: cell
        !
        CHARACTER(LEN=80) :: sub_name = 'generate_gradexponential'
        !
        !--------------------------------------------------------------------------------
        !
        cell => gradient%cell
        !
        IF (axis < 1 .OR. axis > 3) CALL env_errore(sub_name, 'Wrong value of axis', 1)
        !
        ALLOCATE (gradlocal(3, cell%nnr))
        gradlocal = 0.D0
        !
        DO ir = 1, cell%ir_end
            !
            CALL cell%get_min_distance(ir, dim, axis, pos, r, r2, physical)
            ! compute minimum distance using minimum image convention
            !
            IF (.NOT. physical) CYCLE
            !
            dist = SQRT(r2)
            arg = (dist - width) / spread
            !
            ! compute exponentially decaying function
            IF (r2 > tol .AND. ABS(arg) <= exp_tol) &
                gradlocal(:, ir) = r / SQRT(r2) / spread * EXP(-arg)
            !
        END DO
        !
        gradient%of_r = gradient%of_r + gradlocal
        DEALLOCATE (gradlocal)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE generate_gradexponential
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE generate_erfc(dim, axis, charge, width, spread, pos, density)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: dim, axis
        REAL(DP), INTENT(IN) :: charge, width, spread
        REAL(DP), INTENT(IN) :: pos(3)
        !
        CLASS(environ_density), INTENT(INOUT) :: density
        !
        LOGICAL :: physical
        INTEGER :: ir, ir_end, i
        REAL(DP) :: scale, r2, dist, arg, chargeanalytic, chargelocal
        REAL(DP) :: r(3)
        REAL(DP), ALLOCATABLE :: local(:)
        !
        CLASS(environ_cell), POINTER :: cell
        !
        CHARACTER(LEN=80) :: sub_name = 'generate_erfc'
        !
        !--------------------------------------------------------------------------------
        !
        cell => density%cell
        !
        IF (axis < 1 .OR. axis > 3) CALL env_errore(sub_name, 'Wrong value of axis', 1)
        !
        chargeanalytic = erfcvolume(dim, axis, width, spread, cell)
        scale = charge / chargeanalytic * 0.5D0
        !
        ALLOCATE (local(cell%nnr))
        local = 0.D0
        !
        DO ir = 1, cell%ir_end
            !
            CALL cell%get_min_distance(ir, dim, axis, pos, r, r2, physical)
            ! compute minimum distance using minimum image convention
            !
            IF (.NOT. physical) CYCLE
            !
            dist = SQRT(r2)
            arg = (dist - width) / spread
            !
            local(ir) = environ_erfc(arg) ! compute error function
            !
        END DO
        !
        !--------------------------------------------------------------------------------
        ! Check integral of function is consistent with analytic one
        !
        chargelocal = SUM(local) * cell%omega / DBLE(cell%ntot) * 0.5D0
        !
        CALL env_mp_sum(chargelocal, cell%dfft%comm)
        !
        IF (ABS(chargelocal - chargeanalytic) / chargeanalytic > 1.D-4) &
            CALL env_warning('wrong integral of erfc function')
        !
        !--------------------------------------------------------------------------------
        !
        density%of_r = density%of_r + scale * local
        ! rescale generated function to obtain the requested integral
        !
        DEALLOCATE (local)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE generate_erfc
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE generate_graderfc(dim, axis, charge, width, spread, pos, gradient)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: dim, axis
        REAL(DP), INTENT(IN) :: charge, width, spread
        REAL(DP), INTENT(IN) :: pos(3)
        !
        CLASS(environ_gradient), INTENT(INOUT) :: gradient
        !
        LOGICAL :: physical
        INTEGER :: ir
        REAL(DP) :: scale, r2, dist, arg, chargeanalytic
        REAL(DP) :: r(3)
        REAL(DP), ALLOCATABLE :: gradlocal(:, :)
        !
        CLASS(environ_cell), POINTER :: cell
        !
        CHARACTER(LEN=80) :: sub_name = 'generate_graderfc'
        !
        !--------------------------------------------------------------------------------
        !
        cell => gradient%cell
        !
        IF (axis < 1 .OR. axis > 3) CALL env_errore(sub_name, 'Wrong value of axis', 1)
        !
        chargeanalytic = erfcvolume(dim, axis, width, spread, cell)
        !
        scale = charge / chargeanalytic / sqrtpi / spread
        ! scaling factor, take into account rescaling of generated density
        ! to obtain the correct integrated total charge
        !
        ALLOCATE (gradlocal(3, cell%nnr))
        gradlocal = 0.D0
        !
        DO ir = 1, cell%ir_end
            !
            CALL cell%get_min_distance(ir, dim, axis, pos, r, r2, physical)
            ! compute minimum distance using minimum image convention
            !
            IF (.NOT. physical) CYCLE
            !
            dist = SQRT(r2)
            arg = (dist - width) / spread
            !
            IF (dist > tol) gradlocal(:, ir) = EXP(-arg**2) * r / dist
            ! compute gradient of error function
            !
        END DO
        !
        gradient%of_r = gradient%of_r + gradlocal * scale
        DEALLOCATE (gradlocal)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE generate_graderfc
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE generate_laplerfc(dim, axis, charge, width, spread, pos, laplacian)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: dim, axis
        REAL(DP), INTENT(IN) :: charge, width, spread
        REAL(DP), INTENT(IN) :: pos(3)
        !
        CLASS(environ_density), INTENT(INOUT) :: laplacian
        !
        LOGICAL :: physical
        INTEGER :: ir
        REAL(DP) :: scale, r2, dist, arg, chargeanalytic
        REAL(DP) :: r(3)
        REAL(DP), ALLOCATABLE :: lapllocal(:)
        !
        CLASS(environ_cell), POINTER :: cell
        !
        CHARACTER(LEN=80) :: sub_name = 'generate_laplerfc'
        !
        !--------------------------------------------------------------------------------
        !
        cell => laplacian%cell
        !
        IF (axis < 1 .OR. axis > 3) CALL env_errore(sub_name, 'Wrong value of axis', 1)
        !
        chargeanalytic = erfcvolume(dim, axis, width, spread, cell)
        !
        scale = charge / chargeanalytic / sqrtpi / spread
        ! scaling factor, take into account rescaling of generated density
        ! to obtain the correct integrated total charge
        !
        ALLOCATE (lapllocal(cell%nnr))
        lapllocal = 0.D0
        !
        DO ir = 1, cell%ir_end
            !
            CALL cell%get_min_distance(ir, dim, axis, pos, r, r2, physical)
            ! compute minimum distance using minimum image convention
            !
            IF (.NOT. physical) CYCLE
            !
            dist = SQRT(r2)
            !
            arg = (dist - width) / spread
            !
            !----------------------------------------------------------------------------
            ! Compute laplacian of error function
            !
            SELECT CASE (dim)
                !
            CASE (0)
                !
                IF (dist > tol) &
                    lapllocal(ir) = -EXP(-arg**2) * (1.D0 / dist - arg / spread) * 2.D0
                !
            CASE (1)
                !
                IF (dist > tol) &
                    lapllocal(ir) = -EXP(-arg**2) * (1.D0 / dist - 2.D0 * arg / spread)
                !
            CASE (2)
                lapllocal(ir) = EXP(-arg**2) * arg / spread * 2.D0
                !
            END SELECT
            !
        END DO
        !
        laplacian%of_r = laplacian%of_r + lapllocal * scale
        DEALLOCATE (lapllocal)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE generate_laplerfc
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE generate_hesserfc(dim, axis, charge, width, spread, pos, hessian)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: dim, axis
        REAL(DP), INTENT(IN) :: charge, width, spread
        REAL(DP), INTENT(IN) :: pos(3)
        !
        CLASS(environ_hessian), INTENT(INOUT) :: hessian
        !
        LOGICAL :: physical
        INTEGER :: ir, ip, jp
        REAL(DP) :: scale, r2, dist, arg, tmp, chargeanalytic
        REAL(DP) :: r(3)
        REAL(DP), ALLOCATABLE :: hesslocal(:, :, :)
        !
        CLASS(environ_cell), POINTER :: cell
        !
        CHARACTER(LEN=80) :: sub_name = 'generate_hesserfc'
        !
        !--------------------------------------------------------------------------------
        !
        cell => hessian%cell
        !
        IF (axis < 1 .OR. axis > 3) CALL env_errore(sub_name, 'Wrong value of axis', 1)
        !
        chargeanalytic = erfcvolume(dim, axis, width, spread, cell)
        scale = charge / chargeanalytic / sqrtpi / spread
        !
        ALLOCATE (hesslocal(3, 3, cell%nnr))
        hesslocal = 0.D0
        !
        DO ir = 1, cell%ir_end
            !
            CALL cell%get_min_distance(ir, dim, axis, pos, r, r2, physical)
            ! compute minimum distance using minimum image convention
            !
            IF (.NOT. physical) CYCLE
            !
            dist = SQRT(r2)
            arg = (dist - width) / spread
            !
            !----------------------------------------------------------------------------
            ! Compute hessian of error function
            !
            IF (dist > tol) THEN
                !
                DO ip = 1, 3
                    !
                    DO jp = 1, 3
                        tmp = -r(ip) * r(jp) * (1.D0 / dist + 2.D0 * arg / spread)
                        !
                        IF (ip == jp) tmp = tmp + dist
                        !
                        hesslocal(ip, jp, ir) = -EXP(-arg**2) * tmp / dist**2
                    END DO
                    !
                END DO
                !
            END IF
            !
        END DO
        !
        hessian%of_r = hessian%of_r + hesslocal * scale
        DEALLOCATE (hesslocal)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE generate_hesserfc
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE generate_axis(cell, icor, pos, axis)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_cell), INTENT(IN) :: cell
        INTEGER, INTENT(IN) :: icor
        REAL(DP), INTENT(IN) :: pos(3)
        !
        REAL(DP), INTENT(OUT) :: axis(cell%nnr)
        !
        LOGICAL :: physical
        INTEGER :: ir
        REAL(DP) :: r(3), r2
        !
        !--------------------------------------------------------------------------------
        !
        DO ir = 1, cell%ir_end
            !
            CALL cell%get_min_distance(ir, 0, 0, pos, r, r2, physical)
            ! compute minimum distance using minimum image convention
            !
            IF (.NOT. physical) CYCLE
            !
            axis(ir) = -r(icor)
            !
        END DO
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE generate_axis
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE generate_distance(cell, pos, distance)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_cell), INTENT(IN) :: cell
        REAL(DP), INTENT(IN) :: pos(3)
        !
        REAL(DP), INTENT(OUT) :: distance(3, cell%nnr)
        !
        LOGICAL :: physical
        INTEGER :: ir
        REAL(DP) :: r(3), r2
        !
        !--------------------------------------------------------------------------------
        !
        DO ir = 1, cell%ir_end
            !
            CALL cell%get_min_distance(ir, 0, 0, pos, r, r2, physical)
            ! compute minimum distance using minimum image convention
            !
            IF (.NOT. physical) CYCLE
            !
            distance(:, ir) = -r
            !
        END DO
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE generate_distance
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    FUNCTION erfcvolume(dim, axis, width, spread, cell)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP) :: erfcvolume
        !
        INTEGER, INTENT(IN) :: dim, axis
        REAL(DP), INTENT(IN) :: width, spread
        CLASS(environ_cell), INTENT(IN) :: cell
        !
        REAL(DP) :: f1 = 0.0_DP, f2 = 0.0_DP
        REAL(DP) :: t, invt
        !
        CHARACTER(LEN=80) :: fun_name = 'erfcvolume'
        !
        !--------------------------------------------------------------------------------
        !
        IF (spread < tol .OR. width < tol) &
            CALL env_errore(fun_name, 'Wrong parameters of erfc function', 1)
        !
        t = spread / width
        invt = width / spread
        f1 = (1.D0 + environ_erf(invt)) / 2.D0 ! f1 is close to one  for t-->0
        f2 = EXP(-(invt)**2) / 2.D0 / sqrtpi ! f2 is close to zero for t-->0
        !
        SELECT CASE (dim)
            !
        CASE (0)
            !
            !----------------------------------------------------------------------------
            ! Zero-dimensional erfc, volume is approx the one of the
            ! sphere of radius=width
            !
            erfcvolume = fpi / 3.D0 * width**3 * &
                         ((1.D0 + 1.5D0 * t**2) * f1 + (1.D0 + t**2) * t * f2)
            !
        CASE (1)
            !
            !----------------------------------------------------------------------------
            ! One-dimensional erfc, volume is approx the one of the
            ! cylinder of radius=width and length=alat*at(axis,axis)
            !
            erfcvolume = pi * width**2 * cell%at(axis, axis) * &
                         ((1.D0 + 0.5D0 * t**2) * f1 + t * f2)
            !
        CASE (2)
            !
            !----------------------------------------------------------------------------
            ! Two-dimensional erfc, volume is exactly the one of the
            ! box, does not depend on spread
            !
            erfcvolume = 2.D0 * width * cell%omega / cell%at(axis, axis)
            !
        END SELECT
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END FUNCTION erfcvolume
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE generate_functions
!----------------------------------------------------------------------------------------
