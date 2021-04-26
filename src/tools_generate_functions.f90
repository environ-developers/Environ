! Copyright (C) 2018 ENVIRON (www.quantum-environment.org)
! Copyright (C) 2006-2010 Quantum ESPRESSO group
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
!
!----------------------------------------------------------------------------------------
!>
!! Module to generate functions on the real space dense grid
!!
!----------------------------------------------------------------------------------------
MODULE tools_generate_functions
    !------------------------------------------------------------------------------------
    !
    USE modules_constants, ONLY: DP, sqrtpi, pi, fpi
    !
    USE cell_types, ONLY: environ_cell
    USE representation_types, ONLY: environ_density, environ_gradient, environ_hessian
    !
    USE tools_cell, ONLY: ir2ijk, ir2r, displacement, minimum_image
    !
    USE modules_erf, ONLY: environ_erfc, environ_erf
    !
    USE mp, ONLY: mp_sum
    !
    !------------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    REAL(DP), PARAMETER :: tol = 1.D-10
    REAL(DP), PARAMETER :: exp_tol = 4.D1
    !
    !------------------------------------------------------------------------------------
    !
    PRIVATE :: tol, exp_tol
    !
    !------------------------------------------------------------------------------------
CONTAINS
    !------------------------------------------------------------------------------------
    !>
    !! #TODO unused
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE planar_average(cell, nnr, naxis, axis, shift, reverse, f, f1d)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_cell), INTENT(IN) :: cell
        INTEGER, INTENT(IN) :: nnr, naxis, axis, shift
        LOGICAL, INTENT(IN) :: reverse
        !
        REAL(DP), INTENT(INOUT) :: f(nnr)
        REAL(DP), INTENT(INOUT) :: f1d(naxis)
        !
        INTEGER :: i, j, k, ir
        INTEGER :: idx, narea
        LOGICAL :: physical
        !
        !--------------------------------------------------------------------------------
        !
        narea = cell%ntot / naxis
        !
        IF (reverse) THEN
            f = 0.D0
        ELSE
            f1d = 0.D0
        END IF
        !
        DO ir = 1, cell%ir_end
            !
            CALL ir2ijk(cell, ir, i, j, k, physical) ! three dimensional indexes
            !
            IF (.NOT. physical) CYCLE ! do not include points outside the physical range
            !
            SELECT CASE (axis)
            CASE (1)
                idx = i
            CASE (2)
                idx = j
            CASE (3)
                idx = k
            END SELECT
            !
            idx = idx + 1 + shift
            !
            IF (idx > naxis) THEN
                idx = idx - naxis
            ELSE IF (idx <= 0) THEN
                idx = idx + naxis
            END IF
            !
            IF (reverse) THEN
                f(ir) = f1d(idx)
            ELSE
                f1d(idx) = f1d(idx) + f(ir)
            END IF
            !
        END DO
        !
        IF (.NOT. reverse) THEN
            !
            CALL mp_sum(f1d(:), cell%dfft%comm)
            !
            f1d = f1d / DBLE(narea)
        END IF
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE planar_average
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
        TYPE(environ_density), INTENT(INOUT) :: density
        !
        LOGICAL :: physical
        INTEGER :: ir
        REAL(DP) :: scale, spr2, length
        REAL(DP) :: r(3), r2
        REAL(DP), ALLOCATABLE :: local(:)
        !
        TYPE(environ_cell), POINTER :: cell
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
            CALL errore(sub_name, 'Wrong spread for Gaussian function', 1)
        !
        IF (axis < 1 .OR. axis > 3) &
            CALL errore(sub_name, 'Wrong axis in generate_gaussian', 1)
        !
        SELECT CASE (dim)
        CASE (0)
            scale = charge / (sqrtpi * spread)**3
        CASE (1)
            length = ABS(cell%at(axis, axis) * cell%alat)
            scale = charge / length / (sqrtpi * spread)**2
        CASE (2)
            length = ABS(cell%at(axis, axis) * cell%alat)
            scale = charge * length / cell%omega / (sqrtpi * spread)
        CASE default
            CALL errore(sub_name, 'Wrong value of dim', 1)
        END SELECT
        !
        spr2 = (spread / cell%alat)**2
        !
        ALLOCATE (local(cell%nnr))
        local = 0.D0
        !
        DO ir = 1, cell%ir_end
            !
            CALL ir2r(cell, ir, r, physical) ! position in real space grid
            !
            IF (.NOT. physical) CYCLE ! do not include points outside the physical range
            !
            CALL displacement(dim, axis, pos, r, r) ! displacement from origin
            !
            CALL minimum_image(cell, r, r2) ! minimum image convention
            !
            !----------------------------------------------------------------------------
            ! Compute Gaussian function
            !
            r2 = r2 / spr2
            !
            IF (r2 <= exp_tol) local(ir) = EXP(-r2)
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
        TYPE(environ_gradient), INTENT(INOUT) :: gradient
        !
        LOGICAL :: physical
        INTEGER :: ir
        REAL(DP) :: scale, spr2, length
        REAL(DP) :: r(3), r2
        REAL(DP), ALLOCATABLE :: gradlocal(:, :)
        !
        TYPE(environ_cell), POINTER :: cell
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
            CALL errore(sub_name, 'Wrong spread for Gaussian function', 1)
        !
        IF (axis < 1 .OR. axis > 3) &
            CALL errore(sub_name, 'Wrong value of axis', 1)
        !
        SELECT CASE (dim)
        CASE (0)
            scale = charge / (sqrtpi * spread)**3
        CASE (1)
            length = ABS(cell%at(axis, axis) * cell%alat)
            scale = charge / length / (sqrtpi * spread)**2
        CASE (2)
            length = ABS(cell%at(axis, axis) * cell%alat)
            scale = charge * length / cell%omega / (sqrtpi * spread)
        CASE default
            CALL errore(sub_name, 'Wrong value of dim', 1)
        END SELECT
        !
        scale = scale * 2.D0 / spread**2 * cell%alat
        !
        spr2 = (spread / cell%alat)**2
        !
        ALLOCATE (gradlocal(3, cell%nnr))
        gradlocal = 0.D0
        !
        DO ir = 1, cell%ir_end
            !
            CALL ir2r(cell, ir, r, physical) ! position in real space grid
            !
            IF (.NOT. physical) CYCLE ! do not include points outside the physical range
            !
            CALL displacement(dim, axis, pos, r, r) ! displacement from origin
            !
            CALL minimum_image(cell, r, r2) ! minimum image convention
            !
            !----------------------------------------------------------------------------
            ! Compute gradient of Gaussian function
            !
            r2 = r2 / spr2
            !
            IF (r2 <= exp_tol) gradlocal(:, ir) = EXP(-r2) * r(:)
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
        TYPE(environ_density), INTENT(INOUT) :: density
        !
        LOGICAL :: physical
        INTEGER :: ir
        REAL(DP) :: r2, dist, arg
        REAL(DP) :: r(3)
        REAL(DP), ALLOCATABLE :: local(:)
        !
        TYPE(environ_cell), POINTER :: cell
        !
        CHARACTER(LEN=80) :: sub_name = 'generate_exponential'
        !
        !--------------------------------------------------------------------------------
        !
        cell => density%cell
        !
        IF (axis < 1 .OR. axis > 3) &
            CALL errore(sub_name, 'Wrong value of axis', 1)
        !
        ALLOCATE (local(cell%nnr))
        local = 0.D0
        !
        DO ir = 1, cell%ir_end
            !
            CALL ir2r(cell, ir, r, physical) ! position in real space grid
            !
            IF (.NOT. physical) CYCLE ! do not include points outside the physical range
            !
            CALL displacement(dim, axis, pos, r, r) ! displacement from origin
            !
            CALL minimum_image(cell, r, r2) ! minimum image convention
            !
            !----------------------------------------------------------------------------
            ! Compute exponentially decaying function
            !
            dist = SQRT(r2) * cell%alat
            arg = (dist - width) / spread
            !
            IF (ABS(arg) <= exp_tol) local(ir) = EXP(-arg)
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
        TYPE(environ_gradient), INTENT(INOUT) :: gradient
        !
        LOGICAL :: physical
        INTEGER :: ir
        REAL(DP) :: r2, dist, arg
        REAL(DP) :: r(3)
        REAL(DP), ALLOCATABLE :: gradlocal(:, :)
        !
        TYPE(environ_cell), POINTER :: cell
        !
        CHARACTER(LEN=80) :: sub_name = 'generate_gradexponential'
        !
        !--------------------------------------------------------------------------------
        !
        cell => gradient%cell
        !
        IF (axis < 1 .OR. axis > 3) &
            CALL errore(sub_name, 'Wrong value of axis', 1)
        !
        ALLOCATE (gradlocal(3, cell%nnr))
        gradlocal = 0.D0
        !
        DO ir = 1, cell%ir_end
            !
            CALL ir2r(cell, ir, r, physical) ! position in real space grid
            !
            IF (.NOT. physical) CYCLE ! do not include points outside the physical range
            !
            CALL displacement(dim, axis, pos, r, r) ! displacement from origin
            !
            CALL minimum_image(cell, r, r2) ! minimum image convention
            !
            !----------------------------------------------------------------------------
            ! Compute exponentially decaying function
            !
            dist = SQRT(r2) * cell%alat
            arg = (dist - width) / spread
            !
            IF (r2 > tol .AND. ABS(arg) <= exp_tol) &
                gradlocal(:, ir) = r(:) / SQRT(r2) / spread * EXP(-arg)
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
        TYPE(environ_density), INTENT(INOUT) :: density
        !
        LOGICAL :: physical
        INTEGER :: ir, ir_end, i
        REAL(DP) :: scale, r2, dist, arg, chargeanalytic, chargelocal
        REAL(DP) :: r(3)
        REAL(DP), ALLOCATABLE :: local(:)
        !
        TYPE(environ_cell), POINTER :: cell
        !
        CHARACTER(LEN=80) :: sub_name = 'generate_erfc'
        !
        !--------------------------------------------------------------------------------
        !
        cell => density%cell
        !
        IF (axis < 1 .OR. axis > 3) &
            CALL errore(sub_name, 'Wrong value of axis', 1)
        !
        chargeanalytic = erfcvolume(dim, axis, width, spread, cell)
        scale = charge / chargeanalytic * 0.5D0
        !
        ALLOCATE (local(cell%nnr))
        local = 0.D0
        !
        DO ir = 1, cell%ir_end
            !
            CALL ir2r(cell, ir, r, physical) ! position in real space grid
            !
            IF (.NOT. physical) CYCLE ! do not include points outside the physical range
            !
            CALL displacement(dim, axis, pos, r, r) ! displacement from origin
            !
            CALL minimum_image(cell, r, r2) ! minimum image convention
            !
            !----------------------------------------------------------------------------
            ! Compute error function
            !
            dist = SQRT(r2) * cell%alat
            arg = (dist - width) / spread
            local(ir) = environ_erfc(arg)
            !
            !----------------------------------------------------------------------------
            !
        END DO
        !
        !--------------------------------------------------------------------------------
        ! Check integral of function is consistent with analytic one
        !
        chargelocal = SUM(local) * cell%omega / DBLE(cell%ntot) * 0.5D0
        !
        CALL mp_sum(chargelocal, cell%dfft%comm)
        !
        IF (ABS(chargelocal - chargeanalytic) / chargeanalytic > 1.D-4) &
            CALL infomsg(sub_name, 'WARNING: wrong integral of erfc function')
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
        TYPE(environ_gradient), INTENT(INOUT) :: gradient
        !
        LOGICAL :: physical
        INTEGER :: ir
        REAL(DP) :: scale, r2, dist, arg, chargeanalytic
        REAL(DP) :: r(3)
        REAL(DP), ALLOCATABLE :: gradlocal(:, :)
        !
        TYPE(environ_cell), POINTER :: cell
        !
        CHARACTER(LEN=80) :: sub_name = 'generate_graderfc'
        !
        !--------------------------------------------------------------------------------
        !
        cell => gradient%cell
        !
        IF (axis < 1 .OR. axis > 3) &
            CALL errore(sub_name, 'Wrong value of axis', 1)
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
            CALL ir2r(cell, ir, r, physical) ! position in real space grid
            !
            IF (.NOT. physical) CYCLE ! do not include points outside the physical range
            !
            CALL displacement(dim, axis, pos, r, r) ! displacement from origin
            !
            CALL minimum_image(cell, r, r2) ! minimum image convention
            !
            !----------------------------------------------------------------------------
            ! Compute gradient of error function
            !
            r = r * cell%alat
            dist = SQRT(r2) * cell%alat
            arg = (dist - width) / spread
            !
            IF (dist > tol) gradlocal(:, ir) = EXP(-arg**2) * r(:) / dist
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
        TYPE(environ_density), INTENT(INOUT) :: laplacian
        !
        LOGICAL :: physical
        INTEGER :: ir
        REAL(DP) :: scale, r2, dist, arg, chargeanalytic
        REAL(DP) :: r(3)
        REAL(DP), ALLOCATABLE :: lapllocal(:)
        !
        TYPE(environ_cell), POINTER :: cell
        !
        CHARACTER(LEN=80) :: sub_name = 'generate_laplerfc'
        !
        !--------------------------------------------------------------------------------
        !
        cell => laplacian%cell
        !
        IF (axis < 1 .OR. axis > 3) &
            CALL errore(sub_name, 'Wrong value of axis', 1)
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
            CALL ir2r(cell, ir, r, physical) ! position in real space grid
            !
            IF (.NOT. physical) CYCLE ! do not include points outside the physical range
            !
            CALL displacement(dim, axis, pos, r, r) ! displacement from origin
            !
            CALL minimum_image(cell, r, r2) ! minimum image convention
            !
            !----------------------------------------------------------------------------
            ! Compute laplacian of error function
            !
            dist = SQRT(r2) * cell%alat
            arg = (dist - width) / spread
            !
            SELECT CASE (dim)
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
        TYPE(environ_hessian), INTENT(INOUT) :: hessian
        !
        LOGICAL :: physical
        INTEGER :: ir, ip, jp
        REAL(DP) :: scale, r2, dist, arg, tmp, chargeanalytic
        REAL(DP) :: r(3)
        REAL(DP), ALLOCATABLE :: hesslocal(:, :, :)
        !
        TYPE(environ_cell), POINTER :: cell
        !
        CHARACTER(LEN=80) :: sub_name = 'generate_hesserfc'
        !
        !--------------------------------------------------------------------------------
        !
        cell => hessian%cell
        !
        IF (axis < 1 .OR. axis > 3) &
            CALL errore(sub_name, 'Wrong value of axis', 1)
        !
        chargeanalytic = erfcvolume(dim, axis, width, spread, cell)
        scale = charge / chargeanalytic / sqrtpi / spread
        !
        ALLOCATE (hesslocal(3, 3, cell%nnr))
        hesslocal = 0.D0
        !
        DO ir = 1, cell%ir_end
            !
            CALL ir2r(cell, ir, r, physical) ! position in real space grid
            !
            IF (.NOT. physical) CYCLE ! do not include points outside the physical range
            !
            CALL displacement(dim, axis, pos, r, r) ! displacement from origin
            !
            CALL minimum_image(cell, r, r2) ! minimum image convention
            !
            !----------------------------------------------------------------------------
            ! Compute hessian of error function
            !
            r = r * cell%alat
            dist = SQRT(r2) * cell%alat
            arg = (dist - width) / spread
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
    !! #TODO field-aware
    !!
    !------------------------------------------------------------------------------------
    ! SUBROUTINE generate_deriverfc(nnr, dim, axis, charge, width, spread, pos, drho)
    !     !--------------------------------------------------------------------------------
    !     !
    !     USE modules_constants, ONLY: DP, sqrtpi
    !     USE fft_base, ONLY: dfftp
    !     USE mp, ONLY: mp_sum
    !     USE mp_bands, ONLY: me_bgrp, intra_bgrp_comm
    !     !
    !     IMPLICIT NONE
    !     !
    !     REAL(DP), PARAMETER :: tol = 1.D-10
    !     !
    !     INTEGER, INTENT(IN) :: nnr, dim, axis
    !     REAL(DP), INTENT(IN) :: charge, width, spread
    !     REAL(DP), INTENT(IN) :: pos(3)
    !     TYPE(environ_density), TARGET, INTENT(INOUT) :: drho
    !     !
    !     INTEGER :: i, j, j0, k, k0, ir, ir_end, ip
    !     INTEGER :: idx, idx0, ntot
    !     !
    !     REAL(DP) :: inv_nr1, inv_nr2, inv_nr3
    !     REAL(DP) :: scale, dist, arg, chargeanalytic, chargelocal
    !     REAL(DP) :: r(3), s(3)
    !     REAL(DP), ALLOCATABLE :: drholocal(:)
    !     REAL(DP), EXTERNAL :: qe_erfc
    !     REAL(DP) :: alat, at(3, 3), bg(3, 3), omega
    !     !
    !     TYPE(environ_cell) :: cell
    !     alat = drho%cell%alat
    !     at = drho%cell%at
    !     bg = drho%cell%bg
    !     omega = drho%cell%omega
    !     !
    !     IF (dfftp%nr1 == 0 .OR. dfftp%nr2 == 0 .OR. dfftp%nr3 == 0) THEN
    !         WRITE (6, *) 'ERROR: wrong grid dimension', dfftp%nr1, dfftp%nr2, dfftp%nr3
    !         !
    !         STOP
    !         !
    !     END IF
    !     !
    !     inv_nr1 = 1.D0 / DBLE(dfftp%nr1)
    !     inv_nr2 = 1.D0 / DBLE(dfftp%nr2)
    !     inv_nr3 = 1.D0 / DBLE(dfftp%nr3)
    !     !
    !     ntot = dfftp%nr1 * dfftp%nr2 * dfftp%nr3
    !     !
    !     IF (axis < 1 .OR. axis > 3) &
    !         WRITE (6, *) 'WARNING: wrong axis in generate_gaussian'
    !     !
    !     chargeanalytic = erfcvolume(dim, axis, width, spread, cell)
    !     !
    !     !--------------------------------------------------------------------------------
    !     ! scaling factor, take into account rescaling of generated density
    !     ! to obtain the correct integrated total charge
    !     !
    !     scale = charge / chargeanalytic / sqrtpi / spread
    !     !
    !     ALLOCATE (drholocal(nnr))
    !     drholocal = 0.D0
    !     chargelocal = 0.D0
    !     !
    !     ! BACKWARD COMPATIBILITY
    !     ! Compatible with QE-5.X QE-6.1.X
    !     ! idx0 = dfftp%nr1x * dfftp%nr2x * dfftp%ipp(me_bgrp + 1)
    !     ! ir_end = dfftp%nr1x * dfftp%nr2x * dfftp%npl
    !     ! Compatible with QE-6.2, QE-6.2.1 and QE-GIT
    !     #if defined (__MPI)
    !     j0 = dfftp%my_i0r2p; k0 = dfftp%my_i0r3p
    !     ir_end = MIN(nnr, dfftp%nr1x * dfftp%my_nr2p * dfftp%my_nr3p)
    !     #else
    !     j0 = 0; k0 = 0
    !     ir_end = nnr
    !     #endif
    !     ! END BACKWARD COMPATIBILITY
    !     !
    !     DO ir = 1, ir_end
    !         !
    !         !----------------------------------------------------------------------------
    !         ! three dimensional indexes
    !         !
    !         ! BACKWARD COMPATIBILITY
    !         ! Compatible with QE-5.X QE-6.1.X
    !         ! idx = idx0 + ir - 1
    !         ! k = idx / (dfftp%nr1x * dfftp%nr2x)
    !         ! idx = idx - (dfftp%nr1x * dfftp%nr2x) * k
    !         ! j = idx / dfftp%nr1x
    !         ! idx = idx - dfftp%nr1x * j
    !         ! i = idx
    !         ! Compatible with QE-6.2, QE-6.2.1 and QE-GIT
    !         idx = ir - 1
    !         k = idx / (dfftp%nr1x * dfftp%my_nr2p)
    !         idx = idx - (dfftp%nr1x * dfftp%my_nr2p) * k
    !         k = k + k0
    !         j = idx / dfftp%nr1x
    !         idx = idx - dfftp%nr1x * j
    !         j = j + j0
    !         i = idx
    !         ! END BACKWARD COMPATIBILITY
    !         !
    !         !----------------------------------------------------------------------------
    !         !
    !         IF (i >= dfftp%nr1 .OR. j >= dfftp%nr2 .OR. k >= dfftp%nr3) CYCLE
    !         ! do not include points outside the physical range
    !         !
    !         !----------------------------------------------------------------------------
    !         !
    !         DO ip = 1, 3
    !             !
    !             r(ip) = DBLE(i) * inv_nr1 * at(ip, 1) + &
    !                     DBLE(j) * inv_nr2 * at(ip, 2) + &
    !                     DBLE(k) * inv_nr3 * at(ip, 3)
    !             !
    !         END DO
    !         !
    !         r(:) = r(:) - pos(:)
    !         !
    !         !----------------------------------------------------------------------------
    !         ! possibly 2D or 1D erfc
    !         !
    !         IF (dim == 1) THEN
    !             r(axis) = 0.D0
    !         ELSE IF (dim == 2) THEN
    !             !
    !             DO i = 1, 3
    !                 IF (i /= axis) r(i) = 0.D0
    !             END DO
    !             !
    !         END IF
    !         !
    !         !----------------------------------------------------------------------------
    !         ! minimum image convention
    !         !
    !         s(:) = MATMUL(r(:), bg(:, :))
    !         s(:) = s(:) - ANINT(s(:))
    !         r(:) = MATMUL(at(:, :), s(:))
    !         r = r * alat
    !         !
    !         dist = SQRT(SUM(r * r))
    !         arg = (dist - width) / spread
    !         !
    !         IF (dist > tol) drholocal(ir) = -EXP(-arg**2)
    !         chargelocal = chargelocal + qe_erfc(arg)
    !         !
    !     END DO
    !     !
    !     !--------------------------------------------------------------------------------
    !     ! double check that the integral of the generated charge corresponds to
    !     ! what is expected
    !     !
    !     CALL mp_sum(chargelocal, intra_bgrp_comm)
    !     chargelocal = chargelocal * omega / DBLE(ntot) * 0.5D0
    !     !
    !     IF (ABS(chargelocal - chargeanalytic) / chargeanalytic > 1.D-4) &
    !         WRITE (6, *) &
    !         'WARNING: significant discrepancy between &
    !         &the numerical and the !expected erfc charge'
    !     !
    !     drholocal = drholocal * scale
    !     !
    !     drho%of_r = drho%of_r + drholocal
    !     DEALLOCATE (drholocal)
    !     !
    !     RETURN
    !     !
    !     !--------------------------------------------------------------------------------
    ! END SUBROUTINE generate_deriverfc
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE generate_axis(cell, icor, pos, axis)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_cell), INTENT(IN) :: cell
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
            CALL ir2r(cell, ir, r, physical) ! position in real space grid
            !
            IF (.NOT. physical) CYCLE ! do not include points outside the physical range
            !
            r = r - pos ! displacement from origin
            !
            CALL minimum_image(cell, r, r2) ! minimum image convention
            !
            axis(ir) = r(icor)
            !
        END DO
        !
        axis = axis * cell%alat
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
        TYPE(environ_cell), INTENT(IN) :: cell
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
            CALL ir2r(cell, ir, r, physical) ! position in real space grid
            !
            IF (.NOT. physical) CYCLE ! do not include points outside the physical range
            !
            r = r - pos ! displacement from origin
            !
            CALL minimum_image(cell, r, r2) ! minimum image convention
            !
            distance(:, ir) = r(:)
            !
        END DO
        !
        distance = distance * cell%alat
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
        TYPE(environ_cell), INTENT(IN) :: cell
        !
        REAL(DP) :: f1 = 0.0_DP, f2 = 0.0_DP
        REAL(DP) :: t, invt
        !
        CHARACTER(LEN=80) :: fun_name = 'erfcvolume'
        !
        !--------------------------------------------------------------------------------
        !
        IF (spread < tol .OR. width < tol) &
            CALL errore(fun_name, 'Wrong parameters of erfc function', 1)
        !
        t = spread / width
        invt = width / spread
        f1 = (1.D0 + environ_erf(invt)) / 2.D0 ! f1 is close to one  for t-->0
        f2 = EXP(-(invt)**2) / 2.D0 / sqrtpi ! f2 is close to zero for t-->0
        !
        SELECT CASE (dim)
        CASE (0)
            !
            !----------------------------------------------------------------------------
            ! Zero-dimensional erfc, volume is approx the one of the
            ! sphere of radius=width
            !
            erfcvolume = fpi / 3.D0 * width**3 * &
                         ((1.D0 + 1.5D0 * t**2) * f1 + (1.D0 + t**2) * t * f2)
        CASE (1)
            !
            !----------------------------------------------------------------------------
            ! One-dimensional erfc, volume is approx the one of the
            ! cylinder of radius=width and length=alat*at(axis,axis)
            !
            erfcvolume = pi * width**2 * cell%at(axis, axis) * cell%alat * &
                         ((1.D0 + 0.5D0 * t**2) * f1 + t * f2)
        CASE (2)
            !
            !----------------------------------------------------------------------------
            ! Two-dimensional erfc, volume is exactly the one of the
            ! box, does not depend on spread
            !
            erfcvolume = 2.D0 * width * cell%omega / cell%at(axis, axis) / cell%alat
        END SELECT
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END FUNCTION erfcvolume
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE tools_generate_functions
!----------------------------------------------------------------------------------------

