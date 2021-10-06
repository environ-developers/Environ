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
!          Edan Bainglass     (Department of Physics, UNT)
!
!----------------------------------------------------------------------------------------
!>
!!
!----------------------------------------------------------------------------------------
MODULE class_function_erfc
    !------------------------------------------------------------------------------------
    !
    USE class_io, ONLY: io
    USE env_mp, ONLY: env_mp_sum
    !
    USE environ_param, ONLY: DP, sqrtpi, pi, fpi
    !
    USE class_cell
    USE class_density
    USE class_function
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
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    TYPE, EXTENDS(environ_function), PUBLIC :: environ_function_erfc
        !--------------------------------------------------------------------------------
        !
        !--------------------------------------------------------------------------------
    CONTAINS
        !--------------------------------------------------------------------------------
        !
        PROCEDURE :: density => density_of_function_erfc
        PROCEDURE :: gradient => gradient_of_function_erfc
        PROCEDURE :: laplacian => laplacian_of_function_erfc
        PROCEDURE :: hessian => hessian_of_function_erfc
        !
        PROCEDURE, PRIVATE :: get_charge, erfcvolume
        !
        !--------------------------------------------------------------------------------
    END TYPE environ_function_erfc
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
CONTAINS
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
    SUBROUTINE density_of_function_erfc(this, density, zero)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_function_erfc), TARGET, INTENT(IN) :: this
        LOGICAL, INTENT(IN), OPTIONAL :: zero
        !
        TYPE(environ_density), TARGET, INTENT(INOUT) :: density
        !
        LOGICAL :: physical
        INTEGER :: ir
        REAL(DP) :: r(3), r2, scale, dist, arg, chargeanalytic, integral, local_charge
        REAL(DP), ALLOCATABLE :: local(:)
        !
        CHARACTER(LEN=80) :: sub_name = 'density_of_function_erfc'
        !
        !--------------------------------------------------------------------------------
        !
        IF (this%f_type == 5) THEN
            density%of_r = this%volume
        ELSE IF (PRESENT(zero)) THEN
            IF (zero) density%of_r = 0.D0
        END IF
        !
        ASSOCIATE (cell => density%cell, &
                   pos => this%pos, &
                   spread => this%spread, &
                   width => this%width, &
                   dim => this%dim, &
                   axis => this%axis)
            !
            IF (axis < 1 .OR. axis > 3) CALL io%error(sub_name, 'Wrong value of axis', 1)
            !
            local_charge = this%get_charge(cell)
            chargeanalytic = this%erfcvolume(cell)
            !
            scale = local_charge / chargeanalytic * 0.5D0
            ! scaling factor, take into account rescaling of generated density
            ! to obtain the correct integrated total charge
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
            !----------------------------------------------------------------------------
            ! Check integral of function is consistent with analytic one
            !
            integral = SUM(local) * cell%omega / DBLE(cell%ntot) * 0.5D0
            !
            CALL env_mp_sum(integral, cell%dfft%comm)
            !
            IF (ABS(integral - chargeanalytic) / chargeanalytic > 1.D-4) &
                CALL io%warning('wrong integral of erfc function')
            !
            !----------------------------------------------------------------------------
            !
            density%of_r = density%of_r + scale * local
            ! rescale generated function to obtain the requested integral
            !
            DEALLOCATE (local)
            !
        END ASSOCIATE
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE density_of_function_erfc
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE gradient_of_function_erfc(this, gradient, zero)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_function_erfc), TARGET, INTENT(IN) :: this
        LOGICAL, INTENT(IN), OPTIONAL :: zero
        !
        TYPE(environ_gradient), TARGET, INTENT(INOUT) :: gradient
        !
        LOGICAL :: physical
        INTEGER :: ir
        REAL(DP) :: r(3), r2, scale, dist, arg, chargeanalytic, local_charge
        REAL(DP), ALLOCATABLE :: gradlocal(:, :)
        !
        CHARACTER(LEN=80) :: sub_name = 'gradient_of_function_erfc'
        !
        !--------------------------------------------------------------------------------
        !
        IF (PRESENT(zero)) THEN
            IF (zero) gradient%of_r = 0.D0
        END IF
        !
        ASSOCIATE (cell => gradient%cell, &
                   pos => this%pos, &
                   spread => this%spread, &
                   width => this%width, &
                   dim => this%dim, &
                   axis => this%axis)
            !
            IF (axis < 1 .OR. axis > 3) CALL io%error(sub_name, 'Wrong value of axis', 1)
            !
            local_charge = this%get_charge(cell)
            chargeanalytic = this%erfcvolume(cell)
            !
            scale = local_charge / chargeanalytic / sqrtpi / spread
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
                IF (dist > func_tol) gradlocal(:, ir) = -EXP(-arg**2) * r / dist
                ! compute gradient of error function
                !
            END DO
            !
            gradient%of_r = gradient%of_r + gradlocal * scale
            DEALLOCATE (gradlocal)
            !
        END ASSOCIATE
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE gradient_of_function_erfc
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE laplacian_of_function_erfc(this, laplacian, zero)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_function_erfc), TARGET, INTENT(IN) :: this
        LOGICAL, INTENT(IN), OPTIONAL :: zero
        !
        TYPE(environ_density), TARGET, INTENT(INOUT) :: laplacian
        !
        LOGICAL :: physical
        INTEGER :: ir
        REAL(DP) :: r(3), r2, scale, dist, arg, chargeanalytic, local_charge
        REAL(DP), ALLOCATABLE :: lapllocal(:)
        !
        CHARACTER(LEN=80) :: sub_name = 'laplacian_of_function_erfc'
        !
        !--------------------------------------------------------------------------------
        !
        IF (PRESENT(zero)) THEN
            IF (zero) laplacian%of_r = 0.D0
        END IF
        !
        ASSOCIATE (cell => laplacian%cell, &
                   pos => this%pos, &
                   spread => this%spread, &
                   width => this%width, &
                   dim => this%dim, &
                   axis => this%axis)
            !
            IF (axis < 1 .OR. axis > 3) CALL io%error(sub_name, 'Wrong value of axis', 1)
            !
            local_charge = this%get_charge(cell)
            chargeanalytic = this%erfcvolume(cell)
            !
            scale = local_charge / chargeanalytic / sqrtpi / spread
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
                arg = (dist - width) / spread
                !
                !------------------------------------------------------------------------
                ! Compute laplacian of error function
                !
                SELECT CASE (dim)
                    !
                CASE (0)
                    !
                    IF (dist > func_tol) &
                        lapllocal(ir) = -EXP(-arg**2) * &
                                        (1.D0 / dist - arg / spread) * 2.D0
                    !
                CASE (1)
                    !
                    IF (dist > func_tol) &
                        lapllocal(ir) = -EXP(-arg**2) * &
                                        (1.D0 / dist - 2.D0 * arg / spread)
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
        END ASSOCIATE
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE laplacian_of_function_erfc
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE hessian_of_function_erfc(this, hessian, zero)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_function_erfc), TARGET, INTENT(IN) :: this
        LOGICAL, INTENT(IN), OPTIONAL :: zero
        !
        TYPE(environ_hessian), TARGET, INTENT(INOUT) :: hessian
        !
        LOGICAL :: physical
        INTEGER :: ir, ip, jp
        REAL(DP) :: r(3), r2, scale, dist, arg, tmp, chargeanalytic, local_charge
        REAL(DP), ALLOCATABLE :: hesslocal(:, :, :)
        !
        CHARACTER(LEN=80) :: sub_name = 'hessian_of_function_erfc'
        !
        !--------------------------------------------------------------------------------
        !
        IF (PRESENT(zero)) THEN
            IF (zero) hessian%of_r = 0.D0
        END IF
        !
        ASSOCIATE (cell => hessian%cell, &
                   pos => this%pos, &
                   spread => this%spread, &
                   width => this%width, &
                   dim => this%dim, &
                   axis => this%axis)
            !
            IF (axis < 1 .OR. axis > 3) CALL io%error(sub_name, 'Wrong value of axis', 1)
            !
            local_charge = this%get_charge(cell)
            chargeanalytic = this%erfcvolume(cell)
            !
            scale = local_charge / chargeanalytic / sqrtpi / spread
            ! scaling factor, take into account rescaling of generated density
            ! to obtain the correct integrated total charge
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
                !------------------------------------------------------------------------
                ! Compute hessian of error function
                !
                IF (dist > func_tol) THEN
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
        END ASSOCIATE
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE hessian_of_function_erfc
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
    REAL(DP) FUNCTION erfcvolume(this, cell)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_function_erfc), INTENT(IN) :: this
        TYPE(environ_cell), INTENT(IN) :: cell
        !
        REAL(DP) :: f1 = 0.0_DP, f2 = 0.0_DP
        REAL(DP) :: t, invt
        !
        CHARACTER(LEN=80) :: fun_name = 'erfcvolume'
        !
        !--------------------------------------------------------------------------------
        !
        ASSOCIATE (spread => this%spread, &
                   width => this%width, &
                   dim => this%dim, &
                   axis => this%axis)
            !
            IF (spread < func_tol .OR. width < func_tol) &
                CALL io%error(fun_name, 'Wrong parameters of erfc function', 1)
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
                !------------------------------------------------------------------------
                ! Zero-dimensional erfc, volume is approx the one of the
                ! sphere of radius=width
                !
                erfcvolume = fpi / 3.D0 * width**3 * &
                             ((1.D0 + 1.5D0 * t**2) * f1 + (1.D0 + t**2) * t * f2)
                !
            CASE (1)
                !
                !------------------------------------------------------------------------
                ! One-dimensional erfc, volume is approx the one of the
                ! cylinder of radius=width and length=alat*at(axis,axis)
                !
                erfcvolume = pi * width**2 * cell%at(axis, axis) * &
                             ((1.D0 + 0.5D0 * t**2) * f1 + t * f2)
                !
            CASE (2)
                !
                !------------------------------------------------------------------------
                ! Two-dimensional erfc, volume is exactly the one of the
                ! box, does not depend on spread
                !
                erfcvolume = 2.D0 * width * cell%omega / cell%at(axis, axis)
                !
            END SELECT
            !
        END ASSOCIATE
        !
        !--------------------------------------------------------------------------------
    END FUNCTION erfcvolume
    !------------------------------------------------------------------------------------
    !>
    !! type 2 - CHARGE * NORMALIZED_ERFC_HALF(X) ! integrates to charge
    !!
    !! type 4 - HARGE * NORMALIZED_ERFC_HALF(X) * VOLUME_NORMALIZED_ERFC_HALF
    !!          goes from charge to 0
    !!
    !! type 5 - CHARGE * (1 - NORMALIZED_ERFC_HALF(x) * VOLUME_NORMALIZED_ERFC_HALF)
    !!          goes from 0 to charge
    !!
    !------------------------------------------------------------------------------------
    REAL(DP) FUNCTION get_charge(this, cell)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_function_erfc), INTENT(IN) :: this
        TYPE(environ_cell), INTENT(IN) :: cell
        !
        !--------------------------------------------------------------------------------
        !
        ASSOCIATE (f_type => this%f_type, &
                   charge => this%volume)
            !
            SELECT CASE (f_type)
                !
            CASE (2)
                get_charge = charge
                !
            CASE (4)
                get_charge = this%erfcvolume(cell) * charge
                !
            CASE (5)
                get_charge = -this%erfcvolume(cell) * charge
                !
            END SELECT
            !
        END ASSOCIATE
        !
        !--------------------------------------------------------------------------------
    END FUNCTION get_charge
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE class_function_erfc
!----------------------------------------------------------------------------------------
