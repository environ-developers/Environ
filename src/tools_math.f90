!----------------------------------------------------------------------------------------
!
! Copyright (C) 2018-2021 ENVIRON (www.quantum-environ.org)
! Copyright (C) 2002-2009 Quantum ESPRESSO (www.quantum-espresso.org)
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
!!
!----------------------------------------------------------------------------------------
MODULE tools_math
    !------------------------------------------------------------------------------------
    !
    USE env_mp, ONLY: env_mp_sum
    !
    USE environ_param, ONLY: DP
    !
    USE types_representation, ONLY: environ_density, environ_gradient, environ_hessian
    USE types_cell, ONLY: environ_cell
    !
    USE tools_cell, ONLY: ir2r, displacement, minimum_image
    !
    !------------------------------------------------------------------------------------
    !
    PRIVATE
    !
    PUBLIC :: scalar_product_environ_density, scalar_product_environ_gradient, &
              scalar_product_environ_hessian, scalar_product_environ_gradient_density, &
              integrate_environ_density, euclidean_norm_environ_density, &
              quadratic_mean_environ_density, multipoles_environ_density, &
              environ_erf, environ_erfc
    !
    !------------------------------------------------------------------------------------
CONTAINS
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    FUNCTION scalar_product_environ_density(density1, density2) RESULT(scalar_product)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_density), INTENT(IN) :: density1, density2
        !
        INTEGER, POINTER :: ir_end
        REAL(DP) :: scalar_product
        !
        CHARACTER(LEN=80) :: fun_name = 'scalar_product_environ_density'
        !
        !--------------------------------------------------------------------------------
        !
        IF (.NOT. ASSOCIATED(density1%cell, density2%cell)) &
            CALL env_errore(fun_name, 'Operation on fields with inconsistent domains', 1)
        !
        ir_end => density1%cell%ir_end
        scalar_product = DOT_PRODUCT(density1%of_r(1:ir_end), density2%of_r(1:ir_end))
        !
        CALL env_mp_sum(scalar_product, density1%cell%dfft%comm)
        !
        scalar_product = scalar_product * density1%cell%domega
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END FUNCTION scalar_product_environ_density
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE scalar_product_environ_gradient(gradA, gradB, dens)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_gradient), INTENT(IN) :: gradA, gradB
        !
        TYPE(environ_density), INTENT(INOUT) :: dens
        !
        INTEGER :: ir
        !
        CHARACTER(LEN=80) :: sub_name = 'scalar_product_environ_gradient'
        !
        !--------------------------------------------------------------------------------
        !
        dens%of_r = 0.D0
        !
        IF (.NOT. ASSOCIATED(gradA%cell, gradB%cell)) &
            CALL env_errore(sub_name, 'Mismatch in domain of input gradients', 1)
        !
        IF (.NOT. ASSOCIATED(gradA%cell, dens%cell)) &
            CALL env_errore(sub_name, 'Mismatch in domain of input and output', 1)
        !
        DO ir = 1, dens%cell%ir_end
            dens%of_r(ir) = SUM(gradA%of_r(:, ir) * gradB%of_r(:, ir))
        END DO
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE scalar_product_environ_gradient
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE scalar_product_environ_hessian(hess, gradin, gradout)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_hessian), INTENT(IN) :: hess
        TYPE(environ_gradient), INTENT(IN) :: gradin
        !
        TYPE(environ_gradient), INTENT(INOUT) :: gradout
        !
        INTEGER :: ir, ipol
        !
        CHARACTER(LEN=80) :: sub_name = 'scalar_product_environ_hessian'
        !
        !--------------------------------------------------------------------------------
        !
        gradout%of_r = 0.D0
        !
        IF (.NOT. ASSOCIATED(gradin%cell, hess%cell)) &
            CALL env_errore(sub_name, 'Mismatch in domain of input hessian/gradients', 1)
        !
        IF (.NOT. ASSOCIATED(gradin%cell, gradout%cell)) &
            CALL env_errore(sub_name, 'Mismatch in domain of input and output', 1)
        !
        DO ir = 1, hess%cell%ir_end
            !
            DO ipol = 1, 3
                gradout%of_r(ipol, ir) = SUM(hess%of_r(:, ipol, ir) * gradin%of_r(:, ir))
            END DO
            !
        END DO
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE scalar_product_environ_hessian
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    FUNCTION scalar_product_environ_gradient_density(gradient, density) RESULT(res)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_gradient), INTENT(IN) :: gradient
        TYPE(environ_density), INTENT(IN) :: density
        !
        REAL(DP) :: res(3)
        !
        INTEGER, POINTER :: ir_end
        !
        INTEGER :: ipol
        REAL(DP) :: scalar_product
        !
        CHARACTER(LEN=80) :: sub_name = 'scalar_product_environ_gradient_density'
        !
        !--------------------------------------------------------------------------------
        !
        res = 0.D0
        !
        IF (.NOT. ASSOCIATED(gradient%cell, density%cell)) &
            CALL env_errore(sub_name, 'Mismatch in domain of input vectors', 1)
        !
        ir_end => density%cell%ir_end
        !
        DO ipol = 1, 3
            !
            scalar_product = DOT_PRODUCT(gradient%of_r(ipol, 1:ir_end), &
                                         density%of_r(1:ir_end))
            !
            CALL env_mp_sum(scalar_product, density%cell%dfft%comm)
            !
            res(ipol) = scalar_product * density%cell%domega
        END DO
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END FUNCTION scalar_product_environ_gradient_density
    !------------------------------------------------------------------------------------
    !>
    !! #TODO unused
    !!
    !------------------------------------------------------------------------------------
    FUNCTION dipole_of_origin(density, origin) RESULT(dipole)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_density), INTENT(IN) :: density
        !
        REAL(DP), DIMENSION(3) :: origin, dipole
        !
        !--------------------------------------------------------------------------------
        !
        dipole = density%dipole + density%charge * (density%cell%origin - origin)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END FUNCTION dipole_of_origin
    !------------------------------------------------------------------------------------
    !>
    !! #TODO unused
    !!
    !------------------------------------------------------------------------------------
    FUNCTION quadrupole_of_origin(density, origin) RESULT(quadrupole)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_density), INTENT(IN) :: density
        !
        REAL(DP), DIMENSION(3) :: origin, quadrupole
        !
        !--------------------------------------------------------------------------------
        !
        quadrupole = density%quadrupole + &
                     density%charge * (density%cell%origin - origin)**2 + &
                     2.D0 * density%dipole * (density%cell%origin - origin)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END FUNCTION quadrupole_of_origin
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    FUNCTION integrate_environ_density(density) RESULT(integral)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_density), INTENT(IN) :: density
        !
        REAL(DP) :: integral
        !
        !--------------------------------------------------------------------------------
        !
        integral = SUM(density%of_r(1:density%cell%ir_end))
        !
        CALL env_mp_sum(integral, density%cell%dfft%comm)
        !
        integral = integral * density%cell%domega
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END FUNCTION integrate_environ_density
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    FUNCTION euclidean_norm_environ_density(density) RESULT(euclidean_norm)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_density), INTENT(IN) :: density
        !
        INTEGER, POINTER :: ir_end
        !
        REAL(DP) :: euclidean_norm
        !
        !--------------------------------------------------------------------------------
        !
        ir_end => density%cell%ir_end
        euclidean_norm = DOT_PRODUCT(density%of_r(1:ir_end), density%of_r(1:ir_end))
        !
        CALL env_mp_sum(euclidean_norm, density%cell%dfft%comm)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END FUNCTION euclidean_norm_environ_density
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    FUNCTION quadratic_mean_environ_density(density) RESULT(quadratic_mean)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_density), INTENT(IN) :: density
        !
        INTEGER, POINTER :: ir_end
        !
        REAL(DP) :: quadratic_mean
        !
        !--------------------------------------------------------------------------------
        !
        ir_end => density%cell%ir_end
        quadratic_mean = DOT_PRODUCT(density%of_r(1:ir_end), density%of_r(1:ir_end))
        !
        CALL env_mp_sum(quadratic_mean, density%cell%dfft%comm)
        !
        quadratic_mean = SQRT(quadratic_mean / density%cell%ntot)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END FUNCTION quadratic_mean_environ_density
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE multipoles_environ_density(density, origin, monopole, dipole, quadrupole)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_density), TARGET, INTENT(IN) :: density
        REAL(DP), INTENT(IN) :: origin(3)
        !
        REAL(DP), INTENT(OUT) :: monopole
        REAL(DP), DIMENSION(3), INTENT(OUT) :: dipole, quadrupole
        !
        TYPE(environ_cell), POINTER :: cell
        !
        LOGICAL :: physical
        INTEGER :: ir
        REAL(DP) :: r(3), rhoir, r2
        INTEGER :: dim, axis
        !
        !--------------------------------------------------------------------------------
        !
        cell => density%cell
        !
        monopole = 0.D0
        dipole = 0.D0
        quadrupole = 0.D0
        !
        dim = 0
        axis = 3
        !
        DO ir = 1, cell%ir_end
            !
            CALL ir2r(cell, ir, r, physical) ! position in real space grid
            !
            IF (.NOT. physical) CYCLE ! do not include points outside the physical range
            !
            CALL displacement(dim, axis, r, origin, r) ! displacement from origin
            !
            CALL minimum_image(cell, r, r2) ! minimum image convention
            !
            rhoir = density%of_r(ir)
            !
            !----------------------------------------------------------------------------
            ! Multipoles
            !
            monopole = monopole + rhoir
            dipole = dipole + rhoir * r
            quadrupole = quadrupole + rhoir * r**2
            !
            !----------------------------------------------------------------------------
            !
        END DO
        !
        CALL env_mp_sum(monopole, cell%dfft%comm)
        !
        CALL env_mp_sum(dipole, cell%dfft%comm)
        !
        CALL env_mp_sum(quadrupole, cell%dfft%comm)
        !
        monopole = monopole * cell%domega
        dipole = dipole * cell%domega * cell%alat
        quadrupole = quadrupole * cell%domega * cell%alat**2
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE multipoles_environ_density
    !------------------------------------------------------------------------------------
    !>
    !! Error function - computed from the rational approximations of
    !! W. J. Cody, Math. Comp. 22 (1969), pages 631-637.
    !!
    !! for abs(x) le 0.47 erf is calculated directly
    !! for abs(x) gt 0.47 erf is calculated via erf(x)=1-erfc(x)
    !!
    !------------------------------------------------------------------------------------
    FUNCTION environ_erf(x)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP), INTENT(IN) :: x
        !
        REAL(DP) :: x2, p1(4), q1(4)
        REAL(DP) :: environ_erf
        !
        DATA p1/2.426679552305318E2_DP, 2.197926161829415E1_DP, &
            6.996383488619136_DP, -3.560984370181538E-2_DP/
        !
        DATA q1/2.150588758698612E2_DP, 9.116490540451490E1_DP, &
            1.508279763040779E1_DP, 1.000000000000000_DP/
        !
        !--------------------------------------------------------------------------------
        !
        IF (ABS(x) > 6.0_DP) THEN
            !
            environ_erf = SIGN(1.0_DP, x)
            ! erf(6) = 1 - 10^(-17) cannot be distinguished from 1
            !
        ELSE
            !
            IF (ABS(x) <= 0.47_DP) THEN
                x2 = x**2
                !
                environ_erf = x * (p1(1) + x2 * (p1(2) + x2 * (p1(3) + x2 * p1(4)))) / &
                              (q1(1) + x2 * (q1(2) + x2 * (q1(3) + x2 * q1(4))))
                !
            ELSE
                environ_erf = 1.0_DP - environ_erfc(x)
            END IF
            !
        END IF
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END FUNCTION environ_erf
    !------------------------------------------------------------------------------------
    !>
    !! erfc(x) = 1-erf(x) - See comments in erf
    !!
    !------------------------------------------------------------------------------------
    FUNCTION environ_erfc(x)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP), INTENT(IN) :: x
        !
        REAL(DP) :: environ_erfc
        REAL(DP) :: ax, x2, xm2, p2(8), q2(8), p3(5), q3(5), pim1
        !
        DATA p2/3.004592610201616E2_DP, 4.519189537118719E2_DP, &
            3.393208167343437E2_DP, 1.529892850469404E2_DP, &
            4.316222722205674E1_DP, 7.211758250883094_DP, &
            5.641955174789740E-1_DP, -1.368648573827167E-7_DP/
        !
        DATA q2/3.004592609569833E2_DP, 7.909509253278980E2_DP, &
            9.313540948506096E2_DP, 6.389802644656312E2_DP, &
            2.775854447439876E2_DP, 7.700015293522947E1_DP, &
            1.278272731962942E1_DP, 1.000000000000000_DP/
        !
        DATA p3/-2.996107077035422E-3_DP, -4.947309106232507E-2_DP, &
            -2.269565935396869E-1_DP, -2.786613086096478E-1_DP, &
            -2.231924597341847E-2_DP/
        !
        DATA q3/1.062092305284679E-2_DP, 1.913089261078298E-1_DP, &
            1.051675107067932_DP, 1.987332018171353_DP, &
            1.000000000000000_DP/
        !
        DATA pim1/0.56418958354775629_DP/ ! ( pim1= sqrt(1/pi) )
        !
        !--------------------------------------------------------------------------------
        !
        ax = ABS(x)
        !
        IF (ax > 26.0_DP) THEN
            environ_erfc = 0.0_DP ! erfc(26.0) = 10^(-296); erfc(9.0) = 10^(-37)
        ELSE IF (ax > 4.0_DP) THEN
            x2 = x**2
            xm2 = (1.0_DP / ax)**2
            !
            environ_erfc = 1.0_DP / ax * &
                           EXP(-x2) * &
                           (pim1 + xm2 * ( &
                            p3(1) + xm2 * ( &
                            p3(2) + xm2 * ( &
                            p3(3) + xm2 * ( &
                            p3(4) + xm2 * p3(5) &
                            )))) / &
                            (q3(1) + xm2 * ( &
                             q3(2) + xm2 * ( &
                             q3(3) + xm2 * ( &
                             q3(4) + xm2 * q3(5) &
                             )))))
            !
        ELSE IF (ax > 0.47_DP) THEN
            x2 = x**2
            !
            environ_erfc = EXP(-x2) * &
                           (p2(1) + ax * ( &
                            p2(2) + ax * ( &
                            p2(3) + ax * ( &
                            p2(4) + ax * ( &
                            p2(5) + ax * ( &
                            p2(6) + ax * ( &
                            p2(7) + ax * p2(8) &
                            ))))))) / &
                           (q2(1) + ax * ( &
                            q2(2) + ax * ( &
                            q2(3) + ax * ( &
                            q2(4) + ax * ( &
                            q2(5) + ax * ( &
                            q2(6) + ax * ( &
                            q2(7) + ax * q2(8) &
                            )))))))
            !
        ELSE
            environ_erfc = 1.0_DP - environ_erf(ax)
        END IF
        !
        IF (x < 0.0_DP) environ_erfc = 2.0_DP - environ_erfc
        ! erf(-x) = -erf(x) => erfc(-x) = 2 - erfc(x)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END FUNCTION environ_erfc
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE tools_math
!----------------------------------------------------------------------------------------
