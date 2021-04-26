!----------------------------------------------------------------------------------------
!>
!!
!----------------------------------------------------------------------------------------
MODULE tools_math
    !------------------------------------------------------------------------------------
    !
    USE modules_constants, ONLY: DP
    !
    USE representation_types, ONLY: environ_density, environ_gradient, environ_hessian
    USE cell_types, ONLY: environ_cell
    !
    USE tools_cell, ONLY: ir2r, displacement, minimum_image
    !
    USE mp, ONLY: mp_sum
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
            CALL errore(fun_name, 'operation on fields with inconsistent domains', 1)
        !
        ir_end => density1%cell%ir_end
        scalar_product = DOT_PRODUCT(density1%of_r(1:ir_end), density2%of_r(1:ir_end))
        !
        CALL mp_sum(scalar_product, density1%cell%dfft%comm)
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
            CALL errore(sub_name, 'Missmatch in domain of input gradients', 1)
        !
        IF (.NOT. ASSOCIATED(gradA%cell, dens%cell)) &
            CALL errore(sub_name, 'Missmatch in domain of input and output', 1)
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
            CALL errore(sub_name, 'Missmatch in domain of input hessian/gradients', 1)
        !
        IF (.NOT. ASSOCIATED(gradin%cell, gradout%cell)) &
            CALL errore(sub_name, 'Missmatch in domain of input and output', 1)
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
            CALL errore(sub_name, 'Missmatch in domain of input vectors', 1)
        !
        ir_end => density%cell%ir_end
        !
        DO ipol = 1, 3
            !
            scalar_product = DOT_PRODUCT(gradient%of_r(ipol, 1:ir_end), &
                                         density%of_r(1:ir_end))
            !
            CALL mp_sum(scalar_product, density%cell%dfft%comm)
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
        CALL mp_sum(integral, density%cell%dfft%comm)
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
        CALL mp_sum(euclidean_norm, density%cell%dfft%comm)
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
        CALL mp_sum(quadratic_mean, density%cell%dfft%comm)
        !
        quadratic_mean = SQRT(quadratic_mean / density%cell%ntot)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END FUNCTION quadratic_mean_environ_density
    !------------------------------------------------------------------------------------
    !>
    !! #TODO unused
    !!
    !------------------------------------------------------------------------------------
    FUNCTION quadratic_mean_environ_density_old(density) RESULT(quadratic_mean)
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
        CALL mp_sum(quadratic_mean, density%cell%dfft%comm)
        !
        quadratic_mean = SQRT(quadratic_mean) / density%cell%ntot
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END FUNCTION quadratic_mean_environ_density_old
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
        CALL mp_sum(monopole, cell%dfft%comm)
        !
        CALL mp_sum(dipole, cell%dfft%comm)
        !
        CALL mp_sum(quadrupole, cell%dfft%comm)
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
    !
    !------------------------------------------------------------------------------------
END MODULE tools_math
!----------------------------------------------------------------------------------------
