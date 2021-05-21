!
! Copyright (C) 2001-2006 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------------------
!>
!! The variables needed to the Martyna-Tuckerman method for isolated systems
!!
!----------------------------------------------------------------------------------------
MODULE correction_mt
    !------------------------------------------------------------------------------------
    !
    USE env_fft_types, ONLY: env_fft_type_descriptor
    USE env_fft_interfaces, ONLY: env_fwfft
    !
    USE modules_constants, ONLY: DP, pi, tpi, fpi, e2
    !
    USE core_types, ONLY: fft_core
    USE cell_types, ONLY: environ_cell
    USE physical_types, ONLY: environ_ions
    !
    USE tools_cell, ONLY: ir2r, minimum_image
    USE tools_math, ONLY: environ_erf, environ_erfc
    !
    !------------------------------------------------------------------------------------
    !
    PRIVATE
    !
    PUBLIC :: update_mt_correction, calc_vmt, calc_gradvmt, calc_fmt
    !
    !------------------------------------------------------------------------------------
CONTAINS
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE update_mt_correction(fft)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(fft_core), TARGET, INTENT(INOUT) :: fft
        !
        LOGICAL :: physical
        INTEGER :: ir, ig
        REAL(DP) :: r(3), rws, upperbound, ecutrho
        COMPLEX(DP), ALLOCATABLE :: aux(:)
        !
        REAL(DP), POINTER :: alpha, beta, omega, tpiba2
        REAL(DP), POINTER :: mt_corr(:)
        TYPE(env_fft_type_descriptor), POINTER :: dfft
        TYPE(environ_cell), POINTER :: cell
        !
        CHARACTER(LEN=80) :: sub_name = 'update_mt_correction'
        !
        !--------------------------------------------------------------------------------
        !
        alpha => fft%alpha
        beta => fft%beta
        mt_corr => fft%mt_corr
        !
        cell => fft%cell
        dfft => fft%cell%dfft
        omega => fft%cell%omega
        tpiba2 => fft%cell%tpiba2
        !
        ecutrho = fft%gcutm * tpiba2
        !
        !--------------------------------------------------------------------------------
        ! choose alpha in order to have convergence in the sum over G
        ! upperbound is a safe upper bound for the error in the sum over G
        !
        alpha = 2.9D0
        upperbound = 1._DP
        !
        DO WHILE (upperbound > 1.E-7_DP)
            alpha = alpha - 0.1_DP
            !
            IF (alpha <= 0._DP) &
                CALL env_errore(sub_name, 'optimal alpha not found', 1)
            !
            upperbound = e2 * SQRT(2.D0 * alpha / tpi) * &
                         environ_erfc(SQRT(ecutrho / 4.D0 / alpha))
            !
        END DO
        !
        beta = 0.5_DP / alpha ! 1._dp/alpha
        !
        ALLOCATE (aux(dfft%nnr))
        aux = (0._DP, 0._DP)
        !
        DO ir = 1, cell%ir_end
            !
            CALL ir2r(cell, ir, r, physical) ! position in real space grid
            !
            IF (.NOT. physical) CYCLE ! do not include points outside the physical range
            !
            r = r - cell%origin ! the function should be centered on the origin
            !
            CALL minimum_image(cell, r, rws) ! minimum image convention
            !
            aux(ir) = smooth_coulomb_r(alpha, SQRT(rws) * cell%alat)
            !
        END DO
        !
        CALL env_fwfft(aux, dfft)
        !
!$omp parallel do
        DO ig = 1, fft%ngm
            !
            mt_corr(ig) = cell%omega * REAL(aux(dfft%nl(ig))) - &
                          smooth_coulomb_g(alpha, beta, cell%tpiba2 * fft%gg(ig))
            !
        END DO
!$omp end parallel do
        !
        mt_corr(:) = mt_corr(:) * EXP(-cell%tpiba2 * fft%gg(:) * beta / 4._DP)**2
        !
        DEALLOCATE (aux)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE update_mt_correction
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_vmt(fft, rho, v)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(fft_core), INTENT(IN) :: fft
        COMPLEX(DP), INTENT(IN) :: rho(fft%ngm)
        !
        COMPLEX(DP), INTENT(OUT) :: v(fft%ngm)
        !
        INTEGER :: ig
        !
        !--------------------------------------------------------------------------------
        !
        v(:) = (0._DP, 0._DP)
        !
!$omp parallel do
        DO ig = 1, fft%ngm
            v(ig) = fft%mt_corr(ig) * rho(ig)
        END DO
!$omp end parallel do
        !
        v = e2 * v
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_vmt
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_gradvmt(ipol, fft, rho, v)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: ipol
        TYPE(fft_core), INTENT(IN) :: fft
        COMPLEX(DP), INTENT(IN) :: rho(fft%ngm)
        !
        COMPLEX(DP), INTENT(OUT) :: v(fft%ngm)
        !
        INTEGER :: ig
        REAL(DP) :: fac
        !
        !--------------------------------------------------------------------------------
        !
        v(:) = (0._DP, 0._DP)
        !
!$omp parallel do private(fac)
        DO ig = fft%gstart, fft%ngm
            fac = fft%g(ipol, ig) * fft%cell%tpiba
            !
            v(ig) = fft%mt_corr(ig) * &
                    CMPLX(-AIMAG(rho(ig)), REAL(rho(ig)), kind=dp) * fac
            !
        END DO
!$omp end parallel do
        !
        v = e2 * v
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_gradvmt
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_fmt(fft, rho, ions, force)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(fft_core), INTENT(IN) :: fft
        TYPE(environ_ions), INTENT(IN) :: ions
        COMPLEX(DP), INTENT(IN) :: rho(fft%ngm)
        !
        REAL(DP), INTENT(OUT) :: force(3, ions%number)
        !
        INTEGER :: iat, ig
        REAL(DP) :: arg
        !
        !--------------------------------------------------------------------------------
        !
        force = 0._DP
        !
        DO iat = 1, ions%number
            !
            DO ig = 1, fft%ngm
                arg = tpi * SUM(fft%g(:, ig) * ions%tau(:, iat))
                !
                force(:, iat) = &
                    force(:, iat) + fft%g(:, ig) * &
                    (SIN(arg) * REAL(rho(ig)) + COS(arg) * AIMAG(rho(ig))) * &
                    fft%mt_corr(ig)
                !
            END DO
            !
            force(:, iat) = &
                force(:, iat) * ions%iontype(ions%ityp(iat))%zv * fft%cell%tpiba
            !
        END DO
        !
        force = e2 * force
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_fmt
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    REAL(DP) FUNCTION smooth_coulomb_r(alpha, r)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP), INTENT(IN) :: alpha, r
        !
        !--------------------------------------------------------------------------------
        !
        IF (r > 1.E-6_DP) THEN
            smooth_coulomb_r = environ_erf(SQRT(alpha) * r) / r
        ELSE
            smooth_coulomb_r = 2._DP / SQRT(pi) * SQRT(alpha)
        END IF
        !
        !--------------------------------------------------------------------------------
    END FUNCTION smooth_coulomb_r
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    REAL(DP) FUNCTION smooth_coulomb_g(alpha, beta, q2)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP), INTENT(IN) :: alpha, beta, q2
        !
        !--------------------------------------------------------------------------------
        !
        IF (q2 > 1.E-6_DP) THEN
            smooth_coulomb_g = fpi * EXP(-q2 / 4._DP / alpha) / q2
        ELSE
            !
            smooth_coulomb_g = -1._DP * fpi * &
                               (1._DP / 4._DP / alpha + 2._DP * beta / 4._DP)
            !
        END IF
        !
        !--------------------------------------------------------------------------------
    END FUNCTION smooth_coulomb_g
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE correction_mt
!----------------------------------------------------------------------------------------
