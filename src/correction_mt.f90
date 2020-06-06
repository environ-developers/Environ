!
! Copyright (C) 2001-2006 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE correction_mt
  !
  ! ... The variables needed to the Martyna-Tuckerman method for isolated
  !     systems
  !
  USE modules_constants, ONLY: DP, pi, tpi, fpi, e2
  USE cell_types, ONLY : ir2r, minimum_image
  USE core_types
  !
  IMPLICIT NONE
  !
  SAVE
  !
  PRIVATE
  !
  PUBLIC :: update_mt_correction, calc_vmt, calc_gradvmt, calc_fmt
  !
CONTAINS
!----------------------------------------------------------------------------
  SUBROUTINE update_mt_correction( fft )
!----------------------------------------------------------------------------
    !
    USE modules_erf,       ONLY : environ_erfc
    USE fft_interfaces,    ONLY : fwfft
    !
    TYPE( fft_core ), TARGET, INTENT(INOUT) :: fft
    !
    LOGICAL :: physical
    INTEGER :: ir, ig
    REAL(DP) :: r(3), rws, upperbound, ecutrho
    COMPLEX (DP), ALLOCATABLE :: aux(:)
    !
    REAL( DP ), POINTER :: alpha, beta, omega, tpiba2
    REAL( DP ), DIMENSION(:), POINTER :: mt_corr
    TYPE( fft_type_descriptor ), POINTER :: dfft
    TYPE( environ_cell ), POINTER :: cell
    !
    alpha => fft % alpha
    beta => fft % beta
    mt_corr => fft % mt_corr
    !
    cell => fft % cell
    dfft => fft % dfft
    omega => fft % cell % omega
    tpiba2 => fft % cell % tpiba2
    !
    ecutrho = fft % gcutm * tpiba2
    !
    ! choose alpha in order to have convergence in the sum over G
    ! upperbound is a safe upper bound for the error in the sum over G
    !
    alpha = 2.9d0
    upperbound = 1._dp
    DO WHILE ( upperbound > 1.e-7_dp)
       alpha = alpha - 0.1_dp
       IF (alpha<=0._dp) CALL errore('init_mt_correction','optimal alpha not found',1)
       upperbound = e2 * SQRT(2.d0 * alpha / tpi) * &
            environ_erfc( SQRT( ecutrho / 4.d0 / alpha) )
    END DO
    beta = 0.5_dp/alpha ! 1._dp/alpha
    !
    ALLOCATE (aux(dfft%nnr))
    aux = (0._dp,0._dp)
    !
    DO ir = 1, cell%ir_end
       !
       ! ... three dimensional indexes
       !
       CALL ir2r( cell, ir, r, physical )
       !
       ! ... do not include points outside the physical range
       !
       IF ( .NOT. physical ) CYCLE
       !
       ! ... minimum image convention
       !
       CALL minimum_image( cell, r, rws )
       !
       aux(ir) = smooth_coulomb_r( alpha, SQRT(rws)*cell%alat )
       !
    END DO
    !
    CALL fwfft ('Rho', aux, dfft)
    !
!$omp parallel do
    DO ig = 1, fft%ngm
       mt_corr(ig) = cell%omega * REAL(aux(dfft%nl(ig))) - &
            smooth_coulomb_g( alpha, beta, cell%tpiba2*fft%gg(ig))
    ENDDO
!$omp end parallel do
    mt_corr(:) =  mt_corr(:) * exp(-cell%tpiba2*fft%gg(:)*beta/4._dp)**2
    !
    DEALLOCATE (aux)
    !
    RETURN
    !
!----------------------------------------------------------------------------
  END SUBROUTINE update_mt_correction
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
  SUBROUTINE calc_vmt( fft, rho, v )
!----------------------------------------------------------------------------
    !
    TYPE( fft_core ), INTENT(IN) :: fft
    COMPLEX(DP), INTENT(IN)  :: rho(fft%ngm)
    COMPLEX(DP), INTENT(OUT) :: v(fft%ngm)
    !
    INTEGER :: ig
    !
    v(:) = (0._dp,0._dp)
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
!----------------------------------------------------------------------------
  END SUBROUTINE calc_vmt
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
  SUBROUTINE calc_gradvmt( ipol, fft, rho, v )
!----------------------------------------------------------------------------
    !
    INTEGER, INTENT(IN) :: ipol
    TYPE( fft_core ), INTENT(IN) :: fft
    COMPLEX(DP), INTENT(IN)  :: rho(fft%ngm)
    COMPLEX(DP), INTENT(OUT) :: v(fft%ngm)
    !
    INTEGER :: ig
    REAL( DP ) :: fac
    !
    v(:) = (0._dp,0._dp)
    !
!$omp parallel do private(fac)
    DO ig = fft%gstart, fft%ngm
       fac = fft%g(ipol,ig) * fft%cell%tpiba
       v(ig) = fft%mt_corr(ig) * CMPLX(-AIMAG(rho(ig)),REAL(rho(ig)),kind=dp)*fac
    END DO
!$omp end parallel do
    !
    v = e2 * v
    !
    RETURN
    !
!----------------------------------------------------------------------------
  END SUBROUTINE calc_gradvmt
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
  SUBROUTINE calc_fmt( fft, rho, ions, force )
!----------------------------------------------------------------------------
    !
    USE environ_types
    !
    TYPE( fft_core ), INTENT(IN) :: fft
    TYPE( environ_ions ), INTENT(IN) :: ions
    COMPLEX(DP), INTENT(IN) :: rho(fft%ngm)
    REAL(DP), INTENT(OUT) :: force(3,ions%number)
    !
    INTEGER :: iat, ig
    REAL (DP) :: arg
    !
    force = 0._dp
    DO iat = 1, ions%number
       DO ig = 1, fft % ngm
          arg = tpi * SUM ( fft%g(:,ig)*ions%tau(:,iat) )
          force( :, iat ) = force( :, iat ) + fft%g( :, ig ) * &
               ( SIN(arg)*REAL(rho(ig)) + COS(arg)*AIMAG(rho(ig))) * &
               fft%mt_corr(ig)
       END DO
       force( :, iat ) = force( :, iat ) * ions%iontype(ions%ityp(iat))%zv * fft%cell%tpiba
    END DO
    !
    force = e2 * force
    !
    RETURN
    !
!----------------------------------------------------------------------------
  END SUBROUTINE calc_fmt
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
  REAL(DP) FUNCTION smooth_coulomb_r(alpha,r)
 !----------------------------------------------------------------------------
    !
    USE modules_erf, ONLY : environ_erf
    !
    REAL(DP), INTENT(IN) :: alpha, r
    !
    IF (r>1.e-6_dp) THEN
       smooth_coulomb_r = environ_erf(SQRT(alpha)*r)/r
    ELSE
       smooth_coulomb_r = 2._dp/SQRT(pi) * SQRT(alpha)
    END IF
    !
!----------------------------------------------------------------------------
  END FUNCTION smooth_coulomb_r
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
  REAL(DP) FUNCTION smooth_coulomb_g(alpha,beta,q2)
!----------------------------------------------------------------------------
    !
    REAL(DP), INTENT(IN) :: alpha, beta, q2
    !
    IF (q2>1.e-6_dp) THEN
       smooth_coulomb_g = fpi * EXP(-q2/4._dp/alpha)/q2
    ELSE
       smooth_coulomb_g = - 1._dp * fpi * (1._dp/4._dp/alpha + 2._dp*beta/4._dp)
    END IF
    !
!----------------------------------------------------------------------------
  END FUNCTION smooth_coulomb_g
!----------------------------------------------------------------------------
END MODULE correction_mt
