!
! Copyright (C) 2001-2006 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#undef TESTING
MODULE env_martyna_tuckerman
  !
  ! ... The variables needed to the Martyna-Tuckerman method for isolated
  !     systems
  !
  USE env_kinds, ONLY: dp
  USE env_constants, ONLY : e2, pi, tpi, fpi
  USE env_ws_base
  !
  IMPLICIT NONE
  !
  TYPE (ws_type) :: ws
  REAL (DP) :: alpha, beta
  REAL (DP), ALLOCATABLE :: wg_corr(:)
  LOGICAL :: wg_corr_is_updated = .FALSE.
  LOGICAL :: do_comp_mt = .FALSE.
  LOGICAL :: gamma_only = .FALSE.
  integer :: gstart = 1
    !
  SAVE

  PRIVATE

  PUBLIC :: env_tag_wg_corr_as_obsolete, do_comp_mt, &
            env_wg_corr_ewald, env_wg_corr_loc, env_wg_corr_h, env_wg_corr_force

  !PUBLIC :: env_tag_wg_corr_as_obsolete, do_comp_mt, &
  !env_wg_corr_ewald, env_wg_corr_loc, env_wg_corr_h, env_wg_corr_force, &
  !env_init_wg_corr

CONTAINS
!----------------------------------------------------------------------------
  SUBROUTINE env_tag_wg_corr_as_obsolete
!----------------------------------------------------------------------------
     wg_corr_is_updated = .FALSE.
  END SUBROUTINE env_tag_wg_corr_as_obsolete
!----------------------------------------------------------------------------
  SUBROUTINE env_wg_corr_h( omega, ngm, rho, v, eh_corr )
!----------------------------------------------------------------------------
  INTEGER, INTENT(IN) :: ngm
  REAL(DP), INTENT(IN) :: omega
  COMPLEX(DP), INTENT(IN)  :: rho(ngm)
  COMPLEX(DP), INTENT(OUT) :: v(ngm)
  REAL(DP), INTENT(OUT) :: eh_corr

  INTEGER :: ig

  IF (.NOT.wg_corr_is_updated) CALL env_init_wg_corr
!
  v(:) = (0._dp,0._dp)

  eh_corr =  0._dp
  DO ig = 1,ngm
     v(ig) = e2 * wg_corr(ig) * rho(ig) 
     eh_corr = eh_corr + ABS(rho(ig))**2 * wg_corr(ig)
  END DO
  iF (gamma_only) v(gstart:ngm) = 0.5_dp * v(gstart:ngm)

  eh_corr = 0.5_dp * e2 * eh_corr * omega

  RETURN
  END SUBROUTINE env_wg_corr_h
!----------------------------------------------------------------------------
  SUBROUTINE env_wg_corr_loc( omega, ntyp, ngm, zv, strf, v )
!----------------------------------------------------------------------------
  INTEGER, INTENT(IN) :: ntyp, ngm
  REAL(DP), INTENT(IN) :: omega, zv(ntyp)
  COMPLEX(DP), INTENT(IN) :: strf(ngm,ntyp)
  COMPLEX(DP), INTENT(OUT) :: v(ngm)
  INTEGER :: ig

  IF (.NOT.wg_corr_is_updated) CALL env_init_wg_corr
!
  do ig=1,ngm
     v(ig) = - e2 * wg_corr(ig) * SUM(zv(1:ntyp)*strf(ig,1:ntyp)) / omega
  end do
  iF (gamma_only) v(gstart:ngm) = 0.5_dp * v(gstart:ngm)

  RETURN
  END SUBROUTINE env_wg_corr_loc
!----------------------------------------------------------------------------
  SUBROUTINE env_wg_corr_force( lnuclei, omega, nat, ntyp, ityp, ngm, g, tau, zv, strf, &
                            rho, force )
!----------------------------------------------------------------------------
  USE env_cell_base, ONLY : tpiba
  USE env_mp_bands,  ONLY : intra_bgrp_comm
  USE env_mp,        ONLY : env_mp_sum
  INTEGER, INTENT(IN) :: nat, ntyp, ityp(nat), ngm
  REAL(DP), INTENT(IN) :: omega, zv(ntyp), tau(3,nat), g(3,ngm)
  COMPLEX(DP), INTENT(IN) :: strf(ngm,ntyp), rho(ngm)
  LOGICAL, INTENT(IN) :: lnuclei
  ! this variable is used in wg_corr_force to select if
  ! corr should be done on rho and nuclei or only on rho
  REAL(DP), INTENT(OUT) :: force(3,nat)
  INTEGER :: ig, na
  REAL (DP) :: arg
  COMPLEX(DP), ALLOCATABLE :: v(:)
  COMPLEX(DP) :: rho_tot
  !
  IF (.NOT.wg_corr_is_updated) CALL env_init_wg_corr
  !
  allocate ( v(ngm) )
  do ig=1,ngm
     rho_tot = rho(ig)
     if(lnuclei) rho_tot = rho_tot - SUM(zv(1:ntyp)*strf(ig,1:ntyp)) / omega
     v(ig) = e2 * wg_corr(ig) * rho_tot
  end do
  force(:,:) = 0._dp
  do na=1,nat
     do ig=1,ngm
        arg = tpi * SUM ( g(:,ig)*tau(:, na) ) 
        force(:,na) = force(:,na) + g(:,ig) * CMPLX(SIN(arg),-COS(ARG), KIND=dp) * v(ig)
     end do
     force(:,na) = - force(:,na) * zv(ityp(na))  * tpiba
  end do
  deallocate ( v )
  !
  call env_mp_sum(  force, intra_bgrp_comm )
  !
  RETURN
  END SUBROUTINE env_wg_corr_force
!----------------------------------------------------------------------------
  SUBROUTINE env_init_wg_corr
!----------------------------------------------------------------------------
  USE env_mp_bands,      ONLY : me_bgrp
  USE env_fft_base,      ONLY : dfftp
  USE env_fft_interfaces,ONLY : env_fwfft, env_invfft
  USE env_control_flags, ONLY : gamma_only_ => gamma_only
  USE env_gvect,         ONLY : ngm, gg, gstart_ => gstart, ecutrho
  USE env_cell_base,     ONLY : at, alat, tpiba2, omega
  USE env_ws_base

  INTEGER :: idx, ir, i,j,k, j0, k0, ig, nt
  REAL(DP) :: r(3), rws, upperbound, rws2
  COMPLEX (DP), ALLOCATABLE :: aux(:)
  REAL(DP), EXTERNAL :: qe_erfc

  IF ( ALLOCATED(wg_corr) ) DEALLOCATE(wg_corr)
  ALLOCATE(wg_corr(ngm))
  !
  ! choose alpha in order to have convergence in the sum over G
  ! upperbound is a safe upper bound for the error in the sum over G
  !
  alpha = 2.9d0
  upperbound = 1._dp
  DO WHILE ( upperbound > 1.e-7_dp) 
     alpha = alpha - 0.1_dp  
     if (alpha<=0._dp) call env_errore('init_wg_corr','optimal alpha not found',1)
     upperbound = e2 * sqrt (2.d0 * alpha / tpi) * &
                       qe_erfc ( sqrt ( ecutrho / 4.d0 / alpha) )
  END DO
  beta = 0.5_dp/alpha ! 1._dp/alpha
  ! write (*,*) " alpha, beta MT = ", alpha, beta
  !
  call env_ws_init(at,ws)
  !
  gstart = gstart_
  gamma_only = gamma_only_
  !
  ALLOCATE (aux(dfftp%nnr))
  aux = (0._dp,0._dp)
  j0 = dfftp%my_i0r2p ; k0 = dfftp%my_i0r3p
  DO ir = 1, dfftp%nr1x*dfftp%my_nr2p*dfftp%my_nr3p
     !
     ! ... three dimensional indexes
     !
     idx = ir -1
     k   = idx / (dfftp%nr1x*dfftp%my_nr2p)
     idx = idx - (dfftp%nr1x*dfftp%my_nr2p)*k
     k   = k + k0
     j   = idx / dfftp%nr1x
     idx = idx - dfftp%nr1x * j
     j   = j + j0
     i   = idx

     ! ... do not include points outside the physical range

     IF ( i >= dfftp%nr1 .OR. j >= dfftp%nr2 .OR. k >= dfftp%nr3 ) CYCLE

     r(:) = ( at(:,1)/dfftp%nr1*i + at(:,2)/dfftp%nr2*j + at(:,3)/dfftp%nr3*k )

     rws = env_ws_dist(r,ws)

     aux(ir) = env_smooth_coulomb_r( rws*alat )

  END DO

  CALL env_fwfft ('Rho', aux, dfftp)

  do ig =1, ngm
     wg_corr(ig) = omega * REAL(aux(dfftp%nl(ig))) - env_smooth_coulomb_g( tpiba2*gg(ig))
  end do
  wg_corr(:) =  wg_corr(:) * exp(-tpiba2*gg(:)*beta/4._dp)**2
  !
  if (gamma_only) wg_corr(gstart:ngm) = 2.d0 * wg_corr(gstart:ngm)
!
  wg_corr_is_updated = .true.

  DEALLOCATE (aux)

  RETURN

  END SUBROUTINE env_init_wg_corr 
!----------------------------------------------------------------------------
  REAL(DP) FUNCTION env_wg_corr_ewald ( omega, ntyp, ngm, zv, strf )
!----------------------------------------------------------------------------
  INTEGER, INTENT(IN) :: ntyp, ngm
  REAL(DP), INTENT(IN) :: omega, zv(ntyp)
  COMPLEX(DP), INTENT(IN) :: strf(ngm,ntyp)
  INTEGER :: ig
  COMPLEX(DP)  :: rhoion

  IF (.NOT.wg_corr_is_updated) CALL env_init_wg_corr
!
  env_wg_corr_ewald = 0._dp
  DO ig=1,ngm
     rhoion = SUM (zv(1:ntyp)* strf(ig,1:ntyp) ) / omega
     env_wg_corr_ewald = env_wg_corr_ewald + ABS(rhoion)**2 * wg_corr(ig) 
  END DO
  env_wg_corr_ewald = 0.5_dp * e2 * env_wg_corr_ewald * omega
!  write(*,*) "ewald correction   = ", wg_corr_ewald

  END FUNCTION env_wg_corr_ewald
!----------------------------------------------------------------------------
  REAL(DP) FUNCTION env_smooth_coulomb_r(r)
!----------------------------------------------------------------------------
  REAL(DP), INTENT(IN) :: r
  REAL(DP), EXTERNAL :: qe_erf
!  smooth_coulomb_r = sqrt(2._dp*alpha/tpi)**3 * exp(-alpha*r*r) ! to be modified
  IF (r>1.e-6_dp) THEN
    env_smooth_coulomb_r = qe_erf(sqrt(alpha)*r)/r
  ELSE
    env_smooth_coulomb_r = 2._dp/sqrt(pi) * sqrt(alpha)
  END IF

  END FUNCTION env_smooth_coulomb_r
!----------------------------------------------------------------------------
  REAL(DP) FUNCTION env_smooth_coulomb_g(q2)
!----------------------------------------------------------------------------
  REAL(DP), INTENT(IN) :: q2
!  smooth_coulomb_g = exp(-q2/4._dp/alpha) ! to be modified
  IF (q2>1.e-6_dp) THEN
    env_smooth_coulomb_g = fpi * exp(-q2/4._dp/alpha)/q2 ! to be modified
  ELSE 
    env_smooth_coulomb_g = - 1._dp * fpi * (1._dp/4._dp/alpha + 2._dp*beta/4._dp)
  END IF
  END FUNCTION env_smooth_coulomb_g
!----------------------------------------------------------------------------

END MODULE env_martyna_tuckerman
