!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE v_of_rho( rho, rho_core, rhog_core, &
                     ehart, etxc, vtxc, eth, etotefield, charge, v )
  !----------------------------------------------------------------------------
  !
  ! ... This routine computes the Hartree and Exchange and Correlation
  ! ... potential and energies which corresponds to a given charge density
  ! ... The XC potential is computed in real space, while the
  ! ... Hartree potential is computed in reciprocal space.
  !
  USE kinds,            ONLY : DP
  USE fft_base,         ONLY : dfftp
  USE gvect,            ONLY : ngm
  USE ions_base,        ONLY : nat, tau
  USE scf,              ONLY : scf_type
  USE cell_base,        ONLY : alat
  USE control_flags,    ONLY : ts_vdw
  !
  IMPLICIT NONE
  !
  TYPE(scf_type), INTENT(INOUT) :: rho  ! the valence charge
  TYPE(scf_type), INTENT(INOUT) :: v ! the scf (Hxc) potential 
  !!!!!!!!!!!!!!!!! NB: NOTE that in F90 derived data type must be INOUT and 
  !!!!!!!!!!!!!!!!! not just OUT because otherwise their allocatable or pointer
  !!!!!!!!!!!!!!!!! components are NOT defined !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  REAL(DP), INTENT(IN) :: rho_core(dfftp%nnr)
    ! the core charge
  COMPLEX(DP), INTENT(IN) :: rhog_core(ngm)
    ! the core charge in reciprocal space
  REAL(DP), INTENT(OUT) :: vtxc, etxc, ehart, eth, charge
    ! the integral V_xc * rho
    ! the E_xc energy
    ! the hartree energy
    ! the hubbard energy
    ! the integral of the charge
  REAL(DP), INTENT(INOUT) :: etotefield
    ! electric field energy - inout due to the screwed logic of add_efield
  ! ! 
  INTEGER :: is, ir
  print *, "*************************************************"
  print *, "running from Environ!!"
  print *, "*************************************************"
  !
  CALL start_clock( 'v_of_rho' )
  !
  ! ... calculate exchange-correlation potential
  !
  CALL v_xc( rho, rho_core, rhog_core, etxc, vtxc, v%of_r )
  !
  ! ... add a magnetic field  (if any)
  !
  CALL add_bfield( v%of_r, rho%of_r )
  !
  ! ... calculate hartree potential
  !
  CALL v_h( rho%of_g(:,1), ehart, charge, v%of_r )
  !
  CALL stop_clock( 'v_of_rho' )
  !
  RETURN
  !
END SUBROUTINE v_of_rho
!
!----------------------------------------------------------------------------
SUBROUTINE v_h( rhog, ehart, charge, v )
  !----------------------------------------------------------------------------
  !
  ! ... Hartree potential VH(r) from n(G)
  !
  USE constants, ONLY : fpi, e2
  USE kinds,     ONLY : DP
  USE fft_base,  ONLY : dfftp
  USE fft_interfaces,ONLY : invfft
  USE gvect,     ONLY : ngm, gg, gstart
  USE lsda_mod,  ONLY : nspin
  USE cell_base, ONLY : omega, tpiba2
  USE control_flags, ONLY : gamma_only
  USE mp_bands,  ONLY: intra_bgrp_comm
  USE mp,        ONLY: mp_sum
  USE martyna_tuckerman, ONLY : wg_corr_h, do_comp_mt
  USE esm,       ONLY: do_comp_esm, esm_hartree, esm_bc
  !
  IMPLICIT NONE
  !
  COMPLEX(DP), INTENT(IN)  :: rhog(ngm)
  REAL(DP),  INTENT(INOUT) :: v(dfftp%nnr,nspin)
  REAL(DP),    INTENT(OUT) :: ehart, charge
  !
  REAL(DP)              :: fac
  REAL(DP), ALLOCATABLE :: aux1(:,:)
  REAL(DP)              :: rgtot_re, rgtot_im, eh_corr
  INTEGER               :: is, ig
  COMPLEX(DP), ALLOCATABLE :: aux(:), rgtot(:), vaux(:)
  INTEGER               :: nt
  !
  CALL start_clock( 'v_h' )
  !
  ALLOCATE( aux( dfftp%nnr ), aux1( 2, ngm ) )
  charge = 0.D0
  !
  IF ( gstart == 2 ) THEN
     !
     charge = omega*REAL( rhog(1) )
     !
  END IF
  !
  CALL mp_sum(  charge , intra_bgrp_comm )
  !
  ! ... calculate hartree potential in G-space (NB: V(G=0)=0 )
  !
  IF ( do_comp_esm .and. ( esm_bc .ne. 'pbc' ) ) THEN
     !
     ! ... calculate modified Hartree potential for ESM
     !
     CALL esm_hartree (rhog, ehart, aux)
     !
  ELSE
     !
     ehart     = 0.D0
     aux1(:,:) = 0.D0
!$omp parallel do private( fac, rgtot_re, rgtot_im ), reduction(+:ehart)
        DO ig = gstart, ngm
           !
           fac = 1.D0 / gg(ig) 
           !
           rgtot_re = REAL(  rhog(ig) )
           rgtot_im = AIMAG( rhog(ig) )
           !
           ehart = ehart + ( rgtot_re**2 + rgtot_im**2 ) * fac
           !
           aux1(1,ig) = rgtot_re * fac
           aux1(2,ig) = rgtot_im * fac
           !
        ENDDO
!$omp end parallel do
   ENDIF
   !
   fac = e2 * fpi / tpiba2
   !
   ehart = ehart * fac
   !
   aux1 = aux1 * fac
   !
   IF ( gamma_only ) THEN
      !
      ehart = ehart * omega
      !
   ELSE 
      !
      ehart = ehart * 0.5D0 * omega
      !
   END IF
   ! 
   if (do_comp_mt) then
      ALLOCATE( vaux( ngm ), rgtot(ngm) )
      rgtot(:) = rhog(:)
      CALL wg_corr_h (omega, ngm, rgtot, vaux, eh_corr)
      aux1(1,1:ngm) = aux1(1,1:ngm) + REAL( vaux(1:ngm))
      aux1(2,1:ngm) = aux1(2,1:ngm) + AIMAG(vaux(1:ngm))
      ehart = ehart + eh_corr
      DEALLOCATE( rgtot, vaux )
   end if
   !
   CALL mp_sum(  ehart , intra_bgrp_comm )
   ! 
   aux(:) = 0.D0
   !
   aux(dfftp%nl(1:ngm)) = CMPLX ( aux1(1,1:ngm), aux1(2,1:ngm), KIND=dp )
   !
   IF ( gamma_only ) THEN
      !
      aux(dfftp%nlm(1:ngm)) = CMPLX ( aux1(1,1:ngm), -aux1(2,1:ngm), KIND=dp )
      !
   END IF
  !
  ! ... transform hartree potential to real space
  !
  CALL invfft ('Rho', aux, dfftp)
  !
  ! ... add hartree potential to the xc potential
  !
  IF ( nspin == 4 ) THEN
     !
     v(:,1) = v(:,1) + DBLE (aux(:))
     !
  ELSE
     !
     DO is = 1, nspin
        !
        v(:,is) = v(:,is) + DBLE (aux(:))
        !
     END DO
     !
  END IF
  !
  DEALLOCATE( aux, aux1 )
  !
  CALL stop_clock( 'v_h' )
  !
  RETURN
  !
END SUBROUTINE v_h

!----------------------------------------------------------------------------
SUBROUTINE v_h_of_rho_r( rhor, ehart, charge, v )
  !----------------------------------------------------------------------------
  !
  ! ... Hartree potential VH(r) from a density in R space n(r)
  !
  USE kinds,           ONLY : DP
  USE fft_base,        ONLY : dfftp
  USE fft_interfaces,  ONLY : fwfft
  USE lsda_mod,        ONLY : nspin
  !
  IMPLICIT NONE
  !
  ! ... Declares variables
  !
  REAL( DP ), INTENT(IN)     :: rhor( dfftp%nnr )
  REAL( DP ), INTENT(INOUT)  :: v( dfftp%nnr )
  REAL( DP ), INTENT(OUT)    :: ehart, charge
  !
  ! ... Local variables
  !
  COMPLEX( DP ), ALLOCATABLE :: rhog( : )
  COMPLEX( DP ), ALLOCATABLE :: aux( : )
  REAL( DP ), ALLOCATABLE :: vaux(:,:)
  INTEGER :: is
  !
  ! ... bring the (unsymmetrized) rho(r) to G-space (use aux as work array)
  !
  ALLOCATE( rhog( dfftp%ngm ) )
  ALLOCATE( aux( dfftp%nnr ) )
  aux = CMPLX(rhor,0.D0,kind=dp)
  CALL fwfft ('Rho', aux, dfftp)
  rhog(:) = aux(dfftp%nl(:))
  DEALLOCATE( aux )
  !
  ! ... compute VH(r) from n(G)
  !
  ALLOCATE( vaux( dfftp%nnr, nspin ) )
  vaux = 0.D0
  CALL v_h( rhog, ehart, charge, vaux )
  v(:) = v(:) + vaux(:,1)
  !
  DEALLOCATE( rhog )
  DEALLOCATE( vaux )
  !
  RETURN
  !
END SUBROUTINE v_h_of_rho_r
!----------------------------------------------------------------------------
SUBROUTINE gradv_h_of_rho_r( rho, gradv )
  !----------------------------------------------------------------------------
  !
  ! ... Gradient of Hartree potential in R space from a total
  !     (spinless) density in R space n(r)
  !
  USE kinds,           ONLY : DP
  USE fft_base,        ONLY : dfftp
  USE fft_interfaces,  ONLY : fwfft, invfft
  USE constants,       ONLY : fpi, e2
  USE control_flags,   ONLY : gamma_only
  USE cell_base,       ONLY : tpiba, omega
  USE gvect,           ONLY : ngm, gg, gstart, g
  USE martyna_tuckerman, ONLY : wg_corr_h, do_comp_mt
  !
  IMPLICIT NONE
  !
  ! ... Declares variables
  !
  REAL( DP ), INTENT(IN)     :: rho( dfftp%nnr )
  REAL( DP ), INTENT(OUT)    :: gradv( 3, dfftp%nnr )
  !
  ! ... Local variables
  !
  COMPLEX( DP ), ALLOCATABLE :: rhoaux( : )
  COMPLEX( DP ), ALLOCATABLE :: gaux( : )
  COMPLEX( DP ), ALLOCATABLE :: rgtot(:), vaux(:)
  REAL( DP )                 :: fac, eh_corr
  INTEGER                    :: ig, ipol
  !
  ! ... Bring rho to G space
  !
  ALLOCATE( rhoaux( dfftp%nnr ) )
  rhoaux( : ) = CMPLX( rho( : ), 0.D0, KIND=dp ) 
  !
  CALL fwfft('Rho', rhoaux, dfftp)
  !
  ! ... Compute total potential in G space
  !
  ALLOCATE( gaux( dfftp%nnr ) )
  !
  DO ipol = 1, 3
    !
    gaux(:) = (0.0_dp,0.0_dp)
    !
    DO ig = gstart, ngm
      !
      fac = g(ipol,ig) / gg(ig)
      gaux(dfftp%nl(ig)) = CMPLX(-AIMAG(rhoaux(dfftp%nl(ig))),REAL(rhoaux(dfftp%nl(ig))),kind=dp) * fac 
      !
    END DO
    !
    ! ...and add the factor e2*fpi/2\pi/a coming from the missing prefactor of 
    !  V = e2 * fpi divided by the 2\pi/a factor missing in G  
    !
    fac = e2 * fpi / tpiba
    gaux = gaux * fac 
    !
    ! ...add martyna-tuckerman correction, if needed
    ! 
    if (do_comp_mt) then
       ALLOCATE( vaux( ngm ), rgtot(ngm) )
       rgtot(1:ngm) = rhoaux(dfftp%nl(1:ngm))
       CALL wg_corr_h (omega, ngm, rgtot, vaux, eh_corr)
       DO ig = gstart, ngm
         fac = g(ipol,ig) * tpiba
         gaux(dfftp%nl(ig)) = gaux(dfftp%nl(ig)) + CMPLX(-AIMAG(vaux(ig)),REAL(vaux(ig)),kind=dp)*fac 
       END DO
       DEALLOCATE( rgtot, vaux )
    end if
    !
    IF ( gamma_only ) THEN
      !
      gaux(dfftp%nlm(:)) = &
        CMPLX( REAL( gaux(dfftp%nl(:)) ), -AIMAG( gaux(dfftp%nl(:)) ) ,kind=DP)
       !
    END IF
    !
    ! ... bring back to R-space, (\grad_ipol a)(r) ...
    !
    CALL invfft ('Rho', gaux, dfftp)
    !
    gradv(ipol,:) = REAL( gaux(:) )
    !
  ENDDO
  !
  DEALLOCATE(gaux)
  !
  DEALLOCATE(rhoaux)
  !
  RETURN
  !
END SUBROUTINE gradv_h_of_rho_r
