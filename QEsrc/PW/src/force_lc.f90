!
! Copyright (C) 2001-2011 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine force_lc (nat, tau, ityp, alat, omega, ngm, ngl, &
     igtongl, g, rho, nl, nspin, gstart, gamma_only, vloc, forcelc)
  !----------------------------------------------------------------------
  !
  USE kinds
  USE constants, ONLY : tpi
  USE mp_bands,  ONLY : intra_bgrp_comm
  USE mp,        ONLY : mp_sum
  USE fft_base,  ONLY : dfftp
  USE fft_interfaces, ONLY : fwfft
  USE esm,       ONLY : esm_force_lc, do_comp_esm, esm_bc
  implicit none
  !
  !   first the dummy variables
  !
  integer, intent(in) :: nat, ngm, nspin, ngl, gstart, &
                         igtongl (ngm), nl (ngm), ityp (nat)
  ! nat:    number of atoms in the cell
  ! ngm:    number of G vectors
  ! nspin:  number of spin polarizations
  ! ngl:    number of shells
  ! igtongl correspondence G <-> shell of G
  ! nl:     correspondence fft mesh <-> G vec
  ! ityp:   types of atoms

  logical, intent(in) :: gamma_only

  real(DP), intent(in) :: tau (3, nat), g (3, ngm), vloc (ngl, * ), &
       rho (dfftp%nnr, nspin), alat, omega
  ! tau:  coordinates of the atoms
  ! g:    coordinates of G vectors
  ! vloc: local potential
  ! rho:  valence charge
  ! alat: lattice parameter
  ! omega: unit cell volume

  real(DP), intent(out) :: forcelc (3, nat)
  ! the local-potential contribution to forces on atoms

  integer :: ipol, ig, na
  ! counter on polarizations
  ! counter on G vectors
  ! counter on atoms

  complex(DP), allocatable :: aux (:)
  ! auxiliary space for FFT
  real(DP) :: arg, fact
  !
  ! contribution to the force from the local part of the bare potential
  ! F_loc = Omega \Sum_G n*(G) d V_loc(G)/d R_i
  !
  allocate (aux(dfftp%nnr))
  if ( nspin == 2) then
      aux(:) = CMPLX( rho(:,1)+rho(:,2), 0.0_dp, kind=dp )
  else
      aux(:) = CMPLX( rho(:,1), 0.0_dp, kind=dp )
  end if
  CALL fwfft ('Dense', aux, dfftp)
  !
  !    aux contains now  n(G)
  !
  if (gamma_only) then
     fact = 2.d0
  else
     fact = 1.d0
  end if
  do na = 1, nat
     do ipol = 1, 3
        forcelc (ipol, na) = 0.d0
     enddo
     ! contribution from G=0 is zero
     do ig = gstart, ngm
        arg = (g (1, ig) * tau (1, na) + g (2, ig) * tau (2, na) + &
               g (3, ig) * tau (3, na) ) * tpi
        do ipol = 1, 3
           forcelc (ipol, na) = forcelc (ipol, na) + &
                g (ipol, ig) * vloc (igtongl (ig), ityp (na) ) * &
                (sin(arg)*DBLE(aux(nl(ig))) + cos(arg)*AIMAG(aux(nl(ig))) )
        enddo
     enddo
     do ipol = 1, 3
        forcelc (ipol, na) = fact * forcelc (ipol, na) * omega * tpi / alat
     enddo
  enddo
  IF ( do_comp_esm .and. ( esm_bc .ne. 'pbc' ) ) THEN
     !
     ! ... Perform corrections for ESM method (add long-range part)
     !
     CALL esm_force_lc ( aux, forcelc )
  ENDIF
  !
  call mp_sum(  forcelc, intra_bgrp_comm )
  !
  deallocate (aux)
  return
end subroutine force_lc

subroutine external_force_lc( rhor, force )

  use kinds,            only : DP
  use cell_base,        only : at, bg, alat, omega
  use ions_base,        only : nat, ntyp => nsp, ityp, tau, zv, amass
  use fft_base,         only : dfftp
  use fft_interfaces,   only : fwfft
  use gvect,            only : ngm, gstart, ngl, nl, igtongl, g, gg, gcutm
  use lsda_mod,         only : nspin
  use vlocal,           only : strf, vloc
  use control_flags,    only : gamma_only
  use martyna_tuckerman, only: do_comp_mt, wg_corr_force
  implicit none

  real( dp ), intent(in) ::  rhor (dfftp%nnr, nspin)
  real( dp ), intent(out) :: force (3, nat)

  real( dp ), allocatable :: force_tmp(:,:)
  complex( dp ), allocatable :: auxg(:), auxr(:)

  force = 0.0_dp

  allocate(force_tmp(3,nat))

  if ( do_comp_mt) then
     force_tmp = 0.0_dp
     allocate(auxr(dfftp%nnr))
     allocate(auxg(ngm))
     auxg = cmplx(0.0_dp,0.0_dp)
     auxr = cmplx(rhor(:,1),0.0_dp)
     if ( nspin .eq. 2 ) auxr = auxr + cmplx(rhor(:,2),0.0_dp)
     call fwfft ("Dense", auxr, dfftp)
     auxg(:)=auxr(nl(:))
     call wg_corr_force(.false.,omega, nat, ntyp, ityp, ngm, g, tau, zv, strf, &
                        1, auxg, force_tmp)
     deallocate(auxr,auxg)
     force = force + force_tmp
  endif

  force_tmp = 0.0_dp
  call force_lc( nat, tau, ityp, alat, omega, ngm, ngl, igtongl, &
       g, rhor, nl, nspin, gstart, gamma_only, vloc, force_tmp )
  force = force + force_tmp

  return
end subroutine external_force_lc
