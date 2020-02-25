!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
!
MODULE scf
  !  
  !  This module contains variables and auxiliary routines needed for
  !  the self-consistent cycle
  !
  !  ROUTINES: allocate_scf_type
  !
  USE kinds,      ONLY : DP
  USE lsda_mod,     ONLY : nspin
  USE ions_base,    ONLY : nat
  USE fft_base,     ONLY : dfftp
  USE fft_interfaces,ONLY: invfft
  USE gvect,        ONLY : ngm
  USE gvecs,        ONLY : ngms
  USE paw_variables,ONLY : okpaw
  USE uspp_param,   ONLY : nhm
  USE extfield,     ONLY : dipfield, emaxpos, eopreg, edir
  USE control_flags,ONLY : lxdm
  !
  SAVE
  !
! Details of PAW implementation:
! NOTE: scf_type is used for two different quantities: density and potential.
!       These correspond, for PAW, to becsum and D coefficients.
!       Due to interference with the ultrasoft routines only the becsum part
!       is stored in the structure (at the moment).
!       This only holds for scf_type; mix_type is not affected.
! NOTE: rho%bec is different from becsum for two reasons:
!       1. rho%bec is mixed, while becsum is not
!       2. for npool > 1 rho%bec is collected, becsum is not
!          ( this is necessary to make the stress work)

  TYPE scf_type
     REAL(DP),   ALLOCATABLE :: of_r(:,:)  ! the charge density in R-space
     COMPLEX(DP),ALLOCATABLE :: of_g(:,:)  ! the charge density in G-space
     REAL(DP),   ALLOCATABLE :: kin_r(:,:) ! the kinetic energy density in R-space
     COMPLEX(DP),ALLOCATABLE :: kin_g(:,:) ! the kinetic energy density in G-space
     REAL(DP),   ALLOCATABLE :: ns(:,:,:,:)! the LDA+U occupation matrix
     COMPLEX(DP),ALLOCATABLE :: ns_nc(:,:,:,:)!     ---       noncollinear case
     REAL(DP),   ALLOCATABLE :: bec(:,:,:) ! the PAW hamiltonian elements
  END TYPE scf_type
  !
  TYPE mix_type
     COMPLEX(DP), ALLOCATABLE :: of_g(:,:)  ! the charge density in G-space
     COMPLEX(DP), ALLOCATABLE :: kin_g(:,:) ! the charge density in G-space
     REAL(DP),    ALLOCATABLE :: ns(:,:,:,:)! the LDA+U occupation matrix 
     COMPLEX(DP), ALLOCATABLE :: ns_nc(:,:,:,:)!     ---     noncollinear case 
     REAL(DP),    ALLOCATABLE :: bec(:,:,:) ! PAW corrections to hamiltonian
     REAL(DP)                 :: el_dipole  ! electrons dipole
  END TYPE mix_type

  type (scf_type) :: rho  ! the charge density and its other components
  type (scf_type) :: v    ! the scf potential
  type (scf_type) :: vnew ! used to correct the forces

  REAL(DP) :: v_of_0    ! vltot(G=0)      
  REAL(DP), ALLOCATABLE :: &
       vltot(:),       &! the local potential in real space
       vrs(:,:),       &! the total pot. in real space (smooth grid)
       rho_core(:),    &! the core charge in real space
       kedtau(:,:)      ! position dependent kinetic energy enhancement factor
  COMPLEX(DP), ALLOCATABLE :: &
       rhog_core(:)     ! the core charge in reciprocal space

  INTEGER, PRIVATE  :: record_length, &
                       rlen_rho=0,  rlen_kin=0,  rlen_ldaU=0,  rlen_bec=0,&
                       rlen_dip=0, &
                       start_rho=0, start_kin=0, start_ldaU=0, start_bec=0, &
                       start_dipole=0
  COMPLEX(DP), PRIVATE, ALLOCATABLE:: io_buffer(:)
CONTAINS

 SUBROUTINE create_scf_type ( rho, do_not_allocate_becsum )
   IMPLICIT NONE
   TYPE (scf_type) :: rho
   LOGICAL,INTENT(IN),OPTIONAL :: do_not_allocate_becsum ! PAW hack
   LOGICAL                     :: allocate_becsum        ! PAW hack
   allocate ( rho%of_r( dfftp%nnr, nspin) )
   allocate ( rho%of_g( ngm, nspin ) )
   allocate ( rho%kin_r(1,1) )
   allocate ( rho%kin_g(1,1) )

   if (okpaw) then ! See the top of the file for clarification
      if(present(do_not_allocate_becsum)) then
         allocate_becsum = .not. do_not_allocate_becsum
      else
         allocate_becsum = .true.
      endif
      if(allocate_becsum) allocate (rho%bec(nhm*(nhm+1)/2,nat,nspin))
   endif
   
 return
 END SUBROUTINE create_scf_type

 SUBROUTINE destroy_scf_type ( rho )
   IMPLICIT NONE
   TYPE (scf_type) :: rho

   if (ALLOCATED(rho%of_r))  deallocate(rho%of_r)
   if (ALLOCATED(rho%of_g))  deallocate(rho%of_g)
   if (ALLOCATED(rho%kin_r)) deallocate(rho%kin_r)
   if (ALLOCATED(rho%kin_g)) deallocate(rho%kin_g)
   if (ALLOCATED(rho%ns))    deallocate(rho%ns)
   if (ALLOCATED(rho%ns_nc))    deallocate(rho%ns_nc)
   if (ALLOCATED(rho%bec))   deallocate(rho%bec)

   return
 END SUBROUTINE destroy_scf_type
 !
 
 SUBROUTINE create_mix_type ( rho )
   IMPLICIT NONE
   TYPE (mix_type) :: rho
   allocate ( rho%of_g( ngms, nspin ) )
   rho%of_g = 0._dp
   if (okpaw) then
      allocate (rho%bec(nhm*(nhm+1)/2,nat,nspin))
      rho%bec   = 0._dp
   end if
   rho%el_dipole =  0._dp
   
 return
 END SUBROUTINE create_mix_type

 SUBROUTINE destroy_mix_type ( rho )
   IMPLICIT NONE
   TYPE (mix_type) :: rho

   if (ALLOCATED(rho%of_g))  deallocate(rho%of_g)
   if (ALLOCATED(rho%kin_g)) deallocate(rho%kin_g)
   if (ALLOCATED(rho%ns))    deallocate(rho%ns)
   if (ALLOCATED(rho%ns_nc))    deallocate(rho%ns_nc)
   if (ALLOCATED(rho%bec))   deallocate(rho%bec)

   return
 END SUBROUTINE destroy_mix_type
 !
 subroutine assign_scf_to_mix_type(rho_s, rho_m)
   IMPLICIT NONE
   TYPE (scf_type), INTENT(IN)  :: rho_s
   TYPE (mix_type), INTENT(INOUT) :: rho_m
   REAL(DP) :: e_dipole
      
   rho_m%of_g(1:ngms,:) = rho_s%of_g(1:ngms,:)
   if (okpaw)         rho_m%bec = rho_s%bec
   
   if (dipfield) then
      CALL compute_el_dip(emaxpos, eopreg, edir, rho_s%of_r(:,1), e_dipole)
      rho_m%el_dipole = e_dipole
   endif
   
 return
 end subroutine assign_scf_to_mix_type
 !
 !----------------------------------------------------------------------------
 subroutine scf_type_COPY (X,Y)
  !----------------------------------------------------------------------------
  ! works like DCOPY for scf_type copy variables :  Y = X 
  USE kinds, ONLY : DP
  IMPLICIT NONE
  TYPE(scf_type), INTENT(IN)    :: X
  TYPE(scf_type), INTENT(INOUT) :: Y
  Y%of_r  = X%of_r
  Y%of_g  = X%of_g
  if (lxdm) then
     Y%kin_r = X%kin_r
     Y%kin_g = X%kin_g
  end if
  if (okpaw)      Y%bec = X%bec
  !
  RETURN
 end subroutine scf_type_COPY
 !
 !----------------------------------------------------------------------------
 subroutine mix_type_AXPY (A,X,Y)
  !----------------------------------------------------------------------------
  ! works like daxpy for scf_type variables :  Y = A * X + Y
  ! NB: A is a REAL(DP) number
  USE kinds, ONLY : DP
  IMPLICIT NONE
  REAL(DP)                      :: A
  TYPE(mix_type), INTENT(IN)    :: X
  TYPE(mix_type), INTENT(INOUT) :: Y
  Y%of_g  = Y%of_g  + A * X%of_g
  if (okpaw)     Y%bec = Y%bec + A * X%bec
  if (dipfield)  Y%el_dipole =  Y%el_dipole + A * X%el_dipole
  !
  RETURN
 END SUBROUTINE mix_type_AXPY
 !
 !----------------------------------------------------------------------------
 subroutine mix_type_COPY (X,Y)
  !----------------------------------------------------------------------------
  ! works like DCOPY for mix_type copy variables :  Y = X 
  USE kinds, ONLY : DP
  IMPLICIT NONE
  TYPE(mix_type), INTENT(IN)    :: X
  TYPE(mix_type), INTENT(INOUT) :: Y
  Y%of_g  = X%of_g
  if (okpaw)      Y%bec = X%bec
  if (dipfield)   Y%el_dipole =  X%el_dipole
  !
  RETURN
 end subroutine mix_type_COPY
 !
 !----------------------------------------------------------------------------
 subroutine mix_type_SCAL (A,X)
  !----------------------------------------------------------------------------
  ! works like DSCAL for mix_type copy variables :  X = A * X 
  ! NB: A is a REAL(DP) number
  USE kinds, ONLY : DP
  IMPLICIT NONE
  REAL(DP),       INTENT(IN)    :: A
  TYPE(mix_type), INTENT(INOUT) :: X
  X%of_g(:,:)  = A * X%of_g(:,:)
  if (okpaw)      X%bec= A * X%bec
  if (dipfield)   X%el_dipole =  A * X%el_dipole
  !
  RETURN
 end subroutine mix_type_SCAL
 !
 !----------------------------------------------------------------------------
 FUNCTION rho_ddot( rho1, rho2, gf )
  !----------------------------------------------------------------------------
  !
  ! ... calculates 4pi/G^2*rho1(-G)*rho2(G) = V1_Hartree(-G)*rho2(G)
  ! ... used as an estimate of the self-consistency error on the energy
  !
  USE kinds,         ONLY : DP
  USE constants,     ONLY : e2, tpi, fpi
  USE cell_base,     ONLY : omega, tpiba2
  USE gvect,         ONLY : gg, gstart
  USE control_flags, ONLY : gamma_only
  !USE paw_onecenter, ONLY : paw_ddot
  USE mp_bands,      ONLY : intra_bgrp_comm
  USE mp,            ONLY : mp_sum
  !
  IMPLICIT NONE
  !
  type(mix_type), INTENT(IN) :: rho1, rho2
  INTEGER,        INTENT(IN) :: gf
  REAL(DP)                :: rho_ddot
  !
  REAL(DP) :: fac
  INTEGER  :: ig
  !
  fac = e2 * fpi / tpiba2
  !
  rho_ddot = 0.D0
  !
  DO ig = gstart, gf
     !
     rho_ddot = rho_ddot + &
                REAL( CONJG( rho1%of_g(ig,1) )*rho2%of_g(ig,1), DP ) / gg(ig)
     !
  END DO
  !
  rho_ddot = fac*rho_ddot
  !
  IF ( gamma_only ) rho_ddot = 2.D0 * rho_ddot
  !
  IF ( nspin >= 2 )  THEN
     !
     fac = e2*fpi / tpi**2  ! lambda=1 a.u.
     !
     IF ( gstart == 2 ) THEN
        !
        rho_ddot = rho_ddot + &
                fac * SUM( REAL( CONJG( rho1%of_g(1,2:nspin))*(rho2%of_g(1,2:nspin) ), DP ) )
        !
     END IF
     !
     IF ( gamma_only ) fac = 2.D0 * fac
     !
     DO ig = gstart, gf
        !
        rho_ddot = rho_ddot + &
                fac * SUM( REAL( CONJG( rho1%of_g(ig,2:nspin))*(rho2%of_g(ig,2:nspin) ), DP ) )
     !
     END DO
     !
  END IF
  !
  rho_ddot = rho_ddot * omega * 0.5D0
  !
  CALL mp_sum(  rho_ddot , intra_bgrp_comm )
  ! 
  ! Beware: paw_ddot has a hidden parallelization on all processors
  !         it must be called on all processors or else it will hang
  ! Beware: commented out because it yields too often negative values
  ! IF (okpaw)         rho_ddot = rho_ddot + paw_ddot(rho1%bec, rho2%bec)
  IF (dipfield)      rho_ddot = rho_ddot + (e2/2.0_DP)* &
                                    (rho1%el_dipole * rho2%el_dipole)*omega/fpi

  RETURN
  !
 END FUNCTION rho_ddot
 !
!----------------------------------------------------------------------------
FUNCTION tauk_ddot( rho1, rho2, gf )
  !----------------------------------------------------------------------------
  !
  ! ... calculates 4pi/G^2*rho1(-G)*rho2(G) = V1_Hartree(-G)*rho2(G)
  ! ... used as an estimate of the self-consistency error on the energy
  !
  USE kinds,         ONLY : DP
  USE constants,     ONLY : e2, tpi, fpi
  USE cell_base,     ONLY : omega, tpiba2
  USE gvect,         ONLY : gg, gstart
  USE control_flags, ONLY : gamma_only
  USE mp_bands,      ONLY : intra_bgrp_comm
  USE mp,            ONLY : mp_sum
  !
  IMPLICIT NONE
  !
  type(mix_type), INTENT(IN) :: rho1, rho2
  INTEGER,     INTENT(IN) :: gf
  REAL(DP)                :: tauk_ddot
  !
  REAL(DP) :: fac
  INTEGER  :: ig
  !
  tauk_ddot = 0.D0
  !
!  write (*,*) rho1%kin_g(1:4,1)
!  if (.true. ) stop
  !
  DO ig = gstart, gf
     tauk_ddot = tauk_ddot + &
                 REAL( CONJG( rho1%kin_g(ig,1) )*rho2%kin_g(ig,1) ) 
  END DO
  !
  IF ( nspin==1 .and. gamma_only ) tauk_ddot = 2.D0 * tauk_ddot
  !
  ! ... G=0 term
  !
  IF ( gstart == 2 ) THEN
     !
     tauk_ddot = tauk_ddot + &
                 REAL( CONJG( rho1%kin_g(1,1) ) * rho2%kin_g(1,1) )
     !
  END IF
  !
  IF ( nspin >= 2 ) THEN
     !
     DO ig = gstart, gf
        !
        tauk_ddot = tauk_ddot + &
                SUM( REAL( CONJG( rho1%kin_g(1,2:nspin))*(rho2%kin_g(1,2:nspin) ), DP ) )
        !
     END DO
     !
     IF ( gamma_only ) tauk_ddot = 2.D0 * tauk_ddot
     !
     ! ... G=0 term
     !
     IF ( gstart == 2 ) THEN
        !
        tauk_ddot = tauk_ddot + &
                SUM( REAL( CONJG( rho1%kin_g(1,1:nspin))*(rho2%kin_g(1,1:nspin) ), DP ) )
        !
     END IF
     IF ( nspin == 2 ) tauk_ddot = 0.5D0 *  tauk_ddot 
     !
  END IF
  !
  fac = e2 * fpi / tpi**2  ! lambda = 1 a.u.
  !
  tauk_ddot = fac * tauk_ddot * omega * 0.5D0
  !
  CALL mp_sum(  tauk_ddot , intra_bgrp_comm )
  !
  RETURN
  !
END FUNCTION tauk_ddot
 !----------------------------------------------------------------------------
 FUNCTION local_tf_ddot( rho1, rho2, ngm0 )
  !----------------------------------------------------------------------------
  !
  ! ... calculates 4pi/G^2*rho1(-G)*rho2(G) = V1_Hartree(-G)*rho2(G)
  ! ... used as an estimate of the self-consistency error on the energy
  !
  USE kinds,         ONLY : DP
  USE constants,     ONLY : e2, fpi
  USE cell_base,     ONLY : omega, tpiba2
  USE gvect,         ONLY : gg, gstart
  USE control_flags, ONLY : gamma_only
  USE mp_bands,      ONLY : intra_bgrp_comm
  USE mp,            ONLY : mp_sum
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN)     :: ngm0
  COMPLEX(DP), INTENT(IN) :: rho1(ngm0), rho2(ngm0)
  REAL(DP)                :: local_tf_ddot
  !
  REAL(DP) :: fac
  INTEGER  :: ig
  !
  local_tf_ddot = 0.D0
  !
  fac = e2 * fpi / tpiba2
  !
  !$omp parallel do reduction(+:local_tf_ddot)
  DO ig = gstart, ngm0
     local_tf_ddot = local_tf_ddot + REAL( CONJG(rho1(ig))*rho2(ig) ) / gg(ig)
  END DO
  !$omp end parallel do
  !
  local_tf_ddot = fac * local_tf_ddot * omega * 0.5D0
  !
  IF ( gamma_only ) local_tf_ddot = 2.D0 * local_tf_ddot
  !
  CALL mp_sum(  local_tf_ddot , intra_bgrp_comm )
  !
  RETURN
  !
 END FUNCTION local_tf_ddot
 !
 SUBROUTINE bcast_scf_type ( rho, root, comm )
  !----------------------------------------------------------------------------
  ! ... Broadcast all mixed quantities from first pool to all others
  ! ... Needed to prevent divergencies in k-point parallization
  !
  USE mp,            ONLY : mp_bcast
  !
  IMPLICIT NONE
  !
  type(scf_type), INTENT(INOUT) :: rho
  INTEGER, INTENT(IN) :: root, comm
  !
  CALL mp_bcast ( rho%of_g, root, comm )
  CALL mp_bcast ( rho%of_r, root, comm )
  IF (lxdm) THEN
     CALL mp_bcast ( rho%kin_g, root, comm )
     CALL mp_bcast ( rho%kin_r, root, comm )
  END IF
  IF ( okpaw )        CALL mp_bcast ( rho%bec,   root, comm )
  !
  END SUBROUTINE
  !
  SUBROUTINE rhoz_or_updw( rho, sp, dir )
  !--------------------------------------------------------------------------
  ! ... Converts rho(up,dw) into rho(up+dw,up-dw) if dir='->rhoz' and
  ! ... vice versa if dir='->updw'.
  !
  USE gvect,  ONLY : ngm
  !
  IMPLICIT NONE
  !
  type(scf_type), INTENT(INOUT) :: rho
  CHARACTER(len=*), INTENT(IN) :: dir, sp
  INTEGER :: ir
  REAL(DP) :: vi
  !
  IF ( nspin /= 2 ) RETURN
  !
  vi = 0._dp
  IF (dir == '->updw')  vi = 0.5_dp
  IF (dir == '->rhoz')  vi = 1.0_dp
  IF (vi  == 0._dp)  CALL errore( 'rhoz_or_updw', 'wrong input', 1 )
  !
  IF ( sp /= 'only_g' ) THEN
!   !$omp parallel do
     DO ir = 1, dfftp%nnr  
        rho%of_r(ir,1) = ( rho%of_r(ir,1) + rho%of_r(ir,nspin) ) * vi
        rho%of_r(ir,nspin) = rho%of_r(ir,1) - rho%of_r(ir,nspin) * vi * 2._dp
     ENDDO
!   !$omp end parallel do
  ENDIF
  IF ( sp /= 'only_r' ) THEN
!   !$omp parallel do
     DO ir = 1, ngm
        rho%of_g(ir,1) = ( rho%of_g(ir,1) + rho%of_g(ir,nspin) ) * vi
        rho%of_g(ir,nspin) = rho%of_g(ir,1) - rho%of_g(ir,nspin) * vi * 2._dp
     ENDDO
!   !$omp end parallel do  
  ENDIF
  !
  RETURN
  !
  END SUBROUTINE
  !
  !
END MODULE scf
