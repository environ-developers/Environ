!
! Copyright (C) 2007-2011 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Original version by Minoru Otani (AIST), Yoshio Miura (Tohoku U.),
! Nicephore Bonet (MIT), Nicola Marzari (MIT), Brandon Wood (LLNL), 
! Tadashi Ogitsu (LLNL).
! Constant bias potential (constant-mu) method by Minoru Otani (AIST) and
! Nicephore Bonnet (AIST).
!
! Contains SUBROUTINEs for implementation of
! 1) ESM (Effective Screening Medium Method) developed by M. Otani and 
!    O. Sugino (see PRB 73, 115407 [2006])
! 2) Constant-mu method developed by N. Bonnet, T. Morishita, O. Sugino, 
!    and M. Otani (see PRL 109, 266101 [2012]).
!
! ESM enables description of a surface slab sandwiched between two 
! semi-infinite media, making it possible to deal with polarized surfaces 
! without using dipole corrections. It is useful for simulating interfaces 
! with vacuum, one or more electrodes, or an electrolyte.
!
! Constant-mu scheme with the boundary condition 'bc2' and 'bc3' enables
! description of the system is connected to a potentiostat which preserves
! the Fermi energy of the system as the target Fermi energy (mu).
!
! Modified SUBROUTINEs for calculating the Hartree potential, the local 
! potential, and the Ewald sum are contained here, along with SUBROUTINEs for
! calculating force contributions based on the modified local potential and 
! Ewald term. Constant-mu parts are contained in the fcp.f90.
!
!----------------------------------------------------------------------------
MODULE env_esm
  !--------------------------------------------------------------------------
  !
  ! ... this module contains the variables and SUBROUTINEs needed for the 
  ! ... EFFECTIVE SCREENING MEDIUM (ESM) METHOD 
  !
  USE env_kinds, ONLY : DP
  !
  IMPLICIT NONE
  !
  SAVE
  !
  PRIVATE
  !
  PUBLIC :: do_comp_esm, esm_nfit,  esm_w, esm_a, esm_bc, &
            mill_2d, imill_2d, ngm_2d, env_esm_hartree
  !
  LOGICAL              :: do_comp_esm=.FALSE.
  INTEGER              :: esm_nfit
  REAL(KIND=DP)        :: esm_efield, esm_w, esm_a
  CHARACTER (LEN=3)    :: esm_bc

  INTEGER, ALLOCATABLE :: mill_2d(:,:), imill_2d(:,:)
  INTEGER              :: ngm_2d = 0
  real(DP), external   :: qe_erf, qe_erfc
  !
  CONTAINS
  !-----------------------------------------------------------------------
  !--------------ESM ENERGY AND POTENTIAL SUBROUTINE----------------------
  !-----------------------------------------------------------------------
     subroutine env_esm_hartree( rhog, ehart, aux )
        USE env_kinds,    ONLY : DP
        USE env_gvect,    ONLY : ngm
        USE env_fft_base, ONLY : dfftp
        IMPLICIT NONE
        real(DP)    :: ehart             !  Hartree energy
        complex(DP) :: rhog(ngm)         !  n(G)
        complex(DP) :: aux(dfftp%nnr)    !  v_h(G)

        if( esm_bc == 'pbc' ) then
           call env_esm_hartree_pbc ( rhog, ehart, aux )
        else if ( esm_bc == 'bc1' ) then
           call env_esm_hartree_bc1 ( rhog, ehart, aux )
        else if ( esm_bc == 'bc2' ) then
           call env_esm_hartree_bc2 ( rhog, ehart, aux )
        else if ( esm_bc == 'bc3' ) then
           call env_esm_hartree_bc3 ( rhog, ehart, aux )
        else if ( esm_bc == 'bc4' ) then
           call env_esm_hartree_bc4 ( rhog, ehart, aux )
        end if

 
     end subroutine env_esm_hartree
!
!-----------------------------------------------------------------------
!--------------ESM HARTREE SUBROUTINE-----------------------------------
!-----------------------------------------------------------------------
SUBROUTINE env_esm_hartree_pbc (rhog, ehart, aux)
  USE env_gvect,    ONLY : ngm
  USE env_fft_base, ONLY : dfftp
  IMPLICIT NONE
  real(DP)    :: ehart             !  Hartree energy
  complex(DP) :: rhog(ngm)   !  n(G)
  complex(DP) :: aux(dfftp%nnr)    !  v_h(G)
  
  stop 'esm_hartree must not be called for esm_bc = pbc'

END SUBROUTINE env_esm_hartree_pbc

SUBROUTINE env_esm_hartree_bc1(rhog, ehart, aux)

  USE env_constants,        ONLY : tpi, fpi, e2
  USE env_gvect,            ONLY : ngm, mill
  USE env_cell_base,        ONLY : omega, alat, tpiba2, at, bg
  USE env_control_flags,    ONLY : gamma_only
  USE env_mp_bands,         ONLY : intra_bgrp_comm
  USE env_mp,               ONLY : env_mp_sum
  USE env_fft_base,         ONLY : dfftp
  USE env_fft_scalar,       ONLY : env_cft_1z

  !
  IMPLICIT NONE
  !
  real(DP)                 :: ehart             !  Hartree energy
  complex(DP)              :: rhog(ngm)         !  n(G)
  complex(DP)              :: aux(dfftp%nnr)    !  v_h(G)   
  !
  !    here the local variables
  !
  integer                  :: k1, k2, k3, iz, ng, n1, n2, n3, nz_r, nz_l, &
                              ng_2d
  real(DP)                 :: t(2), z, z0, gp, gp2, kn, cc0, ss0, L, &
                              z_l, z_r, eh, arg1, arg2
  complex(DP)              :: ci, f1, f2, f3, f4, a0, a1, a2, a3, c_r, c_l, &
                              s_r, s_l, rg3, tmp1, tmp2, tmp3
  complex(DP),allocatable  :: rhog3(:,:), vg(:), vg_r(:), vg3(:,:)

  allocate(rhog3(dfftp%nr3,ngm_2d), vg3(dfftp%nr3,ngm_2d))
!
! Map to FFT mesh (dfftp%nr3,ngm_2d)
  rhog3(:,:)=(0.d0,0.d0)
  do ng=1,ngm
     n1 = mill(1,ng)
     n2 = mill(2,ng)
     ng_2d = imill_2d(n1,n2)
     n3 = mill(3,ng)+1
     IF (n3<1) n3 = n3 + dfftp%nr3
     rg3 = rhog(ng)
     rhog3(n3,ng_2d)=rg3
     if ( gamma_only .and. n1==0 .and. n2==0 ) then
        n3 = -mill(3,ng)+1
        IF (n3<1) n3 = n3 + dfftp%nr3
        rhog3(n3,ng_2d)=CONJG(rg3)
     endif
  enddo
! End mapping
!
  vg3(:,:)=(0.d0,0.d0)
  L=at(3,3)*alat
  z0=L/2.d0
  ci=(0.d0,1.d0)

!****For gp!=0 case ********************
  allocate(vg(dfftp%nr3),vg_r(dfftp%nr3))
  do ng_2d = 1, ngm_2d
     k1 = mill_2d(1,ng_2d)
     k2 = mill_2d(2,ng_2d)
     if(k1==0.and.k2==0) cycle
     t(1:2) = k1 * bg (1:2, 1) + k2 * bg (1:2, 2)
     gp2 = sum( t(:) * t(:) ) * tpiba2
     gp=sqrt(gp2)
     tmp1=(0.d0,0.d0); tmp2=(0.d0,0.d0);
     vg(:)=(0.d0,0.d0)
     do iz=1, dfftp%nr3
        if(iz<=dfftp%nr3/2) kn=dble(iz-1)*tpi/L
        if(iz> dfftp%nr3/2) kn=dble(iz-1-dfftp%nr3)*tpi/L
        cc0=cos(kn*z0)
        ss0=sin(kn*z0)
        rg3=rhog3(iz,ng_2d)
        ! bc1
        vg(iz)=fpi*rg3/(gp**2+kn**2)
        tmp1=tmp1+rg3*(cc0+ci*ss0)/(gp-ci*kn)
        tmp2=tmp2+rg3*(cc0-ci*ss0)/(gp+ci*kn)
     enddo
     vg3(:,ng_2d)=vg(:)

     ! real part
     vg_r(:) = (0.d0,0.d0)
     do iz=1,dfftp%nr3
        k3=iz-1
        if (k3.gt.dfftp%nr3/2) k3=iz-dfftp%nr3-1
        z=dble(k3)/dble(dfftp%nr3)*L
        ! bc1
        arg1= gp*(z-z0)
        arg2=-gp*(z+z0)
        vg_r(iz)=-tpi/gp*(exp(arg1)*tmp1+exp(arg2)*tmp2)
     enddo
     
     call env_cft_1z(vg_r,1,dfftp%nr3,dfftp%nr3,-1,vg)
     vg3(:,ng_2d)=(vg3(:,ng_2d)+vg(:))*e2 ! factor e2: hartree -> Ry.
  enddo
  deallocate(vg,vg_r)

!****For gp=0 case ********************
  ng_2d = imill_2d(0,0)
  if( ng_2d > 0 ) then
     allocate(vg(dfftp%nr3),vg_r(dfftp%nr3))
     tmp1=(0.d0,0.d0); tmp2=(0.d0,0.d0); tmp3=(0.d0,0.d0);
     vg(:)=(0.d0,0.d0);
     rg3=rhog3(1,ng_2d)
     vg(1)=-tpi*z0**2*rg3
     do iz=2,dfftp%nr3
        if(iz<=dfftp%nr3/2) kn=dble(iz-1)*tpi/L
        if(iz> dfftp%nr3/2) kn=dble(iz-1-dfftp%nr3)*tpi/L
        rg3=rhog3(iz,ng_2d)
        cc0=cos(kn*z0)
        ss0=sin(kn*z0)
        tmp1=tmp1+rg3*ci*(cc0+ci*ss0)/kn
        tmp2=tmp2+rg3*ci*(cc0-ci*ss0)/kn
        tmp3=tmp3+rg3*cc0/kn**2
        vg(iz)=fpi*rg3/(kn**2)
     enddo
     vg3(:,ng_2d)=vg(:)
     
     ! real part
     vg_r(:) = (0.d0,0.d0)
     rg3=rhog3(1,ng_2d)
     do iz=1,dfftp%nr3
        k3=iz-1
        if (k3.gt.dfftp%nr3/2) k3=iz-dfftp%nr3-1
        z=dble(k3)/dble(dfftp%nr3)*L
        ! bc1
        vg_r(iz)=-tpi*z**2*rg3 &
                 -tpi*(z-z0)*tmp1 &
                 -tpi*(z+z0)*tmp2 &
                 -fpi*tmp3
     enddo

     ! start smoothing
     nz_l=dfftp%nr3/2+1+esm_nfit
     nz_r=dfftp%nr3/2+1-esm_nfit
     z_l=dble(nz_l-1)*L/dble(dfftp%nr3)-L
     z_r=dble(nz_r-1)*L/dble(dfftp%nr3)
     f1=-tpi*z_r**2*rg3 &
        -tpi*(z_r-z0)*tmp1 &
        -tpi*(z_r+z0)*tmp2 &
        -fpi*tmp3
     f2=-tpi*z_l**2*rg3 &
        -tpi*(z_l-z0)*tmp1 &
        -tpi*(z_l+z0)*tmp2 &
        -fpi*tmp3
     f3=-fpi*z_r*rg3 &
        -tpi*tmp1 &
        -tpi*tmp2
     f4=-fpi*z_l*rg3 &
        -tpi*tmp1 &
        -tpi*tmp2
     z_r=z_r
     z_l=z_l+L
     a0=(f1*z_l**2*(z_l-3.d0*z_r)+z_r*(f3*z_l**2*(-z_l+z_r) &
          +z_r*(f2*(3.d0*z_l-z_r)+f4*z_l*(-z_l+z_r))))/(z_l-z_r)**3
     a1=(f3*z_l**3+z_l*(6.d0*f1-6.d0*f2+(f3+2.d0*f4)*z_l)*z_r &
          -(2*f3+f4)*z_l*z_r**2-f4*z_r**3)/(z_l-z_r)**3
     a2=(-3*f1*(z_l+z_r)+3.d0*f2*(z_l+z_r)-(z_l-z_r)*(2*f3*z_l &
          +f4*z_l+f3*z_r+2*f4*z_r))/(z_l-z_r)**3
     a3=(2.d0*f1-2.d0*f2+(f3+f4)*(z_l-z_r))/(z_l-z_r)**3
     do iz=nz_r,nz_l
        z=dble(iz-1)/dble(dfftp%nr3)*L
        vg_r(iz)=(a0+a1*z+a2*z**2+a3*z**3)
     enddo
     ! end smoothing

     call env_cft_1z(vg_r,1,dfftp%nr3,dfftp%nr3,-1,vg)
     vg3(:,ng_2d)=(vg3(:,ng_2d)+vg(:))*e2 ! factor e2: hartree -> Ry.

     deallocate (vg,vg_r)
  endif ! if( ng_2d > 0 )

! Hartree Energy
  ehart=0.d0
  eh = 0d0
  do ng_2d = 1, ngm_2d
     k1 = mill_2d(1,ng_2d)
     k2 = mill_2d(2,ng_2d)
     eh = eh + sum( vg3(:,ng_2d)*conjg(rhog3(:,ng_2d)) )
  enddo
  ehart=ehart+eh
  if( gamma_only ) then
     ehart = ehart * 2d0
     ng_2d = imill_2d(0,0)
     if( ng_2d > 0 ) then
        ehart = ehart - sum( vg3(:,ng_2d)*conjg(rhog3(:,ng_2d)) )
     endif
  endif
  ehart = ehart*omega*0.5d0
  !
  call env_mp_sum( ehart, intra_bgrp_comm )
  !
! Map to FFT mesh (dfftp%nrx)
  aux=0.0d0
  do ng=1,ngm
     n1 = mill(1,ng)
     n2 = mill(2,ng)
     ng_2d = imill_2d(n1,n2)
     n3 = mill(3,ng) + 1
     if (n3<1) n3 = n3 + dfftp%nr3
     aux(dfftp%nl(ng))= aux(dfftp%nl(ng)) + vg3(n3,ng_2d)
     if (gamma_only) then
        aux(dfftp%nlm(ng))=CONJG(aux(dfftp%nl(ng)))
     endif
  enddo

  deallocate (rhog3,vg3)

  RETURN
END SUBROUTINE env_esm_hartree_bc1

SUBROUTINE env_esm_hartree_bc2 (rhog, ehart, aux)

  USE env_constants,        ONLY : tpi, fpi, e2
  USE env_gvect,            ONLY : ngm, mill
  USE env_cell_base,        ONLY : omega, alat, tpiba2, at, bg
  USE env_control_flags,    ONLY : gamma_only
  USE env_mp_bands,         ONLY : intra_bgrp_comm
  USE env_mp,               ONLY : env_mp_sum
  USE env_fft_base,         ONLY : dfftp
  USE env_fft_scalar,       ONLY : env_cft_1z
  !
  IMPLICIT NONE
  !
  real(DP)                 :: ehart             !  Hartree energy
  complex(DP)              :: rhog(ngm)         !  n(G)
  complex(DP)              :: aux(dfftp%nnr)    !  v_h(G)
  !
  !    here the local variables
  !
  integer                  :: k1, k2, k3, iz, ng, n1, n2, n3, nz_r, nz_l, &
                              ng_2d
  real(DP)                 :: t(2), z, z0, z1, gp, gp2, kn, cc0, ss0, L, &
                              z_l, z_r, eh, arg1, arg2, arg3, arg4, arg5
  complex(DP)              :: ci, f1, f2, f3, f4, a0, a1, a2, a3, c_r, c_l, &
                              s_r, s_l, rg3, tmp, tmp1, tmp2, tmp3, tmp4
  complex(DP),allocatable  :: rhog3(:,:), vg(:), vg_r(:), vg3(:,:)

  allocate(rhog3(dfftp%nr3,ngm_2d), vg3(dfftp%nr3,ngm_2d))
!
! Map to FFT mesh (dfftp%nr3,ngm_2d)
  rhog3(:,:)=(0.d0,0.d0)
  do ng=1,ngm
     n1 = mill(1,ng)
     n2 = mill(2,ng)
     ng_2d = imill_2d(n1,n2)
     n3 = mill(3,ng)+1
     IF (n3<1) n3 = n3 + dfftp%nr3    
     rg3 = rhog(ng)
     rhog3(n3,ng_2d)=rg3
     if ( gamma_only .and. n1==0 .and. n2==0 ) then
        n3 = -mill(3,ng)+1
        IF (n3<1) n3 = n3 + dfftp%nr3
        rhog3(n3,ng_2d)=CONJG(rg3)
     endif
  enddo
! End mapping
!
  vg3(:,:)=(0.d0,0.d0)
  L=at(3,3)*alat
  z0=L/2.d0
  z1=z0+esm_w
  ci=(0.d0,1.d0)

!****For gp!=0 case ********************
  allocate(vg(dfftp%nr3),vg_r(dfftp%nr3))
  do ng_2d = 1, ngm_2d
     k1 = mill_2d(1,ng_2d)
     k2 = mill_2d(2,ng_2d)
     if(k1==0.and.k2==0) cycle
     t(1:2) = k1 * bg (1:2, 1) + k2 * bg (1:2, 2)
     gp2 = sum( t(:) * t(:) ) * tpiba2
     gp=sqrt(gp2)
     tmp1=(0.d0,0.d0); tmp2=(0.d0,0.d0);
     vg(:)=(0.d0,0.d0)
     do iz=1, dfftp%nr3
        if(iz<=dfftp%nr3/2) kn=dble(iz-1)*tpi/L
        if(iz> dfftp%nr3/2) kn=dble(iz-1-dfftp%nr3)*tpi/L
        cc0=cos(kn*z0)
        ss0=sin(kn*z0)
        rg3=rhog3(iz,ng_2d)
        ! bc2
        arg1= gp*(z1-z0)
        arg2=-gp*(z1-z0)
        vg(iz)=fpi*rg3/(gp**2+kn**2)
        tmp=((gp+ci*kn)*exp(arg1)+(gp-ci*kn)*exp(arg2))/(2.d0*gp)
        tmp1=tmp1+rg3*(cc0+ci*ss0)/(gp**2+kn**2)*tmp
        tmp=((gp-ci*kn)*exp(arg1)+(gp+ci*kn)*exp(arg2))/(2.d0*gp)
        tmp2=tmp2+rg3*(cc0-ci*ss0)/(gp**2+kn**2)*tmp
     enddo
     vg3(:,ng_2d)=vg(:)
     
     ! real part
     vg_r(:)=(0.d0,0.d0)
     do iz=1,dfftp%nr3
        k3=iz-1
        if (k3.gt.dfftp%nr3/2) k3=iz-dfftp%nr3-1
        z=dble(k3)/dble(dfftp%nr3)*L
        ! bc2
        arg1= gp*(z-z1)
        arg2=-gp*(z+z1)
        arg3= gp*(z-3.d0*z1)
        arg4=-gp*(z+3.d0*z1)
        arg5=-4.d0*gp*z1
        vg_r(iz)=-fpi*(exp(arg1)-exp(arg4))*tmp1/(1.d0-exp(arg5)) &
           +fpi*(exp(arg3)-exp(arg2))*tmp2/(1.d0-exp(arg5))
     enddo
     call env_cft_1z(vg_r,1,dfftp%nr3,dfftp%nr3,-1,vg)
     vg3(:,ng_2d)=(vg3(:,ng_2d)+vg(:))*e2 ! factor e2: hartree -> Ry.
  enddo
  deallocate(vg,vg_r)

!****For gp=0 case ********************
  ng_2d = imill_2d(0,0)
  if( ng_2d > 0 ) then
     allocate(vg(dfftp%nr3),vg_r(dfftp%nr3))
     tmp1=(0.d0,0.d0); tmp2=(0.d0,0.d0); tmp3=(0.d0,0.d0); tmp4=(0.d0,0.d0)
     vg(:)=(0.d0,0.d0);
     rg3=rhog3(1,ng_2d)
     vg(1)= tpi*(2.d0*z1-z0)*z0*rg3
     do iz=2,dfftp%nr3
        if(iz<=dfftp%nr3/2) kn=dble(iz-1)*tpi/L
        if(iz> dfftp%nr3/2) kn=dble(iz-1-dfftp%nr3)*tpi/L
        rg3=rhog3(iz,ng_2d)
        cc0=cos(kn*z0)
        ss0=sin(kn*z0)
        tmp1=tmp1+rg3*(cc0+ci*ss0)/kn**2
        tmp2=tmp2+rg3*(cc0-ci*ss0)/kn**2
        tmp3=tmp3+rg3*ci*cc0/kn
        tmp4=tmp4+rg3*ss0/kn
        vg(iz)=fpi*rg3/(kn**2)
     enddo
     vg3(:,ng_2d)=vg(:)
     
     ! real part
     vg_r(:) = (0.d0,0.d0)
     rg3=rhog3(1,ng_2d)
     do iz=1,dfftp%nr3
        k3=iz-1
        if (k3.gt.dfftp%nr3/2) k3=iz-dfftp%nr3-1
        z=dble(k3)/dble(dfftp%nr3)*L
        vg_r(iz)=-tpi*z**2*rg3 &
                 -tpi*(z+z1)*tmp1/z1 &
                 +tpi*(z-z1)*tmp2/z1 &
                 -fpi*z*(z1-z0)/z1*tmp3 &
                 +fpi*(z1-z0)*tmp4
     enddo

     ! start smoothing
     nz_l=dfftp%nr3/2+1+esm_nfit
     nz_r=dfftp%nr3/2+1-esm_nfit
     z_l=dble(nz_l-1)*L/dble(dfftp%nr3)-L
     z_r=dble(nz_r-1)*L/dble(dfftp%nr3)
     f1=-tpi*z_r**2*rg3 &
        -tpi*(z_r+z1)*tmp1/z1 &
        +tpi*(z_r-z1)*tmp2/z1 &
        -fpi*z_r*(z1-z0)/z1*tmp3 &
        +fpi*(z1-z0)*tmp4
     f2=-tpi*z_l**2*rg3 &
        -tpi*(z_l+z1)*tmp1/z1 &
        +tpi*(z_l-z1)*tmp2/z1 &
        -fpi*z_l*(z1-z0)/z1*tmp3 &
        +fpi*(z1-z0)*tmp4
     f3=-fpi*z_r*rg3 &
        -tpi*tmp1/z1 &
        +tpi*tmp2/z1 &
        -fpi*(z1-z0)/z1*tmp3
     f4=-fpi*z_l*rg3 &
        -tpi*tmp1/z1 &
        +tpi*tmp2/z1 &
        -fpi*(z1-z0)/z1*tmp3
     z_r=z_r
     z_l=z_l+L
     a0=(f1*z_l**2*(z_l-3.d0*z_r)+z_r*(f3*z_l**2*(-z_l+z_r) &
          +z_r*(f2*(3.d0*z_l-z_r)+f4*z_l*(-z_l+z_r))))/(z_l-z_r)**3
     a1=(f3*z_l**3+z_l*(6.d0*f1-6.d0*f2+(f3+2.d0*f4)*z_l)*z_r &
          -(2*f3+f4)*z_l*z_r**2-f4*z_r**3)/(z_l-z_r)**3
     a2=(-3*f1*(z_l+z_r)+3.d0*f2*(z_l+z_r)-(z_l-z_r)*(2*f3*z_l &
          +f4*z_l+f3*z_r+2*f4*z_r))/(z_l-z_r)**3
     a3=(2.d0*f1-2.d0*f2+(f3+f4)*(z_l-z_r))/(z_l-z_r)**3
     do iz=nz_r,nz_l
        z=dble(iz-1)/dble(dfftp%nr3)*L
        vg_r(iz)=(a0+a1*z+a2*z**2+a3*z**3)
     enddo
     ! end smoothing
     
     call env_cft_1z(vg_r,1,dfftp%nr3,dfftp%nr3,-1,vg)
     vg3(:,ng_2d)=(vg3(:,ng_2d)+vg(:))*e2 ! factor e2: hartree -> Ry.

     deallocate (vg,vg_r)
  endif ! if( ng_2d > 0 )

! Hartree Energy
  ehart=0.d0
  eh = 0d0
  do ng_2d = 1, ngm_2d
     k1 = mill_2d(1,ng_2d)
     k2 = mill_2d(2,ng_2d)
     eh = eh + sum( vg3(:,ng_2d)*conjg(rhog3(:,ng_2d)) )
  enddo
  ehart=ehart+eh
  if( gamma_only ) then
     ehart = ehart * 2d0
     ng_2d = imill_2d(0,0)
     if( ng_2d > 0 ) then
        ehart = ehart - sum( vg3(:,ng_2d)*conjg(rhog3(:,ng_2d)) )
     endif
  endif
  ehart = ehart*omega*0.5d0
  !
  call env_mp_sum( ehart, intra_bgrp_comm )
  !
! Map to FFT mesh (dfftp%nrx)
  aux=0.0d0
  do ng=1,ngm
     n1 = mill(1,ng)
     n2 = mill(2,ng)
     ng_2d = imill_2d(n1,n2)
     n3 = mill(3,ng) + 1
     if (n3<1) n3 = n3 + dfftp%nr3
     aux(dfftp%nl(ng))= aux(dfftp%nl(ng)) + vg3(n3,ng_2d)
     if (gamma_only) then
        aux(dfftp%nlm(ng))=CONJG(aux(dfftp%nl(ng)))
     endif
  enddo

  deallocate (rhog3,vg3)

  RETURN
END SUBROUTINE env_esm_hartree_bc2

SUBROUTINE env_esm_hartree_bc3 (rhog, ehart, aux)

  USE env_constants,        ONLY : tpi, fpi, e2
  USE env_gvect,            ONLY : ngm, mill
  USE env_cell_base,        ONLY : omega, alat, tpiba2, at, bg
  USE env_control_flags,    ONLY : gamma_only
  USE env_mp_bands,         ONLY : intra_bgrp_comm
  USE env_mp,               ONLY : env_mp_sum
  USE env_fft_base,         ONLY : dfftp
  USE env_fft_scalar,       ONLY : env_cft_1z
  !
  IMPLICIT NONE
  !
  real(DP)                 :: ehart             !  Hartree energy
  complex(DP)              :: rhog(ngm)         !  n(G)
  complex(DP)              :: aux(dfftp%nnr)    !  v_h(G)
  !
  !    here the local variables
  !
  integer                  :: k1, k2, k3, iz, ng, n1, n2, n3, nz_r, nz_l, &
                              ng_2d
  real(DP)                 :: t(2), z, z0, z1, gp, gp2, kn, cc0, ss0, L, &
                              z_l, z_r, eh, arg1, arg2, arg3
  complex(DP)              :: ci, f1, f2, f3, f4, a0, a1, a2, a3, c_r, c_l, &
                              s_r, s_l, rg3, tmp, tmp1, tmp2, tmp3
  complex(DP),allocatable  :: rhog3(:,:), vg(:), vg_r(:), vg3(:,:)

  allocate(rhog3(dfftp%nr3,ngm_2d), vg3(dfftp%nr3,ngm_2d))
!
! Map to FFT mesh (dfftp%nr3,ngm_2d)
  rhog3(:,:)=(0.d0,0.d0)
  do ng=1,ngm
     n1 = mill(1,ng)
     n2 = mill(2,ng)
     ng_2d = imill_2d(n1,n2)
     n3 = mill(3,ng)+1
     IF (n3<1) n3 = n3 + dfftp%nr3    
     rg3 = rhog(ng)
     rhog3(n3,ng_2d)=rg3
     if ( gamma_only .and. n1==0 .and. n2==0 ) then
        n3 = -mill(3,ng)+1
        IF (n3<1) n3 = n3 + dfftp%nr3
        rhog3(n3,ng_2d)=CONJG(rg3)
     endif
  enddo
! End mapping
!
  vg3(:,:)=(0.d0,0.d0)
  L=at(3,3)*alat
  z0=L/2.d0
  z1=z0+esm_w
  ci=(0.d0,1.d0)

!****For gp!=0 case ********************
  allocate(vg(dfftp%nr3),vg_r(dfftp%nr3))
  do ng_2d = 1, ngm_2d
     k1 = mill_2d(1,ng_2d)
     k2 = mill_2d(2,ng_2d)
     if(k1==0.and.k2==0) cycle
     t(1:2) = k1 * bg (1:2, 1) + k2 * bg (1:2, 2)
     gp2 = sum( t(:) * t(:) ) * tpiba2
     gp=sqrt(gp2)
     tmp1=(0.d0,0.d0); tmp2=(0.d0,0.d0);
     vg(:)=(0.d0,0.d0)
     do iz=1, dfftp%nr3
        if(iz<=dfftp%nr3/2) kn=dble(iz-1)*tpi/L
        if(iz> dfftp%nr3/2) kn=dble(iz-1-dfftp%nr3)*tpi/L
        cc0=cos(kn*z0)
        ss0=sin(kn*z0)
        rg3=rhog3(iz,ng_2d)
        ! bc3
        arg1= gp*(z1-z0)
        arg2=-gp*(z1-z0)
        vg(iz)=fpi*rg3/(gp**2+kn**2)
        tmp=((gp+ci*kn)*exp(arg1)+(gp-ci*kn)*exp(arg2))/(2.d0*gp)
        tmp1=tmp1+rg3*(cc0+ci*ss0)/(gp**2+kn**2)*tmp
        tmp=(gp-ci*kn)/gp
        tmp2=tmp2+rg3*(cc0-ci*ss0)/(gp**2+kn**2)*tmp
     enddo
     vg3(:,ng_2d)=vg(:)

     ! real part
     vg_r(:)=(0.d0,0.d0)
     do iz=1,dfftp%nr3
        k3=iz-1
        if (k3.gt.dfftp%nr3/2) k3=iz-dfftp%nr3-1
        z=dble(k3)/dble(dfftp%nr3)*L
        ! bc3
        arg1= gp*(z-z1)
        arg2=-gp*(z+z0)
        arg3= gp*(z-z0-2.d0*z1)
        vg_r(iz)=-fpi*exp(arg1)*tmp1+tpi*(exp(arg3)-exp(arg2))*tmp2
     enddo
     
     call env_cft_1z(vg_r,1,dfftp%nr3,dfftp%nr3,-1,vg)
     vg3(:,ng_2d)=(vg3(:,ng_2d)+vg(:))*e2 ! factor e2: hartree -> Ry.
  enddo
  deallocate(vg,vg_r)

!****For gp=0 case ********************
  ng_2d = imill_2d(0,0)
  if( ng_2d > 0 ) then
     allocate(vg(dfftp%nr3),vg_r(dfftp%nr3))
     tmp1=(0.d0,0.d0); tmp2=(0.d0,0.d0); tmp3=(0.d0,0.d0);
     vg(:)=(0.d0,0.d0);
     rg3=rhog3(1,ng_2d)
     vg(1)= tpi*(4.d0*z1-z0)*z0*rg3
     do iz=2,dfftp%nr3
        if(iz<=dfftp%nr3/2) kn=dble(iz-1)*tpi/L
        if(iz> dfftp%nr3/2) kn=dble(iz-1-dfftp%nr3)*tpi/L
        rg3=rhog3(iz,ng_2d)
        cc0=cos(kn*z0)
        ss0=sin(kn*z0)
        tmp1=tmp1+rg3*(cc0+ci*ss0)/kn**2
        tmp2=tmp2+rg3*(cc0-ci*ss0)/kn
        tmp3=tmp3+rg3*(cc0+ci*ss0)/kn
        vg(iz)=fpi*rg3/(kn**2)
     enddo
     vg3(:,ng_2d)=vg(:)
     
     ! real part
     vg_r(:) = (0.d0,0.d0) 
     rg3=rhog3(1,ng_2d)
     do iz=1,dfftp%nr3
        k3=iz-1
        if (k3.gt.dfftp%nr3/2) k3=iz-dfftp%nr3-1
        z=dble(k3)/dble(dfftp%nr3)*L
        vg_r(iz)=-tpi*(z**2+2.d0*z*z0)*rg3 &
                 -fpi*tmp1 &
                 -fpi*ci*(z-z1)*tmp2 &
                 -fpi*ci*(z1-z0)*tmp3
     enddo

     ! start smoothing
     nz_l=dfftp%nr3/2+1+esm_nfit
     nz_r=dfftp%nr3/2+1-esm_nfit
     z_l=dble(nz_l-1)*L/dble(dfftp%nr3)-L
     z_r=dble(nz_r-1)*L/dble(dfftp%nr3)
     f1=-tpi*(z_r**2+2.d0*z_r*z0)*rg3 &
        -fpi*tmp1 &
        -fpi*ci*(z_r-z1)*tmp2 &
        -fpi*ci*(z1 -z0)*tmp3
     f2=-tpi*(z_l**2+2.d0*z_l*z0)*rg3 &
        -fpi*tmp1 &
        -fpi*ci*(z_l-z1)*tmp2 &
        -fpi*ci*(z1 -z0)*tmp3
     f3=-tpi*(2.d0*z_r+2.d0*z0)*rg3 &
        -fpi*ci*tmp2
     f4=-tpi*(2.d0*z_l+2.d0*z0)*rg3 &
        -fpi*ci*tmp2
     z_r=z_r
     z_l=z_l+L
     a0=(f1*z_l**2*(z_l-3.d0*z_r)+z_r*(f3*z_l**2*(-z_l+z_r) &
          +z_r*(f2*(3.d0*z_l-z_r)+f4*z_l*(-z_l+z_r))))/(z_l-z_r)**3
     a1=(f3*z_l**3+z_l*(6.d0*f1-6.d0*f2+(f3+2.d0*f4)*z_l)*z_r &
          -(2*f3+f4)*z_l*z_r**2-f4*z_r**3)/(z_l-z_r)**3
     a2=(-3*f1*(z_l+z_r)+3.d0*f2*(z_l+z_r)-(z_l-z_r)*(2*f3*z_l &
          +f4*z_l+f3*z_r+2*f4*z_r))/(z_l-z_r)**3
     a3=(2.d0*f1-2.d0*f2+(f3+f4)*(z_l-z_r))/(z_l-z_r)**3
     do iz=nz_r,nz_l
        z=dble(iz-1)/dble(dfftp%nr3)*L
        vg_r(iz)=(a0+a1*z+a2*z**2+a3*z**3)
     enddo
     ! end smoothing

     call env_cft_1z(vg_r,1,dfftp%nr3,dfftp%nr3,-1,vg)
     vg3(:,ng_2d)=(vg3(:,ng_2d)+vg(:))*e2 ! factor e2: hartree -> Ry.

     deallocate (vg,vg_r)
  endif ! if( ng_2d > 0 )

! Hartree Energy
  ehart=0.d0
  eh = 0d0
  do ng_2d = 1, ngm_2d
     k1 = mill_2d(1,ng_2d)
     k2 = mill_2d(2,ng_2d)
     eh = eh + sum( vg3(:,ng_2d)*conjg(rhog3(:,ng_2d)) )
  enddo
  ehart=ehart+eh
  if( gamma_only ) then
     ehart = ehart * 2d0
     ng_2d = imill_2d(0,0)
     if( ng_2d > 0 ) then
        ehart = ehart - sum( vg3(:,ng_2d)*conjg(rhog3(:,ng_2d)) )
     endif
  endif
  ehart = ehart*omega*0.5d0
  !
  call env_mp_sum( ehart, intra_bgrp_comm )
  !
! Map to FFT mesh (dfftp%nrx)
  aux=0.0d0
  do ng=1,ngm
     n1 = mill(1,ng)
     n2 = mill(2,ng)
     ng_2d = imill_2d(n1,n2)
     n3 = mill(3,ng) + 1
     if (n3<1) n3 = n3 + dfftp%nr3
     aux(dfftp%nl(ng))= aux(dfftp%nl(ng)) + vg3(n3,ng_2d)
     if (gamma_only) then
        aux(dfftp%nlm(ng))=CONJG(aux(dfftp%nl(ng)))
     endif
  enddo

  deallocate (rhog3,vg3)

  RETURN
END SUBROUTINE env_esm_hartree_bc3

SUBROUTINE env_esm_hartree_bc4 (rhog, ehart, aux)

  USE env_constants,        ONLY : pi, tpi, fpi, e2
  USE env_gvect,            ONLY : ngm, mill
  USE env_cell_base,        ONLY : omega, alat, tpiba2, at, bg
  USE env_control_flags,    ONLY : gamma_only
  USE env_mp_bands,         ONLY : intra_bgrp_comm
  USE env_mp,               ONLY : env_mp_sum
  USE env_fft_base,         ONLY : dfftp
  USE env_fft_scalar,       ONLY : env_cft_1z
  !
  IMPLICIT NONE
  !
  real(DP)                 :: ehart             !  Hartree energy
  complex(DP)              :: rhog(ngm)         !  n(G)
  complex(DP)              :: aux(dfftp%nnr)    !  v_h(G)
  !
  !    here the local variables
  !
  integer                  :: k1, k2, k3, iz, ng, n1, n2, n3, nz_r, nz_l, &
                              ng_2d
  real(DP)                 :: t(2), z, z0, z1, gp, gp2, kn, cc0, ss0, L, &
                              z_l, z_r, eh, aaa, cc1, ss1, alpha, beta, &
                              chi, xi, kappa, lambda, arg1, arg2, arg3, &
                              arg4, argr1, argr2, argr3, argr4, argr5
  complex(DP)              :: ci, f1, f2, f3, f4, a0, a1, a2, a3, c_r, c_l, &
                              s_r, s_l, rg3, tmp, tmp1, tmp2, tmp3, tmp4, &
                              tmpr1, tmpr2, tmpr3, tmpr4
  complex(DP),allocatable  :: rhog3(:,:), vg(:), vg_r(:), vg3(:,:), vr(:), &
                              vr_r(:)

  allocate(rhog3(dfftp%nr3,ngm_2d), vg3(dfftp%nr3,ngm_2d))
!
! Map to FFT mesh (dfftp%nr3,ngm_2d)
  rhog3(:,:)=(0.d0,0.d0)
  do ng=1,ngm
     n1 = mill(1,ng)
     n2 = mill(2,ng)
     ng_2d = imill_2d(n1,n2)
     n3 = mill(3,ng)+1
     IF (n3<1) n3 = n3 + dfftp%nr3    
     rg3 = rhog(ng)
     rhog3(n3,ng_2d)=rg3
     if ( gamma_only .and. n1==0 .and. n2==0 ) then
        n3 = -mill(3,ng)+1
        IF (n3<1) n3 = n3 + dfftp%nr3
        rhog3(n3,ng_2d)=CONJG(rg3)
     endif
  enddo
! End mapping
!
  vg3(:,:)=(0.d0,0.d0)
  L=at(3,3)*alat
  z0=L/2.d0
  z1=z0+esm_w
  aaa=esm_a
  ci=(0.d0,1.d0)

!****For gp!=0 case ********************
  allocate(vg(dfftp%nr3),vg_r(dfftp%nr3),vr(dfftp%nr3),vr_r(dfftp%nr3))
  do ng_2d = 1, ngm_2d
     k1 = mill_2d(1,ng_2d)
     k2 = mill_2d(2,ng_2d)
     if(k1==0.and.k2==0) cycle
     t(1:2) = k1 * bg (1:2, 1) + k2 * bg (1:2, 2)
     gp2 = sum( t(:) * t(:) ) * tpiba2
     gp=sqrt(gp2)
     tmp1=(0.d0,0.d0); tmp2=(0.d0,0.d0); tmp3=(0.d0,0.d0)
     tmp4=(0.d0,0.d0); tmpr1=(0.d0,0.d0); tmpr2=(0.d0,0.d0)
     tmpr3=(0.d0,0.d0); tmpr4=(0.d0,0.d0)
     vr(:)=(0.d0,0.d0); vg(:)=(0.d0,0.d0)
     do iz=1, dfftp%nr3
        if(iz<=dfftp%nr3/2) kn=dble(iz-1)*tpi/L
        if(iz> dfftp%nr3/2) kn=dble(iz-1-dfftp%nr3)*tpi/L
        cc0=cos(kn*z0)
        ss0=sin(kn*z0)
        rg3=rhog3(iz,ng_2d)
        ! bc4
        vg(iz)=fpi*rg3/(gp**2+kn**2)
        vr(iz)=fpi*rg3/(gp**2+kn**2+ci*aaa*kn)
        cc1=cos(kn*z1)
        ss1=sin(kn*z1)
        alpha=aaa+gp+sqrt(aaa**2+gp**2)
        beta =aaa+gp-sqrt(aaa**2+gp**2)
        kappa=aaa-gp+sqrt(aaa**2+gp**2)
        xi   =aaa   +sqrt(aaa**2+gp**2)
        chi  =aaa   -sqrt(aaa**2+gp**2)
        lambda=      sqrt(aaa**2+gp**2)
        tmp1=tmp1+rg3*(cc0+ci*ss0)/(xi-ci*kn)/alpha
        tmp2=tmp2+rg3*(cc0-ci*ss0)/(gp+ci*kn)/gp
        tmp3=tmp3+rg3*kappa/alpha*(cc0-ci*ss0)/(gp+ci*kn)/gp
        tmp4=tmp4+rg3*kappa*(cc1+ci*ss1)/(xi-ci*kn)/(gp**2+kn**2)
        tmpr1=tmpr1+rg3*(cc0-ci*ss0)/(gp+ci*kn)/alpha
        tmpr2=tmpr2+rg3*(cc0+ci*ss0)/(xi-ci*kn)/lambda
        tmpr3=tmpr3+rg3*beta/alpha*(cc0+ci*ss0)/(xi-ci*kn)/lambda
        tmpr4=tmpr4+rg3*beta*(cc1+ci*ss1)/(gp+ci*kn) &
           /(gp**2+kn**2+ci*2.d0*aaa*kn)
     enddo
     
     call env_cft_1z(vg,1,dfftp%nr3,dfftp%nr3,1,vg_r)
     ! bc4
     CALL env_cft_1z(vr,1,dfftp%nr3,dfftp%nr3,1,vr_r)
     
     do iz=1,dfftp%nr3
        k3=iz-1
        if (k3.gt.dfftp%nr3/2) k3=iz-dfftp%nr3-1
        z=dble(k3)/dble(dfftp%nr3)*L
        ! bc4
        arg1= gp*(z-z1)-xi*(z0-z1)
        arg2=-gp*(z+z0)
        arg3=-gp*(z0+z1)+gp*(z-z1)
        arg4= gp*(z-z1)
        argr1=-gp*(z0+z1)-xi*(z-z1)
        argr2=-xi*(z0-z1)-chi*(z-z1)
        argr3=-xi*(z-z1)-xi*(z0-z1)
        argr4=-xi*(z-z1)
        argr5=-2.d0*aaa*(z-z1)
        if (z < z1) then
           vg_r(iz) = vg_r(iz)-fpi*exp(arg1)*tmp1-tpi*exp(arg2)*tmp2 &
              +tpi*exp(arg3)*tmp3-fpi*exp(arg4)*tmp4
        else
           vg_r(iz) = vr_r(iz)*exp(argr5) &
              -fpi*exp(argr1)*tmpr1-tpi*exp(argr2)*tmpr2 &
              +tpi*exp(argr3)*tmpr3-fpi*exp(argr4)*tmpr4
        endif
     enddo
     
     call env_cft_1z(vg_r,1,dfftp%nr3,dfftp%nr3,-1,vg)
     vg3(:,ng_2d)=vg(:)*e2 ! factor e2: hartree -> Ry.
  enddo
  deallocate(vg,vg_r,vr,vr_r)

!****For gp=0 case ********************
  ng_2d = imill_2d(0,0)
  if( ng_2d > 0 ) then
     allocate(vg(dfftp%nr3),vg_r(dfftp%nr3),vr(dfftp%nr3),vr_r(dfftp%nr3))
     tmp1=(0.d0,0.d0); tmp2=(0.d0,0.d0); tmp3=(0.d0,0.d0); tmp4=(0.d0,0.d0)
     vg(:)=(0.d0,0.d0); vr(:)=(0.d0,0.d0)
     !for smoothing
     f1=(0.d0,0.d0); f2=(0.d0,0.d0); f3=(0.d0,0.d0); f4=(0.d0,0.d0)
     nz_l=dfftp%nr3/2+1+esm_nfit
     nz_r=dfftp%nr3/2+1-esm_nfit
     z_l=dble(nz_l-1)*L/dble(dfftp%nr3)-L
     z_r=dble(nz_r-1)*L/dble(dfftp%nr3)
     !
     rg3=rhog3(1,ng_2d)
     ! bc4
     arg1=-2.d0*aaa*(z0-z1)
     vg(1)= tpi*((z0+z1)/aaa+2.d0*z0*z1+z1**2)*rg3 &
        -pi *(exp(arg1)-1.d0)/aaa**2*rg3
     vr(1)= tpi*(z0+0.5d0/aaa)/aaa*rg3

     do iz=2,dfftp%nr3
        if(iz<=dfftp%nr3/2) kn=dble(iz-1)*tpi/L
        if(iz> dfftp%nr3/2) kn=dble(iz-1-dfftp%nr3)*tpi/L
        rg3=rhog3(iz,ng_2d)
        ! bc4
        cc0= cos(kn*z0)
        ss0= sin(kn*z0)
        cc1=cos(kn*z1)
        ss1=sin(kn*z1)
        tmp1=tmp1+rg3*(cc1+ci*ss1)/(2.d0*aaa-ci*kn)/kn**2
        tmp2=tmp2+rg3*(cc0-ci*ss0)/kn
        tmp3=tmp3+rg3*(cc0+ci*ss0)/(2.d0*aaa-ci*kn)
        tmp4=tmp4+(0.d0,0.d0)

        vg(iz)=fpi*rg3/(kn**2)
        ! bc4
        vr(iz)=fpi*rg3/(kn**2+ci*2.d0*aaa*kn)

        !for smoothing
        c_r=cos(kn*z_r)
        s_r=sin(kn*z_r)
        c_l=cos(kn*z_l)
        s_l=sin(kn*z_l)
        ! bc4
        f1=f1+fpi*   rg3*(c_r+ci*s_r)/(kn**2+ci*2.d0*aaa*kn)
        f2=f2+fpi*   rg3*(c_l+ci*s_l)/kn**2
        f3=f3+fpi*ci*rg3*(c_r+ci*s_r)/kn
        f4=f4+fpi*ci*rg3*(c_l+ci*s_l)/kn
        !
     enddo
     
     call env_cft_1z(vg,1,dfftp%nr3,dfftp%nr3,1,vg_r)
     ! bc4
     call env_cft_1z(vr,1,dfftp%nr3,dfftp%nr3,1,vr_r)
     
     rg3=rhog3(1,ng_2d)
     do iz=1,dfftp%nr3
        k3=iz-1
        if (k3.gt.dfftp%nr3/2) k3=iz-dfftp%nr3-1
        z=dble(k3)/dble(dfftp%nr3)*L
        ! bc4
        arg1=-2.d0*aaa*(z0-z1)
        arg2=-2.d0*aaa*(z-z1)
        if (z < z1) then
           vg_r(iz)=vg_r(iz) &
              -fpi*2.d0*aaa*tmp1 &
              -tpi*ci*(2.d0*(z -z1)-1.d0/aaa)*tmp2 &
              -tpi*exp(arg1)/aaa*tmp3 &
              -tpi*z*(z+2.d0*z0)*rg3
        else
           vg_r(iz)=vr_r(iz)*exp(arg2) &
              +tpi*ci*exp(arg2)/aaa*tmp2 &
              -tpi*exp(arg1)/aaa*tmp3 &
              +tpi*exp(arg2)*z/aaa*rg3 &
                    -pi *exp(arg1)/aaa**2*rg3
        endif
     enddo

     !for smoothing
     ! bc4
     arg1=-2.d0*aaa*(z0-z1)
     arg2=-2.d0*aaa*(z_r-z1)
     f1=f1+tpi*(z0+0.5d0/aaa)/aaa*rg3
     f1=f1*exp(arg2) &
        +tpi*ci*exp(arg2)/aaa*tmp2 &
        -tpi*exp(arg1)/aaa*tmp3 &
        +tpi*exp(arg2)*z_r/aaa*rg3 &
        -pi *exp(arg1)/aaa**2*rg3
     f2=f2+tpi*((z0+z1)/aaa+2.d0*z0*z1+z1**2)*rg3 &
        -pi *(exp(arg1)-1.d0)/aaa**2*rg3
     f2=f2 &
        -fpi*2.d0*aaa*tmp1 &
        -tpi*ci*(2.d0*(z_l-z1)-1.d0/aaa)*tmp2 &
        -tpi*exp(arg1)/aaa*tmp3 &
        -tpi*z_l*(z_l+2.d0*z0)*rg3
     f3=f3*exp(arg2) &
        -fpi*ci*exp(arg2)*tmp2 &
        -fpi*(z_r+z0)*exp(arg2)*rg3 
     f4=f4-fpi*ci*tmp2-fpi*(z_l+z0)*rg3   
     ! for smoothing
     !factor e2 will be multiplied later (at vg3 <= vg)
     !f1=f1*e2; f2=f2*e2; f3=f3*e2; f4=f4*e2
     z_r=z_r
     z_l=z_l+L
     a0=(f1*z_l**2*(z_l-3.d0*z_r)+z_r*(f3*z_l**2*(-z_l+z_r) &
          +z_r*(f2*(3.d0*z_l-z_r)+f4*z_l*(-z_l+z_r))))/(z_l-z_r)**3
     a1=(f3*z_l**3+z_l*(6.d0*f1-6.d0*f2+(f3+2.d0*f4)*z_l)*z_r &
          -(2*f3+f4)*z_l*z_r**2-f4*z_r**3)/(z_l-z_r)**3
     a2=(-3*f1*(z_l+z_r)+3.d0*f2*(z_l+z_r)-(z_l-z_r)*(2*f3*z_l &
          +f4*z_l+f3*z_r+2*f4*z_r))/(z_l-z_r)**3
     a3=(2.d0*f1-2.d0*f2+(f3+f4)*(z_l-z_r))/(z_l-z_r)**3
     do iz=nz_r,nz_l
        z=dble(iz-1)/dble(dfftp%nr3)*L
        vg_r(iz)=(a0+a1*z+a2*z**2+a3*z**3)
     enddo
     
     call env_cft_1z(vg_r,1,dfftp%nr3,dfftp%nr3,-1,vg)
     
     vg3(:,ng_2d)=vg(:)*e2 ! factor e2: hartree -> Ry.

     deallocate (vg,vg_r,vr,vr_r)
  endif ! if( ng_2d > 0 )

! Hartree Energy
  ehart=0.d0
  eh = 0d0
  do ng_2d = 1, ngm_2d
     k1 = mill_2d(1,ng_2d)
     k2 = mill_2d(2,ng_2d)
     eh = eh + sum( vg3(:,ng_2d)*conjg(rhog3(:,ng_2d)) )
  enddo
  ehart=ehart+eh
  if( gamma_only ) then
     ehart = ehart * 2d0
     ng_2d = imill_2d(0,0)
     if( ng_2d > 0 ) then
        ehart = ehart - sum( vg3(:,ng_2d)*conjg(rhog3(:,ng_2d)) )
     endif
  endif
  ehart = ehart*omega*0.5d0
  !
  call env_mp_sum( ehart, intra_bgrp_comm )
  !
! Map to FFT mesh (dfftp%nrx)
  aux=0.0d0
  do ng=1,ngm
     n1 = mill(1,ng)
     n2 = mill(2,ng)
     ng_2d = imill_2d(n1,n2)
     n3 = mill(3,ng) + 1
     if (n3<1) n3 = n3 + dfftp%nr3
     aux(dfftp%nl(ng))= aux(dfftp%nl(ng)) + vg3(n3,ng_2d)
     if (gamma_only) then
        aux(dfftp%nlm(ng))=CONJG(aux(dfftp%nl(ng)))
     endif
  enddo

  deallocate (rhog3,vg3)

  RETURN
END SUBROUTINE env_esm_hartree_bc4

END MODULE env_esm