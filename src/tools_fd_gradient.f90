! Copyright (C) 2018 ENVIRON (www.quantum-environment.org)
! Copyright (C) 2006-2010 Quantum ESPRESSO group
!
!    This file is part of Environ version 1.1
!
!    Environ 1.1 is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 2 of the License, or
!    (at your option) any later version.
!
!    Environ 1.1 is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more detail, either the file
!    `License' in the root directory of the present distribution, or
!    online at <http://www.gnu.org/licenses/>.
!
!> Module to compute finite-differences gradients on dense real space grid
!
! Authors: Oliviero Andreussi (Department of Physics, UNT)
!
!
!----------------------------------------------------------------------------
MODULE tools_fd_gradient
!----------------------------------------------------------------------------
  !
  USE env_kinds, ONLY: DP
  !
  IMPLICIT NONE
  !
CONTAINS
!----------------------------------------------------------------------
  SUBROUTINE calc_fd_gradient( nfdpoint, icfd, ncfd, nnr, f, grad )
!----------------------------------------------------------------------
    USE env_kinds,         ONLY : DP
    USE env_cell_base,     ONLY : at, bg, alat
    USE env_fft_base,      ONLY : dfftp
    USE env_scatter_mod,   ONLY : env_scatter_grid
    USE env_mp,            ONLY : env_mp_sum
    USE env_mp_bands,      ONLY : me_bgrp, intra_bgrp_comm
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN)  :: nfdpoint
    INTEGER, INTENT(IN)  :: ncfd
    INTEGER, INTENT(IN)  :: icfd(-nfdpoint:nfdpoint)
    INTEGER, INTENT(IN)  :: nnr
    REAL( DP ), DIMENSION( nnr ), INTENT(IN) :: f
    REAL( DP ), DIMENSION( 3, nnr ), INTENT(OUT) :: grad
    !
    INTEGER :: idx, idx0, j0, k0, i, ir, ir_end, ipol, in
    INTEGER :: ix(-nfdpoint:nfdpoint),iy(-nfdpoint:nfdpoint),iz(-nfdpoint:nfdpoint)
    INTEGER :: ixc, iyc, izc, ixp, ixm, iyp, iym, izp, izm
    REAL( DP ), DIMENSION( :, : ), ALLOCATABLE :: gradtmp, gradaux
    !
    grad = 0.D0
    !
    ALLOCATE( gradtmp( dfftp%nr1x*dfftp%nr2x*dfftp%nr3x, 3 ) )
    gradtmp = 0.D0
    !
! BACKWARD COMPATIBILITY
! Compatible with QE-5.X QE-6.1.X
!    idx0 = dfftp%nr1x*dfftp%nr2x*dfftp%ipp(me_bgrp+1)
!    ir_end = dfftp%nr1x*dfftp%nr2x*dfftp%npl
! Compatible with QE-6.2, QE-6.2.1 and QE-GIT
#if defined (__MPI)
    j0 = dfftp%my_i0r2p ; k0 = dfftp%my_i0r3p
    ir_end = MIN(nnr,dfftp%nr1x*dfftp%my_nr2p*dfftp%my_nr3p)
#else
    j0 = 0 ; k0 = 0
    ir_end = nnr
#endif
! END BACKWARD COMPATIBILITY
    !
    DO ir = 1, ir_end
       !
! BACKWARD COMPATIBILITY
! Compatible with QE-5.X QE-6.1.X
!       i = idx0 + ir - 1
!       iz(0) = i / (dfftp%nr1x*dfftp%nr2x)
!       i     = i - (dfftp%nr1x*dfftp%nr2x)*iz(0)
!       iy(0) = i / dfftp%nr1x
!       ix(0) = i - dfftp%nr1x*iy(0)
! Compatible with QE-6.2, QE-6.2.1 and QE-GIT
       idx   = ir - 1
       iz(0) = idx / (dfftp%nr1x*dfftp%my_nr2p)
       idx   = idx - (dfftp%nr1x*dfftp%my_nr2p)*iz(0)
       iz(0) = iz(0) + k0
       iy(0) = idx / dfftp%nr1x
       idx   = idx - dfftp%nr1x*iy(0)
       iy(0) = iy(0) + j0
       ix(0) = idx
! END BACKWARD COMPATIBILITY
       !
       ! ... do not include points outside the physical range
       !
       IF ( ix(0) >= dfftp%nr1 .OR. iy(0) >= dfftp%nr2 .OR. iz(0) >= dfftp%nr3 ) CYCLE
       !
       DO in = 1, nfdpoint
          ix(in) = ix(in-1) + 1
          IF( ix(in) .GT. dfftp%nr1x-1 ) ix(in) = 0
          ix(-in) = ix(-in+1) - 1
          IF( ix(-in) .LT. 0 ) ix(-in) = dfftp%nr1x-1
          iy(in) = iy(in-1) + 1
          IF( iy(in) .GT. dfftp%nr2x-1 ) iy(in) = 0
          iy(-in) = iy(-in+1) - 1
          IF( iy(-in) .LT. 0 ) iy(-in) = dfftp%nr2x-1
          iz(in) = iz(in-1) + 1
          IF( iz(in) .GT. dfftp%nr3x-1 ) iz(in) = 0
          iz(-in) = iz(-in+1) - 1
          IF( iz(-in) .LT. 0 ) iz(-in) = dfftp%nr3x-1
       ENDDO
       !
       DO in = -nfdpoint, nfdpoint
          i = ix(in) + iy(0) * dfftp%nr1x + iz(0) * dfftp%nr1x * dfftp%nr2x + 1
          gradtmp(i,1) = gradtmp(i,1) - icfd(in)*f(ir)*dfftp%nr1
          i = ix(0) + iy(in) * dfftp%nr1x + iz(0) * dfftp%nr1x * dfftp%nr2x + 1
          gradtmp(i,2) = gradtmp(i,2) - icfd(in)*f(ir)*dfftp%nr2
          i = ix(0) + iy(0) * dfftp%nr1x + iz(in) * dfftp%nr1x * dfftp%nr2x + 1
          gradtmp(i,3) = gradtmp(i,3) - icfd(in)*f(ir)*dfftp%nr3
       ENDDO
       !
    ENDDO
    !
    ALLOCATE( gradaux(nnr,3) )
#if defined (__MPI)
    DO ipol = 1, 3
       CALL env_mp_sum( gradtmp(:,ipol), intra_bgrp_comm )
       CALL env_scatter_grid ( dfftp, gradtmp(:,ipol), gradaux(:,ipol) )
    ENDDO
#else
    gradaux(1:nnr,:) = gradtmp(1:nnr,:)
#endif
    !
    DEALLOCATE( gradtmp )
    !
    DO ir = 1,nnr
       grad(:,ir) = MATMUL( bg, gradaux(ir,:) )
    ENDDO
    DEALLOCATE( gradaux )
    !
    grad = grad / DBLE(ncfd) / alat
    !
    RETURN
    !
!----------------------------------------------------------------------
  END SUBROUTINE calc_fd_gradient
!----------------------------------------------------------------------
!----------------------------------------------------------------------
  SUBROUTINE init_fd_gradient( ifdtype, nfdpoint, ncfd, icfd )
!----------------------------------------------------------------------
    !
    USE env_kinds,         ONLY : DP
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN)  :: ifdtype, nfdpoint
    INTEGER, INTENT(OUT) :: ncfd
    INTEGER, INTENT(OUT) :: icfd(-nfdpoint:nfdpoint)
    !
    INTEGER :: in
    !
    ncfd = 0
    icfd = 0
    !
    SELECT CASE ( ifdtype )
       !
    CASE ( 1 )
       ! (2N+1)-point Central Differences
       IF ( nfdpoint .EQ. 1 ) THEN
          ncfd = 2
          icfd(  1 ) =   1
       ELSE IF ( nfdpoint .EQ. 2 ) THEN
          ncfd = 12
          icfd(  2 ) =  -1
          icfd(  1 ) =   8
       ELSE IF ( nfdpoint .EQ. 3 ) THEN
          ncfd = 60
          icfd(  3 ) =   1
          icfd(  2 ) =  -9
          icfd(  1 ) =  45
       ELSE IF ( nfdpoint .EQ. 4 ) THEN
          ncfd = 840
          icfd(  4 ) =  -3
          icfd(  3 ) =  32
          icfd(  2 ) =-168
          icfd(  1 ) = 672
       ELSE
          WRITE(*,*)'ERROR: wrong number of points',nfdpoint,&
               &' for finite difference type ',ifdtype
          STOP
       ENDIF
       !
    CASE ( 2 )
       ! Low-Noise Lanczos Differentiators ( M = 2 )
       IF ( nfdpoint .GE. 2 ) THEN
          ncfd = (nfdpoint)*(nfdpoint+1)*(2*nfdpoint+1)/3
          DO in = 1,nfdpoint
             icfd( in ) = in
          ENDDO
       ELSE
          WRITE(*,*)'ERROR: wrong number of points',nfdpoint,&
               &' for finite difference type ',ifdtype
          STOP
       END IF
       !
    CASE ( 3 )
       ! Super Lanczos Low-Noise Differentiators ( M = 4 )
       IF ( nfdpoint .EQ. 3 ) THEN
          ncfd = 252
          icfd(  3 ) = -22
          icfd(  2 ) =  67
          icfd(  1 ) =  58
       ELSE IF ( nfdpoint .EQ. 4 ) THEN
          ncfd = 1188
          icfd(  4 ) = -86
          icfd(  3 ) = 142
          icfd(  2 ) = 193
          icfd(  1 ) = 126
       ELSE IF ( nfdpoint .EQ. 5 ) THEN
          ncfd = 5148
          icfd(  5 ) =-300
          icfd(  4 ) = 294
          icfd(  3 ) = 532
          icfd(  2 ) = 503
          icfd(  1 ) = 296
       ELSE
          WRITE(*,*)'ERROR: wrong number of points',nfdpoint,&
               &' for finite difference type ',ifdtype
          STOP
       ENDIF
       !
    CASE ( 4 )
       ! Smooth Noise-Robust Differentiators  ( n = 2 )
       IF ( nfdpoint .EQ. 2 ) THEN
          ncfd = 8
          icfd(  2 ) =   1
          icfd(  1 ) =   2
       ELSE IF ( nfdpoint .EQ. 3 ) THEN
          ncfd = 32
          icfd(  3 ) =   1
          icfd(  2 ) =   4
          icfd(  1 ) =   5
       ELSE IF ( nfdpoint .EQ. 4 ) THEN
          ncfd = 128
          icfd(  4 ) =   1
          icfd(  3 ) =   6
          icfd(  2 ) =  14
          icfd(  1 ) =  14
       ELSE IF ( nfdpoint .EQ. 5 ) THEN
          ncfd = 512
          icfd(  5 ) =   1
          icfd(  4 ) =   8
          icfd(  3 ) =  27
          icfd(  2 ) =  48
          icfd(  1 ) =  42
       ELSE
          WRITE(*,*)'ERROR: wrong number of points',nfdpoint,&
               &' for finite difference type ',ifdtype
          STOP
       ENDIF
       !
    CASE ( 5 )
       ! Smooth Noise-Robust Differentiators  ( n = 4 )
       IF ( nfdpoint .EQ. 3 ) THEN
          ncfd = 96
          icfd(  3 ) =  -5
          icfd(  2 ) =  12
          icfd(  1 ) =  39
       ELSE IF ( nfdpoint .EQ. 4 ) THEN
          ncfd = 96
          icfd(  4 ) =  -2
          icfd(  3 ) =  -1
          icfd(  2 ) =  16
          icfd(  1 ) =  27
       ELSE IF ( nfdpoint .EQ. 5 ) THEN
          ncfd = 1536
          icfd(  5 ) = -11
          icfd(  4 ) = -32
          icfd(  3 ) =  39
          icfd(  2 ) = 256
          icfd(  1 ) = 322
       ELSE
          WRITE(*,*)'ERROR: wrong number of points',nfdpoint,&
               &' for finite difference type ',ifdtype
          STOP
       ENDIF
       !
    CASE DEFAULT
       !
       WRITE(*,*)'ERROR: finite difference type unknown, ifdtype=',ifdtype
       STOP
       !
    END SELECT
    !
    DO in = 1,nfdpoint
       icfd( -in ) = - icfd( in )
    ENDDO
    !
    RETURN
    !
!----------------------------------------------------------------------
  END SUBROUTINE init_fd_gradient
!----------------------------------------------------------------------
!----------------------------------------------------------------------------
END MODULE tools_fd_gradient
!----------------------------------------------------------------------------
