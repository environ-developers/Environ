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
MODULE core_fd
!----------------------------------------------------------------------------
  !
  USE environ_types
  USE core_types
  !
  PRIVATE
  !
  PUBLIC :: gradient_fd
  !
CONTAINS
!----------------------------------------------------------------------
  SUBROUTINE gradient_fd( fd, f, grad )
!----------------------------------------------------------------------
    !
    USE scatter_mod,       ONLY : scatter_grid
    USE mp,                ONLY : mp_sum
    !
    IMPLICIT NONE
    !
    TYPE( fd_core ), INTENT(IN) :: fd
    TYPE( environ_density ), INTENT(IN) :: f
    TYPE( environ_gradient ), INTENT(INOUT) :: grad
    !
    INTEGER :: idx, idx0, j0, k0, i, ir, ir_end, ipol, in
    INTEGER :: ixc, iyc, izc, ixp, ixm, iyp, iym, izp, izm
    INTEGER, DIMENSION(:), ALLOCATABLE :: ix, iy, iz
    REAL( DP ), DIMENSION( :, : ), ALLOCATABLE :: gradtmp, gradaux
    !
    dfft => fd % dfft
    nnr => cell%nnr
    nfdpoint => fd % nfdpoint
    !
    ALLOCATE( ix(-nfdpoint:nfdpoint), iy(-nfdpoint:nfdpoint), iz(-nfdpoint:nfdpoint) )
    !
    ALLOCATE( gradtmp( dfft%nr1x*dfft%nr2x*dfft%nr3x, 3 ) )
    gradtmp = 0.D0
    !
! BACKWARD COMPATIBILITY
! Compatible with QE-5.X QE-6.1.X
!    idx0 = dfftp%nr1x*dfftp%nr2x*dfftp%ipp(me_bgrp+1)
!    ir_end = dfftp%nr1x*dfftp%nr2x*dfftp%npl
! Compatible with QE-6.2, QE-6.2.1 and QE-GIT
#if defined (__MPI)
    j0 = dfft%my_i0r2p ; k0 = dfft%my_i0r3p
    ir_end = MIN(nnr,dfft%nr1x*dfft%my_nr2p*dfft%my_nr3p)
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
       iz(0) = idx / (dfftp%nr1x*dfft%my_nr2p)
       idx   = idx - (dfft%nr1x*dfft%my_nr2p)*iz(0)
       iz(0) = iz(0) + k0
       iy(0) = idx / dfft%nr1x
       idx   = idx - dfft%nr1x*iy(0)
       iy(0) = iy(0) + j0
       ix(0) = idx
! END BACKWARD COMPATIBILITY
       !
       ! ... do not include points outside the physical range
       !
       IF ( ix(0) >= dfft%nr1 .OR. iy(0) >= dfft%nr2 .OR. iz(0) >= dfft%nr3 ) CYCLE
       !
       DO in = 1, fd%nfdpoint
          ix(in) = ix(in-1) + 1
          IF( ix(in) .GT. dfft%nr1x-1 ) ix(in) = 0
          ix(-in) = ix(-in+1) - 1
          IF( ix(-in) .LT. 0 ) ix(-in) = dfft%nr1x-1
          iy(in) = iy(in-1) + 1
          IF( iy(in) .GT. dfft%nr2x-1 ) iy(in) = 0
          iy(-in) = iy(-in+1) - 1
          IF( iy(-in) .LT. 0 ) iy(-in) = dfft%nr2x-1
          iz(in) = iz(in-1) + 1
          IF( iz(in) .GT. dfft%nr3x-1 ) iz(in) = 0
          iz(-in) = iz(-in+1) - 1
          IF( iz(-in) .LT. 0 ) iz(-in) = dfft%nr3x-1
       ENDDO
       !
       DO in = -fd%nfdpoint, fd%nfdpoint
          i = ix(in) + iy(0) * dfft%nr1x + iz(0) * dfft%nr1x * dfft%nr2x + 1
          gradtmp(i,1) = gradtmp(i,1) - fd%icfd(in)*f%of_r(ir)*dfft%nr1
          i = ix(0) + iy(in) * dfft%nr1x + iz(0) * dfft%nr1x * dfft%nr2x + 1
          gradtmp(i,2) = gradtmp(i,2) - fd%icfd(in)*f%of_r(ir)*dfft%nr2
          i = ix(0) + iy(0) * dfft%nr1x + iz(in) * dfft%nr1x * dfft%nr2x + 1
          gradtmp(i,3) = gradtmp(i,3) - fd%icfd(in)*f%of_r(ir)*dfft%nr3
       ENDDO
       !
    ENDDO
    !
    DEALLOCATE( ix, iy, iz )
    !
    ALLOCATE( gradaux(nnr,3) )
#if defined (__MPI)
    DO ipol = 1, 3
       CALL mp_sum( gradtmp(:,ipol), cell%comm )
       CALL scatter_grid ( dfft, gradtmp(:,ipol), gradaux(:,ipol) )
    ENDDO
#else
    gradaux(1:nnr,:) = gradtmp(1:nnr,:)
#endif
    !
    DEALLOCATE( gradtmp )
    !
    DO ir = 1,nnr
       grad%of_r(:,ir) = MATMUL( cell%bg, gradaux(ir,:) )
    ENDDO
    DEALLOCATE( gradaux )
    !
    grad%of_r = grad%of_r / DBLE(fd%ncfd) / cell%alat
    !
    RETURN
    !
!----------------------------------------------------------------------
  END SUBROUTINE gradient_fd
!----------------------------------------------------------------------
!----------------------------------------------------------------------------
END MODULE core_fd
!----------------------------------------------------------------------------
