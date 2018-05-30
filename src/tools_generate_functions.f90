! Copyright (C) 2018 ENVIRON (www.quantum-environment.org)
! Copyright (C) 2006-2010 Quantum ESPRESSO group
!
!    This file is part of Environ version 1.0
!
!    Environ 1.0 is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 2 of the License, or
!    (at your option) any later version.
!
!    Environ 1.0 is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more detail, either the file
!    `License' in the root directory of the present distribution, or
!    online at <http://www.gnu.org/licenses/>.
!
! Module to generate functions on the real space dense grid
!
! Authors: Oliviero Andreussi (Department of Physics, UNT)
!
!----------------------------------------------------------------------------
MODULE tools_generate_functions
!----------------------------------------------------------------------------
  !
  USE kinds, ONLY: DP
  !
  IMPLICIT NONE
  !
CONTAINS
!--------------------------------------------------------------------
  SUBROUTINE compute_convolution_fft( nnr, fa, fb, fc )
!--------------------------------------------------------------------
    !
    ! ... Calculates using reciprocal-space techniques the
    ! ... convolution of function fa with function fb and put
    ! ... the result in function fc
    !
    USE cell_base,      ONLY : omega
    USE fft_base,       ONLY : dfftp
    USE fft_interfaces, ONLY : fwfft, invfft
    USE control_flags,  ONLY : gamma_only
    USE gvect,          ONLY : nl, nlm, ngm, gg, gstart
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN)         :: nnr
    REAL( DP ), INTENT(IN)      :: fa( nnr ), fb( nnr )
    REAL( DP ), INTENT(OUT)     :: fc( nnr )
    !
    COMPLEX( DP ), DIMENSION( : ), ALLOCATABLE :: auxr, auxg
    !
    ALLOCATE( auxr( nnr ) )
    auxr(:) = CMPLX( fa(:), 0.D0, kind=DP )
    CALL fwfft('Dense', auxr, dfftp)
    !
    ALLOCATE( auxg( nnr ) ) 
    auxg = 0.D0
    auxg(nl(1:ngm)) = auxr(nl(1:ngm))
    !
    auxr(:) = CMPLX( fb(:), 0.D0, kind=DP )
    CALL fwfft('Dense', auxr, dfftp)
    !
    auxg(nl(1:ngm)) = auxg(nl(1:ngm)) * auxr(nl(1:ngm))
    !
    DEALLOCATE( auxr )
    !
    IF ( gamma_only ) auxg(nlm(1:ngm)) = CMPLX( REAL( auxg(nl(1:ngm)) ), -AIMAG( auxg(nl(1:ngm)) ) ,kind=DP)
    !
    CALL invfft('Dense',auxg, dfftp)
    !
    fc(:) = REAL( auxg(:) ) * omega
    !
    DEALLOCATE( auxg )
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE compute_convolution_fft
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE planar_average( nnr, naxis, axis, shift, reverse, f, f1d )
!--------------------------------------------------------------------
    !
    USE kinds,            ONLY : DP
    USE fft_base,         ONLY : dfftp
    USE mp,               ONLY : mp_sum
    USE mp_bands,         ONLY : me_bgrp, intra_bgrp_comm
    !
    IMPLICIT NONE
    !
    ! ... Declares variables
    !
    INTEGER, INTENT(IN)       :: nnr, naxis, axis, shift
    LOGICAL, INTENT(IN)       :: reverse
    REAL( DP ), INTENT(INOUT) :: f( nnr )
    REAL( DP ), INTENT(INOUT) :: f1d( naxis )
    !
    ! ... Local variables
    !
    INTEGER                   :: i, j, k, ir, ir_end
    INTEGER                   :: idx, idx0, narea
    !
    INTEGER                   :: j0, k0
    !
    narea = dfftp%nr1*dfftp%nr2*dfftp%nr3 / naxis
    !
    IF ( reverse ) THEN
       f = 0.D0
    ELSE
       f1d = 0.D0
    END IF
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
       ! ... three dimensional indexes
       !
! BACKWARD COMPATIBILITY
! Compatible with QE-5.X QE-6.1.X
!       idx = idx0 + ir - 1
!       k   = idx / (dfftp%nr1x*dfftp%nr2x)
!       idx = idx - (dfftp%nr1x*dfftp%nr2x)*k
!       j   = idx / dfftp%nr1x
!       idx = idx - dfftp%nr1x*j
!       i   = idx
! Compatible with QE-6.2, QE-6.2.1 and QE-GIT
       idx = ir - 1
       k   = idx / (dfftp%nr1x*dfftp%my_nr2p)
       idx = idx - (dfftp%nr1x*dfftp%my_nr2p)*k
       k   = k + k0
       j   = idx / dfftp%nr1x
       idx = idx - dfftp%nr1x * j
       j   = j + j0
       i   = idx
! END BACKWARD COMPATIBILITY
       !
       ! ... do not include points outside the physical range
       !
       IF ( i >= dfftp%nr1 .OR. j >= dfftp%nr2 .OR. k >= dfftp%nr3 ) CYCLE
       !
       SELECT CASE ( axis )
       CASE ( 1 )
          idx = i
       CASE ( 2 )
          idx = j
       CASE ( 3 )
          idx = k
       END SELECT
       !
       idx = idx + 1 + shift
       !
       IF ( idx .GT. naxis ) THEN
          idx = idx - naxis
       ELSE IF (idx .LE. 0 ) THEN
          idx = idx + naxis
       ENDIF
       !
       IF ( reverse ) THEN
          f(ir) = f1d(idx)
       ELSE
          f1d(idx) = f1d(idx) + f(ir)
       END IF
       !
    END DO
    !
    IF ( .NOT. reverse ) THEN
       CALL mp_sum( f1d(:), intra_bgrp_comm )
       f1d = f1d / DBLE(narea)
    END IF
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE planar_average
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE generate_gaussian( nnr, dim, axis, charge, spread, pos, rho )
!--------------------------------------------------------------------
    !
    USE kinds,            ONLY : DP
    USE constants,        ONLY : sqrtpi
    USE io_global,        ONLY : stdout
    USE cell_base,        ONLY : at, bg, alat, omega
    USE fft_base,         ONLY : dfftp
    USE mp,               ONLY : mp_sum
    USE mp_bands,         ONLY : me_bgrp, intra_bgrp_comm
    !
    IMPLICIT NONE
    !
    REAL(DP), PARAMETER       :: tol = 1.D-6
    REAL(DP), PARAMETER       :: exp_tol = 3.6D1
    !
    ! ... Declares variables
    !
    INTEGER, INTENT(IN)       :: nnr, dim, axis
    REAL( DP ), INTENT(IN)    :: charge, spread
    REAL( DP ), INTENT(IN)    :: pos( 3 )
    REAL( DP ), INTENT(INOUT) :: rho( nnr )
    !
    ! ... Local variables
    !
    INTEGER                   :: i, j, j0, k, k0, ir, ip, ir_end
    INTEGER                   :: idx, idx0
    !
    REAL( DP )                :: inv_nr1, inv_nr2, inv_nr3
    REAL( DP )                :: scale, spr2, dist, length
    REAL( DP )                :: r( 3 ), s( 3 )
    REAL( DP ), ALLOCATABLE   :: rholocal ( : )
    !
    IF ( dfftp%nr1 .EQ. 0 .OR. dfftp%nr2 .EQ. 0 .OR. dfftp%nr3 .EQ. 0 ) THEN
       WRITE(stdout,*)'ERROR: wrong grid dimension',dfftp%nr1,dfftp%nr2,dfftp%nr3
       STOP
    ENDIF
    inv_nr1 = 1.D0 / DBLE( dfftp%nr1 )
    inv_nr2 = 1.D0 / DBLE( dfftp%nr2 )
    inv_nr3 = 1.D0 / DBLE( dfftp%nr3 )
    !
    IF ( ABS( spread ) .LT. tol ) THEN
       WRITE(stdout,*)'ERROR: wrong spread for Gaussian function',spread
       STOP
    ENDIF
    !
    IF (axis.LT.1.OR.axis.GT.3) &
         WRITE(stdout,*)'WARNING: wrong axis in generate_gaussian'
    IF ( dim .EQ. 0 ) THEN
       scale = charge / ( sqrtpi * spread )**3
    ELSE IF ( dim .EQ. 1 ) THEN
       length = ABS( at(axis,axis) * alat )
       IF ( length .LT. tol ) THEN
          WRITE(stdout,*)'ERROR: unphysically small dimension of cell',length
          STOP
       ENDIF
       scale = charge / length / ( sqrtpi * spread )**2
    ELSE IF ( dim .EQ. 2 ) THEN
       length = ABS( at(axis,axis) * alat )
       IF ( length .LT. tol .OR. omega .LT. tol ) THEN
          WRITE(stdout,*)'ERROR: unphysically small dimensions of cell',length, omega
          STOP
       ENDIF
       scale = charge * length / omega / ( sqrtpi * spread )
    ELSE
       WRITE(stdout,*)'WARNING: wrong dim in generate_gaussian'
    ENDIF
    IF ( ABS( alat ) .LT. tol ) THEN
       WRITE(stdout,*)'ERROR: unphysically small alat',alat
       STOP
    ENDIF
    spr2 = ( spread / alat )**2
    ALLOCATE( rholocal( nnr ) )
    rholocal = 0.D0
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
       ! ... three dimensional indexes
       !
! BACKWARD COMPATIBILITY
! Compatible with QE-5.X QE-6.1.X
!       idx = idx0 + ir - 1
!       k   = idx / (dfftp%nr1x*dfftp%nr2x)
!       idx = idx - (dfftp%nr1x*dfftp%nr2x)*k
!       j   = idx / dfftp%nr1x
!       idx = idx - dfftp%nr1x*j
!       i   = idx
! Compatible with QE-6.2, QE-6.2.1 and QE-GIT
       idx = ir - 1
       k   = idx / (dfftp%nr1x*dfftp%my_nr2p)
       idx = idx - (dfftp%nr1x*dfftp%my_nr2p)*k
       k   = k + k0
       j   = idx / dfftp%nr1x
       idx = idx - dfftp%nr1x * j
       j   = j + j0
       i   = idx
! END BACKWARD COMPATIBILITY
       !
       ! ... do not include points outside the physical range
       !
       IF ( i >= dfftp%nr1 .OR. j >= dfftp%nr2 .OR. k >= dfftp%nr3 ) CYCLE
       !
       DO ip = 1, 3
          r(ip) = DBLE( i )*inv_nr1*at(ip,1) + &
               DBLE( j )*inv_nr2*at(ip,2) + &
               DBLE( k )*inv_nr3*at(ip,3)
       END DO
       !
       r(:) = pos(:) - r(:)
       !
       !  ... possibly 2D or 1D gaussians
       !
       IF ( dim .EQ. 1) THEN
          r(axis) = 0.D0
       ELSE IF ( dim .EQ. 2 ) THEN
          DO i = 1, 3
             IF ( i .NE. axis ) r(i) = 0.D0
          ENDDO
       END IF
       !
       ! ... minimum image convention
       !
       s(:) = MATMUL( r(:), bg(:,:) )
       s(:) = s(:) - ANINT(s(:))
       r(:) = MATMUL( at(:,:), s(:) )
       !
       dist = SUM( r * r ) / spr2
       !
       IF ( dist .GT. exp_tol ) THEN
          rholocal( ir ) = 0.D0
       ELSE
          rholocal( ir ) = scale * EXP(-dist)
       ENDIF
       !
    END DO
    !
    rho = rho + rholocal
    DEALLOCATE( rholocal )
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE generate_gaussian
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE generate_gradgaussian( nnr, dim, axis, charge, spread, pos, gradrho )
!--------------------------------------------------------------------
    !
    USE kinds,            ONLY : DP
    USE constants,        ONLY : sqrtpi
    USE io_global,        ONLY : stdout
    USE cell_base,        ONLY : at, bg, alat, omega
    USE fft_base,         ONLY : dfftp
    USE mp_bands,         ONLY : me_bgrp, intra_bgrp_comm
    !
    IMPLICIT NONE
    !
    REAL(DP), PARAMETER       :: tol = 1.D-6
    REAL(DP), PARAMETER       :: exp_tol = 3.6D1
    !
    ! ... Declares variables
    !
    INTEGER, INTENT(IN)       :: nnr, dim, axis
    REAL( DP ), INTENT(IN)    :: charge, spread
    REAL( DP ), INTENT(IN)    :: pos( 3 )
    REAL( DP ), INTENT(INOUT) :: gradrho( 3, nnr )
    !
    ! ... Local variables
    !
    INTEGER                   :: i, j, j0, k, k0, ir, ip, ir_end
    INTEGER                   :: idx, idx0
    !
    REAL( DP )                :: inv_nr1, inv_nr2, inv_nr3
    REAL( DP )                :: scale, spr2, dist, length
    REAL( DP )                :: r( 3 ), s( 3 )
    REAL( DP ), ALLOCATABLE   :: gradrholocal ( :, : )
    !
    IF ( dfftp%nr1 .EQ. 0 .OR. dfftp%nr2 .EQ. 0 .OR. dfftp%nr3 .EQ. 0 ) THEN
       WRITE(stdout,*)'ERROR: wrong grid dimension',dfftp%nr1,dfftp%nr2,dfftp%nr3
       STOP
    ENDIF
    inv_nr1 = 1.D0 / DBLE( dfftp%nr1 )
    inv_nr2 = 1.D0 / DBLE( dfftp%nr2 )
    inv_nr3 = 1.D0 / DBLE( dfftp%nr3 )
    !
    IF ( ABS( spread ) .LT. tol ) THEN
       WRITE(stdout,*)'ERROR: wrong spread for Gaussian function',spread
       STOP
    ENDIF
    !
    IF (axis.LT.1.OR.axis.GT.3) &
         WRITE(stdout,*)'WARNING: wrong axis in generate_gradgaussian'
    IF ( dim .EQ. 0 ) THEN
       scale = charge / ( sqrtpi * spread )**3
    ELSE IF ( dim .EQ. 1 ) THEN
       length = ABS( at(axis,axis) * alat )
       IF ( length .LT. tol ) THEN
          WRITE(stdout,*)'ERROR: unphysically small dimension of cell',length
          STOP
       ENDIF
       scale = charge / length / ( sqrtpi * spread )**2
    ELSE IF ( dim .EQ. 2 ) THEN
       length = ABS( at(axis,axis) * alat )
       IF ( length .LT. tol .OR. omega .LT. tol ) THEN
          WRITE(stdout,*)'ERROR: unphysically small dimensions of cell',length, omega
          STOP
       ENDIF
       scale = charge * length / omega / ( sqrtpi * spread )
    ELSE
       WRITE(stdout,*)'WARNING: wrong dim in generate_gradgaussian'
    ENDIF
    IF ( ABS( alat ) .LT. tol ) THEN
       WRITE(stdout,*)'ERROR: unphysically small alat',alat
       STOP
    ENDIF
    scale = scale * 2.D0 / spread**2 * alat
    spr2 = ( spread / alat )**2
    ALLOCATE( gradrholocal( 3, nnr ) )
    gradrholocal = 0.D0
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
       ! ... three dimensional indexes
       !
! BACKWARD COMPATIBILITY
! Compatible with QE-5.X QE-6.1.X
!       idx = idx0 + ir - 1
!       k   = idx / (dfftp%nr1x*dfftp%nr2x)
!       idx = idx - (dfftp%nr1x*dfftp%nr2x)*k
!       j   = idx / dfftp%nr1x
!       idx = idx - dfftp%nr1x*j
!       i   = idx
! Compatible with QE-6.2, QE-6.2.1 and QE-GIT
       idx = ir - 1
       k   = idx / (dfftp%nr1x*dfftp%my_nr2p)
       idx = idx - (dfftp%nr1x*dfftp%my_nr2p)*k
       k   = k + k0
       j   = idx / dfftp%nr1x
       idx = idx - dfftp%nr1x * j
       j   = j + j0
       i   = idx
! END BACKWARD COMPATIBILITY
       !
       ! ... do not include points outside the physical range
       !
       IF ( i >= dfftp%nr1 .OR. j >= dfftp%nr2 .OR. k >= dfftp%nr3 ) CYCLE
       !
       DO ip = 1, 3
          r(ip) = DBLE( i )*inv_nr1*at(ip,1) + &
               DBLE( j )*inv_nr2*at(ip,2) + &
               DBLE( k )*inv_nr3*at(ip,3)
       END DO
       !
       r(:) = pos(:) - r(:)
       !
       !  ... possibly 2D or 1D gaussians
       !
       IF ( dim .EQ. 1) THEN
          r(axis) = 0.D0
       ELSE IF ( dim .EQ. 2 ) THEN
          DO i = 1, 3
             IF ( i .NE. axis ) r(i) = 0.D0
          ENDDO
       END IF
       !
       ! ... minimum image convention
       !
       s(:) = MATMUL( r(:), bg(:,:) )
       s(:) = s(:) - ANINT(s(:))
       r(:) = MATMUL( at(:,:), s(:) )
       !
       dist = SUM( r * r ) / spr2
       !
       IF ( dist .GT. exp_tol ) THEN
          gradrholocal( :, ir ) = 0.D0
       ELSE
          gradrholocal( :, ir ) = EXP(-dist) * r(:)
       ENDIF
       !
    END DO
    !
    gradrho = gradrho + gradrholocal * scale
    DEALLOCATE( gradrholocal )
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE generate_gradgaussian
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE generate_exponential(nnr, dim, axis, width, spread, pos, rho )
!--------------------------------------------------------------------
    !
    USE kinds,            ONLY : DP
    USE io_global,        ONLY : stdout
    USE cell_base,        ONLY : at, bg, alat
    USE fft_base,         ONLY : dfftp
    USE mp_bands,         ONLY : me_bgrp, intra_bgrp_comm
    !
    IMPLICIT NONE
    !
    ! ... Declares variables
    !
    INTEGER, INTENT(IN)       :: nnr, dim, axis
    REAL( DP ), INTENT(IN)    :: width, spread
    REAL( DP ), INTENT(IN)    :: pos( 3 )
    REAL( DP ), INTENT(INOUT) :: rho( nnr )
    !
    ! ... Local variables
    !
    INTEGER                   :: i, j, j0, k, k0, ir, ir_end, ip
    INTEGER                   :: idx, idx0
    !
    REAL( DP )                :: inv_nr1, inv_nr2, inv_nr3
    REAL( DP )                :: dist, arg
    REAL( DP )                :: r( 3 ), s( 3 )
    REAL( DP ), ALLOCATABLE   :: rholocal ( : )
    REAL( DP ), PARAMETER     :: exp_arg_limit = 100.D0
    !
    IF ( dfftp%nr1 .EQ. 0 .OR. dfftp%nr2 .EQ. 0 .OR. dfftp%nr3 .EQ. 0 ) THEN
       WRITE(stdout,*)'ERROR: wrong grid dimension',dfftp%nr1,dfftp%nr2,dfftp%nr3
       STOP
    ENDIF
    inv_nr1 = 1.D0 / DBLE( dfftp%nr1 )
    inv_nr2 = 1.D0 / DBLE( dfftp%nr2 )
    inv_nr3 = 1.D0 / DBLE( dfftp%nr3 )
    !
    IF (axis.LT.1.OR.axis.GT.3) &
         WRITE(stdout,*)'WARNING: wrong axis in generate_exponential'
    !
    ALLOCATE( rholocal( nnr ) )
    rholocal = 0.D0
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
       ! ... three dimensional indexes
       !
! BACKWARD COMPATIBILITY
! Compatible with QE-5.X QE-6.1.X
!       idx = idx0 + ir - 1
!       k   = idx / (dfftp%nr1x*dfftp%nr2x)
!       idx = idx - (dfftp%nr1x*dfftp%nr2x)*k
!       j   = idx / dfftp%nr1x
!       idx = idx - dfftp%nr1x*j
!       i   = idx
! Compatible with QE-6.2, QE-6.2.1 and QE-GIT
       idx = ir - 1
       k   = idx / (dfftp%nr1x*dfftp%my_nr2p)
       idx = idx - (dfftp%nr1x*dfftp%my_nr2p)*k
       k   = k + k0
       j   = idx / dfftp%nr1x
       idx = idx - dfftp%nr1x * j
       j   = j + j0
       i   = idx
! END BACKWARD COMPATIBILITY
       !
       ! ... do not include points outside the physical range
       !
       IF ( i >= dfftp%nr1 .OR. j >= dfftp%nr2 .OR. k >= dfftp%nr3 ) CYCLE
       !
       DO ip = 1, 3
          r(ip) = DBLE( i )*inv_nr1*at(ip,1) + &
               DBLE( j )*inv_nr2*at(ip,2) + &
               DBLE( k )*inv_nr3*at(ip,3)
       END DO
       !
       r(:) = pos(:) - r(:)
       !
       !  ... possibly 2D or 1D gaussians
       !
       IF ( dim .EQ. 1) THEN
          r(axis) = 0.D0
       ELSE IF ( dim .EQ. 2 ) THEN
          DO i = 1, 3
             IF ( i .NE. axis ) r(i) = 0.D0
          ENDDO
       END IF
       !
       ! ... minimum image convention
       !
       s(:) = MATMUL( r(:), bg(:,:) )
       s(:) = s(:) - ANINT(s(:))
       r(:) = MATMUL( at(:,:), s(:) )
       !
       dist = SQRT(SUM( r * r )) * alat
       arg = ( dist - width ) / spread
       !
       IF( ABS( arg ) .LT. exp_arg_limit ) THEN
          rholocal( ir ) = EXP( - arg )
       ELSE
          rholocal( ir ) = 0.D0
       END IF
       !
    END DO
    !
    rho = rho + rholocal
    DEALLOCATE( rholocal )
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE generate_exponential
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE generate_gradexponential(nnr, dim, axis, width, spread, pos, gradrho )
!--------------------------------------------------------------------
    !
    USE kinds,            ONLY : DP
    USE io_global,        ONLY : stdout
    USE cell_base,        ONLY : at, bg, alat
    USE fft_base,         ONLY : dfftp
    USE mp_bands,         ONLY : me_bgrp, intra_bgrp_comm
    !
    IMPLICIT NONE
    !
    ! ... Declares variables
    !
    INTEGER, INTENT(IN)       :: nnr, dim, axis
    REAL( DP ), INTENT(IN)    :: width, spread
    REAL( DP ), INTENT(IN)    :: pos( 3 )
    REAL( DP ), INTENT(INOUT) :: gradrho( 3, nnr )
    !
    ! ... Local variables
    !
    INTEGER                   :: i, j, j0, k, k0, ir, ir_end, ip
    INTEGER                   :: idx, idx0
    !
    REAL( DP )                :: inv_nr1, inv_nr2, inv_nr3
    REAL( DP )                :: dist, arg
    REAL( DP )                :: r( 3 ), s( 3 )
    REAL( DP ), ALLOCATABLE   :: gradrholocal ( :, : )
    REAL( DP ), PARAMETER     :: exp_arg_limit = 100.D0, tol = 1.D-10
    !
    IF ( dfftp%nr1 .EQ. 0 .OR. dfftp%nr2 .EQ. 0 .OR. dfftp%nr3 .EQ. 0 ) THEN
       WRITE(stdout,*)'ERROR: wrong grid dimension',dfftp%nr1,dfftp%nr2,dfftp%nr3
       STOP
    ENDIF
    inv_nr1 = 1.D0 / DBLE( dfftp%nr1 )
    inv_nr2 = 1.D0 / DBLE( dfftp%nr2 )
    inv_nr3 = 1.D0 / DBLE( dfftp%nr3 )
    !
    IF (axis.LT.1.OR.axis.GT.3) &
         WRITE(stdout,*)'WARNING: wrong axis in generate_gradexponential'
    !
    ALLOCATE( gradrholocal( 3, nnr ) )
    gradrholocal = 0.D0
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
       ! ... three dimensional indexes
       !
! BACKWARD COMPATIBILITY
! Compatible with QE-5.X QE-6.1.X
!       idx = idx0 + ir - 1
!       k   = idx / (dfftp%nr1x*dfftp%nr2x)
!       idx = idx - (dfftp%nr1x*dfftp%nr2x)*k
!       j   = idx / dfftp%nr1x
!       idx = idx - dfftp%nr1x*j
!       i   = idx
! Compatible with QE-6.2, QE-6.2.1 and QE-GIT
       idx = ir - 1
       k   = idx / (dfftp%nr1x*dfftp%my_nr2p)
       idx = idx - (dfftp%nr1x*dfftp%my_nr2p)*k
       k   = k + k0
       j   = idx / dfftp%nr1x
       idx = idx - dfftp%nr1x * j
       j   = j + j0
       i   = idx
! END BACKWARD COMPATIBILITY
       !
       ! ... do not include points outside the physical range
       !
       IF ( i >= dfftp%nr1 .OR. j >= dfftp%nr2 .OR. k >= dfftp%nr3 ) CYCLE
       !
       DO ip = 1, 3
          r(ip) = DBLE( i )*inv_nr1*at(ip,1) + &
               DBLE( j )*inv_nr2*at(ip,2) + &
               DBLE( k )*inv_nr3*at(ip,3)
       END DO
       !
       r(:) = pos(:) - r(:)
       !
       !  ... possibly 2D or 1D erfc
       !
       IF ( dim .EQ. 1) THEN
          r(axis) = 0.D0
       ELSE IF ( dim .EQ. 2 ) THEN
          DO i = 1, 3
             IF ( i .NE. axis ) r(i) = 0.D0
          ENDDO
       END IF
       !
       ! ... minimum image convention
       !
       s(:) = MATMUL( r(:), bg(:,:) )
       s(:) = s(:) - ANINT(s(:))
       r(:) = MATMUL( at(:,:), s(:) )
       !
       dist = SQRT(SUM( r * r )) * alat
       arg = ( dist - width ) / spread
       !
       IF ( dist .GT. tol .AND. ABS( arg ) .LT. exp_arg_limit ) THEN
          gradrholocal( :, ir ) = r(:) * alat / dist / spread * EXP( - arg )
       ELSE
          gradrholocal( :, ir ) = 0.D0
       ENDIF
       !
    END DO
    !
    gradrho = gradrho + gradrholocal
    DEALLOCATE( gradrholocal )
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE generate_gradexponential
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE generate_erfc( nnr, dim, axis, charge, width, spread, pos, rho )
!--------------------------------------------------------------------
    !
    USE kinds,            ONLY : DP
    USE io_global,        ONLY : stdout
    USE cell_base,        ONLY : at, bg, alat, omega
    USE fft_base,         ONLY : dfftp
    USE mp,               ONLY : mp_sum
    USE mp_bands,         ONLY : me_bgrp, intra_bgrp_comm
    !
    IMPLICIT NONE
    !
    ! ... Declares variables
    !
    INTEGER, INTENT(IN)       :: nnr, dim, axis
    REAL( DP ), INTENT(IN)    :: charge, width, spread
    REAL( DP ), INTENT(IN)    :: pos( 3 )
    REAL( DP ), INTENT(INOUT) :: rho( nnr )
    !
    ! ... Local variables
    !
    INTEGER                   :: i, j, j0, k, k0, ir, ir_end, ip
    INTEGER                   :: idx, idx0, ntot
    !
    REAL( DP )                :: inv_nr1, inv_nr2, inv_nr3
    REAL( DP )                :: scale, dist, arg, chargeanalytic, chargelocal
    REAL( DP )                :: r( 3 ), s( 3 )
    REAL( DP ), ALLOCATABLE   :: rholocal ( : )
    REAL( DP ), EXTERNAL      :: qe_erfc
    !
    IF ( dfftp%nr1 .EQ. 0 .OR. dfftp%nr2 .EQ. 0 .OR. dfftp%nr3 .EQ. 0 ) THEN
       WRITE(stdout,*)'ERROR: wrong grid dimension',dfftp%nr1,dfftp%nr2,dfftp%nr3
       STOP
    ENDIF
    inv_nr1 = 1.D0 / DBLE( dfftp%nr1 )
    inv_nr2 = 1.D0 / DBLE( dfftp%nr2 )
    inv_nr3 = 1.D0 / DBLE( dfftp%nr3 )
    !
    ntot = dfftp%nr1 * dfftp%nr2 * dfftp%nr3
    !
    IF (axis.LT.1.OR.axis.GT.3) &
         WRITE(stdout,*)'WARNING: wrong axis in generate_erfc'
    chargeanalytic = erfcvolume(dim,axis,width,spread,alat,omega,at)
    scale = charge / chargeanalytic * 0.5D0
    !
    ALLOCATE( rholocal( nnr ) )
    rholocal = 0.D0
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
       ! ... three dimensional indexes
       !
! BACKWARD COMPATIBILITY
! Compatible with QE-5.X QE-6.1.X
!       idx = idx0 + ir - 1
!       k   = idx / (dfftp%nr1x*dfftp%nr2x)
!       idx = idx - (dfftp%nr1x*dfftp%nr2x)*k
!       j   = idx / dfftp%nr1x
!       idx = idx - dfftp%nr1x*j
!       i   = idx
! Compatible with QE-6.2, QE-6.2.1 and QE-GIT
       idx = ir - 1
       k   = idx / (dfftp%nr1x*dfftp%my_nr2p)
       idx = idx - (dfftp%nr1x*dfftp%my_nr2p)*k
       k   = k + k0
       j   = idx / dfftp%nr1x
       idx = idx - dfftp%nr1x * j
       j   = j + j0
       i   = idx
! END BACKWARD COMPATIBILITY
       !
       ! ... do not include points outside the physical range
       !
       IF ( i >= dfftp%nr1 .OR. j >= dfftp%nr2 .OR. k >= dfftp%nr3 ) CYCLE
       !
       DO ip = 1, 3
          r(ip) = DBLE( i )*inv_nr1*at(ip,1) + &
               DBLE( j )*inv_nr2*at(ip,2) + &
               DBLE( k )*inv_nr3*at(ip,3)
       END DO
       !
       r(:) = r(:) - pos(:)
       !
       !  ... possibly 2D or 1D gaussians
       !
       IF ( dim .EQ. 1) THEN
          r(axis) = 0.D0
       ELSE IF ( dim .EQ. 2 ) THEN
          DO i = 1, 3
             IF ( i .NE. axis ) r(i) = 0.D0
          ENDDO
       END IF
       !
       ! ... minimum image convention
       !
       s(:) = MATMUL( r(:), bg(:,:) )
       s(:) = s(:) - ANINT(s(:))
       r(:) = MATMUL( at(:,:), s(:) )
       !
       dist = SQRT(SUM( r * r ))
       arg = ( dist * alat - width ) / spread
       !
       rholocal( ir ) = qe_erfc(arg)
       !
    END DO
    !
    ! ... double check that the integral of the generated charge corresponds to
    !     what is expected
    !
    chargelocal = SUM(rholocal)*omega/DBLE(ntot)*0.5D0
    CALL mp_sum(chargelocal,intra_bgrp_comm)
    IF ( ABS(chargelocal-chargeanalytic)/chargeanalytic .GT. 1.D-4 ) &
         WRITE(stdout,*)'WARNING: significant discrepancy between the numerical and the expected erfc charge'
    !
    ! ... rescale generated charge to obtain the correct integrated total charge
    !
    rholocal = rholocal * scale
    !
    rho = rho + rholocal
    DEALLOCATE( rholocal )
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE generate_erfc
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE generate_graderfc( nnr, dim, axis, charge, width, spread, pos, gradrho )
!--------------------------------------------------------------------
    !
    USE kinds,            ONLY : DP
    USE constants,        ONLY : sqrtpi
    USE io_global,        ONLY : stdout
    USE cell_base,        ONLY : at, bg, alat, omega
    USE fft_base,         ONLY : dfftp
    USE mp,               ONLY : mp_sum
    USE mp_bands,         ONLY : me_bgrp, intra_bgrp_comm
    !
    IMPLICIT NONE
    !
    REAL( DP ), PARAMETER :: tol = 1.D-10
    !
    ! ... Declares variables
    !
    INTEGER, INTENT(IN)       :: nnr, dim, axis
    REAL( DP ), INTENT(IN)    :: charge, width, spread
    REAL( DP ), INTENT(IN)    :: pos( 3 )
    REAL( DP ), INTENT(INOUT) :: gradrho( 3, nnr )
    !
    ! ... Local variables
    !
    INTEGER                   :: i, j, j0, k, k0, ir, ir_end, ip
    INTEGER                   :: idx, idx0, ntot
    !
    REAL( DP )                :: inv_nr1, inv_nr2, inv_nr3
    REAL( DP )                :: scale, dist, arg, chargeanalytic, chargelocal
    REAL( DP )                :: r( 3 ), s( 3 )
    REAL( DP ), ALLOCATABLE   :: gradrholocal ( :, : )
    REAL( DP ), EXTERNAL      :: qe_erfc
    !
    IF ( dfftp%nr1 .EQ. 0 .OR. dfftp%nr2 .EQ. 0 .OR. dfftp%nr3 .EQ. 0 ) THEN
       WRITE(stdout,*)'ERROR: wrong grid dimension',dfftp%nr1,dfftp%nr2,dfftp%nr3
       STOP
    ENDIF
    inv_nr1 = 1.D0 / DBLE( dfftp%nr1 )
    inv_nr2 = 1.D0 / DBLE( dfftp%nr2 )
    inv_nr3 = 1.D0 / DBLE( dfftp%nr3 )
    !
    ntot = dfftp%nr1 * dfftp%nr2 * dfftp%nr3
    !
    IF (axis.LT.1.OR.axis.GT.3) &
         WRITE(stdout,*)'WARNING: wrong axis in generate_gaussian'
    chargeanalytic = erfcvolume(dim,axis,width,spread,alat,omega,at)
    !
    ! ... scaling factor, take into account rescaling of generated density
    !     to obtain the correct integrated total charge
    !
    scale = charge / chargeanalytic / sqrtpi / spread
    !
    ALLOCATE( gradrholocal( 3, nnr ) )
    gradrholocal = 0.D0
    chargelocal = 0.D0
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
       ! ... three dimensional indexes
       !
! BACKWARD COMPATIBILITY
! Compatible with QE-5.X QE-6.1.X
!       idx = idx0 + ir - 1
!       k   = idx / (dfftp%nr1x*dfftp%nr2x)
!       idx = idx - (dfftp%nr1x*dfftp%nr2x)*k
!       j   = idx / dfftp%nr1x
!       idx = idx - dfftp%nr1x*j
!       i   = idx
! Compatible with QE-6.2, QE-6.2.1 and QE-GIT
       idx = ir - 1
       k   = idx / (dfftp%nr1x*dfftp%my_nr2p)
       idx = idx - (dfftp%nr1x*dfftp%my_nr2p)*k
       k   = k + k0
       j   = idx / dfftp%nr1x
       idx = idx - dfftp%nr1x * j
       j   = j + j0
       i   = idx
! END BACKWARD COMPATIBILITY
       !
       ! ... do not include points outside the physical range
       !
       IF ( i >= dfftp%nr1 .OR. j >= dfftp%nr2 .OR. k >= dfftp%nr3 ) CYCLE
       !
       DO ip = 1, 3
          r(ip) = DBLE( i )*inv_nr1*at(ip,1) + &
               DBLE( j )*inv_nr2*at(ip,2) + &
               DBLE( k )*inv_nr3*at(ip,3)
       END DO
       !
       r(:) = r(:) - pos(:)
       !
       !  ... possibly 2D or 1D erfc
       !
       IF ( dim .EQ. 1) THEN
          r(axis) = 0.D0
       ELSE IF ( dim .EQ. 2 ) THEN
          DO i = 1, 3
             IF ( i .NE. axis ) r(i) = 0.D0
          ENDDO
       END IF
       !
       ! ... minimum image convention
       !
       s(:) = MATMUL( r(:), bg(:,:) )
       s(:) = s(:) - ANINT(s(:))
       r(:) = MATMUL( at(:,:), s(:) )
       r = r * alat
       !
       dist = SQRT(SUM( r * r ))
       arg = ( dist - width ) / spread
       !
       IF ( dist .GT. tol ) gradrholocal( :, ir ) = - EXP( - arg**2 ) * r(:) / dist
       chargelocal = chargelocal + qe_erfc(arg)
       !
    END DO
    !
    ! ... double check that the integral of the generated charge corresponds to
    !     what is expected
    !
    CALL mp_sum( chargelocal, intra_bgrp_comm )
    chargelocal = chargelocal*omega/DBLE(ntot)*0.5D0
    IF ( ABS(chargelocal-chargeanalytic)/chargeanalytic .GT. 1.D-4 ) &
         WRITE(stdout,*)'WARNING: significant discrepancy between the numerical and the expected erfc charge'
    !
    gradrholocal = gradrholocal * scale
    !
    gradrho = gradrho + gradrholocal
    DEALLOCATE( gradrholocal )
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE generate_graderfc
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE generate_laplerfc( nnr, dim, axis, charge, width, spread, pos, laplrho )
!--------------------------------------------------------------------
    !
    USE kinds,            ONLY : DP
    USE constants,        ONLY : sqrtpi
    USE io_global,        ONLY : stdout
    USE cell_base,        ONLY : at, bg, alat, omega
    USE fft_base,         ONLY : dfftp
    USE mp,               ONLY : mp_sum
    USE mp_bands,         ONLY : me_bgrp, intra_bgrp_comm
    !
    IMPLICIT NONE
    !
    REAL( DP ), PARAMETER :: tol = 1.D-10
    !
    ! ... Declares variables
    !
    INTEGER, INTENT(IN)       :: nnr, dim, axis
    REAL( DP ), INTENT(IN)    :: charge, width, spread
    REAL( DP ), INTENT(IN)    :: pos( 3 )
    REAL( DP ), INTENT(INOUT) :: laplrho( nnr )
    !
    ! ... Local variables
    !
    INTEGER                   :: i, j, j0, k, k0, ir, ir_end, ip
    INTEGER                   :: idx, idx0, ntot
    !
    REAL( DP )                :: inv_nr1, inv_nr2, inv_nr3
    REAL( DP )                :: scale, dist, arg, chargeanalytic, chargelocal
    REAL( DP )                :: r( 3 ), s( 3 )
    REAL( DP ), ALLOCATABLE   :: laplrholocal ( : )
    REAL( DP ), EXTERNAL      :: qe_erfc
    !
    IF ( dfftp%nr1 .EQ. 0 .OR. dfftp%nr2 .EQ. 0 .OR. dfftp%nr3 .EQ. 0 ) THEN
       WRITE(stdout,*)'ERROR: wrong grid dimension',dfftp%nr1,dfftp%nr2,dfftp%nr3
       STOP
    ENDIF
    inv_nr1 = 1.D0 / DBLE( dfftp%nr1 )
    inv_nr2 = 1.D0 / DBLE( dfftp%nr2 )
    inv_nr3 = 1.D0 / DBLE( dfftp%nr3 )
    !
    ntot = dfftp%nr1 * dfftp%nr2 * dfftp%nr3
    !
    IF (axis.LT.1.OR.axis.GT.3) &
         WRITE(stdout,*)'WARNING: wrong axis in generate_gaussian'
    chargeanalytic = erfcvolume(dim,axis,width,spread,alat,omega,at)
    !
    ! ... scaling factor, take into account rescaling of generated density
    !     to obtain the correct integrated total charge
    !
    scale = charge / chargeanalytic / sqrtpi / spread
    !
    ALLOCATE( laplrholocal( nnr ) )
    laplrholocal = 0.D0
    chargelocal = 0.D0
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
       ! ... three dimensional indexes
       !
! BACKWARD COMPATIBILITY
! Compatible with QE-5.X QE-6.1.X
!       idx = idx0 + ir - 1
!       k   = idx / (dfftp%nr1x*dfftp%nr2x)
!       idx = idx - (dfftp%nr1x*dfftp%nr2x)*k
!       j   = idx / dfftp%nr1x
!       idx = idx - dfftp%nr1x*j
!       i   = idx
! Compatible with QE-6.2, QE-6.2.1 and QE-GIT
       idx = ir - 1
       k   = idx / (dfftp%nr1x*dfftp%my_nr2p)
       idx = idx - (dfftp%nr1x*dfftp%my_nr2p)*k
       k   = k + k0
       j   = idx / dfftp%nr1x
       idx = idx - dfftp%nr1x * j
       j   = j + j0
       i   = idx
! END BACKWARD COMPATIBILITY
       !
       ! ... do not include points outside the physical range
       !
       IF ( i >= dfftp%nr1 .OR. j >= dfftp%nr2 .OR. k >= dfftp%nr3 ) CYCLE
       !
       DO ip = 1, 3
          r(ip) = DBLE( i )*inv_nr1*at(ip,1) + &
               DBLE( j )*inv_nr2*at(ip,2) + &
               DBLE( k )*inv_nr3*at(ip,3)
       END DO
       !
       r(:) = r(:) - pos(:)
       !
       !  ... possibly 2D or 1D erfc
       !
       IF ( dim .EQ. 1) THEN
          r(axis) = 0.D0
       ELSE IF ( dim .EQ. 2 ) THEN
          DO i = 1, 3
             IF ( i .NE. axis ) r(i) = 0.D0
          ENDDO
       END IF
       !
       ! ... minimum image convention
       !
       s(:) = MATMUL( r(:), bg(:,:) )
       s(:) = s(:) - ANINT(s(:))
       r(:) = MATMUL( at(:,:), s(:) )
       r = r * alat
       !
       dist = SQRT(SUM( r * r ))
       arg = ( dist - width ) / spread
       !
       IF ( dist .GT. tol ) laplrholocal( ir ) = - EXP( - arg**2 ) * ( 1.D0 / dist - arg / spread ) * 2.D0
       chargelocal = chargelocal + qe_erfc(arg)
       !
    END DO
    !
    ! ... double check that the integral of the generated charge corresponds to
    !     what is expected
    !
    CALL mp_sum( chargelocal, intra_bgrp_comm )
    chargelocal = chargelocal*omega/DBLE(ntot)*0.5D0
    IF ( ABS(chargelocal-chargeanalytic)/chargeanalytic .GT. 1.D-4 ) &
         WRITE(stdout,*)'WARNING: significant discrepancy between the numerical and the expected erfc charge'
    !
    laplrholocal = laplrholocal * scale
    !
    laplrho = laplrho + laplrholocal
    DEALLOCATE( laplrholocal )
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE generate_laplerfc
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE generate_hesserfc( nnr, dim, axis, charge, width, spread, pos, hessrho )
!--------------------------------------------------------------------
    !
    USE kinds,            ONLY : DP
    USE constants,        ONLY : sqrtpi
    USE io_global,        ONLY : stdout
    USE cell_base,        ONLY : at, bg, alat, omega
    USE fft_base,         ONLY : dfftp
    USE mp,               ONLY : mp_sum
    USE mp_bands,         ONLY : me_bgrp, intra_bgrp_comm
    !
    IMPLICIT NONE
    !
    REAL( DP ), PARAMETER :: tol = 1.D-10
    !
    ! ... Declares variables
    !
    INTEGER, INTENT(IN)       :: nnr, dim, axis
    REAL( DP ), INTENT(IN)    :: charge, width, spread
    REAL( DP ), INTENT(IN)    :: pos( 3 )
    REAL( DP ), INTENT(INOUT) :: hessrho( 3, 3, nnr )
    !
    ! ... Local variables
    !
    INTEGER                   :: i, j, j0, k, k0, ir, ir_end, ip, jp
    INTEGER                   :: idx, idx0, ntot
    !
    REAL( DP )                :: inv_nr1, inv_nr2, inv_nr3
    REAL( DP )                :: scale, dist, arg, chargeanalytic, chargelocal, tmp
    REAL( DP )                :: r( 3 ), s( 3 )
    REAL( DP ), ALLOCATABLE   :: hessrholocal ( :, :, : )
    REAL( DP ), EXTERNAL      :: qe_erfc
    !
    IF ( dfftp%nr1 .EQ. 0 .OR. dfftp%nr2 .EQ. 0 .OR. dfftp%nr3 .EQ. 0 ) THEN
       WRITE(stdout,*)'ERROR: wrong grid dimension',dfftp%nr1,dfftp%nr2,dfftp%nr3
       STOP
    ENDIF
    inv_nr1 = 1.D0 / DBLE( dfftp%nr1 )
    inv_nr2 = 1.D0 / DBLE( dfftp%nr2 )
    inv_nr3 = 1.D0 / DBLE( dfftp%nr3 )
    !
    ntot = dfftp%nr1 * dfftp%nr2 * dfftp%nr3
    !
    IF (axis.LT.1.OR.axis.GT.3) &
         WRITE(stdout,*)'WARNING: wrong axis in generate_gaussian'
    chargeanalytic = erfcvolume(dim,axis,width,spread,alat,omega,at)
    !
    ! ... scaling factor, take into account rescaling of generated density
    !     to obtain the correct integrated total charge
    !
    scale = charge / chargeanalytic / sqrtpi / spread
    !
    ALLOCATE( hessrholocal( 3, 3, nnr ) )
    hessrholocal = 0.D0
    chargelocal = 0.D0
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
       ! ... three dimensional indexes
       !
! BACKWARD COMPATIBILITY
! Compatible with QE-5.X QE-6.1.X
!       idx = idx0 + ir - 1
!       k   = idx / (dfftp%nr1x*dfftp%nr2x)
!       idx = idx - (dfftp%nr1x*dfftp%nr2x)*k
!       j   = idx / dfftp%nr1x
!       idx = idx - dfftp%nr1x*j
!       i   = idx
! Compatible with QE-6.2, QE-6.2.1 and QE-GIT
       idx = ir - 1
       k   = idx / (dfftp%nr1x*dfftp%my_nr2p)
       idx = idx - (dfftp%nr1x*dfftp%my_nr2p)*k
       k   = k + k0
       j   = idx / dfftp%nr1x
       idx = idx - dfftp%nr1x * j
       j   = j + j0
       i   = idx
! END BACKWARD COMPATIBILITY
       !
       ! ... do not include points outside the physical range
       !
       IF ( i >= dfftp%nr1 .OR. j >= dfftp%nr2 .OR. k >= dfftp%nr3 ) CYCLE
       !
       DO ip = 1, 3
          r(ip) = DBLE( i )*inv_nr1*at(ip,1) + &
               DBLE( j )*inv_nr2*at(ip,2) + &
               DBLE( k )*inv_nr3*at(ip,3)
       END DO
       !
       r(:) = r(:) - pos(:)
       !
       !  ... possibly 2D or 1D erfc
       !
       IF ( dim .EQ. 1) THEN
          r(axis) = 0.D0
       ELSE IF ( dim .EQ. 2 ) THEN
          DO i = 1, 3
             IF ( i .NE. axis ) r(i) = 0.D0
          ENDDO
       END IF
       !
       ! ... minimum image convention
       !
       s(:) = MATMUL( r(:), bg(:,:) )
       s(:) = s(:) - ANINT(s(:))
       r(:) = MATMUL( at(:,:), s(:) )
       r = r * alat
       !
       dist = SQRT(SUM( r * r ))
       arg = ( dist - width ) / spread
       !
       IF ( dist .GT. tol ) THEN
          DO ip = 1, 3
             DO jp = 1, 3
                tmp = - r(ip) * r(jp) * ( 1.D0 / dist + 2.D0 * arg / spread )
                IF ( ip .EQ. jp ) tmp = tmp + dist
                hessrholocal( ip, jp, ir ) = - EXP( - arg**2 ) * tmp / dist**2
             ENDDO
          ENDDO
       END IF
       chargelocal = chargelocal + qe_erfc(arg)
       !
    END DO
    !
    ! ... double check that the integral of the generated charge corresponds to
    !     what is expected
    !
    CALL mp_sum( chargelocal, intra_bgrp_comm )
    chargelocal = chargelocal*omega/DBLE(ntot)*0.5D0
    IF ( ABS(chargelocal-chargeanalytic)/chargeanalytic .GT. 1.D-4 ) &
         WRITE(stdout,*)'WARNING: significant discrepancy between the numerical and the expected erfc charge'
    !
    hessrholocal = hessrholocal * scale
    !
    hessrho = hessrho + hessrholocal
    DEALLOCATE( hessrholocal )
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE generate_hesserfc
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE generate_axis( nnr, icor, pos, axis )
!--------------------------------------------------------------------
    USE kinds,            ONLY : DP
    USE io_global,        ONLY : stdout
    USE cell_base,        ONLY : at, bg, alat
    USE fft_base,         ONLY : dfftp
    USE mp_bands,         ONLY : me_bgrp, intra_bgrp_comm
    !
    INTEGER, INTENT(IN) :: nnr
    INTEGER, INTENT(IN) :: icor
    REAL(DP), INTENT(IN) :: pos(3)
    REAL(DP), INTENT(OUT) :: axis( dfftp%nnr )
    !
    INTEGER  :: i, j, j0, k, k0, ir, ir_end, ip, idx, idx0
    REAL(DP) :: inv_nr1, inv_nr2, inv_nr3
    REAL(DP) :: r(3)
    !
    IF ( dfftp%nr1 .EQ. 0 .OR. dfftp%nr2 .EQ. 0 .OR. dfftp%nr3 .EQ. 0 ) THEN
       WRITE(stdout,*)'ERROR: wrong grid dimension',dfftp%nr1,dfftp%nr2,dfftp%nr3
       STOP
    ENDIF
    inv_nr1 = 1.D0 / DBLE( dfftp%nr1 )
    inv_nr2 = 1.D0 / DBLE( dfftp%nr2 )
    inv_nr3 = 1.D0 / DBLE( dfftp%nr3 )
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
       ! ... three dimensional indexes
       !
! BACKWARD COMPATIBILITY
! Compatible with QE-5.X QE-6.1.X
!       idx = idx0 + ir - 1
!       k   = idx / (dfftp%nr1x*dfftp%nr2x)
!       idx = idx - (dfftp%nr1x*dfftp%nr2x)*k
!       j   = idx / dfftp%nr1x
!       idx = idx - dfftp%nr1x*j
!       i   = idx
! Compatible with QE-6.2, QE-6.2.1 and QE-GIT
       idx = ir - 1
       k   = idx / (dfftp%nr1x*dfftp%my_nr2p)
       idx = idx - (dfftp%nr1x*dfftp%my_nr2p)*k
       k   = k + k0
       j   = idx / dfftp%nr1x
       idx = idx - dfftp%nr1x * j
       j   = j + j0
       i   = idx
! END BACKWARD COMPATIBILITY
       !
       ! ... do not include points outside the physical range
       !
       IF ( i >= dfftp%nr1 .OR. j >= dfftp%nr2 .OR. k >= dfftp%nr3 ) CYCLE
       !
       DO ip = 1, 3
          r(ip) = DBLE( i )*inv_nr1*at(ip,1) + &
               DBLE( j )*inv_nr2*at(ip,2) + &
               DBLE( k )*inv_nr3*at(ip,3)
       END DO
       !
       r(:) = r(:) - pos(:)
       !
       ! ... minimum image convention
       !
       CALL cryst_to_cart( 1, r, bg, -1 )
       !
       r(:) = r(:) - ANINT( r(:) )
       !
       CALL cryst_to_cart( 1, r, at, 1 )
       !
       axis(ir) = r(icor)
       !
    END DO
    !
    axis = axis * alat
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE generate_axis
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE generate_distance( nnr, pos, distance )
!--------------------------------------------------------------------
    USE kinds,            ONLY : DP
    USE io_global,        ONLY : stdout
    USE cell_base,        ONLY : at, bg, alat
    USE fft_base,         ONLY : dfftp
    USE mp_bands,         ONLY : me_bgrp, intra_bgrp_comm
    !
    INTEGER, INTENT(IN) :: nnr
    REAL(DP), INTENT(IN) :: pos(3)
    REAL(DP), INTENT(OUT) :: distance( 3, dfftp%nnr )
    !
    INTEGER  :: i, j, j0, k, k0, ir, ir_end, ip, idx, idx0
    REAL(DP) :: inv_nr1, inv_nr2, inv_nr3
    REAL(DP) :: r(3), s(3)
    !
    IF ( dfftp%nr1 .EQ. 0 .OR. dfftp%nr2 .EQ. 0 .OR. dfftp%nr3 .EQ. 0 ) THEN
       WRITE(stdout,*)'ERROR: wrong grid dimension',dfftp%nr1,dfftp%nr2,dfftp%nr3
       STOP
    ENDIF
    inv_nr1 = 1.D0 / DBLE( dfftp%nr1 )
    inv_nr2 = 1.D0 / DBLE( dfftp%nr2 )
    inv_nr3 = 1.D0 / DBLE( dfftp%nr3 )
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
       ! ... three dimensional indexes
       !
! BACKWARD COMPATIBILITY
! Compatible with QE-5.X QE-6.1.X
!       idx = idx0 + ir - 1
!       k   = idx / (dfftp%nr1x*dfftp%nr2x)
!       idx = idx - (dfftp%nr1x*dfftp%nr2x)*k
!       j   = idx / dfftp%nr1x
!       idx = idx - dfftp%nr1x*j
!       i   = idx
! Compatible with QE-6.2, QE-6.2.1 and QE-GIT
       idx = ir - 1
       k   = idx / (dfftp%nr1x*dfftp%my_nr2p)
       idx = idx - (dfftp%nr1x*dfftp%my_nr2p)*k
       k   = k + k0
       j   = idx / dfftp%nr1x
       idx = idx - dfftp%nr1x * j
       j   = j + j0
       i   = idx
! END BACKWARD COMPATIBILITY
       !
       ! ... do not include points outside the physical range
       !
       IF ( i >= dfftp%nr1 .OR. j >= dfftp%nr2 .OR. k >= dfftp%nr3 ) CYCLE
       !
       DO ip = 1, 3
          r(ip) = DBLE( i )*inv_nr1*at(ip,1) + &
               DBLE( j )*inv_nr2*at(ip,2) + &
               DBLE( k )*inv_nr3*at(ip,3)
       END DO
       !
       r(:) = r(:) - pos(:)
       !
       ! ... minimum image convention
       !
       CALL cryst_to_cart( 1, r, bg, -1 )
       !
       r(:) = r(:) - ANINT( r(:) )
       !
       CALL cryst_to_cart( 1, r, at, 1 )
       !
       distance(:,ir) = r(:)
       !
    END DO
    !
    distance = distance * alat
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE generate_distance
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  FUNCTION erfcvolume(dim,axis,width,spread,alat,omega,at)
!--------------------------------------------------------------------
    !
    USE constants,        ONLY : sqrtpi, fpi, pi
    USE io_global,        ONLY : stdout
    !
    REAL(DP), PARAMETER :: tol = 1.D-6
    !
    REAL(DP) :: erfcvolume
    !
    INTEGER, INTENT(IN) :: dim, axis
    REAL(DP), INTENT(IN) :: width, spread, alat, omega
    REAL(DP), DIMENSION(3,3), INTENT(IN) :: at
    !
    REAL(DP) :: f1 = 0.0_DP , f2 = 0.0_DP
    REAL(DP) :: t, invt
    REAL( DP ), EXTERNAL      :: qe_erf
    !
    IF ( spread .LT. tol .OR. width .LT. tol ) THEN
       WRITE(stdout,*)'ERROR: wrong parameters of erfc function',spread,width
       STOP
    ENDIF
    t = spread / width
    invt = width / spread
    f1 = ( 1.D0 + qe_erf(invt) ) / 2.D0 ! f1 is close to one  for t-->0
    f2 = exp(-(invt)**2) / 2.D0 / sqrtpi ! f2 is close to zero for t-->0
    SELECT CASE ( dim )
    CASE ( 0 )
       ! zero-dimensional erfc, volume is approx the one of the
       ! sphere of radius=width
       erfcvolume = fpi / 3.D0 * width**3 * &
            ( ( 1.D0 + 1.5D0 * t**2 ) * f1 + ( 1.D0 + t**2 ) * t * f2 )
    CASE ( 1 )
       ! one-dimensional erfc, volume is approx the one of the
       ! cylinder of radius=width and length=alat*at(axis,axis)
       erfcvolume = pi * width**2 * at(axis,axis) * alat * &
            ( ( 1.D0 + 0.5D0 * t**2 ) * f1  + t * f2 )
    CASE ( 2 )
       ! two-dimensional erfc, volume is exactly the one of the
       ! box, does not depend on spread
       erfcvolume = 2.D0 * width * omega / at(axis,axis) / alat
    END SELECT
    !
    RETURN
    !
!--------------------------------------------------------------------
  END FUNCTION erfcvolume
!--------------------------------------------------------------------
!----------------------------------------------------------------------------
END MODULE tools_generate_functions
!----------------------------------------------------------------------------


