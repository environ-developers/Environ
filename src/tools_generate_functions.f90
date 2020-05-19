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
!> Module to generate functions on the real space dense grid
!
! Authors: Oliviero Andreussi (Department of Physics, UNT)
!
!----------------------------------------------------------------------------
MODULE tools_generate_functions
!----------------------------------------------------------------------------
  !
  USE environ_types
  !
  IMPLICIT NONE
  !
  REAL(DP), PARAMETER       :: tol = 1.D-10
  REAL(DP), PARAMETER       :: exp_tol = 4.D1
  !
CONTAINS
!--------------------------------------------------------------------
  SUBROUTINE planar_average( cell, nnr, naxis, axis, shift, reverse, f, f1d )
!--------------------------------------------------------------------
    !
    USE env_mp,   ONLY : env_mp_sum
    !
    IMPLICIT NONE
    !
    ! ... Declares variables
    !
    TYPE( environ_cell ), INTENT( IN ) :: cell
    INTEGER, INTENT(IN)       :: nnr, naxis, axis, shift
    LOGICAL, INTENT(IN)       :: reverse
    REAL( DP ), INTENT(INOUT) :: f( nnr )
    REAL( DP ), INTENT(INOUT) :: f1d( naxis )
    !
    ! ... Local variables
    !
    INTEGER                   :: i, j, k, ir
    INTEGER                   :: idx, narea
    LOGICAL                   :: physical
    !
    narea = cell % ntot / naxis
    !
    IF ( reverse ) THEN
       f = 0.D0
    ELSE
       f1d = 0.D0
    END IF
    !
    DO ir = 1, cell%ir_end
       !
       ! ... three dimensional indexes
       !
       CALL ir2ijk( cell, ir, i, j, k, physical )
       !
       ! ... do not include points outside the physical range
       !
       IF ( .NOT. physical ) CYCLE
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
       CALL env_mp_sum( f1d(:), cell%comm )
       f1d = f1d / DBLE(narea)
    END IF
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE planar_average
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE generate_gaussian( dim, axis, charge, spread, pos, density )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    ! ... Declares variables
    !
    INTEGER, INTENT(IN)       :: dim, axis
    REAL( DP ), INTENT(IN)    :: charge, spread
    REAL( DP ), INTENT(IN)    :: pos( 3 )
    TYPE( environ_density), INTENT(INOUT) :: density
    !
    ! ... Local variables
    !
    LOGICAL                   :: physical
    INTEGER                   :: ir
    REAL( DP )                :: scale, spr2, length
    REAL( DP )                :: r( 3 ), r2
    REAL( DP ), ALLOCATABLE   :: local ( : )
    CHARACTER( LEN=80 )       :: sub_name = 'generate_gaussian'
    !
    ! ... Aliases
    !
    TYPE( environ_cell ), POINTER :: cell
    !
    ! ... sanity checks and initial setup
    !
    cell => density % cell
    !
    IF ( ABS( charge ) .LT. tol ) RETURN
    IF ( ABS( spread ) .LT. tol ) &
         & CALL env_errore(sub_name,'Wrong spread for Gaussian function',1)
    IF ( axis .LT. 1 .OR. axis .GT. 3 ) &
         & CALL env_errore(sub_name,'Wrong axis in generate_gaussian',1)
    !
    SELECT CASE ( dim )
    CASE ( 0 )
       scale = charge / ( sqrtpi * spread )**3
    CASE ( 1 )
       length = ABS( cell%at(axis,axis) * cell%alat )
       scale = charge / length / ( sqrtpi * spread )**2
    CASE ( 2 )
       length = ABS( cell%at(axis,axis) * cell%alat )
       scale = charge * length / cell%omega / ( sqrtpi * spread )
    CASE default
       CALL env_errore(sub_name,'Wrong value of dim',1)
    END SELECT
    !
    spr2 = ( spread / cell%alat )**2
    !
    ALLOCATE( local( cell%nnr ) )
    local = 0.D0
    !
    DO ir = 1, cell%ir_end
       !
       ! ... position in real space grid
       !
       CALL ir2r( cell, ir, r, physical )
       !
       ! ... do not include points outside the physical range
       !
       IF ( .NOT. physical ) CYCLE
       !
       ! ... displacement from origin
       !
       CALL displacement( dim, axis, pos, r, r )
       !
       ! ... minimum image convention
       !
       CALL minimum_image( cell, r, r2 )
       !
       ! ... compute Gaussian function
       !
       r2 = r2 / spr2
       !
       IF ( r2 .LE. exp_tol ) local( ir ) = EXP( -r2 )
       !
    END DO
    !
    density%of_r = density%of_r + scale * local
    DEALLOCATE( local )
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE generate_gaussian
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE generate_gradgaussian( dim, axis, charge, spread, pos, gradient )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    ! ... Declares variables
    !
    INTEGER, INTENT(IN)       :: dim, axis
    REAL( DP ), INTENT(IN)    :: charge, spread
    REAL( DP ), INTENT(IN)    :: pos( 3 )
    TYPE( environ_gradient ), INTENT(INOUT) :: gradient
    !
    ! ... Local variables
    !
    LOGICAL                   :: physical
    INTEGER                   :: ir
    REAL( DP )                :: scale, spr2, length
    REAL( DP )                :: r( 3 ), r2
    REAL( DP ), ALLOCATABLE   :: gradlocal ( :, : )
    CHARACTER( LEN=80 )       :: sub_name = 'generate_gradgaussian'
    !
    ! ... Aliases
    !
    TYPE( environ_cell ), POINTER :: cell
    !
    ! ... sanity checks and initial setup
    !
    cell => gradient % cell
    !
    IF ( ABS( charge ) .LT. tol ) RETURN
    IF ( ABS( spread ) .LT. tol ) &
         & CALL env_errore(sub_name,'Wrong spread for Gaussian function',1)
    IF ( axis .LT. 1 .OR. axis .GT. 3 ) &
         & CALL env_errore(sub_name,'Wrong value of axis',1)
    !
    SELECT CASE ( dim )
    CASE ( 0 )
       scale = charge / ( sqrtpi * spread )**3
    CASE ( 1 )
       length = ABS( cell%at(axis,axis) * cell%alat )
       scale = charge / length / ( sqrtpi * spread )**2
    CASE ( 2 )
       length = ABS( cell%at(axis,axis) * cell%alat )
       scale = charge * length / cell%omega / ( sqrtpi * spread )
    CASE default
       CALL env_errore(sub_name,'Wrong value of dim',1)
    END SELECT
    scale = scale * 2.D0 / spread**2 * cell%alat
    !
    spr2 = ( spread / cell%alat )**2
    !
    ALLOCATE( gradlocal( 3, cell%nnr ) )
    gradlocal = 0.D0
    !
    DO ir = 1, cell%ir_end
       !
       ! ... position in real space grid
       !
       CALL ir2r( cell, ir, r, physical )
       !
       ! ... do not include points outside the physical range
       !
       IF ( .NOT. physical ) CYCLE
       !
       ! ... displacement from origin
       !
       CALL displacement( dim, axis, pos, r, r )
       !
       ! ... minimum image convention
       !
       CALL minimum_image( cell, r, r2 )
       !
       ! ... compute gradient of Gaussian function
       !
       r2 = r2 / spr2
       !
       IF ( r2 .LE. exp_tol ) gradlocal( :, ir ) = EXP(-r2) * r(:)
       !
    END DO
    !
    gradient%of_r = gradient%of_r + scale * gradlocal
    DEALLOCATE( gradlocal )
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE generate_gradgaussian
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE generate_exponential( dim, axis, width, spread, pos, density )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    ! ... Declares variables
    !
    INTEGER, INTENT(IN)       :: dim, axis
    REAL( DP ), INTENT(IN)    :: width, spread
    REAL( DP ), INTENT(IN)    :: pos( 3 )
    TYPE( environ_density ), INTENT(INOUT) :: density
    !
    ! ... Local variables
    !
    LOGICAL                   :: physical
    INTEGER                   :: ir
    REAL( DP )                :: r2, dist, arg
    REAL( DP )                :: r( 3 )
    REAL( DP ), ALLOCATABLE   :: local ( : )
    CHARACTER( LEN=80 )       :: sub_name = 'generate_exponential'
    !
    ! ... Aliases
    !
    TYPE( environ_cell ), POINTER :: cell
    !
    cell => density % cell
    !
    IF ( axis .LT. 1 .OR. axis .GT. 3 ) &
         & CALL env_errore(sub_name,'Wrong value of axis',1)
    !
    ALLOCATE( local( cell%nnr ) )
    local = 0.D0
    !
    DO ir = 1, cell%ir_end
       !
       ! ... position in real space grid
       !
       CALL ir2r( cell, ir, r, physical )
       !
       ! ... do not include points outside the physical range
       !
       IF ( .NOT. physical ) CYCLE
       !
       ! ... displacement from origin
       !
       CALL displacement( dim, axis, pos, r, r )
       !
       ! ... minimum image convention
       !
       CALL minimum_image( cell, r, r2 )
       !
       ! ... compute exponentially decaying function
       !
       dist = SQRT(r2) * cell % alat
       !
       arg = ( dist - width ) / spread
       !
       IF( ABS( arg ) .LE. exp_tol ) local( ir ) = EXP( - arg )
       !
    END DO
    !
    density%of_r = density%of_r + local
    DEALLOCATE( local )
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE generate_exponential
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE generate_gradexponential( dim, axis, width, spread, pos, gradient )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    ! ... Declares variables
    !
    INTEGER, INTENT(IN)       :: dim, axis
    REAL( DP ), INTENT(IN)    :: width, spread
    REAL( DP ), INTENT(IN)    :: pos( 3 )
    TYPE( environ_gradient ), INTENT(INOUT) :: gradient
    !
    ! ... Local variables
    !
    LOGICAL                   :: physical
    INTEGER                   :: ir
    REAL( DP )                :: r2, dist, arg
    REAL( DP )                :: r( 3 )
    REAL( DP ), ALLOCATABLE   :: gradlocal ( :, : )
    CHARACTER( LEN=80 )       :: sub_name = 'generate_gradexponential'
    !
    ! ... Aliases
    !
    TYPE( environ_cell ), POINTER :: cell
    !
    cell => gradient % cell
    !
    IF ( axis .LT. 1 .OR. axis .GT. 3 ) &
         & CALL env_errore(sub_name,'Wrong value of axis',1)
    !
    ALLOCATE( gradlocal( 3, cell%nnr ) )
    gradlocal = 0.D0
    !
    DO ir = 1, cell%ir_end
       !
       ! ... position in real space grid
       !
       CALL ir2r( cell, ir, r, physical )
       !
       ! ... do not include points outside the physical range
       !
       IF ( .NOT. physical ) CYCLE
       !
       ! ... displacement from origin
       !
       CALL displacement( dim, axis, pos, r, r )
       !
       ! ... minimum image convention
       !
       CALL minimum_image( cell, r, r2 )
       !
       ! ... compute exponentially decaying function
       !
       dist = SQRT( r2 ) * cell % alat
       !
       arg = ( dist - width ) / spread
       !
       IF ( r2 .GT. tol .AND. ABS( arg ) .LE. exp_tol ) &
          gradlocal( :, ir ) = r(:) / SQRT(r2) / spread * EXP( - arg )
       !
    END DO
    !
    gradient%of_r = gradient%of_r + gradlocal
    DEALLOCATE( gradlocal )
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE generate_gradexponential
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE generate_erfc( dim, axis, charge, width, spread, pos, density )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    ! ... Declares variables
    !
    INTEGER, INTENT(IN)       :: dim, axis
    REAL( DP ), INTENT(IN)    :: charge, width, spread
    REAL( DP ), INTENT(IN)    :: pos( 3 )
    TYPE( environ_density ), INTENT(INOUT) :: density
    !
    ! ... Local variables
    !
    LOGICAL                   :: physical
    INTEGER                   :: ir, ir_end, i
    REAL( DP )                :: scale, r2, dist, arg, chargeanalytic, chargelocal
    REAL( DP )                :: r( 3 )
    REAL( DP ), ALLOCATABLE   :: local ( : )
    CHARACTER( LEN=80 )       :: sub_name = 'generate_erfc'
    REAL( DP ), EXTERNAL      :: qe_erfc
    !
    ! ... Aliases
    !
    TYPE( environ_cell ), POINTER :: cell
    !
    cell => density % cell
    !
    IF ( axis .LT. 1 .OR. axis .GT. 3 ) &
         & CALL env_errore(sub_name,'Wrong value of axis',1)
    !
    chargeanalytic = erfcvolume(dim,axis,width,spread,cell)
    scale = charge / chargeanalytic * 0.5D0
    !
    ALLOCATE( local( cell%nnr ) )
    local = 0.D0
    !
    DO ir = 1, cell%ir_end
       !
       ! ... position in real space grid
       !
       CALL ir2r( cell, ir, r, physical )
       !
       ! ... do not include points outside the physical range
       !
       IF ( .NOT. physical ) CYCLE
       !
       ! ... displacement from origin
       !
       CALL displacement( dim, axis, pos, r, r )
       !
       ! ... minimum image convention
       !
       CALL minimum_image( cell, r, r2 )
       !
       ! ... compute error function
       !
       dist = SQRT(r2) * cell%alat
       arg = ( dist  - width ) / spread
       !
       local( ir ) = qe_erfc(arg)
       !
    END DO
    !
    ! ... check integral of function is consistent with analytic one
    !
    chargelocal = SUM(local)*cell%omega/DBLE(cell%ntot)*0.5D0
    CALL env_mp_sum(chargelocal,cell%comm)
    IF ( ABS(chargelocal-chargeanalytic)/chargeanalytic .GT. 1.D-4 ) &
         CALL env_infomsg(sub_name,'WARNING: wrong integral of erfc function')
    !
    ! ... rescale generated function to obtain the requested integral
    !
    density%of_r = density%of_r + scale * local
    DEALLOCATE( local )
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE generate_erfc
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE generate_graderfc( dim, axis, charge, width, spread, pos, gradient )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    ! ... Declares variables
    !
    INTEGER, INTENT(IN)       :: dim, axis
    REAL( DP ), INTENT(IN)    :: charge, width, spread
    REAL( DP ), INTENT(IN)    :: pos( 3 )
    TYPE( environ_gradient ), INTENT(INOUT) :: gradient
    !
    ! ... Local variables
    !
    LOGICAL                   :: physical
    INTEGER                   :: ir
    REAL( DP )                :: scale, r2, dist, arg, chargeanalytic
    REAL( DP )                :: r( 3 )
    REAL( DP ), ALLOCATABLE   :: gradlocal ( :, : )
    CHARACTER( LEN=80 )       :: sub_name = 'generate_graderfc'
    !
    ! ... Aliases
    !
    TYPE( environ_cell ), POINTER :: cell
    !
    cell => gradient % cell
    !
    IF ( axis .LT. 1 .OR. axis .GT. 3 ) &
         & CALL env_errore(sub_name,'Wrong value of axis',1)
    !
    chargeanalytic = erfcvolume(dim,axis,width,spread,cell)
    !
    ! ... scaling factor, take into account rescaling of generated density
    !     to obtain the correct integrated total charge
    !
    scale = charge / chargeanalytic / sqrtpi / spread
    !
    ALLOCATE( gradlocal( 3, cell%nnr ) )
    gradlocal = 0.D0
    !
    DO ir = 1, cell%ir_end
       !
       ! ... position in real space grid
       !
       CALL ir2r( cell, ir, r, physical )
       !
       ! ... do not include points outside the physical range
       !
       IF ( .NOT. physical ) CYCLE
       !
       ! ... displacement from origin
       !
       CALL displacement( dim, axis, pos, r, r )
       !
       ! ... minimum image convention
       !
       CALL minimum_image( cell, r, r2 )
       !
       ! ... compute gradient of error function
       !
       r = r * cell % alat
       dist = SQRT(r2) * cell % alat
       arg = ( dist - width ) / spread
       !
       IF ( dist .GT. tol ) gradlocal( :, ir ) = EXP( - arg**2 ) * r(:) / dist
       !
    END DO
    !
    gradient%of_r = gradient%of_r + gradlocal * scale
    DEALLOCATE( gradlocal )
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE generate_graderfc
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE generate_laplerfc( dim, axis, charge, width, spread, pos, laplacian )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    ! ... Declares variables
    !
    INTEGER, INTENT(IN)       :: dim, axis
    REAL( DP ), INTENT(IN)    :: charge, width, spread
    REAL( DP ), INTENT(IN)    :: pos( 3 )
    TYPE( environ_density ), INTENT(INOUT) :: laplacian
    !
    ! ... Local variables
    !
    LOGICAL                   :: physical
    INTEGER                   :: ir
    REAL( DP )                :: scale, r2, dist, arg, chargeanalytic
    REAL( DP )                :: r( 3 )
    REAL( DP ), ALLOCATABLE   :: lapllocal ( : )
    CHARACTER( LEN=80 )       :: sub_name = 'generate_laplerfc'
    !
    ! ... Aliases
    !
    TYPE( environ_cell ), POINTER :: cell
    !
    cell => laplacian % cell
    !
    IF ( axis .LT. 1 .OR. axis .GT. 3 ) &
         & CALL env_errore(sub_name,'Wrong value of axis',1)
    !
    chargeanalytic = erfcvolume(dim,axis,width,spread,cell)
    !
    ! ... scaling factor, take into account rescaling of generated density
    !     to obtain the correct integrated total charge
    !
    scale = charge / chargeanalytic / sqrtpi / spread
    !
    ALLOCATE( lapllocal( cell%nnr ) )
    lapllocal = 0.D0
    !
    DO ir = 1, cell%ir_end
       !
       ! ... position in real space grid
       !
       CALL ir2r( cell, ir, r, physical )
       !
       ! ... do not include points outside the physical range
       !
       IF ( .NOT. physical ) CYCLE
       !
       ! ... displacement from origin
       !
       CALL displacement( dim, axis, pos, r, r )
       !
       ! ... minimum image convention
       !
       CALL minimum_image( cell, r, r2 )
       !
       ! ... compute laplacian of error function
       !
       dist = SQRT(r2) * cell%alat
       !
       arg = ( dist - width ) / spread
       !
       SELECT CASE ( dim )
       CASE( 0 )
         IF ( dist .GT. tol ) lapllocal( ir ) = - EXP( - arg**2 ) * ( 1.D0 / dist - arg / spread ) * 2.D0
       CASE( 1 )
         IF ( dist .GT. tol ) lapllocal( ir ) = - EXP( - arg**2 ) * ( 1.D0 / dist - 2.D0 * arg / spread )
       CASE( 2 )
         lapllocal( ir ) = EXP( - arg**2 ) * arg / spread  * 2.D0
       END SELECT
       !
    END DO
    !
    laplacian%of_r = laplacian%of_r + lapllocal * scale
    DEALLOCATE( lapllocal )
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE generate_laplerfc
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE generate_hesserfc( dim, axis, charge, width, spread, pos, hessian )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    ! ... Declares variables
    !
    INTEGER, INTENT(IN)       :: dim, axis
    REAL( DP ), INTENT(IN)    :: charge, width, spread
    REAL( DP ), INTENT(IN)    :: pos( 3 )
    TYPE( environ_hessian ), INTENT(INOUT) :: hessian
    !
    ! ... Local variables
    !
    LOGICAL                   :: physical
    INTEGER                   :: ir, ip, jp
    REAL( DP )                :: scale, r2, dist, arg, tmp, chargeanalytic
    REAL( DP )                :: r( 3 )
    REAL( DP ), ALLOCATABLE   :: hesslocal ( :, :, : )
    CHARACTER( LEN=80 )       :: sub_name = 'generate_hesserfc'
    REAL( DP ), EXTERNAL      :: qe_erfc
    !
    ! ... Aliases
    !
    TYPE( environ_cell ), POINTER :: cell
    !
    cell => hessian % cell
    !
    IF ( axis .LT. 1 .OR. axis .GT. 3 ) &
         & CALL env_errore(sub_name,'Wrong value of axis',1)
    !
    chargeanalytic = erfcvolume(dim,axis,width,spread,cell)
    scale = charge / chargeanalytic / sqrtpi / spread
    !
    ALLOCATE( hesslocal( 3, 3, cell%nnr ) )
    hesslocal = 0.D0
    !
    DO ir = 1, cell%ir_end
       !
       ! ... position in real space grid
       !
       CALL ir2r( cell, ir, r, physical )
       !
       ! ... do not include points outside the physical range
       !
       IF ( .NOT. physical ) CYCLE
       !
       ! ... displacement from origin
       !
       CALL displacement( dim, axis, pos, r, r )
       !
       ! ... minimum image convention
       !
       CALL minimum_image( cell, r, r2 )
       !
       ! ... compute hessian of error function
       !
       r = r * cell%alat
       dist = SQRT(r2) * cell%alat
       !
       arg = ( dist - width ) / spread
       !
       IF ( dist .GT. tol ) THEN
          DO ip = 1, 3
             DO jp = 1, 3
                tmp = - r(ip) * r(jp) * ( 1.D0 / dist + 2.D0 * arg / spread )
                IF ( ip .EQ. jp ) tmp = tmp + dist
                hesslocal( ip, jp, ir ) = - EXP( - arg**2 ) * tmp / dist**2
             ENDDO
          ENDDO
       END IF
       !
    END DO
    !
    hessian%of_r = hessian%of_r + hesslocal * scale
    DEALLOCATE( hesslocal )
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE generate_hesserfc
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE generate_axis( cell, icor, pos, axis )
!--------------------------------------------------------------------
    !
    TYPE( environ_cell ), INTENT(IN) :: cell
    INTEGER, INTENT(IN) :: icor
    REAL(DP), INTENT(IN) :: pos(3)
    REAL(DP), INTENT(OUT) :: axis( cell%nnr )
    !
    LOGICAL  :: physical
    INTEGER  :: ir
    REAL(DP) :: r(3), r2
    !
    DO ir = 1, cell%ir_end
       !
       ! ... position in real space grid
       !
       CALL ir2r( cell, ir, r, physical )
       !
       ! ... do not include points outside the physical range
       !
       IF ( .NOT. physical ) CYCLE
       !
       ! ... displacement from origin
       !
       r = r - pos
       !
       ! ... minimum image convention
       !
       CALL minimum_image( cell, r, r2 )
       !
       axis(ir) = r(icor)
       !
    END DO
    !
    axis = axis * cell%alat
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE generate_axis
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE generate_distance( cell, pos, distance )
!--------------------------------------------------------------------
    !
    TYPE(environ_cell), INTENT(IN) :: cell
    REAL(DP), INTENT(IN) :: pos(3)
    REAL(DP), INTENT(OUT) :: distance( 3, cell%nnr )
    !
    LOGICAL  :: physical
    INTEGER  :: ir
    REAL(DP) :: r(3), r2
    !
    DO ir = 1, cell%ir_end
       !
       ! ... position in real space grid
       !
       CALL ir2r( cell, ir, r, physical )
       !
       ! ... do not include points outside the physical range
       !
       IF ( .NOT. physical ) CYCLE
       !
       ! ... displacement from origin
       !
       r = r - pos
       !
       ! ... minimum image convention
       !
       CALL minimum_image( cell, r, r2 )
       !
       distance(:,ir) = r(:)
       !
    END DO
    !
    distance = distance * cell%alat
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE generate_distance
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  FUNCTION erfcvolume(dim,axis,width,spread,cell)
!--------------------------------------------------------------------
    !
    REAL(DP) :: erfcvolume
    !
    INTEGER, INTENT(IN) :: dim, axis
    REAL(DP), INTENT(IN) :: width, spread
    TYPE( environ_cell ), INTENT(IN) :: cell
    !
    REAL(DP) :: f1 = 0.0_DP , f2 = 0.0_DP
    REAL(DP) :: t, invt
    REAL( DP ), EXTERNAL      :: qe_erf
    !
    CHARACTER( LEN=80 ) :: fun_name = 'erfcvolume'
    !
    IF ( spread .LT. tol .OR. width .LT. tol ) &
         & CALL env_errore(fun_name,'Wrong parameters of erfc function',1)
    !
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
       erfcvolume = pi * width**2 * cell%at(axis,axis) * cell%alat * &
            ( ( 1.D0 + 0.5D0 * t**2 ) * f1  + t * f2 )
    CASE ( 2 )
       ! two-dimensional erfc, volume is exactly the one of the
       ! box, does not depend on spread
       erfcvolume = 2.D0 * width * cell%omega / cell%at(axis,axis) / cell%alat
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


