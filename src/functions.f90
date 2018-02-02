! Copyright (C) 2018 ENVIRON (www.quantum-environment.org)
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
! Module containing the main routines to handle
!
!              environ_functions
!
! derived data types.
!
! Environ_functions contains all the details of analytic functions needed
! by Environ modules and defined on the three-dimensional real-space
! domain, together with the routines to handle the derived data type and
! to generate the functions from their parameters.
!
!----------------------------------------------------------------------------
!  TYPE environ_functions
!----------------------------------------------------------------------------
!
!     INTEGER :: type
!     INTEGER :: axis, dim
!     REAL( DP ) :: width, spread, volume
!     REAL( DP ), DIMENSION(:), POINTER :: pos
!     ! environ_functions are not designed to be mobile,
!     ! thus position can be included in the definition
!     ! of the type
!
!----------------------------------------------------------------------------
!  END TYPE environ_functions
!----------------------------------------------------------------------------
!
! Authors: Oliviero Andreussi (Department of Physics, UNT)
!
!----------------------------------------------------------------------------
MODULE functions
  !
  USE environ_types
  USE environ_generate_function
  !
  PRIVATE
  !
  PUBLIC :: create_environ_functions, destroy_environ_functions, copy_environ_functions, &
       & density_of_functions, gradient_of_functions, laplacian_of_functions, &
       & hessian_of_functions
  !
  INTERFACE create_environ_functions
     MODULE PROCEDURE create_environ_functions_scalar, create_environ_functions_array
  END INTERFACE create_environ_functions
  !
  INTERFACE density_of_functions
     MODULE PROCEDURE density_of_functions_scalar, density_of_functions_array
  END INTERFACE density_of_functions
  !
  INTERFACE gradient_of_functions
     MODULE PROCEDURE gradient_of_functions_scalar, gradient_of_functions_array
  END INTERFACE gradient_of_functions
  !
  INTERFACE laplacian_of_functions
     MODULE PROCEDURE laplacian_of_functions_scalar, laplacian_of_functions_array
  END INTERFACE laplacian_of_functions
  !
  INTERFACE hessian_of_functions
     MODULE PROCEDURE hessian_of_functions_scalar, hessian_of_functions_array
  END INTERFACE hessian_of_functions
  !
CONTAINS
!--------------------------------------------------------------------
  SUBROUTINE create_environ_functions_scalar( type, dimm, axis, pos, width, spreadd, volume, f )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: type
    INTEGER, INTENT(IN) :: dimm, axis
    REAL( DP ), INTENT(IN) :: width, spreadd, volume
    REAL( DP ), DIMENSION(3), TARGET, INTENT(IN) :: pos
    TYPE( environ_functions ), INTENT(INOUT) :: f
    !
    INTEGER :: i
    !
    f%type   = type
    f%dim    = dimm
    f%axis   = axis
    f%spread = spreadd
    f%width  = width
    f%volume = volume
    !
    f%pos => pos
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE create_environ_functions_scalar
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE create_environ_functions_array( n, type, dimm, axis, pos, width, spreadd, volume, f )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: n
    INTEGER, INTENT(IN) :: type
    INTEGER, DIMENSION(n), INTENT(IN) :: dimm, axis
    REAL( DP ), DIMENSION(n), INTENT(IN) :: width, spreadd, volume
    REAL( DP ), DIMENSION(3,n), TARGET, INTENT(IN) :: pos
    TYPE( environ_functions ), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: f
    !
    INTEGER :: i
    !
    ALLOCATE( f(n) )
    DO i = 1, n
       CALL create_environ_functions_scalar( type, dimm(i), axis(i), pos(:,i), width(i), spreadd(i), volume(i), f(i) )
    ENDDO
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE create_environ_functions_array
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE copy_environ_functions( foriginal, fcopy )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    TYPE( environ_functions ), INTENT(IN) :: foriginal
    TYPE( environ_functions ), INTENT(OUT) :: fcopy
    !
    fcopy % pos   => foriginal % pos
    !
    fcopy % type   = foriginal % type
    fcopy % dim    = foriginal % dim
    fcopy % axis   = foriginal % axis
    fcopy % spread = foriginal % spread
    fcopy % width  = foriginal % width
    fcopy % volume = foriginal % volume
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE copy_environ_functions
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE destroy_environ_functions( n, f )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: n
    TYPE( environ_functions ), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: f
    !
    CHARACTER( LEN=80 ) :: sub_name = 'destroy_environ_functions'
    !
    INTEGER :: i
    !
    IF ( .NOT. ALLOCATED( f ) ) &
         & CALL errore(sub_name,'Trying to destroy a non allocated object',1)
    IF ( SIZE(f) .NE. n ) &
         & CALL errore(sub_name,'Inconsistent size of allocated object',1)
    !
    DO i = 1, n
       IF ( .NOT. ASSOCIATED( f(i)%pos ) ) &
            & CALL errore(sub_name,'Trying to destroy a non allocated object',1)
       NULLIFY( f(i)%pos )
    ENDDO
    !
    DEALLOCATE( f )
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE destroy_environ_functions
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE density_of_functions_scalar( functions, density, zero )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    TYPE( environ_functions ), TARGET, INTENT(IN) :: functions
    TYPE( environ_density ), TARGET, INTENT(INOUT) :: density
    LOGICAL, INTENT(IN), OPTIONAL :: zero
    !
    INTEGER :: i
    REAL( DP ) :: local_charge
    CHARACTER( LEN=80 ) :: sub_name = 'density_of_functions'
    !
    INTEGER, POINTER :: type, dim, axis, nnr
    REAL( DP ), POINTER :: alat, omega, at(:,:)
    REAL( DP ), POINTER :: charge, spread, width
    REAL( DP ), DIMENSION(:), POINTER :: pos
    !
    IF ( PRESENT( zero ) .AND. zero ) density%of_r = 0.D0
    !
    alat => density%cell%alat
    omega => density%cell%omega
    at => density%cell%at
    !
    type   => functions%type
    pos    => functions%pos
    spread => functions%spread
    charge => functions%volume
    width  => functions%width
    dim    => functions%dim
    axis   => functions%axis
    nnr    => density%cell%nnr
    SELECT CASE ( type )
    CASE ( 1 ) ! Gaussian
       CALL generate_gaussian(nnr, dim, axis, charge, spread, pos, density%of_r)
    CASE ( 2 ) ! CHARGE * NORMALIZED_ERFC_HALF(X) ! integrates to charge
       CALL generate_erfc(nnr, dim, axis, charge, width, spread, pos, density%of_r)
    CASE ( 3 ) ! Exponential
       CALL generate_exponential(nnr, spread, pos, density%of_r)
    CASE ( 4 ) ! CHARGE * NORMALIZED_ERFC_HALF(X) * VOLUME_NORMALIZED_ERFC_HALF ! goes from charge to 0
       local_charge = erfcvolume(dim,axis,width,spread,alat,omega,at) * charge
       CALL generate_erfc(nnr, dim, axis, local_charge, width, spread, pos, density%of_r)
    CASE ( 5 ) ! CHARGE * ( 1 - NORMALIZED_ERFC_HALF(x) * VOLUME_NORMALIZED_ERFC_HALF ) ! goes from 0 to charge
       local_charge = - erfcvolume(dim,axis,width,spread,alat,omega,at) * charge
       CALL generate_erfc(nnr, dim, axis, local_charge, width, spread, pos, density%of_r)
       density % of_r = density % of_r + charge
    END SELECT
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE density_of_functions_scalar
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE density_of_functions_array( nfunctions, functions, density, zero )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: nfunctions
    TYPE( environ_functions ), DIMENSION(nfunctions), TARGET, INTENT(IN) :: functions
    TYPE( environ_density ), TARGET, INTENT(INOUT) :: density
    LOGICAL, INTENT(IN), OPTIONAL :: zero
    !
    INTEGER :: i
    !
    IF ( PRESENT( zero ) .AND. zero ) density%of_r = 0.D0
    !
    DO i = 1, nfunctions
       CALL density_of_functions_scalar(functions(i),density)
    END DO
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE density_of_functions_array
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE gradient_of_functions_scalar( functions, gradient, zero )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    TYPE( environ_functions ), TARGET, INTENT(IN) :: functions
    TYPE( environ_gradient ), TARGET, INTENT(INOUT) :: gradient
    LOGICAL, INTENT(IN), OPTIONAL :: zero
    !
    INTEGER :: i
    REAL( DP ) :: local_charge
    CHARACTER( LEN=80 ) :: sub_name = 'gradient_of_functions'
    !
    INTEGER, POINTER :: type, dim, axis, nnr
    REAL( DP ), POINTER :: alat, omega, at(:,:)
    REAL( DP ), POINTER :: charge, spread, width
    REAL( DP ), DIMENSION(:), POINTER :: pos
    !
    alat => gradient%cell%alat
    omega => gradient%cell%omega
    at => gradient%cell%at
    !
    IF ( PRESENT(zero) .AND. zero ) gradient%of_r = 0.D0
    !
    type   => functions%type
    pos    => functions%pos
    spread => functions%spread
    charge => functions%volume
    width  => functions%width
    dim    => functions%dim
    axis   => functions%axis
    nnr    => gradient%cell%nnr
    SELECT CASE ( type )
    CASE ( 1 ) ! Gaussian
       CALL generate_gradgaussian(nnr, dim, axis, charge, spread, pos, gradient%of_r)
    CASE ( 2 ) ! CHARGE * NORMALIZED_ERFC_HALF(X) ! integrates to charge
       CALL generate_graderfc(nnr, dim, axis, charge, width, spread, pos, gradient%of_r)
    CASE ( 3 ) ! Exponential
       CALL generate_gradexponential(nnr, spread, pos, gradient%of_r)
    CASE ( 4 ) ! CHARGE * NORMALIZED_ERFC_HALF(X) * VOLUME_NORMALIZED_ERFC_HALF ! goes from charge to 0
       local_charge = erfcvolume(dim,axis,width,spread,alat,omega,at) * charge
       CALL generate_graderfc(nnr, dim, axis, local_charge, width, spread, pos, gradient%of_r)
    CASE ( 5 ) ! CHARGE * ( 1 - NORMALIZED_ERFC_HALF(x) * VOLUME_NORMALIZED_ERFC_HALF ) ! goes from 0 to charge
       local_charge = - erfcvolume(dim,axis,width,spread,alat,omega,at) * charge
       CALL generate_graderfc(nnr, dim, axis, local_charge, width, spread, pos, gradient%of_r)
    END SELECT
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE gradient_of_functions_scalar
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE gradient_of_functions_array( nfunctions, functions, gradient, zero )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: nfunctions
    TYPE( environ_functions ), DIMENSION(nfunctions), TARGET, INTENT(IN) :: functions
    TYPE( environ_gradient ), TARGET, INTENT(INOUT) :: gradient
    LOGICAL, INTENT(IN), OPTIONAL :: zero
    !
    INTEGER :: i
    !
    IF ( PRESENT(zero) .AND. zero ) gradient%of_r = 0.D0
    !
    DO i = 1, nfunctions
       CALL gradient_of_functions_scalar(functions(i),gradient,.FALSE.)
    END DO
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE gradient_of_functions_array
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE laplacian_of_functions_scalar( functions, laplacian, zero )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    TYPE( environ_functions ), TARGET, INTENT(IN) :: functions
    TYPE( environ_density ), TARGET, INTENT(INOUT) :: laplacian
    LOGICAL, INTENT(IN), OPTIONAL :: zero
    !
    INTEGER :: i
    REAL( DP ) :: local_charge
    CHARACTER( LEN=80 ) :: sub_name = 'laplacian_of_functions'
    !
    INTEGER, POINTER :: type, dim, axis, nnr
    REAL( DP ), POINTER :: alat, omega, at(:,:)
    REAL( DP ), POINTER :: charge, spread, width
    REAL( DP ), DIMENSION(:), POINTER :: pos
    !
    alat => laplacian%cell%alat
    omega => laplacian%cell%omega
    at => laplacian%cell%at
    !
    IF ( PRESENT(zero) .AND. zero ) laplacian%of_r = 0.D0
    !
    type   => functions%type
    pos    => functions%pos
    spread => functions%spread
    charge => functions%volume
    width  => functions%width
    dim    => functions%dim
    axis   => functions%axis
    nnr    => laplacian%cell%nnr
    SELECT CASE ( type )
    CASE ( 1 ) ! Gaussian
       CALL errore(sub_name,'Options not yet implemented',1)
    CASE ( 2 ) ! CHARGE * NORMALIZED_ERFC_HALF(X) ! integrates to charge
       CALL generate_laplerfc(nnr, dim, axis, charge, width, spread, pos, laplacian%of_r)
    CASE ( 3 ) ! Exponential
       CALL errore(sub_name,'Options not yet implemented',1)
    CASE ( 4 ) ! CHARGE * NORMALIZED_ERFC_HALF(X) * VOLUME_NORMALIZED_ERFC_HALF ! goes from charge to 0
       local_charge = erfcvolume(dim,axis,width,spread,alat,omega,at) * charge
       CALL generate_laplerfc(nnr, dim, axis, local_charge, width, spread, pos, laplacian%of_r)
    CASE ( 5 ) ! CHARGE * ( 1 - NORMALIZED_ERFC_HALF(x) * VOLUME_NORMALIZED_ERFC_HALF ) ! goes from 0 to charge
       local_charge = - erfcvolume(dim,axis,width,spread,alat,omega,at) * charge
       CALL generate_laplerfc(nnr, dim, axis, local_charge, width, spread, pos, laplacian%of_r)
    END SELECT
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE laplacian_of_functions_scalar
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE laplacian_of_functions_array( nfunctions, functions, laplacian, zero )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: nfunctions
    TYPE( environ_functions ), DIMENSION(nfunctions), TARGET, INTENT(IN) :: functions
    TYPE( environ_density ), TARGET, INTENT(INOUT) :: laplacian
    LOGICAL, INTENT(IN), OPTIONAL :: zero
    !
    INTEGER :: i
    !
    IF ( PRESENT(zero) .AND. zero ) laplacian%of_r = 0.D0
    !
    DO i = 1, nfunctions
       CALL laplacian_of_functions_scalar( functions(i), laplacian, .FALSE. )
    END DO
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE laplacian_of_functions_array
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE hessian_of_functions_scalar( functions, hessian, zero )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    TYPE( environ_functions ), TARGET, INTENT(IN) :: functions
    TYPE( environ_hessian ), TARGET, INTENT(INOUT) :: hessian
    LOGICAL, INTENT(IN), OPTIONAL :: zero
    !
    INTEGER :: i
    REAL( DP ) :: local_charge
    CHARACTER( LEN=80 ) :: sub_name = 'hessian_of_functions'
    !
    INTEGER, POINTER :: type, dim, axis, nnr
    REAL( DP ), POINTER :: alat, omega, at(:,:)
    REAL( DP ), POINTER :: charge, spread, width
    REAL( DP ), DIMENSION(:), POINTER :: pos
    !
    alat => hessian%cell%alat
    omega => hessian%cell%omega
    at => hessian%cell%at
    !
    IF ( PRESENT(zero) .AND. zero ) hessian%of_r = 0.D0
    !
    type   => functions%type
    pos    => functions%pos
    spread => functions%spread
    charge => functions%volume
    width  => functions%width
    dim    => functions%dim
    axis   => functions%axis
    nnr    => hessian%cell%nnr
    SELECT CASE ( type )
    CASE ( 1 ) ! Gaussian
       CALL errore(sub_name,'Options not yet implemented',1)
    CASE ( 2 ) ! CHARGE * NORMALIZED_ERFC_HALF(X) ! integrates to charge
       CALL generate_hesserfc(nnr, dim, axis, charge, width, spread, pos, hessian%of_r)
    CASE ( 3 ) ! Exponential
       CALL errore(sub_name,'Options not yet implemented',1)
    CASE ( 4 ) ! CHARGE * NORMALIZED_ERFC_HALF(X) * VOLUME_NORMALIZED_ERFC_HALF ! goes from charge to 0
       local_charge = erfcvolume(dim,axis,width,spread,alat,omega,at) * charge
       CALL generate_hesserfc(nnr, dim, axis, local_charge, width, spread, pos, hessian%of_r)
    CASE ( 5 ) ! CHARGE * ( 1 - NORMALIZED_ERFC_HALF(x) * VOLUME_NORMALIZED_ERFC_HALF ) ! goes from 0 to charge
       local_charge = - erfcvolume(dim,axis,width,spread,alat,omega,at) * charge
       CALL generate_hesserfc(nnr, dim, axis, local_charge, width, spread, pos, hessian%of_r)
    END SELECT
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE hessian_of_functions_scalar
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE hessian_of_functions_array( nfunctions, functions, hessian, zero )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: nfunctions
    TYPE( environ_functions ), DIMENSION(nfunctions), TARGET, INTENT(IN) :: functions
    TYPE( environ_hessian ), TARGET, INTENT(INOUT) :: hessian
    LOGICAL, INTENT(IN), OPTIONAL :: zero
    !
    INTEGER :: i
    !
    IF ( PRESENT(zero) .AND. zero ) hessian%of_r = 0.D0
    !
    DO i = 1, nfunctions
       CALL hessian_of_functions_scalar( functions(i), hessian, .FALSE. )
    END DO
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE hessian_of_functions_array
!--------------------------------------------------------------------
!----------------------------------------------------------------------------
END MODULE functions
!----------------------------------------------------------------------------
