! Copyright (C) 2018 ENVIRON (www.quantum-environment.org)
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
!> Module to compute an confinement functionas, defined on a smooth
!> continuum interface function
!
! Authors: Oliviero Andreussi (Department of Physics, UNT)
!
!----------------------------------------------------------------------------
MODULE embedding_confine
!----------------------------------------------------------------------------
  !
  USE environ_types
  USE environ_output
  USE environ_base,  ONLY : e2
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC :: calc_deconfine_dboundary, calc_vconfine
  !
CONTAINS
!  Subroutine: calc_vconfine
!
!> Calculates the confine contribution to the potential
!--------------------------------------------------------------------
  SUBROUTINE calc_vconfine( confine, boundary, vconfine )
!--------------------------------------------------------------------
    IMPLICIT NONE
    !
    ! ... Declares variables
    !
    REAL( DP ), INTENT(IN) :: confine
    TYPE( environ_boundary ), INTENT(IN) :: boundary
    TYPE( environ_density ), INTENT(INOUT) :: vconfine
    !
    CHARACTER( LEN=80 ) :: sub_name = 'calc_vconfine'
    !
    ! ... The confine potetial is defined as confine * ( 1 - s(r) )
    !
    vconfine%of_r = 0.D0
    vconfine%of_r = confine * ( 1.D0 - boundary % scaled % of_r )
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE calc_vconfine
!--------------------------------------------------------------------
!  Subroutine: calc_deconfine_dboundary
!
!> Calculates the confine contribution to the interface potential
!--------------------------------------------------------------------
  SUBROUTINE calc_deconfine_dboundary( confine, rhoelec, de_dboundary )
!--------------------------------------------------------------------
    IMPLICIT NONE
    !
    ! ... Declares variables
    !
    REAL( DP ), INTENT(IN) :: confine
    TYPE( environ_density ), INTENT(IN) :: rhoelec
    TYPE( environ_density ), INTENT(INOUT) :: de_dboundary
    !
    CHARACTER( LEN=80 ) :: sub_name = 'calc_deconfine_dboundary'
    !
    ! ... The functional derivative of the confine term is
    !     - confine * rho^elec(r)
    !
    de_dboundary%of_r = de_dboundary%of_r - confine * rhoelec % of_r
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE calc_deconfine_dboundary
!--------------------------------------------------------------------
!--------------------------------------------------------------------
END MODULE embedding_confine
!--------------------------------------------------------------------
