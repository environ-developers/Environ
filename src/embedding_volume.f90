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
!> Module to compute an enthalpy functiona, defined as the quantum volume
!! of the system times the external pressure of the environment.
!! Original method developed in Cococcioni et al, PRL (2005)
!
! Authors: Oliviero Andreussi (Department of Physics, UNT)
!          Nicola Marzari     (THEOS and NCCR-MARVEL, EPFL)
!
! Original version by M. Cococcioni and N. Marzari (MIT)
!
!----------------------------------------------------------------------------
MODULE embedding_volume
!----------------------------------------------------------------------------
  !
  USE environ_types
  USE environ_output
  USE modules_constants,  ONLY : e2
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC :: calc_devolume_dboundary, calc_evolume
  !
CONTAINS
!  Subroutine: calc_devolume_dboundary
!
!> Calculates the PV contribution to the potential
!--------------------------------------------------------------------
  SUBROUTINE calc_devolume_dboundary( pressure, boundary, de_dboundary )
!--------------------------------------------------------------------
    IMPLICIT NONE
    !
    ! ... Declares variables
    !
    REAL( DP ), INTENT(IN) :: pressure
    TYPE( environ_boundary ), INTENT(IN) :: boundary
    TYPE( environ_density ), INTENT(INOUT) :: de_dboundary
    !
    CHARACTER( LEN=80 ) :: sub_name = 'calc_devolume_dboundary'
    !
    ! ... The functional derivative of the volume term is just unity
    !
    de_dboundary%of_r = de_dboundary%of_r + pressure
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE calc_devolume_dboundary
!--------------------------------------------------------------------
!  Subroutine: calc_evolume
!
!> Calculates the PV contribution to the energy
!--------------------------------------------------------------------
  SUBROUTINE calc_evolume( pressure, boundary, evolume )
!--------------------------------------------------------------------
    IMPLICIT NONE
    !
    ! ... Declares variables
    !
    REAL( DP ), INTENT(IN) :: pressure
    TYPE( environ_boundary ), TARGET, INTENT(IN) :: boundary
    REAL( DP ), INTENT(OUT) :: evolume
    !
    ! ... Computes the PV energy
    !
    evolume = pressure * boundary%volume * e2 / 2.D0
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE calc_evolume
!--------------------------------------------------------------------
!--------------------------------------------------------------------
END MODULE embedding_volume
!--------------------------------------------------------------------
