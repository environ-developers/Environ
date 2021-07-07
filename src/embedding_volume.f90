!----------------------------------------------------------------------------------------
!
! Copyright (C) 2018-2021 ENVIRON (www.quantum-environ.org)
!
!----------------------------------------------------------------------------------------
!
!     This file is part of Environ version 2.0
!     
!     Environ 2.0 is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation, either version 2 of the License, or
!     (at your option) any later version.
!     
!     Environ 2.0 is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more detail, either the file
!     `License' in the root directory of the present distribution, or
!     online at <http://www.gnu.org/licenses/>.
!
!----------------------------------------------------------------------------------------
!
! Authors: Oliviero Andreussi (Department of Physics, UNT)
!          Nicola Marzari     (THEOS and NCCR-MARVEL, EPFL)
!          Edan Bainglass     (Department of Physics, UNT)
!
! Original version by M. Cococcioni and N. Marzari (MIT)
!
!----------------------------------------------------------------------------------------
!>
!! Module to compute an enthalpy functiona, defined as the quantum volume
!! of the system times the external pressure of the environment.
!! Original method developed in Cococcioni et al, PRL (2005)
!!
!----------------------------------------------------------------------------------------
MODULE embedding_volume
    !------------------------------------------------------------------------------------
    !
    USE environ_param, ONLY: DP, e2
    !
    USE types_physical, ONLY: environ_boundary
    USE types_representation, ONLY: environ_density
    !
    !------------------------------------------------------------------------------------
    !
    PRIVATE
    !
    PUBLIC :: calc_devolume_dboundary, calc_evolume
    !
    !------------------------------------------------------------------------------------
CONTAINS
    !------------------------------------------------------------------------------------
    !>
    !! Calculates the PV contribution to the potential
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_devolume_dboundary(pressure, boundary, de_dboundary)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP), INTENT(IN) :: pressure
        TYPE(environ_boundary), INTENT(IN) :: boundary
        !
        TYPE(environ_density), INTENT(INOUT) :: de_dboundary
        !
        CHARACTER(LEN=80) :: sub_name = 'calc_devolume_dboundary'
        !
        !--------------------------------------------------------------------------------
        !
        de_dboundary%of_r = de_dboundary%of_r + pressure
        ! the functional derivative of the volume term is just unity
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_devolume_dboundary
    !------------------------------------------------------------------------------------
    !>
    !! Calculates the PV contribution to the energy
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_evolume(pressure, boundary, evolume)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP), INTENT(IN) :: pressure
        TYPE(environ_boundary), TARGET, INTENT(IN) :: boundary
        !
        REAL(DP), INTENT(OUT) :: evolume
        !
        !--------------------------------------------------------------------------------
        !
        evolume = pressure * boundary%volume * e2 / 2.D0 ! computes the PV energy
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_evolume
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE embedding_volume
!----------------------------------------------------------------------------------------
