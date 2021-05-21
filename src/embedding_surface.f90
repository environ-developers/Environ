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
!----------------------------------------------------------------------------------------
!
! Authors: Oliviero Andreussi (Department of Physics, UNT)
!          Nicola Marzari     (THEOS and NCCR-MARVEL, EPFL)
!
! Original version by D. Scherlis, M. Cococcioni, and N. Marzari (MIT)
!
!----------------------------------------------------------------------------------------
!>
!! Module to compute a cavitation functional, defined as the quantum surface
!! of the system times the surface tension of the environment.
!! Original method developed in Scherlis et al, J. Chem. Phys. (2006),
!! but formulas for the quantum surface were not correct and have been
!! rederived by Andreussi et al., J. Chem. Phys. 136, 064102 (2012)
!!
!! The variables needed to introduce an environment with a surface tension
!! around the system.
!!
!----------------------------------------------------------------------------------------
MODULE embedding_surface
    !------------------------------------------------------------------------------------
    !
    USE modules_constants, ONLY: DP, e2
    !
    USE physical_types, ONLY: environ_boundary
    USE representation_types, ONLY: environ_density
    !
    !------------------------------------------------------------------------------------
    !
    PRIVATE
    !
    PUBLIC :: calc_desurface_dboundary, calc_esurface
    !
    !------------------------------------------------------------------------------------
CONTAINS
    !------------------------------------------------------------------------------------
    !>
    !! Calculates the cavitation contribution to the potential
    !
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_desurface_dboundary(surface_tension, boundary, de_dboundary)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP), INTENT(IN) :: surface_tension
        TYPE(environ_boundary), TARGET, INTENT(IN) :: boundary
        !
        TYPE(environ_density), TARGET, INTENT(INOUT) :: de_dboundary
        !
        CHARACTER(LEN=80) :: sub_name = 'calc_desurface_dboundary'
        !
        !--------------------------------------------------------------------------------
        !
        de_dboundary%of_r = de_dboundary%of_r + surface_tension * boundary%dsurface%of_r
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_desurface_dboundary
    !------------------------------------------------------------------------------------
    !>
    !! Calculates the cavitation contribution to the energy
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_esurface(surface_tension, boundary, esurface)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP), INTENT(IN) :: surface_tension
        TYPE(environ_boundary), TARGET, INTENT(IN) :: boundary
        !
        REAL(DP), INTENT(OUT) :: esurface
        !
        CHARACTER(LEN=80) :: sub_name = 'calc_esurface'
        !
        !--------------------------------------------------------------------------------
        !
        esurface = surface_tension * boundary%surface * e2 / 2.D0
        ! computes the cavitation energy
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_esurface
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE embedding_surface
!----------------------------------------------------------------------------------------
