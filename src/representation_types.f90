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
!          Francesco Nattino  (THEOS and NCCR-MARVEL, EPFL)
!          Nicola Marzari     (THEOS and NCCR-MARVEL, EPFL)
!
!----------------------------------------------------------------------------------------
!>
!! This module contains the main control and parameter variables from QE Modules,
!! the definitions of Environ derived data types and the routines to handle the
!! basic derived data types (cell, density, gradient, hessian, electrons, system)
!!
!----------------------------------------------------------------------------------------
MODULE representation_types
    !------------------------------------------------------------------------------------
    !
    USE modules_constants, ONLY: DP
    !
    USE cell_types, ONLY: environ_cell
    !
    !------------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    TYPE environ_density
        !--------------------------------------------------------------------------------
        !
        LOGICAL :: update = .FALSE. ! optionally have an associated logical status
        !
        CHARACTER(LEN=80) :: label = ' '
        ! optionally have an associated label, used for printout and debugs
        !
        TYPE(environ_cell), POINTER :: cell => NULL()
        ! each quantity in real-space is associated with its definition domain
        !
        REAL(DP), ALLOCATABLE :: of_r(:)
        ! the quantity in real-space, local to each processor
        !
        !--------------------------------------------------------------------------------
        ! Multipole moments of the quantity
        !
        REAL(DP) :: charge
        REAL(DP) :: dipole(3)
        REAL(DP) :: quadrupole(3)
        !
        !--------------------------------------------------------------------------------
    END TYPE environ_density
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    TYPE environ_gradient
        !--------------------------------------------------------------------------------
        !
        LOGICAL :: update = .FALSE. ! optionally have an associated logical status
        !
        CHARACTER(LEN=80) :: label = ' '
        ! optionally have an associated label, used for printout and debugs
        !
        TYPE(environ_cell), POINTER :: cell => NULL()
        ! each quantity in real-space is associated with its definition domain
        !
        REAL(DP), ALLOCATABLE :: of_r(:, :)
        ! the quantity in real-space, local to each processor
        !
        TYPE(environ_density) :: modulus
        !
        !--------------------------------------------------------------------------------
    END TYPE environ_gradient
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    TYPE environ_hessian
        !--------------------------------------------------------------------------------
        !
        LOGICAL :: update = .FALSE. ! optionally have an associated logical status
        !
        CHARACTER(LEN=80) :: label = ' '
        ! optionally have an associated label, used for printout and debugs
        !
        TYPE(environ_cell), POINTER :: cell => NULL()
        ! each quantity in real-space is associated with its definition domain
        !
        REAL(DP), ALLOCATABLE :: of_r(:, :, :)
        ! the quantity in real-space, local to each processor
        !
        TYPE(environ_density) :: laplacian
        !
        !--------------------------------------------------------------------------------
    END TYPE environ_hessian
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    TYPE environ_functions
        !--------------------------------------------------------------------------------
        !
        INTEGER :: type_
        INTEGER :: axis, dim
        REAL(DP) :: width, spread, volume
        REAL(DP), POINTER :: pos(:)
        ! environ_functions are not designed to be mobile, thus position
        ! can be included in the definition of the type
        !
        !--------------------------------------------------------------------------------
    END TYPE environ_functions
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
    ! Keeping imports private
    !
    PRIVATE :: DP, environ_cell
    !
    !------------------------------------------------------------------------------------
END MODULE representation_types
!----------------------------------------------------------------------------------------
