!----------------------------------------------------------------------------------------
!
! Copyright (C) 2021 ENVIRON (www.quantum-environ.org)
!
!----------------------------------------------------------------------------------------
!
!     This file is part of Environ version 2.1
!
!     Environ 2.1 is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation, either version 2 of the License, or
!     (at your option) any later version.
!
!     Environ 2.1 is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more detail, either the file
!     `License' in the root directory of the present distribution, or
!     online at <http://www.gnu.org/licenses/>.
!
!----------------------------------------------------------------------------------------
!
! Authors: Edan Bainglass (Department of Physics, UNT)
!
!----------------------------------------------------------------------------------------
!>
!!
!----------------------------------------------------------------------------------------
MODULE prog_utils
    !------------------------------------------------------------------------------------
    !
    USE class_io, ONLY: io
    !
    USE environ_param, ONLY: DP
    !
    USE environ_api, ONLY: environ_interface
    !
    USE cmdline_args
    !
    USE parsers
    !
    !------------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    PRIVATE
    !
    PUBLIC :: init_environ_from_cube
    !
    !------------------------------------------------------------------------------------
CONTAINS
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !                                  READING ROUTINES
    !
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_environ_from_cube(environ, nelec, rho)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_interface), INTENT(INOUT) :: environ
        REAL(DP), INTENT(OUT) :: nelec
        REAL(DP), ALLOCATABLE, OPTIONAL, INTENT(OUT) :: rho(:)
        !
        INTEGER :: nat
        INTEGER :: ntyp
        INTEGER, ALLOCATABLE :: ityp(:)
        CHARACTER(LEN=2), ALLOCATABLE :: label(:)
        REAL(DP), ALLOCATABLE :: zv(:)
        REAL(DP), ALLOCATABLE :: tau(:, :)
        REAL(DP) :: origin(3)
        INTEGER :: nr(3)
        REAL(DP) :: at(3, 3)
        !
        CHARACTER(LEN=80) :: sub_name = 'init_environ_from_cube'
        !
        !--------------------------------------------------------------------------------
        !
        CALL read_cube(nat, ntyp, ityp, label, zv, nelec, tau, origin, nr, at, rho)
        !
        !--------------------------------------------------------------------------------
        ! Initialize Environ
        !
        CALL environ%read_input(inputfile, SIZE(label))
        !
        CALL environ%setup%init(use_pbc_corr)
        !
        IF (ANY(nr == 1)) THEN
            CALL environ%setup%init_cell(io%comm, at)
        ELSE
            CALL environ%setup%init_cell(io%comm, at, nr=nr)
        END IF
        !
        CALL environ%setup%init_cores()
        !
        CALL environ%main%init(nat, ntyp, label, ityp, zv)
        !
        CALL environ%main%update_ions(nat, tau, origin)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_environ_from_cube
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE prog_utils
!----------------------------------------------------------------------------------------
