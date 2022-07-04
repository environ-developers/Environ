!----------------------------------------------------------------------------------------
!
! Copyright (C) 2018-2022 ENVIRON (www.quantum-environ.org)
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
    SUBROUTINE init_environ_from_cube(environ, rho, reduce_cell)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        LOGICAL, OPTIONAL, INTENT(IN) :: reduce_cell
        !
        TYPE(environ_interface), INTENT(INOUT) :: environ
        REAL(DP), ALLOCATABLE, OPTIONAL, INTENT(OUT) :: rho(:)
        !
        INTEGER :: nat
        INTEGER :: ntyp
        INTEGER, ALLOCATABLE :: ityp(:)
        INTEGER, ALLOCATABLE :: species(:)
        REAL(DP), ALLOCATABLE :: zv(:)
        REAL(DP), ALLOCATABLE :: tau(:, :)
        REAL(DP) :: origin(3)
        INTEGER :: nr(3)
        REAL(DP) :: at(3, 3)
        !
        CHARACTER(LEN=80) :: routine = 'init_environ_from_cube'
        !
        !--------------------------------------------------------------------------------
        !
        CALL read_cube(nat, ntyp, ityp, species, zv, tau, origin, nr, at, rho)
        !
        !--------------------------------------------------------------------------------
        ! Initialize Environ
        !
        CALL environ%read_input(inputfile, ntyp)
        !
        CALL environ%setup%init()
        !
        IF (PRESENT(reduce_cell)) THEN
            IF (reduce_cell) CALL get_reduced_cell(nat, at, tau)
        END IF
        !
        IF (ANY(ABS(nr) == 1)) THEN
            CALL environ%setup%init_cell(io%comm, at)
        ELSE
            CALL environ%setup%init_cell(io%comm, at, nr=nr)
        END IF
        !
        CALL environ%setup%init_numerical(use_internal_pbc_corr)
        !
        CALL environ%main%init(nat, ntyp, ityp, zv, number=species)
        !
        CALL environ%main%update_ions(nat, tau, origin)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_environ_from_cube
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !                              PRIVATE HELPER ROUTINES
    !
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE get_reduced_cell(nat, at, tau)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: nat
        !
        REAL(DP), INTENT(INOUT) :: at(3, 3)
        REAL(DP), INTENT(INOUT) :: tau(:, :)
        !
        INTEGER :: i
        !
        REAL(DP), PARAMETER :: fluff = 7.D0
        REAL(DP) :: min_vec(3), max_vec(3), shift(3)
        !
        !--------------------------------------------------------------------------------
        !
        at = 0.D0
        !
        min_vec = MINVAL(tau, DIM=2) - fluff
        max_vec = MAXVAL(tau, DIM=2) + fluff
        !
        shift = max_vec - min_vec
        !
        DO i = 1, nat
            tau(:, i) = tau(:, i) - min_vec
        END DO
        !
        DO i = 1, 3
            at(i, i) = shift(i)
        END DO
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE get_reduced_cell
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE prog_utils
!----------------------------------------------------------------------------------------
