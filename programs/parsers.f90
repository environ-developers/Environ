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
MODULE parsers
    !------------------------------------------------------------------------------------
    !
    USE class_io, ONLY: io
    USE env_mp, ONLY: env_mp_bcast, env_mp_abort
    USE env_array_ops, ONLY: env_get_index
    !
    USE environ_param, ONLY: DP
    !
    USE environ_api, ONLY: get_atom_labels
    !
    USE cmdline_args
    !
    !------------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    PRIVATE
    !
    PUBLIC :: read_cube, char2real
    !
    !------------------------------------------------------------------------------------
CONTAINS
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE read_cube(nat, ntyp, ityp, atom_label, zv, tau, origin, nr, at, rho)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(OUT) :: nat, ntyp
        INTEGER, ALLOCATABLE, INTENT(OUT) :: ityp(:)
        CHARACTER(LEN=2), ALLOCATABLE, INTENT(OUT) :: atom_label(:)
        REAL(DP), INTENT(OUT) :: origin(3)
        INTEGER, INTENT(OUT) :: nr(3)
        REAL(DP), INTENT(OUT) :: at(3, 3)
        REAL(DP), ALLOCATABLE, INTENT(OUT) :: zv(:)
        REAL(DP), ALLOCATABLE, INTENT(OUT) :: tau(:, :)
        REAL(DP), ALLOCATABLE, OPTIONAL, INTENT(OUT) :: rho(:)
        !
        INTEGER :: ios
        LOGICAL :: ext
        !
        INTEGER :: i, j, k, l, count, unit, current, nnt
        !
        REAL(DP) :: fact ! angstrom to bohr conversion factor
        !
        INTEGER, ALLOCATABLE :: atomic_number(:)
        INTEGER, ALLOCATABLE :: charge_index(:) ! index of unique charges
        REAL(DP), ALLOCATABLE :: charge(:) ! unprocessed charges
        INTEGER, ALLOCATABLE :: species(:) ! register for checked atomic numbers
        REAL(DP), ALLOCATABLE :: unsorted(:) ! unsorted density read from file
        !
        CHARACTER(LEN=14) :: current_combo ! current number/charge
        CHARACTER(LEN=14), ALLOCATABLE :: combo(:) ! unique set of number/charge
        !
        CHARACTER(LEN=80) :: sub_name = 'read_cube'
        !
        !--------------------------------------------------------------------------------
        ! Open cube file
        !
        unit = io%find_free_unit()
        INQUIRE (file=TRIM(cubefile), exist=ext)
        !
        IF (.NOT. ext) CALL io%error(sub_name, TRIM(cubefile)//' not found', 1)
        !
        OPEN (unit=unit, file=TRIM(cubefile), status='old')
        !
        !--------------------------------------------------------------------------------
        !
        CALL io%writer('') ! blank line
        !
        CALL io%writer('Reading cube data from '//TRIM(cubefile))
        !
        !--------------------------------------------------------------------------------
        ! Discard comment lines
        !
        IF (io%lnode) THEN
            !
            DO i = 1, 2
                READ (unit, *, iostat=ios)
                !
                CALL check_ios(ios)
                !
            END DO
            !
        END IF
        !
        !--------------------------------------------------------------------------------
        ! Read number of atoms and origin
        !
        IF (io%lnode) THEN
            READ (unit, *, iostat=ios) nat, origin
            !
            CALL check_ios(ios)
            !
        END IF
        !
        CALL env_mp_bcast(nat, io%node, io%comm)
        !
        CALL env_mp_bcast(origin, io%node, io%comm)
        !
        !--------------------------------------------------------------------------------
        ! Read grid size and lattice vectors
        ! Note: vector(i) in units of bohr/grid(i)
        !
        IF (io%lnode) THEN
            !
            DO i = 1, 3
                READ (unit, *, iostat=ios) nr(i), at(i, :)
                !
                CALL check_ios(ios)
                !
            END DO
            !
        END IF
        !
        CALL env_mp_bcast(nr, io%node, io%comm)
        !
        CALL env_mp_bcast(at, io%node, io%comm)
        !
        IF (ANY(nr < 0)) THEN
            fact = 1.8897259886D0
            nr = -nr
        ELSE
            fact = 1.D0
        END IF
        !
        at = at * SPREAD(nr, 1, 3) * fact
        !
        !--------------------------------------------------------------------------------
        ! Read atomic numbers, charges, and positions
        !
        ALLOCATE (atomic_number(nat))
        ALLOCATE (charge(nat))
        ALLOCATE (tau(3, nat))
        !
        IF (io%lnode) THEN
            !
            DO i = 1, nat
                READ (unit, *, iostat=ios) atomic_number(i), charge(i), tau(:, i)
                !
                CALL check_ios(ios)
                !
            END DO
            !
        END IF
        !
        CALL env_mp_bcast(atomic_number, io%node, io%comm)
        !
        CALL env_mp_bcast(charge, io%node, io%comm)
        !
        CALL env_mp_bcast(tau, io%node, io%comm)
        !
        tau = tau * fact
        !
        !--------------------------------------------------------------------------------
        ! Determine number of species and assign species index per atom
        !
        ALLOCATE (species(nat))
        ALLOCATE (ityp(nat))
        ALLOCATE (combo(nat))
        ALLOCATE (charge_index(nat))
        !
        ntyp = 0
        species = 0
        combo = ''
        charge_index = 0
        !
        DO i = 1, nat
            WRITE (current_combo, '(I3, 1X, F10.6)') atomic_number(i), charge(i)
            !
            IF (.NOT. ANY(combo == current_combo)) THEN
                ntyp = ntyp + 1
                species(ntyp) = atomic_number(i)
                combo(ntyp) = current_combo
                charge_index(ntyp) = i
            END IF
            !
            ityp(i) = env_get_index(current_combo, combo)
        END DO
        !
        !--------------------------------------------------------------------------------
        ! Reduce charges to unique array
        !
        ALLOCATE (zv(ntyp))
        zv = 0.D0
        !
        DO i = 1, ntyp
            zv(i) = charge(charge_index(i))
        END DO
        !
        !--------------------------------------------------------------------------------
        ! Convert atomic numbers to atom labels
        !
        ALLOCATE (atom_label(ntyp))
        !
        CALL get_atom_labels(species, atom_label)
        !
        !--------------------------------------------------------------------------------
        ! Read density
        !
        IF (PRESENT(rho)) THEN
            !
            nnt = PRODUCT(nr)
            ALLOCATE (unsorted(nnt))
            ALLOCATE (rho(nnt))
            !
            IF (io%lnode) THEN
                READ (unit, *, iostat=ios) unsorted ! read unsorted cube density
                !
                CALL check_ios(ios)
                !
            END IF
            !
            CALL env_mp_bcast(unsorted, io%node, io%comm)
            !
            count = 0
            !
            DO i = 1, nr(1)
                !
                DO j = 1, nr(2)
                    !
                    DO k = 1, nr(3)
                        count = count + 1
                        l = i + (j - 1) * nr(1) + (k - 1) * nr(1) * nr(2)
                        rho(l) = unsorted(count)
                    END DO
                    !
                END DO
                !
            END DO
            !
        END IF
        !
        !--------------------------------------------------------------------------------
        !
        CLOSE (unit)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE read_cube
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE char2real(number,real_)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CHARACTER(LEN=80), INTENT(IN) :: number
        REAL(DP), INTENT(OUT) :: real_
        !
        CHARACTER(LEN=80) :: sub_name = 'check_EOF'
        !
        INTEGER :: i
        !
        !--------------------------------------------------------------------------------
        !
        READ(number,*) real_
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE char2real
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
    SUBROUTINE check_ios(ios)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: ios
        !
        CHARACTER(LEN=80) :: sub_name = 'check_EOF'
        !
        !--------------------------------------------------------------------------------
        !
        IF (ios == 0) THEN
            RETURN
        ELSE IF (ios > 0) THEN
            WRITE (io%unit, 1) "Error encountered while reading. Check input file."
        ELSE IF (ios < 0) THEN
            WRITE (io%unit, 1) "End-of-file encountered while reading. Check input file."
        END IF
        !
        CALL env_mp_abort(ios, io%comm)
        !
        !--------------------------------------------------------------------------------
        !
1       FORMAT(/, 5X, A,/)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE check_ios
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE parsers
!----------------------------------------------------------------------------------------
