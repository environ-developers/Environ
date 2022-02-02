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
    USE env_mp, ONLY: env_mp_bcast
    USE env_array_ops, ONLY: env_get_index
    !
    USE environ_param, ONLY: DP
    !
    USE environ_api, ONLY: get_atom_labels
    !
    !------------------------------------------------------------------------------------
CONTAINS
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE read_cube(nat, ntyp, ityp, atom_label, zv, nelec, tau, origin, nr, at, &
                         rho, cubefile)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: cubefile
        !
        INTEGER, INTENT(OUT) :: nat, ntyp
        REAL(DP), INTENT(OUT) :: nelec
        INTEGER, ALLOCATABLE, INTENT(OUT) :: ityp(:)
        CHARACTER(LEN=2), ALLOCATABLE, INTENT(OUT) :: atom_label(:)
        REAL(DP), INTENT(OUT) :: origin(3)
        INTEGER, INTENT(OUT) :: nr(3)
        REAL(DP), INTENT(OUT) :: at(3, 3)
        REAL(DP), ALLOCATABLE, INTENT(OUT) :: zv(:)
        REAL(DP), ALLOCATABLE, INTENT(OUT) :: tau(:, :)
        REAL(DP), ALLOCATABLE, OPTIONAL, INTENT(OUT) :: rho(:)
        !
        INTEGER :: i, j, k, l, count, unit, current, nnt
        LOGICAL :: ext
        !
        REAL(DP) :: current_charge
        !
        INTEGER, ALLOCATABLE :: atomic_number(:)
        REAL(DP), ALLOCATABLE :: charge(:) ! unprocessed charges
        INTEGER, ALLOCATABLE :: species(:) ! register for checked atomic numbers
        REAL(DP), ALLOCATABLE :: unsorted(:) ! unsorted density read from file
        !
        CHARACTER(LEN=80) :: filename = 'density.cube'
        !
        CHARACTER(LEN=80) :: sub_name = 'read_cube'
        !
        !--------------------------------------------------------------------------------
        ! Open cube file
        !
        IF (PRESENT(cubefile)) filename = cubefile
        !
        unit = io%find_free_unit()
        INQUIRE (file=TRIM(filename), exist=ext)
        !
        IF (.NOT. ext) CALL io%error(sub_name, TRIM(filename)//' not found', 1)
        !
        OPEN (unit=unit, file=TRIM(filename), status='old')
        !
        !--------------------------------------------------------------------------------
        !
        CALL io%writer('Reading cube data from '//TRIM(filename))
        !
        CALL io%writer('') ! blank line
        !
        !--------------------------------------------------------------------------------
        ! Discard comment lines
        !
        IF (io%lnode) THEN
            !
            DO i = 1, 2
                READ (unit, *)
            END DO
            !
        END IF
        !
        !--------------------------------------------------------------------------------
        ! Read number of atoms and origin
        !
        IF (io%lnode) READ (unit, *) nat, origin
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
                READ (unit, *) nr(i), at(i, :)
            END DO
            !
        END IF
        !
        CALL env_mp_bcast(nr, io%node, io%comm)
        !
        CALL env_mp_bcast(at, io%node, io%comm)
        !
        at = at * SPREAD(nr, 1, 3)
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
                READ (unit, *) atomic_number(i), charge(i), tau(:, i)
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
        nelec = SUM(charge)
        !
        !--------------------------------------------------------------------------------
        ! Determine number of species and assign species index per atom
        !
        ALLOCATE (species(nat))
        ALLOCATE (ityp(nat))
        !
        ntyp = 0
        species = 0
        !
        DO i = 1, nat
            !
            IF (.NOT. ANY(species == atomic_number(i))) THEN
                ntyp = ntyp + 1
                species(i) = atomic_number(i)
            END IF
            !
            ityp(i) = env_get_index(atomic_number(i), species)
        END DO
        !
        !--------------------------------------------------------------------------------
        ! Reduce charges to unique array
        !
        ALLOCATE (zv(ntyp))
        current_charge = 0.D0
        count = 0
        !
        DO i = 1, nat
            !
            IF (charge(i) /= current_charge) THEN
                count = count + 1
                zv(count) = charge(i)
                current_charge = charge(i)
            END IF
            !
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
            IF (io%lnode) READ (unit, *) unsorted ! read unsorted cube density
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
    !
    !------------------------------------------------------------------------------------
END MODULE parsers
!----------------------------------------------------------------------------------------
