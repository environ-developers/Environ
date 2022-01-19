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
MODULE tester_utils
    !------------------------------------------------------------------------------------
    !
    USE class_io, ONLY: io
    USE env_mp, ONLY: env_mp_bcast
    USE env_array_ops, ONLY: env_get_index
    !
    USE environ_param, ONLY: DP, tpi2
    !
    USE environ_api, ONLY: environ_interface, get_atom_labels
    !
    !------------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    PUBLIC :: from_input, from_cube, from_xml
    !
    !------------------------------------------------------------------------------------
CONTAINS
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE from_input(filename)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CHARACTER(LEN=*), INTENT(IN) :: filename
        !
        CHARACTER(LEN=80) :: sub_name = 'from_input'
        !
        !--------------------------------------------------------------------------------
        !

        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE from_input
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE from_cube(natoms, nspecies, ispecies, label, ionic_charge, nelectrons, &
                         atom_position, system_origin, grid, lattice, density, cubefile)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: cubefile
        !
        INTEGER, INTENT(OUT) :: natoms, nspecies
        REAL(DP), INTENT(OUT) :: nelectrons
        INTEGER, ALLOCATABLE, INTENT(OUT) :: ispecies(:)
        CHARACTER(LEN=2), ALLOCATABLE, INTENT(OUT) :: label(:)
        REAL(DP), INTENT(OUT) :: system_origin(3)
        INTEGER, INTENT(OUT) :: grid(3)
        REAL(DP), INTENT(OUT) :: lattice(3, 3)
        REAL(DP), ALLOCATABLE, INTENT(OUT) :: ionic_charge(:)
        REAL(DP), ALLOCATABLE, INTENT(OUT) :: atom_position(:, :)
        REAL(DP), ALLOCATABLE, INTENT(OUT) :: density(:)
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
        CHARACTER(LEN=80) :: sub_name = 'from_cube'
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
        CALL io%header('Reading cube data from '//TRIM(filename))
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
        IF (io%lnode) READ (unit, *) natoms, system_origin
        !
        CALL env_mp_bcast(natoms, io%node, io%comm)
        !
        CALL env_mp_bcast(system_origin, io%node, io%comm)
        !
        !--------------------------------------------------------------------------------
        ! Read grid size and lattice vectors
        ! Note: vector(i) in units of bohr/grid(i)
        !
        IF (io%lnode) THEN
            !
            DO i = 1, 3
                READ (unit, *) grid(i), lattice(i, :)
            END DO
            !
        END IF
        !
        CALL env_mp_bcast(grid, io%node, io%comm)
        !
        CALL env_mp_bcast(lattice, io%node, io%comm)
        !
        lattice = lattice * SPREAD(grid, 1, 3)
        !
        !--------------------------------------------------------------------------------
        ! Read atomic numbers, charges, and positions
        !
        ALLOCATE (atomic_number(natoms))
        ALLOCATE (charge(natoms))
        ALLOCATE (atom_position(3, natoms))
        !
        IF (io%lnode) THEN
            !
            DO i = 1, natoms
                READ (unit, *) atomic_number(i), charge(i), atom_position(:, i)
            END DO
            !
        END IF
        !
        CALL env_mp_bcast(atomic_number, io%node, io%comm)
        !
        CALL env_mp_bcast(charge, io%node, io%comm)
        !
        CALL env_mp_bcast(atom_position, io%node, io%comm)
        !
        nelectrons = SUM(charge)
        !
        !--------------------------------------------------------------------------------
        ! Determine number of species and assign species index per atom
        !
        ALLOCATE (species(natoms))
        ALLOCATE (ispecies(natoms))
        !
        nspecies = 0
        species = 0
        !
        DO i = 1, natoms
            !
            IF (.NOT. ANY(species == atomic_number(i))) THEN
                nspecies = nspecies + 1
                species(i) = atomic_number(i)
            END IF
            !
            ispecies(i) = env_get_index(atomic_number(i), species)
        END DO
        !
        !--------------------------------------------------------------------------------
        ! Reduce charges to unique array
        !
        ALLOCATE (ionic_charge(nspecies))
        current_charge = 0.D0
        count = 0
        !
        DO i = 1, natoms
            !
            IF (charge(i) /= current_charge) THEN
                count = count + 1
                ionic_charge(count) = charge(i)
                current_charge = charge(i)
            END IF
            !
        END DO
        !
        !--------------------------------------------------------------------------------
        ! Convert atomic numbers to atom labels
        !
        ALLOCATE (label(nspecies))
        !
        CALL get_atom_labels(species, label)
        !
        !--------------------------------------------------------------------------------
        ! Read density
        !
        nnt = PRODUCT(grid)
        ALLOCATE (unsorted(nnt))
        ALLOCATE (density(nnt))
        !
        IF (io%lnode) READ (unit, *) unsorted ! read unsorted cube density
        !
        CALL env_mp_bcast(unsorted, io%node, io%comm)
        !
        count = 0
        !
        DO i = 1, grid(1)
            !
            DO j = 1, grid(2)
                !
                DO k = 1, grid(3)
                    count = count + 1
                    l = i + (j - 1) * grid(1) + (k - 1) * grid(1) * grid(2)
                    density(l) = unsorted(count)
                END DO
                !
            END DO
            !
        END DO
        !
        !--------------------------------------------------------------------------------
        !
        CLOSE (unit)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE from_cube
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE from_xml(filename)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CHARACTER(LEN=*), INTENT(IN) :: filename
        !
        CHARACTER(LEN=80) :: sub_name = 'from_xml'
        !
        !--------------------------------------------------------------------------------
        !

        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE from_xml
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE tester_utils
!----------------------------------------------------------------------------------------
