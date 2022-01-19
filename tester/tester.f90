!----------------------------------------------------------------------------------------
!
! Copyright (C) 2021 ENVIRON (www.quantum-environ.org)
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
! Authors: Edan Bainglass (Department of Physics, UNT)
!
!----------------------------------------------------------------------------------------
!>
!!
!----------------------------------------------------------------------------------------
PROGRAM tester
    !------------------------------------------------------------------------------------
    !
    USE env_parallel_include
    USE class_io, ONLY: io
    !
    USE environ_param, ONLY: DP, tpi2
    !
    USE environ_api, ONLY: environ_interface, get_atom_labels
    !
    USE tester_utils, ONLY: from_cube
    !
    !------------------------------------------------------------------------------------
    !
    !                                     PARAMETERS
    !
    !------------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    !------------------------------------------------------------------------------------
    ! Declare interface
    !
    TYPE(environ_interface) :: environ
    !
    !------------------------------------------------------------------------------------
    ! Declare MPI parameters
    !
    INTEGER :: comm, ionode
    LOGICAL :: lnode
    LOGICAL :: initialized
    INTEGER :: ierr ! mpi error
    !
    INTEGER, EXTERNAL :: env_mp_rank
    !
    !------------------------------------------------------------------------------------
    ! Initialize MPI
    !
    comm = MPI_COMM_WORLD
    !
    CALL MPI_Init(ierr)
    !
    IF (ierr /= 0) CALL env_mp_stop(8000)
    !
    !------------------------------------------------------------------------------------
    ! Initialize I/O
    !
    ionode = 0
    lnode = env_mp_rank(comm) == ionode
    !
    CALL environ%init_interface()
    !
    CALL environ%init_io(lnode, ionode, comm, 6, .FALSE.)
    !
    !------------------------------------------------------------------------------------
    !
    !                                   SELECT PROGRAM
    !
    !------------------------------------------------------------------------------------
    !

    !
    !------------------------------------------------------------------------------------
    ! Clean up
    !
    CALL environ%destroy()
    !
    CALL MPI_Finalize(ierr)
    !
    IF (ierr /= 0) CALL env_mp_stop(8001)
    !
    !------------------------------------------------------------------------------------
CONTAINS
    !------------------------------------------------------------------------------------
    !>
    !! An Environ calculation on a "frozen" density provided in a cube file
    !!
    !! Default values:
    !!
    !! - inputfile          = 'environ.in'
    !! - cubefile           = 'density.cube'
    !! - use_pbc_correction = .FALSE.
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE run_environ_from_cube(inputfile, cubefile, use_pbc_correction)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: inputfile
        CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: cubefile
        LOGICAL, OPTIONAL, INTENT(IN) :: use_pbc_correction
        !
        INTEGER :: nat
        INTEGER :: ntyp
        REAL(DP) :: nelec
        INTEGER, ALLOCATABLE :: ityp(:)
        CHARACTER(LEN=2), ALLOCATABLE :: label(:)
        REAL(DP), ALLOCATABLE :: zv(:)
        REAL(DP), ALLOCATABLE :: tau(:, :)
        REAL(DP) :: origin(3)
        INTEGER :: nr(3)
        REAL(DP) :: at(3, 3)
        REAL(DP), ALLOCATABLE :: rho(:)
        !
        REAL(DP), ALLOCATABLE :: env_potential(:)
        REAL(DP) :: env_energy
        !
        REAL(DP) :: volume, avg_dvtot, avg_velectrostatic
        !
        INTEGER :: i, nrmax
        !
        REAL(DP) :: gcutm, tmp, a1(3), sumat2, est
        !
        CHARACTER(LEN=80) :: sub_name = 'run_environ_from_cube'
        !
        !--------------------------------------------------------------------------------
        ! Get input parameters from file
        !
        CALL from_cube(nat, ntyp, ityp, label, zv, nelec, tau, origin, nr, at, rho, &
                       cubefile)
        !
        !--------------------------------------------------------------------------------
        ! Initialize Environ
        !
        CALL environ%read_input(inputfile)
        !
        CALL environ%setup%init(use_pbc_correction)
        !
        CALL environ%setup%init_cell(io%comm, at, nr=nr)
        !
        CALL environ%setup%init_cores()
        !
        CALL environ%main%init(nat, ntyp, label, ityp, zv)
        !
        !--------------------------------------------------------------------------------
        ! Update Environ parameters
        !
        CALL environ%main%update_ions(nat, tau, origin)
        !
        CALL environ%update_electrons(rho, nelec=nelec, lscatter=.TRUE.)
        !
        !--------------------------------------------------------------------------------
        ! Compute potential
        !
        ALLOCATE (env_potential(SIZE(rho)))
        !
        CALL environ%calc_potential(.TRUE., env_potential, lgather=.TRUE.)
        !
        CALL environ%calc%energy(env_energy)
        !
        !--------------------------------------------------------------------------------
        ! Print results
        !
        volume = environ%setup%system_cell%omega
        avg_dvtot = environ%main%dvtot%integrate() / volume
        !
        IF (environ%setup%has_electrostatics()) &
            avg_velectrostatic = environ%main%velectrostatic%integrate() / volume
        !
        IF (io%lnode) THEN
            !
            CALL environ%main%print_energies('PW', .FALSE.)
            !
            WRITE (io%unit, 1000), env_energy
            WRITE (io%unit, 1001), environ%main%system_charges%charge
            WRITE (io%unit, 1002), avg_dvtot
            !
            IF (environ%setup%has_electrostatics()) &
                WRITE (io%unit, 1003), avg_velectrostatic
            !
            WRITE (io%unit, *) ! final blank line
        END IF
        !
        !--------------------------------------------------------------------------------
        !
1000    FORMAT(5X, "total energy              =", F17.8, " Ry")
        !
1001    FORMAT(/, 5X, "total charge              =", F17.8, " a.u.")
        !
1002    FORMAT(/, 5X, "average total potential           =", F17.8, " Ry")

1003    FORMAT(5X, "average electrostatic potential   =", F17.8, " Ry")
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE run_environ_from_cube
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END PROGRAM tester
!----------------------------------------------------------------------------------------
