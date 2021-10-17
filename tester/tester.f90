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
    USE class_setup
    USE class_environ
    USE class_calculator, ONLY: calc
    USE class_destructor, ONLY: clean
    !
    USE environ_input
    USE env_base_input
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
    ! Global objects
    !
    TYPE(environ_setup), SAVE :: setup1
    TYPE(environ_obj), SAVE :: env1
    !
    !------------------------------------------------------------------------------------
    !
    INTEGER :: ierr ! mpi error
    !
    INTEGER :: i, j ! loop indices
    !
    !------------------------------------------------------------------------------------
    ! MP parameters
    !
    INTEGER :: comm
    LOGICAL :: ionode
    !
    INTEGER, EXTERNAL :: env_mp_rank
    !
    !------------------------------------------------------------------------------------
    ! Cell parameters
    !
    REAL(DP) :: L
    REAL(DP) :: at(3, 3)
    !
    INTEGER :: ecutrho
    !
    !------------------------------------------------------------------------------------
    ! Ion parameters
    !
    INTEGER :: nat
    INTEGER :: ntyp
    INTEGER, ALLOCATABLE :: ityp(:)
    CHARACTER(LEN=3), ALLOCATABLE :: atom_label(:)
    REAL(DP), ALLOCATABLE :: zv(:)
    REAL(DP), ALLOCATABLE :: tau(:, :)
    !
    !------------------------------------------------------------------------------------
    ! PBC Parameters
    !
    LOGICAL :: apply_pbc_correction
    !
    !------------------------------------------------------------------------------------
    ! Results
    !
    REAL(DP) :: energy
    REAL(DP), ALLOCATABLE :: force(:, :)
    !
    !------------------------------------------------------------------------------------
    !
    !                                    MAIN PROGRAM
    !
    !------------------------------------------------------------------------------------
    ! Initialize MPI
    !
    comm = MPI_COMM_WORLD
    !
    CALL MPI_Init(ierr)
    !
    IF (ierr /= 0) CALL env_mp_abort(ierr, comm)
    !
    CALL env_allocate_mp_buffers()
    !
    !------------------------------------------------------------------------------------
    ! Lattice setup
    !
    L = 40.D0
    !
    CALL set_lattice(L) ! assumes cubic for now
    !
    !------------------------------------------------------------------------------------
    ! Ions setup
    !
    nat = 1
    ntyp = 1
    !
    ALLOCATE (ityp(nat))
    ALLOCATE (atom_label(nat))
    ALLOCATE (zv(nat))
    ALLOCATE (tau(3, nat))
    !
    ityp = 1
    atom_label = 'H'
    zv = 1.D0
    !
    tau = 0.D0
    !
    DO i = 1, nat
        !
        DO j = 1, 3
            tau(j, i) = L / 2.D0
        END DO
        !
    END DO
    !
    !------------------------------------------------------------------------------------
    ! Initialize I/O
    !
    ionode = env_mp_rank(comm) == 0
    !
    CALL io%init(ionode, 0, comm, 6, .TRUE.)
    !
    !------------------------------------------------------------------------------------
    ! Environ input
    !
    CALL read_environ_input()
    !
    io%verbosity = verbose
    !
    !------------------------------------------------------------------------------------
    ! Initialize Environ setup and base
    !
    ecutrho = 300
    apply_pbc_correction = .TRUE.
    !
    CALL init(setup1, env1, ecutrho)
    !
    !------------------------------------------------------------------------------------
    !
    CALL remove_self_interaction(env1)
    !
    !------------------------------------------------------------------------------------
    ! Test
    !
    CALL calc%potential(env1, .TRUE.)
    !
    CALL calc%energy(env1, energy)
    !
    ALLOCATE (force(3, nat))
    !
    CALL calc%force(env1, nat, force)
    !
    CALL print_results()
    !
    !------------------------------------------------------------------------------------
CONTAINS
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE set_lattice(size_in)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP), INTENT(IN) :: size_in
        !
        !--------------------------------------------------------------------------------
        !
        at = 0.D0
        !
        DO i = 1, 3
            at(i, i) = size_in
        END DO
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE set_lattice
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init(setup, env, ecut)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: ecut
        !
        TYPE(environ_setup), INTENT(INOUT) :: setup
        TYPE(environ_obj), INTENT(INOUT) :: env
        !
        REAL(DP) :: gcutm
        !
        !--------------------------------------------------------------------------------
        !
        gcutm = ecut / tpi2
        !
        CALL setup%init()
        !
        CALL setup%init_cell(gcutm, io%comm, at)
        ! !
        CALL setup%init_cores(gcutm, apply_pbc_correction)
        !
        CALL env%init(setup, 0, nat, ntyp, atom_label, ityp, zv)
        !
        CALL env%update_ions(nat, tau)
        !
        CALL setup%update_cell(at)
        !
        CALL env%update_cell_dependent_quantities()
        !
        CALL setup%end_cell_update()
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE remove_self_interaction(env)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_obj), INTENT(INOUT) :: env
        !
        CHARACTER(LEN=80) :: sub_name = 'remove_self_interaction'
        !
        !--------------------------------------------------------------------------------
        !
        CALL env%system_charges%add(externals=env%externals)
        !
        env%system_charges%include_ions = .FALSE.
        !
        CALL env%system_charges%update()
        !
        env%environment_charges%include_ions = .FALSE.
        !
        CALL env%environment_charges%update()
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE remove_self_interaction
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE print_results()
        !--------------------------------------------------------------------------------
        !
        IF (io%lnode) THEN
            !
            PRINT '(5X, A15, F14.8, /)', 'total energy = ', energy
            !
            PRINT '(5X, A6, /)', 'forces'
            !
            DO i = 1, nat
                PRINT '(5X, 3F14.8)', force(:, i)
            END DO
            !
            PRINT *, ''
            !
        END IF
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE print_results
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END PROGRAM tester
!----------------------------------------------------------------------------------------
