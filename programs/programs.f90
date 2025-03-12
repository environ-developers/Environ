!----------------------------------------------------------------------------------------
!
! Copyright (C) 2018-2022 ENVIRON (www.quantum-environ.org)
!
!----------------------------------------------------------------------------------------
!
!     This file is part of Environ version 3.0
!
!     Environ 3.0 is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation, either version 2 of the License, or
!     (at your option) any later version.
!
!     Environ 3.0 is distributed in the hope that it will be useful,
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
MODULE programs
    !------------------------------------------------------------------------------------
    !
    USE env_parallel_include
    USE env_mp, ONLY: env_mp_rank, env_mp_stop, env_mp_sum, env_mp_size
    !
    USE env_mytime, ONLY: env_f_wall
    !
    USE class_io, ONLY: io
    !
    USE environ_param, ONLY: DP
    !
    USE environ_api, ONLY: environ_interface
    !
    USE class_cell, ONLY: environ_cell
    USE class_density, ONLY: environ_density
    USE class_gradient, ONLY: environ_gradient
    !
    USE class_boundary_ionic, ONLY: environ_boundary_ionic
    USE class_boundary_electronic, ONLY: environ_boundary_electronic
    !
    USE env_write_cube, ONLY: write_cube
    !
    USE cmdline_args
    !
    USE prog_utils
    !
    USE parsers
    !
    USE environ_param, only: BOHR_RADIUS_ANGS
    !
    USE tools_math, only: environ_erfc
    !
    USE env_base_input, ONLY: solvent_mode, alpha, softness
    !
    !------------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    PRIVATE
    !
    PUBLIC :: run_tester, run_environ_from_cube, run_descriptors_generator
    PUBLIC :: run_environ_with_aims
    !
    PUBLIC :: initial_setup, clean_up, print_available_programs
    !
    !------------------------------------------------------------------------------------
    ! Declare interface
    !
    TYPE(environ_interface), TARGET :: environ
    !
    !------------------------------------------------------------------------------------
CONTAINS
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !                                  PROGRAM ROUTINES
    !
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE print_available_programs()
        !--------------------------------------------------------------------------------
        !
        IF (io%lnode) &
            PRINT '(4(/, 5X, A), /)', &
            'Available calculations:', &
            '- tester', &
            '- from_cube', &
            '- descriptors', &
            '- with_aims'
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE print_available_programs
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE run_tester()
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !

        !
        !--------------------------------------------------------------------------------
        !

        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE run_tester
    !------------------------------------------------------------------------------------
    !>
    !! An Environ calculation communicating with an FHI-aims calculation
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE run_environ_with_aims()
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER :: comm, myid, mpi_size, aims_to_env_stat, env_to_aims_stat, ierr, &
                   crash_file_size
        CHARACTER(34) :: fname_aims_env, fname_env_aims
        LOGICAL :: exist, finished, initialized, updated
        INTEGER, ALLOCATABLE :: stat_all_tasks(:)
        !
        !--------------------------------------------------------------------------------
        ! This is copied over from aims. Not all values are actually used on this side.
        ! See aims code for meaning of inividual values.
        !
        INTEGER, PARAMETER :: stat_env_err              = -101
        INTEGER, PARAMETER :: stat_env_err_not_ini      = -102
        INTEGER, PARAMETER :: stat_env_err_not_updt     = -103
        INTEGER, PARAMETER :: stat_env_err_hand_ini     = -104
        INTEGER, PARAMETER :: stat_env_err_hand_not_ini = -105
        INTEGER, PARAMETER :: stat_env_err_inval        = -106
        INTEGER, PARAMETER :: stat_env_stopped          = -107
        INTEGER, PARAMETER :: stat_env_undef_env        = -1
        INTEGER, PARAMETER :: stat_env_undef_aims       = -2
        INTEGER, PARAMETER :: stat_env_ini              = 1
        INTEGER, PARAMETER :: stat_env_updates          = 2
        INTEGER, PARAMETER :: stat_env_grid_getters     = 3
        INTEGER, PARAMETER :: stat_env_calc             = 4
        INTEGER, PARAMETER :: stat_env_cleanup          = 5
        INTEGER, PARAMETER :: stat_env_read_success     = 101
        !
        CHARACTER(LEN=80) :: routine = 'run_environ_with_aims'
        !
        !--------------------------------------------------------------------------------
        !
        initialized = .false.
        updated = .false.
        !
        !--------------------------------------------------------------------------------
        !
        comm = MPI_COMM_WORLD
        myid = env_mp_rank(comm)
        mpi_size = env_mp_size(comm)
        !
        ALLOCATE(stat_all_tasks(mpi_size))
        !
        !--------------------------------------------------------------------------------
        ! Filehandles to send and receive status to / from aims
        !
        WRITE(fname_aims_env, '(A,I0.10)') 'aims_to_env_stat_info_id', myid
        WRITE(fname_env_aims, '(A,I0.10)') 'env_to_aims_stat_info_id', myid
        !
        !--------------------------------------------------------------------------------
        ! Keep waiting for status signals from aims until cleanup has been called
        !
        finished = .false.
        !
        steps: DO
            !
            !----------------------------------------------------------------------------
            ! Set status that will be returned to aims to undefined. If work subroutines
            ! do not overwrite it, aims will know that work has not been completed
            ! correctly, and will terminate. A signal to terminate Environ will be sent
            ! from aims, no need to check for the return status here.
            !
            env_to_aims_stat = stat_env_undef_env
            !
            !----------------------------------------------------------------------------
            ! Idle loop, wait until status file exists
            !
            check_stat: DO
                !
                !------------------------------------------------------------------------
                ! Check if aims or a different process wrote a crash file. In this case,
                ! we can come to a controlled termination where aims writes some
                ! diagnostic output.
                !
                INQUIRE(FILE='AIMS_ENV_CRASH', EXIST=exist)
                !
                IF (exist) THEN
                    !
                    aims_to_env_stat = stat_env_cleanup
                    EXIT check_stat
                    !
                ENDIF
                !
                !------------------------------------------------------------------------
                ! Check if an error happened in aims, outside of the Environ interface,
                ! and was redirected to AIMS_CRASH. In this case, we cannot write any
                ! meaningful diagnostic output, but we can still terminate Environ. If
                ! user did not redirect stderr, it is up to them to kill Environ after
                ! aims has crashed.
                !
                INQUIRE(FILE='AIMS_CRASH', EXIST=exist)
                !
                IF (exist) then
                    !
                    INQUIRE(FILE='AIMS_CRASH', SIZE=crash_file_size)
                    !
                    IF (crash_file_size .gt. 0) THEN
                        !
                        aims_to_env_stat = stat_env_cleanup
                        exit check_stat
                        !
                    ENDIF
                    !
                ENDIF
                !
                !------------------------------------------------------------------------
                ! If no process sent a crash file, check for regular communication file
                ! from aims
                !
                INQUIRE(FILE=fname_aims_env, EXIST=exist)
                !
                IF (exist) THEN
                    !
                    OPEN(UNIT=146, FILE=fname_aims_env, STATUS='unknown', &
                         ACTION='read', FORM='unformatted', ACCESS='stream')
                    READ(146) aims_to_env_stat
                    CLOSE(146, STATUS='delete')
                    !
                    EXIT check_stat
                    !
                ELSE
                    !
                    CALL sleep(1)
                    !
                ENDIF
                !
            ENDDO check_stat
            !
            !----------------------------------------------------------------------------
            ! Status read, check if all processes were sent the same status
            !
            stat_all_tasks(:) = 0
            stat_all_tasks(myid+1) = aims_to_env_stat
            CALL MPI_ALLREDUCE(MPI_IN_PLACE, stat_all_tasks, mpi_size, MPI_INTEGER, &
                               MPI_SUM, comm, ierr)
            !
            !----------------------------------------------------------------------------
            ! If different processes exited idle loop with different status, this is
            ! probably because some of them caught a crash file. In this case, terminate
            ! Environ. Aims will terminate on its own; if not because it sent the crash
            ! file, then because Environ will return a different status than it was sent
            ! on at least one process. No need to do anything fancy here.
            !
            IF (any(stat_all_tasks(:mpi_size) .ne. stat_all_tasks(1))) &
                aims_to_env_stat = stat_env_cleanup
            !
            !----------------------------------------------------------------------------
            ! Run the program step that was requested by aims
            !
            SELECT CASE (aims_to_env_stat)
                !
            CASE (stat_env_ini)
                IF (io%lnode) WRITE(io%unit, '(A)') &
                        '  [Environ side interface for FHI-aims] run initialization'
                CALL environ_aims_initializations(env_to_aims_stat)
                !
            CASE (stat_env_updates)
                IF (io%lnode) WRITE(io%unit, '(A)') &
                        '  [Environ side interface for FHI-aims] run updates'
                CALL environ_aims_updates(env_to_aims_stat)
                !
            CASE (stat_env_grid_getters)
                IF (io%lnode) WRITE(io%unit, '(A)') &
                        '  [Environ side interface for FHI-aims] run grid getters'
                CALL environ_aims_grid_getters(env_to_aims_stat)
                !
            CASE (stat_env_calc)
                IF (io%lnode) WRITE(io%unit, '(A)') &
                        '  [Environ side interface for FHI-aims] run calculators'
                CALL environ_aims_calculators(env_to_aims_stat)
                !
            CASE (stat_env_cleanup)
                IF (io%lnode) WRITE(io%unit, '(A)') &
                        '  [Environ side interface for FHI-aims] run cleanup'
                CALL environ_aims_cleanup(env_to_aims_stat)
                finished = .true.
                !
            CASE DEFAULT
                !
                !------------------------------------------------------------------------
                ! In principle, this should already have been caught by aims before
                ! any signal is sent to Environ. The only case when this error should
                ! occur is when a new status is implemented in the version of aims that
                ! was used in the run, but is not implemented in this Environ version.
                !
                ! Generally, erroneous status should be handled by aims.
                !
                IF (io%lnode) WRITE(io%unit, '(A)') &
                        '  [Environ side interface for FHI-aims] unknown'//&
                        ' input status. Check version compatibility.'
                env_to_aims_stat = stat_env_err
                !
            END SELECT
            !
            !----------------------------------------------------------------------------
            ! Pass status back to aims
            !
            OPEN(UNIT=146, FILE=fname_env_aims, STATUS='unknown', ACTION='write', &
                 FORM='unformatted', ACCESS='stream')
            WRITE(146) env_to_aims_stat
            CLOSE(146)
            !
            !----------------------------------------------------------------------------
            !
            IF (finished) EXIT steps
            !
        ENDDO steps
        !
        DEALLOCATE(stat_all_tasks)
        !
        !--------------------------------------------------------------------------------
    CONTAINS
        !--------------------------------------------------------------------------------
        !>
        !!
        !--------------------------------------------------------------------------------
        SUBROUTINE environ_aims_initializations(stat)
            !----------------------------------------------------------------------------
            !
            IMPLICIT NONE
            !
            INTEGER, INTENT(INOUT) :: stat
            !         
            CHARACTER(34) :: fname
            !
            INTEGER :: n_atoms, n_species
            REAL(DP) :: lattice_vector(3,3)
            INTEGER, ALLOCATABLE :: species(:)
            REAL(DP), ALLOCATABLE :: species_z(:)
            !
            !----------------------------------------------------------------------------
            ! Read data from aims
            !
            WRITE(fname, '(A,I0.10)') 'aims_to_env_init_info_id', myid
            OPEN(UNIT=147, FILE=fname, STATUS='unknown', ACTION='read', &
                 FORM='unformatted', ACCESS='stream')
            !
            READ(147) n_atoms
            READ(147) n_species
            !
            ALLOCATE(species(n_atoms))
            ALLOCATE(species_z(n_species))
            !
            READ(147) species
            READ(147) species_z
            READ(147) lattice_vector
            !
            CLOSE(147, STATUS='delete')
            !
            !----------------------------------------------------------------------------
            ! Run Environ initializations
            !
            CALL environ%read_input(nsx=n_species)
            CALL environ%setup%init()
            CALL environ%setup%print_summary()
            CALL environ%setup%init_cell(comm, lattice_vector)
            CALL environ%setup%init_numerical()
            CALL environ%main%init(n_atoms, n_species, species, species_z, &
                                   number=NINT(species_z))
            !
            DEALLOCATE(species)
            DEALLOCATE(species_z)
            !
            initialized = .true.
            !
            !----------------------------------------------------------------------------
            ! Pass 'safe' region of each species to FHI-aims
            !
            WRITE(fname, '(A,I0.10)') 'env_to_aims_init_info_id', myid
            OPEN(UNIT=147, FILE=fname, STATUS='unknown', ACTION='write', &
                 FORM='unformatted', ACCESS='stream')
            !
            SELECT CASE (solvent_mode)
                !
            CASE ('electronic', 'full')
                !
                !------------------------------------------------------------------------
                ! Flag for aims to get radii from isodensity
                !
                WRITE(147) .TRUE.
                !
                SELECT TYPE (slvnt => environ%main%solvent)
                    !
                TYPE IS (environ_boundary_electronic)
                    WRITE(147) slvnt%rhomax
                    !
                CLASS DEFAULT
                    stat = stat_env_err
                    CLOSE(147)
                    RETURN
                    !
                END SELECT
                !
            CASE ('ionic', 'system')
                !
                !------------------------------------------------------------------------
                ! Flag for aims to not get radii from isodensity
                !
                WRITE(147) .FALSE.
                !
                IF (ANY(environ%main%system_ions%iontype(:)%solvationrad .le. &
                        2.D0*softness/alpha)) THEN
                    stat = stat_env_err
                    CLOSE(147)
                    RETURN
                ENDIF
                !
                !------------------------------------------------------------------------
                ! Get radius where (1-boundary function) = 0.5*(1+erf(-2)) < 1.e-2, then
                ! take half of that radius
                !
                WRITE(147) (environ%main%system_ions%iontype(:)%solvationrad * alpha &
                            - 2.D0 * softness) * 0.5D0
                !
            CASE DEFAULT
                stat = stat_env_err
                CLOSE(147)
                RETURN
                !
            END SELECT
            !
            CLOSE(147)
            !
            stat = stat_env_ini
            !
            !----------------------------------------------------------------------------
        END SUBROUTINE environ_aims_initializations
        !--------------------------------------------------------------------------------
        !>
        !!
        !--------------------------------------------------------------------------------
        SUBROUTINE environ_aims_updates(stat)
            !----------------------------------------------------------------------------
            !
            IMPLICIT NONE
            !
            INTEGER, INTENT(INOUT) :: stat
            !
            CHARACTER(34) :: fname
            !
            INTEGER :: n_atoms, ind
            REAL(DP), ALLOCATABLE :: coords_pass(:,:), lattice_vector(:,:)
            !
            !----------------------------------------------------------------------------
            !
            if (.not. initialized) then
                stat = stat_env_err_not_ini
                return
            endif
            !
            !----------------------------------------------------------------------------
            ! Read data from aims
            !
            WRITE(fname, '(A,I0.10)') 'aims_to_env_updt_info_id', myid
            OPEN(UNIT=147, FILE=fname, STATUS='unknown', ACTION='read', &
                 FORM='unformatted', ACCESS='stream')
            !
            READ(147) n_atoms
            !
            ALLOCATE(coords_pass(3,n_atoms))
            ALLOCATE(lattice_vector(3,3))
            !
            READ(147) coords_pass
            READ(147) lattice_vector
            !
            CLOSE(147, STATUS='delete')
            !
            !----------------------------------------------------------------------------
            ! Run Environ updates
            !
            CALL environ%main%update_ions(n_atoms, coords_pass)
            CALL environ%setup%update_cell(lattice_vector)
            !
            !----------------------------------------------------------------------------
            ! Nothing to pass back to aims at this point, just deallocate and return
            !
            DEALLOCATE(coords_pass)
            DEALLOCATE(lattice_vector)
            !
            !----------------------------------------------------------------------------
            ! Caution is advised here. This flag only checks if Environ has ever been
            ! updated at all. It does not check if it was updated since the last
            ! geometry change in aims.
            !
            updated = .true.
            stat = stat_env_updates
            !
            !----------------------------------------------------------------------------
        END SUBROUTINE environ_aims_updates
        !--------------------------------------------------------------------------------
        !>
        !!
        !--------------------------------------------------------------------------------
        SUBROUTINE environ_aims_grid_getters(stat)
            !----------------------------------------------------------------------------
            !
            IMPLICIT NONE
            !
            INTEGER, INTENT(INOUT) :: stat
            !
            CHARACTER(34) :: fname
            !
            INTEGER :: nnr, ir_end, nnt, ind, nr(3), nrx(3), nproc2, nproc3
            REAL(DP), ALLOCATABLE :: points(:,:)
            INTEGER, ALLOCATABLE :: i0r2p(:), i0r3p(:), nr2p(:), nr3p(:)
            !
            !----------------------------------------------------------------------------
            !
            if (.not. initialized) then
                stat = stat_env_err_not_ini
                return
            endif
            if (.not. updated) then
                stat = stat_env_err_not_updt
                return
            endif
            !
            !----------------------------------------------------------------------------
            ! Gather data
            !
            nnr = environ%setup%get_nnr()
            ir_end = environ%setup%get_ir_end()
            nnt = environ%setup%get_nnt()
            DO ind = 1, 3
                nr(ind) = environ%setup%get_nr(ind)
            ENDDO
            nrx = environ%setup%system_cell%nrx
            nproc2 = environ%setup%system_cell%dfft%nproc2
            nproc3 = environ%setup%system_cell%dfft%nproc3
            !
            ALLOCATE(i0r2p(nproc2))
            ALLOCATE(i0r3p(nproc3))
            ALLOCATE(nr2p(nproc2))
            ALLOCATE(nr3p(nproc3))
            ALLOCATE(points(3,nnr))
            !
            i0r2p = environ%setup%system_cell%dfft%i0r2p
            i0r3p = environ%setup%system_cell%dfft%i0r3p
            nr2p = environ%setup%system_cell%dfft%nr2p
            nr3p = environ%setup%system_cell%dfft%nr3p
            points = environ%setup%get_coords(nnr)
            !
            !----------------------------------------------------------------------------
            ! Write data for aims
            !
            WRITE(fname, '(A,I0.10)') 'env_to_aims_grid_info_id', myid
            OPEN(UNIT=147, FILE=fname, STATUS='unknown', ACTION='write', &
                 FORM='unformatted', ACCESS='stream')
            !
            WRITE(147) nnr
            WRITE(147) ir_end
            WRITE(147) nnt
            WRITE(147) nr
            WRITE(147) nrx
            WRITE(147) nproc2
            WRITE(147) nproc3
            WRITE(147) i0r2p
            WRITE(147) nr2p
            WRITE(147) i0r3p
            WRITE(147) nr3p
            WRITE(147) points
            !
            CLOSE(147)
            !
            !----------------------------------------------------------------------------
            !
            DEALLOCATE(points)
            DEALLOCATE(i0r2p)
            DEALLOCATE(i0r3p)
            DEALLOCATE(nr2p)
            DEALLOCATE(nr3p)
            !
            stat = stat_env_grid_getters
            !
            !----------------------------------------------------------------------------
        END SUBROUTINE environ_aims_grid_getters
        !--------------------------------------------------------------------------------
        !>
        !!
        !--------------------------------------------------------------------------------
        SUBROUTINE environ_aims_calculators(stat)
            !----------------------------------------------------------------------------
            !
            IMPLICIT NONE
            !
            INTEGER, INTENT(INOUT) :: stat
            !
            CHARACTER(34) :: fname
            !
            INTEGER :: ind, n_atoms, array_size, mype2, mype3
            LOGICAL :: forces_on
            REAL(DP) :: n_electrons, energy
            REAL(DP), ALLOCATABLE :: rho(:), gradrho(:,:), forces(:,:), dvtot(:)
            !
            !----------------------------------------------------------------------------
            !
            if (.not. initialized) then
                stat = stat_env_err_not_ini
                return
            endif
            if (.not. updated) then
                stat = stat_env_err_not_updt
                return
            endif
            !
            !----------------------------------------------------------------------------
            !
            energy = 0.
            !
            mype2 = environ%setup%system_cell%dfft%mype2
            mype3 = environ%setup%system_cell%dfft%mype3
            !
            !----------------------------------------------------------------------------
            ! Read data from aims
            !
            WRITE(fname, '(A,I0.10)') 'aims_to_env_elec_info_id', myid
            OPEN(UNIT=147, FILE=fname, STATUS='unknown', ACTION='read', &
                 FORM='unformatted', ACCESS='stream')
            !
            READ(147) forces_on
            READ(147) n_atoms
            READ(147) n_electrons
            READ(147) array_size
            !
            ALLOCATE(rho(array_size))
            ALLOCATE(gradrho(3,array_size))
            !
            READ(147) rho
            READ(147) gradrho
            !
            CLOSE(147, STATUS='delete')
            !
            !----------------------------------------------------------------------------
            ! Update environ and deallocate temporary arrays
            !
            CALL environ%update_electrons(rho, nelec=n_electrons, lscatter=.false., &
                                          gradrho_in=gradrho)
            !
            DEALLOCATE(rho)
            DEALLOCATE(gradrho)
            !
            !----------------------------------------------------------------------------
            ! Run calculators
            !
            CALL environ%calc%potential(update=.true.)
            CALL environ%calc%energy(energy)
            !
            ALLOCATE(dvtot(array_size))
            !
            dvtot = environ%main%get_dvtot(array_size)
            !
            !----------------------------------------------------------------------------
            ! Write energy, as well as part of potential owned by this processor
            !
            WRITE(fname, '(A,I0.10)') 'env_to_aims_ener_info_id', myid
            OPEN(UNIT=147, FILE=fname, STATUS='unknown', ACTION='write', &
                 FORM='unformatted', ACCESS='stream')
            !
            WRITE(147) energy
            !
            CLOSE(147)
            !
            !----------------------------------------------------------------------------
            !
            WRITE(fname, '(A,I0.10,A,I0.10)') 'env_pot_proc_', mype2, '_', mype3
            OPEN(UNIT=147, FILE=fname, STATUS='unknown', ACTION='write', &
                 FORM='unformatted', ACCESS='stream')
            !
            WRITE(147) dvtot
            !
            CLOSE(147)
            !
            DEALLOCATE(dvtot)
            !
            !----------------------------------------------------------------------------
            ! Calculate and write forces, if requested
            !
            IF (forces_on) THEN
                ALLOCATE(forces(3,n_atoms))
                CALL environ%calc%force(n_atoms, forces)
                !
                WRITE(fname, '(A,I0.10)') 'env_to_aims_forc_info_id', myid
                OPEN(UNIT=147, FILE=fname, STATUS='unknown', ACTION='write', &
                     FORM='unformatted', ACCESS='stream')
                !
                WRITE(147) forces
                !
                CLOSE(147)
                !
                DEALLOCATE(forces)
            ENDIF
            !
            stat = stat_env_calc
            !
            !----------------------------------------------------------------------------
        END SUBROUTINE environ_aims_calculators
        !--------------------------------------------------------------------------------
        !>
        !!
        !--------------------------------------------------------------------------------
        SUBROUTINE environ_aims_cleanup(stat)
            !----------------------------------------------------------------------------
            !
            IMPLICIT NONE
            !
            INTEGER, INTENT(INOUT) :: stat
            !
            INTEGER :: mype2, mype3
            CHARACTER(34) :: fname
            !
            !----------------------------------------------------------------------------
            ! Clean up potential files
            !
            mype2 = environ%setup%system_cell%dfft%mype2
            mype3 = environ%setup%system_cell%dfft%mype3
            !
            WRITE(fname, '(A,I0.10,A,I0.10)') 'env_pot_proc_', mype2, '_', mype3
            OPEN(UNIT=147, FILE=fname, STATUS='unknown', ACTION='write', &
                 FORM='unformatted', ACCESS='stream')
            !
            CLOSE(147, STATUS='delete')
            !
            !----------------------------------------------------------------------------
            ! clean_up() will be called in driver.f90, we're done here.
            !
            stat = stat_env_cleanup
            !
            !----------------------------------------------------------------------------
        END SUBROUTINE environ_aims_cleanup
        !--------------------------------------------------------------------------------
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE run_environ_with_aims
    !------------------------------------------------------------------------------------
    !>
    !! An Environ calculation on a "frozen" density provided in a cube file
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE run_environ_from_cube()
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP), ALLOCATABLE :: rho(:)
        !
        REAL(DP), ALLOCATABLE :: env_potential(:)
        REAL(DP), ALLOCATABLE :: env_force(:, :)
        REAL(DP) :: env_energy
        !
        REAL(DP) :: volume, avg_dvtot, avg_velectrostatic
        !
        INTEGER :: i, j, nat
        !
        REAL(DP) :: gcutm, tmp, a1(3), sumat2, est
        !
        CHARACTER(LEN=80) :: routine = 'run_environ_from_cube'
        !
        !--------------------------------------------------------------------------------
        ! Initialize Environ
        !
        IF (no_density) THEN
            CALL init_environ_from_cube(environ)
        ELSE
            !
            CALL init_environ_from_cube(environ, rho)
            !
            CALL environ%update_electrons(rho, lscatter=.TRUE.)
            !
        END IF
        !
        !--------------------------------------------------------------------------------
        ! Compute potential
        !
        env_energy = 0.D0
        !
        nat = environ%main%system_ions%number
        !
        ALLOCATE (env_potential(environ%setup%get_nnt()))
        ALLOCATE (env_force(3, nat))
        !
        CALL environ%calc_potential(.TRUE., env_potential, lgather=.TRUE.)
        !
        CALL environ%calc%energy(env_energy)
        !
        IF (calc_force) CALL environ%calc%force(nat, env_force)
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
            WRITE (io%unit, 1000), environ%main%system_charges%charge
            !
            IF (environ%setup%has_electrostatics()) &
                WRITE (io%unit, 1001), avg_velectrostatic
            !
            WRITE (io%unit, 1002), avg_dvtot
            !
            CALL environ%main%print_energies('PW', .FALSE.)
            !
            WRITE (io%unit, 1003), env_energy
            !
            IF (calc_force) THEN
                !
                WRITE (io%unit, 1004), SUM(env_force)
                !
                DO i = 1, nat
                    WRITE (io%unit, 1005) i, (env_force(j, i), j=1, 3)
                END DO
                !
            END IF
            !
            WRITE (io%unit, *) ! final blank line
        END IF
        !
        !--------------------------------------------------------------------------------
        !
1000    FORMAT(5X, "total charge              =", F17.8, " a.u.",/)
!
1001    FORMAT(5X, "electrostatic potential   =", F17.8, " Ry")
        !
1002    FORMAT(5X, "total potential           =", F17.8, " Ry",/)
        !
1003    FORMAT(5X, "total energy              =", F17.8, " Ry",/)
!
1004    FORMAT(5X, "total force               =", F17.8, " Ry/bohr",/)
        !
1005    FORMAT(5X, "force on atom ", I4, "        =", 3F17.8)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE run_environ_from_cube
    !------------------------------------------------------------------------------------
    !>
    !! An Environ calculation of the Soft-Sphere interface function surface,
    !! volume, and local surface and volume. Additionally, the forces due to
    !! the interface function and energy contributions can be calculated.
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE run_descriptors_generator()
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        LOGICAL :: reduce_cell = .FALSE.
        !
        TYPE(environ_cell), POINTER :: cell
        !
        ! TODO: If we want to integrate over cylinders or slabs, we should make these
        !       input parameters or command line arguments. For now, integrate over
        !       spheres.
        INTEGER, PARAMETER :: descriptor_dim = 0
        INTEGER, PARAMETER :: descriptor_axis = 1
        !
        CHARACTER(LEN=80) :: routine = 'run_descriptors_generator'
        !
        !--------------------------------------------------------------------------------
        ! Validate input parameters
        !
        CALL check_input()
        !
        !--------------------------------------------------------------------------------
        !
        IF (.NOT. calc_energy) reduce_cell = .TRUE.
        !
        CALL init_environ_from_cube(environ, reduce_cell=reduce_cell)
        !
        !--------------------------------------------------------------------------------
        !
        cell => environ%setup%environment_cell
        !
        !--------------------------------------------------------------------------------
        ! Get the energy of the system
        !
        IF (calc_energy) CALL get_energy()
        !
        !--------------------------------------------------------------------------------
        ! Get the forces due to the surface and volume
        !
        IF (calc_force) CALL get_forces()
        !
        !--------------------------------------------------------------------------------
        ! Compute descriptors
        !
        CALL calc_descriptors(descriptor_dim, descriptor_axis)
        !
        !--------------------------------------------------------------------------------
        ! Write cube files corresponding to descriptors
        !
        CALL write_cube_descriptors()
        !
        !--------------------------------------------------------------------------------
    CONTAINS
        !--------------------------------------------------------------------------------
        !>
        !!
        !--------------------------------------------------------------------------------
        SUBROUTINE check_input()
            !----------------------------------------------------------------------------
            !
            IF (alpha_max + alpha_min + alpha_step /= -3.D0) THEN
                !
                IF (alpha_max == -1.D0) &
                    CALL io%error(routine, 'Missing maximum alpha value', 1)
                !
                IF (alpha_step == -1.D0) THEN
                    !
                    CALL io%warning("No alpha step value given. Using 1.0 as step value", 1)
                    !
                    alpha_step = 1.D0
                END IF
                !
                IF (alpha_min == -1.D0) THEN
                    !
                    CALL io%warning("No minimum alpha value given. Using step value as minimum", 1)
                    !
                    alpha_min = alpha_step
                END IF
                !
            END IF
            !
            !----------------------------------------------------------------------------
        END SUBROUTINE check_input
        !--------------------------------------------------------------------------------
        !>
        !!
        !--------------------------------------------------------------------------------
        SUBROUTINE calc_descriptors(dim, axis)
            !----------------------------------------------------------------------------
            !
            IMPLICIT NONE
            !
            INTEGER, INTENT(IN) :: dim, axis
            !
            INTEGER :: i, j, ir
            REAL(DP) :: r(3), r2, dist, weight
            LOGICAL :: physical
            REAL(DP), ALLOCATABLE :: alpha(:), pvol(:, :), psurf(:, :), ppol(:, :)
            !
            TYPE(environ_boundary_ionic), POINTER :: solvent
            !
            !----------------------------------------------------------------------------
            !
            SELECT TYPE (main_solvent => environ%main%solvent)
                !
            TYPE IS (environ_boundary_ionic)
                solvent => main_solvent
                !
            END SELECT
            !
            !----------------------------------------------------------------------------
            !
            ASSOCIATE (ions => environ%main%system_ions, &
                       scal => solvent%scaled, &
                       mod_grad => solvent%gradient%modulus, &
                       num => environ%main%system_ions%number, &
                       a_num => NINT(alpha_max / alpha_step), &
                       vol => solvent%volume, &
                       surf => solvent%surface, &
                       pol_den => environ%main%static%density)
                !
                !------------------------------------------------------------------------
                ! Print out surface and volume
                !
                IF (io%lnode) THEN
                    WRITE (io%unit, 2000) vol
                    WRITE (io%unit, 2001) surf
                END IF
                !
                IF (alpha_max + alpha_min + alpha_step == -3.D0) RETURN
                !
                !------------------------------------------------------------------------
                !
                ALLOCATE (alpha(a_num))
                ALLOCATE (pvol(num, a_num))
                ALLOCATE (psurf(num, a_num))
                ALLOCATE (ppol(num, a_num))
                !
                pvol = 0.D0
                psurf = 0.D0
                ppol = 0.D0
                !
                DO i = 1, a_num
                    alpha(i) = alpha_step * i
                END DO
                !
                !------------------------------------------------------------------------
                !
                DO ir = 1, cell%nnr
                    !
                    DO i = 1, num
                        !
                        CALL cell%get_min_distance(ir, dim, axis, &
                                                   ions%tau(:,i), r, r2, physical)
                        !
                        IF (.NOT. physical) CYCLE
                        !
                        dist = SQRT(r2)
                        !
                        DO j = 1, a_num
                            !
                            weight = 0.5D0 * environ_erfc(2*(dist-alpha(j))/alpha_step) &
                                     * cell%domega
                            pvol(i, j) = pvol(i, j) + scal%of_r(ir) * weight
                            psurf(i, j) = psurf(i, j) + mod_grad%of_r(ir) * weight
                            if (environ%main%setup%lelectrostatic) ppol(i, j) = &
                                ppol(i, j) + pol_den%of_r(ir) * weight
                            !
                        END DO
                        !
                    END DO
                    !
                END DO
                !
#if defined (__MPI)
                CALL env_mp_sum(pvol, io%comm)
                !
                CALL env_mp_sum(psurf, io%comm)
                !
                CALL env_mp_sum(ppol, io%comm)
#endif
                !
                IF (io%lnode) THEN
                    !
                    DO i = 1, num
                        WRITE (io%unit, 2002) i, ions%iontype(ions%ityp(i))%atmnum
                        !
                        WRITE (io%unit, 2008) "Bohr", ions%tau(:,i)
                        !
                        WRITE (io%unit, 2008) "Angstrom", ions%tau(:,i) * BOHR_RADIUS_ANGS
                        !
                        WRITE (io%unit, '(A)') ""
                        !
                        WRITE (io%unit, 2003)
                        !
                        DO j = 1, a_num
                            WRITE (io%unit, 2004) alpha(j), pvol(i, j), psurf(i, j), ppol(i, j)
                        END DO
                        !
                        WRITE (io%unit, *)
                    END DO
                    !
                END IF
                !
            END ASSOCIATE
            !
            FLUSH (io%unit)
            !
            !----------------------------------------------------------------------------
            !
2000        FORMAT(5X, "Total volume of the QM region: ", F18.8)
            !
2001        FORMAT(5X, "Total surface of the QM region: ", F17.8,/)
            !
2002        FORMAT(5X, "Atom number: ", I4, "; Atom type: ", I3,/)
            !
2008        FORMAT(5X, "Position (", A8, "): ", F18.8, 1X, F18.8, 1X, F18.8)
            !
2003        FORMAT(10X, "alpha  |    Partial Volume    |     Partial Surface  |     Partial Pol. Density", /, 9X, 81('-'))
            !
2004        FORMAT(10X, F6.3, ' |', F17.8, '     |', f17.8, '     |', f17.8)
            !
            !----------------------------------------------------------------------------
        END SUBROUTINE calc_descriptors
        !--------------------------------------------------------------------------------
        !>
        !!
        !--------------------------------------------------------------------------------
        SUBROUTINE write_cube_descriptors()
            !----------------------------------------------------------------------------
            !
            IMPLICIT NONE
            !
            !----------------------------------------------------------------------------
            !
            CALL write_cube(environ%main%solvent%scaled, &
                            environ%main%system_ions, 3)
            CALL write_cube(environ%main%solvent%gradient%modulus, &
                            environ%main%system_ions, 3)
            IF (environ%main%setup%lelectrostatic) &
                CALL write_cube(environ%main%static%density, &
                                environ%main%system_ions, 3, label='static_density')
            !
            !----------------------------------------------------------------------------
        END SUBROUTINE write_cube_descriptors
        !--------------------------------------------------------------------------------
        !>
        !!
        !--------------------------------------------------------------------------------
        SUBROUTINE get_energy()
            !----------------------------------------------------------------------------
            !
            IMPLICIT NONE
            !
            REAL(DP), ALLOCATABLE :: env_potential(:)
            REAL(DP) :: env_energy
            !
            !----------------------------------------------------------------------------
            !
            env_energy = 0.D0
            !
            ALLOCATE (env_potential(environ%setup%get_nnt()))
            !
            CALL environ%calc_potential(.TRUE., env_potential, lgather=.TRUE.)
            !
            CALL environ%calc%energy(env_energy)
            !
            !----------------------------------------------------------------------------
            !
            IF (io%lnode) THEN
                !
                CALL environ%main%print_energies('PW', .FALSE.)
                !
                WRITE (io%unit, 2005) env_energy
            END IF
            !
            FLUSH (io%unit)
            !
            !----------------------------------------------------------------------------
            !
2005        FORMAT(/, 5X, "total energy ", F17.8,/)
            !
            !----------------------------------------------------------------------------
        END SUBROUTINE get_energy
        !--------------------------------------------------------------------------------
        !>
        !!
        !--------------------------------------------------------------------------------
        SUBROUTINE get_forces()
            !----------------------------------------------------------------------------
            !
            IMPLICIT NONE
            !
            INTEGER :: k
            !
            INTEGER, POINTER :: nat
            !
            REAL(DP), DIMENSION(:, :), ALLOCATABLE :: surf_force, vol_force
            !
            TYPE(environ_density) :: surf_den, vol_den
            TYPE(environ_gradient) :: grad
            !
            !----------------------------------------------------------------------------
            !
            nat => environ%main%environment_ions%number
            !
            ALLOCATE (surf_force(3, nat), vol_force(3, nat))
            !
            surf_force = 0.D0
            vol_force = 0.D0
            !
            !----------------------------------------------------------------------------
            ! Initialize densities
            !
            CALL surf_den%init(cell)
            !
            CALL vol_den%init(cell)
            !
            CALL grad%init(cell)
            !
            !----------------------------------------------------------------------------
            ! Calculate forces
            !
            CALL environ%main%solvent%desurface_dboundary(1.D0, surf_den)
            !
            CALL environ%main%solvent%devolume_dboundary(1.D0, vol_den)
            !
            DO k = 1, nat
                !
                CALL environ%main%solvent%dboundary_dions(k, grad)
                !
                surf_force(:, k) = surf_force(:, k) - grad%scalar_product_density(surf_den)
                vol_force(:, k) = vol_force(:, k) - grad%scalar_product_density(vol_den)
            END DO
            !
            IF (io%lnode) THEN
                !
                WRITE (io%unit, 2006) &
                    'idx', &
                    'S_FORCE_x', 'S_FORCE_y', 'S_FORCE_z', &
                    'V_FORCE_x', 'V_FORCE_y', 'V_FORCE_z'
                !
                DO k = 1, nat
                    WRITE (io%unit, 2007) k, surf_force(:, k), vol_force(:, k)
                END DO
                !
                WRITE (io%unit, *)
            END IF
            !
            CALL grad%destroy()
            !
            CALL surf_den%destroy()
            !
            CALL vol_den%destroy()
            !
            FLUSH (io%unit)
            !
            !----------------------------------------------------------------------------
            !
2006        FORMAT(5X, A3, 7X, A9, 8X, A9, 8X, A9, 8X, A9, 8X, A9, 8X, A9)
            !
2007        FORMAT(5X, I3, 3F17.8, 3F17.8)
            !
            !----------------------------------------------------------------------------
        END SUBROUTINE get_forces
        !--------------------------------------------------------------------------------
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE run_descriptors_generator
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !                                   SETUP ROUTINES
    !
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE initial_setup()
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER :: comm
        LOGICAL :: lnode
        INTEGER :: ionode = 0
        !
        !--------------------------------------------------------------------------------
        !
        comm = get_comm()
        lnode = env_mp_rank(comm) == ionode
        !
        CALL environ%init_interface()
        !
        CALL environ%init_io(lnode, ionode, comm, 6, .FALSE.)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE initial_setup
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE clean_up()
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER :: ierr
        !
        CHARACTER(LEN=80) :: routine = 'clean_up'
        !
        !--------------------------------------------------------------------------------
        !
        IF (environ%main%initialized) CALL environ%destroy()
        !
#if defined(__MPI)
        CALL MPI_Finalize(ierr)
        !
        IF (ierr /= 0) CALL env_mp_stop(8001)
#endif
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE clean_up
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
    INTEGER FUNCTION get_comm() RESULT(comm)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER :: ierr
        !
        !--------------------------------------------------------------------------------
        !
#if defined(__MPI)
        comm = MPI_COMM_WORLD
        !
        CALL MPI_Init(ierr)
        !
        IF (ierr /= 0) CALL env_mp_stop(8000)
#else
        comm = 0
#endif
        !
        !--------------------------------------------------------------------------------
    END FUNCTION get_comm
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE programs
!----------------------------------------------------------------------------------------
