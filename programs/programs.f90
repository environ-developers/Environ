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
    USE env_mp, ONLY: env_mp_rank, env_mp_stop, env_mp_sum
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
    !------------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    PRIVATE
    !
    PUBLIC :: run_tester, run_environ_from_cube, run_descriptors_generator
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
            '- descriptors'
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
        CALL environ%calc%force(nat, env_force)
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
            WRITE (io%unit, 1004), SUM(env_force)
            !
            DO i = 1, nat
                WRITE (io%unit, 1005) i, (env_force(j, i), j=1, 3)
            END DO
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
            REAL(DP) :: r(3), r2, dist
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
                            IF (dist <= alpha(j)) THEN
                                pvol(i, j) = pvol(i, j) + scal%of_r(ir) * cell%domega
                                psurf(i, j) = psurf(i, j) + mod_grad%of_r(ir) * cell%domega
                                if (environ%main%setup%lelectrostatic) &
                                    ppol(i, j) = ppol(i, j) + pol_den%of_r(ir) * cell%domega
                            END IF
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
                        WRITE (io%unit, 2002) i, ions%ityp(i)
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
2002        FORMAT(5X, "Atom number: ", I4, "; Atom type: ", I4,/)
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
