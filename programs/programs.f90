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
MODULE programs
    !------------------------------------------------------------------------------------
    !
    USE env_parallel_include
    USE env_mp, ONLY: env_mp_rank, env_mp_stop
    !
    USE env_mytime, ONLY: env_f_wall
    !
    USE class_io, ONLY: io
    !
    USE environ_param, ONLY: DP
    !
    USE environ_api, ONLY: environ_interface
    !
    USE class_density
    USE class_gradient
    !
    USE env_write_cube, ONLY: write_cube
    !
    USE cmdline_args
    !
    USE prog_utils
    !
    USE parsers
    !
    !------------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    PRIVATE
    !
    PUBLIC :: run_tester, run_environ_from_cube, run_desc
    !
    PUBLIC :: initial_setup, clean_up, print_available_programs
    !
    !------------------------------------------------------------------------------------
    ! Declare interface
    !
    TYPE(environ_interface) :: environ
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
            '- desc'
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
    !! An Environ calculation of the Soft-Sphere interface function surface,
    !! volume, and local surface and volume. Additionally, the forces due to
    !! the interface function and energy contributions can be calculated.
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE run_desc()
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        !------------------------------------------------------------------------------------
        ! Output file parameters and conversions
        !
        REAL(DP) :: time1, time2, ttime
        !
        CHARACTER(LEN=80) :: tlabel
        !
        !------------------------------------------------------------------------------------
        !
        !                                    MAIN PROGRAM
        !
        !------------------------------------------------------------------------------------
        ! Initialize interface, I/O, and get initial time
        !
        time1 = env_f_wall()
        !
        CALL init()
        !
        !------------------------------------------------------------------------------------
        ! Get the surface and volume descriptors
        !
        CALL get_desc()
        !
        !------------------------------------------------------------------------------------
        ! Get the forces due to the surface and volume
        !
        IF (calc_force) CALL get_forces()
        !
        !------------------------------------------------------------------------------------
        ! Get the energy of the system
        !
        IF (calc_energy) CALL get_energy()
        !
        !------------------------------------------------------------------------------------
        ! Get final time and run time of code
        !
        time2 = env_f_wall()
        IF (io%lnode) THEN
            !
            CALL get_time(time1, time2, tlabel)
            WRITE(io%unit,"(A,F17.8,1X,A7)") 'Total time: ', time2, tlabel
            FLUSH(io%unit)
            !
        END IF
        !
        !--------------------------------------------------------------------------------
        !
    CONTAINS
        !------------------------------------------------------------------------------------
        !>
        !!
        !------------------------------------------------------------------------------------
        SUBROUTINE init()
            !--------------------------------------------------------------------------------
            !
            IMPLICIT NONE
            !
            REAL(DP) :: at(3,3), new_at(3,3), origin(3)
            INTEGER :: nr(3), ntyp, nat
            REAL(DP), ALLOCATABLE :: zv(:)
            CHARACTER(LEN=2), ALLOCATABLE :: label(:)
            INTEGER, ALLOCATABLE :: ityp(:), atomic_number(:)
            REAL(DP), ALLOCATABLE :: tau(:, :)
            !
            !--------------------------------------------------------------------------------
            !
            io%unit = 6
            !
            CALL read_cube(nat, ntyp, ityp, label, zv, tau, origin, nr, at)
            !
            new_at = at
            IF (.NOT. calc_energy) THEN
                !
                CALL get_new_cell(at, tau, nat, new_at)
                !
                IF (io%lnode) THEN
                    !
                    WRITE(io%unit,"(5X,A)") 'Resizing cell for calculation of descriptors'
                    FLUSH(io%unit)
                    !
                END IF
                !
            END IF
            !
            CALL environ%read_input(inputfile,nsx=SIZE(label))
            !
            CALL environ%setup%init()
            !
            CALL environ%setup%init_cell(io%comm, new_at)
            !
            CALL environ%setup%init_numerical(use_internal_pbc_corr)
            !
            CALL environ%main%init(nat, ntyp, label, ityp, zv)
            !
            CALL environ%main%update_ions(nat, tau)
            !
            IF (io%lnode) THEN
                !
                WRITE(io%unit,"(5X,A)") 'Summary of system information'
                !
                WRITE(io%unit,"(10X,A)") 'Values will be printed using ENVIRON internal units, e.g. length is in Bohr'
                WRITE(io%unit,"(10X,A,3F17.8)") 'Cell parameters: ', new_at(1,1), new_at(2,2), new_at(3,3)
                WRITE(io%unit,"(10X,A,I4)") 'Number of atoms: ', nat
                WRITE(io%unit,"(10X,A,I4)") 'Number of atomic types: ', ntyp
                FLUSH(io%unit)
                !
            END IF
            !
            IF (alpha_max+alpha_min+alpha_step /= -3.D0) THEN
                !
                IF (alpha_max == -1.D0) CALL io%error('get_desc', 'Missing maximum alpha value', 1)
                !
                IF (alpha_step == -1.D0) THEN
                    !
                    IF (io%lnode) WRITE(io%unit,"(5X,A)") 'No alpha step value given. Using 1.0 as step value.'
                    alpha_step = 1.D0
                    !
                END IF
                !
                IF (alpha_min == -1.D0) THEN
                    !
                    IF (io%lnode) WRITE(io%unit,"(5X,A)") 'No minimum alpha value given. Using step value as minimum.'
                    alpha_min = alpha_step
                    !
                END IF
                !
            END IF
            !
        END SUBROUTINE init
        !--------------------------------------------------------------------------------
        !>
        !!
        !------------------------------------------------------------------------------------
        SUBROUTINE get_new_cell( at, tau, nat, new_at )
            !--------------------------------------------------------------------------------
            !
            IMPLICIT NONE
            !
            REAL(DP), INTENT(IN) :: at(3,3)
            REAL(DP), INTENT(INOUT) :: tau(:,:)
            INTEGER, INTENT(IN) :: nat
            REAL(DP), INTENT(OUT) :: new_at(3,3)
            !
            REAL(DP), PARAMETER :: fluff = 7.D0
            REAL(DP) :: min_vec(3), max_vec(3), shift(3)
            !
            INTEGER :: i
            !
            !--------------------------------------------------------------------------------
            !
            new_at = 0.D0
            !
            min_vec = MINVAL( tau, DIM=2 ) - fluff
            max_vec = MAXVAL( tau, DIM=2 ) + fluff
            shift = max_vec - min_vec
            !
            DO i=1,nat
                tau(:,i) = tau(:,i) - min_vec
            END DO
            !
            DO i=1,3
                new_at(i,i) = shift(i)
            END DO
            !
        END SUBROUTINE get_new_cell
        !------------------------------------------------------------------------------------
        !>
        !!
        !------------------------------------------------------------------------------------
        SUBROUTINE get_desc()
            !--------------------------------------------------------------------------------
            !
            IMPLICIT NONE
            !
            REAL(DP) :: volu, surf
            !
            REAL(DP) :: r(3), r2, dist
            REAL(DP), ALLOCATABLE :: alpha(:), pvolu(:,:), psurf(:,:)
            INTEGER :: a_num, num, nnr, i, j, ir
            LOGICAL :: physical
            !
            !--------------------------------------------------------------------------------
            ! Print out surface and volume
            !
            volu = environ%main%solvent%volume
            surf = environ%main%solvent%surface
            !
            IF (io%lnode) THEN
                !
                WRITE(io%unit,"(/,5X,A,F17.8)") 'Total Volume of the QM region: ', volu
                WRITE(io%unit,"(5X,A,F17.8,/)") 'Total Surface of the QM region: ', surf
                !
            END IF
            !
            IF (alpha_max+alpha_min+alpha_step == -3.D0) RETURN
            !
            num = environ%main%solvent%soft_spheres%number
            nnr = environ%setup%environment_cell%nnr
            !
            !--------------------------------------------------------------------------------
            ! Beginning of local surface and volume calculations
            !
            a_num = NINT(alpha_max/alpha_step)
            ALLOCATE(alpha(a_num))
            ALLOCATE(pvolu(num,a_num))
            ALLOCATE(psurf(num,a_num))
            !
            pvolu = 0.D0
            psurf = 0.D0
            !
            DO i=1,a_num
                alpha(i) = alpha_step*i
            END DO
            !
            ASSOCIATE(solvent => environ%main%solvent, &
                ss => environ%main%solvent%soft_spheres%array, &
                scal => environ%main%solvent%scaled, &
                mod_grad => environ%main%solvent%gradient%modulus)
                !
                DO ir=1,nnr
                    !
                    DO i=1,solvent%soft_spheres%number
                        !
                        ASSOCIATE (cell => environ%setup%environment_cell, &
                            pos => ss(i)%pos, &
                            spread => ss(i)%spread, &
                            width => ss(i)%width, &
                            dim => ss(i)%dim, &
                            axis => ss(i)%axis)
                            !
                            CALL cell%get_min_distance(ir, dim, axis, pos, r, r2, physical)
                            !
                            IF (.NOT. physical) CYCLE
                            dist = SQRT(r2)
                            !
                            DO j=1,a_num
                                !
                                IF (dist <= alpha(j)) THEN
                                    !
                                    pvolu(i,j) = pvolu(i,j) + scal%of_r(ir)*cell%domega
                                    psurf(i,j) = psurf(i,j) + mod_grad%of_r(ir)*cell%domega
                                    !
                                END IF
                                !
                            END DO
                            !
                        END ASSOCIATE
                        !
                    END DO
                    !
                END DO
                !
            END ASSOCIATE
            !
            IF (io%lnode) THEN
                !
                DO i=1,num
                    !
                    WRITE(io%unit,"(5X,A,I4,A,I4,/)") 'Atom number: ', i, &
                        '; Atom type: ', environ%main%solvent%ions%ityp(i)
                    !
                    WRITE(io%unit,"(10X,A,4X,A,5X,A)") 'alpha  |', 'Partial Volume    |', 'Partial Surface'
                    WRITE(io%unit,"(8X,52('-'))")
                    !
                    DO j=1,a_num
                        !
                        WRITE(io%unit,"(10X,F6.3,' |',F17.8,'     |',f17.8)") alpha(j), pvolu(i,j), psurf(i,j)
                        !
                    END DO
                    !
                    WRITE(io%unit,*) ''
                    !
                END DO
                !
                FLUSH(io%unit)
                !
            END IF
            !
            !--------------------------------------------------------------------------------
        END SUBROUTINE get_desc
        !------------------------------------------------------------------------------------
        !>
        !!
        !------------------------------------------------------------------------------------
        SUBROUTINE get_forces()
            !--------------------------------------------------------------------------------
            !
            IMPLICIT NONE
            !
            REAL(DP), ALLOCATABLE :: surf_for(:,:), vol_for(:,:)
            Integer :: k, nat
            !
            TYPE(environ_density) :: surf_den, vol_den
            TYPE(environ_gradient) :: grad
            !
            nat = environ%main%environment_ions%number
            ALLOCATE( surf_for(3,nat), vol_for(3,nat) )
            surf_for = 0.D0
            vol_for = 0.D0
            !
            !Initialize densities
            CALL surf_den%init( environ%setup%environment_cell )
            CALL vol_den%init( environ%setup%environment_cell )
            CALL grad%init( environ%setup%environment_cell )
            !
            !Calculate forces
            CALL environ%main%solvent%desurface_dboundary(1.D0,surf_den)
            CALL environ%main%solvent%devolume_dboundary(1.D0,vol_den)
            !
            DO k=1,nat
                !
                CALL environ%main%solvent%dboundary_dions(k,grad)
                !
                surf_for(:,k) = surf_for(:,k) - grad%scalar_product_density(surf_den)
                !
                vol_for(:,k) = vol_for(:,k) - grad%scalar_product_density(vol_den)
                !
            END DO
            !
            IF (io%lnode) THEN
                !
                WRITE(io%unit, '(5X,A3,7X,A9,8X,A9,8X,A9,8X,A9,8X,A9,8X,A9)') 'idx',&
                    &'S_FORCE_x','S_FORCE_y','S_FORCE_z','V_FORCE_x','V_FORCE_y',&
                    &'V_FORCE_z'
                !
                DO k=1,nat
                    !
                    WRITE(io%unit, '(5X,I3,3F17.8,3F17.8)') k, surf_for(:,k), vol_for(:,k)
                    !
                END DO
                !
                WRITE(io%unit,*) ''
                FLUSH(io%unit)
                !
            END IF
            !
            CALL grad%destroy()
            CALL surf_den%destroy()
            CALL vol_den%destroy()
            !
            !--------------------------------------------------------------------------------
        END SUBROUTINE get_forces
        !------------------------------------------------------------------------------------
        !>
        !!
        !------------------------------------------------------------------------------------
        SUBROUTINE get_energy()
            !--------------------------------------------------------------------------------
            !
            IMPLICIT NONE
            !
            REAL(DP), ALLOCATABLE :: env_potential(:)
            REAL(DP) :: env_energy = 0.D0
            !
            !------------------------------------------------------------------------------------
            !
            ALLOCATE(env_potential(environ%setup%get_nnt()))
            !
            CALL environ%calc_potential(.TRUE., env_potential, lgather=.TRUE.)
            !
            CALL environ%calc%energy(env_energy)
            !
            !------------------------------------------------------------------------------------
            !
            IF (io%lnode) THEN
                !
                CALL environ%main%print_energies('PW', .FALSE.)
                !
                WRITE(io%unit,"(/,5X,A,F17.8,/)") 'total energy ', env_energy
                FLUSH(io%unit)
                !
            END IF
            !
            !--------------------------------------------------------------------------------
        END SUBROUTINE get_energy
        !------------------------------------------------------------------------------------
        !>
        !!
        !------------------------------------------------------------------------------------
        SUBROUTINE get_time(t1, t2, label)
            !--------------------------------------------------------------------------------
            !
            IMPLICIT NONE
            !
            REAL(DP), INTENT(IN) :: t1
            REAL(DP), INTENT(INOUT) :: t2
            CHARACTER(LEN=*), INTENT(OUT) :: label
            !
            !------------------------------------------------------------------------------------
            !
            t2 = t2 - t1
            label = 'seconds'
            !
            IF (t2>60.D0 .AND. t2<3600.D0) THEN
                !
                t2 = t2/60.D0
                label = 'minutes'
                !
            ELSE IF (t2>3600.D0) THEN
                !
                t2 = t2/3600.D0
                label = 'hours'
                !
            END IF
            !
            !--------------------------------------------------------------------------------
        END SUBROUTINE get_time
        !------------------------------------------------------------------------------------
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE run_desc
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
        CHARACTER(LEN=80) :: sub_name = 'run_environ_from_cube'
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
        CHARACTER(LEN=80) :: sub_name = 'clean_up'
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
