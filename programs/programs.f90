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
    USE class_io, ONLY: io
    !
    USE environ_param, ONLY: DP
    !
    USE environ_api, ONLY: environ_interface
    !
    USE cmdline_args
    !
    USE prog_utils
    !
    !------------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    PRIVATE
    !
    PUBLIC :: run_tester, run_environ_from_cube
    !
    PUBLIC :: general_setup, clean_up, print_available_programs
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
            PRINT '(3(/, A))', &
            'Available calculations:', &
            '- tester', &
            '- from_cube'
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
        REAL(DP) :: nelec
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
        ! Initialize Environ
        !
        CALL init_environ_from_cube(environ, cubefile, rho, nelec)
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
    !------------------------------------------------------------------------------------
    !
    !                                   SETUP ROUTINES
    !
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE general_setup()
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER :: comm
        LOGICAL :: lnode
        INTEGER :: ionode = 0
        !
        INTEGER, EXTERNAL :: env_mp_rank
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
        CALL environ%read_input(inputfile)
        !
        CALL environ%setup%init(use_pbc_corr)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE general_setup
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
