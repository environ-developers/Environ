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
! Authors: Edan Bainglass   (Department of Physics, UNT)
!          Matthew Truscott (Department of Physics, UNT)
!
!----------------------------------------------------------------------------------------
!>
!!
!----------------------------------------------------------------------------------------
MODULE environ_api
    !------------------------------------------------------------------------------------
    !
    USE class_io, ONLY: io
    !
    USE env_base_scatter, ONLY: env_scatter_grid, env_gather_grid
    !
    USE environ_param, ONLY: DP
    !
    USE class_calculator
    USE class_destructor
    USE class_environ
    USE class_setup
    !
    USE environ_input, ONLY: read_environ_input
    !
    !------------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    PRIVATE
    !
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    TYPE :: environ_interface
        !--------------------------------------------------------------------------------
        !
        TYPE(environ_setup) :: setup
        TYPE(environ_main) :: main
        TYPE(environ_calculator) :: calc
        TYPE(environ_destructor) :: clean
        !
        !--------------------------------------------------------------------------------
    CONTAINS
        !--------------------------------------------------------------------------------
        !
        PROCEDURE :: init => init_interface
        !
        PROCEDURE, NOPASS :: init_io
        PROCEDURE, NOPASS :: read_input
        !
        PROCEDURE :: update_cell
        PROCEDURE :: update_electrons
        !
        PROCEDURE :: add_charges
        !
        PROCEDURE, NOPASS :: get_verbosity
        !
        PROCEDURE :: calc_potential
        !
        !--------------------------------------------------------------------------------
    END TYPE environ_interface
    !------------------------------------------------------------------------------------
    !
    TYPE(environ_interface), PUBLIC :: environ
    !
    !------------------------------------------------------------------------------------
CONTAINS
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !                                   ADMIN METHODS
    !
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_interface(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_interface), TARGET, INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: sub_name = 'init_interface'
        !
        !--------------------------------------------------------------------------------
        !
        this%calc%main => this%main
        this%clean%main => this%main
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_interface
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_io(ionode, ionode_id, comm, program_unit, lstdout)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        LOGICAL, INTENT(IN) :: ionode
        INTEGER, INTENT(IN) :: ionode_id
        INTEGER, INTENT(IN) :: comm
        INTEGER, INTENT(IN) :: program_unit
        LOGICAL, INTENT(IN) :: lstdout
        !
        CHARACTER(LEN=80) :: sub_name = 'init_io'
        !
        !--------------------------------------------------------------------------------
        !
        CALL io%init(ionode, ionode_id, comm, program_unit, lstdout)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_io
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE read_input(filename, nsx)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: filename
        INTEGER, OPTIONAL, INTENT(IN) :: nsx
        !
        CHARACTER(LEN=80) :: sub_name = 'read_input'
        !
        !--------------------------------------------------------------------------------
        !
        CALL read_environ_input(filename, nsx)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE read_input
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !                                   UPDATE METHODS
    !
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE update_cell(this, at)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP), INTENT(IN) :: at(3, 3)
        !
        CLASS(environ_interface), INTENT(INOUT) :: this
        !
        !--------------------------------------------------------------------------------
        !
        CALL this%setup%update_cell(at)
        !
        CALL this%main%update_cell_dependent_quantities()
        !
        CALL this%setup%end_cell_update()
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE update_cell
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE update_electrons(this, rho, lscatter)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP), INTENT(IN) :: rho(:)
        LOGICAL, OPTIONAL, INTENT(IN) :: lscatter
        !
        CLASS(environ_interface), INTENT(INOUT) :: this
        !
        REAL(DP) :: aux(this%setup%system_cell%dfft%nnr)
        REAL(DP) :: nelec
        !
        !--------------------------------------------------------------------------------
        !
#if defined(__MPI)
        IF (PRESENT(lscatter)) THEN
            !
            IF (lscatter) THEN
                CALL env_scatter_grid(this%setup%system_cell%dfft, rho, aux)
            ELSE
                aux = rho
            END IF
            !
        ELSE
            aux = rho
        END IF
        !
#else
        aux = rho
#endif
        nelec = REAL(this%main%system_electrons%number, DP)
        !
        CALL this%main%update_electrons(this%setup%system_cell%dfft%nnr, aux, nelec)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE update_electrons
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE add_charges(this, rho, lscatter)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP), INTENT(IN) :: rho(:)
        LOGICAL, OPTIONAL, INTENT(IN) :: lscatter
        !
        CLASS(environ_interface), INTENT(INOUT) :: this
        !
        REAL(DP) :: aux(this%setup%system_cell%dfft%nnr)
        !
        CHARACTER(LEN=80) :: local_label = 'mbx_charges'
        !
        !--------------------------------------------------------------------------------
        !
#if defined(__MPI)
        IF (PRESENT(lscatter)) THEN
            !
            IF (lscatter) THEN
                CALL env_scatter_grid(this%setup%system_cell%dfft, rho, aux)
            ELSE
                aux = rho
            END IF
            !
        ELSE
            aux = rho
        END IF
        !
#else
        aux = rho
        !
#endif
        CALL this%main%add_charges(this%setup%system_cell%dfft%nnr, aux, local_label)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE add_charges
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !                                   ACCESS METHODS
    !
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    INTEGER FUNCTION get_verbosity()
        !--------------------------------------------------------------------------------
        !
        get_verbosity = io%verbosity
        !
        !--------------------------------------------------------------------------------
    END FUNCTION get_verbosity
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !                                COMPUTATION METHODS
    !
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_potential(this, update, potential, local_verbose, lgather)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        LOGICAL, INTENT(IN) :: update
        INTEGER, OPTIONAL, INTENT(IN) :: local_verbose
        LOGICAL, OPTIONAL, INTENT(IN) :: lgather
        !
        CLASS(environ_interface), INTENT(INOUT) :: this
        !
        REAL(DP), INTENT(OUT) :: potential(this%setup%system_cell%dfft%nnt)
        !
        REAL(DP) :: aux(this%setup%system_cell%dfft%nnr)
        !
        !--------------------------------------------------------------------------------
        !
        CALL this%calc%potential(update, local_verbose)
        !
        aux = this%main%dvtot%of_r
        !
#if defined(__MPI)
        IF (PRESENT(lgather)) THEN
            !
            IF (lgather) THEN
                CALL env_gather_grid(this%setup%system_cell%dfft, aux, potential)
            ELSE
                potential = aux
            END IF
            !
        ELSE
            potential = aux
        END IF
#else
        !
        potential = aux
#endif
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_potential
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE environ_api
!----------------------------------------------------------------------------------------
