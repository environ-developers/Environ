!----------------------------------------------------------------------------------------
!
! Copyright (C) 2021 ENVIRON (www.quantum-environment.org)
!
!----------------------------------------------------------------------------------------
!
! This file is part of Environ version 2.0
!
! Environ 2.0 is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 2 of the License, or
! (at your option) any later version.
!
! Environ 2.0 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more detail, either the file
! `License' in the root directory of the present distribution, or
! online at <http://www.gnu.org/licenses/>.
!
!----------------------------------------------------------------------------------------
!
! Authors: Compiled by Edan Bainglass (Department of Physics, UNT)
!
!----------------------------------------------------------------------------------------
!>
!!
!----------------------------------------------------------------------------------------
MODULE env_base_io
    !------------------------------------------------------------------------------------
    !
    USE env_char_ops, ONLY: env_uppercase
    USE env_utils_param
    USE env_mp, ONLY: env_mp_bcast, env_mp_abort, env_mp_rank
    !
#if defined(__PTRACE)&&defined(__INTEL_COMPILER)
    USE ifcore, ONLY: tracebackqq
#endif
    !
    IMPLICIT NONE
    !
    PRIVATE
    !
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    TYPE environ_io
        !--------------------------------------------------------------------------------
        !
        LOGICAL :: lnode = .TRUE.
        INTEGER :: node = 0
        !
        LOGICAL :: lstdout = .FALSE. ! whether environ can print on standard output
        !
        INTEGER :: comm ! WE MAY NEED A SECOND COMMUNICATOR FOR IMAGE PARALLELIZATION
        !
        INTEGER :: unit
        INTEGER :: debug_unit
        INTEGER :: verbosity
        !
        CHARACTER(LEN=2) :: prog ! the calling program
        !
        !--------------------------------------------------------------------------------
    CONTAINS
        !--------------------------------------------------------------------------------
        !
        PROCEDURE :: init => init_base_io
        PROCEDURE :: update_unit => update_output_program_unit
        !
        PROCEDURE, NOPASS :: find_free_unit => env_find_free_unit
        !
        !--------------------------------------------------------------------------------
    END TYPE environ_io
    !------------------------------------------------------------------------------------
    !
    TYPE(environ_io), PUBLIC, SAVE :: io
    !
    !------------------------------------------------------------------------------------
CONTAINS
    !------------------------------------------------------------------------------------
    !>
    !! Set global I/O constants
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_base_io(this, can_write, id, comm, prog_unit, lstdout)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_io), INTENT(INOUT) :: this
        !
        LOGICAL, INTENT(IN) :: can_write
        INTEGER, INTENT(IN) :: id
        INTEGER, INTENT(IN) :: comm
        INTEGER, INTENT(IN) :: prog_unit
        LOGICAL, INTENT(IN) :: lstdout
        !
        !--------------------------------------------------------------------------------
        !
        this%lnode = can_write
        this%node = id
        this%comm = comm
        this%unit = prog_unit
        this%lstdout = lstdout
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_base_io
    !------------------------------------------------------------------------------------
    !>
    !! Sets the output file target
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE update_output_program_unit(this, program_unit_in)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: program_unit_in
        !
        CLASS(environ_io), INTENT(INOUT) :: this
        !
        !--------------------------------------------------------------------------------
        !
        this%unit = program_unit_in
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE update_output_program_unit
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    INTEGER FUNCTION env_find_free_unit()
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER :: iunit
        LOGICAL :: opnd
        !
        CHARACTER(LEN=80) :: fun_name = 'env_find_free_unit'
        !
        !--------------------------------------------------------------------------------
        !
        env_find_free_unit = -1
        !
        unit_loop: DO iunit = 99, 1, -1
            !
            INQUIRE (UNIT=iunit, OPENED=opnd)
            !
            IF (.NOT. opnd) THEN
                env_find_free_unit = iunit
                !
                RETURN
                !
            END IF
            !
        END DO unit_loop
        !
        CALL env_warning('free unit not found?!?')
        !
        !--------------------------------------------------------------------------------
    END FUNCTION env_find_free_unit
    !------------------------------------------------------------------------------------

    !
    !------------------------------------------------------------------------------------
END MODULE env_base_io
!----------------------------------------------------------------------------------------
