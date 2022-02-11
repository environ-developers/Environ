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
! Authors: Oliviero Andreussi (Department of Physics, UNT)
!          Gabriel Medrano    (Department of Physics, UNT)
!          Edan Bainglass     (Department of Physics, UNT)
!
!----------------------------------------------------------------------------------------
!>
!!
!----------------------------------------------------------------------------------------
MODULE class_mapping
    !------------------------------------------------------------------------------------
    !
    USE class_io, ONLY: io
    USE env_mp, ONLY: env_mp_sum
    !
    USE env_scatter_mod, ONLY: env_scatter_grid, env_gather_grid
    !
    USE environ_param, ONLY: DP
    !
    USE class_cell
    USE class_density
    !
    !------------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    PRIVATE
    !
    !------------------------------------------------------------------------------------
    !>
    !! The mapping between system and environment simulation cells
    !!
    !------------------------------------------------------------------------------------
    TYPE, PUBLIC :: environ_mapping
        !--------------------------------------------------------------------------------
        !
        INTEGER :: nrep(3) = 0
        ! number of system cells in environment cell along a_i = 2 * nrep_i + 1
        !
        TYPE(environ_cell), POINTER :: small => NULL() ! system cell
        TYPE(environ_cell), POINTER :: large => NULL() ! environment cell
        !
        INTEGER, ALLOCATABLE :: map(:)
        !
        !--------------------------------------------------------------------------------
    CONTAINS
        !--------------------------------------------------------------------------------
        !
        PROCEDURE, PRIVATE :: create => create_environ_mapping
        PROCEDURE :: init => init_environ_mapping
        PROCEDURE :: update => update_environ_mapping
        PROCEDURE :: destroy => destroy_environ_mapping
        !
        PROCEDURE, PRIVATE :: map_small_to_large_real
        PROCEDURE, PRIVATE :: map_small_to_large_density
        PROCEDURE, PRIVATE :: map_large_to_small_real
        PROCEDURE, PRIVATE :: map_large_to_small_density
        !
        GENERIC :: to_large => map_small_to_large_real, map_small_to_large_density
        GENERIC :: to_small => map_large_to_small_real, map_large_to_small_density
        !
        PROCEDURE :: printout => print_environ_mapping
        !
        !--------------------------------------------------------------------------------
    END TYPE environ_mapping
    !------------------------------------------------------------------------------------
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
    SUBROUTINE create_environ_mapping(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_mapping), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: sub_name = 'create_environ_mapping'
        !
        !--------------------------------------------------------------------------------
        !
        IF (ASSOCIATED(this%large)) CALL io%create_error(sub_name)
        !
        IF (ASSOCIATED(this%small)) CALL io%create_error(sub_name)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE create_environ_mapping
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_environ_mapping(this, nrep, small, large)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: nrep(3)
        TYPE(environ_cell), TARGET, INTENT(IN) :: small, large
        !
        CLASS(environ_mapping), INTENT(INOUT) :: this
        !
        !--------------------------------------------------------------------------------
        !
        CALL this%create()
        !
        this%nrep = nrep
        !
        !--------------------------------------------------------------------------------
        ! Check that small%at and large%at are compatible with this%nrep
        !
        this%small => small
        this%large => large
        !
        IF (.NOT. ASSOCIATED(this%small, this%large)) THEN
            ALLOCATE (this%map(this%small%nnr))
            this%map = 0
        END IF
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_environ_mapping
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE update_environ_mapping(this, pos)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP), OPTIONAL, INTENT(IN) :: pos(3)
        !
        CLASS(environ_mapping), INTENT(INOUT) :: this
        !
        INTEGER :: ir
        LOGICAL :: physical
        INTEGER, DIMENSION(3) :: small_n, large_n, center, origin, shift, ijk
        REAL(DP) :: tmp(3)
        !
        !--------------------------------------------------------------------------------
        !
        IF (ASSOCIATED(this%small, this%large)) RETURN ! identical cells
        !
        !--------------------------------------------------------------------------------
        ! Set grid points
        !
        small_n(1) = this%small%dfft%nr1
        small_n(2) = this%small%dfft%nr2
        small_n(3) = this%small%dfft%nr3
        !
        large_n(1) = this%large%dfft%nr1
        large_n(2) = this%large%dfft%nr2
        large_n(3) = this%large%dfft%nr3
        !
        !--------------------------------------------------------------------------------
        ! Define origin of small cell (in units of nnr)
        !
        IF (PRESENT(pos)) THEN
            tmp = MATMUL(this%small%bg, pos)
            origin = NINT(tmp * small_n) ! center of charge
        ELSE
            origin = 0
        END IF
        !
        !--------------------------------------------------------------------------------
        ! Compute shift placing center of charge at center of small cell
        ! (Minimizes potential cutting of DFT densities)
        !
        center = NINT(small_n / 2.D0)
        shift = center - origin
        !
        !--------------------------------------------------------------------------------
        ! Shift large cell origin in internal length units
        !
        this%large%origin = -MATMUL(this%large%at, (0.5 - origin / DBLE(large_n)))
        !
        !--------------------------------------------------------------------------------
        ! Generate mapping
        !
        DO ir = 1, this%small%ir_end
            !
            CALL this%small%ir2ijk(ir, ijk(1), ijk(2), ijk(3), physical)
            !
            IF (.NOT. physical) CYCLE
            !
            !----------------------------------------------------------------------------
            ! Shift center of charge to center of small cell
            !
            ijk = ijk + shift
            ijk = ijk - FLOOR(DBLE(ijk) / small_n) * small_n ! enforce periodicity
            ijk = ijk + small_n * this%nrep
            !
            !----------------------------------------------------------------------------
            ! Map small cell to large cell
            !
            this%map(ir) = 1 + ijk(1) & ! x-point
                           & + ijk(2) * large_n(1) & ! y-row
                           & + ijk(3) * large_n(1) * large_n(2) ! z-plane
            !
        END DO
        !
        !--------------------------------------------------------------------------------
        ! Output current state
        !
        CALL this%printout()
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE update_environ_mapping
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE destroy_environ_mapping(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_mapping), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: sub_name = 'destroy_environ_mapping'
        !
        !--------------------------------------------------------------------------------
        !
        IF (.NOT. ASSOCIATED(this%small)) CALL io%destroy_error(sub_name)
        !
        IF (.NOT. ASSOCIATED(this%large)) CALL io%destroy_error(sub_name)
        !
        !--------------------------------------------------------------------------------
        !
        NULLIFY (this%small)
        NULLIFY (this%large)
        !
        IF (ALLOCATED(this%map)) DEALLOCATE (this%map)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE destroy_environ_mapping
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !                                  GENERAL METHODS
    !
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE map_small_to_large_real(this, nsmall, nlarge, fsmall, flarge)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_mapping), INTENT(IN) :: this
        INTEGER, INTENT(IN) :: nsmall, nlarge
        REAL(DP), INTENT(IN) :: fsmall(nsmall)
        !
        REAL(DP), INTENT(INOUT) :: flarge(nlarge)
        !
        INTEGER :: i
        REAL(DP), ALLOCATABLE :: auxlarge(:)
        !
        CHARACTER(LEN=80) :: sub_name = 'map_small_to_large_real'
        !
        !--------------------------------------------------------------------------------
        ! Check if input/output dimensions match mapping cells
        !
        IF (nsmall /= this%small%nnr) &
            CALL io%error(sub_name, "Wrong dimension of small cell", 1)
        !
        IF (nlarge /= this%large%nnr) &
            CALL io%error(sub_name, "Wrong dimension of large cell", 1)
        !
        !--------------------------------------------------------------------------------
        ! If the cells are the same, just copy
        !
        IF (nsmall == nlarge) THEN
            flarge = fsmall
        ELSE
            !
            !----------------------------------------------------------------------------
            ! Copy small cell to corresponding gridpoints in the full large cell
            !
            ALLOCATE (auxlarge(this%large%nnt))
            auxlarge = 0.D0
            !
            DO i = 1, this%small%ir_end
                IF (this%map(i) > 0) auxlarge(this%map(i)) = fsmall(i)
            END DO
            !
#if defined(__MPI)
            CALL env_mp_sum(auxlarge, this%large%dfft%comm)
            !
            CALL env_scatter_grid(this%large%dfft, auxlarge, flarge)
            !
#else
            flarge = auxlarge
#endif
            DEALLOCATE (auxlarge)
        END IF
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE map_small_to_large_real
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE map_small_to_large_density(this, fsmall, flarge)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_mapping), INTENT(IN) :: this
        TYPE(environ_density), INTENT(IN) :: fsmall
        !
        TYPE(environ_density), INTENT(INOUT) :: flarge
        !
        INTEGER :: i
        REAL(DP), ALLOCATABLE :: auxlarge(:)
        !
        CHARACTER(LEN=80) :: sub_name = 'map_small_to_large_density'
        !
        !--------------------------------------------------------------------------------
        ! Check if input/output dimensions match mapping cells
        !
        IF (.NOT. ASSOCIATED(fsmall%cell, this%small)) &
            CALL io%error(sub_name, "Mismatch of small cell", 1)
        !
        IF (.NOT. ASSOCIATED(flarge%cell, this%large)) &
            CALL io%error(sub_name, "Mismatch of large cell", 1)
        !
        !--------------------------------------------------------------------------------
        ! If the cells are the same, just copy
        !
        IF (ASSOCIATED(this%large, this%small)) THEN
            flarge%of_r = fsmall%of_r
        ELSE
            !
            !----------------------------------------------------------------------------
            ! Copy small cell to corresponding gridpoints in the full large cell
            !
            ALLOCATE (auxlarge(this%large%nnt))
            auxlarge = 0.D0
            !
            DO i = 1, this%small%ir_end
                IF (this%map(i) > 0) auxlarge(this%map(i)) = fsmall%of_r(i)
            END DO
            !
#if defined(__MPI)
            CALL env_mp_sum(auxlarge, this%large%dfft%comm)
            !
            CALL env_scatter_grid(this%large%dfft, auxlarge, flarge%of_r)
            !
#else
            flarge%of_r = auxlarge
#endif
            DEALLOCATE (auxlarge)
        END IF
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE map_small_to_large_density
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE map_large_to_small_real(this, nlarge, nsmall, flarge, fsmall)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_mapping), INTENT(IN) :: this
        INTEGER, INTENT(IN) :: nsmall, nlarge
        REAL(DP), INTENT(IN) :: flarge(nlarge)
        !
        REAL(DP), INTENT(INOUT) :: fsmall(nsmall)
        !
        INTEGER :: i
        REAL(DP), ALLOCATABLE :: auxlarge(:)
        !
        CHARACTER(LEN=80) :: sub_name = 'map_large_to_small_real'
        !
        !--------------------------------------------------------------------------------
        ! Check if input/output dimensions match mapping cells
        !
        IF (nsmall /= this%small%nnr) &
            CALL io%error(sub_name, "Wrong dimension of small cell", 1)
        !
        IF (nlarge /= this%large%nnr) &
            CALL io%error(sub_name, "Wrong dimension of large cell", 1)
        !
        !--------------------------------------------------------------------------------
        ! If the cells are the same, just copy
        !
        IF (nsmall == nlarge) THEN
            fsmall = flarge
        ELSE
            !
            !----------------------------------------------------------------------------
            ! Copy portion of large cell to corresponding gridpoints in the small cell
            !
            ALLOCATE (auxlarge(this%large%nnt))
            auxlarge = 0.D0
#if defined(__MPI)
            CALL env_gather_grid(this%large%dfft, flarge, auxlarge)
            !
            CALL env_mp_sum(auxlarge, this%large%dfft%comm)
            !
#else
            auxlarge = flarge
#endif
            fsmall = 0.D0
            !
            DO i = 1, this%small%ir_end
                IF (this%map(i) > 0) fsmall(i) = auxlarge(this%map(i))
            END DO
            !
            DEALLOCATE (auxlarge)
        END IF
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE map_large_to_small_real
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE map_large_to_small_density(this, flarge, fsmall)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_mapping), INTENT(IN) :: this
        TYPE(environ_density), INTENT(IN) :: flarge
        !
        TYPE(environ_density), INTENT(INOUT) :: fsmall
        !
        INTEGER :: i
        REAL(DP), ALLOCATABLE :: auxlarge(:)
        !
        CHARACTER(LEN=80) :: sub_name = 'map_large_to_small_density'
        !
        !--------------------------------------------------------------------------------
        ! Check if input/output dimensions match mapping cells
        !
        IF (.NOT. ASSOCIATED(fsmall%cell, this%small)) &
            CALL io%error(sub_name, "Mismatch of small cell", 1)
        !
        IF (.NOT. ASSOCIATED(flarge%cell, this%large)) &
            CALL io%error(sub_name, "Mismatch of large cell", 1)
        !
        !--------------------------------------------------------------------------------
        ! If the cells are the same, just copy
        !
        IF (ASSOCIATED(this%large, this%small)) THEN
            fsmall%of_r = flarge%of_r
        ELSE
            !
            !----------------------------------------------------------------------------
            ! Copy portion of large cell to corresponding gridpoints in the small cell
            !
            ALLOCATE (auxlarge(this%large%nnt))
            auxlarge = 0.D0
            !
#if defined(__MPI)
            CALL env_gather_grid(this%large%dfft, flarge%of_r, auxlarge)
            !
            CALL env_mp_sum(auxlarge, this%large%dfft%comm)
            !
#else
            auxlarge = flarge%of_r
#endif
            fsmall%of_r = 0.D0
            !
            DO i = 1, this%small%ir_end
                IF (this%map(i) > 0) fsmall%of_r(i) = auxlarge(this%map(i))
            END DO
            !
            DEALLOCATE (auxlarge)
        END IF
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE map_large_to_small_density
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !                                   OUTPUT METHODS
    !
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !>
    !! Prints the details of the mapping
    !!
    !! @param verbose       : (INTEGER) adds verbosity to global verbose
    !! @param debug_verbose : (INTEGER) replaces global verbose for debugging
    !! @param unit          : (INTEGER) output target (default = io%debug_unit)
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE print_environ_mapping(this, verbose, debug_verbose, unit)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_mapping), INTENT(IN) :: this
        INTEGER, OPTIONAL, INTENT(IN) :: verbose, debug_verbose, unit
        !
        INTEGER :: base_verbose, local_verbose, local_unit
        !
        CHARACTER(LEN=80) :: sub_name = 'print_environ_mapping'
        !
        !--------------------------------------------------------------------------------
        !
        IF (PRESENT(debug_verbose)) THEN
            base_verbose = debug_verbose
            !
            IF (PRESENT(verbose)) THEN
                local_verbose = verbose
            ELSE
                local_verbose = debug_verbose
            END IF
            !
        ELSE IF (io%verbosity > 0) THEN
            base_verbose = io%verbosity
            !
            IF (PRESENT(verbose)) THEN
                local_verbose = base_verbose + verbose
            ELSE
                local_verbose = base_verbose
            END IF
            !
        ELSE
            RETURN
        END IF
        !
        IF (PRESENT(unit)) THEN
            local_unit = unit
        ELSE
            local_unit = io%debug_unit
        END IF
        !
        IF (local_verbose >= 1) THEN
            !
            IF (io%lnode) THEN
                WRITE (local_unit, 1000)
                WRITE (local_unit, 1001) this%nrep
                WRITE (local_unit, 1002) this%small%origin, this%large%origin
            END IF
            !
        END IF
        !
        FLUSH (local_unit)
        !
        !--------------------------------------------------------------------------------
        !
1000    FORMAT(/, 4('%'), " MAPPING ", 67('%'))
!
1001    FORMAT(/, " number of replicas (x,y,z) = ", 3I14)
!
1002    FORMAT(/, " cell origins:", /, &
                " system                     = ", 3F14.7, /, &
                " environment                = ", 3F14.7)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE print_environ_mapping
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE class_mapping
!----------------------------------------------------------------------------------------
