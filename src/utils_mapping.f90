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
    USE env_mp, ONLY: env_mp_sum
    !
    USE env_base_scatter, ONLY: env_scatter_grid, env_gather_grid
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
        INTEGER :: nrep(3)
        ! number of system cells in environment cell along a_i = 2 * nrep_i + 1
        !
        TYPE(environ_cell), POINTER :: small ! system cell
        TYPE(environ_cell), POINTER :: large ! environment cell
        !
        INTEGER, ALLOCATABLE :: map(:)
        !
        !--------------------------------------------------------------------------------
    CONTAINS
        !--------------------------------------------------------------------------------
        !
        PROCEDURE :: create => create_environ_mapping
        PROCEDURE :: init_first => init_environ_mapping_first
        PROCEDURE :: init_second => init_environ_mapping_second
        PROCEDURE :: update => update_environ_mapping
        PROCEDURE :: destroy => destroy_environ_mapping
        !
        PROCEDURE, PRIVATE :: &
            map_small_to_large_real, map_small_to_large_density, &
            map_large_to_small_real, map_large_to_small_density
        !
        GENERIC :: to_large => map_small_to_large_real, map_small_to_large_density
        GENERIC :: to_small => map_large_to_small_real, map_large_to_small_density
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
        !--------------------------------------------------------------------------------
        !
        this%nrep = 1
        this%large => NULL()
        this%small => NULL()
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE create_environ_mapping
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_environ_mapping_first(this, nrep)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: nrep(3)
        !
        CLASS(environ_mapping), INTENT(INOUT) :: this
        !
        !--------------------------------------------------------------------------------
        !
        this%nrep = nrep
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_environ_mapping_first
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_environ_mapping_second(this, small, large)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_cell), TARGET, INTENT(IN) :: small, large
        !
        CLASS(environ_mapping), INTENT(INOUT) :: this
        !
        !--------------------------------------------------------------------------------
        ! Check that small%at and large%at are compatible with this%nrep
        !
        this%small => small
        this%large => large
        !
        IF (.NOT. ASSOCIATED(this%small, this%large)) &
            ALLOCATE (this%map(this%small%nnr))
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_environ_mapping_second
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE update_environ_mapping(this, pos)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP), INTENT(IN), OPTIONAL :: pos(3)
        !
        CLASS(environ_mapping), INTENT(INOUT) :: this
        !
        LOGICAL :: physical
        INTEGER :: ir, ipol
        INTEGER, DIMENSION(3) :: small_n, large_n, center, origin, shift, ijk
        REAL(DP) :: tmp(3)
        !
        !--------------------------------------------------------------------------------
        !
        IF (ASSOCIATED(this%small, this%large)) RETURN
        ! environment cell is identical to system cell (env_nrep = 0 in environ.in)
        !
        !--------------------------------------------------------------------------------
        ! Compute mapping
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
        !
        center = NINT(small_n / 2.D0) ! Indexes of center of small cell
        !
        !--------------------------------------------------------------------------------
        ! Indexes of origin of small cell
        !
        IF (PRESENT(pos)) THEN
            tmp = MATMUL(this%small%bg, pos) ! #TODO center of charge (think molecule)
            origin = NINT(tmp * small_n)
        ELSE
            origin = 0 ! center of charge
        END IF
        !
        shift = center - origin
        !
        !--------------------------------------------------------------------------------
        ! Shift origin of large cell
        !
        this%large%origin = -MATMUL(this%large%at, (0.5 - origin / DBLE(large_n)))
        !
        this%map = 0
        !
        DO ir = 1, this%small%ir_end
            !
            CALL this%small%ir2ijk(ir, ijk(1), ijk(2), ijk(3), physical)
            !
            IF (.NOT. physical) CYCLE
            !
            !----------------------------------------------------------------------------
            ! Shift to center small cell
            !
            ijk = ijk + shift
            ijk = ijk - FLOOR(DBLE(ijk) / small_n) * small_n ! enforce periodicity
            !
            !----------------------------------------------------------------------------
            ! Map small cell to large cell #TODO check if this works in parallel
            !
            ijk = ijk + small_n * this%nrep
            !
            this%map(ir) = 1 + ijk(1) & ! x-point
                           & + ijk(2) * large_n(1) & ! y-row
                           & + ijk(3) * large_n(1) * large_n(2) ! z-plane
            !
        END DO
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE update_environ_mapping
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE destroy_environ_mapping(this, lflag)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        LOGICAL, INTENT(IN) :: lflag
        !
        CLASS(environ_mapping), INTENT(INOUT) :: this
        !
        !--------------------------------------------------------------------------------
        !
        this%nrep = 1
        NULLIFY (this%small)
        NULLIFY (this%large)
        !
        IF (ALLOCATED(this%map)) DEALLOCATE (this%map)
        !
        RETURN
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
        INTEGER :: ir
        REAL(DP), ALLOCATABLE :: auxlarge(:)
        !
        CHARACTER(LEN=80) :: sub_name = 'map_small_to_large_real'
        !
        !--------------------------------------------------------------------------------
        ! Check if input/output dimensions match mapping cells
        !
        IF (nsmall /= this%small%nnr) &
            CALL env_errore(sub_name, 'Wrong dimension of small cell', 1)
        !
        IF (nlarge /= this%large%nnr) &
            CALL env_errore(sub_name, 'Wrong dimension of large cell', 1)
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
            ALLOCATE (auxlarge(this%large%ntot))
            auxlarge = 0.D0
            !
            DO ir = 1, this%small%ir_end
                IF (this%map(ir) > 0) auxlarge(this%map(ir)) = fsmall(ir)
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
        RETURN
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
        INTEGER :: ir
        REAL(DP), ALLOCATABLE :: auxlarge(:)
        !
        CHARACTER(LEN=80) :: sub_name = 'map_small_to_large_density'
        !
        !--------------------------------------------------------------------------------
        ! Check if input/output dimensions match mapping cells
        !
        IF (.NOT. ASSOCIATED(fsmall%cell, this%small)) &
            CALL env_errore(sub_name, 'Mismatch of small cell', 1)
        !
        IF (.NOT. ASSOCIATED(flarge%cell, this%large)) &
            CALL env_errore(sub_name, 'Mismatch of large cell', 1)
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
            ALLOCATE (auxlarge(this%large%ntot))
            auxlarge = 0.D0
            !
            DO ir = 1, this%small%ir_end
                IF (this%map(ir) > 0) auxlarge(this%map(ir)) = fsmall%of_r(ir)
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
        RETURN
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
        INTEGER :: ir
        REAL(DP), ALLOCATABLE :: auxlarge(:)
        !
        CHARACTER(LEN=80) :: sub_name = 'map_large_to_small_real'
        !
        !--------------------------------------------------------------------------------
        ! Check if input/output dimensions match mapping cells
        !
        IF (nsmall /= this%small%nnr) &
            CALL env_errore(sub_name, 'Wrong dimension of small cell', 1)
        !
        IF (nlarge /= this%large%nnr) &
            CALL env_errore(sub_name, 'Wrong dimension of large cell', 1)
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
            ALLOCATE (auxlarge(this%large%ntot))
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
            DO ir = 1, this%small%ir_end
                IF (this%map(ir) > 0) fsmall(ir) = auxlarge(this%map(ir))
            END DO
            !
            DEALLOCATE (auxlarge)
        END IF
        !
        RETURN
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
        INTEGER :: ir
        REAL(DP), ALLOCATABLE :: auxlarge(:)
        !
        CHARACTER(LEN=80) :: sub_name = 'map_large_to_small_density'
        !
        !--------------------------------------------------------------------------------
        ! Check if input/output dimensions match mapping cells
        !
        IF (.NOT. ASSOCIATED(fsmall%cell, this%small)) &
            CALL env_errore(sub_name, 'Mismatch of small cell', 1)
        !
        IF (.NOT. ASSOCIATED(flarge%cell, this%large)) &
            CALL env_errore(sub_name, 'Mismatch of large cell', 1)
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
            ALLOCATE (auxlarge(this%large%ntot))
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
            DO ir = 1, this%small%ir_end
                IF (this%map(ir) > 0) fsmall%of_r(ir) = auxlarge(this%map(ir))
            END DO
            !
            DEALLOCATE (auxlarge)
        END IF
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE map_large_to_small_density
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE class_mapping
!----------------------------------------------------------------------------------------
