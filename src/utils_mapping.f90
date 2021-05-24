!----------------------------------------------------------------------------------------
!>
!!
!----------------------------------------------------------------------------------------
MODULE utils_mapping
    !------------------------------------------------------------------------------------
    !
    USE environ_param, ONLY: DP
    !
    USE types_cell
    !
    USE tools_cell, ONLY: ir2ijk
    !
    !------------------------------------------------------------------------------------
    !
    PRIVATE
    !
    PUBLIC :: create_environ_mapping, init_environ_mapping_first, &
              init_environ_mapping_second, update_environ_mapping, &
              destroy_environ_mapping
    !
    !------------------------------------------------------------------------------------
CONTAINS
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE create_environ_mapping(mapping)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_mapping), INTENT(INOUT) :: mapping
        !
        !--------------------------------------------------------------------------------
        !
        mapping%nrep = 1
        mapping%large => NULL()
        mapping%small => NULL()
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE create_environ_mapping
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_environ_mapping_first(nrep, mapping)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: nrep(3)
        !
        TYPE(environ_mapping), INTENT(INOUT) :: mapping
        !
        !--------------------------------------------------------------------------------
        !
        mapping%nrep = nrep
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_environ_mapping_first
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_environ_mapping_second(small, large, mapping)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_cell), TARGET, INTENT(IN) :: small, large
        !
        TYPE(environ_mapping), INTENT(INOUT) :: mapping
        !
        !--------------------------------------------------------------------------------
        ! Check that small%at and large%at are compatible with mapping%nrep
        !
        mapping%small => small
        mapping%large => large
        !
        IF (.NOT. ASSOCIATED(mapping%small, mapping%large)) &
            ALLOCATE (mapping%map(mapping%small%nnr))
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_environ_mapping_second
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE update_environ_mapping(mapping, pos)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP), INTENT(IN), OPTIONAL :: pos(3)
        !
        TYPE(environ_mapping), INTENT(INOUT) :: mapping
        !
        LOGICAL :: physical
        INTEGER :: ir, ipol
        INTEGER, DIMENSION(3) :: small_n, large_n, center, origin, shift, ijk
        REAL(DP) :: tmp(3)
        !
        !--------------------------------------------------------------------------------
        !
        IF (ASSOCIATED(mapping%small, mapping%large)) RETURN
        ! environment cell is identical to system cell (env_nrep = 0 in environ.in)
        !
        !--------------------------------------------------------------------------------
        ! Compute mapping
        !
        small_n(1) = mapping%small%dfft%nr1
        small_n(2) = mapping%small%dfft%nr2
        small_n(3) = mapping%small%dfft%nr3
        !
        large_n(1) = mapping%large%dfft%nr1
        large_n(2) = mapping%large%dfft%nr2
        large_n(3) = mapping%large%dfft%nr3
        !
        !--------------------------------------------------------------------------------
        !
        center = NINT(small_n / 2.D0) ! Indexes of center of small cell
        !
        !--------------------------------------------------------------------------------
        ! Indexes of origin of small cell
        !
        IF (PRESENT(pos)) THEN
            tmp = MATMUL(mapping%small%bg, pos) ! #TODO center of charge (think molecule)
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
        mapping%large%origin = -MATMUL(mapping%large%at, (0.5 - origin / DBLE(large_n)))
        !
        mapping%map = 0
        !
        DO ir = 1, mapping%small%ir_end
            !
            CALL ir2ijk(mapping%small, ir, ijk(1), ijk(2), ijk(3), physical)
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
            ijk = ijk + small_n * mapping%nrep
            !
            mapping%map(ir) = 1 + ijk(1) & ! x-point
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
    SUBROUTINE destroy_environ_mapping(lflag, mapping)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        LOGICAL, INTENT(IN) :: lflag
        !
        TYPE(environ_mapping), INTENT(INOUT) :: mapping
        !
        !--------------------------------------------------------------------------------
        !
        mapping%nrep = 1
        NULLIFY (mapping%small)
        NULLIFY (mapping%large)
        !
        IF (ALLOCATED(mapping%map)) DEALLOCATE (mapping%map)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE destroy_environ_mapping
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE utils_mapping
!----------------------------------------------------------------------------------------
