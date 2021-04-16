!----------------------------------------------------------------------------------------
!>
!!
!----------------------------------------------------------------------------------------
MODULE utils_mapping
    !------------------------------------------------------------------------------------
    !
    USE cell_types
    USE environ_types
    USE scatter_mod, ONLY: scatter_grid, gather_grid
    USE mp, ONLY: mp_sum
    !
    PRIVATE
    PUBLIC :: map_small_to_large, map_large_to_small
    !
    INTERFACE map_small_to_large
        MODULE PROCEDURE map_small_to_large_real, map_small_to_large_density
    END INTERFACE map_small_to_large
    !
    INTERFACE map_large_to_small
        MODULE PROCEDURE map_large_to_small_real, map_large_to_small_density
    END INTERFACE map_large_to_small
    !
    !------------------------------------------------------------------------------------
CONTAINS
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE map_small_to_large_real(mapping, nsmall, nlarge, fsmall, flarge)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_mapping), INTENT(IN) :: mapping
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
        IF (nsmall /= mapping%small%nnr) &
            CALL errore(sub_name, 'Wrong dimension of small cell', 1)
        !
        IF (nlarge /= mapping%large%nnr) &
            CALL errore(sub_name, 'Wrong dimension of large cell', 1)
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
            ALLOCATE (auxlarge(mapping%large%ntot))
            auxlarge = 0.D0
            !
            DO ir = 1, mapping%small%ir_end
                !
                IF (mapping%map(ir) > 0) &
                    auxlarge(mapping%map(ir)) = fsmall(ir)
                !
            END DO
            !
#if defined (__MPI)
            CALL mp_sum(auxlarge, mapping%large%dfft%comm)
            !
            CALL scatter_grid(mapping%large%dfft, auxlarge, flarge)
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
    SUBROUTINE map_small_to_large_density(mapping, fsmall, flarge)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_mapping), INTENT(IN) :: mapping
        TYPE(environ_density), INTENT(IN) :: fsmall
        !
        TYPE(environ_density), INTENT(INOUT) :: flarge
        !
        INTEGER :: ir
        REAL(DP), ALLOCATABLE :: auxlarge(:)
        !
        CHARACTER(LEN=80) :: sub_name = 'map_small_to_large'
        !
        !--------------------------------------------------------------------------------
        ! Check if input/output dimensions match mapping cells
        !
        IF (.NOT. ASSOCIATED(fsmall%cell, mapping%small)) &
            CALL errore(sub_name, 'Mismatch of small cell', 1)
        !
        IF (.NOT. ASSOCIATED(flarge%cell, mapping%large)) &
            CALL errore(sub_name, 'Mismatch of large cell', 1)
        !
        !--------------------------------------------------------------------------------
        ! If the cells are the same, just copy
        !
        IF (ASSOCIATED(mapping%large, mapping%small)) THEN
            flarge%of_r = fsmall%of_r
        ELSE
            !
            !----------------------------------------------------------------------------
            ! Copy small cell to corresponding gridpoints in the full large cell
            !
            ALLOCATE (auxlarge(mapping%large%ntot))
            auxlarge = 0.D0
            !
            DO ir = 1, mapping%small%ir_end
                !
                IF (mapping%map(ir) > 0) & ! This test may be redundant
                    auxlarge(mapping%map(ir)) = fsmall%of_r(ir)
                !
            END DO
            !
#if defined (__MPI)
            CALL mp_sum(auxlarge, mapping%large%dfft%comm)
            !
            CALL scatter_grid(mapping%large%dfft, auxlarge, flarge%of_r)
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
    SUBROUTINE map_large_to_small_real(mapping, nlarge, nsmall, flarge, fsmall)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_mapping), INTENT(IN) :: mapping
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
        IF (nsmall /= mapping%small%nnr) &
            CALL errore(sub_name, 'Wrong dimension of small cell', 1)
        !
        IF (nlarge /= mapping%large%nnr) &
            CALL errore(sub_name, 'Wrong dimension of large cell', 1)
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
            ALLOCATE (auxlarge(mapping%large%ntot))
            auxlarge = 0.D0
#if defined (__MPI)
            CALL gather_grid(mapping%large%dfft, flarge, auxlarge)
            !
            CALL mp_sum(auxlarge, mapping%large%dfft%comm)
            !
#else
            auxlarge = flarge
#endif
            fsmall = 0.D0
            !
            DO ir = 1, mapping%small%ir_end
                !
                IF (mapping%map(ir) > 0) & ! #TODO This test may be redundant
                    fsmall(ir) = auxlarge(mapping%map(ir))
                !
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
    SUBROUTINE map_large_to_small_density(mapping, flarge, fsmall)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_mapping), INTENT(IN) :: mapping
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
        IF (.NOT. ASSOCIATED(fsmall%cell, mapping%small)) &
            CALL errore(sub_name, 'Mismatch of small cell', 1)
        !
        IF (.NOT. ASSOCIATED(flarge%cell, mapping%large)) &
            CALL errore(sub_name, 'Mismatch of large cell', 1)
        !
        !--------------------------------------------------------------------------------
        ! If the cells are the same, just copy
        !
        IF (ASSOCIATED(mapping%large, mapping%small)) THEN
            fsmall%of_r = flarge%of_r
        ELSE
            !
            !----------------------------------------------------------------------------
            ! Copy portion of large cell to corresponding gridpoints in the small cell
            !
            ALLOCATE (auxlarge(mapping%large%ntot))
            auxlarge = 0.D0
            !
#if defined(__MPI)
            CALL gather_grid(mapping%large%dfft, flarge%of_r, auxlarge)
            !
            CALL mp_sum(auxlarge, mapping%large%dfft%comm)
            !
#else
            auxlarge = flarge%of_r
#endif
            fsmall%of_r = 0.D0
            !
            DO ir = 1, mapping%small%ir_end
                !
                IF (mapping%map(ir) > 0) & ! This test may be redundant
                    fsmall%of_r(ir) = auxlarge(mapping%map(ir))
                !
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
END MODULE utils_mapping
!----------------------------------------------------------------------------------------
