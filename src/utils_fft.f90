!----------------------------------------------------------------------------------------
!>
!!
!----------------------------------------------------------------------------------------
MODULE utils_fft
    !------------------------------------------------------------------------------------
    !
    USE environ_param, ONLY: DP
    !
    USE types_core, ONLY: fft_core
    USE types_cell, ONLY: environ_cell
    !
    USE generate_gvectors, ONLY: env_gvect_init, env_ggen
    !
    USE correction_mt, ONLY: update_mt_correction
    !
    !------------------------------------------------------------------------------------
    !
    PRIVATE
    !
    PUBLIC :: init_fft_core_first, init_fft_core_second, update_fft_core_cell, &
              destroy_fft_core
    !
    !------------------------------------------------------------------------------------
CONTAINS
    !------------------------------------------------------------------------------------
    !>
    !! #TODO unused
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE create_fft_core(fft)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(fft_core), INTENT(INOUT) :: fft
        !
        !--------------------------------------------------------------------------------
        !
        NULLIFY (fft%cell) ! create empty fft core
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE create_fft_core
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_fft_core_first(fft, use_internal_pbc_corr)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        LOGICAL, INTENT(IN), OPTIONAL :: use_internal_pbc_corr
        !
        TYPE(fft_core), INTENT(INOUT) :: fft
        !
        !--------------------------------------------------------------------------------
        !
        fft%index = 1
        !
        IF (PRESENT(use_internal_pbc_corr)) THEN
            fft%use_internal_pbc_corr = use_internal_pbc_corr
        ELSE
            fft%use_internal_pbc_corr = .FALSE.
        END IF
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_fft_core_first
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_fft_core_second(gcutm, cell, fft)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP), INTENT(IN) :: gcutm
        TYPE(environ_cell), TARGET, INTENT(IN) :: cell
        !
        TYPE(fft_core), INTENT(INOUT) :: fft
        !
        INTEGER :: ngm_g ! global number of G vectors (summed on all procs)
        ! in serial execution, ngm_g = ngm
        !
        !--------------------------------------------------------------------------------
        !
        fft%gcutm = gcutm
        fft%cell => cell
        !
        fft%ngm = cell%dfft%ngm
        !
        IF (fft%use_internal_pbc_corr) ALLOCATE (fft%mt_corr(fft%ngm))
        !
        !--------------------------------------------------------------------------------
        ! #TODO The following routines are in generate_gvectors
        ! and may need to be simplified
        !
        CALL env_gvect_init(fft, ngm_g, cell%dfft%comm)
        !
        CALL env_ggen(fft%cell%dfft, cell%dfft%comm, cell%at, cell%bg, fft%gcutm, &
                      ngm_g, fft%ngm, fft%g, fft%gg, fft%gstart, .TRUE.)
        !
        IF (fft%use_internal_pbc_corr) CALL update_mt_correction(fft)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_fft_core_second
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE update_fft_core_cell(cell, fft)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_cell), TARGET, INTENT(IN) :: cell
        !
        TYPE(fft_core), INTENT(INOUT) :: fft
        !
        !--------------------------------------------------------------------------------
        !
        fft%cell => cell
        !
        IF (fft%use_internal_pbc_corr) CALL update_mt_correction(fft)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE update_fft_core_cell
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE destroy_fft_core(lflag, fft)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        LOGICAL, INTENT(IN) :: lflag
        !
        TYPE(fft_core), INTENT(INOUT) :: fft
        !
        CHARACTER(LEN=80) :: sub_name = 'destroy_fft_core'
        !
        !--------------------------------------------------------------------------------
        !
        NULLIFY (fft%cell)
        DEALLOCATE (fft%gg)
        DEALLOCATE (fft%g)
        !
        IF (fft%use_internal_pbc_corr) DEALLOCATE (fft%mt_corr)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE destroy_fft_core
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE utils_fft
!----------------------------------------------------------------------------------------
