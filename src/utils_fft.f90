!----------------------------------------------------------------------------------------
!>
!!
!----------------------------------------------------------------------------------------
MODULE utils_fft
    !------------------------------------------------------------------------------------
    !
    USE modules_constants, ONLY: DP
    !
    USE core_types, ONLY: fft_core
    USE cell_types, ONLY: environ_cell
    !
    USE tools_generate_gvectors, ONLY: env_gvect_init, env_ggen
    !
    USE correction_mt, ONLY: update_mt_correction
    !
    !------------------------------------------------------------------------------------
CONTAINS
    !------------------------------------------------------------------------------------
    !>
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
    ! BACKWARD COMPATIBILITY
    ! Compatible with QE-6.0 QE-6.1.X QE-6.2.X QE-6.3
    ! SUBROUTINE init_fft_core( fft, use_internal_pbc_corr, nspin )
    ! Compatible with QE-6.4.X QE-GIT
    SUBROUTINE init_fft_core_first(fft, use_internal_pbc_corr)
        ! END BACKWARD COMPATIBILITY
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        LOGICAL, INTENT(IN), OPTIONAL :: use_internal_pbc_corr
        ! BACKWARD COMPATIBILITY
        ! Compatible with QE-6.0 QE-6.1.X QE-6.2.X QE-6.3.X
        ! INTEGER, INTENT(IN), OPTIONAL :: nspin
        ! Compatible with QE-6.4.X QE-GIT
        ! END BACKWARD COMPATIBILITY
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
        ! BACKWARD COMPATIBILITY
        ! Compatible with QE-6.0 QE-6.1.X QE-6.2.X QE-6.3.X
        ! IF (PRESENT(nspin)) THEN
        !     fft%nspin = nspin
        ! ELSE
        !     fft%nspin = 1
        ! END IF
        ! Compatible with QE-6.4.X QE-GIT
        !
        ! END BACKWARD COMPATIBILITY
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
        ! #TODO The following routines are in tools_generate_gvectors
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
