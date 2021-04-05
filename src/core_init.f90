!----------------------------------------------------------------------------------------
!>
!!
!----------------------------------------------------------------------------------------
MODULE core_init
    !------------------------------------------------------------------------------------
    !
    USE core_base
    !
    !------------------------------------------------------------------------------------
CONTAINS
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE set_core_base(ifdtype, nfdpoint, use_internal_pbc_corr, dim, axis)
        !--------------------------------------------------------------------------------
        !
        USE utils_fd, ONLY: init_fd_core_first
        USE utils_oned_analytic, ONLY: init_oned_analytic_core_first
        USE utils_fft, ONLY: init_fft_core_first
        !
        LOGICAL, INTENT(IN) :: use_internal_pbc_corr
        INTEGER, INTENT(IN) :: ifdtype, nfdpoint, dim, axis
        !
        !--------------------------------------------------------------------------------
        ! Set up active numerical cores
        !
        IF (lfd) CALL init_fd_core_first(ifdtype, nfdpoint, fd)
        !
        IF (loned_analytic) CALL init_oned_analytic_core_first(dim, axis, oned_analytic)
        !
        IF (lfft_system) CALL init_fft_core_first(system_fft, use_internal_pbc_corr)
        !
        IF (lfft_environment) &
            CALL init_fft_core_first(environment_fft, use_internal_pbc_corr)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE set_core_base
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE core_initbase(gcutm, environment_cell, system_cell)
        !--------------------------------------------------------------------------------
        !
        USE utils_oned_analytic, ONLY: init_oned_analytic_core_second
        USE utils_fd, ONLY: init_fd_core_second
        USE utils_fft, ONLY: init_fft_core_second
        !
        IMPLICIT NONE
        !
        REAL(DP), INTENT(IN) :: gcutm
        TYPE(environ_cell), INTENT(IN) :: system_cell, environment_cell
        !
        !--------------------------------------------------------------------------------
        ! calc is used so oned_analytic isn't calculated twice #TODO What?
        !
        IF (loned_analytic) &
            CALL init_oned_analytic_core_second(environment_cell, oned_analytic)
        !
        IF (lfd) CALL init_fd_core_second(environment_cell, fd)
        !
        IF (lfft_environment) &
            CALL init_fft_core_second(gcutm, environment_cell, environment_fft)
        !
        IF (lfft_system) CALL init_fft_core_second(gcutm, system_cell, system_fft)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE core_initbase
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE core_initcell(system_cell, environment_cell)
        !--------------------------------------------------------------------------------
        !
        USE utils_oned_analytic, ONLY: update_oned_analytic_core_cell
        USE utils_fft, ONLY: update_fft_core_cell
        !
        IMPLICIT NONE
        !
        TYPE(environ_cell), INTENT(IN) :: system_cell, environment_cell
        !
        !--------------------------------------------------------------------------------
        !
        IF (loned_analytic) &
            CALL update_oned_analytic_core_cell(environment_cell, oned_analytic)
        !
        IF (lfft_environment) &
            CALL update_fft_core_cell(environment_cell, environment_fft)
        !
        IF (lfft_system) CALL update_fft_core_cell(system_cell, system_fft)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE core_initcell
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE core_initions(pos)
        !--------------------------------------------------------------------------------
        !
        USE utils_oned_analytic, ONLY: update_oned_analytic_core_origin
        !
        IMPLICIT NONE
        !
        REAL(DP), INTENT(IN) :: pos(3)
        !
        !--------------------------------------------------------------------------------
        !
        IF (loned_analytic) CALL update_oned_analytic_core_origin(pos, oned_analytic)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE core_initions
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE core_clean(lflag)
        !--------------------------------------------------------------------------------
        !
        USE utils_oned_analytic, ONLY: destroy_oned_analytic_core
        USE utils_fd, ONLY: destroy_fd_core
        USE utils_fft, ONLY: destroy_fft_core
        !
        IMPLICIT NONE
        !
        LOGICAL, INTENT(IN) :: lflag
        !
        !--------------------------------------------------------------------------------
        !
        IF (lfft_environment) CALL destroy_fft_core(lflag, environment_fft)
        !
        IF (lfft_system) CALL destroy_fft_core(lflag, system_fft)
        !
        IF (lfd) CALL destroy_fd_core(lflag, fd)
        !
        IF (loned_analytic) CALL destroy_oned_analytic_core(lflag, oned_analytic)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE core_clean
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE core_init
!----------------------------------------------------------------------------------------
