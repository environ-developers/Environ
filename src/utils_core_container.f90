!----------------------------------------------------------------------------------------
!>
!!
!----------------------------------------------------------------------------------------
MODULE utils_core_container
    !------------------------------------------------------------------------------------
    !
    USE environ_param, ONLY: DP, tpi
    !
    USE types_core
    !
    !------------------------------------------------------------------------------------
    !
    PRIVATE
    !
    PUBLIC :: create_core_container, init_core_container, destroy_core_container
    !
    !------------------------------------------------------------------------------------
CONTAINS
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE create_core_container(core)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(core_container), INTENT(INOUT) :: core
        !
        !--------------------------------------------------------------------------------
        !
        core%type_ = 'default'
        !
        core%use_fd = .FALSE.
        NULLIFY (core%fd)
        !
        core%use_fft = .FALSE.
        NULLIFY (core%fft)
        !
        core%use_oned_analytic = .FALSE.
        NULLIFY (core%oned_analytic)
        !
        core%need_correction = .FALSE.
        NULLIFY (core%correction)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE create_core_container
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_core_container(type_, core, fd, fft, oned_analytic)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CHARACTER(LEN=80), INTENT(IN) :: type_
        TYPE(fd_core), INTENT(IN), TARGET, OPTIONAL :: fd
        TYPE(fft_core), INTENT(IN), TARGET, OPTIONAL :: fft
        TYPE(oned_analytic_core), INTENT(IN), TARGET, OPTIONAL :: oned_analytic
        !
        TYPE(core_container), INTENT(INOUT) :: core
        !
        INTEGER :: number
        !
        CHARACTER(LEN=80) :: sub_name = 'init_core_container'
        !
        !--------------------------------------------------------------------------------
        ! Assign the selected numerical core
        !
        core%type_ = type_
        !
        SELECT CASE (TRIM(ADJUSTL(type_)))
            !
        CASE ('fd')
            !
            !----------------------------------------------------------------------------
            ! Note: finite differences core only works for gradient,
            ! other derivatives still require fft core
            !
            IF (.NOT. PRESENT(fd)) &
                CALL env_errore(sub_name, 'Missing specified core type', 1)
            !
            IF (.NOT. PRESENT(fft)) &
                CALL env_errore(sub_name, 'Missing specified core type', 1)
            !
            core%use_fd = .TRUE.
            core%use_fft = .TRUE.
            core%fft => fft
            core%fd => fd
            !
        CASE ('chain', 'fft', 'highmem', 'lowmem')
            !
            IF (.NOT. PRESENT(fft)) &
                CALL env_errore(sub_name, 'Missing specified core type', 1)
            !
            core%use_fft = .TRUE.
            core%fft => fft
            !
        CASE ('1da', '1d-analytic', 'oned_analytic', 'gcs', 'gouy-chapman', &
              'gouy-chapman-stern', 'ms', 'mott-schottky', 'ms-gcs', &
              'mott-schottky-gouy-chapman-stern')
            !
            IF (.NOT. PRESENT(oned_analytic)) &
                CALL env_errore(sub_name, 'Missing specified core type', 1)
            !
            core%use_oned_analytic = .TRUE.
            core%oned_analytic => oned_analytic
            !
        CASE DEFAULT
            CALL env_errore(sub_name, 'Unexpected keyword for core_container type', 1)
            !
        END SELECT
        !
        !--------------------------------------------------------------------------------
        ! Double check number of active cores
        !
        number = 0
        !
        IF (core%use_fd) number = number + 1
        !
        IF (core%use_fft) number = number + 1
        !
        IF (core%use_oned_analytic) number = number + 1
        !
        IF (number /= 1) &
            CALL env_errore(sub_name, 'Incorrect number of active cores', 1)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_core_container
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE destroy_core_container(lflag, core)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        LOGICAL, INTENT(IN) :: lflag
        !
        TYPE(core_container), INTENT(INOUT) :: core
        !
        !--------------------------------------------------------------------------------
        !
        IF (lflag) THEN
            core%use_fd = .FALSE.
            NULLIFY (core%fd)
            core%use_fft = .FALSE.
            NULLIFY (core%fft)
            core%use_oned_analytic = .FALSE.
            NULLIFY (core%oned_analytic)
            core%need_correction = .FALSE.
            NULLIFY (core%correction)
        END IF
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE destroy_core_container
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE utils_core_container
!----------------------------------------------------------------------------------------
