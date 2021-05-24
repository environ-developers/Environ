!----------------------------------------------------------------------------------------
!>
!! Internal setup of numerical cores
!!
!----------------------------------------------------------------------------------------
MODULE base_core
    !------------------------------------------------------------------------------------
    !
    USE types_core, ONLY: fd_core, fft_core, oned_analytic_core
    !
    !------------------------------------------------------------------------------------
    !
    SAVE
    !
    LOGICAL :: lfd
    TYPE(fd_core) :: fd
    !
    LOGICAL :: lfft_system
    TYPE(fft_core) :: system_fft
    !
    LOGICAL :: lfft_environment
    TYPE(fft_core) :: environment_fft
    !
    LOGICAL :: loned_analytic
    TYPE(oned_analytic_core) :: oned_analytic
    !
    !------------------------------------------------------------------------------------
    ! Keeping imports private
    !
    PRIVATE :: fd_core, fft_core, oned_analytic_core
    !
    !------------------------------------------------------------------------------------
END MODULE base_core
!----------------------------------------------------------------------------------------
