!----------------------------------------------------------------------------------------
!>
!! Internal setup of numerical cores
!!
!----------------------------------------------------------------------------------------
MODULE core_base
    !------------------------------------------------------------------------------------
    !
    USE core_types
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
END MODULE core_base
!----------------------------------------------------------------------------------------
