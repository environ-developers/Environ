!----------------------------------------------------------------------------------------
!>
!!
!----------------------------------------------------------------------------------------
MODULE types_core
    !------------------------------------------------------------------------------------
    !
    USE environ_param, ONLY: DP
    !
    USE types_cell, ONLY: environ_cell
    !
    !------------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    TYPE fd_core
        !--------------------------------------------------------------------------------
        !
        INTEGER :: ifdtype
        INTEGER :: nfdpoint
        INTEGER, ALLOCATABLE :: icfd(:)
        INTEGER :: ncfd
        !
        TYPE(environ_cell), POINTER :: cell
        !
        !--------------------------------------------------------------------------------
    END TYPE fd_core
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    TYPE fft_core
        !--------------------------------------------------------------------------------
        !
        INTEGER :: index
        LOGICAL :: use_internal_pbc_corr = .FALSE.
        !
        TYPE(environ_cell), POINTER :: cell
        !
        INTEGER :: ngm = 0 ! local  number of G vectors (on this processor)
        ! with gamma tricks, only vectors in G>
        !
        REAL(DP) :: gcutm = 0.0_DP ! ecutrho/(2 pi/a)^2, cut-off for |G|^2
        !
        INTEGER :: gstart = 2 ! index of the first G vector whose module is > 0
        ! needed in parallel execution:
        ! gstart=2 for the proc that holds G=0
        ! gstart=1 for all others
        !
        REAL(DP), ALLOCATABLE :: gg(:)
        ! G^2 in increasing order (in units of tpiba2=(2pi/a)^2)
        !
        REAL(DP), ALLOCATABLE :: g(:, :)
        ! G-vectors cartesian components ( in units tpiba =(2pi/a) )
        !
        !--------------------------------------------------------------------------------
        ! Martyna-Tuckerman correction
        !
        REAL(DP) :: alpha, beta
        REAL(DP), ALLOCATABLE :: mt_corr(:)
        !
        !--------------------------------------------------------------------------------
    END TYPE fft_core
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    TYPE oned_analytic_core
        !--------------------------------------------------------------------------------
        !
        LOGICAL :: initialized = .FALSE.
        INTEGER :: n, d, p, axis
        REAL(DP) :: size, origin(3)
        REAL(DP), ALLOCATABLE :: x(:, :)
        !
        TYPE(environ_cell), POINTER :: cell
        !
        !--------------------------------------------------------------------------------
    END TYPE oned_analytic_core
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    TYPE core_container
        !--------------------------------------------------------------------------------
        !
        CHARACTER(LEN=80) :: type_
        !
        LOGICAL :: use_fd
        TYPE(fd_core), POINTER :: fd => NULL()
        !
        LOGICAL :: use_fft
        TYPE(fft_core), POINTER :: fft => NULL()
        !
        LOGICAL :: use_oned_analytic
        TYPE(oned_analytic_core), POINTER :: oned_analytic => NULL()
        !
        ! #TODO future work
        !
        ! LOGICAL :: use_oned_numeric
        ! TYPE(oned_numeric_core), POINTER :: oned_numeric => NULL()
        !
        ! LOGICAL :: use_multigrid
        ! TYPE(multigrid_core), POINTER :: multigrid => NULL()
        !
        ! LOGICAL :: use_bigdft
        ! TYPE(bigdft_core), POINTER :: bigdft => NULL()
        !
        LOGICAL :: need_correction
        TYPE(core_container), POINTER :: correction => NULL()
        !
        !--------------------------------------------------------------------------------
    END TYPE core_container
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
    ! Keeping imports private
    !
    PRIVATE :: DP, environ_cell
    !
    !------------------------------------------------------------------------------------
END MODULE types_core
!----------------------------------------------------------------------------------------