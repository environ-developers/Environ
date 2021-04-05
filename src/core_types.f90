!----------------------------------------------------------------------------------------
!>
!!
!----------------------------------------------------------------------------------------
MODULE core_types
    !------------------------------------------------------------------------------------
    !
    USE modules_constants, ONLY: DP
    USE cell_types, ONLY: environ_cell
    !
    !------------------------------------------------------------------------------------
    !>
    !! # TODO: add descriptive comments to variables
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
        ! BACKWARD COMPATIBILITY
        ! Compatible with QE-6.0 QE-6.1.X QE-6.2.X QE-6.3.X
        ! INTEGER :: nspin
        ! Compatible with QE-6.4.X QE-GIT
        !
        ! END BACKWARD COMPATIBILITY
        LOGICAL :: use_internal_pbc_corr = .FALSE.
        !
        TYPE(environ_cell), POINTER :: cell
        !
        INTEGER :: ngm = 0  ! local  number of G vectors (on this processor)
        ! with gamma tricks, only vectors in G>
        !
        REAL(DP) :: gcutm = 0.0_DP   ! ecutrho/(2 pi/a)^2, cut-off for |G|^2
        !
        INTEGER :: gstart = 2 ! index of the first G vector whose module is > 0
        ! Needed in parallel execution:
        ! gstart=2 for the proc that holds G=0
        ! gstart=1 for all others
        !
        REAL(DP), ALLOCATABLE :: gg(:)
        ! G^2 in increasing order (in units of tpiba2=(2pi/a)^2)
        !
        REAL(DP), ALLOCATABLE :: g(:, :)
        ! G-vectors cartesian components ( in units tpiba =(2pi/a)  )
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
        INTEGER :: n, d, p, axis ! #TODO: add comments
        REAL(DP) :: size, origin(3)
        REAL(DP), ALLOCATABLE :: x(:, :) ! #TODO: add comment
        !
        TYPE(environ_cell), POINTER :: cell
        !
        !--------------------------------------------------------------------------------
    END TYPE oned_analytic_core
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE core_types
!----------------------------------------------------------------------------------------
