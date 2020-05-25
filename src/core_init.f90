MODULE core_init
  !
  USE core_types
  USE core_base
  !
CONTAINS
!--------------------------------------------------------------------
  SUBROUTINE set_core_base( ifdtype, nfdpoint, use_internal_pbc_corr, dim, axis )
!--------------------------------------------------------------------
    !
    LOGICAL, INTENT(IN) :: use_internal_pbc_corr
    INTEGER, INTENT(IN) :: ifdtype, nfdpoint, dim, axis
    !
    ! Set up active numerical cores
    !
    IF ( lfd ) CALL init_fd_core_first( ifdtype, nfdpoint, fd )
! BACKWARD COMPATIBILITY
! Compatible with QE-6.0 QE-6.1.X QE-6.2.X QE-6.3.X
!    IF ( fft ) CALL init_fft_core( fft, use_internal_pbc_corr, nspin )
! Compatible with QE-6.4.X QE-GIT
    IF ( lfft ) CALL init_fft_core_first( fft, use_internal_pbc_corr )
! END BACKWARD COMPATIBILITY
    IF ( loned_analytic ) CALL init_oned_analytic_core_first( dim, axis, oned_analytic )
    !
    RETURN
!--------------------------------------------------------------------
  END SUBROUTINE set_core_base
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE core_initbase( cell, dfft, tpiba, tpiba2, ngm, gcutm, gstart, g, gg  )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    TYPE( environ_cell ), INTENT(IN) :: cell
    TYPE( fft_type_descriptor ), INTENT(IN) :: dfft
    INTEGER, INTENT(IN) :: ngm, gstart
    REAL( DP ), INTENT(IN) :: gcutm, tpiba, tpiba2
    REAL( DP ), DIMENSION(3,ngm) :: g
    REAL( DP ), DIMENSION(ngm) :: gg
    !
    IF ( loned_analytic ) CALL init_oned_analytic_core_second( cell, oned_analytic )
    !
    IF ( lfd ) CALL init_fd_core_secdon( cell, dfft, fd )
    !
    IF ( lfft ) CALL init_fft_core_second( dfft, cell%omega, tpiba, tpiba2, ngm, gcutm, gstart, g, gg, fft )
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE core_initbase
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE core_initcell( cell )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    TYPE( environ_cell ), INTENT(IN) :: cell
    !
    IF ( loned_analytic ) CALL update_oned_analytic_core_cell( cell, oned_analytic )
    !
    IF ( lfft ) CALL update_fft_core_cell( cell%omega, fft ) ! THIS SHOULD NOT BE USED AND NEEDS TO BE FIXED
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE core_initcell
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE core_initions( pos )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    REAL( DP ), DIMENSION(3), INTENT(IN) :: pos
    !
    IF ( loned_analytic ) CALL update_oned_analytic_core_origin( pos, oned_analytic )
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE core_initions
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE core_clean( lflag )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    LOGICAL, INTENT(IN) :: lflag
    !
    IF ( lfft ) CALL destroy_fft_core( lflag, fft )
    !
    IF ( lfd ) CALL destroy_fd_core( lflag, fd )
    !
    IF ( loned_analytic ) CALL destroy_oned_analytic_core( lflag, oned_analytic )
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE core_clean
!--------------------------------------------------------------------
END MODULE core_init
