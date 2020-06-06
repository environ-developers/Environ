MODULE utils_fft
  !
  USE core_types
  !
  PRIVATE
  !
  PUBLIC :: create_fft_core, init_fft_core_first, init_fft_core_second, &
       update_fft_core_cell, destroy_fft_core
  !
CONTAINS
  !--------------------------------------------------------------------
  SUBROUTINE create_fft_core( fft )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    TYPE( fft_core ), INTENT(INOUT) :: fft
    !
    ! Create empty fft core
    !
    NULLIFY(fft%dfft)
    NULLIFY(fft%cell)
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE create_fft_core
!--------------------------------------------------------------------
!--------------------------------------------------------------------
! BACKWARD COMPATIBILITY
! Compatible with QE-6.0 QE-6.1.X QE-6.2.X QE-6.3
!  SUBROUTINE init_fft_core( fft, use_internal_pbc_corr, nspin )
! Compatible with QE-6.4.X QE-GIT
  SUBROUTINE init_fft_core_first( fft, use_internal_pbc_corr )
! END BACKWARD COMPATIBILITY
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    TYPE( fft_core ), INTENT(INOUT) :: fft
    LOGICAL, INTENT(IN), OPTIONAL :: use_internal_pbc_corr
! BACKWARD COMPATIBILITY
! Compatible with QE-6.0 QE-6.1.X QE-6.2.X QE-6.3.X
!    INTEGER, INTENT(IN), OPTIONAL :: nspin
! Compatible with QE-6.4.X QE-GIT
!
! END BACKWARD COMPATIBILITY
    !
    fft % index = 1
    !
    IF ( PRESENT ( use_internal_pbc_corr ) ) THEN
       fft % use_internal_pbc_corr = use_internal_pbc_corr
    ELSE
       fft % use_internal_pbc_corr = .FALSE.
    END IF
    !
! BACKWARD COMPATIBILITY
! Compatible with QE-6.0 QE-6.1.X QE-6.2.X QE-6.3.X
!    IF ( PRESENT( nspin ) ) THEN
!       fft % nspin = nspin
!    ELSE
!       fft % nspin = 1
!    ENDIF
! Compatible with QE-6.4.X QE-GIT
!
! END BACKWARD COMPATIBILITY
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE init_fft_core_first
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE init_fft_core_second( dfft, cell, ngm, &
       & gcutm, gstart, g, gg, fft )
!--------------------------------------------------------------------
    !
    USE correction_mt, ONLY : update_mt_correction
    !
    IMPLICIT NONE
    !
    TYPE( fft_type_descriptor ), TARGET, INTENT(IN) :: dfft
    TYPE( environ_cell ), TARGET, INTENT(IN) :: cell
    INTEGER, INTENT(IN) :: ngm, gstart
    REAL( DP ), INTENT(IN) :: gcutm
    REAL( DP ), DIMENSION( ngm ), INTENT(IN) :: gg
    REAL( DP ), DIMENSION( 3, ngm ), INTENT(IN) :: g
    TYPE( fft_core ), INTENT(INOUT) :: fft
    !
    ! COPYING ALL THE COMPONENTS OF DFFT WOULD BE WORSE THAN BUILDING IT FROM SCRATCH
    ! for the time being we use a pointer
    !
    fft % dfft => dfft
    !
    fft % cell => cell
    fft % gcutm = gcutm
    fft % ngm = ngm
    fft % gstart = gstart
    ALLOCATE( fft % gg( ngm ) )
    fft % gg = gg
    ALLOCATE( fft % g( 3, ngm ) )
    fft % g = g
    !
    IF ( fft % use_internal_pbc_corr ) ALLOCATE( fft % mt_corr( ngm ) )
    !
    ! In the future we will need to initialize things from scratch
    !
    ! CALL fft_type_init( dfft, smap, "environ", .TRUE., .TRUE., comm, cell%at, cell%bg, gcutm, 4.D0, fft_fact, nyfft )
    !
    ! The following routines are in tools_generate_gvect and may need to be simplified
    !
    ! CALL gvect_init( ngm, comm )
    ! CALL ggen( dfft, .TRUE. , cell%at, cell%bg,  gcutm, ngm_g, ngm, g, gg, mill, ig_l2g, gstart, .TRUE. )
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE init_fft_core_second
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE update_fft_core_cell( cell, fft )
!--------------------------------------------------------------------
    !
    USE correction_mt, ONLY : update_mt_correction
    !
    IMPLICIT NONE
    !
    TYPE( environ_cell ), TARGET, INTENT(IN) :: cell
    TYPE( fft_core ), INTENT(INOUT) :: fft
    !
    fft % cell => cell
    IF ( fft % use_internal_pbc_corr) CALL update_mt_correction( fft )
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE update_fft_core_cell
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE destroy_fft_core( lflag, fft )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    LOGICAL, INTENT(IN) :: lflag
    TYPE( fft_core ), INTENT(INOUT) :: fft
    !
    CHARACTER( LEN = 80 ) :: sub_name = 'destroy_fft_core'
    !
    NULLIFY( fft%dfft )
    NULLIFY( fft%cell )
    DEALLOCATE( fft % gg )
    DEALLOCATE( fft % g )
    IF ( fft % use_internal_pbc_corr ) DEALLOCATE( fft % mt_corr )
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE destroy_fft_core
!--------------------------------------------------------------------
END MODULE utils_fft
