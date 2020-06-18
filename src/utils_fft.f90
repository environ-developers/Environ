MODULE utils_fft
  !
  USE modules_constants, ONLY : DP, tpi
  USE core_types
  !
  PRIVATE
  !
  PUBLIC :: create_fft_core, init_fft_core_first, init_fft_core_second, &
       update_fft_core_cell, destroy_fft_core, init_dfft_core, init_dfft_core_second
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
  SUBROUTINE init_fft_core_second( cell, gcutm, ngm, dfft, fft )
!--------------------------------------------------------------------
    !
    USE fft_types,               ONLY : fft_type_init
    USE mp_bands,                ONLY : nyfft
    USE tools_generate_gvectors, ONLY : env_gvect_init, env_ggen
    USE correction_mt,           ONLY : update_mt_correction
    IMPLICIT NONE
    !
    TYPE( environ_cell ), TARGET, INTENT(IN) :: cell
    TYPE( fft_core ), INTENT(INOUT) :: fft
    TYPE( fft_type_descriptor ), TARGET, INTENT(IN) :: dfft
    INTEGER :: fft_fact(3)
    INTEGER :: i, ngm_g
    INTEGER, INTENT(IN) :: ngm
    REAL(DP) :: gcutm
    !
    fft%cell => cell
    !
    fft%gcutm = gcutm
    fft%ngm = ngm
    fft%dfft => dfft
    !
    IF ( fft % use_internal_pbc_corr ) ALLOCATE( fft % mt_corr( ngm ) )
    !
    ! The following routines are in tools_generate_gvect and may need to be simplified
    !
    CALL env_gvect_init( fft, cell%comm )
    CALL env_ggen( fft%dfft, cell%comm, cell%at, cell%bg, fft%gcutm, &
     & ngm_g, fft%ngm, fft%g, fft%gg, fft%gstart, .TRUE. )
    !
    IF ( fft % use_internal_pbc_corr) CALL update_mt_correction( fft )
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE init_fft_core_second
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE init_dfft_core( cell, gcutm, comm, at )
!--------------------------------------------------------------------
    !
    USE modules_constants, ONLY: pi
    USE stick_base,        ONLY: sticks_map
    USE fft_types,         ONLY: fft_type_init
    USE mp_bands,          ONLY: nyfft
    USE core_base,         ONLY: sys_dfft
    IMPLICIT NONE
    !
    TYPE( environ_cell ), INTENT(INOUT) :: cell
    TYPE( sticks_map ) :: smap
    INTEGER, INTENT(IN) :: comm
    REAL(DP), INTENT(IN) :: gcutm, at(3,3)
    !
    ! Needed some values from environ_cell type for fft_type_init
    !
    cell % comm = comm
    cell % at = at
    !
    ! Calculate the reciprocal lattice vectors
    !
    CALL recips( cell%at(1,1), cell%at(1,2), cell%at(1,3), cell%bg(1,1), &
      & cell%bg(1,2), cell%bg(1,3))
    !
    sys_dfft%rho_clock_label='fft'
    sys_dfft%lgamma = .TRUE.
    CALL fft_type_init( sys_dfft, smap, "rho", sys_dfft%lgamma, .TRUE., cell%comm, cell%at, &
         & cell%bg, gcutm, 4.D0, nyfft=nyfft )
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE init_dfft_core
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE init_dfft_core_second( cell, gcutm, comm, at, m )
!--------------------------------------------------------------------
    !
    USE modules_constants, ONLY: pi
    USE stick_base,        ONLY: sticks_map
    USE fft_types,         ONLY: fft_type_init
    USE mp_bands,          ONLY: nyfft
    USE core_base,         ONLY: environment_dfft!, sys_dfft
    !USE environ_base,      ONLY: sys_cell
    IMPLICIT NONE
    !
    TYPE( environ_cell ), INTENT(INOUT) :: cell
    TYPE( sticks_map ) :: smap
    INTEGER, INTENT(IN) :: comm
    INTEGER, INTENT(IN) :: m(3)
    REAL(DP), INTENT(IN) :: gcutm, at(3,3)
    INTEGER :: i
    !
    ! Needed some values from environ_cell type for fft_type_init
    !
    ! If m values equals 1 then need to point to sys_dfft and sys_cell
    !IF ( m(1)*m(2)*m(3) == 1) THEN
    !  environment_dfft => sys_dfft
    !  cell => sys_cell
    !  RETURN
    !ENDIF
    cell % comm = comm
    DO i=1,3
      cell % at(:,i) = at(:,i)*m(i)
    END DO
    !
    ! Calculate the reciprocal lattice vectors
    !
    CALL recips( cell%at(1,1), cell%at(1,2), cell%at(1,3), cell%bg(1,1), &
      & cell%bg(1,2), cell%bg(1,3))
    !
    environment_dfft%rho_clock_label='fft'
    environment_dfft%lgamma = .TRUE.
    CALL fft_type_init( environment_dfft, smap, "rho", environment_dfft%lgamma, .TRUE., &
         & cell%comm, cell%at, cell%bg, gcutm, 4.D0, nyfft=nyfft )
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE init_dfft_core_second
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
