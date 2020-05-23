MODULE core_types
  !
  USE environ_types
  !
  USE fft_types, ONLY : fft_type_descriptor
  !
  TYPE fd_core
     !
     INTEGER :: ifdtype
     INTEGER :: nfdpoint
     INTEGER, ALLOCATABLE :: icfd(:)
     INTEGER :: ncfd
     !
  END TYPE fd_core
  !
  TYPE fft_core
     !
     INTEGER :: index
! BACKWARD COMPATIBILITY
! Compatible with QE-6.0 QE-6.1.X QE-6.2.X QE-6.3.X
!     INTEGER :: nspin
! Compatible with QE-6.4.X QE-GIT
!
! END BACKWARD COMPATIBILITY
     LOGICAL :: use_internal_pbc_corr = .FALSE.
     !
     TYPE(fft_descriptor_type), POINTER :: dfft
     !
     REAL(DP) :: omega, tpiba, tpiba2
     !
     INTEGER :: ngm  = 0  ! local  number of G vectors (on this processor)
                          ! with gamma tricks, only vectors in G>
     !
     REAL(DP) :: gcutm = 0.0_DP   ! ecutrho/(2 pi/a)^2, cut-off for |G|^2

     INTEGER :: gstart = 2 ! index of the first G vector whose module is > 0
                           ! Needed in parallel execution: gstart=2 for the
                           ! proc that holds G=0, gstart=1 for all others
     !
     !     G^2 in increasing order (in units of tpiba2=(2pi/a)^2)
     !
     REAL(DP), ALLOCATABLE, TARGET :: gg(:)
     !
     !     G-vectors cartesian components ( in units tpiba =(2pi/a)  )
     !
     REAL(DP), ALLOCATABLE, TARGET :: g(:,:)
     !
  END TYPE fft_core
  !
  TYPE oned_analytic_core
     !
     INTEGER :: n, d, p, axis
     LOGICAL :: initialized = .FALSE.
     REAL( DP ) :: size
     TYPE( environ_cell ), POINTER :: cell
     REAL( DP ), DIMENSION(3) :: origin
     REAL( DP ), DIMENSION(:,:), ALLOCATABLE :: x
     !
  END TYPE oned_analytic_core
  !
CONTAINS
!--------------------------------------------------------------------
  SUBROUTINE init_fd_core_first( ifdtype, nfdpoint, fd )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: ifdtype, nfdpoint
    TYPE( fd_core ), INTENT(OUT) :: fd
    !
    !
    ! Set finite differences tools
    !
    fd % ifdtype = ifdtype
    fd % nfdpoint = nfdpoint
    ALLOCATE( fd%icfd(-nfdpoint:nfdpoint) )
    !
    CALL set_fd_coefficients( fd )
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE init_fd_core_first
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE set_fd_coefficients( fd )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    TYPE( fd_core ), INTENT(INOUT) :: fd
    !
    INTEGER, POINTER :: ifdtype, nfdpoint, ncfd
    INTEGER, DIMENSION(:), POINTER :: icfd
    !
    INTEGER :: in
    !
    ifdtype => fd % ifdtype
    nfdpoint => fd % nfdpoint
    ncfd => fd % ncfd
    icfd => fd % icfd
    !
    ncfd = 0
    icfd = 0
    !
    SELECT CASE ( ifdtype )
       !
    CASE ( 1 )
       ! (2N+1)-point Central Differences
       IF ( nfdpoint .EQ. 1 ) THEN
          ncfd = 2
          icfd(  1 ) =   1
       ELSE IF ( nfdpoint .EQ. 2 ) THEN
          ncfd = 12
          icfd(  2 ) =  -1
          icfd(  1 ) =   8
       ELSE IF ( nfdpoint .EQ. 3 ) THEN
          ncfd = 60
          icfd(  3 ) =   1
          icfd(  2 ) =  -9
          icfd(  1 ) =  45
       ELSE IF ( nfdpoint .EQ. 4 ) THEN
          ncfd = 840
          icfd(  4 ) =  -3
          icfd(  3 ) =  32
          icfd(  2 ) =-168
          icfd(  1 ) = 672
       ELSE
          WRITE(*,*)'ERROR: wrong number of points',nfdpoint,&
               &' for finite difference type ',ifdtype
          STOP
       ENDIF
       !
    CASE ( 2 )
       ! Low-Noise Lanczos Differentiators ( M = 2 )
       IF ( nfdpoint .GE. 2 ) THEN
          ncfd = (nfdpoint)*(nfdpoint+1)*(2*nfdpoint+1)/3
          DO in = 1,nfdpoint
             icfd( in ) = in
          ENDDO
       ELSE
          WRITE(*,*)'ERROR: wrong number of points',nfdpoint,&
               &' for finite difference type ',ifdtype
          STOP
       END IF
       !
    CASE ( 3 )
       ! Super Lanczos Low-Noise Differentiators ( M = 4 )
       IF ( nfdpoint .EQ. 3 ) THEN
          ncfd = 252
          icfd(  3 ) = -22
          icfd(  2 ) =  67
          icfd(  1 ) =  58
       ELSE IF ( nfdpoint .EQ. 4 ) THEN
          ncfd = 1188
          icfd(  4 ) = -86
          icfd(  3 ) = 142
          icfd(  2 ) = 193
          icfd(  1 ) = 126
       ELSE IF ( nfdpoint .EQ. 5 ) THEN
          ncfd = 5148
          icfd(  5 ) =-300
          icfd(  4 ) = 294
          icfd(  3 ) = 532
          icfd(  2 ) = 503
          icfd(  1 ) = 296
       ELSE
          WRITE(*,*)'ERROR: wrong number of points',nfdpoint,&
               &' for finite difference type ',ifdtype
          STOP
       ENDIF
       !
    CASE ( 4 )
       ! Smooth Noise-Robust Differentiators  ( n = 2 )
       IF ( nfdpoint .EQ. 2 ) THEN
          ncfd = 8
          icfd(  2 ) =   1
          icfd(  1 ) =   2
       ELSE IF ( nfdpoint .EQ. 3 ) THEN
          ncfd = 32
          icfd(  3 ) =   1
          icfd(  2 ) =   4
          icfd(  1 ) =   5
       ELSE IF ( nfdpoint .EQ. 4 ) THEN
          ncfd = 128
          icfd(  4 ) =   1
          icfd(  3 ) =   6
          icfd(  2 ) =  14
          icfd(  1 ) =  14
       ELSE IF ( nfdpoint .EQ. 5 ) THEN
          ncfd = 512
          icfd(  5 ) =   1
          icfd(  4 ) =   8
          icfd(  3 ) =  27
          icfd(  2 ) =  48
          icfd(  1 ) =  42
       ELSE
          WRITE(*,*)'ERROR: wrong number of points',nfdpoint,&
               &' for finite difference type ',ifdtype
          STOP
       ENDIF
       !
    CASE ( 5 )
       ! Smooth Noise-Robust Differentiators  ( n = 4 )
       IF ( nfdpoint .EQ. 3 ) THEN
          ncfd = 96
          icfd(  3 ) =  -5
          icfd(  2 ) =  12
          icfd(  1 ) =  39
       ELSE IF ( nfdpoint .EQ. 4 ) THEN
          ncfd = 96
          icfd(  4 ) =  -2
          icfd(  3 ) =  -1
          icfd(  2 ) =  16
          icfd(  1 ) =  27
       ELSE IF ( nfdpoint .EQ. 5 ) THEN
          ncfd = 1536
          icfd(  5 ) = -11
          icfd(  4 ) = -32
          icfd(  3 ) =  39
          icfd(  2 ) = 256
          icfd(  1 ) = 322
       ELSE
          WRITE(*,*)'ERROR: wrong number of points',nfdpoint,&
               &' for finite difference type ',ifdtype
          STOP
       ENDIF
       !
    CASE DEFAULT
       !
       WRITE(*,*)'ERROR: finite difference type unknown, ifdtype=',ifdtype
       STOP
       !
    END SELECT
    !
    DO in = 1,nfdpoint
       icfd( -in ) = - icfd( in )
    ENDDO
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE set_fd_coefficients
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE init_fd_core_second( dfft, cell, fd )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    TYPE( fft_dlay_descriptor ), TARGET, INTENT(IN) :: dfft
    TYPE( environ_cell ), TARGET, INTENT(IN) :: cell
    TYPE( fd_core ), INTENT(INOUT) :: fd
    !
    fd % dfft => dfft
    fd % cell => cell
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE init_fd_core_second
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE destroy_fd_core( lflag, fd )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    LOGICAL, INTENT(IN) :: lflag
    TYPE( fd_core ), INTENT(INOUT) :: fd
    !
    CHARACTER( LEN = 80 ) :: sub_name = 'destroy_fd_core'
    !
    NULLIFY( fd%dfft )
    NULLIFY( fd%cell )
    !
    IF ( lflag ) THEN
       IF ( .NOT. ALLOCATED( fd%icfd ) ) &
            & CALL errore( sub_name, 'Trying to deallocate a non-allocated object', 1 )
       DEALLOCATE( fd % icfd )
    END IF
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE destroy_fd_core
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
  SUBROUTINE init_fft_core_second( dfft, omega, tpiba, tpiba2, ngm, &
       & gcutm, gstart, g, gg, fft )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    TYPE( fft_dlay_descriptor ), TARGET, INTENT(IN) :: dfft
    REAL( DP ), INTENT(IN) :: omega, gcutm, tpiba, tpiba2
    INTEGER, INTENT(IN) :: ngm, gstart
    REAL( DP ), DIMENSION( ngm ), INTENT(IN) :: gg
    REAL( DP ), DIMENSION( 3, ngm ), INTENT(IN) :: g
    TYPE( fft_core ), INTENT(INOUT) :: fft
    !
    ! COPYING ALL THE COMPONENTS OF DFFT WOULD BE WORSE THAN BUILDING IT FROM SCRATCH
    ! for the time being we use a pointer
    !
    fft % dfft => dfft
    !
    fft % omega = omega
    fft % gcutm = gcutm
    fft % tpiba = tpiba
    fft % tpiba2 = tpiba2
    fft % ngm = ngm
    fft % gstart = gstart
    ALLOCATE( fft % gg( ngm ) )
    fft % gg = gg
    ALLOCATE( fft % g( 3, ngm ) )
    fft % g = g
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
  SUBROUTINE update_fft_core_cell( omega, fft )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    REAL( DP ), INTENT(IN) :: omega
    TYPE( fft_core ), INTENT(INOUT) :: fft
    !
    fft % omega = omega
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
    DEALLOCATE( fft % gg )
    DEALLOCATE( fft % g )
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE destroy_fft_core
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE init_oned_analytic_core_first( dim, axis, oned_analytic )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: dim, axis
    TYPE( oned_analytic_core ), INTENT(OUT) :: oned_analytic
    !
    CHARACTER( LEN = 80 ) :: sub_name = 'init_oned_analytic_core_first'
    !
    IF ( dim .EQ. 3 .OR. dim .LT. 0 ) &
         & CALL errore(sub_name,'Wrong dimensions for analytic one dimensional core',1)
    oned_analytic % d = dim
    oned_analytic % p = 3 - dim
    !
    IF ( ( dim .EQ. 1 .OR. dim .EQ. 2 ) .AND. ( axis .GT. 3 .OR. axis .LT. 1 ) ) &
         & CALL errore(sub_name,'Wrong choice of axis for analytic one dimensional core',1)
    oned_analytic % axis = axis
    !
    oned_analytic % initialized = .FALSE.
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE init_oned_analytic_core_first
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE init_oned_analytic_core_second( cell, oned_analytic )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    TYPE( environ_cell ), INTENT(IN) :: cell
    TYPE( oned_analytic_core ), INTENT(INOUT) :: oned_analytic
    !
    CHARACTER( LEN = 80 ) :: sub_name = 'init_oned_analytic_core_second'
    !
    CALL update_oned_analytic_core_cell( cell, oned_analytic )
    !
    oned_analytic % n = cell % nnr
    ALLOCATE( oned_analytic % x( oned_analytic % p , oned_analytic % n ) )
    !
    CALL update_oned_analytic_core_origin( cell%origin, oned_analytic )
    !
    oned_analytic % initialized = .TRUE.
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE init_oned_analytic_core_second
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE update_oned_analytic_core_cell( cell, oned_analytic )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    TYPE( environ_cell ), TARGET, INTENT(IN) :: cell
    TYPE( oned_analytic_core ), INTENT(INOUT) :: oned_analytic
    !
    oned_analytic % cell => cell
    IF ( oned_analytic % d .EQ. 0 ) THEN
       oned_analytic % size = cell % omega
    ELSE IF ( oned_analytic % d .EQ. 1 ) THEN
       oned_analytic % size = cell % omega / cell % at( oned_analytic % axis, oned_analytic % axis ) / cell % alat
    ELSE IF ( oned_analytic % d .EQ. 2 ) THEN
       oned_analytic % size = cell % at( oned_analytic % axis, oned_analytic % axis ) * cell % alat
    ENDIF
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE update_oned_analytic_core_cell
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE update_oned_analytic_core_origin( origin, oned_analytic )
!--------------------------------------------------------------------
    !
    USE tools_generate_functions, ONLY : generate_axis, generate_distance
    !
    IMPLICIT NONE
    !
    REAL( DP ), DIMENSION( 3 ), INTENT(IN) :: origin
    TYPE( oned_analytic_core ), INTENT(INOUT) :: oned_analytic
    !
    CHARACTER( LEN = 80 ) :: sub_name = 'update_oned_analytic_core_origin'
    !
    oned_analytic % origin = origin
    IF ( oned_analytic % d .EQ. 0 ) THEN
       CALL generate_distance( oned_analytic % cell, oned_analytic % origin, oned_analytic % x )
    ELSE IF ( oned_analytic % d .EQ. 1 ) THEN
       CALL errore( sub_name, 'Option not yet implemented', 1 )
    ELSE IF ( oned_analytic % d .EQ. 2 ) THEN
       CALL generate_axis( oned_analytic % cell, oned_analytic % axis, oned_analytic % origin, oned_analytic % x(1,:) )
    ENDIF
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE update_oned_analytic_core_origin
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE destroy_oned_analytic_core( lflag, oned_analytic )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    LOGICAL, INTENT(IN) :: lflag
    TYPE( oned_analytic_core ), INTENT(INOUT) :: oned_analytic
    !
    CHARACTER( LEN = 80 ) :: sub_name = 'destroy_oned_analytic_core'
    !
    IF (oned_analytic % initialized ) THEN
      IF ( .NOT. ALLOCATED( oned_analytic % x ) ) &
           & CALL errore(sub_name,'Trying to destroy a non-allocated component',1)
      DEALLOCATE( oned_analytic % x )
      IF ( .NOT. ASSOCIATED( oned_analytic % cell ) ) &
           & CALL errore(sub_name,'Trying to nullify a non-associated pointer',1)
      NULLIFY(oned_analytic%cell)
      oned_analytic % initialized = .FALSE.
    END IF
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE destroy_oned_analytic_core
!--------------------------------------------------------------------
  !
END MODULE core_types
