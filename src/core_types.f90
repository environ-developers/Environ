MODULE core_types
  !
  USE modules_constants, ONLY : DP
  USE fft_types, ONLY : fft_type_descriptor
  !
  TYPE environ_cell
     !
     ! Global properties of the simulation cell
     !
     LOGICAL :: update = .FALSE.
     INTEGER :: ibrav
     REAL( DP ) :: alat
     REAL( DP ) :: omega
     REAL( DP ) :: domega
     REAL( DP ) :: origin( 3 )
     REAL( DP ), DIMENSION( 3, 3 ) :: at
     REAL( DP ), DIMENSION( 3, 3 ) :: bg
     REAL( DP ), DIMENSION( 3, 8 ) :: corners
     !
     ! Properties of the grid
     !
     INTEGER :: ntot, n1, n2, n3
     REAL( DP ) :: in1, in2, in3
     INTEGER :: n1x, n2x, n3x
! BACKWARD COMPATIBILITY
! Compatible with QE-5.X QE-6.0.X QE-6.1.X
!     INTEGER :: idx0
! Compatible with QE-6.2.X QE-6.3.X QE-6.4.X QE-GIT
     INTEGER :: j0, k0 , n2p, n3p
! END BACKWARD COMPATIBILITY
     !
     ! Properties of the processor-specific partition
     !
     INTEGER :: nnr    ! size of processor-specific allocated fields
     INTEGER :: ir_end ! actual physical size of processor-specific allocated field
     INTEGER :: comm   ! parallel communicator
     INTEGER :: me     ! index of processor
     INTEGER :: root   ! index of root
     !
  END TYPE environ_cell
  !
  TYPE fd_core
     !
     INTEGER :: ifdtype
     INTEGER :: nfdpoint
     INTEGER, ALLOCATABLE :: icfd(:)
     INTEGER :: ncfd
     !
     TYPE(fft_type_descriptor), POINTER :: dfft
     !
     TYPE(environ_cell), POINTER :: cell
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
     TYPE(fft_type_descriptor), POINTER :: dfft
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
     REAL(DP), ALLOCATABLE :: gg(:)
     !
     !     G-vectors cartesian components ( in units tpiba =(2pi/a)  )
     !
     REAL(DP), ALLOCATABLE :: g(:,:)
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
  SUBROUTINE init_environ_cell( n1, n2, n3, ibrav, alat, omega, at, bg, &
       & nnr, ir_end, n1x, n2x, n3x, &
! BACKWARD COMPATIBILITY
! Compatible with QE-5.X QE-6.0.X QE-6.1.X
!         & idx0, &
! Compatible with QE-6.2.X QE-6.3.X QE-6.4.X QE-GIT
       & j0, k0, n2p, n3p, &
! END BACKWARD COMPATIBILITY
       & comm, me, root, cell )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: n1, n2, n3, ibrav
    INTEGER, INTENT(IN) :: n1x, n2x, n3x
! BACKWARD COMPATIBILITY
! Compatible with QE-5.X QE-6.0.X QE-6.1.X
!    INTEGER, INTENT(IN) :: idx0
! Compatible with QE-6.2.X QE-6.3.X QE-6.4.X QE-GIT
    INTEGER, INTENT(IN) :: j0, k0, n2p, n3p
! END BACKWARD COMPATIBILITY
    INTEGER, INTENT(IN) :: nnr, ir_end, comm, me, root
    REAL( DP ), INTENT(IN) :: alat, omega, at(3,3), bg(3,3)
    TYPE( environ_cell ), INTENT(INOUT) :: cell
    !
    INTEGER :: ic, ix, iy, iz
    REAL( DP ) :: dx, dy, dz
    !
    CHARACTER( LEN=80 ) :: sub_name = 'init_environ_cell'
    !
    IF ( n1 .EQ. 0 .OR. n2 .EQ. 0 .OR. n3 .EQ. 0 ) &
         & CALL errore(sub_name,'Wrong grid dimension',1)
    !
    cell % n1 = n1
    cell % n2 = n2
    cell % n3 = n3
    cell % ibrav = ibrav
    cell % alat = alat
    cell % omega = omega
    cell % at = at
    cell % bg = bg
    !
    cell % in1 = 1.D0 / DBLE(n1)
    cell % in2 = 1.D0 / DBLE(n2)
    cell % in3 = 1.D0 / DBLE(n3)
    cell % n1x = n1x
    cell % n2x = n2x
    cell % n3x = n3x
! BACKWARD COMPATIBILITY
! Compatible with QE-5.X QE-6.0.X QE-6.1.X
!    cell % idx0 = idx0
! Compatible with QE-6.2.X QE-6.3.X QE-6.4.X QE-GIT
    cell % j0 = j0
    cell % k0 = k0
    cell % n2p = n2p
    cell % n3p = n3p
! END BACKWARD COMPATIBILITY
    !
    cell % nnr = nnr
    cell % ir_end = ir_end
    cell % comm = comm
    cell % me   = me
    cell % root = root
    !
    cell % ntot = cell % n1 * cell % n2 * cell % n3
    cell % domega = cell % omega / cell % ntot
    !
    cell % origin = 0.D0
    !
    ic = 0
    DO ix = 0,1
       !
       dx = DBLE(-ix)
       !
       DO iy = 0,1
          !
          dy = DBLE(-iy)
          !
          DO iz = 0,1
             !
             dz = DBLE(-iz)
             !
             ic = ic + 1
             cell%corners(1,ic) = dx*at(1,1) + dy*at(1,2) + dz*at(1,3)
             cell%corners(2,ic) = dx*at(2,1) + dy*at(2,2) + dz*at(2,3)
             cell%corners(3,ic) = dx*at(3,1) + dy*at(3,2) + dz*at(3,3)
             !
          ENDDO
          !
       ENDDO
       !
    ENDDO
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE init_environ_cell
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE update_environ_cell( omega, at, cell )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    REAL( DP ), INTENT(IN) :: omega, at(3,3)
    TYPE( environ_cell ), INTENT(INOUT) :: cell
    CHARACTER( LEN=80 ) :: sub_name = 'update_environ_cell'
    !
    cell % omega = omega
    cell % at = at
    !
    cell % domega = cell % omega / cell % ntot
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE update_environ_cell
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE ir2ijk( cell, ir, i, j, k, physical )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    TYPE( environ_cell ), INTENT(IN) :: cell
    INTEGER, INTENT(IN) :: ir
    INTEGER, INTENT(OUT) :: i, j, k
    LOGICAL, INTENT(OUT) :: physical
    !
    INTEGER :: idx
    !
! BACKWARD COMPATIBILITY
! Compatible with QE-5.X QE-6.0.X QE-6.1.X
!     idx = cell%idx0 + ir - 1
!     k   = idx / (cell%n1x*cell%n2x)
!     idx = idx - (cell%n1x*cell%n2x)*k
!     j   = idx / cell%n1x
!     idx = idx - cell%n1x*j
!     i   = idx
! Compatible with QE-6.2.X QE-6.3.X QE-6.4.X QE-GIT
    idx = ir - 1
    k   = idx / (cell%n1x*cell%n2p)
    idx = idx - (cell%n1x*cell%n2p)*k
    k   = k + cell%k0
    j   = idx / cell%n1x
    idx = idx - cell%n1x * j
    j   = j + cell%j0
    i   = idx
! END BACKWARD COMPATIBILITY
    !
    physical = i < cell%n1 .AND. j < cell%n2 .AND. k < cell%n3
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE ir2ijk
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE ir2r( cell, ir, r, physical )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    TYPE( environ_cell ), INTENT(IN) :: cell
    INTEGER, INTENT(IN) :: ir
    REAL( DP ), DIMENSION(3), INTENT(OUT) :: r
    LOGICAL, INTENT(OUT) :: physical
    !
    INTEGER :: idx, i, j, k, ip
    !
    r = 0.D0
    !
    CALL ir2ijk( cell, ir, i, j, k, physical )
    !
    IF ( .NOT. physical ) RETURN
    !
    DO ip = 1, 3
       r(ip) = DBLE( i ) * cell%in1 * cell%at(ip,1) + &
               DBLE( j ) * cell%in2 * cell%at(ip,2) + &
               DBLE( k ) * cell%in3 * cell%at(ip,3)
    END DO
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE ir2r
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE displacement( dim, axis, r1, r2, dr )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: dim, axis
    REAL( DP ), DIMENSION(3), INTENT(IN) :: r1, r2
    REAL( DP ), DIMENSION(3), INTENT(OUT) :: dr
    !
    INTEGER :: i
    !
    dr(:) = r1(:) - r2(:)
    !
    !  ... possibly only in 1D or 2D
    !
    SELECT CASE ( dim )
    CASE ( 1 )
       dr(axis) = 0.D0
    CASE ( 2 )
       DO i = 1, 3
          IF ( i .NE. axis ) dr(i) = 0.D0
       ENDDO
    END SELECT
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE displacement
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE minimum_image( cell, r, r2 )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    TYPE( environ_cell ), INTENT(IN) :: cell
    REAL( DP ), DIMENSION(3), INTENT(INOUT) :: r
    REAL( DP ), INTENT(OUT) :: r2
    !
    INTEGER :: ic
    REAL( DP ), DIMENSION(3) :: s
    REAL( DP ), DIMENSION(3) :: rmin
    REAL( DP ) :: r2min
    !
    s(:) = MATMUL( r(:), cell%bg(:,:) )
    s(:) = s(:) - FLOOR(s(:))
    r(:) = MATMUL( cell%at(:,:), s(:) )
    !
    rmin = r
    r2min = SUM( r * r )
    DO ic = 2,8
       s(1) = r(1) + cell%corners(1,ic)
       s(2) = r(2) + cell%corners(2,ic)
       s(3) = r(3) + cell%corners(3,ic)
       r2 = SUM( s * s )
       IF (r2<r2min) THEN
          rmin = s
          r2min = r2
       END IF
    ENDDO
    r = rmin
    r2 = r2min
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE minimum_image
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE create_fd_core( fd )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    TYPE( fd_core ), INTENT(INOUT) :: fd
    !
    ! Create empty finite difference core
    !
    NULLIFY( fd%dfft )
    NULLIFY( fd%cell )
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE create_fd_core
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE init_fd_core_first( ifdtype, nfdpoint, fd )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: ifdtype, nfdpoint
    TYPE( fd_core ), INTENT(INOUT) :: fd
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
    TYPE( fd_core ), TARGET, INTENT(INOUT) :: fd
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
  SUBROUTINE init_fd_core_second( cell, dfft, fd )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    TYPE( environ_cell ), TARGET, INTENT(IN) :: cell
    TYPE( fft_type_descriptor ), TARGET, INTENT(IN) :: dfft
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
  SUBROUTINE init_fft_core_second( dfft, omega, tpiba, tpiba2, ngm, &
       & gcutm, gstart, g, gg, fft )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    TYPE( fft_type_descriptor ), TARGET, INTENT(IN) :: dfft
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
    fft % omega = omega ! Calculated
    fft % gcutm = gcutm !Calculated
    fft % tpiba = tpiba !Calculated
    fft % tpiba2 = tpiba2 !Calculated
    fft % ngm = ngm !Passed from PW
    fft % gstart = gstart !Passed from PW
    ALLOCATE( fft % gg( ngm ) )
    fft % gg = gg
    ALLOCATE( fft % g( 3, ngm ) )
    fft % g = g
    print *, 'Printing information'
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
  SUBROUTINE init_fft_core_second_new( cell, ecutrho, ngm, gstart, dfft, g, gg, fft )
!--------------------------------------------------------------------
    !
    USE modules_constants, ONLY: pi
    USE stick_base,        ONLY: sticks_map
    USE fft_types,         ONLY: fft_type_init
    IMPLICIT NONE
    !
    TYPE( environ_cell ), INTENT(IN) :: cell
    TYPE( fft_core ), INTENT(INOUT) :: fft
    TYPE( fft_type_descriptor ), INTENT(IN) :: dfft
    REAL( DP ), DIMENSION( ngm ), INTENT(IN) :: gg
    REAL( DP ), DIMENSION( 3, ngm ), INTENT(IN) :: g
    TYPE( sticks_map ) :: smap
    INTEGER :: fft_fact(3)
    INTEGER :: nyfft
    INTEGER, INTENT(IN) :: ngm, gstart
    REAL(DP) :: ecutrho
    !
    fft%tpiba = 2.D0 * pi / cell%alat
    fft%tpiba2 = fft%tpiba**2.D0
    fft%gcutm = ecutrho / fft%tpiba2
    fft%ngm = ngm
    fft%gstart = gstart
    !
    ! Should fft%omega just be cell%omega? Do we need both environ_cell and fft_core to have
    ! omega?
    CALL volume( cell%alat, cell%at(1,1), cell%at(1,2), cell%at(1,3), fft%omega )
    !
    ! recips calculates the reciprocal lattice vectors
    !
    CALL recips( cell%at(1,1), cell%at(1,2), cell%at(1,3), cell%bg(1,1), &
      & cell%bg(1,2), cell%bg(1,3))
    ALLOCATE( fft % gg( ngm ) )
    fft % gg = gg
    ALLOCATE( fft % g( 3, ngm ) )
    fft % g = g
    print *, 'Printing information for new'
    !
    !
    !CALL fft_type_init( fft%dfft, smap, "environ", .TRUE., .TRUE., cell%comm, cell%at, &
    !  & cell%bg, fft%gcutm, 4.D0, fft_fact, nyfft )
    !
    ! The following routines are in tools_generate_gvect and may need to be simplified
    !
    !CALL env_gvect_init( ngm, cell%comm )
    !CALL ggen( fft, .TRUE. , cell%at, cell%bg, fft%gcutm, fft%ngm, fft%g, fft%gg, &
    !  & fft%gstart )
    !
    RETURN
    !
!--------------------------------------------------------------------
END SUBROUTINE init_fft_core_second_new
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
  !
!---------------------------------------------------------------------
SUBROUTINE volume (alat, a1, a2, a3, omega)
   !---------------------------------------------------------------------
   !
   !     Compute the volume of the unit cell defined by 3 vectors
   !     a1, a2, a3, given in units of "alat" (alat may be 1):
   !        omega = alat^3 * [ a1 . (a2 x a3) ]
   !     ( . = scalar product, x = vector product )
   !
   USE modules_constants, ONLY: dp
   IMPLICIT NONE
   !
   REAL(dp), INTENT(IN) :: alat, a1(3), a2(3), a3(3)
   REAL(dp), INTENT(OUT) :: omega
   !
   omega = a1(1) * ( a2(2)*a3(3)-a2(3)*a3(2) ) - &
           a1(2) * ( a2(1)*a3(3)-a2(3)*a3(1) ) + &
           a1(3) * ( a2(1)*a3(2)-a2(2)*a3(1) )
   !
   IF ( omega < 0.0_dp) THEN
      call infomsg('volume','axis vectors are left-handed')
      omega = ABS (omega)
   END IF
   !
   IF ( alat < 1.0_dp) call infomsg('volume','strange lattice parameter')
   omega = omega * alat**3
   !
   RETURN
   !
 END SUBROUTINE volume
 !---------------------------------------------------------------------
 subroutine recips (a1, a2, a3, b1, b2, b3)
   !---------------------------------------------------------------------
   !
   !   This routine generates the reciprocal lattice vectors b1,b2,b3
   !   given the real space vectors a1,a2,a3. The b's are units of 2 pi/a.
   !
   !     first the input variables
   !
   use modules_constants, ONLY: DP
   implicit none
   real(DP) :: a1 (3), a2 (3), a3 (3), b1 (3), b2 (3), b3 (3)
   ! input: first direct lattice vector
   ! input: second direct lattice vector
   ! input: third direct lattice vector
   ! output: first reciprocal lattice vector
   ! output: second reciprocal lattice vector
   ! output: third reciprocal lattice vector
   !
   !   then the local variables
   !
   real(DP) :: den, s
   ! the denominator
   ! the sign of the permutations
   integer :: iperm, i, j, k, l, ipol
   ! counter on the permutations
   !\
   !  Auxiliary variables
   !/
   !
   ! Counter on the polarizations
   !
   !    first we compute the denominator
   !
   den = 0
   i = 1
   j = 2
   k = 3
   s = 1.d0
 100 do iperm = 1, 3
      den = den + s * a1 (i) * a2 (j) * a3 (k)
      l = i
      i = j
      j = k
      k = l
   enddo
   i = 2
   j = 1
   k = 3
   s = - s
   if (s.lt.0.d0) goto 100
   !
   !    here we compute the reciprocal vectors
   !
   i = 1
   j = 2
   k = 3
   do ipol = 1, 3
      b1 (ipol) = (a2 (j) * a3 (k) - a2 (k) * a3 (j) ) / den
      b2 (ipol) = (a3 (j) * a1 (k) - a3 (k) * a1 (j) ) / den
      b3 (ipol) = (a1 (j) * a2 (k) - a1 (k) * a2 (j) ) / den
      l = i
      i = j
      j = k
      k = l
   enddo
   return
 end subroutine recips
 
END MODULE core_types
