MODULE cell_types
  !
  USE modules_constants, ONLY : DP, tpi
  !
  TYPE environ_cell
     !
     ! Global properties of the simulation cell
     !
     LOGICAL :: update = .FALSE.
     LOGICAL :: cubic = .FALSE.
     REAL( DP ) :: alat
     REAL( DP ) :: omega
     REAL( DP ) :: domega
     REAL( DP ) :: tpiba, tpiba2
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
     !
     ! Properties of the processor-specific partition
     !
     INTEGER :: j0, k0 , n2p, n3p
     INTEGER :: nnr    ! size of processor-specific allocated fields
     INTEGER :: ir_end ! actual physical size of processor-specific allocated field
     INTEGER :: comm   ! parallel communicator
     INTEGER :: me     ! index of processor
     INTEGER :: root   ! index of root
     !
  END TYPE environ_cell
  !
  TYPE environ_mapping
     !
     INTEGER, DIMENSION(3) :: nrep ! number of replicas of smaller cell in larger cell
     !
     TYPE( environ_cell ), POINTER :: small ! small cell
     !
     TYPE( environ_cell ), POINTER :: large ! large cell
     !
     INTEGER, DIMENSION(:), ALLOCATABLE :: map
     !
  END TYPE environ_mapping
  !
  PRIVATE
  !
  PUBLIC :: environ_cell, init_environ_cell, update_environ_cell, ir2ijk, ir2r, &
       displacement, minimum_image, volume, recips, &
       environ_mapping, init_environ_mapping_first, init_environ_mapping_second, &
       destroy_environ_mapping
  !
CONTAINS
!--------------------------------------------------------------------
  SUBROUTINE init_environ_cell( alat, at, comm, me, root, n1, n2, n3, &
       & n1x, n2x, n3x, n2p, n3p, j0, k0, nnr, cell )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: comm, me, root
    INTEGER, INTENT(IN) :: n1, n2, n3, n1x, n2x, n3x
    INTEGER, INTENT(IN) :: n2p, n3p, j0, k0, nnr
    REAL( DP ), INTENT(IN) :: alat, at(3,3)
    TYPE( environ_cell ), INTENT(INOUT) :: cell
    !
    CHARACTER( LEN=80 ) :: sub_name = 'init_environ_cell'
    !
    IF ( n1 .EQ. 0 .OR. n2 .EQ. 0 .OR. n3 .EQ. 0 ) &
         & CALL errore(sub_name,'Wrong grid dimension',1)
    !
    IF ( alat .LT. 1.D-8 ) &
         & CALL errore(sub_name,'Wrong alat',1)
    !
    ! Store cell units
    !
    cell % alat = alat
    !
    ! Calculate units of reciprocal space
    !
    cell % tpiba = tpi/alat
    cell % tpiba2 = cell % tpiba**2
    !
    ! Real space grid, global dimensions
    !
    cell % n1 = n1
    cell % n2 = n2
    cell % n3 = n3
    !
    cell % in1 = 1.D0 / DBLE(n1)
    cell % in2 = 1.D0 / DBLE(n2)
    cell % in3 = 1.D0 / DBLE(n3)
    !
    cell % n1x = n1x
    cell % n2x = n2x
    cell % n3x = n3x
    !
    ! Real space grid, local dimensions (processor-specific)
    !
    cell % n2p = n2p
    cell % n3p = n3p
    !
#if defined (__MPI)
    cell % j0 = j0 ; cell % k0 = k0
    cell % ir_end = MIN(nnr,n1x*n2p*n3p)
#else
    cell % j0 = 0; k0 = 0;
    cell % ir_end = nnr
#endif
    !
    cell % nnr = nnr
    !
    ! Cell parallelization
    !
    cell % comm = comm
    cell % me   = me
    cell % root = root
    !
    ! Total number of physical points
    !
    cell % ntot = cell % n1 * cell % n2 * cell % n3
    !
    ! Set basic cell properties
    !
    CALL update_environ_cell( at, cell )
    !
    ! Cell origin
    !
    cell % origin = 0.D0
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE init_environ_cell
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE update_environ_cell( at, cell )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    REAL( DP ), INTENT(IN) :: at(3,3)
    TYPE( environ_cell ), INTENT(INOUT) :: cell
    CHARACTER( LEN=80 ) :: sub_name = 'update_environ_cell'
    !
    INTEGER :: ic, ix, iy, iz
    REAL( DP ) :: dx, dy, dz
    !
    cell % at = at
    !
    ! Calculate cell volume
    !
    CALL volume( cell%alat, cell%at(1,1), cell%at(1,2), cell%at(1,3), cell%omega )
    !
    ! Check if the cell is cubic
    !
    cell % cubic = iscubic( cell%at )
    !
    ! Calculate reciprocal cell
    !
    CALL recips( cell%at(1,1), cell%at(1,2), cell%at(1,3), &
         & cell%bg(1,1), cell%bg(1,2), cell%bg(1,3) )
    !
    ! Set volume element
    !
    cell % domega = cell % omega / cell % ntot
    !
    ! Calcualte corners for minimum image convetion
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
    ! The following is an alternative way of implementing
    ! minimum image distance, we may want to check if
    ! it is safer/more efficient
    !
!    x = MATMUL(ws%b,r)
!    x(:) = x(:) - NINT(x(:))
!    c  = SUM(x*MATMUL(ws%aa,x))
!    m = 0
!    !
!    lb(:) =  NINT ( x(:) - DSQRT (c) * ws%norm_b(:) )
!    ! CEILING should be enough for lb but NINT might be safer
!    ub(:) =  NINT ( x(:) + DSQRT (c) * ws%norm_b(:) )
!    ! FLOOR should be enough for ub but NINT might be safer
!    !
!    DO i1 = lb(1), ub(1)
!       DO i2 = lb(2), ub(2)
!          DO i3 = lb(3), ub(3)
!             y = x - (/i1,i2,i3/)
!             ctest = SUM(y*MATMUL(ws%aa,y))
!             IF (ctest < c) THEN
!                c = ctest
!                m = (/i1,i2,i3/)
!             END IF
!          END DO
!       END DO
!    END DO
!    !
!    y = x-m
!    r_ws =  MATMUL(ws%a,y)
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE minimum_image
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE volume( alat, a1, a2, a3, omega )
!--------------------------------------------------------------------
    !
    !     Compute the volume of the unit cell defined by 3 vectors
    !     a1, a2, a3, given in units of "alat" (alat may be 1):
    !        omega = alat^3 * [ a1 . (a2 x a3) ]
    !     ( . = scalar product, x = vector product )
    !
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
!---------------------------------------------------------------------
  END SUBROUTINE
!---------------------------------------------------------------------
!---------------------------------------------------------------------
  FUNCTION iscubic( at )
!---------------------------------------------------------------------
    !
    IMPLICIT NONE
    REAL(DP), PARAMETER :: tol = 1.D-8
    REAL(DP) :: at(3,3)
    LOGICAL :: iscubic
    !
    INTEGER :: ipol, jpol
    REAL(DP) :: tmp
    !
    iscubic = .FALSE.
    !
    ! If at(3,3) is a cubic cell, at(1,1)=at(2,2)=at(3,3)
    ! and the other elements are equal to 0.D0
    !
    tmp = 0.D0
    DO ipol = 1, 3
       DO jpol = 1, 3
          IF ( ipol .EQ. jpol ) THEN
             tmp = tmp + ABS( at(ipol,ipol) - at(1,1) )
          ELSE
             tmp = tmp + ABS( at(ipol,jpol) )
          ENDIF
       ENDDO
    ENDDO
    !
    iscubic = tmp .LT. tol
    !
    RETURN
    !
!---------------------------------------------------------------------
  END FUNCTION iscubic
!---------------------------------------------------------------------
!---------------------------------------------------------------------
  SUBROUTINE recips( a1, a2, a3, b1, b2, b3 )
!---------------------------------------------------------------------
    !
    !   This routine generates the reciprocal lattice vectors b1,b2,b3
    !   given the real space vectors a1,a2,a3. The b's are units of 2 pi/a.
    !
    !     first the input variables
    !
    IMPLICIT NONE
    REAL( DP ) :: a1 (3), a2 (3), a3 (3), b1 (3), b2 (3), b3 (3)
    ! input: first direct lattice vector
    ! input: second direct lattice vector
    ! input: third direct lattice vector
    ! output: first reciprocal lattice vector
    ! output: second reciprocal lattice vector
    ! output: third reciprocal lattice vector
    !
    !   then the local variables
    !
    REAL( DP ) :: den, s
    ! the denominator
    ! the sign of the permutations
    INTEGER :: iperm, i, j, k, l
    ! counter on the permutations
    !\
    !  Auxiliary variables
    !/
    !
    INTEGER :: ipol
    ! Counter on the polarizations
    !
    !    first we compute the denominator
    !
    den = 0
    i = 1
    j = 2
    k = 3
    s = 1.d0
100 DO iperm = 1, 3
       den = den + s * a1 (i) * a2 (j) * a3 (k)
       l = i
       i = j
       j = k
       k = l
    ENDDO
    i = 2
    j = 1
    k = 3
    s = - s
    IF (s.LT.0.D0) GOTO 100
    !
    !    here we compute the reciprocal vectors
    !
    i = 1
    j = 2
    k = 3
    DO ipol = 1, 3
       b1 (ipol) = (a2 (j) * a3 (k) - a2 (k) * a3 (j) ) / den
       b2 (ipol) = (a3 (j) * a1 (k) - a3 (k) * a1 (j) ) / den
       b3 (ipol) = (a1 (j) * a2 (k) - a1 (k) * a2 (j) ) / den
       l = i
       i = j
       j = k
       k = l
    ENDDO
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE recips
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE create_environ_mapping( mapping )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    TYPE( environ_mapping ), INTENT(INOUT) :: mapping
    !
    mapping % nrep = 1
    mapping % large => null()
    mapping % small => null()
    !
    RETURN
!--------------------------------------------------------------------
  END SUBROUTINE create_environ_mapping
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE init_environ_mapping_first( nrep, mapping )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    INTEGER, DIMENSION(3), INTENT(IN) :: nrep
    TYPE( environ_mapping ), INTENT(INOUT) :: mapping
    !
    mapping % nrep = nrep
    !
    RETURN
!--------------------------------------------------------------------
  END SUBROUTINE init_environ_mapping_first
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE init_environ_mapping_second( small, large, mapping )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    TYPE( environ_cell ), TARGET, INTENT(IN) :: small, large
    TYPE( environ_mapping ), INTENT(INOUT) :: mapping
    !
    ! Check that small%at and large%at are compatible with mapping%nrep
    !
    mapping % small => small
    mapping % large => large
    !
    ALLOCATE( mapping % map( mapping % small % nnr ) )
    !
    CALL update_environ_mapping( mapping )
    !
    RETURN
!--------------------------------------------------------------------
  END SUBROUTINE init_environ_mapping_second
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE update_environ_mapping( mapping )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    TYPE( environ_mapping ), INTENT(INOUT) :: mapping
    !
    LOGICAL :: physical
    INTEGER :: ir, ipol
    INTEGER, DIMENSION(3) :: small_n, large_n, center, origin, shift, ijk
    !
    ! ... Compute mapping
    !
    small_n(1) = mapping % small % n1
    small_n(2) = mapping % small % n2
    small_n(3) = mapping % small % n3
    !
    large_n(1) = mapping % large % n1
    large_n(2) = mapping % large % n2
    large_n(3) = mapping % large % n3
    !
    ! ... Indexes of center of small cell
    !
    center = NINT( small_n / 2 )
    !
    ! ... Indexes of origin of small cell
    !
    origin = MATMUL( mapping%small%bg, mapping%small%origin )
    origin(:) = NINT( origin(:) * small_n(:) )
    !
    shift = center - origin
    !
    mapping%map = 0
    !
    DO ir = 1, mapping%small%ir_end
       !
       CALL ir2ijk( mapping%small, ir, ijk(1), ijk(2), ijk(3), physical )
       !
       IF ( .NOT. physical ) CYCLE
       !
       ! ... Shift to center small cell
       !
       ijk = ijk + shift
       ijk = ijk - FLOOR(ijk/n)*n
       !
       ! ... Map small cell to large cell
       !
       ijk = ijk + small_n * mapping%nrep
       !
       mapping%map(ir) = 1 + ijk(1) + ijk(2) * large_n(1) + &
            & ijk(3) * large_n(1) * large_n(2)
       !
    END DO
    !
    RETURN
!--------------------------------------------------------------------
  END SUBROUTINE update_environ_mapping
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE destroy_environ_mapping( lflag, mapping )
!--------------------------------------------------------------------
    IMPLICIT NONE
    !
    LOGICAL, INTENT(IN) :: lflag
    TYPE( environ_mapping ), INTENT(INOUT) :: mapping
    !
    mapping % nrep = 1
    NULLIFY( mapping%small )
    NULLIFY( mapping%large )
    !
    DEALLOCATE( mapping%map )
    !
    RETURN
!--------------------------------------------------------------------
  END SUBROUTINE destroy_environ_mapping
!--------------------------------------------------------------------
END MODULE cell_types
