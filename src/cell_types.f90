!----------------------------------------------------------------------------------------
!>
!!
!----------------------------------------------------------------------------------------
MODULE cell_types
    !------------------------------------------------------------------------------------
    !
    USE modules_constants, ONLY: DP, tpi
    USE stick_base, ONLY: sticks_map, sticks_map_deallocate
    USE fft_types, ONLY: fft_type_descriptor, fft_type_init, fft_type_deallocate
    !
    !------------------------------------------------------------------------------------
    !>
    !! A simulation cell
    !!
    !! Notes:
    !!
    !! 1. ir_end can be different from nnr due to FFT-grid optimization yielding
    !!    additional, unphysical grid points
    !!
    !! 2. cell corners are utilized in minimum_image()
    !!
    !------------------------------------------------------------------------------------
    TYPE environ_cell
        !--------------------------------------------------------------------------------
        !
        LOGICAL :: update = .FALSE.
        LOGICAL :: cubic = .FALSE.
        REAL(DP) :: omega, domega ! volume quantities
        REAL(DP) :: origin(3)
        REAL(DP) :: at(3, 3) ! real-space lattice vectors
        REAL(DP) :: bg(3, 3) ! reciprocal lattice vectors
        REAL(DP) :: corners(3, 8)
        !
        !--------------------------------------------------------------------------------
        ! Units needed to scale real and reciprocal space cells
        ! #TODO redesign Environ units
        !
        REAL(DP) :: alat ! lattice parameter
        REAL(DP) :: tpiba, tpiba2 ! tpi = 2*pi
        !
        !--------------------------------------------------------------------------------
        ! Properties of the grid
        !
        TYPE(fft_type_descriptor) :: dfft
        INTEGER :: ntot ! total number of grid points
        INTEGER :: nnr ! number of grid points allocated in every processor
        INTEGER :: ir_end ! actual number grid points accessed by each processor
        INTEGER :: j0, k0 ! starting indexes of processor-specific boxes of grid points
        REAL(DP) :: in1, in2, in3 ! inverse number of grid points
        !
        !--------------------------------------------------------------------------------
    END TYPE environ_cell
    !------------------------------------------------------------------------------------
    !>
    !! The mapping between system and environment simulation cells
    !!
    !------------------------------------------------------------------------------------
    TYPE environ_mapping
        !--------------------------------------------------------------------------------
        !
        INTEGER :: nrep(3)
        ! number of system cells in environment cell along a_i = 2 * nrep_i + 1
        !
        TYPE(environ_cell), POINTER :: small ! system cell
        TYPE(environ_cell), POINTER :: large ! environment cell
        !
        INTEGER, ALLOCATABLE :: map(:)
        !
        !--------------------------------------------------------------------------------
    END TYPE environ_mapping
    !------------------------------------------------------------------------------------
    !
    PRIVATE
    !
    PUBLIC :: fft_type_descriptor, environ_cell, create_environ_cell, &
              init_environ_cell, update_environ_cell, destroy_environ_cell, ir2ijk, &
              ir2r, displacement, minimum_image, volume, recips, environ_mapping, &
              init_environ_mapping_first, init_environ_mapping_second, &
              update_environ_mapping, destroy_environ_mapping
    !
    !------------------------------------------------------------------------------------
CONTAINS
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_dfft(gcutm, comm, at, dfft)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: comm
        REAL(DP), INTENT(IN) :: gcutm, at(3, 3)
        !
        TYPE(fft_type_descriptor), INTENT(INOUT) :: dfft
        !
        TYPE(sticks_map) :: smap
        REAL(DP) :: bg(3, 3)
        !
        !--------------------------------------------------------------------------------
        ! Calculate the reciprocal lattice vectors
        !
        CALL recips(at(1, 1), at(1, 2), at(1, 3), bg(1, 1), bg(1, 2), bg(1, 3))
        !
        CALL fft_type_init(dfft, smap, "rho", .TRUE., .TRUE., comm, at, bg, gcutm, &
                           nyfft=1, nmany=1)
        !
        dfft%rho_clock_label = 'fft'
        !
        CALL sticks_map_deallocate(smap)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_dfft
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE destroy_dfft(lflag, dfft)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        LOGICAL, INTENT(IN) :: lflag
        !
        TYPE(fft_type_descriptor), INTENT(INOUT) :: dfft
        !
        !--------------------------------------------------------------------------------
        !
        CALL fft_type_deallocate(dfft)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE destroy_dfft
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_environ_cell(gcutm, comm, alat, at, cell, nr)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: comm
        REAL(DP), INTENT(IN) :: gcutm, alat, at(3, 3)
        INTEGER, OPTIONAL, INTENT(IN) :: nr(3)
        !
        TYPE(environ_cell), INTENT(INOUT) :: cell
        !
        CHARACTER(LEN=80) :: sub_name = 'init_environ_cell'
        !
        !--------------------------------------------------------------------------------
        !
        IF (alat < 1.D-8) CALL errore(sub_name, 'Wrong alat', 1)
        !
        cell%alat = alat
        !
        !--------------------------------------------------------------------------------
        ! Calculate units of reciprocal space
        !
        cell%tpiba = tpi / alat
        cell%tpiba2 = cell%tpiba**2
        !
        !--------------------------------------------------------------------------------
        ! Create fft descriptor for system cell
        !
        IF (PRESENT(nr)) THEN
            cell%dfft%nr1 = nr(1)
            cell%dfft%nr2 = nr(2)
            cell%dfft%nr3 = nr(3)
        END IF
        !
        CALL init_dfft(gcutm, comm, at, cell%dfft)
        !
        cell%in1 = 1.D0 / DBLE(cell%dfft%nr1)
        cell%in2 = 1.D0 / DBLE(cell%dfft%nr2)
        cell%in3 = 1.D0 / DBLE(cell%dfft%nr3)
        !
        !--------------------------------------------------------------------------------
        ! Real space grid, local dimensions (processor-specific)
        !
        cell%nnr = cell%dfft%nnr
#if defined (__MPI)
        cell%j0 = cell%dfft%my_i0r2p
        cell%k0 = cell%dfft%my_i0r3p
        !
        cell%ir_end = MIN(cell%nnr, &
                          cell%dfft%nr1x * cell%dfft%my_nr2p * cell%dfft%my_nr3p)
#else
        cell%j0 = 0
        cell%k0 = 0
        cell%ir_end = cell%nnr
#endif
        !
        cell%ntot = cell%dfft%nr1 * cell%dfft%nr2 * cell%dfft%nr3
        ! Total number of physical points
        !
        !--------------------------------------------------------------------------------
        ! Set basic cell properties
        !
        CALL update_environ_cell(at, cell)
        !
        cell%origin = 0.D0
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_environ_cell
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE update_environ_cell(at, cell)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP), INTENT(IN) :: at(3, 3)
        !
        TYPE(environ_cell), INTENT(INOUT) :: cell
        !
        CHARACTER(LEN=80) :: sub_name = 'update_environ_cell'
        INTEGER :: ic, ix, iy, iz
        REAL(DP) :: dx, dy, dz
        !
        !--------------------------------------------------------------------------------
        !
        cell%at = at
        !
        CALL volume(cell%alat, cell%at(1, 1), cell%at(1, 2), cell%at(1, 3), cell%omega)
        ! Calculate cell volume
        !
        cell%cubic = iscubic(cell%at) ! Check if the cell is cubic
        !
        ! Calculate reciprocal cell
        CALL recips(cell%at(1, 1), cell%at(1, 2), cell%at(1, 3), &
                    cell%bg(1, 1), cell%bg(1, 2), cell%bg(1, 3))
        !
        cell%domega = cell%omega / cell%ntot ! Set volume element
        !
        !--------------------------------------------------------------------------------
        ! Calculate corners for minimum image convention
        !
        ic = 0
        !
        DO ix = 0, 1
            dx = DBLE(-ix)
            !
            DO iy = 0, 1
                dy = DBLE(-iy)
                !
                DO iz = 0, 1
                    dz = DBLE(-iz)
                    ic = ic + 1
                    cell%corners(1, ic) = dx * at(1, 1) + dy * at(1, 2) + dz * at(1, 3)
                    cell%corners(2, ic) = dx * at(2, 1) + dy * at(2, 2) + dz * at(2, 3)
                    cell%corners(3, ic) = dx * at(3, 1) + dy * at(3, 2) + dz * at(3, 3)
                END DO
                !
            END DO
            !
        END DO
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE update_environ_cell
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE destroy_environ_cell(lflag, cell)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        LOGICAL, INTENT(IN) :: lflag
        !
        TYPE(environ_cell), INTENT(INOUT) :: cell
        !
        !--------------------------------------------------------------------------------
        !
        CALL destroy_dfft(lflag, cell%dfft)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE destroy_environ_cell
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE ir2ijk(cell, ir, i, j, k, physical)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_cell), INTENT(IN) :: cell
        INTEGER, INTENT(IN) :: ir
        !
        INTEGER, INTENT(OUT) :: i, j, k
        LOGICAL, INTENT(OUT) :: physical
        !
        INTEGER :: idx
        !
        !--------------------------------------------------------------------------------
        ! Convert single ir index to i, j, k
        !
        idx = ir - 1
        k = idx / (cell%dfft%nr1x * cell%dfft%my_nr2p)
        idx = idx - (cell%dfft%nr1x * cell%dfft%my_nr2p) * k
        k = k + cell%k0
        j = idx / cell%dfft%nr1x
        idx = idx - cell%dfft%nr1x * j
        j = j + cell%j0
        i = idx
        !
        physical = i < cell%dfft%nr1 .AND. j < cell%dfft%nr2 .AND. k < cell%dfft%nr3
        ! check if current point was generated for optimization of fft grids
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE ir2ijk
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE ir2r(cell, ir, r, physical)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_cell), INTENT(IN) :: cell
        INTEGER, INTENT(IN) :: ir
        !
        REAL(DP), INTENT(OUT) :: r(3)
        LOGICAL, INTENT(OUT) :: physical
        !
        INTEGER :: idx, i, j, k, ip
        !
        !--------------------------------------------------------------------------------
        !
        r = 0.D0
        !
        CALL ir2ijk(cell, ir, i, j, k, physical)
        !
        IF (.NOT. physical) RETURN
        !
        DO ip = 1, 3
            !
            r(ip) = DBLE(i) * cell%in1 * cell%at(ip, 1) + &
                    DBLE(j) * cell%in2 * cell%at(ip, 2) + &
                    DBLE(k) * cell%in3 * cell%at(ip, 3)
            !
        END DO
        !
        r = r + cell%origin
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE ir2r
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE displacement(dim, axis, r1, r2, dr)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: dim, axis
        REAL(DP), DIMENSION(3), INTENT(IN) :: r1, r2
        !
        REAL(DP), INTENT(OUT) :: dr(3)
        !
        INTEGER :: i
        !
        !--------------------------------------------------------------------------------
        ! Possibly only in 1D or 2D
        !
        dr(:) = r1(:) - r2(:)
        !
        SELECT CASE (dim)
        CASE (1)
            dr(axis) = 0.D0
        CASE (2)
            !
            DO i = 1, 3
                IF (i /= axis) dr(i) = 0.D0
            END DO
            !
        END SELECT
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE displacement
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE minimum_image(cell, r, r2)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_cell), INTENT(IN) :: cell
        !
        REAL(DP), INTENT(INOUT) :: r(3)
        REAL(DP), INTENT(OUT) :: r2
        !
        INTEGER :: ic
        REAL(DP) :: s(3)
        REAL(DP) :: rmin(3)
        REAL(DP) :: r2min
        !
        s(:) = MATMUL(r(:), cell%bg(:, :))
        s(:) = s(:) - FLOOR(s(:)) ! #TODO NINT?
        r(:) = MATMUL(cell%at(:, :), s(:))
        !
        rmin = r
        r2min = SUM(r * r)
        !
        DO ic = 2, 8
            s(1) = r(1) + cell%corners(1, ic)
            s(2) = r(2) + cell%corners(2, ic)
            s(3) = r(3) + cell%corners(3, ic)
            r2 = SUM(s * s)
            !
            IF (r2 < r2min) THEN
                rmin = s
                r2min = r2
            END IF
            !
        END DO
        !
        r = rmin
        r2 = r2min
        !
        !--------------------------------------------------------------------------------
        ! The following is an alternative way of implementing minimum image distance,
        ! #TODO we may want to check if it is safer/more efficient
        !
        ! x = MATMUL(ws%b, r)
        ! x(:) = x(:) - NINT(x(:))
        ! c = SUM(x * MATMUL(ws%aa, x))
        ! m = 0
        ! !
        ! lb(:) = NINT(x(:) - DSQRT(c) * ws%norm_b(:))
        ! ! CEILING should be enough for lb but NINT might be safer
        ! !
        ! ub(:) = NINT(x(:) + DSQRT(c) * ws%norm_b(:))
        ! ! FLOOR should be enough for ub but NINT might be safer
        ! !
        ! DO i1 = lb(1), ub(1)
        !     !
        !     DO i2 = lb(2), ub(2)
        !         !
        !         DO i3 = lb(3), ub(3)
        !             y = x - (/i1, i2, i3/)
        !             ctest = SUM(y * MATMUL(ws%aa, y))
        !             !
        !             IF (ctest < c) THEN
        !                 c = ctest
        !                 m = (/i1, i2, i3/)
        !             END IF
        !             !
        !         END DO
        !         !
        !     END DO
        !     !
        ! END DO
        ! !
        ! y = x - m
        ! r_ws = MATMUL(ws%a, y)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE minimum_image
    !------------------------------------------------------------------------------------
    !>
    !! Compute the volume of the unit cell defined by 3 vectors
    !! a1, a2, a3, given in units of "alat" (alat may be 1):
    !!
    !! omega = alat^3 * [ a1 . (a2 x a3) ]
    !!
    !! ( . = scalar product, x = vector product )
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE volume(alat, a1, a2, a3, omega)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(dp), INTENT(IN) :: alat, a1(3), a2(3), a3(3)
        !
        REAL(dp), INTENT(OUT) :: omega
        !
        !--------------------------------------------------------------------------------
        !
        omega = a1(1) * (a2(2) * a3(3) - a2(3) * a3(2)) - &
                a1(2) * (a2(1) * a3(3) - a2(3) * a3(1)) + &
                a1(3) * (a2(1) * a3(2) - a2(2) * a3(1))
        !
        IF (omega < 0.0_DP) THEN
            !
            CALL infomsg('volume', 'axis vectors are left-handed')
            !
            omega = ABS(omega)
        END IF
        !
        IF (alat < 1.0_DP) CALL infomsg('volume', 'strange lattice parameter')
        !
        omega = omega * alat**3
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    FUNCTION iscubic(at)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP), PARAMETER :: tol = 1.D-8
        REAL(DP) :: at(3, 3)
        LOGICAL :: iscubic
        !
        INTEGER :: ipol, jpol
        REAL(DP) :: tmp
        !
        !--------------------------------------------------------------------------------
        ! If at(3, 3) is a cubic cell, at(1, 1) = at(2, 2) = at(3, 3) and
        ! the other elements are equal to 0.D0
        !
        iscubic = .FALSE.
        !
        tmp = 0.D0
        !
        DO ipol = 1, 3
            !
            DO jpol = 1, 3
                !
                IF (ipol == jpol) THEN
                    tmp = tmp + ABS(at(ipol, ipol) - at(1, 1))
                ELSE
                    tmp = tmp + ABS(at(ipol, jpol))
                END IF
                !
            END DO
            !
        END DO
        !
        iscubic = tmp < tol
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END FUNCTION iscubic
    !------------------------------------------------------------------------------------
    !>
    !! This routine generates the reciprocal lattice vectors b1, b2, b3
    !! given the real space vectors a1, a2, a3. The b vectors are in units of 2*pi/a.
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE recips(a1, a2, a3, b1, b2, b3)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP), DIMENSION(3), INTENT(IN) :: a1, a2, a3
        !
        REAL(DP), DIMENSION(3), INTENT(OUT) :: b1, b2, b3
        !
        REAL(DP) :: den ! the denominator
        REAL(DP) :: s ! the sign of the permutations
        !
        INTEGER :: iperm ! counter on the permutations
        INTEGER :: i, j, k, l
        !
        INTEGER :: ipol ! counter on the polarizations
        !
        !--------------------------------------------------------------------------------
        ! Compute the denominator
        !
        den = 0
        i = 1
        j = 2
        k = 3
        s = 1.D0
        !
100     DO iperm = 1, 3
            den = den + s * a1(i) * a2(j) * a3(k)
            l = i
            i = j
            j = k
            k = l
        END DO
        !
        i = 2
        j = 1
        k = 3
        s = -s
        !
        IF (s < 0.D0) GOTO 100
        !
        !--------------------------------------------------------------------------------
        ! Compute the reciprocal vectors
        !
        i = 1
        j = 2
        k = 3
        !
        DO ipol = 1, 3
            b1(ipol) = (a2(j) * a3(k) - a2(k) * a3(j)) / den
            b2(ipol) = (a3(j) * a1(k) - a3(k) * a1(j)) / den
            b3(ipol) = (a1(j) * a2(k) - a1(k) * a2(j)) / den
            l = i
            i = j
            j = k
            k = l
        END DO
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE recips
    !------------------------------------------------------------------------------------
    !>
    !! never used (Leave it)
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE create_environ_mapping(mapping)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_mapping), INTENT(INOUT) :: mapping
        !
        !--------------------------------------------------------------------------------
        !
        mapping%nrep = 1
        mapping%large => NULL()
        mapping%small => NULL()
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE create_environ_mapping
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_environ_mapping_first(nrep, mapping)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: nrep(3)
        !
        TYPE(environ_mapping), INTENT(INOUT) :: mapping
        !
        !--------------------------------------------------------------------------------
        !
        mapping%nrep = nrep
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_environ_mapping_first
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_environ_mapping_second(small, large, mapping)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_cell), TARGET, INTENT(IN) :: small, large
        !
        TYPE(environ_mapping), INTENT(INOUT) :: mapping
        !
        !--------------------------------------------------------------------------------
        ! Check that small%at and large%at are compatible with mapping%nrep
        !
        mapping%small => small
        mapping%large => large
        !
        IF (.NOT. ASSOCIATED(mapping%small, mapping%large)) &
            ALLOCATE (mapping%map(mapping%small%nnr))
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_environ_mapping_second
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE update_environ_mapping(mapping, pos)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP), INTENT(IN), OPTIONAL :: pos(3)
        !
        TYPE(environ_mapping), INTENT(INOUT) :: mapping
        !
        LOGICAL :: physical
        INTEGER :: ir, ipol
        INTEGER, DIMENSION(3) :: small_n, large_n, center, origin, shift, ijk
        REAL(DP) :: tmp(3)
        !
        !--------------------------------------------------------------------------------
        !
        IF (ASSOCIATED(mapping%small, mapping%large)) RETURN
        ! environment cell is identical to system cell (env_nrep = 0 in environ.in)
        !
        !--------------------------------------------------------------------------------
        ! Compute mapping
        !
        small_n(1) = mapping%small%dfft%nr1
        small_n(2) = mapping%small%dfft%nr2
        small_n(3) = mapping%small%dfft%nr3
        !
        large_n(1) = mapping%large%dfft%nr1
        large_n(2) = mapping%large%dfft%nr2
        large_n(3) = mapping%large%dfft%nr3
        !
        !--------------------------------------------------------------------------------
        !
        center = NINT(small_n / 2.D0) ! Indexes of center of small cell
        !
        !--------------------------------------------------------------------------------
        ! Indexes of origin of small cell
        !
        IF (PRESENT(pos)) THEN
            tmp = MATMUL(mapping%small%bg, pos) ! #TODO center of charge (think molecule)
            origin = NINT(tmp * small_n)
        ELSE
            origin = 0 ! center of charge
        END IF
        !
        shift = center - origin
        !
        !--------------------------------------------------------------------------------
        ! Shift origin of large cell
        !
        mapping%large%origin = -MATMUL(mapping%large%at, (0.5 - origin / DBLE(large_n)))
        !
        mapping%map = 0
        !
        DO ir = 1, mapping%small%ir_end
            !
            CALL ir2ijk(mapping%small, ir, ijk(1), ijk(2), ijk(3), physical)
            !
            IF (.NOT. physical) CYCLE
            !
            !----------------------------------------------------------------------------
            ! Shift to center small cell
            !
            ijk = ijk + shift
            ijk = ijk - FLOOR(DBLE(ijk) / small_n) * small_n ! enforce periodicity
            !
            !----------------------------------------------------------------------------
            ! Map small cell to large cell #TODO check if this works in parallel
            !
            ijk = ijk + small_n * mapping%nrep
            !
            mapping%map(ir) = 1 + ijk(1) &                         ! x-point
                              & + ijk(2) * large_n(1) &            ! y-row
                              & + ijk(3) * large_n(1) * large_n(2) ! z-plane
            !
        END DO
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE update_environ_mapping
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE destroy_environ_mapping(lflag, mapping)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        LOGICAL, INTENT(IN) :: lflag
        !
        TYPE(environ_mapping), INTENT(INOUT) :: mapping
        !
        !--------------------------------------------------------------------------------
        !
        mapping%nrep = 1
        NULLIFY (mapping%small)
        NULLIFY (mapping%large)
        !
        IF (ALLOCATED(mapping%map)) DEALLOCATE (mapping%map)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE destroy_environ_mapping
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE cell_types
!----------------------------------------------------------------------------------------
