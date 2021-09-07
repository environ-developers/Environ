!----------------------------------------------------------------------------------------
!
! Copyright (C) 2021 ENVIRON (www.quantum-environ.org)
!
!----------------------------------------------------------------------------------------
!
!     This file is part of Environ version 2.0
!
!     Environ 2.0 is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation, either version 2 of the License, or
!     (at your option) any later version.
!
!     Environ 2.0 is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more detail, either the file
!     `License' in the root directory of the present distribution, or
!     online at <http://www.gnu.org/licenses/>.
!
!----------------------------------------------------------------------------------------
!
! Authors: Oliviero Andreussi (Department of Physics, UNT)
!          Gabriel Medrano    (Department of Physics, UNT)
!          Edan Bainglass     (Department of Physics, UNT)
!
!----------------------------------------------------------------------------------------
!>
!!
!----------------------------------------------------------------------------------------
MODULE class_cell
    !------------------------------------------------------------------------------------
    !
    USE env_base_io, ONLY: ionode, environ_unit, global_verbose
    USE env_mp, ONLY: env_mp_sum
    !
    USE environ_param, ONLY: DP, tpi
    !
    USE env_base_stick, ONLY: env_sticks_map, env_sticks_map_deallocate
    !
    USE env_types_fft, ONLY: env_fft_type_descriptor, env_fft_type_init, &
                             env_fft_type_deallocate
    !
    !------------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    PRIVATE
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
    TYPE, PUBLIC :: environ_cell
        !--------------------------------------------------------------------------------
        !
        LOGICAL :: lupdate
        LOGICAL :: cubic
        REAL(DP) :: omega, domega ! volume quantities
        REAL(DP) :: origin(3)
        REAL(DP) :: at(3, 3) ! real-space lattice vectors
        REAL(DP) :: bg(3, 3) ! reciprocal lattice vectors
        REAL(DP) :: corners(3, 8)
        !
        !--------------------------------------------------------------------------------
        ! Properties of the grid
        !
        TYPE(env_fft_type_descriptor) :: dfft
        INTEGER :: ntot ! total number of grid points
        INTEGER :: nnr ! number of grid points allocated in every processor
        INTEGER :: ir_end ! actual number grid points accessed by each processor
        INTEGER :: j0, k0 ! starting indexes of processor-specific boxes of grid points
        REAL(DP) :: in1, in2, in3 ! inverse number of grid points
        !
        !--------------------------------------------------------------------------------
    CONTAINS
        !--------------------------------------------------------------------------------
        !
        PROCEDURE, PRIVATE :: create => create_environ_cell
        PROCEDURE :: init => init_environ_cell
        PROCEDURE :: update => update_environ_cell
        PROCEDURE :: destroy => destroy_environ_cell
        PROCEDURE, PRIVATE :: init_dfft, destroy_dfft
        !
        PROCEDURE :: volume, get_min_distance, ir2ijk, planar_average
        PROCEDURE, PRIVATE :: ir2r, minimum_image, is_cubic
        !
        PROCEDURE :: printout => print_environ_cell
        PROCEDURE :: write_cube => write_cube_cell
        !
        !--------------------------------------------------------------------------------
    END TYPE environ_cell
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
CONTAINS
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !                                   ADMIN METHODS
    !
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE create_environ_cell(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_cell), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: sub_name = 'create_environ_cell'
        !
        !--------------------------------------------------------------------------------
        !
        this%lupdate = .FALSE.
        this%cubic = .FALSE.
        !
        this%origin = 0.D0
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE create_environ_cell
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_environ_cell(this, gcutm, comm, at, nr)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: comm
        REAL(DP), INTENT(IN) :: gcutm, at(3, 3)
        INTEGER, OPTIONAL, INTENT(IN) :: nr(3)
        !
        CLASS(environ_cell), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: sub_name = 'init_environ_cell'
        !
        !--------------------------------------------------------------------------------
        ! Create fft descriptor for system cell
        !
        CALL this%create()
        !
        IF (PRESENT(nr)) THEN
            this%dfft%nr1 = nr(1)
            this%dfft%nr2 = nr(2)
            this%dfft%nr3 = nr(3)
        END IF
        !
        CALL this%init_dfft(gcutm, comm, at)
        !
        this%in1 = 1.D0 / DBLE(this%dfft%nr1)
        this%in2 = 1.D0 / DBLE(this%dfft%nr2)
        this%in3 = 1.D0 / DBLE(this%dfft%nr3)
        !
        !--------------------------------------------------------------------------------
        ! Real space grid, local dimensions (processor-specific)
        !
        this%nnr = this%dfft%nnr
#if defined(__MPI)
        this%j0 = this%dfft%my_i0r2p
        this%k0 = this%dfft%my_i0r3p
        !
        this%ir_end = MIN(this%nnr, &
                          this%dfft%nr1x * this%dfft%my_nr2p * this%dfft%my_nr3p)
#else
        this%j0 = 0
        this%k0 = 0
        this%ir_end = this%nnr
#endif
        !
        this%ntot = this%dfft%nr1 * this%dfft%nr2 * this%dfft%nr3
        ! total number of physical points
        !
        !--------------------------------------------------------------------------------
        ! Set basic cell properties
        !
        CALL this%update(at)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_environ_cell
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE update_environ_cell(this, at)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP), INTENT(IN) :: at(3, 3)
        !
        CLASS(environ_cell), INTENT(INOUT) :: this
        !
        INTEGER :: ic, ix, iy, iz
        REAL(DP) :: dx, dy, dz
        !
        CHARACTER(LEN=80) :: sub_name = 'update_environ_cell'
        !
        !--------------------------------------------------------------------------------
        !
        this%at = at
        !
        CALL this%volume(at(1, 1), at(1, 2), at(1, 3)) ! calculate cell volume
        !
        this%cubic = this%is_cubic() ! check if the cell is cubic (stored once)
        !
        !--------------------------------------------------------------------------------
        ! Calculate reciprocal cell
        !
        CALL recips(this%at(1, 1), this%at(1, 2), this%at(1, 3), &
                    this%bg(1, 1), this%bg(1, 2), this%bg(1, 3))
        !
        this%domega = this%omega / this%ntot ! Set volume element
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
                    this%corners(1, ic) = dx * at(1, 1) + dy * at(1, 2) + dz * at(1, 3)
                    this%corners(2, ic) = dx * at(2, 1) + dy * at(2, 2) + dz * at(2, 3)
                    this%corners(3, ic) = dx * at(3, 1) + dy * at(3, 2) + dz * at(3, 3)
                END DO
                !
            END DO
            !
        END DO
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE update_environ_cell
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE destroy_environ_cell(this, lflag)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        LOGICAL, INTENT(IN) :: lflag
        !
        CLASS(environ_cell), INTENT(INOUT) :: this
        !
        !--------------------------------------------------------------------------------
        !
        CALL this%destroy_dfft(lflag)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE destroy_environ_cell
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !                               PRIVATE ADMIN METHODS
    !
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !>
    !! Local fft-descriptor constructor
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_dfft(this, gcutm, comm, at)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: comm
        REAL(DP), INTENT(IN) :: gcutm, at(3, 3)
        !
        CLASS(environ_cell), INTENT(INOUT) :: this
        !
        TYPE(env_sticks_map) :: smap
        !
        REAL(DP) :: bg(3, 3)
        !
        !--------------------------------------------------------------------------------
        !
        CALL recips(at(1, 1), at(1, 2), at(1, 3), bg(1, 1), bg(1, 2), bg(1, 3))
        ! calculate the reciprocal lattice vectors
        !
        CALL env_fft_type_init(this%dfft, smap, .TRUE., .TRUE., comm, at, bg, gcutm, &
                               nyfft=1, nmany=1)
        !
        this%dfft%rho_clock_label = 'fft'
        !
        CALL env_sticks_map_deallocate(smap)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_dfft
    !------------------------------------------------------------------------------------
    !>
    !! Local fft-descriptor destructor
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE destroy_dfft(this, lflag)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        LOGICAL, INTENT(IN) :: lflag
        !
        CLASS(environ_cell), INTENT(INOUT) :: this
        !
        !--------------------------------------------------------------------------------
        !
        CALL env_fft_type_deallocate(this%dfft)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE destroy_dfft
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !                                  GENERAL METHODS
    !
    !------------------------------------------------------------------------------------
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
    SUBROUTINE volume(this, a1, a2, a3)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP), DIMENSION(3), INTENT(IN) :: a1, a2, a3
        !
        CLASS(environ_cell), TARGET, INTENT(INOUT) :: this
        !
        !--------------------------------------------------------------------------------
        !
        this%omega = a1(1) * (a2(2) * a3(3) - a2(3) * a3(2)) - &
                     a1(2) * (a2(1) * a3(3) - a2(3) * a3(1)) + &
                     a1(3) * (a2(1) * a3(2) - a2(2) * a3(1))
        !
        IF (this%omega < 0.0_DP) THEN
            this%omega = ABS(this%omega)
            !
            CALL env_warning('axis vectors are left-handed')
            !
        END IF
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE
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
        !--------------------------------------------------------------------------------
    END SUBROUTINE recips
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE get_min_distance(this, ir, dim, axis, pos, r, r2, physical)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_cell), INTENT(IN) :: this
        INTEGER, INTENT(IN) :: ir, dim, axis
        REAL(DP), INTENT(IN) :: pos(3)
        !
        REAL(DP), INTENT(OUT) :: r(3), r2
        !
        LOGICAL, INTENT(INOUT) :: physical
        !
        !--------------------------------------------------------------------------------
        !
        CALL this%ir2r(ir, r, physical) ! position in real space grid
        !
        IF (.NOT. physical) RETURN ! do not include points outside the physical range
        !
        CALL displacement(dim, axis, pos, r, r) ! displacement from origin
        !
        CALL this%minimum_image(r, r2) ! minimum image convention
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE get_min_distance
    !------------------------------------------------------------------------------------
    !>
    !! Returns indices i, j, k yielding the position of grid point ir in the real-space
    !! FFT grid described by descriptor dfft:
    !!
    !! r(:, ir) = i * tau(:, 1) / n1 + j * tau(:, 2) / n2 + k * tau(:, 3) / n3
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE ir2ijk(this, ir, i, j, k, physical)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_cell), INTENT(IN) :: this
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
        k = idx / (this%dfft%nr1x * this%dfft%my_nr2p)
        idx = idx - (this%dfft%nr1x * this%dfft%my_nr2p) * k
        k = k + this%k0
        j = idx / this%dfft%nr1x
        idx = idx - this%dfft%nr1x * j
        j = j + this%j0
        i = idx
        !
        physical = i < this%dfft%nr1 .AND. j < this%dfft%nr2 .AND. k < this%dfft%nr3
        ! check if current point was generated for optimization of fft grids
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE ir2ijk
    !------------------------------------------------------------------------------------
    !>
    !! #TODO unused
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE planar_average(this, nnr, naxis, axis, shift, reverse, f, f1d)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_cell), INTENT(IN) :: this
        INTEGER, INTENT(IN) :: nnr, naxis, axis, shift
        LOGICAL, INTENT(IN) :: reverse
        !
        REAL(DP), INTENT(INOUT) :: f(nnr)
        REAL(DP), INTENT(INOUT) :: f1d(naxis)
        !
        INTEGER :: i, j, k, ir
        INTEGER :: idx, narea
        LOGICAL :: physical
        !
        !--------------------------------------------------------------------------------
        !
        narea = this%ntot / naxis
        !
        IF (reverse) THEN
            f = 0.D0
        ELSE
            f1d = 0.D0
        END IF
        !
        DO ir = 1, this%ir_end
            !
            CALL this%ir2ijk(ir, i, j, k, physical) ! three dimensional indexes
            !
            IF (.NOT. physical) CYCLE ! do not include points outside the physical range
            !
            SELECT CASE (axis)
                !
            CASE (1)
                idx = i
                !
            CASE (2)
                idx = j
                !
            CASE (3)
                idx = k
                !
            END SELECT
            !
            idx = idx + 1 + shift
            !
            IF (idx > naxis) THEN
                idx = idx - naxis
            ELSE IF (idx <= 0) THEN
                idx = idx + naxis
            END IF
            !
            IF (reverse) THEN
                f(ir) = f1d(idx)
            ELSE
                f1d(idx) = f1d(idx) + f(ir)
            END IF
            !
        END DO
        !
        IF (.NOT. reverse) THEN
            !
            CALL env_mp_sum(f1d(:), this%dfft%comm)
            !
            f1d = f1d / DBLE(narea)
        END IF
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE planar_average
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !                               PRIVATE HELPER METHODS
    !
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE ir2r(this, ir, r, physical)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_cell), INTENT(IN) :: this
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
        CALL this%ir2ijk(ir, i, j, k, physical)
        !
        IF (.NOT. physical) RETURN
        !
        DO ip = 1, 3
            !
            r(ip) = DBLE(i) * this%in1 * this%at(ip, 1) + &
                    DBLE(j) * this%in2 * this%at(ip, 2) + &
                    DBLE(k) * this%in3 * this%at(ip, 3)
            !
        END DO
        !
        r = r + this%origin
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
        !
        dr = r1 - r2
        !
        SELECT CASE (dim)
            !
        CASE (1)
            dr(axis) = 0.D0
            !
        CASE (2)
            !
            DO i = 1, 3
                IF (i /= axis) dr(i) = 0.D0
            END DO
            !
        END SELECT
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE displacement
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE minimum_image(this, r, r2)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_cell), INTENT(IN) :: this
        !
        REAL(DP), INTENT(INOUT) :: r(3)
        REAL(DP), INTENT(OUT) :: r2
        !
        INTEGER :: ic
        REAL(DP) :: s(3), rmin(3), r2min
        !
        !--------------------------------------------------------------------------------
        !
        s = MATMUL(r, this%bg)
        s = s - FLOOR(s)
        r = MATMUL(this%at, s)
        !
        rmin = r
        r2min = SUM(r * r)
        !
        DO ic = 2, 8
            s = r + this%corners(:, ic)
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
        !--------------------------------------------------------------------------------
    END SUBROUTINE minimum_image
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    LOGICAL FUNCTION is_cubic(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_cell), INTENT(IN) :: this
        !
        REAL(DP), PARAMETER :: tol = 1.D-8
        !
        INTEGER :: ipol, jpol
        REAL(DP) :: tmp
        !
        !--------------------------------------------------------------------------------
        ! If at(3, 3) is a cubic cell, at(1, 1) = at(2, 2) = at(3, 3) and
        ! the other elements are equal to 0.D0
        !
        tmp = 0.D0
        !
        DO ipol = 1, 3
            !
            DO jpol = 1, 3
                !
                IF (ipol == jpol) THEN
                    tmp = tmp + ABS(this%at(ipol, ipol) - this%at(1, 1))
                ELSE
                    tmp = tmp + ABS(this%at(ipol, jpol))
                END IF
                !
            END DO
            !
        END DO
        !
        is_cubic = tmp < tol
        !
        !--------------------------------------------------------------------------------
    END FUNCTION is_cubic
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !                                   OUTPUT METHODS
    !
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !>
    !! Prints the details of the cell
    !!
    !! If called by a parent object, prints details in block format
    !!
    !! @param verbose       : (INTEGER) adds verbosity to global verbose
    !! @param debug_verbose : (INTEGER) replaces global verbose for debugging
    !! @param unit          : (INTEGER) output target (default = environ_unit)
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE print_environ_cell(this, verbose, debug_verbose, unit)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_cell), INTENT(IN) :: this
        INTEGER, INTENT(IN), OPTIONAL :: verbose, debug_verbose, unit
        !
        INTEGER :: base_verbose, local_verbose, local_unit
        !
        CHARACTER(LEN=80) :: sub_name = 'print_environ_cell'
        !
        !--------------------------------------------------------------------------------
        !
        IF (.NOT. ionode) RETURN
        !
        IF (PRESENT(debug_verbose)) THEN
            base_verbose = debug_verbose
            !
            IF (PRESENT(verbose)) THEN
                local_verbose = verbose
            ELSE
                local_verbose = debug_verbose
            END IF
            !
        ELSE IF (global_verbose > 0) THEN
            base_verbose = global_verbose
            !
            IF (PRESENT(verbose)) THEN
                local_verbose = base_verbose + verbose
            ELSE
                local_verbose = base_verbose
            END IF
            !
        ELSE
            RETURN
        END IF
        !
        IF (PRESENT(unit)) THEN
            local_unit = unit
        ELSE
            local_unit = environ_unit
        END IF
        !
        IF (local_verbose >= 1) THEN
            !
            IF (local_verbose >= base_verbose) THEN ! header
                WRITE (local_unit, 1000)
            ELSE
                !
                CALL env_block_divider(ionode, local_verbose, base_verbose, local_unit)
                !
                WRITE (local_unit, 1001)
            END IF
            !
            WRITE (local_unit, 1002) this%omega
            !
            IF (local_verbose >= 3) THEN
                WRITE (local_unit, 1003) this%at
                WRITE (local_unit, 1004) this%dfft%nr1, this%dfft%nr2, this%dfft%nr3
                WRITE (local_unit, 1005) this%ntot, this%nnr, this%domega
            END IF
            !
            IF (local_verbose < base_verbose) &
                CALL env_block_divider(ionode, local_verbose, base_verbose, local_unit)
            !
        END IF
        !
        FLUSH (local_unit)
        !
        !--------------------------------------------------------------------------------
        !
1000    FORMAT(/, 4('%'), ' CELL ', 70('%'))
1001    FORMAT(/, ' CELL', /, ' ====')
        !
1002    FORMAT(/, ' cell volume                = ', F12.6)
        !
1003    FORMAT(/, ' simulation cell axes       = ', 3F12.6, /, &
                '                              ', 3F12.6, /, &
                '                              ', 3F12.6)
        !
1004    FORMAT(/, ' r-space grid dim           = ', 3I4)
        !
1005    FORMAT(/, ' total size of grid         = ', I12, /, &
                ' r-space size per proc.     = ', I12, /, &
                ' finite element volume      = ', F12.6)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE print_environ_cell
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE write_cube_cell(this, natoms)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_cell), TARGET, INTENT(IN) :: this
        INTEGER, INTENT(IN) :: natoms
        !
        INTEGER :: ipol
        INTEGER :: nr1, nr2, nr3
        !
        REAL(DP), POINTER :: at(:, :), origin(:)
        !
        !--------------------------------------------------------------------------------
        !
        nr1 = this%dfft%nr1
        nr2 = this%dfft%nr2
        nr3 = this%dfft%nr3
        !
        at => this%at
        origin => this%origin
        !
        !--------------------------------------------------------------------------------
        ! Write cube cell data
        !
        WRITE (300, *) 'CUBE FILE GENERATED BY PW.X'
        WRITE (300, *) 'OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z'
        WRITE (300, '(i5,3f12.6)') natoms, origin
        WRITE (300, '(i5,3f12.6)') nr1, (at(ipol, 1) / DBLE(nr1), ipol=1, 3)
        WRITE (300, '(i5,3f12.6)') nr2, (at(ipol, 2) / DBLE(nr2), ipol=1, 3)
        WRITE (300, '(i5,3f12.6)') nr3, (at(ipol, 3) / DBLE(nr3), ipol=1, 3)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE write_cube_cell
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE class_cell
!----------------------------------------------------------------------------------------
