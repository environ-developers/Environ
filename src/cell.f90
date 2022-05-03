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
    USE class_io, ONLY: io
    USE env_mp, ONLY: env_mp_sum, env_mp_rank, env_mp_size
    !
    USE environ_param, ONLY: DP, tpi, eps8
    !
    USE env_stick_base, ONLY: env_sticks_map, env_sticks_map_deallocate
    !
    USE env_fft_types, ONLY: env_fft_type_descriptor, env_fft_type_init, &
                             env_fft_stick_index, env_fft_type_deallocate
    !
    USE env_fft_ggen, ONLY: env_fft_set_nl
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
        LOGICAL :: lupdate = .FALSE.
        LOGICAL :: cubic = .FALSE.
        !
        CHARACTER(LEN=80) :: label = 'system'
        !
        REAL(DP) :: at(3, 3) ! real-space lattice vectors
        REAL(DP) :: bg(3, 3) ! reciprocal lattice vectors
        !
        REAL(DP) :: origin(3) = 0.D0
        REAL(DP) :: corners(3, 8)
        !
        REAL(DP) :: omega, domega ! volume quantities
        !
        !--------------------------------------------------------------------------------
        ! Properties of the grid
        !
        TYPE(env_fft_type_descriptor) :: dfft
        !
        INTEGER :: nr(3) ! number of grid points along each direction
        INTEGER :: nrx(3) ! number of grid points along each direction (optimized)
        INTEGER :: nnt ! total number of grid points
        INTEGER :: nntx ! total number of grid points (optimized)
        INTEGER :: nnr ! number of grid points allocated in every processor
        INTEGER :: ir_end ! actual number grid points accessed by each processor
        INTEGER :: j0, k0 ! starting indexes of processor-specific boxes of grid points
        REAL(DP) :: in1, in2, in3 ! inverse number of grid points
        !
        REAL(DP) :: gcutm = 0.0_DP ! cutoff for |G|^2
        !
        INTEGER :: gstart = 2 ! index of the first G vector whose module is > 0
        ! needed in parallel execution:
        ! gstart=2 for the proc that holds G=0
        ! gstart=1 for all others
        !
        REAL(DP), ALLOCATABLE :: gg(:)
        ! G^2 in increasing order (in units of (2pi/a)^2)
        !
        REAL(DP), ALLOCATABLE :: g(:, :)
        ! G-vectors cartesian components (in units 2pi/a)
        !
        !--------------------------------------------------------------------------------
        ! Reduced arrays for optimization
        !
        REAL(DP), ALLOCATABLE :: coords(:, :) ! position vectors
        INTEGER, ALLOCATABLE :: ir(:) ! indices of physical grid points
        !
        !--------------------------------------------------------------------------------
    CONTAINS
        !--------------------------------------------------------------------------------
        !
        PROCEDURE, PRIVATE :: create => create_environ_cell
        PROCEDURE :: init => init_environ_cell
        PROCEDURE :: update => update_environ_cell
        PROCEDURE :: destroy => destroy_environ_cell
        !
        PROCEDURE, PRIVATE :: init_dfft
        PROCEDURE, PRIVATE :: destroy_dfft
        !
        PROCEDURE :: volume
        PROCEDURE :: get_min_distance
        PROCEDURE :: ir2ijk
        PROCEDURE :: ir2coords
        PROCEDURE :: map_to_gridx
        PROCEDURE :: planar_average
        PROCEDURE :: running_average
        !
        PROCEDURE, PRIVATE :: minimum_image
        PROCEDURE, PRIVATE :: is_cubic
        !
        PROCEDURE, PRIVATE :: init_gvect => env_gvect_init
        PROCEDURE, PRIVATE :: deallocate_gvect => env_deallocate_gvect
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
        CHARACTER(LEN=80) :: routine = 'create_environ_cell'
        !
        !--------------------------------------------------------------------------------
        !
        IF (ALLOCATED(this%g)) CALL io%create_error(routine)
        !
        IF (ALLOCATED(this%gg)) CALL io%create_error(routine)
        !
        IF (ALLOCATED(this%coords)) CALL io%create_error(routine)
        !
        IF (ALLOCATED(this%ir)) CALL io%create_error(routine)
        !
        !--------------------------------------------------------------------------------
        !
        this%lupdate = .FALSE.
        this%cubic = .FALSE.
        this%label = 'system'
        this%at = 0.D0
        this%bg = 0.D0
        this%origin = 0.D0
        this%corners = 0.D0
        this%omega = 0.D0
        this%domega = 0.D0
        this%nr = 0
        this%nrx = 0
        this%nnt = 0
        this%nntx = 0
        this%nnr = 0
        this%ir_end = 0
        this%j0 = 0
        this%k0 = 0
        this%in1 = 0
        this%in2 = 0
        this%in3 = 0
        this%gcutm = 0.D0
        this%gstart = 2
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE create_environ_cell
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_environ_cell(this, comm, at, gcutm, nr, label)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: comm
        REAL(DP), INTENT(IN) :: at(3, 3), gcutm
        INTEGER, OPTIONAL, INTENT(IN) :: nr(3)
        CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: label
        !
        CLASS(environ_cell), INTENT(INOUT) :: this
        !
        INTEGER :: ngm_g ! global number of G vectors (summed on all procs)
        ! in serial execution, ngm_g = ngm
        !
        INTEGER :: i
        LOGICAL :: physical
        REAL(DP) :: coords(3)
        !
        CHARACTER(LEN=80) :: routine = 'init_environ_cell'
        !
        !--------------------------------------------------------------------------------
        !
        CALL this%create()
        !
        IF (PRESENT(label)) this%label = label
        !
        !--------------------------------------------------------------------------------
        ! Create fft descriptor for system cell
        !
        IF (PRESENT(nr)) THEN
            this%dfft%nr1 = nr(1)
            this%dfft%nr2 = nr(2)
            this%dfft%nr3 = nr(3)
        END IF
        !
        CALL this%init_dfft(gcutm, comm, at)
        !
        !--------------------------------------------------------------------------------
        ! Real-space grid details
        !
        this%nr = (/this%dfft%nr1, this%dfft%nr2, this%dfft%nr3/)
        this%nrx = (/this%dfft%nr1x, this%dfft%nr2x, this%dfft%nr3x/)
        !
        this%nnt = PRODUCT(this%nr)
        this%nntx = PRODUCT(this%nrx)
        !
        this%nnr = this%dfft%nnr
        !
#if defined(__MPI)
        this%j0 = this%dfft%my_i0r2p
        this%k0 = this%dfft%my_i0r3p
        !
        this%ir_end = MIN(this%nnr, this%nrx(1) * this%dfft%my_nr2p * this%dfft%my_nr3p)
#else
        this%j0 = 0
        this%k0 = 0
        this%ir_end = this%nnr
#endif
        !
        this%in1 = 1.D0 / DBLE(this%nr(1))
        this%in2 = 1.D0 / DBLE(this%nr(2))
        this%in3 = 1.D0 / DBLE(this%nr(3))
        !
        !--------------------------------------------------------------------------------
        ! Set basic cell properties
        !
        CALL this%update(at)
        !
        !--------------------------------------------------------------------------------
        !
        CALL this%init_gvect(ngm_g, this%dfft%comm)
        !
        CALL env_ggen(this%dfft, this%dfft%comm, this%at, this%bg, this%gcutm, &
                      ngm_g, this%dfft%ngm, this%g, this%gg, this%gstart, .TRUE.)
        !
        !--------------------------------------------------------------------------------
        ! Storing position vectors for optimization
        !
        ALLOCATE (this%coords(3, this%nnr))
        ALLOCATE (this%ir(this%nnr))
        this%coords = 0.D0
        this%ir = 0
        !
        DO i = 1, this%nnr
            !
            CALL this%ir2coords(i, coords, physical)
            !
            IF (.NOT. physical) CYCLE
            !
            this%coords(:, i) = coords
            this%ir(i) = i
        END DO
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
        INTEGER :: i, j, k, l
        REAL(DP) :: dx, dy, dz
        !
        CHARACTER(LEN=80) :: routine = 'update_environ_cell'
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
        this%domega = this%omega / this%nnt ! set volume element
        !
        !--------------------------------------------------------------------------------
        ! Calculate corners for minimum image convention
        !
        l = 0
        !
        DO i = 0, 1
            dx = DBLE(-i)
            !
            DO j = 0, 1
                dy = DBLE(-j)
                !
                DO k = 0, 1
                    dz = DBLE(-k)
                    l = l + 1
                    this%corners(1, l) = dx * at(1, 1) + dy * at(1, 2) + dz * at(1, 3)
                    this%corners(2, l) = dx * at(2, 1) + dy * at(2, 2) + dz * at(2, 3)
                    this%corners(3, l) = dx * at(3, 1) + dy * at(3, 2) + dz * at(3, 3)
                END DO
                !
            END DO
            !
        END DO
        !
        !--------------------------------------------------------------------------------
        ! Output current state
        !
        CALL this%printout()
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE update_environ_cell
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE destroy_environ_cell(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_cell), INTENT(INOUT) :: this
        !
        !--------------------------------------------------------------------------------
        !
        CALL this%destroy_dfft()
        !
        CALL this%deallocate_gvect()
        !
        IF (ALLOCATED(this%coords)) DEALLOCATE (this%coords)
        IF (ALLOCATED(this%ir)) DEALLOCATE (this%ir)
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
        this%gcutm = gcutm
        !
        CALL recips(at(1, 1), at(1, 2), at(1, 3), bg(1, 1), bg(1, 2), bg(1, 3))
        ! calculate the reciprocal lattice vectors
        !
        CALL env_fft_type_init(this%dfft, smap, 'rho', .TRUE., .TRUE., comm, at, bg, &
                               gcutm, nyfft=1, nmany=1)
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
    SUBROUTINE destroy_dfft(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
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
        CLASS(environ_cell), INTENT(INOUT) :: this
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
            CALL io%warning("axis vectors are left-handed", 1004)
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
        INTEGER :: i, j, k, l
        !
        INTEGER :: m ! counter on the permutations
        INTEGER :: n ! counter on the polarizations
        !
        REAL(DP) :: denominator
        REAL(DP) :: sign
        !
        !--------------------------------------------------------------------------------
        ! Compute the denominator
        !
        denominator = 0
        i = 1
        j = 2
        k = 3
        sign = 1.D0
        !
        DO WHILE (sign >= 0.D0)
            !
            DO m = 1, 3
                denominator = denominator + sign * a1(i) * a2(j) * a3(k)
                l = i
                i = j
                j = k
                k = l
            END DO
            !
            i = 2
            j = 1
            k = 3
            sign = -sign
            !
        END DO
        !
        !--------------------------------------------------------------------------------
        ! Compute the reciprocal vectors
        !
        i = 1
        j = 2
        k = 3
        !
        DO n = 1, 3
            b1(n) = (a2(j) * a3(k) - a2(k) * a3(j)) / denominator
            b2(n) = (a3(j) * a1(k) - a3(k) * a1(j)) / denominator
            b3(n) = (a1(j) * a2(k) - a1(k) * a2(j)) / denominator
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
    SUBROUTINE get_min_distance(this, ir, dim, axis, origin, r, r2, physical)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_cell), INTENT(IN) :: this
        INTEGER, INTENT(IN) :: ir, dim, axis
        REAL(DP), INTENT(IN) :: origin(3)
        !
        LOGICAL, INTENT(INOUT) :: physical
        REAL(DP), INTENT(OUT) :: r(3), r2
        !
        REAL(DP) :: coords(3)
        !
        !--------------------------------------------------------------------------------
        !
        physical = .TRUE.
        !
        IF (this%ir(ir) == 0) THEN
            !
            physical = .FALSE.
            !
            RETURN
            !
        END IF
        !
        coords = this%coords(:, ir)
        !
        CALL displacement(dim, axis, coords, origin, r) ! displacement from origin
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
        k = idx / (this%nrx(1) * this%dfft%my_nr2p)
        idx = idx - (this%nrx(1) * this%dfft%my_nr2p) * k
        k = k + this%k0
        j = idx / this%nrx(1)
        idx = idx - this%nrx(1) * j
        j = j + this%j0
        i = idx
        !
        physical = i < this%nr(1) .AND. j < this%nr(2) .AND. k < this%nr(3)
        ! check if current point was generated for optimization of fft grids
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE ir2ijk
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE ir2coords(this, ir, coords, physical)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_cell), INTENT(IN) :: this
        INTEGER, INTENT(IN) :: ir
        !
        REAL(DP), INTENT(OUT) :: coords(3)
        LOGICAL, INTENT(OUT) :: physical
        !
        INTEGER :: i, j, k, l
        !
        !--------------------------------------------------------------------------------
        !
        coords = 0.D0
        !
        CALL this%ir2ijk(ir, i, j, k, physical)
        !
        IF (.NOT. physical) RETURN
        !
        DO l = 1, 3
            !
            coords(l) = DBLE(i) * this%in1 * this%at(l, 1) + &
                        DBLE(j) * this%in2 * this%at(l, 2) + &
                        DBLE(k) * this%in3 * this%at(l, 3)
            !
        END DO
        !
        coords = coords + this%origin
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE ir2coords
    !------------------------------------------------------------------------------------
    !>
    !! Map array onto parallelization-optimized grid
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE map_to_gridx(this, array_in, array_out)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_cell), TARGET, INTENT(IN) :: this
        REAL(DP), INTENT(IN) :: array_in(:)
        !
        REAL(DP), ALLOCATABLE, INTENT(OUT) :: array_out(:)
        !
        INTEGER :: ir, i, j, k
        LOGICAL :: physical
        !
        INTEGER, POINTER :: nntx
        !
        CHARACTER(LEN=80) :: routine = 'map_to_gridx'
        !
        !--------------------------------------------------------------------------------
        !
        nntx => this%nntx
        !
        IF (SIZE(array_in) /= nntx) CALL io%error(routine, "wrong input array size", 1)
        !
        !--------------------------------------------------------------------------------
        !
        ALLOCATE (array_out(nntx))
        !
        array_out = 0.D0
        !
        DO ir = 1, nntx
            !
            CALL this%ir2ijk(ir, i, j, k, physical)
            !
            IF (.NOT. physical) CYCLE
            !
            array_out(ir) = array_in(ir)
        END DO
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE map_to_gridx
    !------------------------------------------------------------------------------------
    !>
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
        CHARACTER(LEN=80) :: routine = 'planar_average'
        !
        !--------------------------------------------------------------------------------
        !
        narea = this%nnt / naxis
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
            CASE DEFAULT
                CALL io%error(routine, "Unexpected axis value", 1)
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
            CALL env_mp_sum(f1d, this%dfft%comm)
            !
            f1d = f1d / DBLE(narea)
        END IF
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE planar_average
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE running_average(this, naxis, window, pot, averaged_pot)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_cell), INTENT(IN) :: this
        INTEGER, INTENT(IN) :: naxis, window
        REAL(DP), INTENT(IN) :: pot(naxis)
        !
        REAL(DP), INTENT(OUT) :: averaged_pot(naxis)
        !
        INTEGER :: i, start_idx, stop_idx
        !
        !--------------------------------------------------------------------------------
        ! Averaging bb
        !
        DO i = 1, naxis
            start_idx = i - window
            stop_idx = i + window
            !
            IF (start_idx < 1) start_idx = 1
            !
            IF (stop_idx > naxis) stop_idx = naxis
            !
            averaged_pot(i) = SUM(pot(start_idx:stop_idx)) / FLOAT(stop_idx - start_idx)
        END DO
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE running_average
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
    SUBROUTINE displacement(dim, axis, r1, r2, r)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: dim, axis
        REAL(DP), DIMENSION(3), INTENT(IN) :: r1, r2
        !
        REAL(DP), INTENT(OUT) :: r(3)
        !
        INTEGER :: i
        !
        CHARACTER(LEN=80) :: routine = 'displacement'
        !
        !--------------------------------------------------------------------------------
        !
        r = r1 - r2
        !
        SELECT CASE (dim)
            !
        CASE (0)
            !
        CASE (1)
            r(axis) = 0.D0
            !
        CASE (2)
            !
            DO i = 1, 3
                IF (i /= axis) r(i) = 0.D0
            END DO
            !
        CASE DEFAULT
            CALL io%error(routine, "Unexpected system dimensions", 1)
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
        INTEGER :: i
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
        DO i = 2, 8
            s = r + this%corners(:, i)
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
        ! x = x - NINT(x)
        ! c = SUM(x * MATMUL(ws%aa, x))
        ! m = 0
        ! !
        ! lb = NINT(x - DSQRT(c) * ws%norm_b)
        ! ! CEILING should be enough for lb but NINT might be safer
        ! !
        ! ub = NINT(x + DSQRT(c) * ws%norm_b)
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
        INTEGER :: i, j
        REAL(DP) :: tmp
        !
        !--------------------------------------------------------------------------------
        ! If at(3, 3) is a cubic cell, at(1, 1) = at(2, 2) = at(3, 3) and
        ! the other elements are equal to 0.D0
        !
        tmp = 0.D0
        !
        DO i = 1, 3
            !
            DO j = 1, 3
                !
                IF (i == j) THEN
                    tmp = tmp + ABS(this%at(i, i) - this%at(1, 1))
                ELSE
                    tmp = tmp + ABS(this%at(i, j))
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
    !                                  G-SPACE METHODS
    !
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !>
    !! This routine generates all the reciprocal lattice vectors
    !! contained in the sphere of radius gcutm. Furthermore it
    !! computes the indices nl which give the correspondence
    !! between the fft mesh points and the array of g vectors.
    !
    !------------------------------------------------------------------------------------
    SUBROUTINE env_ggen(dfftp, comm, at, bg, gcutm, ngm_g, ngm, g, gg, gstart, &
                        no_global_sort)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP), INTENT(IN) :: at(3, 3), bg(3, 3), gcutm
        INTEGER, INTENT(IN) :: ngm_g, comm
        LOGICAL, OPTIONAL, INTENT(IN) :: no_global_sort
        ! if no_global_sort is present (and it is true) G vectors are sorted only
        ! locally and not globally. In this case no global array needs to be
        ! allocated and sorted: saves memory and a lot of time for large systems
        !
        TYPE(env_fft_type_descriptor), INTENT(INOUT) :: dfftp
        INTEGER, INTENT(INOUT) :: ngm
        REAL(DP), INTENT(OUT) :: g(:, :), gg(:)
        INTEGER, INTENT(OUT) :: gstart
        !
        REAL(DP) :: tx(3), ty(3), t(3)
        REAL(DP), ALLOCATABLE :: tt(:)
        INTEGER :: ngm_save, n1, n2, n3, ngm_offset, ngm_max, ngm_local
        !
        REAL(DP), ALLOCATABLE :: g2sort_g(:)
        ! array containing only g vectors for the current processor
        !
        INTEGER, ALLOCATABLE :: mill_unsorted(:, :)
        ! array containing all g vectors generators, on all processors
        ! (replicated data). When no_global_sort is present and .true.,
        ! only g-vectors for the current processor are stored
        !
        INTEGER, DIMENSION(:), ALLOCATABLE :: igsrt, g2l
        !
        INTEGER :: ni, nj, nk, i, j, k, ng
        INTEGER :: istart, jstart, kstart
        INTEGER :: mype, npe
        LOGICAL :: global_sort, is_local
        INTEGER, ALLOCATABLE :: ngmpe(:)
        !
        CHARACTER(LEN=80) :: routine = 'env_ggen'
        !
        !--------------------------------------------------------------------------------
        !
        global_sort = .TRUE.
        !
        IF (PRESENT(no_global_sort)) global_sort = .NOT. no_global_sort
        !
        IF (.NOT. global_sort) THEN
            ngm_max = ngm
        ELSE
            ngm_max = ngm_g
        END IF
        !
        ngm_save = ngm ! save current value of ngm
        !
        ngm = 0
        ngm_local = 0
        !
        gg = gcutm + 1.D0
        ! set the total number of fft mesh points and and initial value of gg
        ! The choice of gcutm is due to the fact that we have to order the
        ! vectors after computing them.
        !
        !--------------------------------------------------------------------------------
        ! Computes all the g vectors inside a sphere
        !
        ALLOCATE (mill_unsorted(3, ngm_save))
        ALLOCATE (igsrt(ngm_max))
        ALLOCATE (g2l(ngm_max))
        ALLOCATE (g2sort_g(ngm_max))
        !
        g2sort_g = 1.0D20
        !
        !--------------------------------------------------------------------------------
        !
        ALLOCATE (tt(dfftp%nr3)) ! allocate temporal array
        !
        !--------------------------------------------------------------------------------
        ! Max miller indices (same convention as in module stick_set)
        !
        ni = (dfftp%nr1 - 1) / 2
        nj = (dfftp%nr2 - 1) / 2
        nk = (dfftp%nr3 - 1) / 2
        !
        !--------------------------------------------------------------------------------
        ! Gamma-only: exclude space with x < 0
        !
        istart = 0
        !
        iloop: DO i = istart, ni
            !
            !----------------------------------------------------------------------------
            ! Gamma-only: exclude plane with x = 0, y < 0
            !
            IF (i == 0) THEN
                jstart = 0
            ELSE
                jstart = -nj
            END IF
            !
            tx(1:3) = i * bg(1:3, 1)
            !
            jloop: DO j = jstart, nj
                !
                IF (.NOT. global_sort) THEN
                    !
                    IF (env_fft_stick_index(dfftp, i, j) == 0) CYCLE jloop
                    !
                    is_local = .TRUE.
                ELSE
                    !
                    IF (dfftp%lpara .AND. env_fft_stick_index(dfftp, i, j) == 0) THEN
                        is_local = .FALSE.
                    ELSE
                        is_local = .TRUE.
                    END IF
                    !
                END IF
                !
                !------------------------------------------------------------------------
                ! Gamma-only: exclude line with x = 0, y = 0, z < 0
                !
                IF (i == 0 .AND. j == 0) THEN
                    kstart = 0
                ELSE
                    kstart = -nk
                END IF
                !
                ty(1:3) = tx(1:3) + j * bg(1:3, 2)
                !
                !------------------------------------------------------------------------
                ! Compute all the norm square
                !
                DO k = kstart, nk
                    t(1) = ty(1) + k * bg(1, 3)
                    t(2) = ty(2) + k * bg(2, 3)
                    t(3) = ty(3) + k * bg(3, 3)
                    tt(k - kstart + 1) = t(1)**2 + t(2)**2 + t(3)**2
                END DO
                !
                !------------------------------------------------------------------------
                ! Save all the norm square within cutoff
                !
                DO k = kstart, nk
                    !
                    IF (tt(k - kstart + 1) <= gcutm) THEN
                        ngm = ngm + 1
                        !
                        IF (ngm > ngm_max) &
                            CALL io%error(routine, "Too many g-vectors", ngm)
                        !
                        IF (tt(k - kstart + 1) > eps8) THEN
                            g2sort_g(ngm) = tt(k - kstart + 1)
                        ELSE
                            g2sort_g(ngm) = 0.D0
                        END IF
                        !
                        IF (is_local) THEN
                            ngm_local = ngm_local + 1
                            mill_unsorted(:, ngm_local) = (/i, j, k/)
                            g2l(ngm) = ngm_local
                        ELSE
                            g2l(ngm) = 0
                        END IF
                        !
                    END IF
                    !
                END DO
                !
            END DO jloop
            !
        END DO iloop
        !
        IF (ngm /= ngm_max) &
            CALL io%error(routine, "G-vectors missing!", ABS(ngm - ngm_max))
        !
        igsrt(1) = 0
        !
        IF (.NOT. global_sort) THEN
            CALL env_hpsort_eps(ngm, g2sort_g, igsrt, eps8)
        ELSE
            CALL env_hpsort_eps(ngm_g, g2sort_g, igsrt, eps8)
        END IF
        !
        DEALLOCATE (g2sort_g, tt)
        !
        IF (.NOT. global_sort) THEN
            !
            !----------------------------------------------------------------------------
            ! Compute adequate offsets in order to avoid overlap between
            ! g vectors once they are gathered on a single (global) array
            !
            mype = env_mp_rank(comm)
            npe = env_mp_size(comm)
            ALLOCATE (ngmpe(npe))
            ngmpe = 0
            ngmpe(mype + 1) = ngm
            !
            CALL env_mp_sum(ngmpe, comm)
            !
            ngm_offset = 0
            !
            DO ng = 1, mype
                ngm_offset = ngm_offset + ngmpe(ng)
            END DO
            !
            DEALLOCATE (ngmpe)
            !
        END IF
        !
        ngm = 0
        !
        ngloop: DO ng = 1, ngm_max
            !
            IF (g2l(igsrt(ng)) > 0) THEN
                !
                !------------------------------------------------------------------------
                ! Fetch the indices
                !
                i = mill_unsorted(1, g2l(igsrt(ng)))
                j = mill_unsorted(2, g2l(igsrt(ng)))
                k = mill_unsorted(3, g2l(igsrt(ng)))
                !
                ngm = ngm + 1
                !
                !------------------------------------------------------------------------
                ! Map local and global g index
                ! N.B: the global G vectors arrangement depends on the number of processors
                !
                g(1:3, ngm) = i * bg(:, 1) + j * bg(:, 2) + k * bg(:, 3)
                gg(ngm) = SUM(g(1:3, ngm)**2)
            END IF
            !
        END DO ngloop
        !
        DEALLOCATE (igsrt, g2l)
        !
        IF (ngm /= ngm_save) &
            CALL io%error(routine, "G-vectors (ngm) missing!", ABS(ngm - ngm_save))
        !
        !--------------------------------------------------------------------------------
        ! Determine first nonzero g vector
        !
        IF (gg(1) <= eps8) THEN
            gstart = 2
        ELSE
            gstart = 1
        END IF
        !
        !--------------------------------------------------------------------------------
        !
        CALL env_fft_set_nl(dfftp, at, g)
        ! set nl and nls with the correct fft correspondence
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_ggen
    !------------------------------------------------------------------------------------
    !>
    !! Set local and global dimensions, allocate arrays
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_gvect_init(this, ngm_g, comm)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_cell), INTENT(INOUT) :: this
        INTEGER, INTENT(INOUT) :: ngm_g
        !
        INTEGER :: ngm
        !
        INTEGER, INTENT(IN) :: comm
        ! communicator of the group on which g-vecs are distributed
        !
        !--------------------------------------------------------------------------------
        ! Calculate sum over all processors
        !
        ngm = this%dfft%ngm ! local
        ngm_g = ngm ! global
        !
        CALL env_mp_sum(ngm_g, comm)
        !
        !--------------------------------------------------------------------------------
        ! Allocate arrays - only those that are always kept until the end
        !
        ALLOCATE (this%gg(ngm))
        ALLOCATE (this%g(3, ngm))
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_gvect_init
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_deallocate_gvect(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_cell), INTENT(INOUT) :: this
        !
        !--------------------------------------------------------------------------------
        !
        IF (ALLOCATED(this%gg)) DEALLOCATE (this%gg)
        !
        IF (ALLOCATED(this%g)) DEALLOCATE (this%g)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_deallocate_gvect
    !------------------------------------------------------------------------------------
    !>
    !! Sort an array ra(1:n) into ascending order using heapsort algorithm,
    !! and considering two elements being equal if their values differ
    !! for less than "eps".
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_hpsort_eps(n, ra, ind, eps)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: n
        REAL(DP), INTENT(IN) :: eps
        !
        INTEGER, INTENT(INOUT) :: ind(*)
        REAL(DP), INTENT(INOUT) :: ra(*)
        !
        INTEGER :: i, ir, j, l, iind
        REAL(DP) :: rra
        !
        !--------------------------------------------------------------------------------
        ! Initialize index array
        !
        IF (ind(1) == 0) THEN
            !
            DO i = 1, n
                ind(i) = i
            END DO
            !
        END IF
        !
        !--------------------------------------------------------------------------------
        !
        IF (n < 2) RETURN ! nothing to order
        !
        l = n / 2 + 1
        ir = n
        !
        sorting: DO
            !
            IF (l > 1) THEN ! hiring phase
                !
                l = l - 1
                rra = ra(l)
                iind = ind(l)
                !
            ELSE ! retirement-promotion phase
                !
                rra = ra(ir)
                iind = ind(ir)
                ! clear a space at the end of the array
                !
                ra(ir) = ra(1)
                ind(ir) = ind(1)
                ! retire the top of the heap into it
                !
                ir = ir - 1 ! decrease the size of the corporation
                !
                IF (ir == 1) THEN ! done with the last promotion
                    !
                    ra(1) = rra
                    ind(1) = iind
                    ! the least competent worker
                    !
                    EXIT sorting
                    !
                END IF
                !
            END IF
            !
            !----------------------------------------------------------------------------
            ! Regardless of phase, we prepare to place rra in its proper level
            !
            i = l
            j = l + l
            !
            DO WHILE (j <= ir)
                !
                IF (j < ir) THEN
                    !
                    IF (ABS(ra(j) - ra(j + 1)) >= eps) THEN
                        !
                        IF (ra(j) < ra(j + 1)) j = j + 1
                        ! compare to better underling
                        !
                    ELSE ! this means ra(j) == ra(j+1) within tolerance
                        !
                        IF (ind(j) < ind(j + 1)) j = j + 1
                        !
                    END IF
                    !
                END IF
                !
                IF (ABS(rra - ra(j)) >= eps) THEN ! demote rra
                    !
                    IF (rra < ra(j)) THEN
                        ra(i) = ra(j)
                        ind(i) = ind(j)
                        i = j
                        j = j + j
                    ELSE
                        j = ir + 1 ! set j to terminate do-while loop
                    END IF
                    !
                ELSE ! this means rra == ra(j) within tolerance
                    !
                    IF (iind < ind(j)) THEN ! demote rra
                        ra(i) = ra(j)
                        ind(i) = ind(j)
                        i = j
                        j = j + j
                    ELSE
                        j = ir + 1 ! set j to terminate do-while loop
                    END IF
                    !
                END IF
                !
            END DO
            !
            ra(i) = rra
            ind(i) = iind
            !
        END DO sorting
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_hpsort_eps
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
    !! @param unit          : (INTEGER) output target (default = io%debug_unit)
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE print_environ_cell(this, verbose, debug_verbose, unit)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_cell), INTENT(IN) :: this
        INTEGER, OPTIONAL, INTENT(IN) :: verbose, debug_verbose, unit
        !
        INTEGER :: base_verbose, local_verbose, local_unit
        !
        CHARACTER(LEN=80) :: routine = 'print_environ_cell'
        !
        !--------------------------------------------------------------------------------
        !
        IF (.NOT. io%lnode) RETURN
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
        ELSE IF (io%verbosity > 0) THEN
            base_verbose = io%verbosity
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
            local_unit = io%debug_unit
        END IF
        !
        IF (local_verbose >= 1) THEN
            !
            IF (local_verbose >= base_verbose) THEN ! header
                WRITE (local_unit, 1000)
            ELSE
                !
                CALL io%block_divider(local_verbose, base_verbose, local_unit)
                !
                WRITE (local_unit, 1001)
            END IF
            !
            WRITE (local_unit, 1002) this%label
            WRITE (local_unit, 1003) this%omega
            !
            IF (local_verbose >= 3) THEN
                WRITE (local_unit, 1004) this%at
                WRITE (local_unit, 1005) this%nr
                WRITE (local_unit, 1006) this%nnt, this%nnr, this%domega
            END IF
            !
            IF (local_verbose < base_verbose) &
                CALL io%block_divider(local_verbose, base_verbose, local_unit)
            !
        END IF
        !
        FLUSH (local_unit)
        !
        !--------------------------------------------------------------------------------
        !
1000    FORMAT(/, 4('%'), " CELL ", 70('%'))
1001    FORMAT(/, " CELL", /, " ====")
        !
1002    FORMAT(/, " cell label                 = ", A15)
        !
1003    FORMAT(/, " cell volume                = ", F14.6)
        !
1004    FORMAT(/, " simulation cell axes       = ", 3F14.6, /, &
                "                              ", 3F14.6, /, &
                "                              ", 3F14.6)
        !
1005    FORMAT(/, " r-space grid dim           = ", 3I14)
        !
1006    FORMAT(/, " total size of grid         = ", I14, /, &
                " r-space size per proc.     = ", I14, /, &
                " finite element volume      = ", F14.6)
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
        INTEGER :: i
        !
        !--------------------------------------------------------------------------------
        !
        ASSOCIATE (nr => this%nr, &
                   at => this%at, &
                   origin => this%origin)
            !
            !----------------------------------------------------------------------------
            ! Write cube cell data
            !
            WRITE (300, *) "CUBE FILE GENERATED BY PW.X"
            WRITE (300, *) "OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z"
            WRITE (300, "(i5,3f12.6)") natoms, origin
            WRITE (300, "(i5,3f12.6)") nr(1), (at(i, 1) / DBLE(nr(1)), i=1, 3)
            WRITE (300, "(i5,3f12.6)") nr(2), (at(i, 2) / DBLE(nr(2)), i=1, 3)
            WRITE (300, "(i5,3f12.6)") nr(3), (at(i, 3) / DBLE(nr(3)), i=1, 3)
            !
        END ASSOCIATE
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE write_cube_cell
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE class_cell
!----------------------------------------------------------------------------------------
