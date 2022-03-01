!----------------------------------------------------------------------------------------
!
! Copyright (C) 2018-2021 ENVIRON (www.quantum-environ.org)
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
!          Francesco Nattino  (THEOS and NCCR-MARVEL, EPFL)
!          Nicola Marzari     (THEOS and NCCR-MARVEL, EPFL)
!          Edan Bainglass     (Department of Physics, UNT)
!
!----------------------------------------------------------------------------------------
!>
!!
!----------------------------------------------------------------------------------------
MODULE class_core_fft
    !------------------------------------------------------------------------------------
    !
    USE class_io, ONLY: io
    USE env_sorting, ONLY: env_hpsort_eps
    USE env_mp, ONLY: env_mp_sum
    !
    USE env_types_fft, ONLY: env_fft_type_descriptor, env_fft_stick_index
    USE env_fft_ggen, ONLY: env_fft_set_nl
    !
    USE environ_param, ONLY: DP, eps8
    !
    USE class_cell
    !
    USE class_core_numerical
    !
    !------------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    PRIVATE
    !
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    TYPE, EXTENDS(numerical_core), PUBLIC :: core_fft
        !--------------------------------------------------------------------------------
        !
        INTEGER :: ngm = 0 ! local  number of G vectors (on this processor)
        ! with gamma tricks, only vectors in G>
        !
        REAL(DP) :: gcutm = 0.0_DP ! ecutrho/(2 pi/a)^2, cut-off for |G|^2
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
    CONTAINS
        !--------------------------------------------------------------------------------
        !
        PROCEDURE :: create => create_core_fft
        PROCEDURE :: init => init_core_fft
        PROCEDURE :: update_cell => update_core_fft_cell
        PROCEDURE :: destroy => destroy_core_fft
        !
        PROCEDURE, PRIVATE :: init_gvect => env_gvect_init
        PROCEDURE, PRIVATE :: deallocate_gvect => env_deallocate_gvect
        !
        !--------------------------------------------------------------------------------
    END TYPE core_fft
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
    SUBROUTINE create_core_fft(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(core_fft), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: sub_name = 'create_core_fft'
        !
        !--------------------------------------------------------------------------------
        !
        IF (ASSOCIATED(this%cell)) CALL io%create_error(sub_name)
        !
        IF (ALLOCATED(this%g)) CALL io%create_error(sub_name)
        !
        IF (ALLOCATED(this%gg)) CALL io%create_error(sub_name)
        !
        !--------------------------------------------------------------------------------
        !
        this%core_type = 'fft'
        this%ngm = 0
        this%gcutm = 0.D0
        this%gstart = 2
        !
        NULLIFY (this%cell)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE create_core_fft
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_core_fft(this, gcutm, cell, use_internal_pbc_corr)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP), INTENT(IN) :: gcutm
        TYPE(environ_cell), TARGET, INTENT(IN) :: cell
        LOGICAL, INTENT(IN), OPTIONAL :: use_internal_pbc_corr
        !
        CLASS(core_fft), INTENT(INOUT) :: this
        !
        INTEGER :: ngm_g ! global number of G vectors (summed on all procs)
        ! in serial execution, ngm_g = ngm
        !
        !--------------------------------------------------------------------------------
        !
        CALL this%create()
        !
        this%gcutm = gcutm
        this%cell => cell
        this%ngm = cell%dfft%ngm
        !
        !--------------------------------------------------------------------------------
        !
        CALL this%init_gvect(ngm_g, cell%dfft%comm)
        !
        CALL env_ggen(this%cell%dfft, cell%dfft%comm, cell%at, cell%bg, this%gcutm, &
                      ngm_g, this%ngm, this%g, this%gg, this%gstart, .TRUE.)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_core_fft
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE update_core_fft_cell(this, cell)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_cell), TARGET, INTENT(IN) :: cell
        !
        CLASS(core_fft), INTENT(INOUT) :: this
        !
        !--------------------------------------------------------------------------------
        !
        this%cell => cell
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE update_core_fft_cell
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE destroy_core_fft(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(core_fft), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: sub_name = 'destroy_core_fft'
        !
        !--------------------------------------------------------------------------------
        !
        IF (.NOT. ASSOCIATED(this%cell)) CALL io%destroy_error(sub_name)
        !
        !--------------------------------------------------------------------------------
        !
        NULLIFY (this%cell)
        !
        CALL this%deallocate_gvect()
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE destroy_core_fft
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !                               PRIVATE HELPER METHODS
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
        INTEGER, ALLOCATABLE :: igsrt(:), g2l(:)
        !
        INTEGER :: ni, nj, nk, i, j, k, ipol, ng, igl, indsw
        INTEGER :: istart, jstart, kstart
        INTEGER :: mype, npe
        LOGICAL :: global_sort, is_local
        INTEGER, ALLOCATABLE :: ngmpe(:)
        !
        INTEGER, EXTERNAL :: env_mp_rank, env_mp_size
        !
        CHARACTER(LEN=80) :: sub_name = 'env_ggen'
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
        gg(:) = gcutm + 1.D0
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
        g2sort_g(:) = 1.0D20
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
                            CALL io%error(sub_name, 'Too many g-vectors', ngm)
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
            CALL io%error(sub_name, 'G-vectors missing!', ABS(ngm - ngm_max))
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
            CALL io%error(sub_name, 'G-vectors (ngm) missing!', ABS(ngm - ngm_save))
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
        CLASS(core_fft), INTENT(INOUT) :: this
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
        ngm = this%ngm
        ngm_g = ngm
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
        CLASS(core_fft), INTENT(INOUT) :: this
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
    !
    !------------------------------------------------------------------------------------
END MODULE class_core_fft
!----------------------------------------------------------------------------------------
