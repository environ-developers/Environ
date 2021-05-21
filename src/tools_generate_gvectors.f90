!
! Copyright (C) 2011 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------------------
!
!> subroutines generating G-vectors and variables nl* needed to map
!! G-vector components onto the FFT grid(s) in reciprocal space
!!
!----------------------------------------------------------------------------------------
MODULE tools_generate_gvectors
    !------------------------------------------------------------------------------------
    !
    USE env_sorting, ONLY: env_hpsort_eps
    USE env_mp, ONLY: env_mp_sum, env_mp_rank, env_mp_size
    !
    USE env_fft_types, ONLY: env_fft_type_descriptor, env_fft_stick_index
    USE env_fft_ggen, ONLY: env_fft_set_nl
    !
    USE modules_constants, ONLY: DP, eps8
    !
    USE core_types, ONLY: fft_core
    !
    !------------------------------------------------------------------------------------
    !
    PRIVATE
    !
    PUBLIC :: env_gvect_init, env_ggen
    !
    !------------------------------------------------------------------------------------
CONTAINS
    !------------------------------------------------------------------------------------
    !>
    !! Set local and global dimensions, allocate arrays
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_gvect_init(fft, ngm_g, comm)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(fft_core), INTENT(INOUT) :: fft
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
        ngm = fft%ngm
        ngm_g = ngm
        !
        CALL env_mp_sum(ngm_g, comm)
        !
        !--------------------------------------------------------------------------------
        ! Allocate arrays - only those that are always kept until the end
        !
        ALLOCATE (fft%gg(ngm))
        ALLOCATE (fft%g(3, ngm))
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_gvect_init
    !------------------------------------------------------------------------------------
    !>
    !! #TODO unused
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_deallocate_gvect(fft)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(fft_core), INTENT(INOUT) :: fft
        !
        !--------------------------------------------------------------------------------
        !
        IF (ALLOCATED(fft%gg)) DEALLOCATE (fft%gg)
        !
        IF (ALLOCATED(fft%g)) DEALLOCATE (fft%g)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_deallocate_gvect
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
        CHARACTER(LEN=80) :: sub_name = 'env_ggen'
        !
        !--------------------------------------------------------------------------------
        !
        global_sort = .TRUE.
        !
        IF (PRESENT(no_global_sort)) THEN
            global_sort = .NOT. no_global_sort
        END IF
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
                            CALL env_errore(sub_name, 'too many g-vectors', ngm)
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
            CALL env_errore(sub_name, 'g-vectors missing !', ABS(ngm - ngm_max))
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
            CALL env_errore(sub_name, 'g-vectors (ngm) missing !', ABS(ngm - ngm_save))
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
    !
    !------------------------------------------------------------------------------------
END MODULE tools_generate_gvectors
!----------------------------------------------------------------------------------------
