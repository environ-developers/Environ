!----------------------------------------------------------------------------------------
!
! Copyright (C) 2021 ENVIRON (www.quantum-environ.org)
! Copyright (C) Quantum ESPRESSO (www.quantum-espresso.org)
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
! Authors: Modified by Edan Bainglass
!
!----------------------------------------------------------------------------------------
!>
!!
!----------------------------------------------------------------------------------------
MODULE env_base_stick
    !------------------------------------------------------------------------------------
    !
    USE env_base_io, ONLY: io
    USE env_sorting, ONLY: env_hpsort
    !
    USE env_fft_param
    !
    !------------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    TYPE env_sticks_map
        !--------------------------------------------------------------------------------
        !
        LOGICAL :: lgamma = .FALSE. ! if .true. the map has gamma symmetry
        !
        LOGICAL :: lpara = .FALSE.
        ! if .true. the map is set for parallel and serial, if .false. only serial
        !
        INTEGER :: mype = 0 ! my task id (starting from 0)
        INTEGER :: nproc = 1 ! number of task (as nproc in env_fft_type_descriptor)
        !
        INTEGER :: nyfft = 1
        ! number processors in y-direction (as nproc2 in env_fft_type_descriptor)
        !
        INTEGER, ALLOCATABLE :: iproc(:, :)
        ! the processor index (as in env_fft_type_descriptor)
        !
        INTEGER, ALLOCATABLE :: iproc2(:)
        ! the Y group processor index (as in env_fft_type_descriptor)
        !
#if defined(__MPI)
        INTEGER :: comm = MPI_COMM_NULL
#else
        INTEGER :: comm = 0 ! communicator of the fft gruop
#endif
        !
        INTEGER :: nstx = 0 ! a safe maximum number of sticks on the map
        INTEGER :: lb(3) = 0 ! map's lower bounds
        INTEGER :: ub(3) = 0 ! map's upper bounds
        INTEGER, ALLOCATABLE :: idx(:) ! the index of each stick
        INTEGER, ALLOCATABLE :: ist(:, :) ! the cartesian coordinates of each stick
        INTEGER, ALLOCATABLE :: stown(:, :) ! the owner of each stick
        !
        INTEGER, ALLOCATABLE :: indmap(:, :)
        ! the index of each stick (represented on the map)
        !
        REAL(DP) :: bg(3, 3) ! base vectors, the generators of the mapped space
        !
        !--------------------------------------------------------------------------------
    END TYPE
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
    !
    PRIVATE
    !
    PUBLIC :: env_sticks_map, env_sticks_map_allocate, env_sticks_map_deallocate, &
              env_get_sticks
    !
    !------------------------------------------------------------------------------------
CONTAINS
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_sticks_map_allocate(smap, lgamma, lpara, nyfft, iproc, iproc2, nr1, &
                                       nr2, nr3, bg, comm)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        LOGICAL, INTENT(IN) :: lgamma, lpara
        INTEGER, INTENT(IN) :: comm, nyfft
        INTEGER, INTENT(IN) :: iproc(:, :), iproc2(:)
        INTEGER, INTENT(IN) :: nr1, nr2, nr3
        REAL(DP), INTENT(IN) :: bg(3, 3)
        !
        TYPE(env_sticks_map) :: smap
        INTEGER, DIMENSION(3) :: lb, ub
        INTEGER :: nzfft, nstx, ierr
        INTEGER, DIMENSION(:, :), ALLOCATABLE :: indmap, stown, ist
        INTEGER, ALLOCATABLE :: idx(:)
        !
        CHARACTER(LEN=80) :: sub_name = 'env_sticks_map_allocate'
        !
        !--------------------------------------------------------------------------------
        !
        ub(1) = (nr1 - 1) / 2
        ub(2) = (nr2 - 1) / 2
        ub(3) = (nr3 - 1) / 2
        lb = -ub
        nstx = (ub(1) - lb(1) + 1) * (ub(2) - lb(2) + 1) ! we stay very large indeed
        !
        IF (smap%nstx == 0) THEN
            !
            !----------------------------------------------------------------------------
            ! This map is clean, allocate
            !
            smap%mype = 0
            smap%nproc = 1
            smap%comm = comm
            !
#if defined(__MPI)
            CALL MPI_COMM_RANK(smap%comm, smap%mype, ierr)
            !
            CALL MPI_COMM_SIZE(smap%comm, smap%nproc, ierr)
#endif
            !
            smap%lgamma = lgamma
            smap%lpara = lpara
            smap%comm = comm
            smap%nstx = nstx
            smap%ub = ub
            smap%lb = lb
            smap%bg = bg
            smap%nyfft = nyfft
            nzfft = smap%nproc / nyfft
            ALLOCATE (smap%iproc(nyfft, nzfft), smap%iproc2(smap%nproc))
            smap%iproc = iproc
            smap%iproc2 = iproc2
            !
            IF (ALLOCATED(smap%indmap)) &
                CALL io%error(sub_name, 'indmap already allocated', 1)
            !
            IF (ALLOCATED(smap%stown)) &
                CALL io%error(sub_name, 'stown already allocated', 1)
            !
            IF (ALLOCATED(smap%idx)) CALL io%error(sub_name, 'idx already allocated', 1)
            !
            IF (ALLOCATED(smap%ist)) CALL io%error(sub_name, 'ist already allocated', 1)
            !
            ALLOCATE (smap%indmap(lb(1):ub(1), lb(2):ub(2)))
            ALLOCATE (smap%stown(lb(1):ub(1), lb(2):ub(2)))
            ALLOCATE (smap%idx(nstx))
            ALLOCATE (smap%ist(nstx, 2))
            smap%stown = 0
            smap%indmap = 0
            smap%idx = 0
            smap%ist = 0
            !
        ELSE IF (smap%nstx < nstx .OR. smap%ub(3) < ub(3)) THEN
            !
            !----------------------------------------------------------------------------
            ! Change the size of the map, but keep the data already there
            !
            IF (smap%lgamma .NEQV. lgamma) &
                CALL io%error(sub_name, 'Changing gamma symmetry not allowed', 1)
            !
            IF (smap%comm /= comm) &
                CALL io%error(sub_name, 'Changing communicator not allowed', 1)
            !
            ALLOCATE (indmap(lb(1):ub(1), lb(2):ub(2)))
            ALLOCATE (stown(lb(1):ub(1), lb(2):ub(2)))
            ALLOCATE (idx(nstx))
            ALLOCATE (ist(nstx, 2))
            idx = 0
            ist = 0
            indmap = 0
            stown = 0
            idx(1:smap%nstx) = smap%idx
            ist(1:smap%nstx, :) = smap%ist
            !
            indmap(smap%lb(1):smap%ub(1), smap%lb(2):smap%ub(2)) = &
                smap%indmap(smap%lb(1):smap%ub(1), smap%lb(2):smap%ub(2))
            !
            stown(smap%lb(1):smap%ub(1), smap%lb(2):smap%ub(2)) = &
                smap%stown(smap%lb(1):smap%ub(1), smap%lb(2):smap%ub(2))
            !
            DEALLOCATE (smap%indmap)
            DEALLOCATE (smap%stown)
            DEALLOCATE (smap%idx)
            DEALLOCATE (smap%ist)
            ALLOCATE (smap%indmap(lb(1):ub(1), lb(2):ub(2)))
            ALLOCATE (smap%stown(lb(1):ub(1), lb(2):ub(2)))
            ALLOCATE (smap%idx(nstx))
            ALLOCATE (smap%ist(nstx, 2))
            smap%indmap = indmap
            smap%stown = stown
            smap%idx = idx
            smap%ist = ist
            DEALLOCATE (indmap)
            DEALLOCATE (stown)
            DEALLOCATE (idx)
            DEALLOCATE (ist)
            smap%nstx = nstx
            smap%ub = ub
            smap%lb = lb
            smap%bg = bg
            smap%nyfft = nyfft
            smap%iproc = iproc
            smap%iproc2 = iproc2
        ELSE
            !
            IF (smap%lgamma .NEQV. lgamma) &
                CALL io%error(sub_name, 'Changing gamma symmetry not allowed', 2)
            !
            IF (smap%comm /= comm) &
                CALL io%error(sub_name, 'Changing communicator not allowed', 1)
            !
        END IF
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_sticks_map_allocate
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_sticks_map_deallocate(smap)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(env_sticks_map) :: smap
        !
        !--------------------------------------------------------------------------------
        !
        IF (ALLOCATED(smap%iproc)) DEALLOCATE (smap%iproc)
        !
        IF (ALLOCATED(smap%iproc2)) DEALLOCATE (smap%iproc2)
        !
        IF (ALLOCATED(smap%idx)) DEALLOCATE (smap%idx)
        !
        IF (ALLOCATED(smap%ist)) DEALLOCATE (smap%ist)
        !
        IF (ALLOCATED(smap%stown)) DEALLOCATE (smap%stown)
        !
        IF (ALLOCATED(smap%indmap)) DEALLOCATE (smap%indmap)
        !
        smap%ub = 0
        smap%lb = 0
        smap%nstx = 0
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_sticks_map_deallocate
    !------------------------------------------------------------------------------------
    !>
    !! Compute the basic maps of sticks
    !! st(i,j) will contain the number of G vectors of the stick whose indices are (i,j)
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_sticks_map_set(lgamma, ub, lb, bg, gcut, st, comm)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        LOGICAL, INTENT(IN) :: lgamma ! if true use gamma point symmetry
        INTEGER, INTENT(IN) :: ub(:) ! upper bounds for i-th grid dimension
        INTEGER, INTENT(IN) :: lb(:) ! lower bounds for i-th grid dimension
        REAL(DP), INTENT(IN) :: bg(:, :) ! reciprocal space base vectors
        REAL(DP), INTENT(IN) :: gcut ! cut-off for potentials
        INTEGER, OPTIONAL, INTENT(IN) :: comm ! communicator of the g-vec group
        !
        INTEGER, INTENT(OUT) :: st(lb(1):ub(1), lb(2):ub(2))
        ! stick map for wave functions, note that map is taken in YZ plane
        !
        REAL(DP), DIMENSION(3) :: b1, b2, b3
        INTEGER :: i1, i2, i3, n1, n2, n3, mype, nproc, ierr, ngm
        REAL(DP) :: amod
        !
        !--------------------------------------------------------------------------------
        !
        st = 0
        b1(:) = bg(:, 1)
        b2(:) = bg(:, 2)
        b3(:) = bg(:, 3)
        !
        n1 = MAX(ABS(lb(1)), ABS(ub(1)))
        n2 = MAX(ABS(lb(2)), ABS(ub(2)))
        n3 = MAX(ABS(lb(3)), ABS(ub(3)))
        !
        mype = 0
        nproc = 1
        !
#if defined(__MPI)
        IF (PRESENT(comm)) THEN
            !
            CALL MPI_COMM_RANK(comm, mype, ierr)
            !
            CALL MPI_COMM_SIZE(comm, nproc, ierr)
            !
        END IF
#endif
        !
        ngm = 0
        !
        loop1: DO i1 = -n1, n1
            !
            IF ((lgamma .AND. i1 < 0) .OR. (MOD(i1 + n1, nproc) /= mype)) CYCLE loop1
            ! gamma-only: exclude space with x<0
            !
            loop2: DO i2 = -n2, n2
                !
                IF (lgamma .AND. i1 == 0 .AND. i2 < 0) CYCLE loop2
                ! gamma-only: exclude plane with x=0, y<0
                !
                loop3: DO i3 = -n3, n3
                    !
                    !
                    IF (lgamma .AND. i1 == 0 .AND. i2 == 0 .AND. i3 < 0) CYCLE loop3
                    ! gamma-only: exclude line with x=0, y=0, z<0
                    !
                    amod = (i1 * b1(1) + i2 * b2(1) + i3 * b3(1))**2 + &
                           (i1 * b1(2) + i2 * b2(2) + i3 * b3(2))**2 + &
                           (i1 * b1(3) + i2 * b2(3) + i3 * b3(3))**2
                    !
                    IF (amod <= gcut) THEN
                        st(i1, i2) = st(i1, i2) + 1
                        ngm = ngm + 1
                    END IF
                    !
                END DO loop3
                !
            END DO loop2
            !
        END DO loop1
        !
#if defined(__MPI)
        IF (PRESENT(comm)) THEN
            !
            CALL MPI_ALLREDUCE(MPI_IN_PLACE, st, &
                               SIZE(st), MPI_INTEGER, MPI_SUM, comm, ierr)
            !
            CALL MPI_ALLREDUCE(MPI_IN_PLACE, ngm, 1, MPI_INTEGER, MPI_SUM, comm, ierr)
            !
        END IF
#endif
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_sticks_map_set
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_sticks_map_index(ub, lb, st, in1, in2, ngc, index_map)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: ub(:), lb(:)
        INTEGER, INTENT(IN) :: st(lb(1):ub(1), lb(2):ub(2)) ! stick map for potential
        !
        INTEGER, INTENT(INOUT) :: index_map(lb(1):ub(1), lb(2):ub(2))
        ! keep track of sticks index
        !
        INTEGER, INTENT(OUT) :: in1(:), in2(:)
        INTEGER, INTENT(OUT) :: ngc(:)
        !
        INTEGER :: j1, j2, i1, i2, nct, min_size, ind
        !
        CHARACTER(LEN=80) :: sub_name = 'env_sticks_map_index'
        !
        !--------------------------------------------------------------------------------
        ! Initialize the sticks indexes array list
        ! nct counts columns containing G-vectors for the dense grid
        ! ncts counts columns contaning G-vectors for the smooth grid
        !
        nct = MAXVAL(index_map)
        ngc = 0
        !
        min_size = MIN(SIZE(in1), SIZE(in2), SIZE(ngc))
        !
        DO j2 = 0, (ub(2) - lb(2))
            !
            DO j1 = 0, (ub(1) - lb(1))
                i1 = j1
                !
                IF (i1 > ub(1)) i1 = lb(1) + (i1 - ub(1)) - 1
                !
                i2 = j2
                !
                IF (i2 > ub(2)) i2 = lb(2) + (i2 - ub(2)) - 1
                !
                IF (st(i1, i2) > 0) THEN
                    !
                    IF (index_map(i1, i2) == 0) THEN
                        nct = nct + 1
                        index_map(i1, i2) = nct
                    END IF
                    !
                    ind = index_map(i1, i2)
                    !
                    IF (nct > min_size) CALL io%error(sub_name, 'Too many sticks', nct)
                    !
                    in1(ind) = i1
                    in2(ind) = i2
                    ngc(ind) = st(i1, i2)
                END IF
                !
            END DO
            !
        END DO
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_sticks_map_index
    !------------------------------------------------------------------------------------
    !>
    !! This subroutine sorts the sticks indexes, according to the length and type of
    !! the sticks, wave functions sticks first, then smooth mesh sticks, and finally
    !! potential sticks
    !!
    !! lengths of sticks, ngc for potential mesh, ngcw for wave functions mesh and
    !! ngcs for smooth mesh
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_sticks_sort_new(parallel, ng, nct, idx)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        LOGICAL, INTENT(IN) :: parallel
        INTEGER, INTENT(IN) :: ng(:)
        INTEGER, INTENT(IN) :: nct ! nct, total number of sticks
        !
        INTEGER, INTENT(INOUT) :: idx(:) ! index, on output, new sticks indexes
        !
        INTEGER :: mc, ic, nc
        INTEGER, DIMENSION(:), ALLOCATABLE :: iaux, itmp
        REAL(DP), ALLOCATABLE :: aux(:)
        !
        CHARACTER(LEN=80) :: sub_name = 'sticks_sort'
        !
        !--------------------------------------------------------------------------------
        ! We need to avoid sorting elements already sorted previously build
        ! inverse indexes
        !
        ALLOCATE (iaux(nct))
        iaux = 0
        !
        DO mc = 1, nct
            IF (idx(mc) > 0) iaux(idx(mc)) = mc
        END DO
        !
        !--------------------------------------------------------------------------------
        ! Check idx has no "hole"
        !
        IF (idx(1) == 0) THEN
            ic = 0
            !
            DO mc = 2, nct
                !
                IF (idx(mc) /= 0) &
                    CALL io%error(sub_name, 'Non-contiguous indexes 1', nct)
                !
            END DO
            !
        ELSE
            ic = 1
            !
            DO mc = 2, nct
                !
                IF (idx(mc) == 0) EXIT
                !
                ic = ic + 1
            END DO
            !
            DO mc = ic + 1, nct
                !
                IF (idx(mc) /= 0) &
                    CALL io%error(sub_name, 'Non-contiguous indexes 2', nct)
                !
            END DO
            !
        END IF
        !
        IF (parallel) THEN
            ALLOCATE (aux(nct))
            ALLOCATE (itmp(nct))
            itmp = 0
            nc = 0
            !
            DO mc = 1, nct
                !
                IF (ng(mc) > 0 .AND. iaux(mc) == 0) THEN
                    nc = nc + 1
                    aux(nc) = -ng(mc)
                    itmp(nc) = mc
                END IF
                !
            END DO
            !
            CALL env_hpsort(nc, aux, itmp)
            !
            DO mc = 1, nc
                idx(ic + mc) = itmp(mc)
            END DO
            !
            DEALLOCATE (itmp)
            DEALLOCATE (aux)
        ELSE
            !
            DO mc = 1, nct
                !
                IF (ng(mc) > 0 .AND. iaux(mc) == 0) THEN
                    ic = ic + 1
                    idx(ic) = mc
                END IF
                !
            END DO
            !
        END IF
        !
        DEALLOCATE (iaux)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_sticks_sort_new
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_sticks_dist_new(lgamma, mype, nproc, nyfft, iproc, iproc2, ub, lb, &
                                   idx, in1, in2, ngc, nct, ncp, ngp, stown, ng)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        LOGICAL, INTENT(IN) :: lgamma
        INTEGER, INTENT(IN) :: mype, nproc, nyfft
        INTEGER, INTENT(IN) :: iproc(:, :), iproc2(:)
        INTEGER, INTENT(IN) :: ub(:), lb(:) ! map's limits
        !
        INTEGER, INTENT(IN) :: in1(:), in2(:), idx(:)
        ! cartesian coordinate of each column, ordered according to idx index
        INTEGER, INTENT(IN) :: ngc(:) ! number of G-vector of each column, ordered according to idx index
        INTEGER, INTENT(IN) :: nct ! total number of relevant sticks
        !
        INTEGER, INTENT(INOUT) :: stown(lb(1):ub(1), lb(2):ub(2))
        ! stick owner's map, including previously assigned ones
        !
        INTEGER, INTENT(OUT) :: ncp(:) ! number of sticks (columns) per processor
        INTEGER, INTENT(OUT) :: ngp(:) ! number of G-vectors per processor
        !
        INTEGER, INTENT(OUT) :: ng
        ! number of G-vector of this processor. that is ng = ngp(mype+1)
        !
        INTEGER :: mc, i1, i2, j, jj, icnt, gr, j2, j3
        !
        INTEGER, ALLOCATABLE :: yc(:), yg(:)
        ! number of relevant columns and G-vectors in a given yz plane
        !
        INTEGER, ALLOCATABLE :: ygr(:)
        ! element in the nyfft group to which a YZ plane belong
        !
        INTEGER, DIMENSION(:), ALLOCATABLE :: ygrp, ygrc, ygrg
        ! number of yz planes, relevant columns and G-vectors
        ! belonging to a given element in the nyfft group
        !
        !
        CHARACTER(LEN=80) :: sub_name = 'sticks_dist'
        !
        !--------------------------------------------------------------------------------
        ! Distribute X values first
        !
        ALLOCATE (yc(lb(1):ub(1)), yg(lb(1):ub(1)), ygr(lb(1):ub(1)))
        yc = 0
        yg = 0
        ygr = 0
        ALLOCATE (ygrp(nyfft), ygrc(nyfft), ygrg(nyfft))
        ygrp = 0
        ygrc = 0
        ygrg = 0
        !
        DO mc = 1, nct
            !
            IF (idx(mc) < 1) CYCLE
            !
            i1 = in1(idx(mc))
            i2 = in2(idx(mc))
            !
            IF (ngc(idx(mc)) > 0) THEN
                yc(i1) = yc(i1) + 1
                yg(i1) = yg(i1) + ngc(idx(mc))
            END IF
            !
            IF (stown(i1, i2) > 0) THEN
                gr = iproc2(stown(i1, i2))
                !
                IF (ygr(i1) == 0) ygr(i1) = gr
                !
                IF (ygr(i1) .NE. gr) &
                    CALL io%error(sub_name, 'ygroups are not compatible', 1)
                !
            END IF
            !
        END DO
        !
        DO i1 = lb(1), ub(1)
            !
            IF (ygr(i1) == 0) CYCLE
            !
            ygrp(ygr(i1)) = ygrp(ygr(i1)) + 1
            ygrc(ygr(i1)) = ygrc(ygr(i1)) + yc(i1)
            ygrg(ygr(i1)) = ygrg(ygr(i1)) + yg(i1)
        END DO
        !
        ncp = 0
        ngp = 0
        icnt = 0
        !
        DO mc = 1, nct
            !
            IF (idx(mc) < 1) CYCLE
            !
            i1 = in1(idx(mc))
            i2 = in2(idx(mc))
            !
            IF (lgamma .AND. ((i1 < 0) .OR. ((i1 == 0) .AND. (i2 < 0)))) GOTO 30
            !
            IF (ygr(i1) == 0) THEN
                j2 = 1
                !
                DO j = 1, nyfft
                    !
                    IF (ygrg(j) < ygrg(j2)) THEN
                        j2 = j
                    ELSE IF ((ygrg(j) == ygrg(j2)) .AND. (ygrc(j) < ygrc(j2))) THEN
                        j2 = j
                    END IF
                    !
                END DO
                !
                ygr(i1) = j2
                ygrp(j2) = ygrp(j2) + 1
                ygrc(j2) = ygrc(j2) + yc(i1)
                ygrg(j2) = ygrg(j2) + yg(i1)
            ELSE
                j2 = ygr(i1)
            END IF
            !
            IF (ngc(idx(mc)) > 0 .AND. stown(i1, i2) == 0) THEN
                jj = iproc(j2, 1)
                !
                DO j3 = 1, nproc / nyfft
                    j = iproc(j2, j3)
                    !
                    IF (ngp(j) < ngp(jj)) THEN
                        jj = j
                    ELSE IF ((ngp(j) == ngp(jj)) .AND. (ncp(j) < ncp(jj))) THEN
                        jj = j
                    END IF
                    !
                END DO
                !
                stown(i1, i2) = jj
            END IF
            !
            IF (ngc(idx(mc)) > 0) THEN
                ncp(stown(i1, i2)) = ncp(stown(i1, i2)) + 1
                ngp(stown(i1, i2)) = ngp(stown(i1, i2)) + ngc(idx(mc))
            END IF
            !
30          CONTINUE
            !
        END DO
        !
        ng = ngp(mype + 1)
        !
        IF (lgamma) THEN
            !
            !----------------------------------------------------------------------------
            ! When gamma symmetry is used only the sticks of half reciprocal space
            ! are generated, then here we pair-up the sticks with those of the other
            ! half of the space, using the gamma symmetry relation
            ! Note that the total numero of stick "nct" is not modified
            !
            DO mc = 1, nct
                !
                IF (idx(mc) < 1) CYCLE
                !
                IF (ngc(idx(mc)) < 1) CYCLE
                !
                i1 = in1(idx(mc))
                i2 = in2(idx(mc))
                !
                IF (i1 == 0 .AND. i2 == 0) THEN
                    jj = stown(i1, i2)
                    !
                    IF (jj > 0) ngp(jj) = ngp(jj) + ngc(idx(mc)) - 1
                    !
                ELSE
                    jj = stown(i1, i2)
                    !
                    IF (jj > 0) THEN
                        stown(-i1, -i2) = jj
                        ncp(jj) = ncp(jj) + 1
                        ngp(jj) = ngp(jj) + ngc(idx(mc))
                    END IF
                    !
                END IF
                !
            END DO
            !
        END IF
        !
        DEALLOCATE (ygrp, ygrc, ygrg)
        DEALLOCATE (yc, yg, ygr)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_sticks_dist_new
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_get_sticks(smap, gcut, nstp, sstp, st, nst, ng)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP), INTENT(IN) :: gcut
        ! kinetic energy cut-off for this stick distribution
        !
        TYPE(env_sticks_map), INTENT(INOUT) :: smap
        !
        INTEGER, INTENT(OUT) :: st(smap%lb(1):smap%ub(1), smap%lb(2):smap%ub(2))
        ! in output it contains (if>0) the processor_id+1 of the owner of a given stick
        ! internally it's used to contain the number of G-vectors of a given stick
        !
        INTEGER, INTENT(OUT) :: nstp(:) ! number of sticks per processor
        INTEGER, INTENT(OUT) :: sstp(:) ! number of G-vectors per processor
        INTEGER, INTENT(OUT) :: nst ! total number of sticks
        !
        INTEGER, INTENT(OUT) :: ng
        ! number of G-vector of this processor. that is ng = sstp(mype+1)
        !
        INTEGER, ALLOCATABLE :: ngc(:)
        INTEGER :: ic
        !
        CHARACTER(LEN=80) :: sub_name = 'env_get_sticks'
        !
        !--------------------------------------------------------------------------------
        !
        IF (.NOT. ALLOCATED(smap%stown)) &
            CALL io%error(sub_name, 'Sticks map not allocated', 1)
        !
        st = 0
        !
        CALL env_sticks_map_set(smap%lgamma, smap%ub, smap%lb, smap%bg, gcut, st, &
                                smap%comm)
        !
        ALLOCATE (ngc(SIZE(smap%idx)))
        ngc = 0
        !
        CALL env_sticks_map_index(smap%ub, smap%lb, st, smap%ist(:, 1), smap%ist(:, 2), &
                                  ngc, smap%indmap)
        !
        nst = COUNT(st > 0)
        !
        CALL env_sticks_sort_new(smap%nproc > 1, ngc, SIZE(smap%idx), smap%idx)
        !
        CALL env_sticks_dist_new(smap%lgamma, smap%mype, smap%nproc, smap%nyfft, &
                                 smap%iproc, smap%iproc2, smap%ub, smap%lb, smap%idx, &
                                 smap%ist(:, 1), smap%ist(:, 2), ngc, SIZE(smap%idx), &
                                 nstp, sstp, smap%stown, ng)
        !
        !--------------------------------------------------------------------------------
        ! Assign the owner of each (relavant) stick
        !
        st = 0
        !
        DO ic = 1, SIZE(smap%idx)
            !
            IF (smap%idx(ic) > 0) THEN
                !
                IF (ngc(smap%idx(ic)) > 0) THEN
                    !
                    st(smap%ist(smap%idx(ic), 1), smap%ist(smap%idx(ic), 2)) = &
                        smap%stown(smap%ist(smap%idx(ic), 1), smap%ist(smap%idx(ic), 2))
                    !
                    IF (smap%lgamma) &
                        st(-smap%ist(smap%idx(ic), 1), -smap%ist(smap%idx(ic), 2)) = &
                        smap%stown(smap%ist(smap%idx(ic), 1), smap%ist(smap%idx(ic), 2))
                    !
                END IF
                !
            END IF
            !
        END DO
        !
        DEALLOCATE (ngc)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_get_sticks
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE env_base_stick
!----------------------------------------------------------------------------------------
