!----------------------------------------------------------------------------------------
!
! Copyright (C) 2018 ENVIRON (www.quantum-environment.org)
! Copyright (C) 2011 Quantum ESPRESSO group
!
!----------------------------------------------------------------------------------------
!
! This file is part of Environ version 2.0
!
! Environ 2.0 is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 2 of the License, or
! (at your option) any later version.
!
! Environ 2.0 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more detail, either the file
! `License' in the root directory of the present distribution, or
! online at <http://www.gnu.org/licenses/>.
!
!----------------------------------------------------------------------------------------
!
! Authors:
!
!----------------------------------------------------------------------------------------
!>
!! Sorting algorithms
!!
!! n is input, ra is replaced on output by its sorted rearrangement.
!!
!! Create an index table (ind) by making an exchange in the index array
!! whenever an exchange is made on the sorted data array (ra).
!!
!! In case of equal values in the data array (ra) the values in the
!! index array (ind) are used to order the entries.
!!
!! if on input ind(1)  = 0 then indices are initialized in the routine,
!! if on input ind(1) != 0 then indices are assumed to have been
!!                         initialized before entering the routine and these
!!                         indices are carried around during the sorting process
!!
!! Adapted from Numerical Recipes pg. 329 (new edition)
!!
!----------------------------------------------------------------------------------------
MODULE env_sorting
    !------------------------------------------------------------------------------------
    !
    USE env_util_param, ONLY: DP
    !
    !------------------------------------------------------------------------------------
    !
    PRIVATE
    !
    PUBLIC :: env_hpsort, env_hpsort_eps
    !
    !------------------------------------------------------------------------------------
CONTAINS
    !------------------------------------------------------------------------------------
    !>
    !! Sort an array ra(1:n) into ascending order using heapsort algorithm.
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_hpsort(n, ra, ind)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER :: n
        INTEGER :: ind(n)
        REAL(DP) :: ra(n)
        INTEGER :: i, ir, j, l, iind
        REAL(DP) :: rra
        !
        !--------------------------------------------------------------------------------
        !
        IF (n < 1) RETURN
        !
        !--------------------------------------------------------------------------------
        ! initialize index array
        !
        IF (ind(1) == 0) THEN
            !
            DO i = 1, n
                ind(i) = i
            END DO
            !
        END IF
        !
        IF (n < 2) RETURN ! nothing to order
        !
        l = n / 2 + 1 ! initialize indices for hiring and retirement-promotion phase
        ir = n
        !
10      CONTINUE
        !
        !--------------------------------------------------------------------------------
        ! Still in hiring phase
        !
        IF (l > 1) THEN
            l = l - 1
            rra = ra(l)
            iind = ind(l)
        ELSE
            !
            !----------------------------------------------------------------------------
            ! In retirement-promotion phase.
            !
            rra = ra(ir) ! clear a space at the end of the array
            iind = ind(ir)
            ra(ir) = ra(1) ! retire the top of the heap into it
            ind(ir) = ind(1)
            ir = ir - 1 ! decrease the size of the corporation
            !
            !----------------------------------------------------------------------------
            ! Done with the last promotion
            !
            IF (ir == 1) THEN
                ra(1) = rra ! the least competent worker at all
                ind(1) = iind
                !
                RETURN
                !
            END IF
            !
        END IF
        !
        !--------------------------------------------------------------------------------
        ! Whether in hiring or promotion phase, we set up to place rra in its proper level
        !
        i = l
        j = l + l
        !
        DO WHILE (j <= ir)
            !
            IF (j < ir) THEN
                !
                !------------------------------------------------------------------------
                ! Compare to better underling
                !
                IF (ra(j) < ra(j + 1)) THEN
                    j = j + 1
                ELSE IF (ra(j) == ra(j + 1)) THEN
                    IF (ind(j) < ind(j + 1)) j = j + 1
                END IF
                !
            END IF
            !
            !----------------------------------------------------------------------------
            ! Demote rra
            !
            IF (rra < ra(j)) THEN
                ra(i) = ra(j)
                ind(i) = ind(j)
                i = j
                j = j + j
            ELSE IF (rra == ra(j)) THEN
                !
                !------------------------------------------------------------------------
                ! Demote rra
                !
                IF (iind < ind(j)) THEN
                    ra(i) = ra(j)
                    ind(i) = ind(j)
                    i = j
                    j = j + j
                ELSE
                    j = ir + 1 ! set j to terminate do-while loop
                END IF
                !
                ! this is the right place for rra
                !
            ELSE
                j = ir + 1 ! set j to terminate do-while loop
            END IF
            !
        END DO
        !
        ra(i) = rra
        ind(i) = iind
        !
        GOTO 10
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_hpsort
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
    !
    !------------------------------------------------------------------------------------
END MODULE env_sorting
!----------------------------------------------------------------------------------------
