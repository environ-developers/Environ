!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------------------
!>
!!
!----------------------------------------------------------------------------------------
MODULE modules_sort
    !------------------------------------------------------------------------------------
    !
    PRIVATE
    !
    PUBLIC :: hpsort_eps, hpsort, ihpsort
    !
    !------------------------------------------------------------------------------------
CONTAINS
    !------------------------------------------------------------------------------------
    !>
    !! sort an array ra(1:n) into ascending order using heapsort algorithm,
    !! and considering two elements being equal if their values differ
    !! for less than "eps".
    !! n is input, ra is replaced on output by its sorted rearrangement.
    !! create an index table (ind) by making an exchange in the index array
    !! whenever an exchange is made on the sorted data array (ra).
    !! in case of equal values in the data array (ra) the values in the
    !! index array (ind) are used to order the entries.
    !! if on input ind(1)  = 0 then indices are initialized in the routine,
    !! if on input ind(1) != 0 then indices are assumed to have been
    !!                initialized before entering the routine and these
    !!                indices are carried around during the sorting process
    !!
    !! no work space needed !
    !! free us from machine-dependent sorting-routines !
    !!
    !! adapted from Numerical Recipes pg. 329 (new edition)
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE hpsort_eps(n, ra, ind, eps)
        !--------------------------------------------------------------------------------
        !
        USE modules_constants, ONLY: DP
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
    END SUBROUTINE hpsort_eps
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE modules_sort
!----------------------------------------------------------------------------------------