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
    USE modules_constants, ONLY: DP
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
        IF (ind(1) .EQ. 0) THEN
            !
            DO i = 1, n
                ind(i) = i
            END DO
            !
        END IF
        !
        !--------------------------------------------------------------------------------
        !
        IF (n .LT. 2) RETURN ! nothing to order
        !
        l = n / 2 + 1
        ir = n
        !
        sorting: DO
            !
            IF (l .GT. 1) THEN ! hiring phase
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
                IF (ir .EQ. 1) THEN ! done with the last promotion
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
            DO WHILE (j .LE. ir)
                !
                IF (j .LT. ir) THEN
                    !
                    IF (ABS(ra(j) - ra(j + 1)) .GE. eps) THEN
                        !
                        IF (ra(j) .LT. ra(j + 1)) j = j + 1
                        ! compare to better underling
                        !
                    ELSE ! this means ra(j) == ra(j+1) within tolerance
                        !
                        IF (ind(j) .LT. ind(j + 1)) j = j + 1
                        !
                    END IF
                    !
                END IF
                !
                IF (ABS(rra - ra(j)) .GE. eps) THEN ! demote rra
                    !
                    IF (rra .LT. ra(j)) THEN
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
                    IF (iind .LT. ind(j)) THEN ! demote rra
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
    !>
    !! sort an array ra(1:n) into ascending order using heapsort algorithm.
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
    !! #TODO unused
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE hpsort(n, ra, ind)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER :: n
        INTEGER :: ind(*)
        REAL(DP) :: ra(*)
        !
        INTEGER :: i, ir, j, l, iind
        REAL(DP) :: rra
        !
        !--------------------------------------------------------------------------------
        ! initialize index array
        !
        IF (ind(1) .EQ. 0) THEN
            DO i = 1, n
                ind(i) = i
            END DO
        END IF
        !
        !--------------------------------------------------------------------------------
        !
        IF (n .LT. 2) RETURN ! nothing to order
        !
        !--------------------------------------------------------------------------------
        ! initialize indices for hiring and retirement-promotion phase
        !
        l = n / 2 + 1
        ir = n
10      CONTINUE
        !
        ! still in hiring phase
        !
        IF (l .GT. 1) THEN
            l = l - 1
            rra = ra(l)
            iind = ind(l)
            !
            ! in retirement-promotion phase.
            !
        ELSE
            !
            ! clear a space at the end of the array
            !
            rra = ra(ir)
            iind = ind(ir)
            !
            ! retire the top of the heap into it
            !
            ra(ir) = ra(1)
            ind(ir) = ind(1)
            !
            ! decrease the size of the corporation
            !
            ir = ir - 1
            !
            ! done with the last promotion
            !
            IF (ir .EQ. 1) THEN
                !
                ! the least competent worker at all
                !
                ra(1) = rra
                !
                ind(1) = iind
                !
                RETURN
                !
            END IF
            !
        END IF
        !
        ! wheter in hiring or promotion phase, we
        !
        i = l
        !
        ! set up to place rra in its proper level
        !
        j = l + l
        !
        DO WHILE (j .LE. ir)
            !
            IF (j .LT. ir) THEN
                !
                ! compare to better underling
                !
                IF (ra(j) .LT. ra(j + 1)) THEN
                    j = j + 1
                ELSE IF (ra(j) .EQ. ra(j + 1)) THEN
                    IF (ind(j) .LT. ind(j + 1)) j = j + 1
                END IF
                !
            END IF
            !
            ! demote rra
            !
            IF (rra .LT. ra(j)) THEN
                ra(i) = ra(j)
                ind(i) = ind(j)
                i = j
                j = j + j
            ELSE IF (rra .EQ. ra(j)) THEN
                !
                ! demote rra
                !
                IF (iind .LT. ind(j)) THEN
                    ra(i) = ra(j)
                    ind(i) = ind(j)
                    i = j
                    j = j + j
                ELSE
                    !
                    ! set j to terminate do-while loop
                    !
                    j = ir + 1
                END IF
                !
                ! this is the right place for rra
                !
            ELSE
                !
                ! set j to terminate do-while loop
                !
                j = ir + 1
            END IF
            !
        END DO
        !
        ra(i) = rra
        ind(i) = iind
        GOTO 10
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE hpsort
    !------------------------------------------------------------------------------------
    !>
    !! sort an integer array ia(1:n) into ascending order using heapsort algorithm.
    !! n is input, ia is replaced on output by its sorted rearrangement.
    !! create an index table (ind) by making an exchange in the index array
    !! whenever an exchange is made on the sorted data array (ia).
    !! in case of equal values in the data array (ia) the values in the
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
    !! #TODO unused
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE ihpsort(n, ia, ind)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER :: n
        INTEGER :: ind(*)
        INTEGER :: ia(*)
        !
        INTEGER :: i, ir, j, l, iind
        INTEGER :: iia
        !
        !--------------------------------------------------------------------------------
        ! initialize index array
        !
        IF (ind(1) .EQ. 0) THEN
            DO i = 1, n
                ind(i) = i
            END DO
        END IF
        !
        !--------------------------------------------------------------------------------
        !
        IF (n .LT. 2) RETURN ! nothing to order
        !
        !--------------------------------------------------------------------------------
        ! initialize indices for hiring and retirement-promotion phase
        !
        l = n / 2 + 1
        ir = n
10      CONTINUE
        !
        ! still in hiring phase
        !
        IF (l .GT. 1) THEN
            l = l - 1
            iia = ia(l)
            iind = ind(l)
            !
            ! in retirement-promotion phase.
            !
        ELSE
            !
            ! clear a space at the end of the array
            !
            iia = ia(ir)
            iind = ind(ir)
            !
            ! retire the top of the heap into it
            !
            ia(ir) = ia(1)
            ind(ir) = ind(1)
            !
            ir = ir - 1 ! decrease the size of the corporation
            !
            ! done with the last promotion
            !
            IF (ir .EQ. 1) THEN
                ! the least competent worker at all !
                ia(1) = iia
                !
                ind(1) = iind
                !
                RETURN
                !
            END IF
            !
        END IF
        !
        ! wheter in hiring or promotion phase, we
        !
        i = l
        !
        ! set up to place iia in its proper level
        !
        j = l + l
        !
        DO WHILE (j .LE. ir)
            !
            IF (j .LT. ir) THEN
                !
                ! compare to better underling
                !
                IF (ia(j) .LT. ia(j + 1)) THEN
                    j = j + 1
                ELSE IF (ia(j) .EQ. ia(j + 1)) THEN
                    IF (ind(j) .LT. ind(j + 1)) j = j + 1
                END IF
                !
            END IF
            !
            ! demote iia
            !
            IF (iia .LT. ia(j)) THEN
                ia(i) = ia(j)
                ind(i) = ind(j)
                i = j
                j = j + j
            ELSE IF (iia .EQ. ia(j)) THEN
                !
                ! demote iia
                !
                IF (iind .LT. ind(j)) THEN
                    ia(i) = ia(j)
                    ind(i) = ind(j)
                    i = j
                    j = j + j
                ELSE
                    !
                    ! set j to terminate do-while loop
                    !
                    j = ir + 1
                END IF
                !
                ! this is the right place for iia
                !
            ELSE
                !
                ! set j to terminate do-while loop
                !
                j = ir + 1
            END IF
            !
        END DO
        !
        ia(i) = iia
        ind(i) = iind
        GOTO 10
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE ihpsort
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE modules_sort
!----------------------------------------------------------------------------------------