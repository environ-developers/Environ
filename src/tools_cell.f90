!----------------------------------------------------------------------------------------
!>
!!
!----------------------------------------------------------------------------------------
MODULE tools_cell
    !------------------------------------------------------------------------------------
    !
    USE modules_constants, ONLY: DP
    !
    USE cell_types, ONLY: environ_cell
    !
    !------------------------------------------------------------------------------------
CONTAINS
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
    !
    !------------------------------------------------------------------------------------
END MODULE tools_cell
!----------------------------------------------------------------------------------------
