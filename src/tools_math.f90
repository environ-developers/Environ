!----------------------------------------------------------------------------------------
!
! Copyright (C) 2018-2022 ENVIRON (www.quantum-environ.org)
! Copyright (C) 2002-2009 Quantum ESPRESSO (www.quantum-espresso.org)
!
!----------------------------------------------------------------------------------------
!
!     This file is part of Environ version 3.0
!
!     Environ 3.0 is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation, either version 2 of the License, or
!     (at your option) any later version.
!
!     Environ 3.0 is distributed in the hope that it will be useful,
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
MODULE tools_math
    !------------------------------------------------------------------------------------
    !
    USE env_util_param, ONLY: DP
    !
    !------------------------------------------------------------------------------------
    !
    PRIVATE
    !
    PUBLIC :: environ_erf, environ_erfc
    !
    !------------------------------------------------------------------------------------
CONTAINS
    !------------------------------------------------------------------------------------
    !>
    !! Error function - computed from the rational approximations of
    !! W. J. Cody, Math. Comp. 22 (1969), pages 631-637.
    !!
    !! for abs(x) le 0.47 erf is calculated directly
    !! for abs(x) gt 0.47 erf is calculated via erf(x)=1-erfc(x)
    !!
    !------------------------------------------------------------------------------------
    FUNCTION environ_erf(x)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP), INTENT(IN) :: x
        !
        REAL(DP) :: x2, p1(4), q1(4)
        REAL(DP) :: environ_erf
        !
        DATA p1/2.426679552305318E2_DP, 2.197926161829415E1_DP, &
            6.996383488619136_DP, -3.560984370181538E-2_DP/
        !
        DATA q1/2.150588758698612E2_DP, 9.116490540451490E1_DP, &
            1.508279763040779E1_DP, 1.000000000000000_DP/
        !
        !--------------------------------------------------------------------------------
        !
        IF (ABS(x) > 6.0_DP) THEN
            !
            environ_erf = SIGN(1.0_DP, x)
            ! erf(6) = 1 - 10^(-17) cannot be distinguished from 1
            !
        ELSE
            !
            IF (ABS(x) <= 0.47_DP) THEN
                x2 = x**2
                !
                environ_erf = x * &
                              (p1(1) + x2 * ( &
                               p1(2) + x2 * ( &
                               p1(3) + x2 * p1(4) &
                               ))) / &
                              (q1(1) + x2 * ( &
                               q1(2) + x2 * ( &
                               q1(3) + x2 * q1(4) &
                               )))
                !
            ELSE
                environ_erf = 1.0_DP - environ_erfc(x)
            END IF
            !
        END IF
        !
        !--------------------------------------------------------------------------------
    END FUNCTION environ_erf
    !------------------------------------------------------------------------------------
    !>
    !! erfc(x) = 1-erf(x) - See comments in erf
    !!
    !------------------------------------------------------------------------------------
    FUNCTION environ_erfc(x)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP), INTENT(IN) :: x
        !
        REAL(DP) :: environ_erfc
        REAL(DP) :: ax, x2, xm2, p2(8), q2(8), p3(5), q3(5), pim1
        !
        DATA p2/3.004592610201616E2_DP, 4.519189537118719E2_DP, &
            3.393208167343437E2_DP, 1.529892850469404E2_DP, &
            4.316222722205674E1_DP, 7.211758250883094_DP, &
            5.641955174789740E-1_DP, -1.368648573827167E-7_DP/
        !
        DATA q2/3.004592609569833E2_DP, 7.909509253278980E2_DP, &
            9.313540948506096E2_DP, 6.389802644656312E2_DP, &
            2.775854447439876E2_DP, 7.700015293522947E1_DP, &
            1.278272731962942E1_DP, 1.000000000000000_DP/
        !
        DATA p3/-2.996107077035422E-3_DP, -4.947309106232507E-2_DP, &
            -2.269565935396869E-1_DP, -2.786613086096478E-1_DP, &
            -2.231924597341847E-2_DP/
        !
        DATA q3/1.062092305284679E-2_DP, 1.913089261078298E-1_DP, &
            1.051675107067932_DP, 1.987332018171353_DP, &
            1.000000000000000_DP/
        !
        DATA pim1/0.56418958354775629_DP/ ! ( pim1= sqrt(1/pi) )
        !
        !--------------------------------------------------------------------------------
        !
        ax = ABS(x)
        !
        IF (ax > 26.0_DP) THEN
            environ_erfc = 0.0_DP ! erfc(26.0) = 10^(-296); erfc(9.0) = 10^(-37)
        ELSE IF (ax > 4.0_DP) THEN
            x2 = x**2
            xm2 = (1.0_DP / ax)**2
            !
            environ_erfc = 1.0_DP / ax * &
                           EXP(-x2) * &
                           (pim1 + xm2 * ( &
                            p3(1) + xm2 * ( &
                            p3(2) + xm2 * ( &
                            p3(3) + xm2 * ( &
                            p3(4) + xm2 * p3(5) &
                            )))) / &
                            (q3(1) + xm2 * ( &
                             q3(2) + xm2 * ( &
                             q3(3) + xm2 * ( &
                             q3(4) + xm2 * q3(5) &
                             )))))
            !
        ELSE IF (ax > 0.47_DP) THEN
            x2 = x**2
            !
            environ_erfc = EXP(-x2) * &
                           (p2(1) + ax * ( &
                            p2(2) + ax * ( &
                            p2(3) + ax * ( &
                            p2(4) + ax * ( &
                            p2(5) + ax * ( &
                            p2(6) + ax * ( &
                            p2(7) + ax * p2(8) &
                            ))))))) / &
                           (q2(1) + ax * ( &
                            q2(2) + ax * ( &
                            q2(3) + ax * ( &
                            q2(4) + ax * ( &
                            q2(5) + ax * ( &
                            q2(6) + ax * ( &
                            q2(7) + ax * q2(8) &
                            )))))))
            !
        ELSE
            environ_erfc = 1.0_DP - environ_erf(ax)
        END IF
        !
        IF (x < 0.0_DP) environ_erfc = 2.0_DP - environ_erfc
        ! erf(-x) = -erf(x) => erfc(-x) = 2 - erfc(x)
        !
        !--------------------------------------------------------------------------------
    END FUNCTION environ_erfc
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE tools_math
!----------------------------------------------------------------------------------------
