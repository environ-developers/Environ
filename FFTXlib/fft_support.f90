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
! Authors
!
!----------------------------------------------------------------------------------------
!>
!! FFT support Functions/Subroutines
!!
!----------------------------------------------------------------------------------------
MODULE env_fft_support
    !------------------------------------------------------------------------------------
    !
    USE env_fft_param
    !
    !------------------------------------------------------------------------------------
    !
    PRIVATE
    !
    PUBLIC :: env_good_fft_dimension, env_good_fft_order
    !
    !------------------------------------------------------------------------------------
CONTAINS
    !------------------------------------------------------------------------------------
    !>
    !! Determines the optimal maximum dimensions of fft arrays
    !! Useful on some machines to avoid memory conflicts
    !!
    !------------------------------------------------------------------------------------
    INTEGER FUNCTION env_good_fft_dimension(n)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER :: n, nx
        REAL(DP) :: log2n
        !
        !--------------------------------------------------------------------------------
        !
        nx = n ! this is the default: max dimension = fft dimension
        !
#if defined(__LINUX_ESSL)
        log2n = LOG(DBLE(n)) / LOG(2.0_DP) ! log2n is the logarithm of n in base 2
        !
        IF (ABS(NINT(log2n) - log2n) < 1.0D-8) nx = n + 1
        ! if n is a power of 2 (log2n is integer) increase dimension by 1
        !
#elif defined(__SX6)
        !
        IF (MOD(n, 2) == 0) nx = n + 1
        ! for nec vector machines: if n is even increase dimension by 1
        !
#endif
        !
        env_good_fft_dimension = nx
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END FUNCTION env_good_fft_dimension
    !------------------------------------------------------------------------------------
    !>
    !! This function find a "good" fft order value greater or equal to "nr"
    !!
    !! nr  (input) tentative order n of a fft
    !!
    !! np  (optional input) if present restrict the search of the order
    !!     in the ensemble of multiples of np
    !!
    !! Output: the same if n is a good number
    !!         the closest higher number that is good
    !!         an fft order is not good if not implemented (as on IBM with ESSL)
    !!         or implemented but with awful performances (most other cases)
    !!
    !------------------------------------------------------------------------------------
    INTEGER FUNCTION env_good_fft_order(nr, np)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: nr
        INTEGER, OPTIONAL, INTENT(IN) :: np
        !
        INTEGER :: new
        !
        CHARACTER(LEN=80) :: fun_name = 'env_good_fft_order'
        !
        !--------------------------------------------------------------------------------
        !
        new = nr
        !
        IF (PRESENT(np)) THEN
            !
            IF (np <= 0 .OR. np > nr) &
                CALL env_fft_error(fun_name, ' invalid np ', 1)
            !
            DO WHILE (((.NOT. env_allowed(new)) .OR. (MOD(new, np) /= 0)) .AND. &
                      (new <= nfftx))
                new = new + 1
            END DO
            !
        ELSE
            !
            DO WHILE ((.NOT. env_allowed(new)) .AND. (new <= nfftx))
                new = new + 1
            END DO
            !
        END IF
        !
        IF (new > nfftx) &
            CALL env_fft_error(fun_name, ' fft order too large ', new)
        !
        env_good_fft_order = new
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END FUNCTION env_good_fft_order
    !------------------------------------------------------------------------------------
    !>
    !! Find if the fft dimension is a good one
    !! a "bad one" is either not implemented (as on IBM with ESSL)
    !! or implemented but with awful performances (most other cases)
    !!
    !------------------------------------------------------------------------------------
    LOGICAL FUNCTION env_allowed(nr)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER :: nr
        INTEGER :: pwr(5)
        INTEGER :: mr, i, fac, p, maxpwr
        INTEGER :: factors(5) = (/2, 3, 5, 7, 11/)
        !
        CHARACTER(LEN=80) :: fun_name = 'allowed'
        !
        !--------------------------------------------------------------------------------
        ! Find the factors of the fft dimension
        !
        mr = nr
        pwr = 0
        factors_loop: DO i = 1, 5
            fac = factors(i)
            maxpwr = NINT(LOG(DBLE(mr)) / LOG(DBLE(fac))) + 1
            !
            DO p = 1, maxpwr
                !
                IF (mr == 1) EXIT factors_loop
                !
                IF (MOD(mr, fac) == 0) THEN
                    mr = mr / fac
                    pwr(i) = pwr(i) + 1
                END IF
                !
            END DO
            !
        END DO factors_loop
        !
        IF (nr /= (mr * 2**pwr(1) * 3**pwr(2) * 5**pwr(3) * 7**pwr(4) * 11**pwr(5))) &
            CALL env_fft_error(fun_name, ' what ?!? ', 1)
        !
        IF (mr /= 1) THEN
            !
            env_allowed = .FALSE.
            ! fft dimension contains factors > 11 : no good in any case
            !
        ELSE
            !
#if defined(__LINUX_ESSL)
            !
            !----------------------------------------------------------------------------
            ! IBM machines with essl libraries
            !
            env_allowed = (pwr(1) >= 1) .AND. (pwr(2) <= 2) .AND. (pwr(3) <= 1) .AND. &
                          (pwr(4) <= 1) .AND. (pwr(5) <= 1) .AND. &
                          (((pwr(2) == 0) .AND. (pwr(3) + pwr(4) + pwr(5)) <= 2) .OR. &
                           ((pwr(2) /= 0) .AND. (pwr(3) + pwr(4) + pwr(5)) <= 1))
            !
#else
            !
            !----------------------------------------------------------------------------
            ! fftw and all other cases: no factors 7 and 11
            !
            env_allowed = ((pwr(4) == 0) .AND. (pwr(5) == 0))
            !
#endif
            !
        END IF
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END FUNCTION env_allowed
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE env_fft_support
!----------------------------------------------------------------------------------------
