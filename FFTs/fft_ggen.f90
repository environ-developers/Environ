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
!! Subroutines generating variables nl* needed to map G-vector
!! components onto the FFT grid(s) in reciprocal space
!!
!----------------------------------------------------------------------------------------
MODULE env_fft_ggen
    !------------------------------------------------------------------------------------
    !
    USE env_fft_param
    !
    USE env_types_fft, ONLY: env_fft_type_descriptor
    !
    !------------------------------------------------------------------------------------
    !
    PRIVATE
    !
    PUBLIC :: env_fft_set_nl
    !
    !------------------------------------------------------------------------------------
CONTAINS
    !------------------------------------------------------------------------------------
    !>
    !! Input:  FFT descriptor dfft, lattice vectors at, list of G-vectors g
    !! Output: indices nl such that G_fft(nl(i)) = G(i)
    !!         indices nlm such that G_fft(nlm(i)) = -G(i) only if lgamma=.true.
    !!         optionally, Miller indices: if bg = reciprocal lattice vectors,
    !! G(:,i) = mill(1,i)*bg(:,1) + mill(2,i)*bg(:,2) + mill(3,i)*bg(:,3)
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_fft_set_nl(dfft, at, g, mill)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP), INTENT(IN) :: g(:, :)
        REAL(DP), INTENT(IN) :: at(:, :)
        !
        TYPE(env_fft_type_descriptor), INTENT(INOUT) :: dfft
        INTEGER, OPTIONAL, INTENT(OUT) :: mill(:, :)
        !
        INTEGER :: ng, n1, n2, n3
        !
        CHARACTER(LEN=80) :: sub_name = 'env_fft_set_nl'
        !
        !--------------------------------------------------------------------------------
        !
        IF (ALLOCATED(dfft%nl)) DEALLOCATE (dfft%nl)
        !
        ALLOCATE (dfft%nl(dfft%ngm))
        !
        IF (dfft%lgamma) THEN
            !
            IF (ALLOCATED(dfft%nlm)) DEALLOCATE (dfft%nlm)
            !
            ALLOCATE (dfft%nlm(dfft%ngm))
        END IF
        !
        DO ng = 1, dfft%ngm
            n1 = NINT(SUM(g(:, ng) * at(:, 1)))
            !
            IF (PRESENT(mill)) mill(1, ng) = n1
            !
            IF (n1 < 0) n1 = n1 + dfft%nr1
            !
            n2 = NINT(SUM(g(:, ng) * at(:, 2)))
            !
            IF (PRESENT(mill)) mill(2, ng) = n2
            !
            IF (n2 < 0) n2 = n2 + dfft%nr2
            !
            n3 = NINT(SUM(g(:, ng) * at(:, 3)))
            !
            IF (PRESENT(mill)) mill(3, ng) = n3
            !
            IF (n3 < 0) n3 = n3 + dfft%nr3
            !
            IF (n1 >= dfft%nr1 .OR. n2 >= dfft%nr2 .OR. n3 >= dfft%nr3) &
                CALL env_fft_error(sub_name, 'Mesh too small?', ng)
            !
            IF (dfft%lpara) THEN
                !
                dfft%nl(ng) = 1 + n3 + &
                              (dfft%isind(1 + n1 + n2 * dfft%nr1x) - 1) * dfft%nr3x
                !
            ELSE
                dfft%nl(ng) = 1 + n1 + n2 * dfft%nr1x + n3 * dfft%nr1x * dfft%nr2x
            END IF
            !
            IF (dfft%lgamma) THEN
                n1 = -n1
                !
                IF (n1 < 0) n1 = n1 + dfft%nr1
                !
                n2 = -n2
                !
                IF (n2 < 0) n2 = n2 + dfft%nr2
                !
                n3 = -n3
                !
                IF (n3 < 0) n3 = n3 + dfft%nr3
                !
                IF (dfft%lpara) THEN
                    !
                    dfft%nlm(ng) = 1 + n3 + &
                                   (dfft%isind(1 + n1 + n2 * dfft%nr1x) - 1) * dfft%nr3x
                    !
                ELSE
                    dfft%nlm(ng) = 1 + n1 + n2 * dfft%nr1x + n3 * dfft%nr1x * dfft%nr2x
                END IF
                !
            END IF
            !
        END DO
        !
#if defined(__CUDA)
        IF (ALLOCATED(dfft%nl_d)) DEALLOCATE (dfft%nl_d)
        !
        ALLOCATE (dfft%nl_d, SOURCE=dfft%nl)
        !
        IF (dfft%lgamma) THEN
            !
            IF (ALLOCATED(dfft%nlm_d)) DEALLOCATE (dfft%nlm_d)
            !
            ALLOCATE (dfft%nlm_d, SOURCE=dfft%nlm)
        END IF
        !
#endif
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_fft_set_nl
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE env_fft_ggen
!----------------------------------------------------------------------------------------
