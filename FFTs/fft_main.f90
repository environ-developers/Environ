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
MODULE env_fft_main
    !------------------------------------------------------------------------------------
    !
    USE env_fft_param
    !
    USE env_types_fft, ONLY: env_fft_type_descriptor
    !
#if defined(__CUDA)
    USE cudafor
    !
    USE env_fft_scalar, ONLY: env_cfft3d_gpu
    !
    USE env_fft_parallel_gpu, ONLY: env_tg_cft3s_gpu, env_many_cft3s_gpu, &
                                    env_tg_cft3s_2d_gpu, env_many_cft3s_2d_gpu
#else
    !
    USE env_fft_scalar, ONLY: env_cfft3d
    !
    USE env_fft_parallel, ONLY: env_tg_cft3s, env_many_cft3s, env_tg_cft3s_2d
    !
#endif
    !------------------------------------------------------------------------------------
    !
    INTERFACE env_invfft
        MODULE PROCEDURE env_invfft_y
#if defined(__CUDA)
        MODULE PROCEDURE env_invfft_y_gpu
#endif
    END INTERFACE env_invfft
    !
    INTERFACE env_fwfft
        MODULE PROCEDURE env_fwfft_y
#if defined(__CUDA)
        MODULE PROCEDURE env_fwfft_y_gpu
#endif
    END INTERFACE env_fwfft
    !
    !------------------------------------------------------------------------------------
    !
    PRIVATE
    !
    PUBLIC :: env_invfft, env_fwfft
    !
    !------------------------------------------------------------------------------------
CONTAINS
    !------------------------------------------------------------------------------------
    !>
    !! Compute G-space to R-space for a specific grid type
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_invfft_y(f, dfft, howmany)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, OPTIONAL, INTENT(IN) :: howmany
        !
        TYPE(env_fft_type_descriptor), INTENT(INOUT) :: dfft
        !
        COMPLEX(DP) :: f(:)
        INTEGER :: howmany_ = 1
        CHARACTER(LEN=12) :: clock_label
        !
        CHARACTER(LEN=80) :: sub_name = 'env_invfft'
        !
        !--------------------------------------------------------------------------------
        !
        IF (PRESENT(howmany)) THEN
            howmany_ = howmany
        ELSE
            howmany_ = 1
        END IF
        !
        clock_label = dfft%rho_clock_label
        !
        CALL env_start_clock(clock_label)
        !
        IF (dfft%lpara .AND. dfft%use_pencil_decomposition) THEN
            !
            IF (howmany_ == 1) THEN
                CALL env_tg_cft3s(f, dfft, 1)
            ELSE
                CALL env_many_cft3s(f, dfft, 1, howmany)
            END IF
            !
        ELSE IF (dfft%lpara) THEN
            !
            IF (howmany_ /= 1) THEN
                !
                CALL env_errore(sub_name, &
                                'howmany not yet implemented for parallel driver', 1)
                !
            END IF
            !
            CALL env_tg_cft3s_2d(f, dfft, 1)
            !
        ELSE
            !
            CALL env_cfft3d(f, dfft%nr1, dfft%nr2, dfft%nr3, &
                            dfft%nr1x, dfft%nr2x, dfft%nr3x, howmany_, 1)
            !
        END IF
        !
        CALL env_stop_clock(clock_label)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_invfft_y
    !------------------------------------------------------------------------------------
    !>
    !! Compute R-space to G-space for a specific grid type
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_fwfft_y(f, dfft, howmany)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, OPTIONAL, INTENT(IN) :: howmany
        !
        TYPE(env_fft_type_descriptor), INTENT(INOUT) :: dfft
        !
        COMPLEX(DP) :: f(:)
        INTEGER :: howmany_ = 1
        CHARACTER(LEN=12) :: clock_label
        !
        CHARACTER(LEN=80) :: sub_name = 'env_fwfft'
        !
        !--------------------------------------------------------------------------------
        !
        IF (PRESENT(howmany)) THEN
            howmany_ = howmany
        ELSE
            howmany_ = 1
        END IF
        !
        clock_label = dfft%rho_clock_label
        !
        CALL env_start_clock(clock_label)
        !
        IF (dfft%lpara .AND. dfft%use_pencil_decomposition) THEN
            !
            IF (howmany_ == 1) THEN
                CALL env_tg_cft3s(f, dfft, -1)
            ELSE
                CALL env_many_cft3s(f, dfft, -1, howmany)
            END IF
            !
        ELSE IF (dfft%lpara) THEN
            !
            IF (howmany_ /= 1) THEN
                CALL env_errore(sub_name, &
                                'howmany not yet implemented for parallel driver', 1)
            END IF
            !
            CALL env_tg_cft3s_2d(f, dfft, -1)
        ELSE
            !
            CALL env_cfft3d(f, dfft%nr1, dfft%nr2, dfft%nr3, &
                            dfft%nr1x, dfft%nr2x, dfft%nr3x, howmany_, -1)
            !
        END IF
        !
        CALL env_stop_clock(clock_label)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_fwfft_y
    !------------------------------------------------------------------------------------
#if defined(__CUDA)
    !>
    !! Compute G-space to R-space for a specific grid type
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_invfft_y_gpu(f_d, dfft, howmany, stream)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(env_fft_type_descriptor), INTENT(IN) :: dfft
        INTEGER, OPTIONAL, INTENT(IN) :: howmany
        INTEGER(kind=cuda_stream_kind), OPTIONAL, INTENT(IN) :: stream
        !
        COMPLEX(DP), DEVICE :: f_d(:)
        !
        INTEGER :: howmany_ = 1
        INTEGER(kind=cuda_stream_kind) :: stream_ = 0
        !
        CHARACTER(LEN=12) :: clock_label
        !
        CHARACTER(LEN=80) :: sub_name = 'env_invfft_y_gpu'
        !
        !--------------------------------------------------------------------------------
        !
        IF (PRESENT(howmany)) THEN
            howmany_ = howmany
        ELSE
            howmany_ = 1
        END IF
        !
        IF (PRESENT(stream)) THEN
            stream_ = stream
        ELSE
            stream_ = 0
        END IF
        !
        clock_label = dfft%rho_clock_label
        !
        CALL env_start_clock_gpu(clock_label)
        !
        IF (dfft%lpara .AND. dfft%use_pencil_decomposition) THEN
            !
            IF (stream_ /= 0) THEN
                !
                CALL env_errore(sub_name, &
                                'Stream support not implemented for parallel driver', 1)
                !
            END IF
            !
            IF (howmany_ == 1) THEN
                CALL env_tg_cft3s_gpu(f_d, dfft, 1)
            ELSE
                CALL env_many_cft3s_gpu(f_d, dfft, 1, howmany_)
            END IF
            !
        ELSE IF (dfft%lpara) THEN
            !
            IF (howmany_ /= 1) THEN
                CALL env_many_cft3s_2d_gpu(f_d, dfft, 1, howmany_)
            ELSE
                CALL env_tg_cft3s_2d_gpu(f_d, dfft, 1)
            END IF
            !
        ELSE
            !
            CALL env_cfft3d_gpu(f_d, dfft%nr1, dfft%nr2, dfft%nr3, &
                                dfft%nr1x, dfft%nr2x, dfft%nr3x, howmany_, 1, stream_)
            !
        END IF
        !
        CALL env_stop_clock_gpu(clock_label)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_invfft_y_gpu
    !------------------------------------------------------------------------------------
    !>
    !! Compute R-space to G-space for a specific grid type
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_fwfft_y_gpu(f_d, dfft, howmany, stream)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(env_fft_type_descriptor), INTENT(IN) :: dfft
        INTEGER, OPTIONAL, INTENT(IN) :: howmany
        INTEGER(kind=cuda_stream_kind), OPTIONAL, INTENT(IN) :: stream
        !
        COMPLEX(DP), DEVICE :: f_d(:)
        !
        INTEGER :: howmany_ = 1
        INTEGER(kind=cuda_stream_kind) :: stream_ = 0
        !
        CHARACTER(LEN=12) :: clock_label
        !
        CHARACTER(LEN=80) :: sub_name = 'env_fwfft_y_gpu'
        !
        !--------------------------------------------------------------------------------
        !
        IF (PRESENT(howmany)) THEN
            howmany_ = howmany
        ELSE
            howmany_ = 1
        END IF
        !
        IF (PRESENT(stream)) THEN
            stream_ = stream
        ELSE
            stream_ = 0
        END IF
        !
        clock_label = dfft%rho_clock_label
        !
        CALL env_start_clock_gpu(clock_label)
        !
        IF (dfft%lpara .AND. dfft%use_pencil_decomposition) THEN
            !
            IF (stream_ /= 0) THEN
                !
                CALL env_errore(sub_name, &
                                'Stream support not implemented for parallel driver', 1)
                !
            END IF
            !
            IF (howmany_ == 1) THEN
                CALL env_tg_cft3s_gpu(f_d, dfft, -1)
            ELSE
                CALL env_many_cft3s_gpu(f_d, dfft, -1, howmany_)
            END IF
            !
        ELSE IF (dfft%lpara) THEN
            !
            IF (howmany_ /= 1) THEN
                CALL env_many_cft3s_2d_gpu(f_d, dfft, -1, howmany_)
            ELSE
                CALL env_tg_cft3s_2d_gpu(f_d, dfft, -1)
            END IF
            !
        ELSE
            !
            CALL env_cfft3d_gpu(f_d, dfft%nr1, dfft%nr2, dfft%nr3, &
                                dfft%nr1x, dfft%nr2x, dfft%nr3x, howmany_, -1, stream_)
            !
        END IF
        !
        CALL env_stop_clock_gpu(clock_label)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_fwfft_y_gpu
    !------------------------------------------------------------------------------------
#endif
    !
    !------------------------------------------------------------------------------------
END MODULE env_fft_main
!----------------------------------------------------------------------------------------
