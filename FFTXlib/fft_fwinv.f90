!
! Copyright (C) Quantum ESPRESSO group
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!=---------------------------------------------------------------------------=!
SUBROUTINE env_invfft_y( fft_kind, f, dfft, howmany )
  !! Compute G-space to R-space for a specific grid type
  !! 
  !! **fft_kind = 'Rho'** : 
  !!   inverse (backward) fourier transform of potentials and charge density f
  !!   On output, f is overwritten
  !! 
  !! **fft_kind = 'Wave'** :
  !!   inverse (backward) fourier transform of  wave functions f
  !!   On output, f is overwritten
  !!
  !! **fft_kind = 'tgWave'** :
  !!   inverse (backward) fourier transform of  wave functions f with task group
  !!   On output, f is overwritten
  !!
  !! **dfft = FFT grid descriptor**, IMPORTANT NOTICE: grid is specified only by dfft.
  !!   No check is performed on the correspondence between dfft and fft_kind. 
  !!   from all other cases
  
  USE env_fft_scalar,    ONLY: env_cfft3d, env_cfft3ds
  USE env_fft_parallel,  ONLY: env_tg_cft3s, env_many_cft3s
  USE env_fft_parallel_2d,  ONLY: tg_cft3s_2d => env_tg_cft3s
  USE env_fft_types,     ONLY: env_fft_type_descriptor
  USE env_fft_param,     ONLY: DP

  IMPLICIT NONE

  TYPE(env_fft_type_descriptor), INTENT(INOUT) :: dfft
  CHARACTER(LEN=*), INTENT(IN) :: fft_kind
  COMPLEX(DP) :: f(:)
  INTEGER, OPTIONAL, INTENT(IN) :: howmany
  INTEGER :: howmany_ = 1
  CHARACTER(LEN=12) :: clock_label

  IF(PRESENT(howmany) ) THEN
     howmany_ = howmany
  ELSE
     howmany_ = 1
  END IF
  !
  IF( fft_kind == 'Rho' ) THEN
     clock_label = dfft%rho_clock_label
  ELSE IF( fft_kind == 'Wave' .OR. fft_kind == 'tgWave' ) THEN
     clock_label = dfft%wave_clock_label
  ELSE
     CALL env_fftx_error__( ' env_invfft ', ' unknown fft kind : '//fft_kind , 1 )
  END IF
  IF (clock_label == ' ') CALL env_fftx_error__( ' env_invfft ', ' uninitialized fft kind : '//fft_kind , 1 )

  CALL start_clock(clock_label)

  IF( dfft%lpara .and. dfft%use_pencil_decomposition ) THEN

     IF( howmany_ == 1 ) THEN
     IF( fft_kind == 'Rho' ) THEN
        CALL env_tg_cft3s( f, dfft, 1 )
     ELSE IF( fft_kind == 'Wave' ) THEN
        CALL env_tg_cft3s( f, dfft, 2 )
     ELSE IF( fft_kind == 'tgWave' ) THEN
        CALL env_tg_cft3s( f, dfft, 3 )
        END IF
     ELSE
        IF( fft_kind == 'Rho' ) THEN
           CALL env_many_cft3s( f, dfft, 1, howmany )
        ELSE IF( fft_kind == 'Wave' ) THEN
           CALL env_many_cft3s( f, dfft, 2, howmany )
        ELSE IF( fft_kind == 'tgWave' ) THEN
           CALL env_fftx_error__( ' env_invfft ', ' howmany not yet implemented for parallel driver ', 1 )
        END IF
     END IF

  ELSE IF( dfft%lpara ) THEN

     IF( howmany_ /= 1 ) THEN
        CALL env_fftx_error__( ' env_invfft ', ' howmany not yet implemented for parallel driver ', 1 )
     END IF

     IF( fft_kind == 'Rho' ) THEN
        CALL tg_cft3s_2d(f,dfft, 1)
     ELSE IF( fft_kind == 'Wave' ) THEN
        CALL tg_cft3s_2d(f,dfft, 2 )
     ELSE IF( fft_kind == 'tgWave' ) THEN
        CALL env_fftx_error__( ' env_fwfft ', ' tgWave not implemented  ', 1 )
     END IF

  ELSE

     IF( fft_kind == 'Rho' ) THEN
        CALL env_cfft3d( f, dfft%nr1, dfft%nr2, dfft%nr3, &
                        dfft%nr1x, dfft%nr2x, dfft%nr3x, howmany_ , 1)
     ELSE 
        CALL env_cfft3ds( f, dfft%nr1, dfft%nr2, dfft%nr3, &
                        dfft%nr1x,dfft%nr2x,dfft%nr3x, howmany_ , 1, &
                        dfft%isind, dfft%iplw )
     END IF

  END IF

  CALL stop_clock( clock_label )

  RETURN

END SUBROUTINE env_invfft_y
!
!=---------------------------------------------------------------------------=!
!
SUBROUTINE env_fwfft_y( fft_kind, f, dfft, howmany )
  !! Compute R-space to G-space for a specific grid type
  !! 
  !! **fft_kind = 'Rho'**
  !!   forward fourier transform of potentials and charge density f
  !!   On output, f is overwritten
  !! 
  !! **fft_kind = 'Wave'**
  !!   forward fourier transform of  wave functions f
  !!   On output, f is overwritten
  !!
  !! **fft_kind = 'tgWave'**
  !!   forward fourier transform of wave functions f with task group
  !!   On output, f is overwritten
  !! 
  
  USE env_fft_scalar,    ONLY: env_cfft3d, env_cfft3ds
  USE env_fft_parallel,  ONLY: env_tg_cft3s, env_many_cft3s
  USE env_fft_parallel_2d,  ONLY: tg_cft3s_2d => env_tg_cft3s
  USE env_fft_types,     ONLY: env_fft_type_descriptor
  USE env_fft_param,     ONLY: DP

  IMPLICIT NONE

  TYPE(env_fft_type_descriptor), INTENT(INOUT) :: dfft
  CHARACTER(LEN=*), INTENT(IN) :: fft_kind
  COMPLEX(DP) :: f(:)
  INTEGER, OPTIONAL, INTENT(IN) :: howmany
  INTEGER :: howmany_ = 1
  CHARACTER(LEN=12) :: clock_label

  IF(PRESENT(howmany) ) THEN
     howmany_ = howmany
  ELSE
     howmany_ = 1
  END IF

  IF( fft_kind == 'Rho' ) THEN
     clock_label = dfft%rho_clock_label
  ELSE IF( fft_kind == 'Wave' .OR. fft_kind == 'tgWave' ) THEN
     clock_label = dfft%wave_clock_label
  ELSE
     CALL env_fftx_error__( ' env_fwfft ', ' unknown fft kind: '//fft_kind , 1 )
  END IF
  IF (clock_label == ' ') CALL env_fftx_error__( ' env_fwfft ', ' uninitialized fft kind : '//fft_kind , 1 )

  CALL start_clock(clock_label)

  IF( dfft%lpara .and. dfft%use_pencil_decomposition ) THEN
     
     IF( howmany_ == 1 ) THEN
     IF( fft_kind == 'Rho' ) THEN
        CALL env_tg_cft3s(f,dfft,-1)
     ELSE IF( fft_kind == 'Wave' ) THEN
        CALL env_tg_cft3s(f,dfft,-2 )
     ELSE IF( fft_kind == 'tgWave' ) THEN
        CALL env_tg_cft3s(f,dfft,-3 )
        END IF
     ELSE
        IF( fft_kind == 'Rho' ) THEN
           CALL env_many_cft3s( f, dfft, -1, howmany )
        ELSE IF( fft_kind == 'Wave' ) THEN
           CALL env_many_cft3s( f, dfft, -2, howmany )
        ELSE IF( fft_kind == 'tgWave' ) THEN
           CALL env_many_cft3s( f, dfft, -3, howmany )
        END IF
     END IF

  ELSE IF( dfft%lpara ) THEN

     IF( howmany_ /= 1 ) THEN
        CALL env_fftx_error__( ' env_fwfft ', ' howmany not yet implemented for parallel driver ', 1 )
     END IF

     IF( fft_kind == 'Rho' ) THEN
        CALL tg_cft3s_2d(f,dfft,-1)
     ELSE IF( fft_kind == 'Wave' ) THEN
        CALL tg_cft3s_2d(f,dfft,-2 )
     ELSE IF( fft_kind == 'tgWave' ) THEN
        CALL env_fftx_error__( ' env_fwfft ', ' tgWave not implemented  ', 1 )
     END IF

  ELSE

     IF( fft_kind == 'Rho' ) THEN
        CALL env_cfft3d( f, dfft%nr1, dfft%nr2, dfft%nr3, &
                        dfft%nr1x,dfft%nr2x,dfft%nr3x, howmany_ , -1)
     ELSE 
        CALL env_cfft3ds( f, dfft%nr1, dfft%nr2, dfft%nr3, &
                         dfft%nr1x,dfft%nr2x,dfft%nr3x, howmany_ , -1, &
                         dfft%isind, dfft%iplw )
     END IF

  END IF

  CALL stop_clock( clock_label )
  
  RETURN
  !
END SUBROUTINE env_fwfft_y
!=---------------------------------------------------------------------------=!

!
!=---------------------------------------------------------------------------=!
!=---------------------------------------------------------------------------=!

#if defined(__CUDA)

SUBROUTINE env_invfft_y_gpu( fft_kind, f_d, dfft, howmany, stream )
  !! Compute G-space to R-space for a specific grid type
  !! 
  !! **fft_kind = 'Rho'** : 
  !!   inverse (backward) fourier transform of potentials and charge density f
  !!   On output, f is overwritten
  !! 
  !! **fft_kind = 'Wave'** :
  !!   inverse (backward) fourier transform of  wave functions f
  !!   On output, f is overwritten
  !!
  !! **fft_kind = 'tgWave'** :
  !!   inverse (backward) fourier transform of  wave functions f with task group
  !!   On output, f is overwritten
  !!
  !! **dfft = FFT grid descriptor**, IMPORTANT NOTICE: grid is specified only by dfft.
  !!   No check is performed on the correspondence between dfft and fft_kind. 
  !!   from all other cases
  USE cudafor
  USE env_fft_scalar,    ONLY: env_cfft3d_gpu, env_cfft3ds_gpu
  USE env_fft_parallel,  ONLY: env_tg_cft3s_gpu, env_many_cft3s_gpu
  USE env_fft_parallel_2d,  ONLY: tg_cft3s_2d_gpu => env_tg_cft3s_gpu, &
                            & many_cft3s_2d_gpu => env_many_cft3s_gpu
  USE env_fft_types,     ONLY: env_fft_type_descriptor
  USE env_fft_param,     ONLY: DP

  IMPLICIT NONE

  TYPE(env_fft_type_descriptor), INTENT(IN) :: dfft
  CHARACTER(LEN=*), INTENT(IN) :: fft_kind
  COMPLEX(DP), DEVICE :: f_d(:)
  INTEGER, OPTIONAL, INTENT(IN) :: howmany
  INTEGER(kind = cuda_stream_kind), OPTIONAL, INTENT(IN) :: stream
  !
  INTEGER                          :: howmany_ = 1
  INTEGER(kind = cuda_stream_kind) :: stream_  = 0

  CHARACTER(LEN=12) :: clock_label

  IF(PRESENT(howmany) ) THEN
     howmany_ = howmany
  ELSE
     howmany_ = 1
  END IF
  !
  IF( present( stream ) ) THEN
    stream_ = stream
  ELSE
    stream_ = 0
  END IF
  !
  IF( fft_kind == 'Rho' ) THEN
     clock_label = dfft%rho_clock_label
  ELSE IF( fft_kind == 'Wave' .OR. fft_kind == 'tgWave' ) THEN
     clock_label = dfft%wave_clock_label
  ELSE
     CALL env_fftx_error__( ' env_invfft ', ' unknown fft kind : '//fft_kind , 1 )
  END IF
  IF (clock_label == ' ') CALL env_fftx_error__( ' env_invfft ', ' uninitialized fft kind : '//fft_kind , 1 )

  CALL start_clock_gpu(clock_label)

  IF( dfft%lpara .and. dfft%use_pencil_decomposition ) THEN

     !IF( howmany_ /= 1 .and. fft_kind /= 'Wave' ) THEN
     !   CALL env_fftx_error__( ' env_invfft ', ' howmany not yet implemented for parallel driver ', 1 )
     !END IF
     IF( stream_ /= 0 ) THEN
        CALL env_fftx_error__( ' env_invfft ', ' stream support not implemented for parallel driver ', 1 )
     END IF
     
     IF( fft_kind == 'Rho' ) THEN
        IF( howmany_ == 1 ) THEN
           CALL env_tg_cft3s_gpu( f_d, dfft, 1 )
        ELSE
            CALL env_many_cft3s_gpu( f_d, dfft, 1, howmany_ )
        END IF
     ELSE IF( fft_kind == 'Wave' ) THEN
        IF( howmany_ == 1 ) THEN
            CALL env_tg_cft3s_gpu( f_d, dfft, 2 )
        ELSE
            CALL env_many_cft3s_gpu( f_d, dfft, 2, howmany_ )
        END IF
     ELSE IF( fft_kind == 'tgWave' ) THEN
        CALL env_tg_cft3s_gpu( f_d, dfft, 3 )
     END IF

  ELSE IF( dfft%lpara ) THEN

     IF( howmany_ /= 1 ) THEN
        IF( fft_kind == 'Rho' ) THEN
           CALL many_cft3s_2d_gpu( f_d, dfft, 1,  howmany_)
        ELSE IF( fft_kind == 'Wave' ) THEN
           CALL many_cft3s_2d_gpu( f_d, dfft, 2, howmany_ )
        END IF
     ELSE
        IF( fft_kind == 'Rho' ) THEN
           CALL tg_cft3s_2d_gpu( f_d, dfft, 1 )
        ELSE IF( fft_kind == 'Wave' ) THEN
           CALL tg_cft3s_2d_gpu( f_d, dfft, 2 )
        END IF
     END IF

  ELSE

     IF( fft_kind == 'Rho' ) THEN
        CALL env_cfft3d_gpu( f_d, dfft%nr1, dfft%nr2, dfft%nr3, &
                        dfft%nr1x, dfft%nr2x, dfft%nr3x, howmany_ , 1, stream_)
     ELSE 
        CALL env_cfft3ds_gpu( f_d, dfft%nr1, dfft%nr2, dfft%nr3, &
                        dfft%nr1x,dfft%nr2x,dfft%nr3x, howmany_ , 1, &
                        dfft%isind, dfft%iplw, stream_ )
     END IF

  END IF

  CALL stop_clock_gpu( clock_label )

  RETURN

END SUBROUTINE env_invfft_y_gpu
!
!=---------------------------------------------------------------------------=!
!
SUBROUTINE env_fwfft_y_gpu( fft_kind, f_d, dfft, howmany, stream )
  !! Compute R-space to G-space for a specific grid type
  !! 
  !! **fft_kind = 'Rho'**
  !!   forward fourier transform of potentials and charge density f
  !!   On output, f is overwritten
  !! 
  !! **fft_kind = 'Wave'**
  !!   forward fourier transform of  wave functions f
  !!   On output, f is overwritten
  !!
  !! **fft_kind = 'tgWave'**
  !!   forward fourier transform of wave functions f with task group
  !!   On output, f is overwritten
  !! 
  USE cudafor
  USE env_fft_scalar,    ONLY: env_cfft3d_gpu, env_cfft3ds_gpu
  USE env_fft_parallel,  ONLY: env_tg_cft3s_gpu, env_many_cft3s_gpu
  USE env_fft_parallel_2d,  ONLY: tg_cft3s_2d_gpu => env_tg_cft3s_gpu, &
                              & many_cft3s_2d_gpu => env_many_cft3s_gpu
  USE env_fft_types,     ONLY: env_fft_type_descriptor
  USE env_fft_param,     ONLY: DP

  IMPLICIT NONE

  TYPE(env_fft_type_descriptor), INTENT(IN) :: dfft
  CHARACTER(LEN=*), INTENT(IN) :: fft_kind
  COMPLEX(DP), DEVICE :: f_d(:)
  INTEGER, OPTIONAL, INTENT(IN) :: howmany
  INTEGER(kind = cuda_stream_kind), OPTIONAL, INTENT(IN) :: stream

  INTEGER                          :: howmany_ = 1
  INTEGER(kind = cuda_stream_kind) :: stream_  = 0

  CHARACTER(LEN=12) :: clock_label
  !
  IF(PRESENT(howmany) ) THEN
     howmany_ = howmany
  ELSE
     howmany_ = 1
  END IF
  !
  IF( present( stream ) ) THEN
    stream_ = stream
  ELSE
    stream_ = 0
  END IF
  !
  IF( fft_kind == 'Rho' ) THEN
     clock_label = dfft%rho_clock_label
  ELSE IF( fft_kind == 'Wave' .OR. fft_kind == 'tgWave' ) THEN
     clock_label = dfft%wave_clock_label
  ELSE
     CALL env_fftx_error__( ' env_fwfft ', ' unknown fft kind: '//fft_kind , 1 )
  END IF
  IF (clock_label == ' ') CALL env_fftx_error__( ' env_fwfft ', ' uninitialized fft kind : '//fft_kind , 1 )

  CALL start_clock_gpu(clock_label)

  IF( dfft%lpara .and. dfft%use_pencil_decomposition ) THEN

     !IF( howmany_ /= 1 .and. fft_kind /= 'Wave' ) THEN
     !   CALL env_fftx_error__( ' env_fwfft ', ' howmany not yet implemented for parallel driver ', 1 )
     !END IF
     IF( stream_ /= 0 ) THEN
        CALL env_fftx_error__( ' env_fwfft ', ' stream support not implemented for parallel driver ', 1 )
     END IF
     
     IF( fft_kind == 'Rho' ) THEN
        IF( howmany_ == 1 ) THEN
           CALL env_tg_cft3s_gpu(f_d,dfft,-1)
        ELSE
           CALL env_many_cft3s_gpu(f_d,dfft,-1, howmany_)
        END IF
     ELSE IF( fft_kind == 'Wave' ) THEN
        IF( howmany_ == 1 ) THEN
           CALL env_tg_cft3s_gpu(f_d,dfft,-2)
        ELSE
           CALL env_many_cft3s_gpu(f_d,dfft,-2, howmany_)
        END IF
     ELSE IF( fft_kind == 'tgWave' ) THEN
        CALL env_tg_cft3s_gpu( f_d, dfft, -3 )
     END IF

  ELSE IF( dfft%lpara ) THEN

     IF( howmany_ /= 1 ) THEN
        IF( fft_kind == 'Rho' ) THEN
           CALL many_cft3s_2d_gpu( f_d, dfft, -1, howmany_)
        ELSE IF( fft_kind == 'Wave' ) THEN
           CALL many_cft3s_2d_gpu( f_d, dfft, -2, howmany_ )
        END IF
     ELSE
        IF( fft_kind == 'Rho' ) THEN
           CALL tg_cft3s_2d_gpu( f_d, dfft, -1 )
        ELSE IF( fft_kind == 'Wave' ) THEN
           CALL tg_cft3s_2d_gpu( f_d, dfft, -2 )
        END IF
     END IF

  ELSE

     IF( fft_kind == 'Rho' ) THEN
        CALL env_cfft3d_gpu( f_d, dfft%nr1, dfft%nr2, dfft%nr3, &
                        dfft%nr1x,dfft%nr2x,dfft%nr3x, howmany_ , -1, stream_ )
     ELSE 
        CALL env_cfft3ds_gpu( f_d, dfft%nr1, dfft%nr2, dfft%nr3, &
                         dfft%nr1x,dfft%nr2x,dfft%nr3x, howmany_ , -1, &
                         dfft%isind, dfft%iplw, stream_ )
     END IF

  END IF

  CALL stop_clock_gpu( clock_label )
  
  RETURN
  !
END SUBROUTINE env_fwfft_y_gpu
!=---------------------------------------------------------------------------=!

#endif
