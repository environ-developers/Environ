!
! Copyright (C) Quantum ESPRESSO group
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------!
! FFT scalar drivers Module - contains machine-dependent routines for      !
! internal FFTW, FFTW v.3, IBM ESSL, Intel DFTI
! (both 3d for serial execution and 1d+2d FFTs for parallel execution);    !
! legacy NEC ASL libraries (3d only, no parallel execution)                !
! CUDA FFT for NVidiia GPUs
! Written by Carlo Cavazzoni, modified by P. Giannozzi, contributions      !
! by Martin Hilgemans, Guido Roma, Pascal Thibaudeau, Stephane Lefranc,    !
! Nicolas Lacorne, Filippo Spiga, Nicola Varini, Jason Wood                !
! Last update Feb 2021
!--------------------------------------------------------------------------!

!=----------------------------------------------------------------------=!
   MODULE env_fft_scalar
!=----------------------------------------------------------------------=!

     USE env_fft_param
#if defined(__FFTW3)
     USE env_fft_scalar_fftw3
#elif defined(__DFTI)
     USE env_fft_scalar_dfti
#elif defined(__LINUX_ESSL)
     USE env_fft_scalar_essl
#elif defined(__SX6)
     USE env_fft_scalar_sx6
#else
#error No env_fft_scalar backend selected!
#endif
#if defined(__CUDA)
     USE env_fft_scalar_cuFFT
#endif
     IMPLICIT NONE
     SAVE

     PRIVATE
     PUBLIC :: env_cft_1z, env_cft_2xy, env_cfft3d, env_cfft3ds
#if defined(__CUDA)
     PUBLIC :: env_cft_1z_gpu, env_cft_2xy_gpu, env_cfft3d_gpu, env_cfft3ds_gpu
#endif

   END MODULE env_fft_scalar
