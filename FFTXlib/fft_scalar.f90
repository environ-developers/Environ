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
! internal FFTW, FFTW v.3, IBM ESSL, Intel DFTI, ARMlib                    !
! (both 3d for serial execution and 1d+2d FFTs for parallel execution);    !
! legacy NEC ASL libraries (3d only, no parallel execution)                !
! Written by Carlo Cavazzoni, modified by P. Giannozzi, contributions      !
! by Martin Hilgemans, Guido Roma, Pascal Thibaudeau, Stephane Lefranc,    !
! Nicolas Lacorne, Filippo Spiga, Nicola Varini, Jason Wood                !
! Last update Oct 2017                                                     !
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
#elif defined(__ARM_LIB)
     USE env_fft_scalar_arm
#else
     USE env_fft_scalar_fftw
#endif
       
     IMPLICIT NONE
     SAVE

     PRIVATE
     PUBLIC :: env_cft_1z, env_cft_2xy, env_cfft3d, env_cfft3ds

   END MODULE env_fft_scalar
