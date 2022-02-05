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
! Authors: Carlo Cavazzoni, modified by P. Giannozzi, contributions
!          by Martin Hilgemans, Guido Roma, Pascal Thibaudeau, Stephane Lefranc,
!          Nicolas Lacorne, Filippo Spiga, Nicola Varini, Jason Wood
!          Modified by Edan Bainglass
!
!----------------------------------------------------------------------------------------
!>
!! FFT scalar drivers contains machine-dependent routines for
!! internal FFTW, FFTW v.3, IBM ESSL, Intel DFTI, ARMlib, and NEC ASL
!!
!----------------------------------------------------------------------------------------
MODULE env_fft_scalar
    !------------------------------------------------------------------------------------
    !
#if defined(__FFTW3)
    USE env_fft_scalar_fftw3, ONLY: env_cft_1z, env_cft_2xy, env_cfft3d
#elif defined(__DFTI)
    USE env_fft_scalar_dfti, ONLY: env_cft_1z, env_cft_2xy, env_cfft3d
#elif defined(__LINUX_ESSL)
    USE env_fft_scalar_essl, ONLY: env_cft_1z, env_cft_2xy, env_cfft3d
#elif defined(__SX6)
    USE env_fft_scalar_sx6, ONLY: env_cft_1z, env_cft_2xy, env_cfft3d
#elif defined(__ARM_LIB)
    USE env_fft_scalar_arm, ONLY: env_cft_1z, env_cft_2xy, env_cfft3d
#elif defined(__FFTW)
    USE env_fft_scalar_fftw, ONLY: env_cft_1z, env_cft_2xy, env_cfft3d
#else
#error No fft scalar backend selected!
#endif
#if defined(__CUDA)
    USE env_fft_scalar_cuFFT, ONLY: env_cft_1z_gpu, env_cft_2xy_gpu, env_cfft3d_gpu
#endif
    !
    !------------------------------------------------------------------------------------
END MODULE env_fft_scalar
!----------------------------------------------------------------------------------------
