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
!!
!----------------------------------------------------------------------------------------
MODULE env_fft_param
    !------------------------------------------------------------------------------------
    !
    USE env_utils_param
    !
    !------------------------------------------------------------------------------------
    !
    INTEGER, PARAMETER :: ndims = 10
    ! number of different FFT tables that the module could keep into memory
    ! without reinitialization
    !
    INTEGER, PARAMETER :: nfftx = 2049
    ! max allowed fft dimension
    !
    !------------------------------------------------------------------------------------
END MODULE env_fft_param
!----------------------------------------------------------------------------------------
