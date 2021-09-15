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
MODULE env_data_buffer
    !------------------------------------------------------------------------------------
    !
    USE env_kinds, ONLY: DP
    !
#if defined(__CUDA)
    USE cudafor
#endif
    !
    !------------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    REAL(DP), ALLOCATABLE :: mp_buff_r(:)
    INTEGER, ALLOCATABLE :: mp_buff_i(:)
    !
#if defined(__CUDA)
    REAL(DP), ALLOCATABLE :: mp_buff_r_d(:)
    INTEGER, ALLOCATABLE :: mp_buff_i_d(:)
    ATTRIBUTES(DEVICE) :: mp_buff_r_d, mp_buff_i_d
    !
    PUBLIC :: mp_buff_r_d, mp_buff_i_d
    !
#endif
    !
    !------------------------------------------------------------------------------------
    !
    PRIVATE
    !
    PUBLIC :: mp_buff_r, mp_buff_i
    !
    !------------------------------------------------------------------------------------
END MODULE env_data_buffer
!----------------------------------------------------------------------------------------
