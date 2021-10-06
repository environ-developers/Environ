!----------------------------------------------------------------------------------------
!
! Copyright (C) 2021 ENVIRON (www.quantum-environ.org)
! Copyright (C) 2003-2004 Quantum ESPRESSO (www.quantum-espresso.org)
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
! Authors: Carlo Cavazzoni
!          Modified by Edan Bainglass
!
!----------------------------------------------------------------------------------------
!>
!! SISSA Code Interface
!!
!----------------------------------------------------------------------------------------
MODULE env_parallel_include
    !------------------------------------------------------------------------------------
    !
#if defined(__MPI)
#if defined(__MPI_MODULE)
    USE mpi
#else
    INCLUDE 'mpif.h'
#endif
#else
    INTEGER, PARAMETER :: MPI_COMM_WORLD = 0
    INTEGER, PARAMETER :: MPI_COMM_NULL = -1
    INTEGER, PARAMETER :: MPI_COMM_SELF = -2
#endif
    !
    !------------------------------------------------------------------------------------
END MODULE env_parallel_include
!----------------------------------------------------------------------------------------
