!----------------------------------------------------------------------------------------
!
! Copyright (C) 2018-2021 ENVIRON (www.quantum-environ.org)
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
! Authors: Oliviero Andreussi (Department of Physics, UNT)
!          Francesco Nattino  (THEOS and NCCR-MARVEL, EPFL)
!          Nicola Marzari     (THEOS and NCCR-MARVEL, EPFL)
!          Edan Bainglass     (Department of Physics, UNT)
!
!----------------------------------------------------------------------------------------
!>
!! Internal setup of numerical cores
!!
!----------------------------------------------------------------------------------------
MODULE base_core
    !------------------------------------------------------------------------------------
    !
    USE types_core, ONLY: fd_core, fft_core, oned_analytic_core
    !
    !------------------------------------------------------------------------------------
    !
    SAVE
    !
    LOGICAL :: lfd
    TYPE(fd_core) :: fd
    !
    LOGICAL :: lfft_system
    TYPE(fft_core) :: system_fft
    !
    LOGICAL :: lfft_environment
    TYPE(fft_core) :: environment_fft
    !
    LOGICAL :: loned_analytic
    TYPE(oned_analytic_core) :: oned_analytic
    !
    !------------------------------------------------------------------------------------
    ! Keeping imports private
    !
    PRIVATE :: fd_core, fft_core, oned_analytic_core
    !
    !------------------------------------------------------------------------------------
END MODULE base_core
!----------------------------------------------------------------------------------------
