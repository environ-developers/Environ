!----------------------------------------------------------------------------------------
!
! Copyright (C) 2021 ENVIRON (www.quantum-environ.org)
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
! Authors: Edan Bainglass (Department of Physics, UNT)
!
!----------------------------------------------------------------------------------------
MODULE cmdline_args
    !
    USE environ_param, ONLY: DP
    !------------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=80) :: inputfile = 'environ.in'
    CHARACTER(LEN=80) :: cubefile = 'density.cube'
    !
    LOGICAL :: use_pbc_corr = .FALSE.
    !
    LOGICAL :: no_density = .FALSE.
    !
    LOGICAL :: calc_energy = .FALSE.
    !
    LOGICAL :: calc_force = .FALSE.
    !
    REAL(DP) :: alpha_min = -1.D0
    REAL(DP) :: alpha_max = -1.D0
    REAL(DP) :: alpha_step = -1.D0
    !
    !------------------------------------------------------------------------------------
END MODULE cmdline_args
!----------------------------------------------------------------------------------------
