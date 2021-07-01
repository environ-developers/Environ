!----------------------------------------------------------------------------------------
!
! Copyright (C) 2018 ENVIRON (www.quantum-environment.org)
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
MODULE env_base_io
    !------------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    SAVE
    !
    LOGICAL :: ionode = .TRUE.
    INTEGER :: ionode_id
    !
    LOGICAL :: lstdout ! whether environ can print on standard output
    !
    INTEGER :: comm ! WE MAY NEED A SECOND COMMUNICATOR FOR IMAGE PARALLELIZATION
    !
    INTEGER :: program_unit
    INTEGER :: environ_unit
    INTEGER :: verbose
    INTEGER :: depth = 1
    !
    CHARACTER(LEN=2) :: prog ! the calling program
    !
    !------------------------------------------------------------------------------------
END MODULE env_base_io
!----------------------------------------------------------------------------------------
