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
!!
!----------------------------------------------------------------------------------------
MODULE env_write_cube
    !------------------------------------------------------------------------------------
    !
    USE env_base_io, ONLY: ionode, verbose
    !
    USE class_density
    !
    USE class_ions
    !
    !------------------------------------------------------------------------------------
CONTAINS
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE write_cube(density, ions, local_verbose, local_depth, idx, label)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_density), INTENT(IN) :: density
        TYPE(environ_ions), INTENT(IN) :: ions
        INTEGER, INTENT(IN), OPTIONAL :: local_verbose, local_depth, idx
        CHARACTER(LEN=100), INTENT(IN), OPTIONAL :: label
        !
        CHARACTER(LEN=100) :: filename, filemod, local_label
        !
        INTEGER :: verbosity
        !
        !--------------------------------------------------------------------------------
        !
        IF (verbose == 0) RETURN
        !
        IF (PRESENT(local_verbose)) THEN
            verbosity = verbose + local_verbose
        ELSE
            verbosity = verbose
        END IF
        !
        IF (PRESENT(idx)) THEN
            WRITE (filemod, '(i4.4)') idx
        ELSE
            filemod = ''
        END IF
        !
        IF (PRESENT(label)) THEN
            local_label = label
        ELSE
            local_label = density%label
        END IF
        !
        filename = TRIM(ADJUSTL(local_label))//TRIM(filemod)//'.cube'
        !
        !--------------------------------------------------------------------------------
        ! Write cube
        !
        CALL density%printout(local_verbose, local_depth, .FALSE.)
        !
        IF (verbosity >= 3) THEN
            !
            OPEN (300, file=TRIM(filename), status='unknown')
            !
            CALL density%cell%write_cube(ions%number) ! write cube cell data
            !
            CALL ions%write_cube() ! write cube ion data
            !
            CALL density%write_cube()
            !
            CLOSE (300)
            !
        END IF
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE write_cube
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE env_write_cube
!----------------------------------------------------------------------------------------