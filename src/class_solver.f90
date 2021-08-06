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
! Authors: Francesco Nattino  (THEOS and NCCR-MARVEL, EPFL)
!          Oliviero Andreussi (Department of Physics, UNT)
!          Nicola Marzari     (THEOS and NCCR-MARVEL, EPFL)
!          Edan Bainglass     (Department of Physics, UNT)
!
!----------------------------------------------------------------------------------------
!>
!!
!----------------------------------------------------------------------------------------
MODULE class_solver
    !------------------------------------------------------------------------------------
    !
    USE class_core_container_electrostatics
    !
    !------------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    PRIVATE
    !
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    TYPE, ABSTRACT, PUBLIC :: electrostatic_solver
        !--------------------------------------------------------------------------------
        !
        CHARACTER(LEN=80) :: solver_type
        !
        TYPE(container_electrostatics), POINTER :: cores => NULL()
        !
        !--------------------------------------------------------------------------------
    CONTAINS
        !--------------------------------------------------------------------------------
        !
        PROCEDURE :: set_cores => set_electrostatic_cores
        !
        !--------------------------------------------------------------------------------
    END TYPE electrostatic_solver
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
CONTAINS
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !                                   ADMIN METHODS
    !
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE set_electrostatic_cores(this, cores)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(container_electrostatics), TARGET, INTENT(IN) :: cores
        !
        CLASS(electrostatic_solver), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: sub_name = 'set_electrostatic_cores'
        !
        !--------------------------------------------------------------------------------
        !
        IF (ASSOCIATED(this%cores)) &
            CALL env_errore(sub_name, 'Trying to associate an associated core', 1)
        !
        this%cores => cores
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE set_electrostatic_cores
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE class_solver
!----------------------------------------------------------------------------------------
