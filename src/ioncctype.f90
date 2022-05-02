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
! Authors: Oliviero Andreussi (Department of Physics, UNT)
!          Francesco Nattino  (THEOS and NCCR-MARVEL, EPFL)
!          Nicola Marzari     (THEOS and NCCR-MARVEL, EPFL)
!          Edan Bainglass     (Department of Physics, UNT)
!
!----------------------------------------------------------------------------------------
!>
!!
!----------------------------------------------------------------------------------------
MODULE class_ioncctype
    !------------------------------------------------------------------------------------
    !
    USE environ_param, ONLY: DP, BOHR_RADIUS_SI, AMU_SI
    !
    USE class_cell
    USE class_density
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
    TYPE, PUBLIC :: environ_ioncctype
        !--------------------------------------------------------------------------------
        !
        INTEGER :: index
        REAL(DP) :: cbulk ! bulk concentration
        REAL(DP) :: z ! charge
        !
        TYPE(environ_density) :: c ! local concentration
        TYPE(environ_density) :: cfactor ! exp(-z\phi\beta) or 1 - z\phi\beta
        !
        !--------------------------------------------------------------------------------
    CONTAINS
        !--------------------------------------------------------------------------------
        !
        PROCEDURE, PRIVATE :: create => create_environ_ioncctype
        PROCEDURE :: init => init_environ_ioncctype
        PROCEDURE :: destroy => destroy_environ_ioncctype
        !
        !--------------------------------------------------------------------------------
    END TYPE environ_ioncctype
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
    SUBROUTINE create_environ_ioncctype(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_ioncctype), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: routine = 'create_environ_ioncctype'
        !
        !--------------------------------------------------------------------------------
        !
        this%index = 0
        this%cbulk = 0.D0
        this%z = 0.D0
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE create_environ_ioncctype
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_environ_ioncctype(this, index, cbulk, z, cell)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: index
        REAL(DP), INTENT(IN) :: cbulk, z
        TYPE(environ_cell), INTENT(IN) :: cell
        !
        CLASS(environ_ioncctype), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: index_s
        !
        CHARACTER(LEN=80) :: routine = 'init_environ_ioncctype'
        !
        !--------------------------------------------------------------------------------
        !
        CALL this%create()
        !
        this%index = index
        this%z = -z
        !
        this%cbulk = cbulk * BOHR_RADIUS_SI**3 / AMU_SI
        ! convert bulk concentrations to atomic units
        !
        WRITE (index_s, '(I2.2)') index
        !
        CALL this%c%init(cell, 'c_electrolyte_'//TRIM(index_s))
        !
        CALL this%cfactor%init(cell, 'cfactor_electrolyte_'//TRIM(index_s))
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_environ_ioncctype
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE destroy_environ_ioncctype(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_ioncctype), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: routine = 'destroy_environ_ioncctype'
        !
        !--------------------------------------------------------------------------------
        !
        CALL this%c%destroy()
        !
        CALL this%cfactor%destroy()
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE destroy_environ_ioncctype
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE class_ioncctype
!----------------------------------------------------------------------------------------
