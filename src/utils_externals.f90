! Copyright (C) 2018 ENVIRON (www.quantum-environment.org)
!
!    This file is part of Environ version 1.1
!
!    Environ 1.1 is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 2 of the License, or
!    (at your option) any later version.
!
!    Environ 1.1 is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more detail, either the file
!    `License' in the root directory of the present distribution, or
!    online at <http://www.gnu.org/licenses/>.
!
!----------------------------------------------------------------------------------------
!
! Authors: Oliviero Andreussi (Department of Physics, UNT)
!          Nicola Marzari     (THEOS and NCCR-MARVEL, EPFL)
!
!----------------------------------------------------------------------------------------
!>
!! Module containing the main routines to handle environ_externals
!! derived data types.
!!
!! Environ_externals contains all the specifications and the details of
!! the external user-defined (thus fixed) smooth charge distributions
!! introduced in the Environ simulation cell. Distributions of gaussian
!! shape and different dimensionalities (0D,1D,2D) are available.
!!
!----------------------------------------------------------------------------------------
MODULE utils_externals
    !------------------------------------------------------------------------------------
    !
    USE environ_param, ONLY: DP
    !
    USE types_physical, ONLY: environ_externals
    USE types_cell, ONLY: environ_cell
    !
    USE utils_functions, ONLY: destroy_environ_functions
    !
    USE utils_density, ONLY: create_environ_density, init_environ_density, &
                             destroy_environ_density
    !
    USE tools_functions, ONLY: density_of_functions
    USE tools_math, ONLY: integrate_environ_density
    !
    !------------------------------------------------------------------------------------
    !
    PRIVATE
    !
    PUBLIC :: create_environ_externals, init_environ_externals_first, &
              init_environ_externals_second, update_environ_externals, &
              destroy_environ_externals
    !
    !------------------------------------------------------------------------------------
CONTAINS
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE create_environ_externals(externals)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_externals), INTENT(INOUT) :: externals
        !
        CHARACTER(LEN=80) :: label = 'externals'
        !
        CHARACTER(LEN=80) :: sub_name = 'create_environ_externals'
        !
        !--------------------------------------------------------------------------------
        !
        externals%update = .FALSE.
        externals%number = 0
        !
        IF (ALLOCATED(externals%functions)) &
            CALL env_errore(sub_name, 'Trying to create an already allocated object', 1)
        !
        CALL create_environ_density(externals%density, label)
        !
        externals%charge = 0.D0
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE create_environ_externals
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_environ_externals_first(nexternals, dims, axis, pos, spreads, &
                                            charge, externals)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: nexternals
        INTEGER, DIMENSION(nexternals), INTENT(IN) :: dims, axis
        REAL(DP), INTENT(IN) :: pos(3, nexternals)
        REAL(DP), DIMENSION(nexternals), INTENT(IN) :: spreads, charge
        !
        TYPE(environ_externals), INTENT(INOUT) :: externals
        !
        INTEGER :: i
        !
        !--------------------------------------------------------------------------------
        !
        externals%number = nexternals
        ALLOCATE (externals%functions(externals%number))
        !
        DO i = 1, externals%number
            ALLOCATE (externals%functions(i)%pos(3))
            externals%functions(i)%type_ = 1
            externals%functions(i)%dim = dims(i)
            externals%functions(i)%axis = axis(i)
            externals%functions(i)%pos(:) = pos(:, i)
            externals%functions(i)%spread = spreads(i)
            externals%functions(i)%width = spreads(i)
            externals%functions(i)%volume = -charge(i)
        END DO
        !
        externals%initialized = .FALSE.
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_environ_externals_first
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_environ_externals_second(cell, externals)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_cell), INTENT(IN) :: cell
        !
        TYPE(environ_externals), INTENT(INOUT) :: externals
        !
        INTEGER :: i
        !
        !--------------------------------------------------------------------------------
        !
        IF (externals%number > 0) THEN
            !
            DO i = 1, externals%number
                externals%functions(i)%pos = externals%functions(i)%pos / cell%alat
            END DO
            !
        END IF
        !
        CALL init_environ_density(cell, externals%density)
        !
        externals%initialized = .TRUE.
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_environ_externals_second
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE update_environ_externals(externals)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_externals), INTENT(INOUT) :: externals
        !
        !--------------------------------------------------------------------------------
        !
        CALL density_of_functions(externals%number, externals%functions, &
                                  externals%density, .TRUE.)
        !
        externals%charge = integrate_environ_density(externals%density)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE update_environ_externals
    !------------------------------------------------------------------------------------
    !
    !>
    !------------------------------------------------------------------------------------
    SUBROUTINE destroy_environ_externals(lflag, externals)
        !--------------------------------------------------------------------------------
        !!
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        LOGICAL, INTENT(IN) :: lflag
        !
        TYPE(environ_externals), INTENT(INOUT) :: externals
        !
        !--------------------------------------------------------------------------------
        !
        IF (lflag) CALL destroy_environ_functions(externals%number, externals%functions)
        !
        IF (externals%initialized) THEN
            !
            CALL destroy_environ_density(externals%density)
            !
            externals%initialized = .FALSE.
        END IF
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE destroy_environ_externals
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE utils_externals
!----------------------------------------------------------------------------------------
