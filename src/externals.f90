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
!          Nicola Marzari     (THEOS and NCCR-MARVEL, EPFL)
!          Edan Bainglass     (Department of Physics, UNT)
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
MODULE class_externals
    !------------------------------------------------------------------------------------
    !
    USE class_io, ONLY: io
    !
    USE environ_param, ONLY: DP
    !
    USE class_cell
    USE class_density
    USE class_function
    USE class_function_gaussian
    USE class_functions
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
    TYPE, PUBLIC :: environ_externals
        !--------------------------------------------------------------------------------
        !
        LOGICAL :: lupdate = .FALSE.

        INTEGER :: number = 0
        REAL(DP) :: charge = 0.D0
        !
        TYPE(environ_density) :: density
        !
        TYPE(environ_functions) :: functions
        !
        !--------------------------------------------------------------------------------
    CONTAINS
        !--------------------------------------------------------------------------------
        !
        PROCEDURE :: init => init_environ_externals
        PROCEDURE :: update => update_environ_externals
        PROCEDURE :: destroy => destroy_environ_externals
        !
        PROCEDURE :: printout => print_environ_externals
        !
        !--------------------------------------------------------------------------------
    END TYPE environ_externals
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
    SUBROUTINE init_environ_externals(this, n, dims, axes, pos, spreads, charges, cell)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: n
        INTEGER, DIMENSION(n), INTENT(IN) :: dims, axes
        REAL(DP), DIMENSION(n), INTENT(IN) :: spreads, charges
        REAL(DP), INTENT(IN) :: pos(3, n)
        TYPE(environ_cell), INTENT(IN) :: cell
        !
        CLASS(environ_externals), INTENT(INOUT) :: this
        !
        !--------------------------------------------------------------------------------
        !
        CALL this%density%init(cell, 'externals')
        !
        this%number = n
        !
        IF (n > 0) &
            CALL this%functions%init(n, 1, axes, dims, spreads, spreads, -charges, pos)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_environ_externals
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE update_environ_externals(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_externals), INTENT(INOUT) :: this
        !
        !--------------------------------------------------------------------------------
        !
        CALL this%functions%density(this%density, .TRUE.)
        !
        this%charge = this%density%integrate()
        !
        !--------------------------------------------------------------------------------
        ! Output current state
        !
        IF (.NOT. this%lupdate) CALL this%printout()
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE update_environ_externals
    !------------------------------------------------------------------------------------
    !!
    !>
    !------------------------------------------------------------------------------------
    SUBROUTINE destroy_environ_externals(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_externals), INTENT(INOUT) :: this
        !
        !--------------------------------------------------------------------------------
        !
        CALL this%density%destroy()
        !
        CALL this%functions%destroy()
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE destroy_environ_externals
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !                                   OUTPUT METHODS
    !
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !>
    !! Prints the details of the external charges
    !!
    !! Nested objects receive a decremented passed verbose to trigger block printing
    !!
    !! @param verbose       : (INTEGER) adds verbosity to global verbose
    !! @param debug_verbose : (INTEGER) replaces global verbose for debugging
    !! @param unit          : (INTEGER) output target (default = io%debug_unit)
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE print_environ_externals(this, verbose, debug_verbose, unit)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_externals), INTENT(IN) :: this
        INTEGER, OPTIONAL, INTENT(IN) :: verbose, debug_verbose, unit
        !
        INTEGER :: base_verbose, local_verbose, passed_verbose, local_unit
        !
        CHARACTER(LEN=80) :: sub_name = 'print_environ_externals'
        !
        !--------------------------------------------------------------------------------
        !
        IF (PRESENT(debug_verbose)) THEN
            base_verbose = debug_verbose
            !
            IF (PRESENT(verbose)) THEN
                local_verbose = verbose
            ELSE
                local_verbose = debug_verbose
            END IF
            !
            passed_verbose = verbose - 1
            !
        ELSE IF (io%verbosity > 0) THEN
            base_verbose = io%verbosity
            !
            IF (PRESENT(verbose)) THEN
                local_verbose = base_verbose + verbose
            ELSE
                local_verbose = base_verbose
            END IF
            !
            passed_verbose = local_verbose - base_verbose - 1
            !
        ELSE
            RETURN
        END IF
        !
        IF (PRESENT(unit)) THEN
            local_unit = unit
        ELSE
            local_unit = io%debug_unit
        END IF
        !
        IF (local_verbose >= 1) THEN
            !
            IF (io%lnode) THEN
                WRITE (local_unit, 1000)
                WRITE (local_unit, 1001) this%number
                WRITE (local_unit, 1002) this%charge
            END IF
            !
            IF (local_verbose >= 3) THEN
                !
                CALL this%functions%printout(passed_verbose, debug_verbose, local_unit)
                !
                CALL this%density%printout(passed_verbose, debug_verbose, local_unit)
                !
            END IF
            !
        END IF
        !
        FLUSH (local_unit)
        !
        !--------------------------------------------------------------------------------
        !
1000    FORMAT(/, 4('%'), ' EXTERNALS ', 65('%'))
        !
1001    FORMAT(/, ' number of external charges = ', I14)
        !
1002    FORMAT(/, ' total external charge      = ', F14.7)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE print_environ_externals
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE class_externals
!----------------------------------------------------------------------------------------
