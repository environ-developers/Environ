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
    USE env_base_io, ONLY: ionode, environ_unit, verbose, depth
    !
    USE environ_param, ONLY: DP
    !
    USE class_cell
    USE class_density
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
        LOGICAL :: initialized = .FALSE.
        INTEGER :: number = 0
        !
        TYPE(environ_function), ALLOCATABLE :: functions(:)
        TYPE(environ_density) :: density
        REAL(DP) :: charge = 0.0_DP
        !
        !--------------------------------------------------------------------------------
    CONTAINS
        !--------------------------------------------------------------------------------
        !
        PROCEDURE :: create => create_environ_externals
        PROCEDURE :: init_first => init_environ_externals_first
        PROCEDURE :: init_second => init_environ_externals_second
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
    SUBROUTINE create_environ_externals(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_externals), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: label = 'externals'
        !
        CHARACTER(LEN=80) :: sub_name = 'create_environ_externals'
        !
        !--------------------------------------------------------------------------------
        !
        this%lupdate = .FALSE.
        this%number = 0
        !
        IF (ALLOCATED(this%functions)) &
            CALL env_errore(sub_name, 'Trying to create an already allocated object', 1)
        !
        CALL this%density%create(label)
        !
        this%charge = 0.D0
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE create_environ_externals
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_environ_externals_first(this, nexternals, dims, axes, pos, &
                                            spreads, charges)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: nexternals
        INTEGER, DIMENSION(nexternals), INTENT(IN) :: dims, axes
        REAL(DP), DIMENSION(nexternals), INTENT(IN) :: spreads, charges
        REAL(DP), INTENT(IN) :: pos(3, nexternals)
        !
        CLASS(environ_externals), INTENT(INOUT) :: this
        !
        !--------------------------------------------------------------------------------
        !
        this%number = nexternals
        !
        IF (this%number > 0) THEN
            !
            CALL init_environ_functions(this%functions, nexternals, 1, &
                                        axes, dims, spreads, spreads, -charges, pos)
            !
        END IF
        !
        this%initialized = .FALSE.
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_environ_externals_first
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_environ_externals_second(this, cell)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_cell), INTENT(IN) :: cell
        !
        CLASS(environ_externals), INTENT(INOUT) :: this
        !
        INTEGER :: i
        !
        !--------------------------------------------------------------------------------
        !
        IF (this%number > 0) THEN
            !
            DO i = 1, this%number
                this%functions(i)%pos = this%functions(i)%pos
            END DO
            !
        END IF
        !
        CALL this%density%init(cell)
        !
        this%initialized = .TRUE.
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_environ_externals_second
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
        CALL density_of_functions(this%functions, this%number, this%density, .TRUE.)
        !
        this%charge = this%density%integrate()
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE update_environ_externals
    !------------------------------------------------------------------------------------
    !
    !>
    !------------------------------------------------------------------------------------
    SUBROUTINE destroy_environ_externals(this, lflag)
        !--------------------------------------------------------------------------------
        !!
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        LOGICAL, INTENT(IN) :: lflag
        !
        CLASS(environ_externals), INTENT(INOUT) :: this
        !
        !--------------------------------------------------------------------------------
        !
        IF (lflag) CALL destroy_environ_functions(this%functions, this%number)
        !
        IF (this%initialized) THEN
            !
            CALL this%density%destroy()
            !
            this%initialized = .FALSE.
        END IF
        !
        RETURN
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
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE print_environ_externals(this, local_verbose, local_depth)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_externals), INTENT(IN) :: this
        INTEGER, INTENT(IN), OPTIONAL :: local_verbose
        INTEGER, INTENT(IN), OPTIONAL :: local_depth
        !
        INTEGER :: verbosity, passed_verbosity, passed_depth
        !
        CHARACTER(LEN=80) :: sub_name = 'print_environ_externals'
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
        IF (verbosity == 0) RETURN
        !
        IF (PRESENT(local_depth)) THEN
            passed_verbosity = verbosity - verbose - local_depth
            passed_depth = local_depth
        ELSE
            passed_verbosity = verbosity - verbose - depth
            passed_depth = depth
        END IF
        !
        IF (verbosity >= 1) THEN
            !
            IF (ionode) THEN
                !
                IF (verbosity >= verbose) THEN ! header
                    WRITE (environ_unit, 1000)
                ELSE
                    !
                    CALL env_block_divider(verbosity)
                    !
                    WRITE (environ_unit, 1001)
                END IF
                !
                WRITE (environ_unit, 1002) this%number
                WRITE (environ_unit, 1003) this%charge
                !
            END IF
            !
            IF (verbosity >= 3) THEN
                !
                CALL print_environ_functions(this%functions, this%number, &
                                             passed_verbosity, passed_depth)
                !
                CALL this%density%printout(passed_verbosity, passed_depth)
                !
            END IF
            !
            IF (verbosity < verbose) CALL env_block_divider(verbosity)
            !
        END IF
        !
        FLUSH (environ_unit)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
        !
1000    FORMAT(/, 4('%'), ' EXTERNALS ', 65('%'))
1001    FORMAT(/, ' EXTERNALS', /, ' =========')
        !
1002    FORMAT(/, ' number of external charges = ', I10)
        !
1003    FORMAT(/, ' total external charge      = ', F14.7)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE print_environ_externals
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE class_externals
!----------------------------------------------------------------------------------------
