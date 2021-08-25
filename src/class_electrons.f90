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
MODULE class_electrons
    !------------------------------------------------------------------------------------
    !
    USE env_base_io, ONLY: ionode, environ_unit, verbose, depth
    !
    USE environ_param, ONLY: DP
    !
    USE class_cell
    !
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
    TYPE, PUBLIC :: environ_electrons
        !--------------------------------------------------------------------------------
        !
        LOGICAL :: lupdate
        LOGICAL :: initialized
        INTEGER :: number
        !
        TYPE(environ_density) :: density
        REAL(DP) :: charge
        !
        !--------------------------------------------------------------------------------
    CONTAINS
        !--------------------------------------------------------------------------------
        !
        PROCEDURE, PRIVATE :: create => create_environ_electrons
        PROCEDURE :: init_first => init_environ_electrons_first
        PROCEDURE :: init_second => init_environ_electrons_second
        PROCEDURE :: update => update_environ_electrons
        PROCEDURE :: destroy => destroy_environ_electrons
        !
        PROCEDURE :: printout => print_environ_electrons
        !
        !--------------------------------------------------------------------------------
    END TYPE environ_electrons
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
    SUBROUTINE create_environ_electrons(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_electrons), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: sub_name = 'create_environ_electrons'
        !
        !--------------------------------------------------------------------------------
        !
        this%lupdate = .FALSE.
        this%initialized = .FALSE.
        !
        this%number = 0
        this%charge = 0.D0
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE create_environ_electrons
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_environ_electrons_first(this, nelec)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: nelec
        !
        CLASS(environ_electrons), INTENT(INOUT) :: this
        !
        !--------------------------------------------------------------------------------
        !
        CALL this%create()
        !
        this%number = nelec
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_environ_electrons_first
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_environ_electrons_second(this, cell)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_cell), INTENT(IN) :: cell
        !
        CLASS(environ_electrons), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: local_label = 'electrons'
        !
        !--------------------------------------------------------------------------------
        !
        CALL this%density%init(cell, local_label)
        !
        this%initialized = .TRUE.
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_environ_electrons_second
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE update_environ_electrons(this, nnr, rho, nelec)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: nnr
        REAL(DP), INTENT(IN) :: rho(nnr)
        !
        REAL(DP), INTENT(IN), OPTIONAL :: nelec
        !
        CLASS(environ_electrons), INTENT(INOUT) :: this
        !
        REAL(DP), PARAMETER :: tol = 1.D-4
        REAL(DP) :: charge
        !
        CHARACTER(LEN=80) :: sub_name = 'update_environ_electrons'
        !
        !--------------------------------------------------------------------------------
        !
        ! check on dimensions
        IF (nnr /= this%density%cell%nnr) &
            CALL env_errore(sub_name, 'Mismatch in grid size', 1)
        !
        this%density%of_r = rho
        !
        !--------------------------------------------------------------------------------
        ! Update integral of electronic density and, if provided, check
        ! against input value
        !
        this%charge = this%density%integrate()
        this%number = NINT(this%charge)
        !
        IF (PRESENT(nelec)) THEN
            !
            IF (ABS(this%charge - nelec) > tol) &
                CALL env_errore(sub_name, 'Mismatch in integrated electronic charge', 1)
            !
        END IF
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE update_environ_electrons
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE destroy_environ_electrons(this, lflag)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        LOGICAL, INTENT(IN) :: lflag
        !
        CLASS(environ_electrons), INTENT(INOUT) :: this
        !
        !--------------------------------------------------------------------------------
        !
        IF (this%initialized) THEN
            !
            CALL this%density%destroy()
            !
            this%charge = 0.D0
            this%initialized = .FALSE.
        END IF
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE destroy_environ_electrons
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
    SUBROUTINE print_environ_electrons(this, local_verbose, local_depth)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_electrons), INTENT(IN) :: this
        INTEGER, INTENT(IN), OPTIONAL :: local_verbose
        INTEGER, INTENT(IN), OPTIONAL :: local_depth
        !
        INTEGER :: verbosity, passed_verbosity, passed_depth
        !
        CHARACTER(LEN=80) :: sub_name = 'print_environ_electrons'
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
            IF (verbosity >= 3) &
                CALL this%density%printout(passed_verbosity, passed_depth)
            !
            IF (verbosity < verbose) CALL env_block_divider(verbosity)
            !
        END IF
        !
        FLUSH (environ_unit)
        !
        !--------------------------------------------------------------------------------
        !
1000    FORMAT(/, 4('%'), ' ELECTRONS ', 65('%'))
1001    FORMAT(/, ' ELECTRONS', /, ' =========')
        !
1002    FORMAT(/, ' number of electrons        = ', I10)
        !
1003    FORMAT(/, ' total electronic charge    = ', F14.7)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE print_environ_electrons
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE class_electrons
!----------------------------------------------------------------------------------------
