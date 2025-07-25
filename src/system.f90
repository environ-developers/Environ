!----------------------------------------------------------------------------------------
!
! Copyright (C) 2018-2022 ENVIRON (www.quantum-environ.org)
!
!----------------------------------------------------------------------------------------
!
!     This file is part of Environ version 3.0
!
!     Environ 3.0 is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation, either version 2 of the License, or
!     (at your option) any later version.
!
!     Environ 3.0 is distributed in the hope that it will be useful,
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
MODULE class_system
    !------------------------------------------------------------------------------------
    !
    USE class_io, ONLY: io
    !
    USE environ_param, ONLY: DP
    !
    USE class_ions
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
    TYPE, PUBLIC :: environ_system
        !--------------------------------------------------------------------------------
        !
        LOGICAL :: initialized = .FALSE.
        !
        !--------------------------------------------------------------------------------
        !
        LOGICAL :: lupdate = .FALSE.
        !
        INTEGER :: ntyp = 0
        INTEGER :: dim = 0
        INTEGER :: axis = 1
        REAL(DP) :: width = 0.D0
        REAL(DP) :: com(3) = 0.D0 ! center of mass
        !
        TYPE(environ_ions), POINTER :: ions => NULL()
        !
        !--------------------------------------------------------------------------------
    CONTAINS
        !--------------------------------------------------------------------------------
        !
        PROCEDURE, PRIVATE :: create => create_environ_system
        PROCEDURE :: init => init_environ_system
        PROCEDURE :: update => update_environ_system
        PROCEDURE :: destroy => destroy_environ_system
        !
        PROCEDURE :: printout => print_environ_system
        !
        !--------------------------------------------------------------------------------
    END TYPE environ_system
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
    SUBROUTINE create_environ_system(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_system), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: routine = 'create_environ_system'
        !
        !--------------------------------------------------------------------------------
        !
        IF (ASSOCIATED(this%ions)) CALL io%create_error(routine)
        !
        !--------------------------------------------------------------------------------
        !
        this%lupdate = .FALSE.
        this%ntyp = 0
        this%dim = 0
        this%axis = 0
        this%width = 0.D0
        this%com = 0.D0
        !
        NULLIFY (this%ions)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE create_environ_system
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_environ_system(this, ntyp, dim, axis, ions)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: ntyp, dim, axis
        TYPE(environ_ions), TARGET, INTENT(IN) :: ions
        !
        CLASS(environ_system), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: routine = 'init_environ_system'
        !
        !--------------------------------------------------------------------------------
        !
        CALL this%create()
        !
        this%ntyp = ntyp
        this%dim = dim
        this%axis = axis
        !
        this%ions => ions
        !
        this%initialized = .TRUE.
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_environ_system
    !------------------------------------------------------------------------------------
    !>
    !! Given the system definition compute position (centre of charge)
    !! and width (maximum distance from centre) of the system.
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE update_environ_system(this, pos)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP), OPTIONAL, INTENT(IN) :: pos(3)
        !
        CLASS(environ_system), INTENT(INOUT) :: this
        !
        INTEGER :: i, j
        INTEGER :: max_ntyp
        REAL(DP) :: dist
        REAL(DP) :: tot_weight
        !
        CHARACTER(LEN=80) :: routine = 'update_environ_system'
        !
        !--------------------------------------------------------------------------------
        !
        this%com = 0.D0
        this%width = 0.D0
        !
        max_ntyp = this%ntyp
        !
        IF (this%ntyp == 0) max_ntyp = this%ions%ntyp
        !
        IF (PRESENT(pos)) THEN
            this%com = pos ! fixed position (debugging with finite-differences)
        ELSE
            !
            !----------------------------------------------------------------------------
            ! Compute center of mass
            !
            tot_weight = 0.D0
            !
            DO i = 1, this%ions%number
                !
                IF (this%ions%ityp(i) > max_ntyp) CYCLE
                !
                this%com = this%com + this%ions%tau(:, i) * &
                           this%ions%iontype(this%ions%ityp(i))%weight
                !
                tot_weight = tot_weight + this%ions%iontype(this%ions%ityp(i))%weight
            END DO
            !
            this%com = this%com / tot_weight
        END IF
        !
        this%width = 0.D0
        !
        DO i = 1, this%ions%number
            !
            IF (this%ions%ityp(i) > max_ntyp) CYCLE
            !
            dist = 0.D0
            !
            DO j = 1, 3
                !
                IF ((this%dim == 1 .AND. j == this%axis) .OR. &
                    (this%dim == 2 .AND. j /= this%axis)) CYCLE
                !
                dist = dist + (this%ions%tau(j, i) - this%com(j))**2
            END DO
            !
            ! need to modify it into a smooth maximum to compute derivatives
            this%width = MAX(this%width, dist)
        END DO
        !
        this%width = SQRT(this%width)
        !
        !--------------------------------------------------------------------------------
        ! Output current state
        !
        CALL this%printout()
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE update_environ_system
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE destroy_environ_system(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_system), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: routine = 'destroy_environ_system'
        !
        !--------------------------------------------------------------------------------
        !
        IF (.NOT. ASSOCIATED(this%ions)) CALL io%destroy_error(routine)
        !
        !--------------------------------------------------------------------------------
        !
        NULLIFY (this%ions)
        !
        this%initialized = .FALSE.
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE destroy_environ_system
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !                                   OUTPUT METHODS
    !
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !>
    !! Prints the details of the system
    !!
    !! @param verbose       : (INTEGER) adds verbosity to global verbose
    !! @param debug_verbose : (INTEGER) replaces global verbose for debugging
    !! @param unit          : (INTEGER) output target (default = io%debug_unit)
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE print_environ_system(this, verbose, debug_verbose, unit)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_system), INTENT(IN) :: this
        INTEGER, OPTIONAL, INTENT(IN) :: verbose, debug_verbose, unit
        !
        INTEGER :: base_verbose, local_verbose, local_unit
        !
        CHARACTER(LEN=80) :: routine = 'print_environ_system'
        !
        !--------------------------------------------------------------------------------
        !
        IF (.NOT. io%lnode) RETURN
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
        ELSE IF (io%verbosity > 0) THEN
            base_verbose = io%verbosity
            !
            IF (PRESENT(verbose)) THEN
                local_verbose = base_verbose + verbose
            ELSE
                local_verbose = base_verbose
            END IF
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
            WRITE (local_unit, 1000)
            !
            IF (this%ntyp == 0) THEN
                WRITE (local_unit, 1001)
            ELSE
                WRITE (local_unit, 1002) this%ntyp
            END IF
            !
            WRITE (local_unit, 1003) this%dim, this%axis
            WRITE (local_unit, 1004) this%com, this%width
            !
        END IF
        !
        FLUSH (local_unit)
        !
        !--------------------------------------------------------------------------------
        !
1000    FORMAT(/, 4('%'), " SYSTEM ", 68('%'))
!
1001    FORMAT(/, " system is built from all present ionic types")
!
1002    FORMAT(/, " system is built from the first ", I3, " ionic types")
        !
1003    FORMAT(/, " system defined dimension   = ", I14, /, &
                " system defined axis        = ", I14)
        !
1004    FORMAT(/, " system center of mass      = ", 3F14.7, /, &
                " system width               = ", F14.7)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE print_environ_system
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE class_system
!----------------------------------------------------------------------------------------
