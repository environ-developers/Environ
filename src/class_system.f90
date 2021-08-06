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
MODULE class_system
    !------------------------------------------------------------------------------------
    !
    USE env_base_io, ONLY: ionode, environ_unit, verbose, depth
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
        LOGICAL :: lupdate = .FALSE.
        INTEGER :: ntyp
        INTEGER :: dim
        INTEGER :: axis
        REAL(DP) :: pos(3)
        REAL(DP) :: width
        !
        TYPE(environ_ions), POINTER :: ions
        !
        !--------------------------------------------------------------------------------
    CONTAINS
        !--------------------------------------------------------------------------------
        !
        PROCEDURE :: create => create_environ_system
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
        !--------------------------------------------------------------------------------
        !
        this%lupdate = .FALSE.
        this%ntyp = 0
        this%dim = 0
        this%axis = 1
        this%pos = 0.D0
        this%width = 0.D0
        NULLIFY (this%ions)
        !
        RETURN
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
        CHARACTER(LEN=80) :: sub_name = 'init_environ_system'
        !
        !--------------------------------------------------------------------------------
        !
        this%ntyp = ntyp
        this%dim = dim
        this%axis = axis
        !
        this%ions => ions
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_environ_system
    !------------------------------------------------------------------------------------
    !>
    !! Given the system definition compute position (centre of charge)
    !! and width (maximum distance from centre) of the system.
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE update_environ_system(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_system), INTENT(INOUT) :: this
        !
        INTEGER :: i, icor, max_ntyp
        REAL(DP) :: charge, dist
        INTEGER, POINTER :: ityp
        REAL(DP), POINTER :: zv
        !
        CHARACTER(LEN=80) :: sub_name = 'update_environ_system'
        !
        !--------------------------------------------------------------------------------
        !
        IF (.NOT. ASSOCIATED(this%ions)) &
            CALL env_errore(sub_name, 'Trying to use a non associated object', 1)
        !
        this%pos = 0.D0
        this%width = 0.D0
        !
        max_ntyp = this%ntyp
        !
        IF (this%ntyp == 0) max_ntyp = this%ions%ntyp
        !
        charge = 0.D0
        !
        DO i = 1, this%ions%number
            ityp => this%ions%ityp(i)
            !
            IF (ityp > max_ntyp) CYCLE
            !
            zv => this%ions%iontype(ityp)%zv
            charge = charge + zv
            this%pos(:) = this%pos(:) + this%ions%tau(:, i) * zv
        END DO
        !
        IF (ABS(charge) < 1.D-8) CALL env_errore(sub_name, 'System charge is zero', 1)
        !
        this%pos(:) = this%pos(:) / charge
        !
        this%width = 0.D0
        !
        DO i = 1, this%ions%number
            ityp => this%ions%ityp(i)
            !
            IF (ityp > max_ntyp) CYCLE
            !
            dist = 0.D0
            !
            DO icor = 1, 3
                !
                IF ((this%dim == 1 .AND. icor == this%axis) .OR. &
                    (this%dim == 2 .AND. icor /= this%axis)) CYCLE
                !
                dist = dist + (this%ions%tau(icor, i) - this%pos(icor))**2
            END DO
            !
            ! need to modify it into a smooth maximum to compute derivatives
            this%width = MAX(this%width, dist)
        END DO
        !
        this%width = SQRT(this%width)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE update_environ_system
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE destroy_environ_system(this, lflag)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        LOGICAL, INTENT(IN) :: lflag
        !
        CLASS(environ_system), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: sub_name = 'destroy_environ_system'
        !
        !--------------------------------------------------------------------------------
        !
        IF (lflag) THEN
            !
            IF (.NOT. ASSOCIATED(this%ions)) &
                CALL env_errore(sub_name, 'Trying to destroy a non associated object', 1)
            !
            NULLIFY (this%ions)
        END IF
        !
        RETURN
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
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE print_environ_system(this, local_verbose, local_depth)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_system), INTENT(IN) :: this
        INTEGER, INTENT(IN), OPTIONAL :: local_verbose
        INTEGER, INTENT(IN), OPTIONAL :: local_depth
        !
        INTEGER :: verbosity, passed_verbosity, passed_depth
        !
        CHARACTER(LEN=80) :: sub_name = 'print_environ_system'
        !
        !--------------------------------------------------------------------------------
        !
        IF (.NOT. ionode .OR. verbose == 0) RETURN
        !
        IF (PRESENT(local_verbose)) THEN
            verbosity = verbose + local_verbose
        ELSE
            verbosity = verbose
        END IF
        !
        IF (verbosity == 0) RETURN
        !
        IF (verbosity >= 1) THEN
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
            IF (this%ntyp == 0) THEN
                WRITE (environ_unit, 1002)
            ELSE
                WRITE (environ_unit, 1003) this%ntyp
            END IF
            !
            WRITE (environ_unit, 1004) this%dim, this%axis
            WRITE (environ_unit, 1005) this%pos, this%width
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
1000    FORMAT(/, 4('%'), ' SYSTEM ', 68('%'))
1001    FORMAT(/, ' SYSTEM', /, ' ======')
!
1002    FORMAT(/, ' system is built from all present ionic types')
!
1003    FORMAT(/, ' system is built from the first ', I3, ' ionic types')
        !
1004    FORMAT(/, ' system defined dimension   = ', I1, /, &
                ' system defined axis        = ', I1)
        !
1005    FORMAT(/, ' system center              =', 3F10.7, /, &
                ' system width               = ', F9.7)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE print_environ_system
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE class_system
!----------------------------------------------------------------------------------------
