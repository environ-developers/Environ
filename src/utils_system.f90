!----------------------------------------------------------------------------------------
!>
!!
!----------------------------------------------------------------------------------------
MODULE utils_system
    !------------------------------------------------------------------------------------
    !
    USE environ_param, ONLY: DP
    !
    USE types_physical, ONLY: environ_system, environ_ions
    !
    !------------------------------------------------------------------------------------
    !
    PRIVATE
    !
    PUBLIC :: create_environ_system, init_environ_system, update_environ_system, &
              destroy_environ_system
    !
    !------------------------------------------------------------------------------------
CONTAINS
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE create_environ_system(system)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_system), INTENT(INOUT) :: system
        !
        !--------------------------------------------------------------------------------
        !
        system%update = .FALSE.
        system%ntyp = 0
        system%dim = 0
        system%axis = 1
        system%pos = 0.D0
        system%width = 0.D0
        NULLIFY (system%ions)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE create_environ_system
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_environ_system(ntyp, dim, axis, ions, system)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: ntyp, dim, axis
        TYPE(environ_ions), TARGET, INTENT(IN) :: ions
        !
        TYPE(environ_system), INTENT(INOUT) :: system
        !
        CHARACTER(LEN=80) :: sub_name = 'init_environ_system'
        !
        !--------------------------------------------------------------------------------
        !
        system%ntyp = ntyp
        system%dim = dim
        system%axis = axis
        !
        system%ions => ions
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
    SUBROUTINE update_environ_system(system)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_system), INTENT(INOUT) :: system
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
        IF (.NOT. ASSOCIATED(system%ions)) &
            CALL env_errore(sub_name, 'Trying to use a non associated object', 1)
        !
        system%pos = 0.D0
        system%width = 0.D0
        !
        max_ntyp = system%ntyp
        !
        IF (system%ntyp == 0) max_ntyp = system%ions%ntyp
        !
        charge = 0.D0
        !
        DO i = 1, system%ions%number
            ityp => system%ions%ityp(i)
            !
            IF (ityp > max_ntyp) CYCLE
            !
            zv => system%ions%iontype(ityp)%zv
            charge = charge + zv
            system%pos(:) = system%pos(:) + system%ions%tau(:, i) * zv
        END DO
        !
        IF (ABS(charge) < 1.D-8) &
            CALL env_errore(sub_name, 'System charge is zero', 1)
        !
        system%pos(:) = system%pos(:) / charge
        !
        system%width = 0.D0
        !
        DO i = 1, system%ions%number
            ityp => system%ions%ityp(i)
            !
            IF (ityp > max_ntyp) CYCLE
            !
            dist = 0.D0
            !
            DO icor = 1, 3
                !
                IF ((system%dim == 1 .AND. icor == system%axis) .OR. &
                    (system%dim == 2 .AND. icor /= system%axis)) CYCLE
                !
                dist = dist + (system%ions%tau(icor, i) - system%pos(icor))**2
            END DO
            !
            ! need to modify it into a smooth maximum to compute derivatives
            system%width = MAX(system%width, dist)
        END DO
        !
        system%width = SQRT(system%width)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE update_environ_system
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE destroy_environ_system(lflag, system)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        LOGICAL, INTENT(IN) :: lflag
        !
        TYPE(environ_system), INTENT(INOUT) :: system
        !
        CHARACTER(LEN=80) :: sub_name = 'destroy_environ_system'
        !
        !--------------------------------------------------------------------------------
        !
        IF (lflag) THEN
            !
            IF (.NOT. ASSOCIATED(system%ions)) &
                CALL env_errore(sub_name, 'Trying to destroy a non associated object', 1)
            !
            NULLIFY (system%ions)
        END IF
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE destroy_environ_system
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE utils_system
!----------------------------------------------------------------------------------------
