!----------------------------------------------------------------------------------------
!>
!!
!----------------------------------------------------------------------------------------
MODULE utils_electrons
    !------------------------------------------------------------------------------------
    !
    USE modules_constants, ONLY: DP
    !
    USE physical_types, ONLY: environ_electrons
    USE cell_types, ONLY: environ_cell
    !
    USE utils_density, ONLY: create_environ_density, init_environ_density, &
                             destroy_environ_density
    !
    USE tools_math, ONLY: integrate_environ_density
    !
    !------------------------------------------------------------------------------------
    !
    PRIVATE
    !
    PUBLIC :: create_environ_electrons, init_environ_electrons_first, &
              init_environ_electrons_second, update_environ_electrons, &
              destroy_environ_electrons
    !
    !------------------------------------------------------------------------------------
CONTAINS
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE create_environ_electrons(electrons)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_electrons), INTENT(INOUT) :: electrons
        !
        CHARACTER(LEN=80) :: label = 'electrons'
        !
        !--------------------------------------------------------------------------------
        !
        electrons%update = .FALSE.
        electrons%number = 0
        electrons%charge = 0.D0
        !
        CALL create_environ_density(electrons%density, label)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE create_environ_electrons
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_environ_electrons_first(nelec, electrons)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: nelec
        !
        TYPE(environ_electrons), INTENT(INOUT) :: electrons
        !
        !--------------------------------------------------------------------------------
        !
        electrons%initialized = .FALSE.
        electrons%number = nelec
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_environ_electrons_first
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_environ_electrons_second(cell, electrons)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_cell), INTENT(IN) :: cell
        !
        TYPE(environ_electrons), INTENT(INOUT) :: electrons
        !
        !--------------------------------------------------------------------------------
        !
        CALL init_environ_density(cell, electrons%density)
        !
        electrons%initialized = .TRUE.
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_environ_electrons_second
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE update_environ_electrons(nnr, rho, electrons, nelec)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: nnr
        REAL(DP), INTENT(IN) :: rho(nnr)
        !
        REAL(DP), INTENT(IN), OPTIONAL :: nelec
        !
        TYPE(environ_electrons), INTENT(INOUT) :: electrons
        !
        REAL(DP), PARAMETER :: tol = 1.D-4
        REAL(DP) :: charge
        !
        CHARACTER(LEN=80) :: sub_name = 'update_environ_electrons'
        !
        !--------------------------------------------------------------------------------
        !
        ! check on dimensions
        IF (nnr /= electrons%density%cell%nnr) &
            CALL env_errore(sub_name, 'Mismatch in grid size', 1)
        !
        !--------------------------------------------------------------------------------
        ! Assign input density to electrons%density%of_r
        !
        electrons%density%of_r = rho ! assign input density to electrons%density%of_r
        !
        !--------------------------------------------------------------------------------
        ! Update integral of electronic density and, if provided, check
        ! against input value
        !
        electrons%charge = integrate_environ_density(electrons%density)
        electrons%number = NINT(electrons%charge)
        !
        IF (PRESENT(nelec)) THEN
            !
            IF (ABS(electrons%charge - nelec) > tol) &
                CALL env_errore(sub_name, 'Mismatch in integrated electronic charge', 1)
            !
        END IF
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE update_environ_electrons
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE destroy_environ_electrons(lflag, electrons)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        LOGICAL, INTENT(IN) :: lflag
        !
        TYPE(environ_electrons), INTENT(INOUT) :: electrons
        !
        !--------------------------------------------------------------------------------
        !
        IF (electrons%initialized) THEN
            !
            CALL destroy_environ_density(electrons%density)
            !
            electrons%charge = 0.D0
            electrons%initialized = .FALSE.
        END IF
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE destroy_environ_electrons
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE utils_electrons
!----------------------------------------------------------------------------------------
