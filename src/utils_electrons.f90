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
        !
        ! BACKWARD COMPATIBILITY
        ! Compatible with QE-6.0 QE-6.1.X QE-6.2.X QE-6.3.X
        ! electrons%nspin = 1
        ! Compatible with QE-6.4.X QE-GIT
        ! END BACKWARD COMPATIBILITY
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
    ! BACKWARD COMPATIBILITY
    ! Compatible with QE-6.0 QE-6.1.X QE-6.2.X QE-6.3.X
    ! SUBROUTINE init_environ_electrons_first(nelec, nspin, electrons)
    ! Compatible with QE-6.4.X QE-GIT
    SUBROUTINE init_environ_electrons_first(nelec, electrons)
        ! END BACKWARD COMPATIBILITY
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: nelec
        ! BACKWARD COMPATIBILITY
        ! Compatible with QE-6.0 QE-6.1.X QE-6.2.X QE-6.3.X
        ! INTEGER, INTENT(IN) :: nspin
        ! Compatible with QE-6.4.X QE-GIT
        !
        ! END BACKWARD COMPATIBILITY
        TYPE(environ_electrons), INTENT(INOUT) :: electrons
        !
        !--------------------------------------------------------------------------------
        !
        electrons%initialized = .FALSE.
        electrons%number = nelec
        !
        ! BACKWARD COMPATIBILITY
        ! Compatible with QE-6.0 QE-6.1.X QE-6.2.X QE-6.3.X
        ! electrons%nspin = nspin
        ! Compatible with QE-6.4.X QE-GIT
        ! END BACKWARD COMPATIBILITY
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
    ! BACKWARD COMPATIBILITY
    ! Compatible with QE-6.0 QE-6.1.X QE-6.2.X QE-6.3.X
    ! SUBROUTINE update_environ_electrons(nspin, nnr, rho, electrons, nelec)
    ! Compatible with QE-6.4.X QE-GIT
    SUBROUTINE update_environ_electrons(nnr, rho, electrons, nelec)
        ! END BACKWARD COMPATIBILITY
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: nnr
        !
        ! BACKWARD COMPATIBILITY
        ! Compatible with QE-6.0 QE-6.1.X QE-6.2.X QE-6.3.X
        ! INTEGER, INTENT(IN) :: nspin
        ! REAL(DP), INTENT(IN) :: rho(nnr, nspin)
        ! Compatible with QE-6.4.X QE-GIT
        REAL(DP), INTENT(IN) :: rho(nnr)
        ! END BACKWARD COMPATIBILITY
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
        ! Check on dimensions
        !
        ! BACKWARD COMPATIBILITY
        ! Compatible with QE-6.0 QE-6.1.X QE-6.2.X QE-6.3.X
        ! IF (nspin .NE. electrons%nspin) &
        !     CALL errore(sub_name, 'Missmatch in spin size', 1)
        ! Compatible with QE-6.4.X QE-GIT
        ! END BACKWARD COMPATIBILITY
        !
        IF (nnr .NE. electrons%density%cell%nnr) &
            CALL errore(sub_name, 'Missmatch in grid size', 1)
        !
        !--------------------------------------------------------------------------------
        ! Assign input density to electrons%density%of_r
        !
        ! BACKWARD COMPATIBILITY
        ! Compatible with QE-5.X QE-6.1.X QE-6.2.X QE-6.3.X
        ! electrons%density%of_r(:) = rho(:, 1)
        ! IF (electrons%nspin .EQ. 2) &
        !     electrons%density%of_r(:) = electrons%density%of_r(:) + rho(:, 2)
        ! Compatible with QE-6.4.X and QE-GIT
        electrons%density%of_r = rho
        ! END BACKWARD COMPATIBILITY
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
            IF (ABS(electrons%charge - nelec) .GT. tol) &
                CALL errore(sub_name, 'Missmatch in integrated electronic charge', 1)
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
        LOGICAL, INTENT(IN) :: lflag ! #TODO unused variable
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
