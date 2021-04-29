!----------------------------------------------------------------------------------------
!>
!!
!----------------------------------------------------------------------------------------
MODULE utils_electrostatics
    !------------------------------------------------------------------------------------
    !
    USE modules_constants, ONLY: DP, tpi
    !
    USE electrostatic_types
    USE core_types, ONLY: fft_core, oned_analytic_core, core_container
    !
    !------------------------------------------------------------------------------------
CONTAINS
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_gradient_solver(lconjugate, tol, step_type, step, maxstep, &
                                    preconditioner, screening_type, screening, gradient)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        LOGICAL, INTENT(IN) :: lconjugate
        INTEGER, INTENT(IN) :: maxstep
        REAL(DP), INTENT(IN) :: tol, step, screening
        CHARACTER(LEN=80), INTENT(IN) :: step_type, preconditioner, screening_type
        !
        TYPE(gradient_solver), INTENT(INOUT) :: gradient
        !
        !--------------------------------------------------------------------------------
        !
        gradient%lconjugate = lconjugate
        gradient%tol = tol
        gradient%step_type = step_type
        gradient%step = step
        gradient%maxstep = maxstep
        gradient%preconditioner = preconditioner
        gradient%screening_type = screening_type
        gradient%screening = screening
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_gradient_solver
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_iterative_solver(tol, mix_type, mix, maxiter, ndiis, iterative)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: maxiter, ndiis
        REAL(DP), INTENT(IN) :: tol, mix
        CHARACTER(LEN=80), INTENT(IN) :: mix_type
        !
        TYPE(iterative_solver), INTENT(INOUT) :: iterative
        !
        !--------------------------------------------------------------------------------
        !
        iterative%tol = tol
        iterative%mix_type = mix_type
        iterative%mix = mix
        iterative%maxiter = maxiter
        iterative%ndiis = ndiis
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_iterative_solver
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_newton_solver(tol, maxiter, newton)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: maxiter
        REAL(DP), INTENT(IN) :: tol
        !
        TYPE(newton_solver), INTENT(INOUT) :: newton
        !
        !--------------------------------------------------------------------------------
        !
        newton%tol = tol
        newton%maxiter = maxiter
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_newton_solver
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE create_electrostatic_solver(solver)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(electrostatic_solver), INTENT(INOUT) :: solver
        !
        !--------------------------------------------------------------------------------
        !
        solver%type_ = 'default'
        solver%auxiliary = 'none'
        solver%use_direct = .FALSE.
        !
        solver%use_gradient = .FALSE.
        NULLIFY (solver%gradient)
        !
        solver%use_iterative = .FALSE.
        NULLIFY (solver%iterative)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE create_electrostatic_solver
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_electrostatic_solver(type_, solver, gradient, iterative, newton, &
                                         auxiliary)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CHARACTER(LEN=80), INTENT(IN) :: type_
        TYPE(gradient_solver), TARGET, INTENT(IN), OPTIONAL :: gradient
        TYPE(iterative_solver), TARGET, INTENT(IN), OPTIONAL :: iterative
        TYPE(newton_solver), TARGET, INTENT(IN), OPTIONAL :: newton
        CHARACTER(LEN=80), INTENT(IN), OPTIONAL :: auxiliary
        !
        TYPE(electrostatic_solver), INTENT(INOUT) :: solver
        !
        INTEGER :: number
        !
        CHARACTER(LEN=80) :: sub_name = 'init_electrostatic_solver'
        !
        !--------------------------------------------------------------------------------
        !
        solver%type_ = type_
        !
        IF (PRESENT(auxiliary)) solver%auxiliary = auxiliary
        !
        SELECT CASE (TRIM(ADJUSTL(solver%type_)))
        CASE ('direct', 'default')
            solver%use_direct = .TRUE.
        CASE ('cg', 'sd', 'gradient')
            !
            IF (.NOT. PRESENT(gradient)) &
                CALL errore(sub_name, 'Missing specified solver type', 1)
            !
            solver%use_gradient = .TRUE.
            solver%gradient => gradient
        CASE ('iterative')
            !
            IF (.NOT. PRESENT(iterative)) &
                CALL errore(sub_name, 'Missing specified solver type', 1)
            !
            solver%use_iterative = .TRUE.
            solver%iterative => iterative
        CASE ('newton')
            !
            IF (.NOT. PRESENT(newton)) &
                CALL errore(sub_name, 'Missing specified solver type', 1)
            !
            solver%use_newton = .TRUE.
            solver%newton => newton
        CASE DEFAULT
            CALL errore(sub_name, 'Unexpected option for electrostatic solver type', 1)
        END SELECT
        !
        !--------------------------------------------------------------------------------
        ! Double check that one and only one solver is specified
        !
        number = 0
        !
        IF (solver%use_direct) number = number + 1
        !
        IF (solver%use_gradient) number = number + 1
        !
        IF (solver%use_iterative) number = number + 1
        !
        IF (solver%use_newton) number = number + 1
        !
        IF (number /= 1) &
            CALL errore(sub_name, 'Too few or too many solvers are active', 1)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_electrostatic_solver
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE destroy_electrostatic_solver(lflag, solver)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        LOGICAL, INTENT(IN) :: lflag
        !
        TYPE(electrostatic_solver), INTENT(INOUT) :: solver
        !
        !--------------------------------------------------------------------------------
        !
        IF (lflag) THEN
            solver%use_direct = .FALSE.
            solver%use_gradient = .FALSE.
            NULLIFY (solver%gradient)
            solver%use_iterative = .FALSE.
            NULLIFY (solver%iterative)
        END IF
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE destroy_electrostatic_solver
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE add_correction(correction, core)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(core_container), TARGET, INTENT(IN) :: correction
        !
        TYPE(core_container), INTENT(INOUT) :: core
        !
        !--------------------------------------------------------------------------------
        !
        core%need_correction = .TRUE.
        core%correction => correction
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE add_correction
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE create_electrostatic_setup(setup)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(electrostatic_setup), INTENT(INOUT) :: setup
        !
        !--------------------------------------------------------------------------------
        !
        setup%problem = 'poisson'
        NULLIFY (setup%solver)
        NULLIFY (setup%core)
        !
        setup%nested_problem = .FALSE.
        NULLIFY (setup%inner)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE create_electrostatic_setup
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_electrostatic_setup(problem, solver, core, setup)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CHARACTER(LEN=80), INTENT(IN) :: problem
        TYPE(electrostatic_solver), TARGET, INTENT(IN) :: solver
        TYPE(core_container), TARGET, INTENT(IN) :: core
        !
        TYPE(electrostatic_setup), INTENT(INOUT) :: setup
        !
        CHARACTER(LEN=80) :: sub_name = 'init_electrostatic_setup'
        !
        !--------------------------------------------------------------------------------
        ! Sanity check on the global setup
        !
        setup%problem = problem
        !
        SELECT CASE (TRIM(ADJUSTL(setup%problem)))
        CASE ('poisson', 'default')
        CASE ('generalized', 'gpe')
            !
            IF (solver%use_direct) &
                CALL errore(sub_name, &
                            'Cannot use a direct solver for &
                            &the Generalized Poisson eq.', 1)
            !
        CASE ('linpb', 'linmodpb', 'linearized-pb')
            !
            IF (solver%use_direct .OR. solver%use_iterative) &
                CALL errore(sub_name, &
                            'Only gradient-based solver for &
                            &the linearized Poisson-Boltzmann eq.', 1)
            !
            IF (core%need_correction) THEN
                !
                IF (core%correction%type_ /= '1da') &
                    CALL errore(sub_name, &
                                'linearized-PB problem requires &
                                &parabolic pbc correction.', 1)
                !
            ELSE
                !
                CALL errore(sub_name, &
                            'linearized-PB problem requires &
                            &parabolic pbc correction.', 1)
                !
            END IF
            !
        CASE ('pb', 'modpb', 'poisson-boltzmann')
            !
            IF (solver%use_direct .OR. solver%use_gradient) &
                CALL errore(sub_name, &
                            'No direct or gradient-based solver for &
                            &the full Poisson-Boltzmann eq.', 1)
            !
            IF (core%need_correction) THEN
                !
                IF (core%correction%type_ /= '1da') &
                    CALL errore(sub_name, &
                                'full-PB problem requires &
                                &parabolic pbc correction.', 1)
                !
            ELSE
                !
                CALL errore(sub_name, &
                            'full-PB problem requires parabolic pbc correction.', 1)
                !
            END IF
            !
        CASE DEFAULT
            CALL errore(sub_name, 'Unexpected keyword for electrostatic problem', 1)
        END SELECT
        !
        setup%solver => solver
        setup%core => core
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_electrostatic_setup
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE add_inner_setup(inner, outer)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(electrostatic_setup), TARGET, INTENT(IN) :: inner
        !
        TYPE(electrostatic_setup), INTENT(INOUT) :: outer
        !
        !--------------------------------------------------------------------------------
        !
        outer%nested_problem = .TRUE.
        outer%inner => inner
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE add_inner_setup
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE destroy_electrostatic_setup(lflag, setup)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        LOGICAL, INTENT(IN) :: lflag
        !
        TYPE(electrostatic_setup), INTENT(INOUT) :: setup
        !
        !--------------------------------------------------------------------------------
        !
        IF (lflag) THEN
            NULLIFY (setup%solver)
            NULLIFY (setup%core)
            setup%nested_problem = .FALSE.
            NULLIFY (setup%inner)
        END IF
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE destroy_electrostatic_setup
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE set_electrostatic_flags(setup, need_auxiliary, need_gradient, &
                                       need_factsqrt)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(electrostatic_setup), INTENT(IN) :: setup
        !
        LOGICAL, INTENT(INOUT) :: need_auxiliary, need_gradient, need_factsqrt
        !
        !--------------------------------------------------------------------------------
        !
        SELECT CASE (setup%problem)
        CASE ('generalized', 'linpb', 'linmodpb', 'pb', 'modpb')
            !
            IF (setup%solver%use_gradient) THEN
                !
                SELECT CASE (setup%solver%gradient%preconditioner)
                CASE ('sqrt')
                    need_factsqrt = .TRUE.
                CASE ('left', 'none')
                    need_gradient = .TRUE.
                END SELECT
                !
            END IF
            !
            IF (setup%solver%auxiliary /= 'none') need_auxiliary = .TRUE.
            !
        END SELECT
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE set_electrostatic_flags
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE utils_electrostatics
!----------------------------------------------------------------------------------------
