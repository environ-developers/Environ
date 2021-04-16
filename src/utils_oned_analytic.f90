!----------------------------------------------------------------------------------------
!>
!!
!----------------------------------------------------------------------------------------
MODULE utils_oned_analytic
    !------------------------------------------------------------------------------------
    !
    USE core_types
    !
CONTAINS
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE create_oned_analytic_core(oned_analytic)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(oned_analytic_core), INTENT(INOUT) :: oned_analytic
        !
        CHARACTER(LEN=80) :: sub_name = 'create_oned_analytic_core'
        !
        !--------------------------------------------------------------------------------
        !
        NULLIFY (oned_analytic%cell)
        oned_analytic%initialized = .FALSE.
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE create_oned_analytic_core
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_oned_analytic_core_first(dim, axis, oned_analytic)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: dim, axis
        !
        TYPE(oned_analytic_core), INTENT(INOUT) :: oned_analytic
        !
        CHARACTER(LEN=80) :: sub_name = 'init_oned_analytic_core_first'
        !
        !--------------------------------------------------------------------------------
        !
        IF (dim == 3 .OR. dim < 0) &
            CALL errore(sub_name, &
                        'Wrong dimensions for analytic one dimensional core', 1)
        !
        oned_analytic%d = dim
        oned_analytic%p = 3 - dim
        !
        IF ((dim == 1 .OR. dim == 2) .AND. (axis > 3 .OR. axis < 1)) &
            CALL errore(sub_name, &
                        'Wrong choice of axis for analytic one dimensional core', 1)
        !
        oned_analytic%axis = axis
        !
        oned_analytic%initialized = .FALSE.
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_oned_analytic_core_first
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_oned_analytic_core_second(cell, oned_analytic)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_cell), INTENT(IN) :: cell
        !
        TYPE(oned_analytic_core), INTENT(INOUT) :: oned_analytic
        !
        CHARACTER(LEN=80) :: sub_name = 'init_oned_analytic_core_second'
        !
        !--------------------------------------------------------------------------------
        !
        CALL update_oned_analytic_core_cell(cell, oned_analytic)
        !
        oned_analytic%n = cell%nnr
        ALLOCATE (oned_analytic%x(oned_analytic%p, oned_analytic%n))
        !
        CALL update_oned_analytic_core_origin(cell%origin, oned_analytic)
        !
        oned_analytic%initialized = .TRUE.
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_oned_analytic_core_second
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE update_oned_analytic_core_cell(cell, oned_analytic)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_cell), TARGET, INTENT(IN) :: cell
        !
        TYPE(oned_analytic_core), INTENT(INOUT) :: oned_analytic
        !
        !--------------------------------------------------------------------------------
        !
        oned_analytic%cell => cell
        !
        IF (oned_analytic%d == 0) THEN
            oned_analytic%size = cell%omega
        ELSE IF (oned_analytic%d == 1) THEN
            !
            oned_analytic%size = &
                cell%omega / cell%at(oned_analytic%axis, oned_analytic%axis) / cell%alat
            !
        ELSE IF (oned_analytic%d == 2) THEN
            !
            oned_analytic%size = &
                cell%at(oned_analytic%axis, oned_analytic%axis) * cell%alat
            !
        END IF
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE update_oned_analytic_core_cell
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE update_oned_analytic_core_origin(origin, oned_analytic)
        !--------------------------------------------------------------------------------
        !
        USE tools_generate_functions, ONLY: generate_axis, generate_distance
        !
        IMPLICIT NONE
        !
        REAL(DP), INTENT(IN) :: origin(3)
        TYPE(oned_analytic_core), INTENT(INOUT) :: oned_analytic
        !
        CHARACTER(LEN=80) :: sub_name = 'update_oned_analytic_core_origin'
        !
        !--------------------------------------------------------------------------------
        !
        oned_analytic%origin = origin
        !
        IF (oned_analytic%d == 0) THEN
            !
            CALL generate_distance(oned_analytic%cell, oned_analytic%origin, &
                                   oned_analytic%x)
            !
        ELSE IF (oned_analytic%d == 1) THEN
            CALL errore(sub_name, 'Option not yet implemented', 1)
        ELSE IF (oned_analytic%d == 2) THEN
            !
            CALL generate_axis(oned_analytic%cell, oned_analytic%axis, &
                               oned_analytic%origin, oned_analytic%x(1, :))
            !
        END IF
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE update_oned_analytic_core_origin
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE destroy_oned_analytic_core(lflag, oned_analytic)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        LOGICAL, INTENT(IN) :: lflag
        !
        TYPE(oned_analytic_core), INTENT(INOUT) :: oned_analytic
        !
        CHARACTER(LEN=80) :: sub_name = 'destroy_oned_analytic_core'
        !
        !--------------------------------------------------------------------------------
        !
        IF (oned_analytic%initialized) THEN
            !
            IF (.NOT. ALLOCATED(oned_analytic%x)) &
                CALL errore(sub_name, 'Trying to destroy a non-allocated component', 1)
            !
            DEALLOCATE (oned_analytic%x)
            !
            IF (.NOT. ASSOCIATED(oned_analytic%cell)) &
                CALL errore(sub_name, 'Trying to nullify a non-associated pointer', 1)
            !
            NULLIFY (oned_analytic%cell)
            oned_analytic%initialized = .FALSE.
        END IF
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE destroy_oned_analytic_core
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE utils_oned_analytic
!----------------------------------------------------------------------------------------
