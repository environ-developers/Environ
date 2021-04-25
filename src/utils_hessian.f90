!----------------------------------------------------------------------------------------
!>
!!
!----------------------------------------------------------------------------------------
MODULE utils_hessian
    !------------------------------------------------------------------------------------
    !
    USE cell_types, ONLY: environ_cell
    USE representation_types, ONLY: environ_density, environ_hessian
    !
    USE utils_density, ONLY: create_environ_density, init_environ_density, &
                             copy_environ_density, destroy_environ_density
    !
    !------------------------------------------------------------------------------------
CONTAINS
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE create_environ_hessian(hessian, label)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CHARACTER(LEN=80), OPTIONAL, INTENT(IN) :: label
        !
        TYPE(environ_hessian), INTENT(INOUT) :: hessian
        !
        CHARACTER(LEN=80) :: laplacian_label
        !
        CHARACTER(LEN=80) :: sub_name = 'destroy_environ_hessian'
        !
        !--------------------------------------------------------------------------------
        !
        IF (PRESENT(label)) THEN
            hessian%label = label
            laplacian_label = TRIM(ADJUSTL(label))//'_laplacian'
        ELSE
            hessian%label = 'hessian'
            laplacian_label = 'hessian_laplacian'
        END IF
        !
        NULLIFY (hessian%cell)
        !
        IF (ALLOCATED(hessian%of_r)) &
            CALL errore(sub_name, 'Trying to create an already allocated object', 1)
        !
        CALL create_environ_density(hessian%laplacian, laplacian_label)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE create_environ_hessian
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_environ_hessian(cell, hessian)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_cell), TARGET, INTENT(IN) :: cell
        !
        TYPE(environ_hessian), INTENT(INOUT) :: hessian
        !
        CHARACTER(LEN=80) :: sub_name = 'init_environ_hessian'
        !
        !--------------------------------------------------------------------------------
        !
        hessian%update = .FALSE.
        !
        IF (ASSOCIATED(hessian%cell)) &
            CALL errore(sub_name, 'Trying to associate an associated object', 1)
        !
        hessian%cell => cell
        !
        IF (ALLOCATED(hessian%of_r)) &
            CALL errore(sub_name, 'Trying to allocate an allocated object', 1)
        !
        ALLOCATE (hessian%of_r(3, 3, hessian%cell%nnr))
        hessian%of_r = 0.D0
        !
        CALL init_environ_density(cell, hessian%laplacian)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_environ_hessian
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE copy_environ_hessian(horiginal, hcopy)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_hessian), INTENT(IN) :: horiginal
        !
        TYPE(environ_hessian), INTENT(OUT) :: hcopy
        !
        INTEGER :: n
        !
        CHARACTER(LEN=80) :: sub_name = 'copy_environ_hessian'
        !
        !--------------------------------------------------------------------------------
        !
        IF (.NOT. ASSOCIATED(horiginal%cell)) &
            CALL errore(sub_name, 'Trying to copy a non associated object', 1)
        !
        hcopy%cell => horiginal%cell
        !
        hcopy%update = horiginal%update
        hcopy%label = horiginal%label
        !
        IF (ALLOCATED(horiginal%of_r)) THEN
            n = SIZE(horiginal%of_r, 3)
            !
            IF (ALLOCATED(hcopy%of_r)) DEALLOCATE (hcopy%of_r)
            !
            ALLOCATE (hcopy%of_r(3, 3, n))
            hcopy%of_r = horiginal%of_r
        END IF
        !
        CALL copy_environ_density(horiginal%laplacian, hcopy%laplacian)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE copy_environ_hessian
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE update_hessian_laplacian(hessian)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_hessian), INTENT(INOUT) :: hessian
        !
        INTEGER, POINTER :: ir_end
        !
        !--------------------------------------------------------------------------------
        !
        ir_end => hessian%cell%ir_end
        !
        hessian%laplacian%of_r(1:ir_end) = hessian%of_r(1, 1, 1:ir_end) + &
                                           hessian%of_r(2, 2, 1:ir_end) + &
                                           hessian%of_r(3, 3, 1:ir_end)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE update_hessian_laplacian
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE destroy_environ_hessian(hessian)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_hessian), INTENT(INOUT) :: hessian
        !
        CHARACTER(LEN=80) :: sub_name = 'destroy_environ_hessian'
        !
        hessian%update = .FALSE.
        !
        !--------------------------------------------------------------------------------
        !
        IF (.NOT. ASSOCIATED(hessian%cell)) &
            CALL errore(sub_name, 'Trying to destroy a non associated object', 1)
        !
        NULLIFY (hessian%cell)
        !
        IF (.NOT. ALLOCATED(hessian%of_r)) &
            CALL errore(sub_name, 'Trying to destroy a non allocated object', 1)
        !
        DEALLOCATE (hessian%of_r)
        !
        CALL destroy_environ_density(hessian%laplacian)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE destroy_environ_hessian
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE utils_hessian
!----------------------------------------------------------------------------------------
