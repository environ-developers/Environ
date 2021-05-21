!----------------------------------------------------------------------------------------
!>
!!
!----------------------------------------------------------------------------------------
MODULE utils_density
    !------------------------------------------------------------------------------------
    !
    USE cell_types, ONLY: environ_cell
    USE representation_types, ONLY: environ_density
    !
    USE tools_math, ONLY: multipoles_environ_density
    !
    !------------------------------------------------------------------------------------
    !
    PRIVATE
    !
    PUBLIC :: create_environ_density, init_environ_density, copy_environ_density, &
              update_environ_density, destroy_environ_density
    !
    !------------------------------------------------------------------------------------
CONTAINS
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE create_environ_density(density, local_label)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CHARACTER(LEN=80), INTENT(IN), OPTIONAL :: local_label
        !
        TYPE(environ_density), INTENT(INOUT) :: density
        !
        CHARACTER(LEN=80) :: sub_name = 'create_environ_density'
        !
        CHARACTER(LEN=80) :: label = 'density'
        !
        !--------------------------------------------------------------------------------
        !
        IF (PRESENT(local_label)) THEN
            density%label = local_label
        ELSE
            density%label = label
        END IF
        !
        NULLIFY (density%cell)
        !
        IF (ALLOCATED(density%of_r)) &
            CALL env_errore(sub_name, 'Trying to create an already allocated object', 1)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE create_environ_density
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_environ_density(cell, density)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_cell), TARGET, INTENT(IN) :: cell
        !
        TYPE(environ_density), INTENT(INOUT) :: density
        !
        CHARACTER(LEN=80) :: sub_name = 'init_environ_density'
        !
        !--------------------------------------------------------------------------------
        !
        density%update = .FALSE.
        !
        IF (ASSOCIATED(density%cell)) &
            CALL env_errore(sub_name, 'Trying to associate an associated object', 1)
        !
        density%cell => cell
        !
        IF (ALLOCATED(density%of_r)) &
            CALL env_errore(sub_name, 'Trying to allocate an allocated object', 1)
        !
        ALLOCATE (density%of_r(density%cell%nnr))
        density%of_r = 0.D0
        !
        density%charge = 0.D0
        density%dipole = 0.D0
        density%quadrupole = 0.D0
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_environ_density
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE copy_environ_density(doriginal, dcopy)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_density), INTENT(IN) :: doriginal
        !
        TYPE(environ_density), INTENT(OUT) :: dcopy
        !
        CHARACTER(LEN=80) :: sub_name = 'copy_environ_density'
        !
        INTEGER :: n
        !
        !--------------------------------------------------------------------------------
        !
        IF (.NOT. ASSOCIATED(doriginal%cell)) &
            CALL env_errore(sub_name, 'Trying to copy a non associated object', 1)
        !
        dcopy%cell => doriginal%cell
        !
        dcopy%update = doriginal%update
        dcopy%label = doriginal%label
        dcopy%charge = doriginal%charge
        dcopy%dipole = doriginal%dipole
        dcopy%quadrupole = doriginal%quadrupole
        !
        IF (ALLOCATED(doriginal%of_r)) THEN
            n = SIZE(doriginal%of_r)
            !
            IF (ALLOCATED(dcopy%of_r)) DEALLOCATE (dcopy%of_r)
            !
            ALLOCATE (dcopy%of_r(n))
            dcopy%of_r = doriginal%of_r
        END IF
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE copy_environ_density
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE update_environ_density(density)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_density), INTENT(INOUT) :: density
        !
        !--------------------------------------------------------------------------------
        !
        CALL multipoles_environ_density(density, density%cell%origin, &
                                        density%charge, density%dipole, &
                                        density%quadrupole)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE update_environ_density
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE destroy_environ_density(density)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_density), INTENT(INOUT) :: density
        !
        CHARACTER(LEN=80) :: sub_name = 'destroy_environ_density'
        !
        !--------------------------------------------------------------------------------
        !
        density%update = .FALSE.
        !
        IF (.NOT. ASSOCIATED(density%cell)) &
            CALL env_errore(sub_name, 'Trying to destroy a non associated object', 1)
        !
        NULLIFY (density%cell)
        !
        IF (.NOT. ALLOCATED(density%of_r)) &
            CALL env_errore(sub_name, 'Trying to destroy a non allocated object', 1)
        !
        DEALLOCATE (density%of_r)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE destroy_environ_density
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE utils_density
!----------------------------------------------------------------------------------------
