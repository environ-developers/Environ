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
MODULE class_hessian
    !------------------------------------------------------------------------------------
    !
    USE env_base_io, ONLY: ionode, environ_unit, verbose, depth
    !
    USE environ_param, ONLY: DP
    !
    USE class_cell
    USE class_density
    USE class_gradient
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
    TYPE, PUBLIC :: environ_hessian
        !--------------------------------------------------------------------------------
        !
        LOGICAL :: lupdate ! optionally have an associated logical status
        !
        CHARACTER(LEN=80) :: label
        ! optionally have an associated label, used for printout and debugs
        !
        TYPE(environ_cell), POINTER :: cell => NULL()
        ! each quantity in real-space is associated with its definition domain
        !
        REAL(DP), ALLOCATABLE :: of_r(:, :, :)
        ! the quantity in real-space, local to each processor
        !
        TYPE(environ_density) :: laplacian
        !
        !--------------------------------------------------------------------------------
    CONTAINS
        !--------------------------------------------------------------------------------
        !
        PROCEDURE, PRIVATE :: create => create_environ_hessian
        PROCEDURE :: init => init_environ_hessian
        PROCEDURE :: copy => copy_environ_hessian
        PROCEDURE :: update_laplacian => update_hessian_laplacian
        PROCEDURE :: destroy => destroy_environ_hessian
        !
        PROCEDURE :: scalar_product => scalar_product_environ_hessian
        !
        PROCEDURE :: printout => print_environ_hessian
        !
        !--------------------------------------------------------------------------------
    END TYPE environ_hessian
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
    SUBROUTINE create_environ_hessian(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_hessian), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: sub_name = 'create_environ_hessian'
        !
        !--------------------------------------------------------------------------------
        !
        IF (ASSOCIATED(this%cell)) &
            CALL env_errore(sub_name, 'Trying to create an existing object', 1)
        !
        IF (ALLOCATED(this%of_r)) &
            CALL env_errore(sub_name, 'Trying to create an existing object', 1)
        !
        !--------------------------------------------------------------------------------
        !
        NULLIFY (this%cell)
        !
        this%label = 'hessian'
        this%lupdate = .FALSE.
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE create_environ_hessian
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_environ_hessian(this, cell, label)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_cell), TARGET, INTENT(IN) :: cell
        CHARACTER(LEN=80), INTENT(IN), OPTIONAL :: label
        !
        CLASS(environ_hessian), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: laplacian_label = 'hessian_laplacian'
        !
        CHARACTER(LEN=80) :: sub_name = 'init_environ_hessian'
        !
        !--------------------------------------------------------------------------------
        !
        CALL this%create()
        !
        IF (PRESENT(label)) THEN
            this%label = label
            laplacian_label = TRIM(ADJUSTL(label))//'_laplacian'
        END IF
        !
        this%cell => cell
        !
        ALLOCATE (this%of_r(3, 3, this%cell%nnr))
        this%of_r = 0.D0
        !
        CALL this%laplacian%init(cell, laplacian_label)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_environ_hessian
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE copy_environ_hessian(this, copy)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_hessian), INTENT(IN) :: this
        !
        TYPE(environ_hessian), INTENT(OUT) :: copy
        !
        INTEGER :: n
        !
        CHARACTER(LEN=80) :: sub_name = 'copy_environ_hessian'
        !
        !--------------------------------------------------------------------------------
        !
        IF (.NOT. ASSOCIATED(this%cell)) &
            CALL env_errore(sub_name, 'Trying to copy a non associated object', 1)
        !
        copy%cell => this%cell
        copy%lupdate = this%lupdate
        copy%label = this%label
        !
        IF (ALLOCATED(this%of_r)) THEN
            n = SIZE(this%of_r, 3)
            !
            IF (ALLOCATED(copy%of_r)) DEALLOCATE (copy%of_r)
            !
            ALLOCATE (copy%of_r(3, 3, n))
            copy%of_r = this%of_r
        END IF
        !
        CALL this%laplacian%copy(copy%laplacian)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE copy_environ_hessian
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE update_hessian_laplacian(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_hessian), INTENT(INOUT) :: this
        !
        INTEGER, POINTER :: ir_end
        !
        !--------------------------------------------------------------------------------
        !
        ir_end => this%cell%ir_end
        !
        this%laplacian%of_r(1:ir_end) = this%of_r(1, 1, 1:ir_end) + &
                                        this%of_r(2, 2, 1:ir_end) + &
                                        this%of_r(3, 3, 1:ir_end)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE update_hessian_laplacian
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE destroy_environ_hessian(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_hessian), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: sub_name = 'destroy_environ_hessian'
        !
        !--------------------------------------------------------------------------------
        !
        IF (.NOT. ASSOCIATED(this%cell)) &
            CALL env_errore(sub_name, 'Trying to destroy an empty object', 1)
        !
        IF (.NOT. ALLOCATED(this%of_r)) &
            CALL env_errore(sub_name, 'Trying to destroy an empty object', 1)
        !
        !--------------------------------------------------------------------------------
        !
        NULLIFY (this%cell)
        !
        DEALLOCATE (this%of_r)
        !
        CALL this%laplacian%destroy()
        !
        this%lupdate = .FALSE.
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE destroy_environ_hessian
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !                                  GENERAL METHODS
    !
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE scalar_product_environ_hessian(this, gradin, gradout)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_hessian), INTENT(IN) :: this
        TYPE(environ_gradient), INTENT(IN) :: gradin
        !
        TYPE(environ_gradient), INTENT(INOUT) :: gradout
        !
        INTEGER :: ir, ipol
        !
        CHARACTER(LEN=80) :: sub_name = 'scalar_product_environ_hessian'
        !
        !--------------------------------------------------------------------------------
        !
        gradout%of_r = 0.D0
        !
        IF (.NOT. ASSOCIATED(gradin%cell, this%cell)) &
            CALL env_errore(sub_name, 'Mismatch in domain of input hessian/gradients', 1)
        !
        IF (.NOT. ASSOCIATED(gradin%cell, gradout%cell)) &
            CALL env_errore(sub_name, 'Mismatch in domain of input and output', 1)
        !
        DO ir = 1, this%cell%ir_end
            !
            DO ipol = 1, 3
                gradout%of_r(ipol, ir) = SUM(this%of_r(:, ipol, ir) * gradin%of_r(:, ir))
            END DO
            !
        END DO
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE scalar_product_environ_hessian
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
    SUBROUTINE print_environ_hessian(this, local_verbose, local_depth)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_hessian), INTENT(IN) :: this
        INTEGER, INTENT(IN), OPTIONAL :: local_verbose
        INTEGER, INTENT(IN), OPTIONAL :: local_depth
        !
        TYPE(environ_cell), POINTER :: cell
        TYPE(environ_density) :: dens
        !
        INTEGER :: verbosity, passed_verbosity, passed_depth
        REAL(DP) :: integral
        !
        CHARACTER(LEN=80) :: sub_name = 'print_environ_hessian'
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
                WRITE (environ_unit, 1002) ADJUSTL(this%label)
                !
            END IF
            !
            ! #TODO ADD MAXVAL AND MINVAL
            !
            IF (verbosity >= 3) CALL this%laplacian%write_cube_no_ions()
            !
            IF (verbosity >= 4) THEN
                cell => this%cell
                !
                CALL dens%init(cell)
                !
                dens%label = TRIM(ADJUSTL(this%label))//'_xx'
                dens%of_r(:) = this%of_r(1, 1, :)
                !
                CALL dens%printout(passed_verbosity, passed_depth)
                !
                dens%label = TRIM(ADJUSTL(this%label))//'_xy'
                dens%of_r(:) = this%of_r(1, 2, :)
                !
                CALL dens%printout(passed_verbosity, passed_depth)
                !
                dens%label = TRIM(ADJUSTL(this%label))//'_xz'
                dens%of_r(:) = this%of_r(1, 3, :)
                !
                CALL dens%printout(passed_verbosity, passed_depth)
                !
                dens%label = TRIM(ADJUSTL(this%label))//'_yy'
                dens%of_r(:) = this%of_r(2, 2, :)
                !
                CALL dens%printout(passed_verbosity, passed_depth)
                !
                dens%label = TRIM(ADJUSTL(this%label))//'_yz'
                dens%of_r(:) = this%of_r(2, 3, :)
                !
                CALL dens%printout(passed_verbosity, passed_depth)
                !
                dens%label = TRIM(ADJUSTL(this%label))//'_zz'
                dens%of_r(:) = this%of_r(3, 3, :)
                !
                CALL dens%printout(passed_verbosity, passed_depth)
                !
                CALL dens%destroy()
                !
            END IF
            !
            IF (verbosity < verbose) CALL env_block_divider(verbosity)
            !
        END IF
        !
        FLUSH (environ_unit)
        !
        !--------------------------------------------------------------------------------
        !
1000    FORMAT(/, 4('%'), ' HESSIAN ', 67('%'))
1001    FORMAT(/, ' HESSIAN', /, '=======')
        !
1002    FORMAT(/, ' hessian label              = ', A50)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE print_environ_hessian
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE class_hessian
!----------------------------------------------------------------------------------------
