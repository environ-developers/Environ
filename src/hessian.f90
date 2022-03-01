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
    USE class_io, ONLY: io
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
        LOGICAL :: lupdate = .FALSE. ! optionally have an associated logical status
        !
        CHARACTER(LEN=80) :: label = 'hessian'
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
        IF (ASSOCIATED(this%cell)) CALL io%create_error(sub_name)
        !
        IF (ALLOCATED(this%of_r)) CALL io%create_error(sub_name)
        !
        !--------------------------------------------------------------------------------
        !
        this%lupdate = .FALSE.
        this%label = 'hessian'
        !
        NULLIFY (this%cell)
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
        CALL this%laplacian%init(cell, laplacian_label)
        !
        this%cell => cell
        !
        ALLOCATE (this%of_r(3, 3, this%cell%nnr))
        this%of_r = 0.D0
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
        IF (.NOT. ASSOCIATED(this%cell)) CALL io%destroy_error(sub_name)
        !
        IF (.NOT. ALLOCATED(this%of_r)) CALL io%destroy_error(sub_name)
        !
        !--------------------------------------------------------------------------------
        !
        CALL this%laplacian%destroy()
        !
        NULLIFY (this%cell)
        !
        DEALLOCATE (this%of_r)
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
            CALL io%error(sub_name, 'Mismatch in domain of input hessian/gradients', 1)
        !
        IF (.NOT. ASSOCIATED(gradin%cell, gradout%cell)) &
            CALL io%error(sub_name, 'Mismatch in domain of input and output', 1)
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
    !! Prints the details of the hessian
    !!
    !! Nested objects receive a decremented passed verbose to trigger block printing
    !! If called by a parent object, prints details in block format
    !!
    !! @param verbose       : (INTEGER) adds verbosity to global verbose
    !! @param debug_verbose : (INTEGER) replaces global verbose for debugging
    !! @param unit          : (INTEGER) output target (default = io%debug_unit)
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE print_environ_hessian(this, verbose, debug_verbose, unit)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_hessian), INTENT(IN) :: this
        INTEGER, INTENT(IN), OPTIONAL :: verbose, debug_verbose, unit
        !
        INTEGER :: base_verbose, local_verbose, passed_verbose, local_unit
        !
        REAL(DP) :: integral
        TYPE(environ_density) :: dens
        !
        CHARACTER(LEN=80) :: sub_name = 'print_environ_hessian'
        !
        !--------------------------------------------------------------------------------
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
            passed_verbose = verbose - 1
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
            passed_verbose = local_verbose - base_verbose - 1
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
            !
            IF (io%lnode) THEN
                !
                IF (local_verbose >= base_verbose) THEN ! header
                    WRITE (local_unit, 1000)
                ELSE
                    !
                    CALL io%block_divider(local_verbose, base_verbose, local_unit)
                    !
                    WRITE (local_unit, 1001)
                END IF
                !
                WRITE (local_unit, 1002) ADJUSTL(this%label)
                !
            END IF
            !
            ! #TODO ADD MAXVAL AND MINVAL
            !
            IF (local_verbose >= 3) CALL this%laplacian%write_cube_no_ions()
            !
            IF (local_verbose >= 4) THEN
                !
                CALL dens%init(this%cell)
                !
                dens%label = TRIM(ADJUSTL(this%label))//'_xx'
                dens%of_r(:) = this%of_r(1, 1, :)
                !
                CALL dens%printout(passed_verbose, debug_verbose, local_unit)
                !
                dens%label = TRIM(ADJUSTL(this%label))//'_xy'
                dens%of_r(:) = this%of_r(1, 2, :)
                !
                CALL dens%printout(passed_verbose, debug_verbose, local_unit)
                !
                dens%label = TRIM(ADJUSTL(this%label))//'_xz'
                dens%of_r(:) = this%of_r(1, 3, :)
                !
                CALL dens%printout(passed_verbose, debug_verbose, local_unit)
                !
                dens%label = TRIM(ADJUSTL(this%label))//'_yy'
                dens%of_r(:) = this%of_r(2, 2, :)
                !
                CALL dens%printout(passed_verbose, debug_verbose, local_unit)
                !
                dens%label = TRIM(ADJUSTL(this%label))//'_yz'
                dens%of_r(:) = this%of_r(2, 3, :)
                !
                CALL dens%printout(passed_verbose, debug_verbose, local_unit)
                !
                dens%label = TRIM(ADJUSTL(this%label))//'_zz'
                dens%of_r(:) = this%of_r(3, 3, :)
                !
                CALL dens%printout(passed_verbose, debug_verbose, local_unit)
                !
                CALL dens%destroy()
                !
            END IF
            !
            IF (local_verbose < base_verbose) &
                CALL io%block_divider(local_verbose, base_verbose, local_unit)
            !
        END IF
        !
        FLUSH (local_unit)
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
