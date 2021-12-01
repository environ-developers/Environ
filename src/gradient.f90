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
MODULE class_gradient
    !------------------------------------------------------------------------------------
    !
    USE class_io, ONLY: io
    USE env_mp, ONLY: env_mp_sum
    !
    USE environ_param, ONLY: DP
    !
    USE class_cell
    USE class_density
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
    TYPE, PUBLIC :: environ_gradient
        !--------------------------------------------------------------------------------
        !
        LOGICAL :: lupdate = .FALSE. ! optionally have an associated logical status
        !
        CHARACTER(LEN=80) :: label = 'gradient'
        ! optionally have an associated label, used for printout and debugs
        !
        TYPE(environ_cell), POINTER :: cell => NULL()
        ! each quantity in real-space is associated with its definition domain
        !
        REAL(DP), ALLOCATABLE :: of_r(:, :)
        ! the quantity in real-space, local to each processor
        !
        TYPE(environ_density) :: modulus
        !
        !--------------------------------------------------------------------------------
    CONTAINS
        !--------------------------------------------------------------------------------
        !
        PROCEDURE, PRIVATE :: create => create_environ_gradient
        PROCEDURE :: init => init_environ_gradient
        PROCEDURE :: copy => copy_environ_gradient
        PROCEDURE :: update_modulus => update_gradient_modulus
        PROCEDURE :: destroy => destroy_environ_gradient
        !
        PROCEDURE :: scalar_product => scalar_product_environ_gradient
        PROCEDURE :: scalar_product_density => scalar_product_environ_gradient_density
        !
        PROCEDURE :: printout => print_environ_gradient
        !
        !--------------------------------------------------------------------------------
    END TYPE environ_gradient
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
    SUBROUTINE create_environ_gradient(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_gradient), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: sub_name = 'create_environ_gradient'
        !
        !--------------------------------------------------------------------------------
        !
        IF (ASSOCIATED(this%cell)) CALL io%create_error(sub_name)
        !
        IF (ALLOCATED(this%of_r)) CALL io%create_error(sub_name)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE create_environ_gradient
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_environ_gradient(this, cell, label)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_cell), TARGET, INTENT(IN) :: cell
        CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: label
        !
        CLASS(environ_gradient), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: modulus_label = 'gradient_modulus'
        !
        CHARACTER(LEN=80) :: sub_name = 'init_environ_gradient'
        !
        !--------------------------------------------------------------------------------
        !
        CALL this%create()
        !
        IF (PRESENT(label)) THEN
            this%label = label
            modulus_label = TRIM(ADJUSTL(label))//'_modulus'
        END IF
        !
        CALL this%modulus%init(cell, modulus_label)
        !
        this%cell => cell
        !
        ALLOCATE (this%of_r(3, this%cell%nnr))
        this%of_r = 0.D0
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_environ_gradient
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE copy_environ_gradient(this, copy)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_gradient), INTENT(IN) :: this
        !
        TYPE(environ_gradient), INTENT(OUT) :: copy
        !
        INTEGER :: n
        !
        CHARACTER(LEN=80) :: sub_name = 'copy_environ_gradient'
        !
        !--------------------------------------------------------------------------------
        !
        copy%cell => this%cell
        !
        copy%lupdate = this%lupdate
        copy%label = this%label
        !
        IF (ALLOCATED(this%of_r)) THEN
            n = SIZE(this%of_r, 2)
            !
            IF (ALLOCATED(copy%of_r)) DEALLOCATE (copy%of_r)
            !
            ALLOCATE (copy%of_r(3, n))
            copy%of_r = this%of_r
        END IF
        !
        CALL this%modulus%copy(copy%modulus)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE copy_environ_gradient
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE update_gradient_modulus(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_gradient), INTENT(INOUT) :: this
        !
        INTEGER, POINTER :: ir_end
        !
        !--------------------------------------------------------------------------------
        !
        ir_end => this%cell%ir_end
        !
        this%modulus%of_r(1:ir_end) = SQRT(this%of_r(1, 1:ir_end)**2 + &
                                           this%of_r(2, 1:ir_end)**2 + &
                                           this%of_r(3, 1:ir_end)**2)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE update_gradient_modulus
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE destroy_environ_gradient(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_gradient), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: sub_name = 'destroy_environ_gradient'
        !
        !--------------------------------------------------------------------------------
        !
        IF (.NOT. ASSOCIATED(this%cell)) CALL io%destroy_error(sub_name)
        !
        IF (.NOT. ALLOCATED(this%of_r)) CALL io%destroy_error(sub_name)
        !
        !--------------------------------------------------------------------------------
        !
        CALL this%modulus%destroy()
        !
        NULLIFY (this%cell)
        !
        DEALLOCATE (this%of_r)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE destroy_environ_gradient
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
    SUBROUTINE scalar_product_environ_gradient(this, gradB, dens)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_gradient), INTENT(IN) :: this, gradB
        !
        TYPE(environ_density), INTENT(INOUT) :: dens
        !
        INTEGER :: ir
        !
        CHARACTER(LEN=80) :: sub_name = 'scalar_product_environ_gradient'
        !
        !--------------------------------------------------------------------------------
        !
        dens%of_r = 0.D0
        !
        IF (.NOT. ASSOCIATED(this%cell, gradB%cell)) &
            CALL io%error(sub_name, 'Mismatch in domain of input gradients', 1)
        !
        IF (.NOT. ASSOCIATED(this%cell, dens%cell)) &
            CALL io%error(sub_name, 'Mismatch in domain of input and output', 1)
        !
        DO ir = 1, dens%cell%ir_end
            dens%of_r(ir) = SUM(this%of_r(:, ir) * gradB%of_r(:, ir))
        END DO
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE scalar_product_environ_gradient
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    FUNCTION scalar_product_environ_gradient_density(this, density) RESULT(res)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_gradient), INTENT(IN) :: this
        TYPE(environ_density), INTENT(IN) :: density
        !
        REAL(DP) :: res(3)
        !
        INTEGER, POINTER :: ir_end
        !
        INTEGER :: ipol
        REAL(DP) :: scalar_product
        !
        CHARACTER(LEN=80) :: sub_name = 'scalar_product_environ_gradient_density'
        !
        !--------------------------------------------------------------------------------
        !
        res = 0.D0
        !
        IF (.NOT. ASSOCIATED(this%cell, density%cell)) &
            CALL io%error(sub_name, 'Mismatch in domain of input vectors', 1)
        !
        ir_end => density%cell%ir_end
        !
        DO ipol = 1, 3
            !
            scalar_product = DOT_PRODUCT(this%of_r(ipol, 1:ir_end), &
                                         density%of_r(1:ir_end))
            !
            CALL env_mp_sum(scalar_product, density%cell%dfft%comm)
            !
            res(ipol) = scalar_product * density%cell%domega
        END DO
        !
        !--------------------------------------------------------------------------------
    END FUNCTION scalar_product_environ_gradient_density
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !                                   OUTPUT METHODS
    !
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !>
    !! Prints the details of the gradient
    !!
    !! Nested objects receive a decremented passed verbose to trigger block printing
    !! If called by a parent object, prints details in block format
    !!
    !! @param verbose       : (INTEGER) adds verbosity to global verbose
    !! @param debug_verbose : (INTEGER) replaces global verbose for debugging
    !! @param unit          : (INTEGER) output target (default = io%debug_unit)
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE print_environ_gradient(this, verbose, debug_verbose, unit)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_gradient), INTENT(IN) :: this
        INTEGER, INTENT(IN), OPTIONAL :: verbose, debug_verbose, unit
        !
        INTEGER :: base_verbose, local_verbose, passed_verbose, local_unit
        !
        TYPE(environ_density) :: dens
        !
        CHARACTER(LEN=80) :: sub_name = 'print_environ_gradient'
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
            END IF
            !
            ! #TODO ADD MAXVAL AND MINVAL
            !
            IF (local_verbose >= 3) CALL this%modulus%write_cube_no_ions()
            !
            IF (local_verbose >= 4) THEN
                !
                CALL dens%init(this%cell)
                !
                dens%label = TRIM(ADJUSTL(this%label))//'_x'
                dens%of_r = this%of_r(1, :)
                !
                CALL dens%printout(passed_verbose, debug_verbose, local_unit)
                !
                dens%label = TRIM(ADJUSTL(this%label))//'_y'
                dens%of_r = this%of_r(2, :)
                !
                CALL dens%printout(passed_verbose, debug_verbose, local_unit)
                !
                dens%label = TRIM(ADJUSTL(this%label))//'_z'
                dens%of_r = this%of_r(3, :)
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
1000    FORMAT(/, 4('%'), ' GRADIENT ', 66('%'))
1001    FORMAT(/, ' GRADIENT', /, '========')
        !
1002    FORMAT(/, ' gradient label             = ', A50)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE print_environ_gradient
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE class_gradient
!----------------------------------------------------------------------------------------
