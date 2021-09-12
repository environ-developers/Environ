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
!          Quinn Campbell     (Sandia National Laboratories)
!          Ismaila Dabo       (DMSE, Penn State)
!          Edan Bainglass     (Department of Physics, UNT)
!
!----------------------------------------------------------------------------------------
!>
!!
!----------------------------------------------------------------------------------------
MODULE class_core_1da
    !------------------------------------------------------------------------------------
    !
    USE environ_param, ONLY: DP
    !
    USE class_cell
    !
    USE class_core_numerical
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
    TYPE, EXTENDS(numerical_core), PUBLIC :: core_1da
        !--------------------------------------------------------------------------------
        !
        INTEGER :: nnr = 0
        INTEGER :: dim = 0
        INTEGER :: pdim = 0
        INTEGER :: axis = 0
        !
        REAL(DP) :: size = 0.D0
        REAL(DP) :: origin(3) = 0.D0
        !
        REAL(DP), ALLOCATABLE :: x(:, :)
        !
        !--------------------------------------------------------------------------------
    CONTAINS
        !--------------------------------------------------------------------------------
        !
        PROCEDURE :: create => create_core_1da
        PROCEDURE :: init => init_core_1da
        PROCEDURE :: update_cell => update_core_1da_cell
        PROCEDURE :: update_origin => update_core_1da_origin
        PROCEDURE :: destroy => destroy_core_1da
        !
        !--------------------------------------------------------------------------------
    END TYPE core_1da
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
    SUBROUTINE create_core_1da(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(core_1da), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: sub_name = 'create_core_1da'
        !
        !--------------------------------------------------------------------------------
        !
        IF (ASSOCIATED(this%cell)) CALL env_create_error(sub_name)
        !
        IF (ALLOCATED(this%x)) CALL env_create_error(sub_name)
        !
        !--------------------------------------------------------------------------------
        !
        this%core_type = '1d-analytic'
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE create_core_1da
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_core_1da(this, dim, axis, cell)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: dim, axis
        TYPE(environ_cell), TARGET, INTENT(IN) :: cell
        !
        CLASS(core_1da), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: sub_name = 'init_core_1da'
        !
        !--------------------------------------------------------------------------------
        !
        CALL this%create()
        !
        IF (dim == 3 .OR. dim < 0) &
            CALL env_errore(sub_name, &
                            'Wrong dimensions for analytic one dimensional core', 1)
        !
        this%dim = dim
        this%pdim = 3 - dim
        !
        IF ((dim == 1 .OR. dim == 2) .AND. (axis > 3 .OR. axis < 1)) &
            CALL env_errore(sub_name, &
                            'Wrong choice of axis for analytic one dimensional core', 1)
        !
        this%axis = axis
        this%nnr = cell%nnr
        this%cell => cell
        ALLOCATE (this%x(this%pdim, this%nnr))
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_core_1da
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE update_core_1da_cell(this, cell)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_cell), TARGET, INTENT(IN) :: cell
        !
        CLASS(core_1da), INTENT(INOUT) :: this
        !
        !--------------------------------------------------------------------------------
        !
        this%cell => cell
        !
        IF (this%dim == 0) THEN
            this%size = cell%omega
        ELSE IF (this%dim == 1) THEN
            this%size = cell%omega / cell%at(this%axis, this%axis)
        ELSE IF (this%dim == 2) THEN
            this%size = cell%at(this%axis, this%axis)
        END IF
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE update_core_1da_cell
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE update_core_1da_origin(this, origin)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP), INTENT(IN) :: origin(3)
        !
        CLASS(core_1da), INTENT(INOUT) :: this
        !
        LOGICAL :: physical
        INTEGER :: ir
        REAL(DP) :: r(3), r2
        !
        CHARACTER(LEN=80) :: sub_name = 'update_core_1da_origin'
        !
        !--------------------------------------------------------------------------------
        !
        this%origin = origin
        !
        IF (this%dim == 0) THEN
            !
            DO ir = 1, this%cell%ir_end
                !
                CALL this%cell%get_min_distance(ir, 0, 0, origin, r, r2, physical)
                ! compute minimum distance using minimum image convention
                !
                IF (.NOT. physical) CYCLE
                !
                this%x(:, ir) = -r
            END DO
            !
        ELSE IF (this%dim == 1) THEN
            CALL env_errore(sub_name, 'Option not yet implemented', 1)
        ELSE IF (this%dim == 2) THEN
            !
            DO ir = 1, this%cell%ir_end
                !
                CALL this%cell%get_min_distance(ir, 0, 0, origin, r, r2, physical)
                ! compute minimum distance using minimum image convention
                !
                IF (.NOT. physical) CYCLE
                !
                this%x(1, ir) = -r(this%axis)
            END DO
            !
        END IF
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE update_core_1da_origin
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE destroy_core_1da(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(core_1da), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: sub_name = 'destroy_core_1da'
        !
        !--------------------------------------------------------------------------------
        !
        IF (.NOT. ALLOCATED(this%x)) CALL env_destroy_error(sub_name)
        !
        IF (.NOT. ASSOCIATED(this%cell)) CALL env_destroy_error(sub_name)
        !
        !--------------------------------------------------------------------------------
        !
        DEALLOCATE (this%x)
        !
        NULLIFY (this%cell)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE destroy_core_1da
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE class_core_1da
!----------------------------------------------------------------------------------------
