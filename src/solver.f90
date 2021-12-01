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
! Authors: Francesco Nattino  (THEOS and NCCR-MARVEL, EPFL)
!          Oliviero Andreussi (Department of Physics, UNT)
!          Nicola Marzari     (THEOS and NCCR-MARVEL, EPFL)
!          Edan Bainglass     (Department of Physics, UNT)
!
!----------------------------------------------------------------------------------------
!>
!!
!----------------------------------------------------------------------------------------
MODULE class_solver
    !------------------------------------------------------------------------------------
    !
    USE class_io, ONLY: io
    !
    USE class_density
    USE class_gradient
    !
    USE class_core_container
    !
    USE class_charges
    USE class_dielectric
    USE class_electrolyte
    USE class_semiconductor
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
    TYPE, ABSTRACT, PUBLIC :: electrostatic_solver
        !--------------------------------------------------------------------------------
        !
        CHARACTER(LEN=80) :: solver_type
        !
        TYPE(core_container), POINTER :: cores => NULL()
        !
        !--------------------------------------------------------------------------------
    CONTAINS
        !--------------------------------------------------------------------------------
        ! Admin
        !
        PROCEDURE, PRIVATE :: create => create_solver
        PROCEDURE :: set_cores => set_solver_cores
        PROCEDURE :: destroy_cores => destroy_solver_cores
        !
        PROCEDURE(destroy_solver), DEFERRED :: destroy
        !
        !--------------------------------------------------------------------------------
        ! Direct
        !
        PROCEDURE :: poisson_charges
        PROCEDURE :: poisson_density
        !
        GENERIC :: poisson => &
            poisson_charges, &
            poisson_density
        !
        PROCEDURE :: grad_poisson_charges
        PROCEDURE :: grad_poisson_density
        !
        GENERIC :: grad_poisson => &
            grad_poisson_charges, &
            grad_poisson_density
        !
        !--------------------------------------------------------------------------------
        ! Iterative
        !
        PROCEDURE :: generalized_charges
        PROCEDURE :: generalized_density
        !
        GENERIC :: generalized => &
            generalized_charges, &
            generalized_density
        !
        PROCEDURE :: linearized_pb_charges
        PROCEDURE :: linearized_pb_density
        !
        GENERIC :: linearized_pb => &
            linearized_pb_charges, &
            linearized_pb_density
        !
        PROCEDURE :: pb_nested_charges
        PROCEDURE :: pb_nested_density
        !
        GENERIC :: pb_nested => &
            pb_nested_charges, &
            pb_nested_density
        !
        !--------------------------------------------------------------------------------
    END TYPE electrostatic_solver
    !------------------------------------------------------------------------------------
    !
    ABSTRACT INTERFACE
        SUBROUTINE destroy_solver(this)
            IMPORT electrostatic_solver
            CLASS(electrostatic_solver), INTENT(INOUT) :: this
        END SUBROUTINE
    END INTERFACE
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
    SUBROUTINE create_solver(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(electrostatic_solver), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: sub_name = 'create_solver'
        !
        !--------------------------------------------------------------------------------
        !
        IF (ASSOCIATED(this%cores)) CALL io%create_error(sub_name)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE create_solver
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE set_solver_cores(this, cores)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(core_container), TARGET, INTENT(IN) :: cores
        !
        CLASS(electrostatic_solver), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: sub_name = 'set_solver_cores'
        !
        !--------------------------------------------------------------------------------
        !
        CALL this%create()
        !
        this%cores => cores
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE set_solver_cores
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE destroy_solver_cores(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(electrostatic_solver), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: sub_name = 'destroy_solver_cores'
        !
        !--------------------------------------------------------------------------------
        !
        IF (.NOT. ASSOCIATED(this%cores)) CALL io%destroy_error(sub_name)
        !
        !--------------------------------------------------------------------------------
        !
        CALL this%cores%destroy()
        !
        NULLIFY (this%cores)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE destroy_solver_cores
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !                               DIRECT SOLVER METHODS
    !
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE poisson_charges(this, charges, v)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(electrostatic_solver), INTENT(IN) :: this
        TYPE(environ_charges), INTENT(IN) :: charges
        !
        TYPE(environ_density), INTENT(INOUT) :: v
        !
        CHARACTER(LEN=80) :: sub_name = 'poisson_charges'
        !
        !--------------------------------------------------------------------------------
        !
        CALL io%error(sub_name, "Not implemented", 1)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE poisson_charges
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE poisson_density(this, charges, v, electrolyte, semiconductor)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(electrostatic_solver), INTENT(IN) :: this
        TYPE(environ_density), INTENT(IN) :: charges
        TYPE(environ_electrolyte), OPTIONAL, INTENT(IN) :: electrolyte
        TYPE(environ_semiconductor), OPTIONAL, INTENT(IN) :: semiconductor
        !
        TYPE(environ_density), INTENT(INOUT) :: v
        !
        TYPE(environ_density) :: local
        !
        CHARACTER(LEN=80) :: sub_name = 'poisson_density'
        !
        !--------------------------------------------------------------------------------
        !
        CALL io%error(sub_name, "Not implemented", 1)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE poisson_density
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE grad_poisson_charges(this, charges, grad_v)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(electrostatic_solver), INTENT(IN) :: this
        TYPE(environ_charges), INTENT(IN) :: charges
        !
        TYPE(environ_gradient), INTENT(INOUT) :: grad_v
        !
        CHARACTER(LEN=80) :: sub_name = 'grad_poisson_charges'
        !
        !--------------------------------------------------------------------------------
        !
        CALL io%error(sub_name, "Not implemented", 1)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE grad_poisson_charges
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE grad_poisson_density(this, charges, grad_v, electrolyte, semiconductor)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(electrostatic_solver), INTENT(IN) :: this
        TYPE(environ_density), INTENT(IN) :: charges
        TYPE(environ_electrolyte), OPTIONAL, INTENT(IN) :: electrolyte
        TYPE(environ_semiconductor), OPTIONAL, INTENT(IN) :: semiconductor
        !
        TYPE(environ_gradient), INTENT(INOUT) :: grad_v
        !
        CHARACTER(LEN=80) :: sub_name = 'grad_poisson_density'
        !
        !--------------------------------------------------------------------------------
        !
        CALL io%error(sub_name, "Not implemented", 1)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE grad_poisson_density
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !                              ITERATIVE SOLVER METHODS
    !
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE generalized_charges(this, charges, v)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(electrostatic_solver), INTENT(IN) :: this
        !
        TYPE(environ_density), INTENT(INOUT) :: v
        TYPE(environ_charges), INTENT(INOUT) :: charges
        !
        CHARACTER(LEN=80) :: sub_name = 'generalized_charges'
        !
        !--------------------------------------------------------------------------------
        !
        CALL io%error(sub_name, "Not implemented", 1)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE generalized_charges
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE generalized_density(this, charges, dielectric, v, electrolyte, &
                                   semiconductor)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(electrostatic_solver), INTENT(IN) :: this
        TYPE(environ_density), INTENT(IN) :: charges
        TYPE(environ_electrolyte), OPTIONAL, INTENT(IN) :: electrolyte
        TYPE(environ_semiconductor), OPTIONAL, INTENT(IN) :: semiconductor
        !
        TYPE(environ_density), INTENT(INOUT) :: v
        TYPE(environ_dielectric), INTENT(INOUT) :: dielectric
        !
        CHARACTER(LEN=80) :: sub_name = 'generalized_density'
        !
        !--------------------------------------------------------------------------------
        !
        CALL io%error(sub_name, "Not implemented", 1)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE generalized_density
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE linearized_pb_charges(this, charges, v, screening)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(electrostatic_solver), INTENT(IN) :: this
        TYPE(environ_density), OPTIONAL, INTENT(IN) :: screening
        !
        TYPE(environ_density), INTENT(INOUT) :: v
        TYPE(environ_charges), INTENT(INOUT) :: charges
        !
        TYPE(environ_density) :: local_screening
        !
        CHARACTER(LEN=80) :: sub_name = 'linearized_pb_charges'
        !
        !--------------------------------------------------------------------------------
        !
        CALL io%error(sub_name, "Not implemented", 1)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE linearized_pb_charges
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE linearized_pb_density(this, charges, electrolyte, v, dielectric, &
                                     screening)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(electrostatic_solver), INTENT(IN) :: this
        TYPE(environ_density), INTENT(IN) :: charges
        TYPE(environ_electrolyte), INTENT(IN) :: electrolyte
        TYPE(environ_density), OPTIONAL, INTENT(IN) :: screening
        !
        TYPE(environ_density), INTENT(INOUT) :: v
        TYPE(environ_dielectric), OPTIONAL, INTENT(INOUT) :: dielectric
        !
        TYPE(environ_density) :: local_screening
        !
        CHARACTER(LEN=80) :: sub_name = 'linearized_pb_density'
        !
        !--------------------------------------------------------------------------------
        !
        CALL io%error(sub_name, "Not implemented", 1)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE linearized_pb_density
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE pb_nested_charges(this, charges, v, inner)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(electrostatic_solver), INTENT(IN) :: this
        CLASS(electrostatic_solver), OPTIONAL, INTENT(IN) :: inner
        !
        TYPE(environ_density), INTENT(INOUT) :: v
        TYPE(environ_charges), INTENT(INOUT) :: charges
        !
        CHARACTER(LEN=80) :: sub_name = 'pb_nested_charges'
        !
        !--------------------------------------------------------------------------------
        !
        CALL io%error(sub_name, "Not implemented", 1)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE pb_nested_charges
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE pb_nested_density(this, charges, v, electrolyte, dielectric, inner)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(electrostatic_solver), INTENT(IN) :: this
        TYPE(environ_density), INTENT(IN) :: charges
        TYPE(environ_electrolyte), INTENT(IN) :: electrolyte
        CLASS(electrostatic_solver), OPTIONAL, INTENT(IN) :: inner
        !
        TYPE(environ_density), INTENT(INOUT) :: v
        TYPE(environ_dielectric), OPTIONAL, INTENT(INOUT) :: dielectric
        !
        CHARACTER(LEN=80) :: sub_name = 'pb_nested_density'
        !
        !--------------------------------------------------------------------------------
        !
        CALL io%error(sub_name, "Not implemented", 1)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE pb_nested_density
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE class_solver
!----------------------------------------------------------------------------------------
