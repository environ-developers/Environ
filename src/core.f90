!----------------------------------------------------------------------------------------
!
! Copyright (C) 2021 ENVIRON (www.quantum-environ.org)
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
! Authors: Edan Bainglass (Department of Physics, UNT)
!
!----------------------------------------------------------------------------------------
!>
!!
!----------------------------------------------------------------------------------------
MODULE class_core
    !------------------------------------------------------------------------------------
    !
    USE class_io, ONLY: io
    !
    USE environ_param, ONLY: DP
    !
    USE class_cell
    USE class_density
    USE class_gradient
    USE class_function
    USE class_functions
    USE class_hessian
    !
    USE class_electrolyte_base
    USE class_semiconductor_base
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
    TYPE, ABSTRACT, PUBLIC :: environ_core
        !--------------------------------------------------------------------------------
        !
        CHARACTER(LEN=80) :: core_type
        !
        TYPE(environ_cell), POINTER :: cell => NULL()
        !
        !--------------------------------------------------------------------------------
    CONTAINS
        !--------------------------------------------------------------------------------
        ! Admin
        !
        PROCEDURE(create_core), DEFERRED :: create
        PROCEDURE(update_core_cell), DEFERRED :: update_cell
        PROCEDURE(destroy_core), DEFERRED :: destroy
        !
        !--------------------------------------------------------------------------------
        ! Derivatives
        !
        PROCEDURE :: gradient => calc_gradient
        PROCEDURE :: divergence => calc_divergence
        PROCEDURE :: laplacian => calc_laplacian
        PROCEDURE :: hessian => calc_hessian
        !
        PROCEDURE :: calc_convolution_density
        PROCEDURE :: calc_convolution_gradient
        PROCEDURE :: calc_convolution_hessian
        !
        GENERIC :: convolution => &
            calc_convolution_density, &
            calc_convolution_gradient, &
            calc_convolution_hessian
        !
        !--------------------------------------------------------------------------------
        ! Electrostatics
        !
        PROCEDURE :: poisson => calc_poisson
        PROCEDURE :: grad_poisson => calc_grad_poisson
        PROCEDURE :: force => calc_force
        !
        PROCEDURE :: hess_v_h_of_rho_r
        PROCEDURE :: field_of_grad_rho
        !
        !--------------------------------------------------------------------------------
        ! Corrections
        !
        PROCEDURE :: calc_vperiodic
        PROCEDURE :: calc_grad_vperiodic
        PROCEDURE :: calc_fperiodic
        !
        PROCEDURE :: calc_vgcs
        PROCEDURE :: calc_grad_vgcs
        !
        PROCEDURE :: calc_vms
        PROCEDURE :: calc_grad_vms
        !
        GENERIC :: potential => calc_vperiodic, calc_vgcs, calc_vms
        GENERIC :: grad_potential => calc_grad_vperiodic, calc_grad_vgcs, calc_grad_vms
        GENERIC :: force_periodic => calc_fperiodic
        !
        !--------------------------------------------------------------------------------
    END TYPE environ_core
    !------------------------------------------------------------------------------------
    !
    ABSTRACT INTERFACE
        SUBROUTINE create_core(this)
            IMPORT environ_core
            CLASS(environ_core), INTENT(INOUT) :: this
        END SUBROUTINE
        SUBROUTINE update_core_cell(this, cell)
            IMPORT environ_core, environ_cell
            TYPE(environ_cell), TARGET, INTENT(IN) :: cell
            CLASS(environ_core), INTENT(INOUT) :: this
        END SUBROUTINE
        SUBROUTINE destroy_core(this)
            IMPORT environ_core
            CLASS(environ_core), INTENT(INOUT) :: this
        END SUBROUTINE
    END INTERFACE
    !
    !------------------------------------------------------------------------------------
CONTAINS
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !                                 DERIVATIVE METHODS
    !
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_gradient(this, f, grad)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_core), TARGET, INTENT(IN) :: this
        TYPE(environ_density), INTENT(IN) :: f
        !
        TYPE(environ_gradient), INTENT(INOUT) :: grad
        !
        CHARACTER(LEN=80) :: sub_name = 'calc_gradient'
        !
        !--------------------------------------------------------------------------------
        !
        CALL io%error(sub_name, "Not implemented", 1)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_gradient
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_divergence(this, grad, div)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_core), TARGET, INTENT(IN) :: this
        TYPE(environ_gradient), INTENT(IN) :: grad
        !
        TYPE(environ_density), INTENT(INOUT) :: div
        !
        CHARACTER(LEN=80) :: sub_name = 'calc_divergence'
        !
        !--------------------------------------------------------------------------------
        !
        CALL io%error(sub_name, "Not implemented", 1)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_divergence
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_laplacian(this, f, lapla)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_core), TARGET, INTENT(IN) :: this
        TYPE(environ_density), INTENT(IN) :: f
        !
        TYPE(environ_density), INTENT(INOUT) :: lapla
        !
        CHARACTER(LEN=80) :: sub_name = 'calc_laplacian'
        !
        !--------------------------------------------------------------------------------
        !
        CALL io%error(sub_name, "Not implemented", 1)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_laplacian
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_hessian(this, f, grad, hess)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_core), TARGET, INTENT(IN) :: this
        TYPE(environ_density), INTENT(IN) :: f
        !
        TYPE(environ_gradient), INTENT(INOUT) :: grad
        TYPE(environ_hessian), INTENT(INOUT) :: hess
        !
        CHARACTER(LEN=80) :: sub_name = 'calc_hessian'
        !
        !--------------------------------------------------------------------------------
        !
        CALL io%error(sub_name, "Not implemented", 1)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_hessian
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_convolution_density(this, f1, f2, f_out)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_core), TARGET, INTENT(IN) :: this
        TYPE(environ_density), INTENT(IN) :: f1
        TYPE(environ_density), INTENT(IN) :: f2
        !
        TYPE(environ_density), INTENT(INOUT) :: f_out
        !
        CHARACTER(LEN=80) :: sub_name = 'calc_convolution_density'
        !
        !--------------------------------------------------------------------------------
        !
        CALL io%error(sub_name, "Not implemented", 1)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_convolution_density
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_convolution_gradient(this, f, grad, grad_out)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_core), TARGET, INTENT(IN) :: this
        TYPE(environ_density), INTENT(IN) :: f
        TYPE(environ_gradient), INTENT(IN) :: grad
        !
        TYPE(environ_gradient), INTENT(INOUT) :: grad_out
        !
        CHARACTER(LEN=80) :: sub_name = 'calc_convolution_gradient'
        !
        !--------------------------------------------------------------------------------
        !
        CALL io%error(sub_name, "Not implemented", 1)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_convolution_gradient
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_convolution_hessian(this, f, hess, hess_out)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_core), TARGET, INTENT(IN) :: this
        TYPE(environ_density), INTENT(IN) :: f
        TYPE(environ_hessian), INTENT(IN) :: hess
        !
        TYPE(environ_hessian), INTENT(INOUT) :: hess_out
        !
        CHARACTER(LEN=80) :: sub_name = 'calc_convolution_hessian'
        !
        !--------------------------------------------------------------------------------
        !
        CALL io%error(sub_name, "Not implemented", 1)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_convolution_hessian
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !                               ELECTROSTATIC METHODS
    !
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_poisson(this, rho, v)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_core), TARGET, INTENT(IN) :: this
        TYPE(environ_density), INTENT(IN) :: rho
        !
        TYPE(environ_density), INTENT(INOUT) :: v
        !
        CHARACTER(LEN=80) :: sub_name = 'calc_poisson'
        !
        !--------------------------------------------------------------------------------
        !
        CALL io%error(sub_name, "Not implemented", 1)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_poisson
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_grad_poisson(this, rho, grad_v)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_core), TARGET, INTENT(IN) :: this
        TYPE(environ_density), INTENT(IN) :: rho
        !
        TYPE(environ_gradient), INTENT(INOUT) :: grad_v
        !
        CHARACTER(LEN=80) :: sub_name = 'calc_grad_poisson'
        !
        !--------------------------------------------------------------------------------
        !
        CALL io%error(sub_name, "Not implemented", 1)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_grad_poisson
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_force(this, nat, rho, ions, force)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_core), TARGET, INTENT(IN) :: this
        INTEGER, INTENT(IN) :: nat
        TYPE(environ_density), INTENT(IN) :: rho
        TYPE(environ_functions), TARGET, INTENT(IN) :: ions
        !
        REAL(DP), INTENT(INOUT) :: force(3, nat)
        !
        CHARACTER(LEN=80) :: sub_name = 'calc_force'
        !
        !--------------------------------------------------------------------------------
        !
        CALL io%error(sub_name, "Not implemented", 1)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_force
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE hess_v_h_of_rho_r(this, nnr, rho, hess_v)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_core), TARGET, INTENT(IN) :: this
        INTEGER, INTENT(IN) :: nnr
        REAL(DP), INTENT(IN) :: rho(nnr)
        !
        REAL(DP), INTENT(OUT) :: hess_v(3, 3, nnr)
        !
        CHARACTER(LEN=80) :: sub_name = 'hess_v_h_of_rho_r'
        !
        !--------------------------------------------------------------------------------
        !
        hess_v = 0.D0
        !
        CALL io%error(sub_name, "Not implemented", 1)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE hess_v_h_of_rho_r
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE field_of_grad_rho(this, nnr, grad_rho, field)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_core), TARGET, INTENT(IN) :: this
        INTEGER, INTENT(IN) :: nnr
        REAL(DP), INTENT(IN) :: grad_rho(3, nnr)
        !
        REAL(DP), INTENT(OUT) :: field(nnr)
        !
        CHARACTER(LEN=80) :: sub_name = 'field_of_grad_rho'
        !
        !--------------------------------------------------------------------------------
        !
        field = 0.D0
        !
        CALL io%error(sub_name, "Not implemented", 1)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE field_of_grad_rho
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !                                 CORRECTION METHODS
    !
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_vperiodic(this, charges, v)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_core), TARGET, INTENT(IN) :: this
        TYPE(environ_density), TARGET, INTENT(IN) :: charges
        !
        TYPE(environ_density), INTENT(INOUT) :: v
        !
        CHARACTER(LEN=80) :: sub_name = 'calc_vperiodic'
        !
        !--------------------------------------------------------------------------------
        !
        CALL io%error(sub_name, "Not implemented", 1)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_vperiodic
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_grad_vperiodic(this, charges, grad_v)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_core), TARGET, INTENT(IN) :: this
        TYPE(environ_density), TARGET, INTENT(IN) :: charges
        !
        TYPE(environ_gradient), INTENT(INOUT) :: grad_v
        !
        CHARACTER(LEN=80) :: sub_name = 'calc_grad_vperiodic'
        !
        !--------------------------------------------------------------------------------
        !
        CALL io%error(sub_name, "Not implemented", 1)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_grad_vperiodic
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_fperiodic(this, nat, ions, auxiliary, force)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_core), TARGET, INTENT(IN) :: this
        INTEGER, INTENT(IN) :: nat
        TYPE(environ_functions), TARGET, INTENT(IN) :: ions
        TYPE(environ_density), INTENT(IN) :: auxiliary
        !
        REAL(DP), INTENT(INOUT) :: force(3, nat)
        !
        CHARACTER(LEN=80) :: sub_name = 'calc_fperiodic'
        !
        !--------------------------------------------------------------------------------
        !
        CALL io%error(sub_name, "Not implemented", 1)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_fperiodic
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_vgcs(this, electrolyte, charges, v)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_core), TARGET, INTENT(IN) :: this
        TYPE(environ_electrolyte_base), TARGET, INTENT(IN) :: electrolyte
        TYPE(environ_density), TARGET, INTENT(IN) :: charges
        !
        TYPE(environ_density), INTENT(INOUT) :: v
        !
        CHARACTER(LEN=80) :: sub_name = 'calc_vgcs'
        !
        !--------------------------------------------------------------------------------
        !
        CALL io%error(sub_name, "Not implemented", 1)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_vgcs
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_grad_vgcs(this, electrolyte, charges, grad_v)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_core), TARGET, INTENT(IN) :: this
        TYPE(environ_electrolyte_base), TARGET, INTENT(IN) :: electrolyte
        TYPE(environ_density), TARGET, INTENT(IN) :: charges
        !
        TYPE(environ_gradient), INTENT(INOUT) :: grad_v
        !
        CHARACTER(LEN=80) :: sub_name = 'calc_grad_vgcs'
        !
        !--------------------------------------------------------------------------------
        !
        CALL io%error(sub_name, "Not implemented", 1)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_grad_vgcs
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_vms(this, semiconductor, charges, v)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_core), TARGET, INTENT(IN) :: this
        TYPE(environ_semiconductor_base), TARGET, INTENT(IN) :: semiconductor
        TYPE(environ_density), TARGET, INTENT(IN) :: charges
        !
        TYPE(environ_density), INTENT(INOUT) :: v
        !
        CHARACTER(LEN=80) :: sub_name = 'calc_vms'
        !
        !--------------------------------------------------------------------------------
        !
        CALL io%error(sub_name, "Not implemented", 1)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_vms
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_grad_vms(this, semiconductor, charges, grad_v)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_core), TARGET, INTENT(IN) :: this
        TYPE(environ_semiconductor_base), TARGET, INTENT(IN) :: semiconductor
        TYPE(environ_density), TARGET, INTENT(IN) :: charges
        !
        TYPE(environ_gradient), INTENT(INOUT) :: grad_v
        !
        CHARACTER(LEN=80) :: sub_name = 'calc_grad_vms'
        !
        !--------------------------------------------------------------------------------
        !
        CALL io%error(sub_name, "Not implemented", 1)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_grad_vms
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE class_core
!----------------------------------------------------------------------------------------
