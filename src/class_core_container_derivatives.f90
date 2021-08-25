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
MODULE class_core_container_derivatives
    !------------------------------------------------------------------------------------
    !
    USE class_density
    USE class_gradient
    USE class_hessian
    !
    USE class_core_container
    USE class_core_fd_derivatives
    USE class_core_fft_derivatives
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
    TYPE, EXTENDS(core_container), PUBLIC :: container_derivatives
        !--------------------------------------------------------------------------------
        !
        !--------------------------------------------------------------------------------
    CONTAINS
        !--------------------------------------------------------------------------------
        !
        PROCEDURE :: gradient => calc_gradient
        PROCEDURE :: graddot => calc_graddot
        PROCEDURE :: hessian => calc_hessian
        PROCEDURE :: laplacian => calc_laplacian
        !
        PROCEDURE, PRIVATE :: &
            calc_convolution_density, &
            calc_convolution_gradient, &
            calc_convolution_hessian
        !
        GENERIC :: convolution => &
            calc_convolution_density, &
            calc_convolution_gradient, &
            calc_convolution_hessian
        !
        !--------------------------------------------------------------------------------
    END TYPE container_derivatives
    !------------------------------------------------------------------------------------
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
    SUBROUTINE calc_gradient(this, density, gradient)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_density), INTENT(IN) :: density
        !
        CLASS(container_derivatives), INTENT(INOUT) :: this
        TYPE(environ_gradient), INTENT(INOUT) :: gradient
        !
        CHARACTER(LEN=80) :: sub_name = 'calc_gradient'
        !
        !--------------------------------------------------------------------------------
        !
        SELECT TYPE (core => this%core)
            !
        TYPE IS (core_fd_derivatives)
            CALL core%gradient(density, gradient)
            !
        TYPE IS (core_fft_derivatives)
            CALL core%gradient(density, gradient)
            !
        CLASS DEFAULT
            CALL env_errore(sub_name, 'Unexpected core', 1)
            !
        END SELECT
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_gradient
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_graddot(this, gradient, density)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_gradient), INTENT(IN) :: gradient
        !
        CLASS(container_derivatives), INTENT(INOUT) :: this
        TYPE(environ_density), INTENT(INOUT) :: density
        !
        CHARACTER(LEN=80) :: sub_name = 'calc_graddot'
        !
        !--------------------------------------------------------------------------------
        !
        SELECT TYPE (core => this%core)
            !
        TYPE IS (core_fft_derivatives)
            CALL core%graddot(gradient, density)
            !
        CLASS DEFAULT
            CALL env_errore(sub_name, 'Unexpected core', 1)
            !
        END SELECT
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_graddot
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_hessian(this, density, gradient, hessian)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_density), INTENT(IN) :: density
        !
        CLASS(container_derivatives), INTENT(INOUT) :: this
        TYPE(environ_gradient), INTENT(INOUT) :: gradient
        TYPE(environ_hessian), INTENT(INOUT) :: hessian
        !
        CHARACTER(LEN=80) :: sub_name = 'calc_hessian'
        !
        !--------------------------------------------------------------------------------
        !
        SELECT TYPE (core => this%core)
            !
        TYPE IS (core_fft_derivatives)
            CALL core%hessian(density, gradient, hessian)
            !
        CLASS DEFAULT
            CALL env_errore(sub_name, 'Unexpected core', 1)
            !
        END SELECT
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_hessian
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_laplacian(this, density, laplacian)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_density), INTENT(IN) :: density
        !
        CLASS(container_derivatives), INTENT(INOUT) :: this
        TYPE(environ_density), INTENT(INOUT) :: laplacian
        !
        CHARACTER(LEN=80) :: sub_name = 'calc_laplacian'
        !
        !--------------------------------------------------------------------------------
        !
        SELECT TYPE (core => this%core)
            !
        TYPE IS (core_fft_derivatives)
            CALL core%laplacian(density, laplacian)
            !
        CLASS DEFAULT
            CALL env_errore(sub_name, 'Unexpected core', 1)
            !
        END SELECT
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_laplacian
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_convolution_density(this, fa, fb, fc)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(container_derivatives), INTENT(IN) :: this
        TYPE(environ_density), INTENT(IN) :: fa
        TYPE(environ_density), INTENT(IN) :: fb
        !
        TYPE(environ_density), INTENT(INOUT) :: fc
        !
        CHARACTER(LEN=80) :: sub_name = 'calc_convolution_density'
        !
        !--------------------------------------------------------------------------------
        !
        SELECT TYPE (core => this%core)
            !
        TYPE IS (core_fft_derivatives)
            CALL core%convolution_density(fa, fb, fc)
            !
        CLASS DEFAULT
            CALL env_errore(sub_name, 'Unexpected core', 1)
            !
        END SELECT
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_convolution_density
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_convolution_gradient(this, fa, gb, gc)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(container_derivatives), INTENT(IN) :: this
        TYPE(environ_density), INTENT(IN) :: fa
        TYPE(environ_gradient), INTENT(IN) :: gb
        !
        TYPE(environ_gradient), INTENT(INOUT) :: gc
        !
        CHARACTER(LEN=80) :: sub_name = 'calc_convolution_gradient'
        !
        !--------------------------------------------------------------------------------
        !
        SELECT TYPE (core => this%core)
            !
        TYPE IS (core_fft_derivatives)
            CALL core%convolution_gradient(fa, gb, gc)
            !
        CLASS DEFAULT
            CALL env_errore(sub_name, 'Unexpected core', 1)
            !
        END SELECT
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_convolution_gradient
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_convolution_hessian(this, fa, hb, hc)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(container_derivatives), INTENT(IN) :: this
        TYPE(environ_density), INTENT(IN) :: fa
        TYPE(environ_hessian), INTENT(IN) :: hb
        !
        TYPE(environ_hessian), INTENT(INOUT) :: hc
        !
        CHARACTER(LEN=80) :: sub_name = 'calc_convolution_hessian'
        !
        !--------------------------------------------------------------------------------
        !
        SELECT TYPE (core => this%core)
            !
        TYPE IS (core_fft_derivatives)
            CALL core%convolution_hessian(fa, hb, hc)
            !
        CLASS DEFAULT
            CALL env_errore(sub_name, 'Unexpected core', 1)
            !
        END SELECT
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_convolution_hessian
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE class_core_container_derivatives
!----------------------------------------------------------------------------------------
