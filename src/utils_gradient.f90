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
MODULE utils_gradient
    !------------------------------------------------------------------------------------
    !
    USE types_cell, ONLY: environ_cell
    USE types_representation, ONLY: environ_gradient, environ_density
    !
    USE utils_density, ONLY: create_environ_density, init_environ_density, &
                             copy_environ_density, destroy_environ_density
    !
    !------------------------------------------------------------------------------------
    !
    PRIVATE
    !
    PUBLIC :: create_environ_gradient, init_environ_gradient, copy_environ_gradient, &
              update_gradient_modulus, destroy_environ_gradient
    !
    !------------------------------------------------------------------------------------
CONTAINS
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE create_environ_gradient(gradient, label)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CHARACTER(LEN=80), INTENT(IN), OPTIONAL :: label
        !
        TYPE(environ_gradient), INTENT(INOUT) :: gradient
        !
        CHARACTER(LEN=80) :: modulus_label
        !
        CHARACTER(LEN=80) :: sub_name = 'destroy_environ_density'
        !
        !--------------------------------------------------------------------------------
        !
        IF (PRESENT(label)) THEN
            gradient%label = label
            modulus_label = TRIM(ADJUSTL(label))//'_modulus'
        ELSE
            gradient%label = 'gradient'
            modulus_label = 'gradient_modulus'
        END IF
        !
        NULLIFY (gradient%cell)
        !
        IF (ALLOCATED(gradient%of_r)) &
            CALL env_errore(sub_name, 'Trying to create an already allocated object', 1)
        !
        CALL create_environ_density(gradient%modulus, modulus_label)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE create_environ_gradient
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_environ_gradient(cell, gradient)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_cell), TARGET, INTENT(IN) :: cell
        !
        TYPE(environ_gradient), INTENT(INOUT) :: gradient
        !
        CHARACTER(LEN=80) :: sub_name = 'init_environ_gradient'
        !
        !--------------------------------------------------------------------------------
        !
        gradient%update = .FALSE.
        !
        IF (ASSOCIATED(gradient%cell)) &
            CALL env_errore(sub_name, 'Trying to associate an associated object', 1)
        !
        gradient%cell => cell
        !
        IF (ALLOCATED(gradient%of_r)) &
            CALL env_errore(sub_name, 'Trying to allocate an allocated object', 1)
        !
        ALLOCATE (gradient%of_r(3, gradient%cell%nnr))
        gradient%of_r = 0.D0
        !
        CALL init_environ_density(cell, gradient%modulus)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_environ_gradient
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE copy_environ_gradient(goriginal, gcopy)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_gradient), INTENT(IN) :: goriginal
        !
        TYPE(environ_gradient), INTENT(OUT) :: gcopy
        !
        INTEGER :: n
        !
        CHARACTER(LEN=80) :: sub_name = 'copy_environ_gradient'
        !
        !--------------------------------------------------------------------------------
        !
        IF (.NOT. ASSOCIATED(goriginal%cell)) &
            CALL env_errore(sub_name, 'Trying to copy a non associated object', 1)
        !
        gcopy%cell => goriginal%cell
        !
        gcopy%update = goriginal%update
        gcopy%label = goriginal%label
        !
        IF (ALLOCATED(goriginal%of_r)) THEN
            n = SIZE(goriginal%of_r, 2)
            !
            IF (ALLOCATED(gcopy%of_r)) DEALLOCATE (gcopy%of_r)
            !
            ALLOCATE (gcopy%of_r(3, n))
            gcopy%of_r = goriginal%of_r
        END IF
        !
        CALL copy_environ_density(goriginal%modulus, gcopy%modulus)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE copy_environ_gradient
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE update_gradient_modulus(gradient)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_gradient), INTENT(INOUT) :: gradient
        !
        INTEGER, POINTER :: ir_end
        !
        !--------------------------------------------------------------------------------
        !
        ir_end => gradient%cell%ir_end
        !
        gradient%modulus%of_r(1:ir_end) = SQRT(gradient%of_r(1, 1:ir_end)**2 + &
                                               gradient%of_r(2, 1:ir_end)**2 + &
                                               gradient%of_r(3, 1:ir_end)**2)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE update_gradient_modulus
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE destroy_environ_gradient(gradient)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_gradient), INTENT(INOUT) :: gradient
        !
        CHARACTER(LEN=80) :: sub_name = 'destroy_environ_gradient'
        !
        !--------------------------------------------------------------------------------
        !
        gradient%update = .FALSE.
        !
        IF (.NOT. ASSOCIATED(gradient%cell)) &
            CALL env_errore(sub_name, 'Trying to destroy a non associated object', 1)
        !
        NULLIFY (gradient%cell)
        !
        IF (.NOT. ALLOCATED(gradient%of_r)) &
            CALL env_errore(sub_name, 'Trying to destroy a non allocated object', 1)
        !
        DEALLOCATE (gradient%of_r)
        !
        CALL destroy_environ_density(gradient%modulus)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE destroy_environ_gradient
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE utils_gradient
!----------------------------------------------------------------------------------------
