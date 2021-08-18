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
MODULE class_solver_direct
    !------------------------------------------------------------------------------------
    !
    USE environ_param, ONLY: DP
    !
    USE class_cell
    USE class_density
    USE class_gradient
    !
    USE class_core_container_electrostatics
    USE class_core_fft_electrostatics
    USE class_core_1da
    !
    USE class_solver
    !
    USE class_charges
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
    TYPE, EXTENDS(electrostatic_solver), PUBLIC :: solver_direct
        !--------------------------------------------------------------------------------
        !
        !--------------------------------------------------------------------------------
    CONTAINS
        !--------------------------------------------------------------------------------
        !
        PROCEDURE :: set_core => create_solver_direct
        PROCEDURE :: destroy => destroy_solver_direct
        !
        PROCEDURE, PRIVATE :: poisson_direct_charges, poisson_direct_density
        GENERIC :: poisson => poisson_direct_charges, poisson_direct_density
        !
        PROCEDURE, PRIVATE :: &
            poisson_gradient_direct_charges, poisson_gradient_direct_density
        !
        GENERIC :: poisson_gradient => &
            poisson_gradient_direct_charges, poisson_gradient_direct_density
        !
        !--------------------------------------------------------------------------------
    END TYPE solver_direct
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
    SUBROUTINE create_solver_direct(this, cores)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(container_electrostatics), TARGET, INTENT(IN) :: cores
        !
        CLASS(solver_direct), INTENT(INOUT) :: this
        !
        !--------------------------------------------------------------------------------
        !
        this%solver_type = 'direct'
        !
        CALL this%set_cores(cores)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE create_solver_direct
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE destroy_solver_direct(this, lflag)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        LOGICAL, INTENT(IN) :: lflag
        !
        CLASS(solver_direct), INTENT(INOUT) :: this
        !
        !--------------------------------------------------------------------------------
        !
        CALL this%cores%destroy(lflag)
        !
        NULLIFY (this%cores)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE destroy_solver_direct
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !                                   SOLVER METHODS
    !
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE poisson_direct_charges(this, charges, potential)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(solver_direct), INTENT(IN) :: this
        !
        TYPE(environ_charges), TARGET, INTENT(INOUT) :: charges
        TYPE(environ_density), INTENT(INOUT) :: potential
        !
        TYPE(environ_cell), POINTER :: cell
        !
        CHARACTER(LEN=80) :: sub_name = 'poisson_direct_charges'
        !
        !--------------------------------------------------------------------------------
        !
        IF (.NOT. ASSOCIATED(charges%density%cell, potential%cell)) &
            CALL env_errore(sub_name, 'Mismatch in domains of charges and potential', 1)
        !
        cell => charges%density%cell
        !
        CALL this%cores%poisson(charges%density, potential)
        !
        !--------------------------------------------------------------------------------
        ! PBC corrections, if needed
        !
        IF (ASSOCIATED(this%cores%correction)) THEN
            !
            SELECT CASE (TRIM(ADJUSTL(this%cores%correction%type_)))
                !
            CASE ('1da', 'oned_analytic')
                CALL this%cores%correction%calc_v(charges%density, potential)
                !
            CASE ('gcs', 'gouy-chapman', 'gouy-chapman-stern')
                !
                IF (.NOT. ASSOCIATED(charges%electrolyte)) &
                    CALL env_errore(sub_name, &
                                    'Missing electrolyte for electrochemical &
                                    &boundary correction', 1)
                !
                CALL this%cores%correction%calc_v(charges%electrolyte, &
                                                  charges%density, potential)
                !
            CASE ('ms', 'mott-schottky')
                !
                IF (.NOT. ASSOCIATED(charges%semiconductor)) &
                    CALL env_errore(sub_name, &
                                    'Missing semiconductor for electrochemical &
                                    &boundary correction', 1)
                !
                CALL this%cores%correction%calc_v(charges%semiconductor, &
                                                  charges%density, potential)
                !
            CASE ('ms-gcs', 'mott-schottky-guoy-chapman-stern') ! #TODO fix when working
                !
                IF (.NOT. ASSOCIATED(charges%semiconductor)) &
                    CALL env_errore(sub_name, &
                                    'Missing semiconductor for &electrochemical &
                                    &boundary correction', 1)
                !
                IF (.NOT. ASSOCIATED(charges%electrolyte)) &
                    CALL env_errore(sub_name, &
                                    'Missing electrolyte for electrochemical &
                                    &boundary correction', 1)
                !
                CALL env_errore(sub_name, 'ms-gcs in development', 1)
                !
                ! CALL calc_vms_gcs(this%core%correction%oned_analytic, charges%electrolyte, &
                !                   charges%semiconductor, charges%density, potential)
                !
            CASE DEFAULT
                CALL env_errore(sub_name, 'Unexpected option for pbc correction core', 1)
                !
            END SELECT
            !
        END IF
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE poisson_direct_charges
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE poisson_direct_density(this, charges, potential, electrolyte, &
                                      semiconductor)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(solver_direct), INTENT(IN) :: this
        TYPE(environ_electrolyte), INTENT(IN), OPTIONAL :: electrolyte
        !
        TYPE(environ_density), INTENT(INOUT) :: charges
        TYPE(environ_density), INTENT(INOUT) :: potential
        TYPE(environ_semiconductor), INTENT(INOUT), OPTIONAL :: semiconductor
        !
        TYPE(environ_cell), POINTER :: cell
        TYPE(environ_density) :: local
        !
        CHARACTER(LEN=80) :: sub_name = 'poisson_direct_density'
        !
        !--------------------------------------------------------------------------------
        !
        IF (.NOT. ASSOCIATED(charges%cell, potential%cell)) &
            CALL env_errore(sub_name, 'Mismatch in domains of charges and potential', 1)
        !
        cell => charges%cell
        !
        !--------------------------------------------------------------------------------
        ! Using a local variable for the potential because the routine may be
        ! called with the same argument for charges and potential
        !
        CALL local%create()
        !
        CALL local%init(cell)
        !
        CALL this%cores%poisson(charges, local)
        !
        !--------------------------------------------------------------------------------
        ! PBC corrections, if needed
        !
        IF (ASSOCIATED(this%cores%correction)) THEN
            !
            SELECT CASE (TRIM(ADJUSTL(this%cores%correction%type_)))
                !
            CASE ('1da', 'oned_analytic')
                CALL this%cores%correction%calc_v(charges, local)
                !
            CASE ('gcs', 'gouy-chapman', 'gouy-chapman-stern')
                !
                IF (.NOT. PRESENT(electrolyte)) &
                    CALL env_errore(sub_name, &
                                    'Missing electrolyte for electrochemical &
                                    &boundary correction', 1)
                !
                CALL this%cores%correction%calc_v(electrolyte, charges, local)
                !
            CASE ('ms', 'mott-schottky')
                !
                IF (.NOT. PRESENT(semiconductor)) &
                    CALL env_errore(sub_name, &
                                    'Missing semiconductor for electrochemical &
                                    &boundary correction', 1)
                !
                CALL this%cores%correction%calc_v(semiconductor, charges, local)
                !
            CASE ('ms-gcs', 'mott-schottky-guoy-chapman-stern') ! #TODO fix when working
                !
                IF (.NOT. PRESENT(semiconductor)) &
                    CALL env_errore(sub_name, &
                                    'Missing semiconductor for electrochemical &
                                    &boundary correction', 1)
                !
                IF (.NOT. PRESENT(electrolyte)) &
                    CALL env_errore(sub_name, &
                                    'Missing electrolyte for electrochemical &
                                    &boundary correction', 1)
                !
                CALL env_errore(sub_name, 'ms-gcs in development', 1)
                !
                ! CALL calc_vms_gcs(this%core%correction%oned_analytic, electrolyte, &
                !                   semiconductor, charges, local)
                !
            CASE DEFAULT
                CALL env_errore(sub_name, 'Unexpected option for pbc correction core', 1)
                !
            END SELECT
            !
        END IF
        !
        potential%of_r = local%of_r ! only update the potential at the end
        !
        CALL local%destroy()
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE poisson_direct_density
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE poisson_gradient_direct_charges(this, charges, gradient)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(solver_direct), INTENT(IN) :: this
        !
        TYPE(environ_charges), TARGET, INTENT(INOUT) :: charges
        TYPE(environ_gradient), INTENT(INOUT) :: gradient
        !
        CHARACTER(LEN=80) :: sub_name = 'poisson_gradient_direct_charges'
        !
        !--------------------------------------------------------------------------------
        !
        CALL this%cores%gradpoisson(charges%density, gradient)
        !
        !--------------------------------------------------------------------------------
        ! PBC corrections, if needed
        !
        IF (ASSOCIATED(this%cores%correction)) THEN
            !
            SELECT CASE (TRIM(ADJUSTL(this%cores%correction%type_)))
                !
            CASE ('1da', 'oned_analytic')
                CALL this%cores%correction%calc_gradv(charges%density, gradient)
                !
            CASE ('gcs', 'gouy-chapman', 'gouy-chapman-stern')
                !
                IF (.NOT. ASSOCIATED(charges%electrolyte)) &
                    CALL env_errore(sub_name, &
                                    'Missing electrolyte for electrochemical &
                                    &boundary correction', 1)
                !
                CALL this%cores%correction%calc_gradv(charges%electrolyte, &
                                                      charges%density, gradient)
                !
            CASE ('ms', 'mott-schottky')
                !
                IF (.NOT. ASSOCIATED(charges%semiconductor)) &
                    CALL env_errore(sub_name, &
                                    'Missing semiconductor for electrochemical &
                                    &boundary correction', 1)
                !
                CALL this%cores%correction%calc_gradv(charges%semiconductor, &
                                                      charges%density, gradient)
                !
            CASE ('ms-gcs', 'mott-schottky-guoy-chapman-stern') ! #TODO fix when working
                !
                IF (.NOT. ASSOCIATED(charges%semiconductor)) &
                    CALL env_errore(sub_name, &
                                    'Missing semiconductor for electrochemical &
                                    &boundary correction', 1)
                !
                IF (.NOT. ASSOCIATED(charges%electrolyte)) &
                    CALL env_errore(sub_name, &
                                    'Missing electrolyte for electrochemical &
                                    &boundary correction', 1)
                !
                CALL env_errore(sub_name, 'ms-gcs in development', 1)
                !
                ! CALL calc_gradvms_gcs(this%core%correction%oned_analytic, &
                !                       charges%electrolyte, charges%semiconductor, &
                !                       charges%density, gradient)
                !
            CASE DEFAULT
                CALL env_errore(sub_name, 'Unexpected option for pbc correction core', 1)
                !
            END SELECT
            !
        END IF
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE poisson_gradient_direct_charges
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE poisson_gradient_direct_density(this, charges, gradient, electrolyte, &
                                               semiconductor)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(solver_direct), INTENT(IN) :: this
        TYPE(environ_electrolyte), INTENT(IN), OPTIONAL :: electrolyte
        !
        TYPE(environ_density), INTENT(INOUT) :: charges
        TYPE(environ_gradient), INTENT(INOUT) :: gradient
        TYPE(environ_semiconductor), INTENT(INOUT), OPTIONAL :: semiconductor
        !
        CHARACTER(LEN=80) :: sub_name = 'poisson_gradient_direct_density'
        !
        !--------------------------------------------------------------------------------
        !
        CALL this%cores%gradpoisson(charges, gradient)
        !
        !--------------------------------------------------------------------------------
        ! PBC corrections, if needed
        !
        IF (ASSOCIATED(this%cores%correction)) THEN
            !
            SELECT CASE (TRIM(ADJUSTL(this%cores%correction%type_)))
                !
            CASE ('1da', 'oned_analytic')
                CALL this%cores%correction%calc_gradv(charges, gradient)
                !
            CASE ('gcs', 'gouy-chapman', 'gouy-chapman-stern')
                !
                IF (.NOT. PRESENT(electrolyte)) &
                    CALL env_errore(sub_name, &
                                    'Missing electrolyte for &
                                    &electrochemical boundary correction', 1)
                !
                CALL this%cores%correction%calc_gradv(electrolyte, charges, gradient)
                !
            CASE ('ms', 'mott-schottky')
                !
                IF (.NOT. PRESENT(semiconductor)) &
                    CALL env_errore(sub_name, &
                                    'Missing semiconductor for &
                                    &electrochemical boundary correction', 1)
                !
                CALL this%cores%correction%calc_gradv(semiconductor, charges, gradient)
                !
            CASE ('ms-gcs', 'mott-schottky-guoy-chapman-stern') ! #TODO fix when working
                !
                IF (.NOT. PRESENT(semiconductor)) &
                    CALL env_errore(sub_name, &
                                    'Missing semiconductor for &
                                    &electrochemical boundary correction', 1)
                !
                IF (.NOT. PRESENT(electrolyte)) &
                    CALL env_errore(sub_name, &
                                    'Missing electrolyte for &
                                    &electrochemical boundary correction', 1)
                !
                CALL env_errore(sub_name, 'ms-gcs in development', 1)
                !
                ! CALL calc_gradvms_gcs(this%core%correction%oned_analytic, electrolyte, &
                !                       semiconductor, charges, gradient)
                !
            CASE DEFAULT
                CALL env_errore(sub_name, 'Unexpected option for pbc correction core', 1)
                !
            END SELECT
            !
        END IF
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE poisson_gradient_direct_density
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE class_solver_direct
!----------------------------------------------------------------------------------------