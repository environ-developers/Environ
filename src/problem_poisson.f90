! Copyright (C) 2018 ENVIRON (www.quantum-environment.org)
!
!    This file is part of Environ version 1.1
!
!    Environ 1.1 is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 2 of the License, or
!    (at your option) any later version.
!
!    Environ 1.1 is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more detail, either the file
!    `License' in the root directory of the present distribution, or
!    online at <http://www.gnu.org/licenses/>.
!
!----------------------------------------------------------------------------------------
!
! Authors: Oliviero Andreussi (Department of Physics, UNT)
!
!----------------------------------------------------------------------------------------
!>
!! This module contains the main drivers and routines to compute the
!! electrostatic potential that is the solution of a Poisson equation:
!!
!! \f[
!!      \nabla ^2 \phi = -4 \pi \rho
!! \f]
!!
!----------------------------------------------------------------------------------------
MODULE problem_poisson
    !------------------------------------------------------------------------------------
    !
    USE modules_constants, ONLY: DP, e2
    !
    USE core_types, ONLY: core_container
    USE cell_types, ONLY: environ_cell
    USE physical_types, ONLY: environ_charges, environ_electrolyte, environ_semiconductor
    USE representation_types, ONLY: environ_density, environ_gradient
    !
    USE utils_density, ONLY: create_environ_density, init_environ_density, &
                             destroy_environ_density
    !
    USE core_fft, ONLY: poisson_fft, gradpoisson_fft
    !
    USE correction_periodic, ONLY: calc_vperiodic, calc_gradvperiodic
    USE correction_gcs
    USE correction_ms
    USE correction_ms_gcs
    !
    USE environ_output, ONLY: verbose, ionode, environ_unit
    !
    !------------------------------------------------------------------------------------
    !
    INTERFACE poisson_direct
        MODULE PROCEDURE poisson_direct_charges, poisson_direct_density
    END INTERFACE poisson_direct
    !
    INTERFACE poisson_gradient_direct
        MODULE PROCEDURE poisson_gradient_direct_charges, poisson_gradient_direct_density
    END INTERFACE poisson_gradient_direct
    !
    !------------------------------------------------------------------------------------
    !
    PRIVATE
    !
    PUBLIC :: poisson_direct, poisson_gradient_direct
    !
    !------------------------------------------------------------------------------------
CONTAINS
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE poisson_direct_charges(core, charges, potential)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(core_container), INTENT(IN) :: core
        !
        TYPE(environ_charges), INTENT(INOUT) :: charges
        TYPE(environ_density), INTENT(INOUT) :: potential
        !
        TYPE(environ_cell), POINTER :: cell
        !
        REAL(DP) :: edummy, cdummy
        ! BACKWARD COMPATIBILITY
        ! Compatible with QE-6.0 QE-6.1.X QE-6.2.X QE-6.3
        ! REAL(DP), DIMENSION(:, :), ALLOCATABLE :: rhoaux, vaux
        ! Compatible with QE-6.4.X QE-GIT
        ! END BACKWARD COMPATIBILITY
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
        IF (core%use_fft) THEN
            ! BACKWARD COMPATIBILITY
            ! Compatible with QE-6.0 QE-6.1.X QE-6.2.X QE-6.3
            ! ALLOCATE (rhoaux(cell%nnr, core%fft%nspin))
            ! rhoaux(:, 1) = charges%density%of_r
            ! IF (core%fft%nspin == 2) rhoaux(:, 2) = 0.D0
            ! ALLOCATE (vaux(cell%nnr, core%fft%nspin))
            ! vaux = 0.D0
            ! CALL v_h_of_rho_r(rhoaux, edummy, cdummy, vaux)
            ! potential%of_r = vaux(:, 1)
            ! DEALLOCATE (rhoaux)
            ! DEALLOCATE (vaux)
            ! Compatible with QE-6.4 and QE-GIT
            ! potential%of_r = 0.D0
            ! CALL v_h_of_rho_r(charges%density%of_r, edummy, cdummy, potential%of_r)
            ! END BACKWARD COMPATIBILITY
            CALL poisson_fft(core%fft, charges%density, potential)
        ELSE IF (core%use_oned_analytic) THEN
            CALL env_errore(sub_name, 'Analytic 1D Poisson kernel is not available', 1)
        ELSE
            CALL env_errore(sub_name, 'Unexpected setup of core container', 1)
        END IF
        !
        !--------------------------------------------------------------------------------
        ! PBC corrections, if needed
        !
        IF (core%need_correction) THEN
            !
            SELECT CASE (TRIM(ADJUSTL(core%correction%type_)))
            CASE ('1da', 'oned_analytic')
                !
                CALL calc_vperiodic(core%correction%oned_analytic, charges%density, &
                                    potential)
                !
            CASE ('gcs', 'gouy-chapman', 'gouy-chapman-stern')
                !
                IF (.NOT. ASSOCIATED(charges%electrolyte)) &
                    CALL env_errore(sub_name, &
                                    'Missing electrolyte for electrochemical &
                                    &boundary correction', 1)
                !
                CALL calc_vgcs(core%correction%oned_analytic, charges%electrolyte, &
                               charges%density, potential)
                !
            CASE ('ms', 'mott-schottky')
                !
                IF (.NOT. ASSOCIATED(charges%semiconductor)) &
                    CALL env_errore(sub_name, &
                                    'Missing semiconductor for electrochemical &
                                    &boundary correction', 1)
                !
                CALL calc_vms(core%correction%oned_analytic, charges%semiconductor, &
                              charges%density, potential)
                !
            CASE ('ms-gcs', 'mott-schottky-guoy-chapman-stern')
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
                CALL calc_vms_gcs(core%correction%oned_analytic, charges%electrolyte, &
                                  charges%semiconductor, charges%density, potential)
                !
            CASE DEFAULT
                CALL env_errore(sub_name, 'Unexpected option for pbc correction core', 1)
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
    SUBROUTINE poisson_direct_density(core, charges, potential, electrolyte, &
                                      semiconductor)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(core_container), INTENT(IN) :: core
        TYPE(environ_electrolyte), INTENT(IN), OPTIONAL :: electrolyte
        !
        TYPE(environ_density), INTENT(INOUT) :: charges
        TYPE(environ_density), INTENT(INOUT) :: potential
        TYPE(environ_semiconductor), INTENT(INOUT), OPTIONAL :: semiconductor
        !
        TYPE(environ_cell), POINTER :: cell
        TYPE(environ_density) :: local
        !
        REAL(DP) :: edummy, cdummy
        ! BACKWARD COMPATIBILITY
        ! Compatible with QE-6.0 QE-6.1.X QE-6.2.X QE-6.3
        ! REAL(DP), DIMENSION(:, :), ALLOCATABLE :: rhoaux, vaux
        ! Compatible with QE-6.4.X QE-GIT
        ! END BACKWARD COMPATIBILITY
        !
        CHARACTER(LEN=80) :: sub_name = 'poisson_direct_density', llab
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
        llab = 'vloc'
        !
        CALL create_environ_density(local, llab)
        !
        CALL init_environ_density(cell, local)
        !
        IF (core%use_fft) THEN
            ! BACKWARD COMPATIBILITY
            ! Compatible with QE-6.0 QE-6.1.X QE-6.2.X QE-6.3
            ! ALLOCATE (rhoaux(cell%nnr, core%fft%nspin))
            ! rhoaux(:, 1) = charges%of_r
            ! IF (core%fft%nspin == 2) rhoaux(:, 2) = 0.D0
            ! ALLOCATE (vaux(cell%nnr, core%fft%nspin))
            ! vaux = 0.D0
            ! CALL v_h_of_rho_r(rhoaux, edummy, cdummy, vaux)
            ! local%of_r = vaux(:, 1)
            ! DEALLOCATE (rhoaux)
            ! DEALLOCATE (vaux)
            ! Compatible with QE-6.4.X QE-GIT
            ! CALL v_h_of_rho_r(charges%of_r, edummy, cdummy, local%of_r)
            ! END BACKWARD COMPATIBILITY
            CALL poisson_fft(core%fft, charges, local)
        ELSE IF (core%use_oned_analytic) THEN
            CALL env_errore(sub_name, 'Analytic 1D Poisson kernel is not available', 1)
        ELSE
            CALL env_errore(sub_name, 'Unexpected setup of core container', 1)
        END IF
        !
        !--------------------------------------------------------------------------------
        ! PBC corrections, if needed
        !
        IF (core%need_correction) THEN
            !
            SELECT CASE (TRIM(ADJUSTL(core%correction%type_)))
            CASE ('1da', 'oned_analytic')
                CALL calc_vperiodic(core%correction%oned_analytic, charges, local)
            CASE ('gcs', 'gouy-chapman', 'gouy-chapman-stern')
                !
                IF (.NOT. PRESENT(electrolyte)) &
                    CALL env_errore(sub_name, &
                                    'Missing electrolyte for electrochemical &
                                    &boundary correction', 1)
                !
                CALL calc_vgcs(core%correction%oned_analytic, electrolyte, charges, &
                               local)
                !
            CASE ('ms', 'mott-schottky')
                !
                IF (.NOT. PRESENT(semiconductor)) &
                    CALL env_errore(sub_name, &
                                    'Missing semiconductor for electrochemical &
                                    &boundary correction', 1)
                !
                CALL calc_vms(core%correction%oned_analytic, semiconductor, charges, &
                              local)
                !
            CASE ('ms-gcs', 'mott-schottky-guoy-chapman-stern')
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
                CALL calc_vms_gcs(core%correction%oned_analytic, electrolyte, &
                                  semiconductor, charges, local)
                !
            CASE DEFAULT
                !
                CALL env_errore(sub_name, 'Unexpected option for pbc correction core', 1)
                !
            END SELECT
            !
        END IF
        !
        potential%of_r = local%of_r ! only update the potential at the end
        !
        CALL destroy_environ_density(local)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE poisson_direct_density
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE poisson_gradient_direct_charges(core, charges, gradient)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(core_container), INTENT(IN) :: core
        !
        TYPE(environ_charges), INTENT(INOUT) :: charges
        TYPE(environ_gradient), INTENT(INOUT) :: gradient
        !
        CHARACTER(LEN=80) :: sub_name = 'poisson_gradient_direct_charges'
        !
        !--------------------------------------------------------------------------------
        !
        IF (core%use_fft) THEN
            !
            CALL gradpoisson_fft(core%fft, charges%density, gradient)
            !
        ELSE IF (core%use_oned_analytic) THEN
            CALL env_errore(sub_name, 'Analytic 1D Poisson kernel is not available', 1)
        ELSE
            CALL env_errore(sub_name, 'Unexpected setup of core container', 1)
        END IF
        !
        !--------------------------------------------------------------------------------
        ! PBC corrections, if needed
        !
        IF (core%need_correction) THEN
            !
            SELECT CASE (TRIM(ADJUSTL(core%correction%type_)))
            CASE ('1da', 'oned_analytic')
                !
                CALL calc_gradvperiodic(core%correction%oned_analytic, &
                                        charges%density, gradient)
                !
            CASE ('gcs', 'gouy-chapman', 'gouy-chapman-stern')
                !
                IF (.NOT. ASSOCIATED(charges%electrolyte)) &
                    CALL env_errore(sub_name, &
                                    'Missing electrolyte for electrochemical &
                                    &boundary correction', 1)
                !
                CALL calc_gradvgcs(core%correction%oned_analytic, charges%electrolyte, &
                                   charges%density, gradient)
                !
            CASE ('ms', 'mott-schottky')
                !
                IF (.NOT. ASSOCIATED(charges%semiconductor)) &
                    CALL env_errore(sub_name, &
                                    'Missing semiconductor for electrochemical &
                                    &boundary correction', 1)
                !
                CALL calc_gradvms(core%correction%oned_analytic, &
                                  charges%semiconductor, charges%density, gradient)
                !
            CASE ('ms-gcs', 'mott-schottky-guoy-chapman-stern')
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
                CALL calc_gradvms_gcs(core%correction%oned_analytic, &
                                      charges%electrolyte, charges%semiconductor, &
                                      charges%density, gradient)
                !
            CASE DEFAULT
                CALL env_errore(sub_name, 'Unexpected option for pbc correction core', 1)
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
    SUBROUTINE poisson_gradient_direct_density(core, charges, gradient, electrolyte, &
                                               semiconductor)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(core_container), INTENT(IN) :: core
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
        IF (core%use_fft) THEN
            !
            CALL gradpoisson_fft(core%fft, charges, gradient)
            !
        ELSE IF (core%use_oned_analytic) THEN
            CALL env_errore(sub_name, 'Analytic 1D Poisson kernel is not available', 1)
        ELSE
            CALL env_errore(sub_name, 'Unexpected setup of core container', 1)
        END IF
        !
        !--------------------------------------------------------------------------------
        ! PBC corrections, if needed
        !
        IF (core%need_correction) THEN
            !
            SELECT CASE (TRIM(ADJUSTL(core%correction%type_)))
            CASE ('1da', 'oned_analytic')
                !
                CALL calc_gradvperiodic(core%correction%oned_analytic, charges, &
                                        gradient)
                !
            CASE ('gcs', 'gouy-chapman', 'gouy-chapman-stern')
                !
                IF (.NOT. PRESENT(electrolyte)) &
                    CALL env_errore(sub_name, &
                                    'Missing electrolyte for &
                                    &electrochemical boundary correction', 1)
                !
                CALL calc_gradvgcs(core%correction%oned_analytic, electrolyte, &
                                   charges, gradient)
                !
            CASE ('ms', 'mott-schottky')
                !
                IF (.NOT. PRESENT(semiconductor)) &
                    CALL env_errore(sub_name, &
                                    'Missing semiconductor for &
                                    &electrochemical boundary correction', 1)
                !
                CALL calc_gradvms(core%correction%oned_analytic, semiconductor, &
                                  charges, gradient)
                !
            CASE ('ms-gcs', 'mott-schottky-guoy-chapman-stern')
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
                CALL calc_gradvms_gcs(core%correction%oned_analytic, electrolyte, &
                                      semiconductor, charges, gradient)
                !
            CASE DEFAULT
                CALL env_errore(sub_name, 'Unexpected option for pbc correction core', 1)
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
END MODULE problem_poisson
!----------------------------------------------------------------------------------------
