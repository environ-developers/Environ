! Copyright (C) 2018 ENVIRON (www.quantum-environment.org)
!
!    This file is part of Environ version 1.0
!
!    Environ 1.0 is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 2 of the License, or
!    (at your option) any later version.
!
!    Environ 1.0 is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more detail, either the file
!    `License' in the root directory of the present distribution, or
!    online at <http://www.gnu.org/licenses/>.
!
! This module contains the main drivers and routines to compute the
! electrostatic potential that is the solution of a Poisson equation:
!
!      \nabla ^2 \phi = -4 \pi \rho
!
! At this time the subroutines are mostly wrappers for the FFT
! solver of Quantum ESPRESSO.
!
! Authors: Oliviero Andreussi (Department of Physics, UNT)
!
!--------------------------------------------------------------------
MODULE problem_poisson
!--------------------------------------------------------------------
  !
  USE environ_types
  USE electrostatic_types
  USE correction_periodic
  USE correction_gcs
  USE environ_base, ONLY : e2, oldenviron
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC :: poisson_direct, poisson_gradient_direct, poisson_energy
  !
  INTERFACE poisson_direct
     MODULE PROCEDURE poisson_direct_charges, poisson_direct_density
  END INTERFACE poisson_direct
  !
  INTERFACE poisson_gradient_direct
     MODULE PROCEDURE poisson_gradient_direct_charges, poisson_gradient_direct_density
  END INTERFACE poisson_gradient_direct
  !
  INTERFACE poisson_energy
     MODULE PROCEDURE poisson_energy_charges, poisson_energy_density
  END INTERFACE poisson_energy
  !
CONTAINS
!--------------------------------------------------------------------
  SUBROUTINE poisson_direct_charges( core, charges, potential )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    TYPE( electrostatic_core ), INTENT(IN) :: core
    TYPE( environ_charges ), INTENT(IN) :: charges
    TYPE( environ_density ), INTENT(INOUT) :: potential
    !
    TYPE( environ_cell ), POINTER :: cell
    !
    REAL( DP ) :: edummy, cdummy
    REAL( DP ), DIMENSION(:,:), ALLOCATABLE :: rhoaux, vaux
    CHARACTER( LEN = 80 ) :: sub_name = 'poisson_direct_charges'
    !
    ! Aliases and sanity checks
    !
    IF ( .NOT. ASSOCIATED( charges % density % cell, potential % cell ) ) &
         & CALL errore(sub_name,'Missmatch in domains of charges and potential',1)
    cell => charges % density % cell
    !
    IF ( core % use_qe_fft ) THEN
       !
       ALLOCATE( rhoaux ( cell % nnr, core % qe_fft % nspin ) )
       rhoaux( :, 1 ) = charges % density % of_r
       IF ( core % qe_fft % nspin .EQ. 2 ) rhoaux( :, 2 ) = 0.D0
       ALLOCATE( vaux ( cell % nnr, core % qe_fft % nspin ) )
       vaux = 0.D0
       CALL v_h_of_rho_r( rhoaux, edummy, cdummy, vaux )
       potential % of_r = vaux( :, 1 )
       DEALLOCATE( rhoaux )
       DEALLOCATE( vaux )
       !
    ELSE IF ( core % use_oned_analytic ) THEN
       !
       CALL errore(sub_name,'Analytic 1D Poisson kernel is not available',1)
       !
    ELSE
       !
       CALL errore(sub_name,'Unexpected setup of electrostatic core',1)
       !
    ENDIF
    !
    IF ( core % need_correction ) THEN
       !
       SELECT CASE ( TRIM( ADJUSTL( core%correction%type ) ) )
          !
       CASE ( '1da', 'oned_analytic' )
          !
          CALL calc_vperiodic( core%correction%oned_analytic, charges%density, potential )
          !
       CASE ( 'gcs', 'gouy-chapman', 'gouy-chapman-stern' )
          !
          IF ( .NOT. ASSOCIATED( charges%electrolyte ) ) &
               & CALL errore(sub_name,'Missing electrolyte for electrochemical boundary correction',1)
          CALL calc_vgcs( core%correction%oned_analytic, charges%electrolyte, charges%density, potential )
          !
       CASE ( 'ms' , 'mott-schottky')
          IF ( .NOT. ASSOCIATED( charges%semiconductor ) ) &
               & CALL errore(sub_name,'Missing semiconductor for electrochemical boundary correction',1)
          CALL calc_vms( core%correction%oned_analytic, charges%semiconductor, charges%density, potential )
       CASE DEFAULT
          !
          CALL errore(sub_name,'Unexpected option for pbc correction core',1)
          !
       END SELECT
       !
    ENDIF
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE poisson_direct_charges
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE poisson_direct_density( core, charges, potential, electrolyte, semiconductor )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    TYPE( electrostatic_core ), INTENT(IN) :: core
    TYPE( environ_density ), INTENT(IN) :: charges
    TYPE( environ_density ), INTENT(INOUT) :: potential
    TYPE( environ_electrolyte), INTENT(IN), OPTIONAL :: electrolyte
    TYPE( environ_semiconductor), INTENT(IN), OPTIONAL :: semiconductor
    !
    TYPE( environ_cell ), POINTER :: cell
    TYPE( environ_density ) :: local
    !
    REAL( DP ) :: edummy, cdummy
    REAL( DP ), DIMENSION( :, : ), ALLOCATABLE :: rhoaux, vaux
    CHARACTER( LEN=80 ) :: sub_name = 'poisson_direct_density'
    !
    IF ( .NOT. ASSOCIATED(charges%cell,potential%cell) ) &
         & CALL errore(sub_name,'Missmatch in domains of charges and potential',1)
    cell => charges % cell
    !
    ! ... Using a local variable for the potential because the
    !     routine may be called with the same argument for charges and potential
    !
    CALL init_environ_density( cell, local )
    !
    IF ( core % use_qe_fft ) THEN
       !
       ALLOCATE( rhoaux ( cell % nnr, core % qe_fft % nspin ) )
       rhoaux( :, 1 ) = charges % of_r
       IF ( core % qe_fft % nspin .EQ. 2 ) rhoaux( :, 2 ) = 0.D0
       ALLOCATE( vaux ( cell % nnr, core % qe_fft % nspin ) )
       vaux = 0.D0
       CALL v_h_of_rho_r( rhoaux, edummy, cdummy, vaux )
       local % of_r = vaux( :, 1 )
       DEALLOCATE( rhoaux )
       DEALLOCATE( vaux )
       !
    ELSE IF ( core % use_oned_analytic ) THEN
       !
       CALL errore(sub_name,'Analytic 1D Poisson kernel is not available',1)
       !
    ELSE
       !
       CALL errore(sub_name,'Unexpected setup of electrostatic core',1)
       !
    ENDIF
    !
    ! ... PBC corrections, if needed
    !
    IF ( core % need_correction ) THEN
       !
       SELECT CASE ( TRIM( ADJUSTL( core%correction%type ) ) )
          !
       CASE ( '1da', 'oned_analytic' )
          !
          CALL calc_vperiodic( core%correction%oned_analytic, charges, local )
          !
       CASE ( 'gcs','gouy-chapman', 'gouy-chapman-stern' )
          !
          IF ( .NOT. PRESENT( electrolyte ) ) &
               & CALL errore(sub_name,'Missing electrolyte for electrochemical boundary correction',1)
          CALL calc_vgcs( core%correction%oned_analytic, electrolyte, charges, local )
          !
       CASE ( 'ms','mott-schottky' )
          !
          IF ( .NOT. PRESENT( semiconductor ) ) &
               & CALL errore(sub_name,'Missing semiconductor for electrochemical boundary correction',1)
          CALL calc_vgcs( core%correction%oned_analytic, semiconductor, charges, local )
       CASE DEFAULT
          !
          CALL errore(sub_name,'Unexpected option for pbc correction core',1)
          !
       END SELECT
       !
    ENDIF
    !
    ! ... Only update the potential at the end
    !
    potential % of_r = local % of_r
    !
    CALL destroy_environ_density( local )
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE poisson_direct_density
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE poisson_gradient_direct_charges( core, charges, gradient )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    TYPE( electrostatic_core ), INTENT(IN) :: core
    TYPE( environ_charges ), INTENT(IN) :: charges
    TYPE( environ_gradient ), INTENT(INOUT) :: gradient
    !
    CHARACTER( LEN = 80 ) :: sub_name = 'poisson_gradient_direct_charges'
    !
    IF ( core % use_qe_fft ) THEN
       !
       CALL gradv_h_of_rho_r( charges%density%of_r, gradient%of_r )
       !
    ELSE IF ( core % use_oned_analytic ) THEN
       !
       CALL errore(sub_name,'Analytic 1D Poisson kernel is not available',1)
       !
    ELSE
       !
       CALL errore(sub_name,'Unexpected setup of electrostatic core',1)
       !
    ENDIF
    !
    IF ( core % need_correction ) THEN
       !
       SELECT CASE ( TRIM( ADJUSTL( core%correction%type ) ) )
          !
       CASE ( '1da', 'oned_analytic' )
          !
          CALL calc_gradvperiodic( core%correction%oned_analytic, charges%density, gradient )
          !
       CASE ( 'gcs','gouy-chapman', 'gouy-chapman-stern' )
          !
          IF ( .NOT. ASSOCIATED( charges%electrolyte ) ) &
               & CALL errore(sub_name,'Missing electrolyte for electrochemical boundary correction',1)
          CALL calc_gradvgcs( core%correction%oned_analytic, charges%electrolyte, charges%density, gradient )
          !
       CASE ( 'ms','mott-schottky')
          IF ( .NOT. ASSOCIATED( charges%semiconductor ) ) &
              & CALL errore(sub_name,'Missing semiconductor for electrochemical boundary correction',1)
          CALL calc_gradvms( core%correction%oned_analytic, charges%semiconductor, charges%density, gradient )          

       CASE DEFAULT
          !
          CALL errore(sub_name,'Unexpected option for pbc correction core',1)
          !
       END SELECT
       !
    ENDIF
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE poisson_gradient_direct_charges
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE poisson_gradient_direct_density( core, charges, gradient, electrolyte,&
            & semiconductor )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    TYPE( electrostatic_core ), INTENT(IN) :: core
    TYPE( environ_density ), INTENT(IN) :: charges
    TYPE( environ_gradient ), INTENT(INOUT) :: gradient
    TYPE( environ_electrolyte ), INTENT(IN), OPTIONAL :: electrolyte
    TYPE( environ_semiconductor ), INTENT(IN), OPTIONAL :: semiconductor
    !
    CHARACTER( LEN = 80 ) :: sub_name = 'poisson_gradient_direct_density'
    !
    IF ( core % use_qe_fft ) THEN
       !
       CALL gradv_h_of_rho_r( charges%of_r, gradient%of_r )
       !
    ELSE IF ( core % use_oned_analytic ) THEN
       !
       CALL errore(sub_name,'Analytic 1D Poisson kernel is not available',1)
       !
    ELSE
       !
       CALL errore(sub_name,'Unexpected setup of electrostatic core',1)
       !
    ENDIF
    !
    IF ( core % need_correction ) THEN
       !
       SELECT CASE ( TRIM( ADJUSTL( core%correction%type ) ) )
          !
       CASE ( '1da', 'oned_analytic' )
          !
          CALL calc_gradvperiodic( core%correction%oned_analytic, charges, gradient )
          !
       CASE ( 'gcs','gouy-chapman', 'gouy-chapman-stern' )
          !
          IF ( .NOT. PRESENT( electrolyte ) ) &
               & CALL errore(sub_name,'Missing electrolyte for electrochemical boundary correction',1)
          CALL calc_gradvgcs( core%correction%oned_analytic, electrolyte, charges, gradient )
          !
       CASE ( 'ms','mott-schottky')
          !
          IF ( .NOT. PRESENT( semiconductor ) ) &
               & CALL errore(sub_name,'Missing semiconductor for electrochemical boundary correction', 1)
          CALL calc_gradvgcs( core%correction%oned_analytic, semiconductor, charges, gradient )

       CASE DEFAULT
          !
          CALL errore(sub_name,'Unexpected option for pbc correction core',1)
          !
       END SELECT
       !
    ENDIF
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE poisson_gradient_direct_density
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE poisson_energy_charges( core, charges, potential, energy )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    TYPE( electrostatic_core ), INTENT(IN) :: core
    TYPE( environ_charges ), INTENT(IN) :: charges
    TYPE( environ_density ), INTENT(IN) :: potential
    REAL( DP ), INTENT(OUT) :: energy
    !
    REAL( DP ) :: degauss, eself
    CHARACTER( LEN = 80 ) :: sub_name = 'poisson_energy_charges'
    !
    energy = 0.d0
    !
    IF ( core % use_qe_fft ) THEN
       !
       energy = 0.5D0 * scalar_product_environ_density( charges%density, potential )
       !
    ELSE IF ( core % use_oned_analytic ) THEN
       !
       CALL errore(sub_name,'Analytic 1D Poisson kernel is not available',1)
       !
    ELSE
       !
       CALL errore(sub_name,'Unexpected setup of electrostatic core',1)
       !
    ENDIF
    !
    ! Remove self-interaction and correct energy of Gaussian ions
    !
    eself = 0.D0
    degauss = 0.D0
    !
    IF ( charges % include_ions .AND. charges % ions % use_smeared_ions ) THEN
       !
       eself = charges % ions % selfenergy_correction * e2
       !
       IF ( core % use_qe_fft ) THEN
          !
          IF ( core % qe_fft % use_internal_pbc_corr .OR. core % need_correction ) THEN
             !
             degauss = 0.D0
             !
          ELSE
             !
             degauss = - charges % ions % quadrupole_correction * charges % charge * e2 * tpi &
                  & / charges % density % cell % omega
             !
          END IF
          !
       ENDIF
       !
       degauss = degauss ! we are missing the difference in interaction of electrons with gaussian vs. pc nuclei
       !
    ENDIF
    !
    energy = energy + eself + degauss
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE poisson_energy_charges
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE poisson_energy_density( core, charges, potential, energy )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    TYPE( electrostatic_core ), INTENT(IN) :: core
    TYPE( environ_density ), INTENT(IN) :: charges
    TYPE( environ_density ), INTENT(IN) :: potential
    REAL( DP ), INTENT(OUT) :: energy
    !
    REAL( DP ) :: degauss, eself
    CHARACTER( LEN = 80 ) :: sub_name = 'poisson_energy_density'
    !
    energy = 0.D0
    !
    IF ( core % use_qe_fft ) THEN
       !
       energy = 0.5D0 * scalar_product_environ_density( charges, potential )
       !
    ELSE IF ( core % use_oned_analytic ) THEN
       !
       CALL errore(sub_name,'Analytic 1D Poisson kernel is not available',1)
       !
    ELSE
       !
       CALL errore(sub_name,'Unexpected setup of electrostatic core',1)
       !
    ENDIF
    !
    ! Remove self-interaction and correct energy of Gaussian ions
    ! WARNING: THESE CORRECTIONS ARE ONLY POSSIBLE WHEN INFORMATIONS
    ! ABOUT THE POINT-LIKE NUCLEI ARE PASSED
    !
    eself = 0.D0
    !
    degauss = 0.D0
    !
    energy = energy + eself + degauss
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE poisson_energy_density
!--------------------------------------------------------------------
!--------------------------------------------------------------------
END MODULE problem_poisson
!--------------------------------------------------------------------
