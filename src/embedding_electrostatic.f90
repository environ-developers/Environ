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
!> Module containing the main drivers to compute electrostatic contributions
!! to Kohn-Sham potential, total Energy and inter-atomic forces.
!!
!! The main inputs provided to the module are:
!!      - the setup of the electrostatic calculation
!!        (in an electrostatic_setup derived data type)
!!      - the charges and environment details
!!        (in an environ_charges derived data type)
!
! Authors: Oliviero Andreussi (Department of Physics, UNT)
!          Francesco Nattino  (THEOS and NCCR-MARVEL, EPFL)
!          Nicola Marzari     (THEOS and NCCR-MARVEL, EPFL)
!
!----------------------------------------------------------------------------
MODULE embedding_electrostatic
!----------------------------------------------------------------------------
  !
  USE electrostatic_types
  USE environ_types
  USE environ_output
  USE environ_base, ONLY : e2, add_jellium
  !
  SAVE
  !
  PRIVATE
  !
  PUBLIC :: calc_velectrostatic, calc_eelectrostatic, &
       & calc_felectrostatic
  !
CONTAINS
!
!  Subroutine: calc_velectrostatic
!
!> Calculates the electrostatic embedding contribution to the
!! Kohn-Sham potential
!--------------------------------------------------------------------
  SUBROUTINE calc_velectrostatic( setup, charges, potential )
!--------------------------------------------------------------------
    !
    USE problem_poisson,       ONLY : poisson_direct
    USE problem_generalized,   ONLY : generalized_gradient
    USE problem_linearized_pb, ONLY : linearized_pb_gradient
    USE problem_pb,            ONLY : pb_nested
    !
    IMPLICIT NONE
    !
    ! ... Declares variables
    !
    TYPE( electrostatic_setup ), INTENT(IN) :: setup
    TYPE( environ_charges ), INTENT(IN) :: charges
    TYPE( environ_density ), INTENT(INOUT) :: potential
    !
    CHARACTER( LEN=80 ) :: sub_name = 'calc_velectrostatic'
    !
    ! ... Local variables
    !
    !
    ! ... Calculates contribution to Hartree and local potentials
    !
    CALL start_clock( 'calc_velect' )
    !
    potential % of_r = 0.D0
    !
    ! ... Select the appropriate combination of problem and solver
    !
    SELECT CASE ( setup % problem )
       !
    CASE ( 'poisson' )
       !
       SELECT CASE ( setup % solver % type )
          !
       CASE ( 'direct', 'default' )
          !
          CALL poisson_direct( setup % core, charges, potential )
          !
       CASE DEFAULT
          !
          CALL errore( sub_name, 'unexpected solver keyword', 1 )
          !
       END SELECT
       !
    CASE ( 'generalized' )
       !
       IF ( .NOT. ASSOCIATED( charges % dielectric ) ) &
            CALL errore( sub_name, 'missing details of dielectric medium', 1 )
       !
       SELECT CASE ( setup % solver % type )
          !
       CASE ( 'direct' )
          !
          CALL errore( sub_name, 'option not yet implemented', 1 )
!          CALL generalized_direct()
          !
       CASE ( 'cg', 'sd', 'iterative', 'default' )
          !
          CALL generalized_gradient( setup % solver, setup % core, charges, potential )
          !
       CASE ( 'lbfgs' )
          !
          CALL errore( sub_name, 'option not implemented', 1 )
!          CALL generalized_lbfgs()
          !
       CASE  DEFAULT
          !
          CALL errore( sub_name, 'unexpected solver keyword', 1 )
          !
       END SELECT
       !
    CASE ( 'linpb', 'linmodpb' )
       !
       IF ( .NOT. ASSOCIATED( charges % electrolyte ) ) &
            CALL errore( sub_name, 'missing details of electrolyte ions', 1 )
       !
       SELECT CASE ( setup % solver % type )
          !
       CASE ( 'direct' )
          !
          CALL errore( sub_name, 'option not yet implemented', 1 )
!          CALL linpb_direct()
          !
       CASE ( 'cg', 'sd' )
          !
          CALL linearized_pb_gradient( setup % solver, setup % core, charges, potential )
          !
       CASE ( 'lbfgs' )
          !
          CALL errore( sub_name, 'option not yet implemented', 1 )
!          CALL linpb_lbfgs()
          !
       CASE DEFAULT
          !
          CALL errore( sub_name, 'unexpected solver keyword', 1 )
          !
       END SELECT
       !
    CASE ( 'pb', 'modpb' )
       !
       IF ( .NOT. ASSOCIATED( charges % electrolyte ) ) &
            CALL errore( sub_name, 'missing details of electrolyte ions', 1 )
       !
       SELECT CASE ( setup % solver % type )
          !
       CASE ( 'direct' )
          !
          CALL errore( sub_name, 'option not yet implemented', 1 )
!          CALL pb_direct()
          !
       CASE ( 'iterative' )
          !
          IF ( ASSOCIATED(setup % inner) ) THEN
             !
             IF ( .NOT. ASSOCIATED( charges % dielectric ) ) &
                CALL errore( sub_name, 'missing details of dielectric medium', 1 )
             !
             CALL pb_nested( setup % solver, setup % core, charges, potential, setup % inner )
             !
          ELSE
             !
             CALL pb_nested( setup % solver, setup % core, charges, potential )
             !
          END IF
          !
       CASE ( 'newton' )
          !
          IF ( .NOT. ASSOCIATED(setup % inner) ) &
              CALL errore( sub_name, 'missing details of inner electrostatic setup', 1 )
          !
          CALL pb_nested( setup % solver, setup % core, charges, potential, setup % inner )
          !
       CASE ( 'lbfgs' )
          !
          CALL errore( sub_name, 'option not yet implemented', 1 )
!          CALL pb_lbfgs()
          !
       CASE DEFAULT
          !
          CALL errore( sub_name, 'unexpected solver keyword', 1 )
          !
       END SELECT
       !
    CASE DEFAULT
       !
       CALL errore( sub_name, 'unexpected problem keyword', 1 )
       !
    END SELECT
    !
    CALL stop_clock( 'calc_velect' )
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE calc_velectrostatic
!--------------------------------------------------------------------
!  Subroutine: calc_eelectrostatic
!
!> Calculates the electrostatic embedding contribution to the
!! energy
!--------------------------------------------------------------------
  SUBROUTINE calc_eelectrostatic( core, charges, potential, energy, add_environment )
!--------------------------------------------------------------------
!    !
!    USE problem_poisson,       ONLY : poisson_energy
!    USE problem_generalized,   ONLY : generalized_energy
!    USE problem_pb,            ONLY : pb_energy
    !
    IMPLICIT NONE
    !
    TYPE( electrostatic_core ), TARGET, INTENT(IN) :: core
    TYPE( environ_charges ), INTENT(IN) :: charges
    TYPE( environ_density ), INTENT(IN) :: potential
    LOGICAL, OPTIONAL, INTENT(IN) :: add_environment
    REAL( DP ), INTENT(OUT) :: energy
    !
    REAL( DP ) :: degauss, eself
    CHARACTER( LEN = 80 ) :: sub_name = 'calc_eelectrostatic'
    !
    CALL start_clock ('calc_eelect')
    !
    ! Aliases and sanity checks
    !
    IF ( .NOT. ASSOCIATED( charges % density % cell, potential % cell ) ) &
         & CALL errore(sub_name,'Missmatch in charges and potential domains',1)
    !
    energy = 0.D0
    eself = 0.D0
    degauss = 0.D0
    !
    ! Electrons and nuclei
    !
    IF ( core % use_fft ) THEN
       !
       energy = energy + 0.5D0 * scalar_product_environ_density( charges%density, potential )
       degauss = degauss + charges % charge
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
    !  Include environment contributions
    !
    IF ( PRESENT(add_environment) .AND. add_environment ) THEN
       !
       ! External charges
       !
       IF ( charges % include_externals ) THEN
          !
          energy = energy + 0.5D0 * scalar_product_environ_density( charges%externals%density, potential )
          degauss = degauss + charges % externals % charge
          !
       END IF
       !
       ! Polarization charge
       !
       IF ( charges % include_dielectric ) THEN
          !
          degauss = degauss + charges % dielectric % charge * 0.5D0
          !
       END IF
       !
       ! Electrolyte charge
       !
       IF ( charges % include_electrolyte ) THEN
          ! NOTE: electrolyte electrostatic interaction should be negative
          energy = energy - 0.5D0 * scalar_product_environ_density( charges%electrolyte%density, potential )
          degauss = degauss + charges % electrolyte % charge
          !
       END IF
       !
    END IF
    !
    ! Adding correction for point-like nuclei: only affects simulations of charged
    ! systems, it does not affect forces, but shift the energy depending on the
    ! fictitious Gaussian spread of the nuclei
    !
    IF ( charges % include_ions .AND. charges % ions % use_smeared_ions ) THEN
       !
       ! Compute spurious self-polarization energy
       !
       eself = charges % ions % selfenergy_correction * e2
       !
       IF ( core % use_fft ) THEN
          !
          IF ( core % fft % use_internal_pbc_corr .OR. core % need_correction ) THEN
             !
             degauss = 0.D0
             !
          ELSE
             !
             degauss = - degauss * charges % ions % quadrupole_correction * e2 * tpi &
                  & / charges % density % cell % omega
             !
          ENDIF
          !
       ENDIF
       !
    ENDIF
    !
    energy = energy + eself + degauss
    !
    CALL stop_clock ('calc_eelect')
    !
    RETURN
!!!!!
!    !
!    ! ... Select the right expression
!    !
!    SELECT CASE ( setup % problem )
!       !
!    CASE ( 'poisson' )
!       !
!       IF ( setup % core % need_correction ) THEN
!          IF ( setup % core % correction % type .EQ. 'gcs' ) THEN
!             IF ( .NOT. ASSOCIATED( charges%electrolyte ) ) &
!                  CALL errore( sub_name, 'missing details of electrolyte ions', 1 )
!             CALL pb_energy( setup % core, charges, potential, energy )
!          ELSE
!             CALL poisson_energy( setup % core, charges, potential, energy )
!          ENDIF
!       ELSE
!          CALL poisson_energy( setup % core, charges, potential, energy )
!       END IF
!       !
!    CASE ( 'generalized' )
!       !
!       IF ( .NOT. ASSOCIATED( charges%dielectric ) ) &
!            CALL errore( sub_name, 'missing details of dielectric medium', 1 )
!       !
!       IF ( setup % core % need_correction ) THEN
!          IF ( setup % core % correction % type .EQ. 'gcs' ) THEN
!             IF ( .NOT. ASSOCIATED( charges%electrolyte ) ) &
!                  CALL errore( sub_name, 'missing details of electrolyte ions', 1 )
!             CALL pb_energy( setup % core, charges, potential, energy )
!          ELSE
!             CALL generalized_energy( setup % core, charges, potential, energy )
!          ENDIF
!       ELSE
!          CALL generalized_energy( setup % core, charges, potential, energy )
!       END IF
!       !
!    CASE ( 'pb', 'modpb', 'linpb', 'linmodpb' )
!       !
!       IF ( .NOT. ASSOCIATED( charges%electrolyte ) ) &
!            CALL errore( sub_name, 'missing details of electrolyte ions', 1 )
!       !
!       CALL pb_energy( setup % core, charges, potential, energy )
!       !
!    CASE DEFAULT
!       !
!       CALL errore( sub_name, 'unexpected problem keyword', 1 )
!       !
!    END SELECT
!    !
!    CALL stop_clock ('calc_eelect')
!    !
!    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE calc_eelectrostatic
!--------------------------------------------------------------------
!  Subroutine: calc_felectrostatic
!
!> Calculates the electrostatic embedding contribution to the
!! interatomic forces
!--------------------------------------------------------------------
  SUBROUTINE calc_felectrostatic( setup, natoms, charges, forces )
!--------------------------------------------------------------------
    !
    USE core_fft,            ONLY : force_fft
    USE correction_periodic, ONLY : calc_fperiodic
    !
    IMPLICIT NONE
    !
    ! ... Declares variables
    !
    TYPE( electrostatic_setup ), INTENT(IN) :: setup
    INTEGER, INTENT(IN) :: natoms
    TYPE( environ_charges ), INTENT(IN) :: charges
    REAL(DP), INTENT(INOUT) :: forces( 3, natoms )
    !
    ! ... Local variables
    !
    TYPE( environ_cell ), POINTER :: cell
    !
    REAL(DP)            :: ftmp( 3, natoms )
! BACKWARD COMPATIBILITY
! Compatible with QE-6.0 QE-6.1.X QE-6.2.X QE-6.3.X
!    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: rhoaux
! Compatible with QE-6.4.X QE-GIT
! END BACKWARD COMPATIBILITY
    TYPE( environ_density ) :: aux
    !
    CHARACTER( LEN=80 ) :: sub_name = 'calc_felectrostatic'
    !
    CALL start_clock ('calc_felect')
    !
    ! ... Sanity checks and aliases
    !
    cell => charges % density % cell
    !
    ! ... Electrostatic contributions to forces are formulated in terms
    !     of the auxiliary charges times the derivatives of the pseudopotentials
    !
    CALL init_environ_density( cell, aux )
    !
    IF ( charges % include_dielectric ) THEN
       IF ( .NOT. ASSOCIATED( charges % dielectric ) ) &
            & CALL errore(sub_name,'Missing expected charge component',1)
       aux % of_r = aux % of_r + charges % dielectric % density % of_r
    ENDIF
    !
    IF ( charges % include_electrolyte ) THEN
       IF ( .NOT. ASSOCIATED( charges%electrolyte ) ) &
            & CALL errore(sub_name,'Missing expected charge component',1)
       aux % of_r = aux % of_r + charges % electrolyte % density % of_r
    END IF
    !
    IF ( charges % include_externals ) THEN
       IF ( .NOT. ASSOCIATED( charges % externals ) ) &
            & CALL errore(sub_name,'Missing expected charge component',1)
       aux % of_r = aux % of_r + charges % externals % density % of_r
    ENDIF
    !
    IF ( setup % core % use_fft ) THEN
       !
       ftmp = 0.D0
! BACKWARD COMPATIBILITY
! Compatible with QE-6.0 QE-6.1.X QE-6.2.X QE-6.3.X
!       ALLOCATE( rhoaux( cell % nnr, setup % core % fft % nspin ) )
!       rhoaux( :, 1 ) = aux % of_r
!       IF ( setup % core % fft % nspin .EQ. 2 ) rhoaux( :, 2 ) = 0.D0
!       CALL external_force_lc(rhoaux,ftmp)
! Compatible with QE-6.4.X and QE-GIT
!       CALL external_force_lc(aux%of_r,ftmp)
! END BACKWARD COMPATIBILITY
       CALL force_fft( setup%core%fft, aux, charges%ions, natoms, ftmp )
       forces = forces + ftmp
       !
! BACKWARD COMPATIBILITY
! Compatible with QE-6.0 QE-6.1.X QE-6.2.X
!
! Compatible with QE-6.3.X
!       IF ( setup % core % fft % use_internal_pbc_corr ) THEN
!          ftmp = 0.D0
!          CALL external_wg_corr_force(rhoaux,ftmp)
!          forces = forces + ftmp
!       END IF
!       DEALLOCATE( rhoaux )
! Compatible with QE-6.4.X and QE-GIT
!       IF ( setup % core % fft % use_internal_pbc_corr ) THEN
!          ftmp = 0.D0
!          CALL external_wg_corr_force(aux%of_r,ftmp)
!          forces = forces + ftmp
!       END IF
! END BACKWARD COMPATIBILITY
       !
    END IF
    !
    IF ( setup % core % need_correction ) THEN
       !
       ftmp = 0.D0
       CALL calc_fperiodic( setup%core%correction%oned_analytic, natoms, charges, aux, ftmp )
       forces = forces + ftmp
       !
    END IF
    !
    CALL destroy_environ_density( aux )
    !
    CALL stop_clock ('calc_felect')
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE calc_felectrostatic
!--------------------------------------------------------------------
!----------------------------------------------------------------------------
END MODULE embedding_electrostatic
!----------------------------------------------------------------------------
