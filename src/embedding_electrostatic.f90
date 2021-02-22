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
! Module containing the main drivers to compute electrostatic contributions
! to Kohn-Sham potential, total Energy and inter-atomic forces.
!
! The main inputs provided to the module are:
!      - the setup of the electrostatic calculation
!        (in an electrostatic_setup derived data type)
!      - the charges and environment details
!        (in an environ_charges derived data type)
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
!--------------------------------------------------------------------
  SUBROUTINE calc_velectrostatic( setup, charges, potential )
!--------------------------------------------------------------------
    !
    ! ... Calculates the electrostatic embedding contribution to the
    !     Kohn-Sham potential
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
    TYPE( environ_charges ), INTENT(INOUT) :: charges
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
!--------------------------------------------------------------------
  SUBROUTINE calc_eelectrostatic( setup, charges, potential, energy )
!--------------------------------------------------------------------
    !
    USE problem_poisson,       ONLY : poisson_energy
    USE problem_generalized,   ONLY : generalized_energy
    USE problem_pb,            ONLY : pb_energy
    !
    IMPLICIT NONE
    !
    TYPE( electrostatic_setup ), INTENT(IN) :: setup
    TYPE( environ_charges ), INTENT(INOUT) :: charges
    TYPE( environ_density ), INTENT(IN) :: potential
    REAL( DP ), INTENT(OUT) :: energy
    !
    CHARACTER( LEN = 80 ) :: sub_name = 'calc_eelectrostatic'
    !
    CALL start_clock ('calc_eelect')
    !
    ! ... Select the right expression
    !
    SELECT CASE ( setup % problem )
       !
    CASE ( 'poisson' )
       !
       IF ( setup % core % need_correction ) THEN
          IF ( setup % core % correction % type .EQ. 'gcs'  .OR. setup % core % correction % type .EQ. 'ms-gcs' ) THEN
             IF ( .NOT. ASSOCIATED( charges%electrolyte ) ) &
                  CALL errore( sub_name, 'missing details of electrolyte ions', 1 )
             CALL pb_energy( setup % core, charges, potential, energy )
          ELSE
             CALL poisson_energy( setup % core, charges, potential, energy )
          ENDIF
       ELSE
          CALL poisson_energy( setup % core, charges, potential, energy )
       END IF
       !
    CASE ( 'generalized' )
       !
       IF ( .NOT. ASSOCIATED( charges%dielectric ) ) &
            CALL errore( sub_name, 'missing details of dielectric medium', 1 )
       !
       IF ( setup % core % need_correction ) THEN
          IF ( setup % core % correction % type .EQ. 'gcs' .OR. setup % core % correction % type .EQ. 'ms-gcs') THEN
             IF ( .NOT. ASSOCIATED( charges%electrolyte ) ) &
                  CALL errore( sub_name, 'missing details of electrolyte ions', 1 )
             CALL pb_energy( setup % core, charges, potential, energy )
          ELSE
             CALL generalized_energy( setup % core, charges, potential, energy )
          ENDIF
       ELSE
          CALL generalized_energy( setup % core, charges, potential, energy )
       END IF
       !
    CASE ( 'pb', 'modpb', 'linpb', 'linmodpb' )
       !
       IF ( .NOT. ASSOCIATED( charges%electrolyte ) ) &
            CALL errore( sub_name, 'missing details of electrolyte ions', 1 )
       !
       CALL pb_energy( setup % core, charges, potential, energy )
       !
    CASE DEFAULT
       !
       CALL errore( sub_name, 'unexpected problem keyword', 1 )
       !
    END SELECT
    !
    CALL stop_clock ('calc_eelect')
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE calc_eelectrostatic
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE calc_felectrostatic( setup, natoms, charges, forces )
!--------------------------------------------------------------------
    !
    USE correction_periodic, ONLY : calc_fperiodic
    !
    ! ... Calculates the electrostatic embedding contribution to
    !     interatomic forces
    !
    IMPLICIT NONE
    !
    ! ... Declares variables
    !
    TYPE( electrostatic_setup ), INTENT(IN) :: setup
    INTEGER, INTENT(IN) :: natoms
    TYPE( environ_charges ), INTENT(INOUT) :: charges
    REAL(DP), INTENT(INOUT) :: forces( 3, natoms )
    !
    ! ... Local variables
    !
    TYPE( environ_cell ), POINTER :: cell
    !
    REAL(DP)            :: ftmp( 3, natoms )
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: rhoaux
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
    IF ( setup % core % use_qe_fft ) THEN
       !
       ftmp = 0.D0
       ALLOCATE( rhoaux( cell % nnr, setup % core % qe_fft % nspin ) )
       rhoaux( :, 1 ) = aux % of_r
       IF ( setup % core % qe_fft % nspin .EQ. 2 ) rhoaux( :, 2 ) = 0.D0
       CALL external_force_lc(rhoaux,ftmp)
       forces = forces + ftmp
       !
! BACKWARD COMPATIBILITY
! Compatible with QE-6.0 QE-6.1.X QE-6.2.X
!
! Compatible with QE-6.3.X and QE-GIT
       IF ( setup % core % qe_fft % use_internal_pbc_corr ) THEN
          ftmp = 0.D0
          CALL external_wg_corr_force(rhoaux,ftmp)
          forces = forces + ftmp
       END IF
! END BACKWARD COMPATIBILITY
       !
       DEALLOCATE( rhoaux )
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
