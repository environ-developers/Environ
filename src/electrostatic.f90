!
! Copyright (C) 2007-2008 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! The different modules in this file contain all the Environ related variables
! that need to be passed in and out. All module-specific variables are declared
! inside the appropriate modules.
!
! original version by O. Andreussi and N. Marzari (MIT)
!
!----------------------------------------------------------------------------
MODULE electrostatic
  !--------------------------------------------------------------------------
  !
  ! ... this module contains all the main variables and routines needed for the
  ! ... electrostatic embedding.
  !
  USE environ_types
  USE electrostatic_types
  USE environ_output
  USE environ_base, ONLY : e2, add_jellium
  !
  SAVE

  PRIVATE

  PUBLIC :: calc_velectrostatic, calc_eelectrostatic, &
            calc_felectrostatic

CONTAINS
!--------------------------------------------------------------------
  SUBROUTINE calc_velectrostatic( setup, charges, potential )
!--------------------------------------------------------------------
    !
    ! ... Calculates the electrostatic embedding contribution to the
    !     Kohn-Sham potential
    !
    USE poisson, ONLY : poisson_direct
    USE generalized, ONLY : generalized_gradient
    USE linearized_pb, ONLY : linearized_pb_gradient
    USE poisson_boltzmann, ONLY : pb_direct
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

    ! ... Calculates contribution to Hartree and local potentials

    CALL start_clock( 'calc_velect' )

    potential % of_r = 0.D0

    ! ... Select the appropriate combination of problem and solver

    SELECT CASE ( setup % problem )

    CASE ( 'poisson' )

       SELECT CASE ( setup % solver % type )

       CASE ( 'direct', 'default' )

          CALL poisson_direct( setup % core, charges, potential )

       CASE DEFAULT

          CALL errore( sub_name, 'unexpected solver keyword', 1 )

       END SELECT

    CASE ( 'generalized' )

       IF ( .NOT. ASSOCIATED( charges % dielectric ) ) &
            CALL errore( sub_name, 'missing details of dielectric medium', 1 )

       SELECT CASE ( setup % solver % type )

       CASE ( 'direct' )

          CALL errore( sub_name, 'option not yet implemented', 1 )
!          CALL generalized_direct()

       CASE ( 'cg', 'sd', 'iterative', 'default' )

          CALL generalized_gradient( setup % solver, setup % core, charges, potential )

       CASE ( 'lbfgs' )

          CALL errore( sub_name, 'option not implemented', 1 )
!          CALL generalized_lbfgs()

       CASE  DEFAULT

          CALL errore( sub_name, 'unexpected solver keyword', 1 )

       END SELECT

    CASE ( 'linpb', 'linmodpb' )

       IF ( .NOT. ASSOCIATED( charges % electrolyte ) ) &
            CALL errore( sub_name, 'missing details of electrolyte ions', 1 )

       SELECT CASE ( setup % solver % type )

       CASE ( 'direct' )

          CALL errore( sub_name, 'option not yet implemented', 1 )
!          CALL linpb_direct()

       CASE ( 'cg', 'sd' )

          CALL linearized_pb_gradient( setup % solver, setup % core, charges, potential )

       CASE ( 'lbfgs' )

          CALL errore( sub_name, 'option not yet implemented', 1 )
!          CALL linpb_lbfgs()

       CASE DEFAULT

          CALL errore( sub_name, 'unexpected solver keyword', 1 )

       END SELECT

    CASE ( 'pb', 'modpb' )

       IF ( .NOT. ASSOCIATED( charges % electrolyte ) ) &
            CALL errore( sub_name, 'missing details of electrolyte ions', 1 )

       SELECT CASE ( setup % solver % type )

       CASE ( 'direct' )

          CALL pb_direct( setup % solver, setup % core, charges, potential )

       CASE ( 'nested' )

          CALL errore( sub_name, 'option not yet implemented', 1 )
!          CALL pb_nested()

       CASE ( 'lbfgs' )

          CALL errore( sub_name, 'option not yet implemented', 1 )
!          CALL pb_lbfgs()

       CASE DEFAULT

          CALL errore( sub_name, 'unexpected solver keyword', 1 )

       END SELECT

    CASE DEFAULT

       CALL errore( sub_name, 'unexpected problem keyword', 1 )

    END SELECT

    CALL stop_clock( 'calc_velect' )

    RETURN

!--------------------------------------------------------------------
  END SUBROUTINE calc_velectrostatic
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE calc_eelectrostatic( setup, charges, potential, energy )
!--------------------------------------------------------------------

    USE poisson, ONLY : poisson_energy
    USE generalized, ONLY: generalized_energy
    USE linearized_pb, ONLY: linearized_pb_energy

    IMPLICIT NONE

    TYPE( electrostatic_setup ), INTENT(IN) :: setup
    TYPE( environ_charges ), INTENT(IN) :: charges
    TYPE( environ_density ), INTENT(IN) :: potential
    REAL( DP ), INTENT(OUT) :: energy

    CHARACTER( LEN = 80 ) :: sub_name = 'calc_eelectrostatic'

    CALL start_clock ('calc_eelect')

    ! ... Select the right expression

    SELECT CASE ( setup % problem )

    CASE ( 'poisson' )

       CALL poisson_energy( setup % core, charges, potential, energy )

    CASE ( 'generalized' )

       IF ( .NOT. ASSOCIATED( charges%dielectric ) ) &
            CALL errore( sub_name, 'missing details of dielectric medium', 1 )

       CALL generalized_energy( setup % core, charges, potential, energy )

    CASE ( 'linpb', 'linmodpb' )

       IF ( .NOT. ASSOCIATED( charges%electrolyte ) ) &
            CALL errore( sub_name, 'missing details of electrolyte ions', 1 )

!       CALL errore( sub_name, 'option not yet implemented', 1 )
       CALL linearized_pb_energy( setup % core, charges, potential, energy )

    CASE ( 'pb', 'modpb' )

       IF ( .NOT. ASSOCIATED( charges%electrolyte ) ) &
            CALL errore( sub_name, 'missing details of electrolyte ions', 1 )

       CALL errore( sub_name, 'option not yet implemented', 1 )
!       CALL pb_energy()

    CASE DEFAULT

       CALL errore( sub_name, 'unexpected problem keyword', 1 )

    END SELECT

    CALL stop_clock ('calc_eelect')

    RETURN
!--------------------------------------------------------------------
  END SUBROUTINE calc_eelectrostatic
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE calc_felectrostatic( setup, natoms, charges, forces )
!--------------------------------------------------------------------
    !
    USE periodic, ONLY : calc_fperiodic
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
    TYPE( environ_charges ), INTENT(IN) :: charges
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
       ftmp = 0.D0
       ALLOCATE( rhoaux( cell % nnr, setup % core % qe_fft % nspin ) )
       rhoaux( :, 1 ) = aux % of_r
       IF ( setup % core % qe_fft % nspin .EQ. 2 ) rhoaux( :, 2 ) = 0.D0
       CALL external_force_lc(rhoaux,ftmp)
       DEALLOCATE( rhoaux )
       forces = forces + ftmp
    END IF
    !
    IF ( setup % core % need_correction ) THEN
       ftmp = 0.D0
       CALL calc_fperiodic( setup%core%correction%oned_analytic, natoms, charges, aux, ftmp )
       forces = forces + ftmp
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
!--------------------------------------------------------------------
END MODULE electrostatic
!--------------------------------------------------------------------
