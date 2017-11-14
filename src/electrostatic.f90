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
            calc_felectrostatic, calc_vreference, calc_ereference

CONTAINS
!--------------------------------------------------------------------
  SUBROUTINE calc_vreference( setup, charges, potential )
!--------------------------------------------------------------------
    IMPLICIT NONE

    TYPE( electrostatic_setup ), INTENT(IN) :: setup
    TYPE( environ_charges ), INTENT(INOUT) :: charges
    TYPE( environ_density ), INTENT(INOUT) :: potential

    LOGICAL :: include_externals

    potential % of_r = 0.D0

    ! ... Remove external charges, if present

    include_externals = charges % include_externals
    charges % include_externals = .FALSE.

    CALL update_environ_charges( charges )

    CALL calc_velectrostatic( setup = setup, charges = charges, potential = potential )

    ! ... Reset flags

    charges % include_externals = include_externals

    CALL update_environ_charges( charges )

    RETURN
!--------------------------------------------------------------------
  END SUBROUTINE calc_vreference
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE calc_velectrostatic( setup, charges, dielectric, electrolyte, &
                                & potential )
!--------------------------------------------------------------------
    !
    ! ... Calculates the electrostatic embedding contribution to the
    !     Kohn-Sham potential
    !
    USE poisson, ONLY : poisson_direct
    USE generalized, ONLY : generalized_gradient
    USE linearized_pb, ONLY : linearized_pb_gradient
    !
    IMPLICIT NONE
    !
    ! ... Declares variables
    !
    TYPE( electrostatic_setup ), INTENT(IN) :: setup
    TYPE( environ_charges ), INTENT(INOUT) :: charges
    TYPE( environ_dielectric ), OPTIONAL, INTENT(IN) :: dielectric
    TYPE( environ_electrolyte ), OPTIONAL, INTENT(IN) :: electrolyte
    TYPE( environ_density ), INTENT(INOUT) :: potential
    !
    CHARACTER( LEN=80 ) :: sub_name = 'calc_velectrostatic'
    !
    ! ... Local variables
    !

    ! ... Calculates contribution to Hartree and local potentials

    CALL start_clock( 'calc_velect' )

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

       IF ( .NOT. PRESENT( dielectric ) ) &
            CALL errore( sub_name, 'missing details of dielectric medium', 1 )

       SELECT CASE ( setup % solver % type )

       CASE ( 'direct' )

          CALL errore( sub_name, 'option not yet implemented', 1 )
!          CALL generalized_direct()

       CASE ( 'cg', 'sd', 'iterative', 'default' )

          CALL generalized_gradient( setup % solver, setup % core, charges, dielectric, potential )

       CASE ( 'lbfgs' )

          CALL errore( sub_name, 'option not implemented', 1 )
!          CALL generalized_lbfgs()

       CASE  DEFAULT

          CALL errore( sub_name, 'unexpected solver keyword', 1 )

       END SELECT

    CASE ( 'linpb', 'linmodpb' )

       IF ( .NOT. PRESENT( electrolyte ) ) &
            CALL errore( sub_name, 'missing details of electrolyte ions', 1 )

       SELECT CASE ( setup % solver % type )

       CASE ( 'direct' )

          CALL errore( sub_name, 'option not yet implemented', 1 )
!          CALL linpb_direct()

       CASE ( 'cg', 'sd' )

!          CALL errore( sub_name, 'option not yet implemented', 1 )
          CALL linearized_pb_gradient( setup % solver, setup % core, charges, &
                                        & dielectric, electrolyte, potential )

       CASE ( 'lbfgs' )

          CALL errore( sub_name, 'option not yet implemented', 1 )
!          CALL linpb_lbfgs()

       CASE DEFAULT

          CALL errore( sub_name, 'unexpected solver keyword', 1 )

       END SELECT

    CASE ( 'pb', 'modpb' )

       IF ( .NOT. PRESENT( electrolyte ) ) &
            CALL errore( sub_name, 'missing details of electrolyte ions', 1 )

       SELECT CASE ( setup % solver % type )

       CASE ( 'direct' )

          CALL errore( sub_name, 'option not yet implemented', 1 )
!          CALL pb_direct()

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
  SUBROUTINE calc_ereference( setup, charges, potential, energy )
!--------------------------------------------------------------------
    IMPLICIT NONE

    TYPE( electrostatic_setup ), INTENT(IN) :: setup
    TYPE( environ_charges ), INTENT(INOUT) :: charges
    TYPE( environ_density ), INTENT(IN) :: potential
    REAL( DP ), INTENT(OUT) :: energy

    LOGICAL :: include_externals

    ! ... Remove external charges, if present

    include_externals = charges % include_externals
    charges % include_externals = .FALSE.

    CALL update_environ_charges( charges )

    CALL calc_eelectrostatic( setup = setup, charges = charges, potential = potential, energy = energy )

    ! ... Reset flags

    charges % include_externals = include_externals

    CALL update_environ_charges( charges )

    RETURN
!--------------------------------------------------------------------
  END SUBROUTINE calc_ereference
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE calc_eelectrostatic( setup, charges, energy, dielectric, potential, electrolyte )
!--------------------------------------------------------------------

    USE poisson, ONLY : poisson_energy
    USE generalized, ONLY: generalized_energy
    USE linearized_pb, ONLY: linearized_pb_energy

    IMPLICIT NONE

    TYPE( electrostatic_setup ), INTENT(IN) :: setup
    TYPE( environ_charges ), INTENT(INOUT) :: charges
    TYPE( environ_density ), INTENT(IN) :: potential
    REAL( DP ), INTENT(OUT) :: energy
    TYPE( environ_dielectric ), OPTIONAL, INTENT(IN) :: dielectric
    TYPE( environ_electrolyte ), OPTIONAL, INTENT(IN) :: electrolyte

    CHARACTER( LEN = 80 ) :: sub_name = 'calc_eelectrostatic'

    CALL start_clock ('calc_eelect')

    ! ... Select the right expression

    SELECT CASE ( setup % problem )

    CASE ( 'poisson' )

       CALL poisson_energy( setup % core, charges, potential, energy )

    CASE ( 'generalized' )

       IF ( .NOT. PRESENT( dielectric ) ) &
            CALL errore( sub_name, 'missing details of dielectric medium', 1 )

       CALL generalized_energy( setup % core, charges, dielectric, potential, energy )

    CASE ( 'linpb', 'linmodpb' )

       IF ( .NOT. PRESENT( electrolyte ) ) &
            CALL errore( sub_name, 'missing details of electrolyte ions', 1 )

!       CALL errore( sub_name, 'option not yet implemented', 1 )
       CALL linearized_pb_energy( setup % core, charges, dielectric, electrolyte, potential, energy )

    CASE ( 'pb', 'modpb' )

       IF ( .NOT. PRESENT( electrolyte ) ) &
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
  SUBROUTINE calc_felectrostatic( setup, natoms, charges, forces, potential, dielectric, electrolyte )
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
    TYPE( environ_charges ), INTENT(INOUT) :: charges
    REAL(DP), INTENT(INOUT) :: forces( 3, natoms )
    TYPE( environ_density ), OPTIONAL, INTENT(IN) :: potential
    TYPE( environ_dielectric ), OPTIONAL, INTENT(IN) :: dielectric
    TYPE( environ_electrolyte ), OPTIONAL, INTENT(IN) :: electrolyte
    !
    ! ... Local variables
    !
    REAL(DP)            :: ftmp( 3, natoms ), jellium
    CHARACTER( LEN=80 ) :: sub_name = 'calc_felectrostatic'
    !
    TYPE( environ_density ) :: aux
    TYPE( environ_gradient ) :: gradaux
    !
    TYPE( environ_cell ), POINTER :: cell
    !
    CALL start_clock ('calc_felect')
    !
    ! ... Sanity checks and aliases
    !
    cell => charges % density % cell
    !
    ! ... Contribution from the interaction of ionic density with
    !     modified electrostatic potential is equivalent to the interaction
    !     of ionic potential with auxiliary charge
    !
    ftmp = 0.D0
    !
    ! ... Select the right expression
    !
    SELECT CASE ( setup % problem )
       !
    CASE ( 'poisson', 'generalized' )
       !
       CALL init_environ_density( cell, aux )
       !
       IF ( charges % include_auxiliary ) THEN
          IF ( .NOT. ASSOCIATED( charges % auxiliary ) ) &
               & CALL errore(sub_name,'Missing expected charge component',1)
          aux % of_r = charges % auxiliary % density % of_r
       ELSE IF ( PRESENT( dielectric ) ) THEN
          IF ( .NOT. ALLOCATED( dielectric%gradlog%of_r ) ) &
               & CALL errore(sub_name,'Missing gradient of logarithm of dielectric',1)
          IF ( .NOT. PRESENT( potential ) ) CALL errore(sub_name,'Missing required input',1)
          CALL init_environ_gradient( cell, gradaux )
          CALL external_gradient( potential%of_r, gradaux%of_r )
          CALL scalar_product_environ_gradient( dielectric%gradlog, gradaux, aux )
          CALL destroy_environ_gradient( gradaux )
          jellium = 0.D0
          IF ( add_jellium ) jellium = charges % charge / cell % omega
          aux % of_r = aux % of_r / fpi / e2 + ( charges % density % of_r - jellium ) * &
               & ( 1.D0 - dielectric % epsilon % of_r ) / dielectric % epsilon % of_r
       ENDIF
       !
       IF ( setup % core % need_correction ) THEN
          CALL calc_fperiodic( setup%core%correction%oned_analytic, natoms, charges, aux, ftmp )
          forces = forces + ftmp
       END IF
       !
       IF ( charges % include_externals ) THEN
          IF ( .NOT. ASSOCIATED( charges % externals ) ) &
               & CALL errore(sub_name,'Missing expected charge component',1)
          aux % of_r = aux % of_r + charges % externals % density % of_r
       ENDIF
       !
       CALL external_force_lc(aux%of_r,ftmp)
       !
       CALL destroy_environ_density( aux )
       !
    CASE ( 'linpb', 'linmodpb' )

       IF ( .NOT. PRESENT( electrolyte ) ) &
            CALL errore( sub_name, 'missing details of electrolyte ions', 1 )

       CALL errore( sub_name, 'option not yet implemented', 1 )
!       CALL linpb_forces()

    CASE ( 'pb', 'modpb' )

       IF ( .NOT. PRESENT( electrolyte ) ) &
            CALL errore( sub_name, 'missing details of electrolyte ions', 1 )

       CALL errore( sub_name, 'option not yet implemented', 1 )
!       CALL pb_forces()

    CASE DEFAULT

       CALL errore( sub_name, 'unexpected problem keyword', 1 )

    END SELECT

    CALL stop_clock ('calc_felect')
    !
    !
    forces = forces + ftmp
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE calc_felectrostatic
!--------------------------------------------------------------------
!--------------------------------------------------------------------
END MODULE electrostatic
!--------------------------------------------------------------------
