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
  USE environ_base, ONLY : e2
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

    LOGICAL :: include_externals, include_auxiliary

    potential % of_r = 0.D0

    ! ... Remove external charges, if present

    include_externals = charges % include_externals
    charges % include_externals = .FALSE.

    ! ... Remove auxiliary charges, if present and activated

    include_auxiliary = charges % include_auxiliary
    charges % include_auxiliary = .FALSE.

    CALL update_environ_charges( charges )

    CALL calc_velectrostatic( setup = setup, charges = charges, potential = potential )

    ! ... Reset flags

    charges % include_externals = include_externals
    charges % include_auxiliary = include_auxiliary

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

       CASE ( 'direct' )

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

       CASE ( 'cg', 'sd', 'iterative' )

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

          CALL errore( sub_name, 'option not yet implemented', 1 )
!          CALL linpb_gradient( )

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
  SUBROUTINE calc_ereference( charges, potential, energy )
!--------------------------------------------------------------------
    IMPLICIT NONE

    TYPE( environ_charges ), INTENT(INOUT) :: charges
    TYPE( environ_density ), INTENT(IN) :: potential
    REAL( DP ), INTENT(OUT) :: energy

    LOGICAL :: include_externals, include_auxiliary

    energy = 0.D0

    ! ... Remove external charges, if present

    include_externals = charges % include_externals
    charges % include_externals = .FALSE.

    ! ... Remove auxiliary charges, if present and activated

    include_auxiliary = charges % include_auxiliary
    charges % include_auxiliary = .FALSE.

    CALL update_environ_charges( charges )

    CALL calc_eelectrostatic( charges = charges, potential = potential, energy = energy )

    ! ... Reset flags

    charges % include_externals = include_externals
    charges % include_auxiliary = include_auxiliary

    CALL update_environ_charges( charges )

    RETURN
!--------------------------------------------------------------------
  END SUBROUTINE calc_ereference
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE calc_eelectrostatic( charges, energy, dielectric, potential, electrolyte )
!--------------------------------------------------------------------
    IMPLICIT NONE

    TYPE( environ_charges ), INTENT(INOUT) :: charges
    TYPE( environ_density ), INTENT(IN) :: potential
    REAL( DP ), INTENT(OUT) :: energy
    TYPE( environ_dielectric ), OPTIONAL, INTENT(IN) :: dielectric
    TYPE( environ_electrolyte ), OPTIONAL, INTENT(IN) :: electrolyte

    LOGICAL :: include_auxiliary

    CALL start_clock ('calc_eelect')

    energy = 0.D0

    ! ... Remove auxiliary charges, if present

    include_auxiliary = charges % include_auxiliary
    charges % include_auxiliary = .FALSE.

    CALL update_environ_charges( charges )

    energy = 0.5D0 * scalar_product_environ_density( charges%density, potential )

    ! ... Reset flags

    charges % include_auxiliary = include_auxiliary

    CALL update_environ_charges( charges )

    CALL stop_clock ('calc_eelect')

    RETURN
!--------------------------------------------------------------------
  END SUBROUTINE calc_eelectrostatic
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE calc_felectrostatic( natoms, charges, forces, potential, dielectric, electrolyte )
!--------------------------------------------------------------------
    !
    ! ... Calculates the electrostatic embedding contribution to
    !     interatomic forces
    !
    IMPLICIT NONE
    !
    ! ... Declares variables
    !
    INTEGER, INTENT(IN) :: natoms
    TYPE( environ_charges ), INTENT(INOUT) :: charges
    REAL(DP), INTENT(INOUT) :: forces( 3, natoms )
    TYPE( environ_density ), OPTIONAL, INTENT(IN) :: potential
    TYPE( environ_dielectric ), OPTIONAL, INTENT(IN) :: dielectric
    TYPE( environ_electrolyte ), OPTIONAL, INTENT(IN) :: electrolyte
    !
    ! ... Local variables
    !
    REAL(DP)            :: ftmp( 3, natoms )
    CHARACTER( LEN=80 ) :: sub_name = 'calc_felectrostatic'
    !
    TYPE( environ_density ) :: aux
    TYPE( environ_gradient ) :: gradaux
    !
    TYPE( environ_cell ), POINTER :: cell
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
       aux % of_r = aux % of_r / fpi / e2 + charges % density % of_r * &
            & ( 1.D0 - dielectric % epsilon % of_r ) / dielectric % epsilon % of_r
    ENDIF
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
