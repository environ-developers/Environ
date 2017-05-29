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
  USE environ_output
  USE electrostatic_base
  !
  SAVE

  PRIVATE

  PUBLIC :: calc_velectrostatic, calc_eelectrostatic, &
            calc_felectrostatic, calc_vreference, calc_ereference

CONTAINS
!--------------------------------------------------------------------
  SUBROUTINE calc_vreference( charges, potential )
!--------------------------------------------------------------------
    IMPLICIT NONE

    TYPE( environ_charges ), INTENT(INOUT) :: charges
    TYPE( environ_density ), INTENT(INOUT) :: potential

    CHARACTER ( LEN = 80 ) :: local_problem
    CHARACTER ( LEN = 80 ) :: local_solver

    LOGICAL :: include_externals, include_auxiliary

    potential % of_r = 0.D0

    ! ... Compute reference potential, in the future need to pass the
    !     same set up used by the calling program

    local_problem = problem
    problem = 'poisson'
    local_solver = solver
    solver = 'direct'

    ! ... Remove external charges, if present

    include_externals = charges % include_externals
    charges % include_externals = .FALSE.

    ! ... Remove auxiliary charges, if present and activated

    include_auxiliary = charges % include_auxiliary
    charges % include_auxiliary = .FALSE.

    CALL update_environ_charges( charges )

    CALL calc_velectrostatic( charges = charges, potential = potential )

    ! ... Reset flags

    charges % include_externals = include_externals
    charges % include_auxiliary = include_auxiliary

    CALL update_environ_charges( charges )

    ! ... Reset electrostatic controls

    problem = local_problem
    solver = local_solver

    RETURN
!--------------------------------------------------------------------
  END SUBROUTINE calc_vreference
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE calc_velectrostatic( charges, dielectric, electrolyte, &
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

    potential%of_r = 0.D0

    ! ... Select the appropriate combination of problem and solver

    SELECT CASE ( problem )

    CASE ( 'poisson' )

       SELECT CASE ( solver )

       CASE ( 'direct' )

          CALL poisson_direct( charges, potential )

       CASE DEFAULT

          CALL errore( sub_name, 'unexpected solver keyword', 1 )

       END SELECT

    CASE ( 'free' )

       SELECT CASE ( solver )

       CASE ( 'direct' )

          CALL errore( sub_name, 'option not implemented', 1 )
!          CALL free_direct( charges, potential )

       CASE ( 'lbfgs' )

          CALL errore( sub_name, 'option not implemented', 1 )
!          CALL free_lbfgs()

       CASE ( 'cg', 'sd' )

          CALL errore( sub_name, 'option not implemented', 1 )
!          CALL free_gradient()

       CASE DEFAULT

          CALL errore( sub_name, 'unexpected solver keyword', 1 )

       END SELECT

    CASE ( 'generalized' )

       IF ( .NOT. PRESENT( dielectric ) ) &
            CALL errore( sub_name, 'missing details of dielectric medium', 1 )

       SELECT CASE ( solver )

       CASE ( 'direct' )

          CALL errore( sub_name, 'option not yet implemented', 1 )
!          CALL generalized_direct()

       CASE ( 'cg', 'sd', 'iterative' )

          CALL generalized_gradient( charges, dielectric, potential )

       CASE ( 'lbfgs' )

          CALL errore( sub_name, 'option not implemented', 1 )
!          CALL generalized_lbfgs()

       CASE  DEFAULT

          CALL errore( sub_name, 'unexpected solver keyword', 1 )

       END SELECT

    CASE ( 'linpb', 'linmodpb' )

       IF ( .NOT. PRESENT( electrolyte ) ) &
            CALL errore( sub_name, 'missing details of electrolyte ions', 1 )

       SELECT CASE ( solver )

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

       SELECT CASE ( solver )

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

    CHARACTER ( LEN = 80 ) :: local_problem
    CHARACTER ( LEN = 80 ) :: local_solver

    LOGICAL :: include_externals, include_auxiliary

    energy = 0.D0

    ! ... Compute reference energy, in the future need to pass the
    !     same set up used by the calling program

    local_problem = problem
    problem = 'poisson'
    local_solver = solver
    solver = 'direct'

    ! ... Remove external charges, if present

    include_externals = charges % include_externals
    charges % include_externals = .FALSE.

    ! ... Remove auxiliary charges, if present and activated

    include_auxiliary = charges % include_auxiliary
    charges % include_auxiliary = .FALSE.

    CALL calc_eelectrostatic( charges = charges, potential = potential, energy = energy )

    ! ... Reset flags

    charges % include_externals = include_externals
    charges % include_auxiliary = include_auxiliary

    ! ... Reset electrostatic controls

    problem = local_problem
    solver = local_solver

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
  SUBROUTINE calc_felectrostatic( natoms, charges, forces, dielectric, electrolyte )
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
    TYPE( environ_dielectric ), OPTIONAL, INTENT(IN) :: dielectric
    TYPE( environ_electrolyte ), OPTIONAL, INTENT(IN) :: electrolyte
    !
    ! ... Local variables
    !
    REAL(DP)            :: ftmp( 3, natoms )
    CHARACTER( LEN=80 ) :: sub_name = 'calc_felectrostatic'
    !
    ! ... Contribution from the interaction of ionic density with
    !     modified electrostatic potential
    !
    ftmp = 0.D0
    SELECT CASE ( auxiliary )
    CASE ('full')
       charges%include_ions = .FALSE.
       charges%include_electrons = .FALSE.
       CALL update_environ_charges( charges )
       CALL external_force_lc(charges%density%of_r,ftmp)
    CASE ('none')
       CALL errore(sub_name,'Option not yet implemented',1)
    CASE DEFAULT
       CALL errore(sub_name,'Option not yet implemented',1)
    END SELECT
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
