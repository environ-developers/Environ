!--------------------------------------------------------------------
MODULE problem_pb
!--------------------------------------------------------------------

  USE environ_types
  USE electrostatic_types
  USE environ_base, ONLY : e2, oldenviron

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: pb_direct, pb_energy

  INTERFACE pb_direct
     MODULE PROCEDURE pb_direct_charges, pb_direct_density
  END INTERFACE pb_direct

CONTAINS
!--------------------------------------------------------------------
  SUBROUTINE pb_direct_charges( solver, core, charges, potential )
!--------------------------------------------------------------------
    !
    USE problem_generalized, ONLY : generalized_gradient
    USE problem_poisson,     ONLY : poisson_direct
    !
    IMPLICIT NONE
    !
    TYPE( electrostatic_solver ), INTENT(IN) :: solver
    TYPE( electrostatic_core ), INTENT(IN) :: core
    TYPE( environ_charges ), INTENT(IN) :: charges
    TYPE( environ_density ), INTENT(INOUT) :: potential
    !
    CHARACTER( LEN = 80 ) :: sub_name = 'pb_direct_charges'
    !
    ! Aliases and sanity checks
    !
    IF ( .NOT. ASSOCIATED( charges % density % cell, potential % cell ) ) &
         & CALL errore(sub_name,'Missmatch in domains of charges and potential',1)
    !
    IF ( .NOT. ASSOCIATED( charges % electrolyte ) ) &
         & CALL errore(sub_name,'Missing electrolyte definition',1)
    !
    IF ( .NOT. ( core%need_correction .AND. core%correction%type .EQ. 'stern' ) ) &
         & CALL errore(sub_name,'Direct solution of PB only through analytic PBC correction',1)
    !
    ! Compute electrostatic potential
    !
    IF ( charges % electrolyte % permittivity .GT. 1.D0 ) THEN
       !
       ! If there is a dielectric embedding, follow the generalized gradient path
       !
       IF ( .NOT. ASSOCIATED( charges % dielectric ) ) &
            & CALL errore(sub_name,'Missing dielectric definition',1)
       SELECT CASE ( solver % type )
       CASE ( 'direct' )
          CALL errore(sub_name,'option not yet implemented',1)
       CASE ( 'cg', 'sd', 'iterative', 'default' )
          CALL generalized_gradient( solver, core, charges, potential, charges%electrolyte )
       CASE ( 'lbfgs' )
          CALL errore(sub_name,'option not yet implemented',1)
       CASE DEFAULT
          CALL errore(sub_name,'unexpected solver keyword',1)
       END SELECT
       !
    ELSE
       !
       ! Otherwise, in vacuum follow the poisson path
       !
       SELECT CASE ( solver % type )
       CASE ( 'direct', 'default' )
          CALL poisson_direct( core, charges, potential, charges%electrolyte )
       CASE DEFAULT
          CALL errore(sub_name,'unexpected solver keyword',1)
       END SELECT
       !
    END IF
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE pb_direct_charges
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE pb_direct_density( solver, core, charges, electrolyte, potential, dielectric )
!--------------------------------------------------------------------
    !
    USE problem_generalized, ONLY : generalized_gradient
    USE problem_poisson,     ONLY : poisson_direct
    !
    IMPLICIT NONE
    !
    TYPE( electrostatic_solver ), INTENT(IN) :: solver
    TYPE( electrostatic_core ), INTENT(IN) :: core
    TYPE( environ_density ), INTENT(IN) :: charges
    TYPE( environ_electrolyte ), INTENT(IN) :: electrolyte
    TYPE( environ_density ), INTENT(INOUT) :: potential
    TYPE( environ_dielectric ), INTENT(IN), OPTIONAL :: dielectric
    !
    TYPE( environ_cell ), POINTER :: cell
    TYPE( environ_density ) :: local
    !
    CHARACTER( LEN=80 ) :: sub_name = 'pb_direct_density'
    !
    IF ( .NOT. ASSOCIATED(charges%cell,potential%cell) ) &
         & CALL errore(sub_name,'Missmatch in domains of charges and potential',1)
    cell => charges % cell
    !
    IF ( .NOT. ( core%need_correction .AND. core%correction%type .EQ. 'stern' ) ) &
         & CALL errore(sub_name,'Direct solution of PB only through analytic PBC correction',1)
    !
    ! ... Using a local variable for the potential because the
    !     routine may be called with the same argument for charges and potential
    !
    CALL init_environ_density( cell, local )
    !
    ! Compute electrostatic potential
    !
    IF ( electrolyte % permittivity .GT. 1.D0 ) THEN
       !
       ! If there is a dielectric embedding, follow the generalized gradient path
       !
       IF ( .NOT. PRESENT( dielectric ) ) &
            & CALL errore(sub_name,'Missing dielectric definition',1)
       SELECT CASE ( solver % type )
       CASE ( 'direct' )
          CALL errore(sub_name,'option not yet implemented',1)
       CASE ( 'cg', 'sd', 'iterative', 'default' )
          CALL generalized_gradient( solver, core, charges, dielectric, local, electrolyte )
       CASE ( 'lbfgs' )
          CALL errore(sub_name,'option not yet implemented',1)
       CASE DEFAULT
          CALL errore(sub_name,'unexpected solver keyword',1)
       END SELECT
       !
    ELSE
       !
       ! Otherwise, in vacuum follow the poisson path
       !
       SELECT CASE ( solver % type )
       CASE ( 'direct', 'default' )
          CALL poisson_direct( core, charges, local, electrolyte )
       CASE DEFAULT
          CALL errore(sub_name,'unexpected solver keyword',1)
       END SELECT
       !
    END IF
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
  END SUBROUTINE pb_direct_density
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE pb_energy( core, charges, potential, energy )
!--------------------------------------------------------------------
    !
    USE problem_generalized, ONLY : generalized_energy
    USE problem_poisson,     ONLY : poisson_energy
    !
    IMPLICIT NONE
    !
    TYPE( electrostatic_core ), INTENT(IN) :: core
    TYPE( environ_charges ), INTENT(IN) :: charges
    TYPE( environ_density ), INTENT(IN) :: potential
    REAL( DP ), INTENT(OUT) :: energy
    !
    CHARACTER( LEN = 80 ) :: sub_name = 'pb_direct_charges'
    !
    ! Aliases and sanity checks
    !
    IF ( .NOT. ASSOCIATED( charges % density % cell, potential % cell ) ) &
         & CALL errore(sub_name,'Missmatch in domains of charges and potential',1)
    !
    IF ( .NOT. ASSOCIATED( charges % electrolyte ) ) &
         & CALL errore(sub_name,'Missing electrolyte definition',1)
    !
    ! Compute electrostatic energy
    !
    IF ( charges % electrolyte % permittivity .GT. 1.D0 ) THEN
       !
       ! If there is a dielectric embedding, follow the generalized gradient path
       !
       IF ( .NOT. ASSOCIATED( charges % dielectric ) ) &
            & CALL errore(sub_name,'Missing dielectric definition',1)
       CALL generalized_energy( core, charges, potential, energy )
       !
    ELSE
       !
       ! Otherwise, in vacuum follow the poisson path
       !
       CALL poisson_energy( core, charges, potential, energy )
       !
    END IF
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE pb_energy
!--------------------------------------------------------------------
!--------------------------------------------------------------------
END MODULE problem_pb
!--------------------------------------------------------------------
