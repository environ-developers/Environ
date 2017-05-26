MODULE softcavity

  USE environ_types
  USE environ_output
  USE environ_base, ONLY : e2

  PRIVATE

  PUBLIC :: calc_vsoftcavity

CONTAINS

  SUBROUTINE calc_vsoftcavity( charges, velectrostatic, vsoftcavity, dielectric, electrolyte )

    IMPLICIT NONE

    TYPE( environ_charges ), INTENT(INOUT) :: charges
    TYPE( environ_density ), INTENT(IN) :: velectrostatic
    TYPE( environ_density ), INTENT(INOUT) :: vsoftcavity
    TYPE( environ_dielectric ), INTENT(IN), OPTIONAL :: dielectric
    TYPE( environ_electrolyte ), INTENT(IN), OPTIONAL :: electrolyte

    TYPE( environ_cell ), POINTER :: cell
    TYPE( environ_gradient ) :: gradient

    vsoftcavity % of_r = 0.D0
    IF ( .NOT. PRESENT(dielectric) .AND. .NOT. PRESENT(electrolyte) ) RETURN

    cell => velectrostatic%cell

    IF ( PRESENT( dielectric ) ) THEN
       IF ( dielectric % boundary % mode .NE. 'ionic' ) THEN
          CALL init_environ_gradient( cell, gradient )
          CALL external_gradient( velectrostatic%of_r, gradient%of_r )
          CALL update_gradient_modulus( gradient )
          vsoftcavity % of_r = - gradient % modulus % of_r * dielectric % depsilon % of_r
          vsoftcavity % of_r = vsoftcavity % of_r * 0.5D0 / fpi / e2
          CALL destroy_environ_gradient( gradient )
       END IF
    END IF

    RETURN

  END SUBROUTINE calc_vsoftcavity

END MODULE softcavity
