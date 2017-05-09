!--------------------------------------------------------------------
MODULE poisson
!--------------------------------------------------------------------

  USE environ_types

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: poisson_direct

  INTERFACE poisson_direct
     MODULE PROCEDURE poisson_direct_charges, poisson_direct_density
  END INTERFACE poisson_direct

CONTAINS
!--------------------------------------------------------------------
  SUBROUTINE poisson_direct_charges( charges, potential )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    TYPE( environ_charges ), INTENT(IN) :: charges
    TYPE( environ_density ), INTENT(OUT) :: potential
    !
    REAL( DP ) :: edummy, cdummy
    !
    ! TO IMPLEMENT THE CASE OF nspin .NE. 1
    !
    CALL v_h_of_rho_r( charges%density%of_r, edummy, cdummy, potential%of_r )
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE poisson_direct_charges
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE poisson_direct_density( charges, potential )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    TYPE( environ_density ), INTENT(IN) :: charges
    TYPE( environ_density ), INTENT(OUT) :: potential
    !
    REAL( DP ) :: edummy, cdummy
    !
    ! TO IMPLEMENT THE CASE OF nspin .NE. 1
    !
    CALL v_h_of_rho_r( charges%of_r, edummy, cdummy, potential%of_r )
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE poisson_direct_density
!--------------------------------------------------------------------
!--------------------------------------------------------------------
END MODULE poisson
!--------------------------------------------------------------------
