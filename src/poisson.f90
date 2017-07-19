!--------------------------------------------------------------------
MODULE poisson
!--------------------------------------------------------------------

  USE environ_types

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: poisson_direct, poisson_gradient_direct

  INTERFACE poisson_direct
     MODULE PROCEDURE poisson_direct_charges, poisson_direct_density
  END INTERFACE poisson_direct

  INTERFACE poisson_gradient_direct
     MODULE PROCEDURE poisson_gradient_direct_charges, poisson_gradient_direct_density
  END INTERFACE poisson_gradient_direct

CONTAINS
!--------------------------------------------------------------------
  SUBROUTINE poisson_direct_charges( charges, potential )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    TYPE( environ_charges ), INTENT(IN) :: charges
    TYPE( environ_density ), INTENT(INOUT) :: potential
    !
    REAL( DP ) :: edummy, cdummy
    !
    ! TO IMPLEMENT THE CASE OF nspin .NE. 1
    !
    potential % of_r = 0.D0
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
    TYPE( environ_density ), INTENT(INOUT) :: potential
    !
    TYPE( environ_cell ), POINTER :: cell
    TYPE( environ_density ) :: local
    !
    REAL( DP ) :: edummy, cdummy
    CHARACTER( LEN=80 ) :: sub_name = 'poisson_direct_density'
    !
    ! TO IMPLEMENT THE CASE OF nspin .NE. 1
    !
    IF ( .NOT. ASSOCIATED(charges%cell,potential%cell) ) &
         & CALL errore(sub_name,'Missmatch in domains of charges and potential',1)
    cell => charges % cell
    !
    CALL init_environ_density( cell, local )
    !
    CALL v_h_of_rho_r( charges%of_r, edummy, cdummy, local%of_r )
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
  SUBROUTINE poisson_gradient_direct_charges( charges, gradient )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    TYPE( environ_charges ), INTENT(IN) :: charges
    TYPE( environ_gradient ), INTENT(INOUT) :: gradient
    !
    ! TO IMPLEMENT THE CASE OF nspin .NE. 1
    !
    CALL gradv_h_of_rho_r( charges%density%of_r, gradient%of_r )
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE poisson_gradient_direct_charges
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE poisson_gradient_direct_density( charges, gradient )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    TYPE( environ_density ), INTENT(IN) :: charges
    TYPE( environ_gradient ), INTENT(INOUT) :: gradient
    !
    ! TO IMPLEMENT THE CASE OF nspin .NE. 1
    !
    CALL gradv_h_of_rho_r( charges%of_r, gradient%of_r )
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE poisson_gradient_direct_density
!--------------------------------------------------------------------
!--------------------------------------------------------------------
END MODULE poisson
!--------------------------------------------------------------------
