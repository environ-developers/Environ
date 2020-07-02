MODULE utils_mapping

  USE cell_types
  USE environ_types

  PRIVATE
  PUBLIC :: map_small_to_large, map_large_to_small

CONTAINS
!--------------------------------------------------------------------
  SUBROUTINE map_small_to_large( mapping, fsmall, flarge )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    TYPE( environ_mapping ), INTENT(IN) :: mapping
    TYPE( environ_density ), INTENT(IN) :: fsmall
    TYPE( environ_density ), INTENT(INOUT) :: flarge
    !
    ! Add checks on correspondence between mapping cells
    ! and input/output cells
    !
    flarge % of_r = fsmall % of_r
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE map_small_to_large
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE map_large_to_small( mapping, flarge, fsmall )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    TYPE( environ_mapping ), INTENT(IN) :: mapping
    TYPE( environ_density ), INTENT(IN) :: flarge
    TYPE( environ_density ), INTENT(INOUT) :: fsmall
    !
    ! Add checks on correspondence between mapping cells
    ! and input/output cells
    !
    fsmall % of_r = flarge % of_r
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE map_large_to_small
!--------------------------------------------------------------------
END MODULE utils_mapping
