!--------------------------------------------------------------------
MODULE ioncc
!--------------------------------------------------------------------

  USE kinds,          ONLY: DP
  USE constants,      ONLY: k_boltzmann_ry, pi, tpi, fpi
  USE io_global,      ONLY: stdout
  USE mp,             ONLY: mp_sum
  USE mp_bands,       ONLY: intra_bgrp_comm
  USE environ_cell,   ONLY: domega, alat, omega, ntot
  USE environ_ions,   ONLY: avg_pos, rhoions
  USE environ_base,   ONLY: verbose, e2, ir_end, environ_unit,           & 
                            env_ioncc_level, nrep, cion, zion,           &
                            solvent_temperature, rhoioncc, rhopolcc,     & 
                            env_static_permittivity, rhopol, slab_axis,  &
                            rhomin, rhopb
  USE generate_function, ONLY: planar_average
  USE periodic,       ONLY: calc_vperiodic

  IMPLICIT NONE

  SAVE

  PRIVATE

  PUBLIC :: ioncc_initbase, ioncc_initcell, ioncc_clean, calc_vioncc, &
            calc_eioncc

CONTAINS

!--------------------------------------------------------------------
  SUBROUTINE ioncc_initbase( nnr )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: nnr
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE ioncc_initbase
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE ioncc_initcell( nnr, n1, n2, n3, at )
!--------------------------------------------------------------------
    !
    USE generate_function, ONLY: generate_axis
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: nnr, n1, n2, n3
    REAL(DP), INTENT(IN) :: at(3,3)
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE ioncc_initcell
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE ioncc_clean( )
!--------------------------------------------------------------------
    !
    ! ... Clean up of local variables
    !
    IMPLICIT NONE
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE ioncc_clean
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE calc_vioncc(  nnr, nspin, rhoelec, vioncc )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: nnr, nspin
    REAL( DP ), INTENT(IN) :: rhoelec(nnr)
    REAL( DP ), INTENT(OUT) :: vioncc(nnr)
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE calc_vioncc
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE calc_eioncc(  nnr, rho, eioncc )
!--------------------------------------------------------------------
    !
    USE environ_base,  ONLY : vioncc
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: nnr
    REAL( DP ), INTENT(IN) :: rho(nnr)
    REAL( DP ), INTENT(OUT) :: eioncc
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE calc_eioncc
!--------------------------------------------------------------------
!--------------------------------------------------------------------
END MODULE ioncc
!--------------------------------------------------------------------
