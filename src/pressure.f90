!
! Copyright (C) 2007-2008 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Module to include an environment pressure potential, proportional to
! the quantum volume of the system. The method was developed by
! Cococcioni et al. PRL 94, 145501 (2005).
!
! original version by O. Andreussi, N. Marzari (MIT)
!
!----------------------------------------------------------------------------
MODULE pressure
!----------------------------------------------------------------------------
  !
  ! ... The variables needed to compute PV energy and potential corrections
  !
  USE environ_types
  USE environ_output
  USE environ_base,  ONLY : env_pressure, e2
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC :: calc_vpressure, calc_epressure, calc_fpressure
  !
CONTAINS
!--------------------------------------------------------------------
  SUBROUTINE calc_vpressure( boundary, potential )
!--------------------------------------------------------------------
    !
    ! ... Calculates the PV contribution to the potential
    !
    IMPLICIT NONE
    !
    ! ... Declares variables
    !
    TYPE( environ_boundary ), INTENT(IN) :: boundary
    TYPE( environ_density ), INTENT(INOUT) :: potential
    !
    CHARACTER( LEN=80 ) :: sub_name = 'calc_vpressure'
    !
    CALL start_clock ('calc_vpre')
    !
    IF ( boundary%need_electrons ) THEN
       !
       ! ... Multiply the derivative of the volume by the external pressure
       !
       potential%of_r = env_pressure * boundary%dscaled%of_r
       !
    ENDIF
    !
    CALL stop_clock ('calc_vpre')
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE calc_vpressure
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE calc_epressure( boundary, epressure )
!--------------------------------------------------------------------
    !
    ! ... Calculates the PV contribution to the energy
    !
    IMPLICIT NONE
    !
    ! ... Declares variables
    !
    TYPE( environ_boundary ), TARGET, INTENT(IN) :: boundary
    REAL( DP ), INTENT(OUT) :: epressure
    !
    REAL( DP )              :: volume
    TYPE( environ_density ) :: volofrho
    !
    CALL start_clock ('calc_epre')
    !
    ! ... Computes step function and volume
    !
    volume = integrate_environ_density( boundary % scaled )
    !
    ! ... Computes the PV energy
    !
    epressure = env_pressure * volume * e2 / 2.D0
    !
    CALL stop_clock ('calc_epre')
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE calc_epressure
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE calc_fpressure( nat, boundary, forces )
!--------------------------------------------------------------------
    !
    ! ... Calculates the cavitation contribution to the energy
    !
    IMPLICIT NONE
    !
    ! ... Declares variables
    !
    INTEGER, INTENT(IN) :: nat
    TYPE( environ_boundary ), INTENT(IN) :: boundary
    REAL( DP ), DIMENSION(3,nat), INTENT(INOUT) :: forces
    !
    REAL( DP ), DIMENSION(3,nat) :: fpressure
    !
    CHARACTER( LEN=80 ) :: sub_name = 'calc_fpressure'
    !
    CALL start_clock ('calc_fpre')
    !
    fpressure = 0.D0
    !
    ! ... Computes
    !
    IF ( boundary%need_ions ) THEN
       !
       ! ... Rigid cases to be implemented
       !
       CALL errore(sub_name,'Option not yet implemented',1)
       !
    ENDIF
    !
    CALL stop_clock ('calc_fpre')
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE calc_fpressure
!--------------------------------------------------------------------
!--------------------------------------------------------------------
END MODULE pressure
!--------------------------------------------------------------------
