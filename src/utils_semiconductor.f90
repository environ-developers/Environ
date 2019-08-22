! Copyright (C) 2018 ENVIRON (www.quantum-environment.org)
!
!    This file is part of Environ version 1.0
!
!    Environ 1.0 is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 2 of the License, or
!    (at your option) any later version.
!
!    Environ 1.0 is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more detail, either the file
!    `License' in the root directory of the present distribution, or
!    online at <http://www.gnu.org/licenses/>.
!
! Module containing the main routines to handle
!
!              environ_semiconductor
!
! derived data types.
!
! Environ_semiconductor contains all the specifications and the details of
! the user defined semiconductor region
!
!  TYPE environ_semiconductor
     !
     ! Update status
     !
!     LOGICAL :: update = .FALSE.
     !
!     LOGICAL :: initialized = .FALSE.
     !
     !
!     REAL( DP ) :: temperature
!     REAL( DP ) :: sc_permittivity
!     REAL( DP ) :: sc_carrier_density
     !
     ! As far as I can tell this is not relevant for the semicondutor
     !TYPE( environ_function ) :: simple
     !
!     TYPE( environ_density ) :: density
     !
     ! The electrolyte switch function and relate quantities
     !
     !TYPE( environ_density ) :: gamma
     !TYPE( environ_density ) :: dgamma
     !
     !TYPE( environ_density ) :: de_dboundary_second_order
     !REAL( DP ) :: energy_second_order
     !
!     REAL( DP ) :: charge = 0.0_DP
     !
!  END TYPE environ_semiconductor
!
!
! Authors: Quinn Campbell (Department of Materials Science and Engineering, Penn State)
!
!----------------------------------------------------------------------------
MODULE utils_semiconductor
!----------------------------------------------------------------------------
  !
  USE environ_types
  USE environ_base, ONLY: semiconductor
  USE environ_output
  USE utils_functions
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC :: create_environ_semiconductor, init_environ_semiconductor_first, &
       & init_environ_semiconductor_second, update_environ_semiconductor, &
       & destroy_environ_semiconductor
  !
CONTAINS
!--------------------------------------------------------------------
  SUBROUTINE create_environ_semiconductor(semiconductor_in)
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    TYPE( environ_semiconductor ), INTENT(INOUT) :: semiconductor_in
    CHARACTER ( LEN=80 ) :: sub_name = 'create_environ_semiconductor'
    CHARACTER ( LEN=80 ) :: label = 'semiconductor'
    !
    semiconductor_in%update = .FALSE.

    !semiconductor%temperature = 300
    !CALL create_environ_function( semiconductor%simple, label )
    CALL create_environ_density( semiconductor_in%density, label )
    semiconductor_in%charge = 0.D0
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE create_environ_semiconductor
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE init_environ_semiconductor_first( temperature, &
       & sc_permittivity,sc_carrier_density, sc_electrode_chg, sc_distance, &
       & sc_spread, sc_chg_thr, system, semiconductor_in )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    REAL( DP ), INTENT(IN) :: temperature, sc_permittivity, sc_electrode_chg
    REAL( DP ), INTENT(IN) :: sc_carrier_density, sc_distance, sc_spread, sc_chg_thr
    TYPE( environ_system ), TARGET, INTENT(IN) :: system
    TYPE( environ_semiconductor ), INTENT(INOUT) :: semiconductor_in
    !
    semiconductor_in%temperature = temperature
    semiconductor_in%permittivity = sc_permittivity
    semiconductor_in%carrier_density = sc_carrier_density

    !   convert carrier density to units of (bohr)^-3
    semiconductor_in%carrier_density = semiconductor_in%carrier_density *1.25D-25


    ! Possibly should have an if statement to check if electrode charge is
    ! greater than tot_charge on DFT, need to consider this option
    semiconductor_in%electrode_charge = sc_electrode_chg
    semiconductor_in%charge_threshold = sc_chg_thr

    semiconductor_in%simple%type = 4
    semiconductor_in%simple%pos => system%pos
    semiconductor_in%simple%volume = 1.D0
    semiconductor_in%simple%dim = system%dim
    semiconductor_in%simple%axis = system%axis
    semiconductor_in%simple%width = sc_distance
    semiconductor_in%simple%spread = sc_spread


    semiconductor_in%initialized = .FALSE.
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE init_environ_semiconductor_first
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE init_environ_semiconductor_second( cell, semiconductor_in )
!--------------------------------------------------------------------
    !
    !IMPLICIT NONE
    !
    TYPE( environ_cell ), INTENT(IN) :: cell
    TYPE( environ_semiconductor ), INTENT(INOUT) :: semiconductor_in
    !
    !INTEGER :: i
    !
    !IF ( externals % number .GT. 0 ) THEN
    !   DO i = 1, externals % number
    !      externals % functions(i) % pos = externals % functions(i) % pos / cell % alat
    !   END DO
    !END IF
    !
    CALL init_environ_density( cell, semiconductor_in%density )
    semiconductor_in%initialized = .TRUE.
    semiconductor = semiconductor_in
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE init_environ_semiconductor_second
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE update_environ_semiconductor( semiconductor_in )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    TYPE( environ_semiconductor ), INTENT(INOUT) :: semiconductor_in
    !
    !CALL density_of_functions( externals%number, externals%functions, externals%density, .TRUE. )
    !
    semiconductor_in % charge = integrate_environ_density( semiconductor_in % density )
    semiconductor = semiconductor_in
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE update_environ_semiconductor
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE destroy_environ_semiconductor( lflag, semiconductor_in )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    LOGICAL, INTENT(IN) :: lflag
    TYPE( environ_semiconductor ), INTENT(INOUT) :: semiconductor_in
    !
    !IF ( lflag ) CALL destroy_environ_functions( externals%number, externals%functions )
    !
    IF ( semiconductor_in%initialized ) THEN
      CALL destroy_environ_density( semiconductor_in%density )
      semiconductor_in%initialized = .FALSE.
    END IF
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE destroy_environ_semiconductor
!--------------------------------------------------------------------
!----------------------------------------------------------------------------
END MODULE utils_semiconductor
!----------------------------------------------------------------------------
