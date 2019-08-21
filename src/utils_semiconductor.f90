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
  SUBROUTINE create_environ_semiconductor(semiconductor)
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    TYPE( environ_semiconductor ), INTENT(INOUT) :: semiconductor
    CHARACTER ( LEN=80 ) :: sub_name = 'create_environ_semiconductor'
    CHARACTER ( LEN=80 ) :: label = 'semiconductor'
    !
    semiconductor%update = .FALSE.

    !semiconductor%temperature = 300
    !CALL create_environ_function( semiconductor%simple, label )
    CALL create_environ_density( semiconductor%density, label )
    semiconductor%charge = 0.D0
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE create_environ_semiconductor
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE init_environ_semiconductor_first( temperature,   &
      & sc_permittivity, sc_carrier_density , sc_electrode_chg, sc_distance, sc_spread, system, &
      & env_dielectric_regions, epsregion_eps, epsregion_dim, epsregion_axis,  &
      & epsregion_pos, epsregion_spread, epsregion_width, nat,env_static_permittivity,semiconductor)
!--------------------------------------------------------------------
    !
    IMPLICIT NONE

    USE environ_base, ONLY: cell, ions
    !
    REAL( DP ), INTENT(IN) :: temperature, sc_permittivity, sc_electrode_chg
    REAL( DP ), INTENT(IN) :: sc_carrier_density, sc_distance, sc_spread
    INTEGER,    INTENT(IN) :: nat
    INTEGER,    INTENT(INOUT) :: env_dielectric_regions
    INTEGER, DIMENSION(env_dielectric_regions+1), INTENT(INOUT) :: epsregion_dim, epsregion_axis
    REAL( DP ), DIMENSION(env_dielectric_regions+1), INTENT(INOUT) :: epsregion_width, epsregion_spread, epsregion_eps
    REAL( DP ), DIMENSION(3,env_dielectric_regions+1), INTENT(INOUT) :: epsregion_pos
    REAL( DP ) :: at_max, at_min, tot_len
    TYPE( environ_system ), TARGET, INTENT(IN) :: system
    TYPE( environ_semiconductor ), INTENT(INOUT) :: semiconductor
    !
    semiconductor%temperature = temperature
    semiconductor%permittivity = sc_permittivity
    semiconductor%carrier_density = sc_carrier_density

    !   convert carrier density to units of (bohr)^-3
    semiconductor%carrier_density = semiconductor%carrier_density *1.25D-25


    ! Possibly should have an if statement to check if electrode charge is
    ! greater than tot_charge on DFT, need to consider this option
    semiconductor%electrode_charge = sc_electrode_chg

    semiconductor%simple%type = 4
    semiconductor%simple%pos => system%pos
    semiconductor%simple%volume = 1.D0
    semiconductor%simple%dim = system%dim
    semiconductor%simple%axis = system%axis
    semiconductor%simple%width = sc_distance
    semiconductor%simple%spread = sc_spread


    semiconductor%initialized = .FALSE.

    ! Adding a dielectric region with the permittivity of the semiconductor
    ! to the side of the slab that wil be experiencing a mott schottky correction
    env_dielectric_regions = env_dielectric_regions + 1
    epsregion_eps(env_dielectric_regions) = semiconductor%permittivity
    epsregion_dim(env_dielectric_regions) = 2
    epsregion_axis(env_dielectric_regions) = semiconductor%simple%axis

    at_max = -0.4
    at_min = 10000.D0
    DO ipol = 1, nat
       IF (system%ions%tau(semiconductor%simple%axis) > at_max ) &
            & at_max = system%ions%tau(semiconductor%simple%axis)
       IF (system%ions%tau(semiconductor%simple%axis) < at_min ) &
            & at_min = system%ions%tau(semiconductor%simple%axis)
    END DO

    tot_len = cell%alat*cell%at(semiconductor%simple%axis,semiconductor%simple%axis)
    epsregion_width(env_dielectric_regions) = ( tot_len- (at_max-at_min))/2.D0

    IF ((at_max + epsregion_width(env_dielectric_regions)) >  tot_len) THEN
        epsregion_center(env_dielectric_regions) = at_min -epsregion_width(env_dielectric_regions)
    ELSE
        epsregion_center(env_dielectric_regions) = at_max +epsregion_width(env_dielectric_regions)
    END IF

    epsregion_spread(env_dielectric_regions) = sc_spread


    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE init_environ_semiconductor_first
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE init_environ_semiconductor_second( cell, semiconductor )
!--------------------------------------------------------------------
    !
    !IMPLICIT NONE
    !
    TYPE( environ_cell ), INTENT(IN) :: cell
    TYPE( environ_semiconductor ), INTENT(INOUT) :: semiconductor
    !
    !INTEGER :: i
    !
    !IF ( externals % number .GT. 0 ) THEN
    !   DO i = 1, externals % number
    !      externals % functions(i) % pos = externals % functions(i) % pos / cell % alat
    !   END DO
    !END IF
    !
    CALL init_environ_density( cell, semiconductor%density )
    semiconductor%initialized = .TRUE.
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE init_environ_semiconductor_second
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE update_environ_semiconductor( semiconductor )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    TYPE( environ_semiconductor ), INTENT(INOUT) :: semiconductor
    !
    !CALL density_of_functions( externals%number, externals%functions, externals%density, .TRUE. )
    !
    semiconductor % charge = integrate_environ_density( semiconductor % density )
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE update_environ_semiconductor
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE destroy_environ_semiconductor( lflag, semiconductor )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    LOGICAL, INTENT(IN) :: lflag
    TYPE( environ_semiconductor ), INTENT(INOUT) :: semiconductor
    !
    !IF ( lflag ) CALL destroy_environ_functions( externals%number, externals%functions )
    !
    IF ( semiconductor%initialized ) THEN
      CALL destroy_environ_density( semiconductor%density )
      semiconductor%initialized = .FALSE.
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
