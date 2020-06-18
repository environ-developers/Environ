! Copyright (C) 2018 ENVIRON (www.quantum-environment.org)
!
!    This file is part of Environ version 1.1
!
!    Environ 1.1 is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 2 of the License, or
!    (at your option) any later version.
!
!    Environ 1.1 is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more detail, either the file
!    `License' in the root directory of the present distribution, or
!    online at <http://www.gnu.org/licenses/>.
!
!> Module for calculation of the response "polarization" and "dielectric"
!! potentials, to be coupled with the TDDFPT algorithms of Quantum ESPRESSO.
!! Original formulas derived in I. Timrov et al, JCP (2015)
!
! Authors: Oliviero Andreussi (Department of Physics, UNT)
!          Iurii Timrov       (THEOS and NCCR-MARVEL, EPFL)
!          Nicola Marzari     (THEOS and NCCR-MARVEL, EPFL)
!
! Original version by I. Timrov and S. Baroni (SISSA), 09/2013
!
!----------------------------------------------------------------------------
MODULE solvent_tddfpt
!----------------------------------------------------------------------------
  !
  USE environ_types
  USE environ_output
  USE modules_constants, ONLY : e2, fpi
  !
  IMPLICIT NONE
  !
  SAVE
  !
  PRIVATE
  !
  PUBLIC :: solvent_clean_tddfpt, calc_vsolvent_tddfpt
  !
CONTAINS
!  Subroutine: solvent_clean_tddfpt
!
!> Local clean up
!--------------------------------------------------------------------
 SUBROUTINE solvent_clean_tddfpt()
!--------------------------------------------------------------------
   USE environ_init, ONLY : environ_clean_tddfpt
   !
   IMPLICIT NONE
   !
   CALL environ_clean_tddfpt(.TRUE.)
   !
   RETURN
   !
!--------------------------------------------------------------------
 END SUBROUTINE solvent_clean_tddfpt
!--------------------------------------------------------------------
!  Subroutine: calc_vsolvent_tddfpt
!
!> This subroutine calculates:
!! -# the response "polarization" from the response polarization
!! density
!! -# the response "dielectric" potential
!--------------------------------------------------------------------
! BACKWARD COMPATIBILITY
! Compatible with QE-6.0 QE-6.1.X QE-6.2.X QE-6.3
! SUBROUTINE calc_vsolvent_tddfpt(nnr, nspin, rho_0, drho_elec, dv_pol, dv_epsilon)
! Compatible with QE-6.4.X QE-GIT
 SUBROUTINE calc_vsolvent_tddfpt(nnr, rho_0, drho_elec, dv_pol, dv_epsilon)
! END BACKWARD COMPATIBILITY
!--------------------------------------------------------------------
   USE environ_base, ONLY: sys_cell, velectrostatic, lsoftsolvent, loptical, optical, ltddfpt
   USE electrostatic_base, ONLY : reference, outer
   USE embedding_electrostatic, ONLY : calc_velectrostatic
   USE core_fft, ONLY : gradient_fft
   USE utils_charges
   !
   IMPLICIT NONE
   !
   INTEGER, INTENT(IN)     :: nnr              ! number of grid points in R-space
! BACKWARD COMPATIBILITY
! Compatible with QE-6.0 QE-6.1.X QE-6.2.X QE-6.3
!   INTEGER, INTENT(IN)     :: nspin            ! if nspin=2 spin-polarized case (not supported)
! Compatible with QE-6.4.X QE-GIT
!
! END BACKWARD COMPATIBILITY
   REAL( DP ), INTENT(IN)  :: rho_0(nnr),    & ! ground-state charge-density
                              drho_elec(nnr)   ! response charge-density
   REAL( DP ), INTENT(OUT) :: dv_pol(nnr),   & ! response polarization potential
                              dv_epsilon(nnr)  ! response dielectric potential
   !
   ! ... Local variables
   !
   TYPE( environ_charges ) :: response_charges
   TYPE( environ_electrons ) :: response_electrons
   TYPE( environ_density ) :: dvreference
   TYPE( environ_density ) :: dvelectrostatic
   !
   INTEGER :: ir
   TYPE( environ_gradient ) :: gvtot0, gdvtot
   CHARACTER( LEN=80 ) :: sub_name = 'calc_vsolvent_tddfpt'
   !
   IF ( .NOT. ltddfpt .OR. .NOT. loptical ) RETURN
   !
   CALL start_clock( 'calc_vsolvent_tddfpt' )
   !
   ! ... Sanity checks and local allocation
   !
   IF ( nnr .NE. sys_cell % nnr ) CALL errore(sub_name,'Missmatch in passed and stored grid dimension',1)
   !
   ! ... Create source response electronic density
   !
   CALL create_environ_electrons( response_electrons )
! BACKWARD COMPATIBILITY
! Compatible with QE-6.0 QE-6.1.X QE-6.2.X QE-6.3.X
!   CALL init_environ_electrons_first( 0, 1, response_electrons )
!   CALL init_environ_electrons_second( cell, response_electrons )
!   CALL update_environ_electrons( 1, nnr, drho_elec, response_electrons, 0.D0 )
! Compatible with QE-6.4.X QE-GIT
   CALL init_environ_electrons_first( 0, response_electrons )
   CALL init_environ_electrons_second( sys_cell, response_electrons )
   CALL update_environ_electrons( nnr, drho_elec, response_electrons, 0.D0 )
! END BACKWARD COMPATIBILITY
   !
   ! ... Link together different sources of electrostatic potential ( charges + dielectric + electrolyte )
   !
   CALL create_environ_charges( response_charges )
   CALL init_environ_charges_first( electrons=response_electrons, dielectric=optical, charges=response_charges )
   CALL init_environ_charges_second( sys_cell, response_charges )
   CALL update_environ_charges( response_charges )
   !
   ! ... Compute reference potential of response density
   !
   CALL init_environ_density( sys_cell, dvreference )
   !
   CALL calc_velectrostatic( reference, response_charges, dvreference )
   !
   ! ... Compute full electrostatic potential of response density
   !
   CALL init_environ_density( sys_cell, dvelectrostatic )
   !
   CALL calc_velectrostatic( outer, response_charges, dvelectrostatic )
   !
   ! ... Response polarization potential
   !
   dv_pol = dvelectrostatic % of_r - dvreference % of_r
   !
   IF ( lsoftsolvent ) THEN
      !
      ! ... Calculate the gradient of the total static potential
      !
      CALL init_environ_gradient( sys_cell, gvtot0 )
      !
      CALL gradient_fft( outer % core % fft, velectrostatic, gvtot0 )
      !
      ! ... Calculate the gradient of the total response potential
      !
      CALL init_environ_gradient( sys_cell, gdvtot )
      !
      CALL gradient_fft( outer % core % fft, dvelectrostatic, gdvtot )
      !
      ! ... Response dielectric potential
      !
      DO ir = 1, nnr
         dv_epsilon(ir) = - SUM( gvtot0%of_r(:,ir) * gdvtot%of_r(:,ir) ) * &
              & optical % depsilon % of_r(ir) * optical % boundary % dscaled % of_r(ir) / (fpi * e2)
      END DO
      !
      CALL destroy_environ_gradient( gdvtot )
      CALL destroy_environ_gradient( gvtot0 )
      !
   ELSE
      !
      dv_epsilon = 0.D0
      !
   END IF
   !
   CALL destroy_environ_density( dvelectrostatic )
   CALL destroy_environ_density( dvreference )
   CALL destroy_environ_charges( .TRUE., response_charges )
   CALL destroy_environ_electrons( .TRUE., response_electrons )
   !
   CALL stop_clock( 'calc_vsolvent_tddfpt' )
   !
   RETURN
   !
!--------------------------------------------------------------------
END SUBROUTINE calc_vsolvent_tddfpt
!--------------------------------------------------------------------
!----------------------------------------------------------------------------
END MODULE solvent_tddfpt
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
MODULE environ_info
!----------------------------------------------------------------------------
  USE environ_output
!----------------------------------------------------------------------------
END MODULE environ_info
!----------------------------------------------------------------------------
