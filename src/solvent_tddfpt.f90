!
! Copyright (C) 2001-2013 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE solvent_tddfpt
!----------------------------------------------------------------------------
  !
  ! ... Module for calculation of the response "polarization" and "dielectric" potentials.
  ! ... Inspired by Environ/src/solvent.f90
  ! ... Written by I. Timrov 09/2013
  ! ... Re-written by O. Andreussi 11/2017
  !
  USE environ_type
  USE environ_output
  USE environ_base,       ONLY : e2
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
!--------------------------------------------------------------------
 SUBROUTINE solvent_clean_tddfpt()
!--------------------------------------------------------------------
   !
   ! ... Local clean up
   !
   IMPLICIT NONE
   !
   CALL environ_clean_tddfpt(.TRUE.)
   !
   RETURN
   !
!---------------------------------------------------------------------------------
 END SUBROUTINE solvent_clean_tddfpt
!---------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
 SUBROUTINE calc_vsolvent_tddfpt(nnr, nspin, rho_0, drho_elec, dv_pol, dv_epsilon)
!-----------------------------------------------------------------------------------
   !
   ! ... This subroutine calculates:
   ! ... 1. the response "polarization" potential from the response polarization density
   ! ... 2. the response "dielectric" potential
   !
   USE environ_base, ONLY : cell, velectrostatic, lsoftcavity
   USE electrostatic_base, ONLY : reference, outer
   !
   IMPLICIT NONE
   !
   INTEGER, INTENT(IN)     :: nnr,           & ! number of grid points in R-space
                              nspin            ! if nspin=2 spin-polarized case (not supported)
   REAL( DP ), INTENT(IN)  :: rho_0(nnr),    & ! ground-state charge-density
                              drho_elec(nnr)   ! response charge-density
   REAL( DP ), INTENT(OUT) :: dv_pol(nnr),   & ! response polarization potential
                              dv_epsilon(nnr)  ! response dielectric potential
   !
   ! ... Local variables
   !
   TYPE( environ_density ) :: response_density
   TYPE( environ_density ) :: dvreference
   TYPE( environ_density ) :: dvelectrostatic
   !
   INTEGER :: ir
   REAL( DP ), DIMENSION( :, : ), ALLOCATABLE :: gvtot0, gdvtot
   !
   IF ( .NOT. tddfpt ) RETURN
   !
   CALL start_clock( 'calc_vsolvent_tddfpt' )
   !
   ! ... Sanity checks and local allocation
   !
   IF ( nnr .NE. cell % nnr ) CALL errore(sub_name,'Missmatch in passed and stored grid dimension',1)
   !
   CALL init_environ_density( cell, response_density )
   !
   response_density % of_r = drho_elec
   !
   ! ... Compute reference potential of response density
   !
   CALL init_environ_density( cell, dvreference )
   !
   CALL calc_vreference( reference, response_density, dvreference )
   !
   ! ... Compute full electrostatic potential of response density
   !
   CALL init_environ_density( cell, dvelectrostatic )
   !
   CALL calc_velectrostatic( setup = outer, charges=response_density, dielectric=optical,  &
        & potential=dvelectrostatic )
   !
   ! ... Response polarization potential
   !
   dvpol = dvelectrostatic % of_r - dvreference % of_r
   !
   IF ( lsoftcavity ) THEN
      !
      ! ... Calculate a gradient of the total static potential
      !
      ALLOCATE( gvtot0( 3, nnr ) )
      !
      CALL external_gradient( velectrostatic % of_r, gvtot0 )
      !
      ! ... Calculate a gradient of the total response potential
      !
      ALLOCATE( gdvtot( 3, nnr ) )
      !
      CALL external_gradient( dvelectrostatic % of_r, gdvtot )
      !
      ! ... Response dielectric potential
      !
      DO ir = 1, nnr
         dv_epsilon(ir) = - SUM( gvtot0(:,ir) * gdvtot(:,ir) ) * &
              & optical % depsilon % of_r(ir) / (fpi * e2)
      END DO
      !
      DEALLOCATE( gdvtot )
      DEALLOCATE( gvtot0 )
      !
   ELSE
      !
      dv_epsilon = 0.D0
      !
   END IF
   !
   CALL destroy_environ_density( dvelectrostatic )
   CALL destroy_environ_density( dvreference )
   CALL destroy_environ_density( response_density )
   !
   CALL stop_clock( 'calc_vsolvent_tddfpt' )
   !
   RETURN
   !
!--------------------------------------------------------------------
END SUBROUTINE calc_vsolvent_tddfpt
!--------------------------------------------------------------------
!--------------------------------------------------------------------
END MODULE solvent_tddfpt
!--------------------------------------------------------------------
