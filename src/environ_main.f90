! Copyright (C) 2018 ENVIRON (www.quantum-environment.org)
! Copyright (C) 2006-2010 Quantum ESPRESSO group
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
! Module containing the main drivers to compute Environ contributions
! to Kohn-Sham potential, total energy and inter-atomic forces
!
! Authors: Oliviero Andreussi (Department of Physics, UNT)
!          Francesco Nattino  (THEOS and NCCR-MARVEL, EPFL)
!          Ismaila Dabo       (DMSE, Penn State)
!          Nicola Marzari     (THEOS and NCCR-MARVEL, EPFL)
!
!----------------------------------------------------------------------------
MODULE environ_main
!----------------------------------------------------------------------------
  !
  USE environ_types
  USE electrostatic_types
  USE environ_output
  !
PRIVATE
!
PUBLIC :: calc_venviron, calc_eenviron, calc_fenviron
!
CONTAINS
!--------------------------------------------------------------------
  SUBROUTINE calc_venviron( update, nnr, vtot )
!--------------------------------------------------------------------
    !
    ! ... Calculates the environ contribution to the local
    !     potential. All the Environ modules need to be called here.
    !     The potentials are all computed on the dense real-space
    !     grid and added to vtot.
    !
    USE environ_base,  ONLY : vzero, solvent,                       &
                              lelectrostatic, velectrostatic,       &
                              vreference, lexternals,               &
                              lsoftcavity, vsoftcavity,             &
                              lsurface, env_surface_tension,        &
                              lvolume, env_pressure,                &
                              charges, lstatic, static,             &
                              lelectrolyte, electrolyte,            &
                              cell, lsoftsolvent, lsoftelectrolyte
    USE electrostatic_base, ONLY : reference, outer
    !
    ! ... Each contribution to the potential is computed in its module
    !
    USE electrostatic, ONLY : calc_velectrostatic
    USE cavity,        ONLY : calc_decavity_dboundary
    USE pressure,      ONLY : calc_depressure_dboundary
    USE dielectric,    ONLY : calc_dedielectric_dboundary
    USE electrolyte_utils, ONLY : calc_deelectrolyte_dboundary
    USE generate_boundary, ONLY : solvent_aware_de_dboundary
    USE charges_utils,     ONLY : update_environ_charges, charges_of_potential
    !
    IMPLICIT NONE
    !
    ! ... Declares variables
    !
    LOGICAL, INTENT(IN)       :: update
    INTEGER, INTENT(IN)       :: nnr
    REAL( DP ), INTENT(OUT)   :: vtot( nnr )
    !
    TYPE( environ_density ) :: de_dboundary
    !
    ! ... If not updating the potentials, add old potentials and exit
    !
    vtot = vzero % of_r
    !
    IF ( .NOT. update ) THEN
       IF ( lelectrostatic ) vtot = vtot + velectrostatic % of_r - vreference % of_r
       IF ( lsoftcavity ) vtot = vtot + vsoftcavity % of_r
       RETURN
    END IF
    !
    ! ... If any form of electrostatic embedding is present, calculate its contribution
    !
    IF ( lelectrostatic ) THEN
       !
       ! ... Electrostatics is also computed inside the calling program, need to remove the reference
       !
       CALL calc_velectrostatic( reference, charges, vreference )
       IF ( verbose .GE. 3 ) CALL print_environ_density( vreference )
       !
       IF ( lexternals ) CALL update_environ_charges( charges, lexternals )
       !
       CALL calc_velectrostatic( outer, charges, velectrostatic )
       IF ( verbose .GE. 3 ) CALL print_environ_density( velectrostatic )
       !
       vtot = vtot + velectrostatic % of_r - vreference % of_r
       !
       CALL charges_of_potential( velectrostatic, charges )
       !
       IF ( lexternals ) CALL update_environ_charges( charges )
       IF ( verbose .GE. 3 ) CALL print_environ_charges( charges )
       !
    END IF
    !
    ! ... Compute the total potential depending on the boundary
    !
    IF ( lsoftcavity ) THEN
       !
       vsoftcavity % of_r = 0.D0
       CALL init_environ_density( cell, de_dboundary )
       !
       IF ( lsoftsolvent ) THEN
          !
          de_dboundary % of_r = 0.D0
          !
          ! ... If surface tension greater than zero, calculates cavity contribution
          !
          IF ( lsurface ) CALL calc_decavity_dboundary( env_surface_tension, solvent, de_dboundary )
          !
          ! ... If external pressure different from zero, calculates PV contribution
          !
          IF ( lvolume ) CALL calc_depressure_dboundary( env_pressure, solvent, de_dboundary )
          !
          ! ... If dielectric embedding, calcultes dielectric contribution
          !
          IF ( lstatic ) CALL calc_dedielectric_dboundary( static, velectrostatic, de_dboundary )
          !
          ! ... If solvent-aware interface correct the potential
          !
          IF ( solvent % solvent_aware ) CALL solvent_aware_de_dboundary( solvent, de_dboundary )
          !
          ! ... Multiply for the derivative of the boundary wrt electronic density
          !
          vsoftcavity % of_r = de_dboundary % of_r * solvent % dscaled % of_r
          !
       END IF
       !
       IF ( lsoftelectrolyte ) THEN
          !
          de_dboundary % of_r = 0.D0
          !
          ! ... If electrolyte is present add its non-electrostatic contribution
          !
          CALL calc_deelectrolyte_dboundary( electrolyte, de_dboundary )
          !
          ! ... If solvent-aware interface correct the potential
          !
          IF ( electrolyte % boundary % solvent_aware ) CALL solvent_aware_de_dboundary( electrolyte % boundary, de_dboundary )
          !
          ! ... Multiply for the derivative of the boundary wrt electronic density
          !
          vsoftcavity % of_r = vsoftcavity % of_r + de_dboundary % of_r * electrolyte % boundary % dscaled % of_r
          !
       END IF
       !
       IF ( verbose .GE. 3 ) CALL print_environ_density( vsoftcavity )
       vtot = vtot + vsoftcavity % of_r
       !
       CALL destroy_environ_density( de_dboundary )
       !
    END IF
    !
    RETURN
!--------------------------------------------------------------------
  END SUBROUTINE calc_venviron
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE calc_eenviron( deenviron, eelectrostatic, ecavity, &
       & epressure, eelectrolyte )
!--------------------------------------------------------------------
    !
    ! ... Calculates the environ contribution to the Energy.
    !     We must remove \int v_environ * rhoelec that is
    !     automatically included in the energy computed as sum of
    !     Kohn-Sham eigenvalues.
    !
    USE environ_base,  ONLY : electrons, solvent,                   &
                              lelectrostatic, velectrostatic,       &
                              vreference, lexternals,               &
                              lsoftcavity, vsoftcavity,             &
                              lsurface, env_surface_tension,        &
                              lvolume, env_pressure,                &
                              charges, lstatic, static,             &
                              lelectrolyte, electrolyte
    USE electrostatic_base, ONLY : reference, outer
    !
    ! ... Each contribution to the energy is computed in its module
    !
    USE electrostatic, ONLY : calc_eelectrostatic
    USE cavity,        ONLY : calc_ecavity
    USE pressure,      ONLY : calc_epressure
    USE electrolyte_utils, ONLY : calc_eelectrolyte
    USE charges_utils,     ONLY : update_environ_charges
    !
    IMPLICIT NONE
    !
    ! ... Declares variables
    !
    REAL( DP ), INTENT(OUT) :: deenviron, eelectrostatic, ecavity, &
         epressure, eelectrolyte
    REAL( DP ) :: ereference
    !
    ! ... Initializes the variables
    !
    deenviron      = 0.D0
    eelectrostatic = 0.D0
    ecavity        = 0.D0
    epressure      = 0.D0
    eelectrolyte   = 0.D0
    !
    ! ... Calculates the energy corrections
    !
    ! if electrostatic is on compute electrostatic energy
    !
    IF ( lelectrostatic ) THEN
       !
       deenviron = deenviron + &
            & scalar_product_environ_density(electrons%density,vreference)
       !
       CALL calc_eelectrostatic( reference, charges, vreference, ereference )
       !
       deenviron = deenviron - &
            & scalar_product_environ_density(electrons%density,velectrostatic)
       !
       IF ( lexternals ) CALL update_environ_charges( charges, add_externals=lexternals )
       !
       CALL calc_eelectrostatic( outer, charges, velectrostatic, eelectrostatic )
       !
       eelectrostatic = eelectrostatic - ereference
       !
       IF ( lexternals ) CALL update_environ_charges( charges )
       !
    END IF
    !
    IF ( lsoftcavity ) deenviron = deenviron - &
         & scalar_product_environ_density(electrons%density,vsoftcavity)
    !
    !  if surface tension different from zero compute cavitation energy
    !
    IF ( lsurface ) CALL calc_ecavity( env_surface_tension, solvent, ecavity )
    !
    !  if pressure different from zero compute PV energy
    !
    IF ( lvolume ) CALL calc_epressure( env_pressure, solvent, epressure )
    !
    !  if electrolyte is present calculate its non-electrostatic contribution
    !
    IF ( lelectrolyte ) CALL calc_eelectrolyte( electrolyte, eelectrolyte )
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE calc_eenviron
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE calc_fenviron( nat, force_environ )
!--------------------------------------------------------------------
    !
    ! ... Calculates the environ contribution to the forces.
    !     Due to Hellman-Feynman only a few of Environ modules
    !     have an effect on atomic forces.
    !
    USE environ_base, ONLY : lelectrostatic, velectrostatic,    &
                             charges, lstatic, static,          &
                             lelectrolyte, electrolyte,         &
                             lrigidcavity, lrigidsolvent,       &
                             lrigidelectrolyte,                 &
                             lsurface, env_surface_tension,     &
                             lvolume, env_pressure,             &
                             lsolvent, solvent, cell
    !
    USE electrostatic_base, ONLY : outer
    !
    ! ... Each contribution to the forces is computed in its module
    !
    USE electrostatic,     ONLY : calc_felectrostatic
    USE cavity,            ONLY : calc_decavity_dboundary
    USE pressure,          ONLY : calc_depressure_dboundary
    USE dielectric,        ONLY : calc_dedielectric_dboundary
    USE electrolyte_utils, ONLY : calc_deelectrolyte_dboundary
    USE generate_boundary, ONLY : calc_dboundary_dions, solvent_aware_de_dboundary
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: nat
    REAL( DP ), INTENT(INOUT) :: force_environ( 3, nat )
    !
    INTEGER :: i
    TYPE( environ_density ) :: de_dboundary
    TYPE( environ_gradient ) :: partial
    !
    force_environ = 0.D0
    !
    IF ( lelectrostatic ) CALL calc_felectrostatic( outer, nat, charges, force_environ )
    !
    ! ... Compute the total forces depending on the boundary
    !
    IF ( lrigidcavity ) THEN
       !
       CALL init_environ_density( cell, de_dboundary )
       !
       CALL init_environ_gradient( cell, partial )
       !
       IF ( lrigidsolvent ) THEN
          !
          de_dboundary % of_r = 0.D0
          !
          ! ... If surface tension greater than zero, calculates cavity contribution
          !
          IF ( lsurface ) CALL calc_decavity_dboundary( env_surface_tension, solvent, de_dboundary )
          !
          ! ... If external pressure different from zero, calculates PV contribution
          !
          IF ( lvolume ) CALL calc_depressure_dboundary( env_pressure, solvent, de_dboundary )
          !
          ! ... If dielectric embedding, calcultes dielectric contribution
          !
          IF ( lstatic ) CALL calc_dedielectric_dboundary( static, velectrostatic, de_dboundary )
          !
          ! ... If solvent-aware interface correct the potential
          !
          IF ( solvent % solvent_aware ) CALL solvent_aware_de_dboundary( solvent, de_dboundary )
          !
          ! ... Multiply for the derivative of the boundary wrt ionic positions
          !
          DO i = 1, nat
             !
             CALL calc_dboundary_dions( i, solvent, partial )
             !
             force_environ( :, i ) = force_environ( :, i ) &
                  & - scalar_product_environ_gradient_density( partial, de_dboundary )
             !
          END DO
          !
       END IF
       !
       IF ( lrigidelectrolyte ) THEN
          !
          de_dboundary % of_r = 0.D0
          !
          ! ... If electrolyte is present, add its non-electrostatic contribution
          !
          CALL calc_deelectrolyte_dboundary( electrolyte, de_dboundary )
          !
          ! ... If solvent-aware interface correct the potential
          !
          IF ( electrolyte % boundary % solvent_aware ) CALL solvent_aware_de_dboundary( electrolyte % boundary, de_dboundary )
          !
          ! ... Multiply for the derivative of the boundary wrt ionic positions
          !
          DO i = 1, nat
             !
             CALL calc_dboundary_dions( i, electrolyte % boundary, partial )
             !
             force_environ( :, i ) = force_environ( :, i ) &
                  & - scalar_product_environ_gradient_density( partial, de_dboundary )
             !
          END DO
          !
       END IF
       !
       CALL destroy_environ_gradient( partial )
       !
       CALL destroy_environ_density( de_dboundary )
       !
    END IF
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE calc_fenviron
!--------------------------------------------------------------------
!----------------------------------------------------------------------------
END MODULE environ_main
!----------------------------------------------------------------------------
