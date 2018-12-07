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
!> Module containing the main drivers to compute Environ contributions
!! to Kohn-Sham potential, total energy and inter-atomic forces
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
PUBLIC :: calc_venviron, calc_eenviron, calc_fenviron, calc_dvenviron
!
CONTAINS
!  Subroutine: calc_venviron
!
!> Calculates the Environ contribution to the local potential. All
!! the Environ modules need to be called here. The potentials are
!! all computed on the dense real-space grid and added to vtot.
!--------------------------------------------------------------------
  SUBROUTINE calc_venviron( update, nnr, vtot )
!--------------------------------------------------------------------
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
    USE embedding_electrostatic, ONLY : calc_velectrostatic
    USE embedding_surface,       ONLY : calc_desurface_dboundary
    USE embedding_volume,        ONLY : calc_devolume_dboundary
    USE utils_dielectric,        ONLY : calc_dedielectric_dboundary
    USE utils_electrolyte,       ONLY : calc_deelectrolyte_dboundary
    USE utils_charges,           ONLY : update_environ_charges, &
                                      & charges_of_potential
    USE tools_generate_boundary, ONLY : solvent_aware_de_dboundary, &
                                      & field_aware_dboundary_drho
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
       IF ( verbose .GE. 2 ) CALL print_environ_density( vreference )
       !
       IF ( lexternals ) CALL update_environ_charges( charges, lexternals )
       !
       CALL calc_velectrostatic( outer, charges, velectrostatic )
       IF ( verbose .GE. 2 ) CALL print_environ_density( velectrostatic )
       !
       vtot = vtot + velectrostatic % of_r - vreference % of_r
       !
       CALL charges_of_potential( velectrostatic, charges )
       !
       IF ( lexternals ) CALL update_environ_charges( charges )
       IF ( verbose .GE. 2 ) CALL print_environ_charges( charges )
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
          IF ( lsurface ) CALL calc_desurface_dboundary( env_surface_tension, solvent, de_dboundary )
          !
          ! ... If external pressure different from zero, calculates PV contribution
          !
          IF ( lvolume ) CALL calc_devolume_dboundary( env_pressure, solvent, de_dboundary )
          !
          ! ... If dielectric embedding, calcultes dielectric contribution
          !
          IF ( lstatic ) CALL calc_dedielectric_dboundary( static, velectrostatic, de_dboundary )
          !
          ! ... If solvent-aware interface correct the potential
          !
          IF ( solvent % solvent_aware ) CALL solvent_aware_de_dboundary( solvent, de_dboundary )
          !
          ! ... If field-aware interface correct the derivative of the interface function
          !
          IF ( solvent % field_aware ) CALL field_aware_dboundary_drho( solvent, solvent%dscaled )
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
       IF ( verbose .GE. 4 ) CALL print_environ_density( vsoftcavity )
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
!  Subroutine: calc_eenviron
!
!> Calculates the Environ contribution to the energy. We must remove
!! int v_environ * rhoelec that is automatically included in the
!! energy computed as the sum of Kohn-Sham eigenvalues.
!--------------------------------------------------------------------
  SUBROUTINE calc_eenviron( deenviron, eelectrostatic, esurface, &
       & evolume, eelectrolyte )
!--------------------------------------------------------------------
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
    USE embedding_electrostatic, ONLY : calc_eelectrostatic
    USE embedding_surface,       ONLY : calc_esurface
    USE embedding_volume,        ONLY : calc_evolume
    USE utils_electrolyte,       ONLY : calc_eelectrolyte
    USE utils_charges,           ONLY : update_environ_charges
    !
    IMPLICIT NONE
    !
    ! ... Declares variables
    !
    REAL( DP ), INTENT(OUT) :: deenviron, eelectrostatic, esurface, &
         evolume, eelectrolyte
    REAL( DP ) :: ereference
    !
    ! ... Initializes the variables
    !
    deenviron      = 0.D0
    eelectrostatic = 0.D0
    esurface       = 0.D0
    evolume        = 0.D0
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
       CALL calc_eelectrostatic( reference%core, charges, vreference, ereference )
       !
       deenviron = deenviron - &
            & scalar_product_environ_density(electrons%density,velectrostatic)
       !
!       IF ( lexternals ) CALL update_environ_charges( charges, add_externals=lexternals )
!       !
       CALL calc_eelectrostatic( outer%core, charges, velectrostatic, eelectrostatic, &
                                 add_environment=.true. )
       !
       eelectrostatic = eelectrostatic - ereference
       !
!       IF ( lexternals ) CALL update_environ_charges( charges )
!       !
    END IF
    !
    IF ( lsoftcavity ) deenviron = deenviron - &
         & scalar_product_environ_density(electrons%density,vsoftcavity)
    !
    !  if surface tension different from zero compute cavitation energy
    !
    IF ( lsurface ) CALL calc_esurface( env_surface_tension, solvent, esurface )
    !
    !  if pressure different from zero compute PV energy
    !
    IF ( lvolume ) CALL calc_evolume( env_pressure, solvent, evolume )
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
!  Subroutine: calc_fenviron
!
!> Calculates the Environ contribution to the forces. Due to
!! Hellman-Feynman only a few of the Environ modules have an
!! effect on the atomic forces.
!--------------------------------------------------------------------
  SUBROUTINE calc_fenviron( nat, force_environ )
!--------------------------------------------------------------------
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
    USE embedding_electrostatic, ONLY : calc_felectrostatic
    USE embedding_surface,       ONLY : calc_desurface_dboundary
    USE embedding_volume,        ONLY : calc_devolume_dboundary
    USE utils_dielectric,        ONLY : calc_dedielectric_dboundary
    USE utils_electrolyte,       ONLY : calc_deelectrolyte_dboundary
    USE tools_generate_boundary, ONLY : calc_dboundary_dions, solvent_aware_de_dboundary, &
         & field_aware_dboundary_dions, compute_ion_field_partial
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
          IF ( lsurface ) CALL calc_desurface_dboundary( env_surface_tension, solvent, de_dboundary )
          !
          ! ... If external pressure different from zero, calculates PV contribution
          !
          IF ( lvolume ) CALL calc_devolume_dboundary( env_pressure, solvent, de_dboundary )
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
          IF ( solvent % field_aware ) CALL compute_ion_field_partial( solvent%ions%number, solvent%soft_spheres, &
               & solvent%ions, solvent%electrons, solvent%ion_field, solvent%partial_of_ion_field )
          !
          DO i = 1, nat
             !
             CALL calc_dboundary_dions( i, solvent, partial )
             !
             ! ... If field-aware interface correct the derivative of the interface function
             !
             IF ( solvent % field_aware ) CALL field_aware_dboundary_dions( i, solvent, partial )
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
!  Subroutine: calc_denviron
!
!> Calculates the Environ contribution to the local potential. All
!! the Environ modules need to be called here. The potentials are
!! all computed on the dense real-space grid and added to vtot.
!--------------------------------------------------------------------
  SUBROUTINE calc_dvenviron( nnr, rho, drho, dvtot )
!--------------------------------------------------------------------
    USE environ_base,  ONLY : vzero, solvent,                       &
                              lelectrostatic, loptical
    !
    ! ... Each contribution to the potential is computed in its module
    !
    USE solvent_tddfpt, ONLY : calc_vsolvent_tddfpt
    !
    IMPLICIT NONE
    !
    ! ... Declares variables
    !
    INTEGER, INTENT(IN)     :: nnr
    REAL( DP ), INTENT(IN)  :: rho( nnr )    !> ground-state charge-density
    REAL( DP ), INTENT(IN)  :: drho( nnr )   !> response charge-density
    REAL( DP ), INTENT(OUT) :: dvtot( nnr )
    !
    ! ... Local variables
    !
    REAL( DP ), DIMENSION( : ), ALLOCATABLE :: dvpol, dvepsilon
    !
    ! ... If any form of electrostatic embedding is present, calculate its contribution
    !
    IF ( lelectrostatic ) THEN
       !
       ! ... Electrostatics is also computed inside the calling program, need to remove the reference
       !
       IF ( loptical ) THEN
          !
          ALLOCATE( dvpol( nnr ) )
          dvpol = 0.D0
          ALLOCATE( dvepsilon( nnr ) )
          dvepsilon = 0.D0
          !
          CALL calc_vsolvent_tddfpt(nnr, 1, rho, drho, dvpol, dvepsilon)
          !
          dvtot = dvtot + dvpol + dvepsilon
          !
          DEALLOCATE( dvpol )
          DEALLOCATE( dvepsilon )
          !
       ENDIF
       !
    END IF
    !
    RETURN
!--------------------------------------------------------------------
  END SUBROUTINE calc_dvenviron
!--------------------------------------------------------------------
!----------------------------------------------------------------------------
END MODULE environ_main
!----------------------------------------------------------------------------
