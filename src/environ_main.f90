!
! Copyright (C) 2007-2008 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! This file contains the main drivers to perform the Environ
! calculation. The three subroutines in this file are designed to
! accomodate all the Environ contributions to Energy, Kohn-Sham
! Potential and Forces. Each subroutine is just the collection of
! the calls to the subroutines of each specific contribution.
!
! original version by O. Andreussi, I. Dabo and N. Marzari
!
MODULE environ_main
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
                                vreference,                           &
                                lsoftcavity, vsoftcavity,             &
                                lsurface, lvolume,                    &
                                charges, lstatic, static,             &
                                lelectrolyte, electrolyte
      USE electrostatic_base, ONLY : reference, outer
      !
      ! ... Each contribution to the potential is computed in its module
      !
      USE electrostatic, ONLY : calc_velectrostatic, calc_vreference
      USE cavity,        ONLY : calc_decavity_dboundary
      USE pressure,      ONLY : calc_depressure_dboundary
      USE dielectric,    ONLY : calc_dedielectric_dboundary
      USE generate_boundary, ONLY : solvent_aware_de_dboundary
      !
      IMPLICIT NONE
      !
      ! ... Declares variables
      !
      LOGICAL, INTENT(IN)       :: update
      INTEGER, INTENT(IN)       :: nnr
      REAL( DP ), INTENT(OUT)   :: vtot( nnr )

      ! ... If not updating the potentials, add old potentials and exit

      vtot = vzero % of_r

      IF ( .NOT. update ) THEN
        IF ( lelectrostatic ) vtot = vtot + velectrostatic % of_r - vreference % of_r
        IF ( lsoftcavity ) vtot = vtot + vsoftcavity % of_r
        RETURN
      END IF

      ! ... If any form of electrostatic embedding is present, calculate its contribution

      IF ( lelectrostatic ) THEN

         ! ... Electrostatics is also computed inside the calling program, need to remove the reference

         CALL calc_vreference( reference, charges, vreference )
         IF ( verbose .GE. 3 ) CALL print_environ_density( vreference )

         IF ( lstatic ) THEN
            IF ( lelectrolyte ) THEN
               CALL calc_velectrostatic( setup = outer, charges=charges, dielectric=static,  &
                    & electrolyte=electrolyte, potential=velectrostatic )
            ELSE
               CALL calc_velectrostatic( setup = outer, charges=charges, dielectric=static, &
                    & potential=velectrostatic )
            ENDIF
         ELSE
            IF ( lelectrolyte ) THEN
               CALL calc_velectrostatic( setup = outer, charges=charges, electrolyte=electrolyte, &
                    & potential=velectrostatic )
            ELSE
               CALL calc_velectrostatic( setup = outer, charges=charges, potential=velectrostatic )
            ENDIF
         ENDIF

         IF ( verbose .GE. 3 ) CALL print_environ_density( velectrostatic )
         vtot = vtot + velectrostatic % of_r - vreference % of_r

      END IF

      ! ... Compute the total potential depending on the boundary

      IF ( lsoftcavity ) THEN

         vsoftcavity % of_r = 0.D0

         ! ... If surface tension greater than zero, calculates cavity contribution

         IF ( lsurface ) CALL calc_decavity_dboundary( solvent, vsoftcavity )

         ! ... If external pressure different from zero, calculates PV contribution

         IF ( lvolume ) CALL calc_depressure_dboundary( solvent, vsoftcavity )

         ! ... If dielectric embedding, calcultes dielectric contribution

         IF ( lstatic ) CALL calc_dedielectric_dboundary( static, velectrostatic, vsoftcavity )

         ! ... If solvent-aware interface correct the potential

         IF ( solvent % solvent_aware ) CALL solvent_aware_de_dboundary( solvent, vsoftcavity )

         ! ... Multiply for the derivative of the boundary wrt electronic density

         vsoftcavity % of_r = vsoftcavity % of_r * solvent % dscaled % of_r
         IF ( verbose .GE. 3 ) CALL print_environ_density( vsoftcavity )
         vtot = vtot + vsoftcavity % of_r

      END IF


      RETURN
!--------------------------------------------------------------------
      END SUBROUTINE calc_venviron
!--------------------------------------------------------------------
!--------------------------------------------------------------------
      SUBROUTINE calc_eenviron( deenviron, eelectrostatic, ecavity, &
                              & epressure )
!--------------------------------------------------------------------
      !
      ! ... Calculates the environ contribution to the Energy.
      !     We must remove \int v_environ * rhoelec that is
      !     automatically included in the energy computed as sum of
      !     Kohn-Sham eigenvalues.
      !
      USE environ_base,  ONLY : electrons, solvent,                   &
                                lelectrostatic, velectrostatic,       &
                                vreference,                           &
                                lsoftcavity, vsoftcavity,             &
                                lsurface, lvolume,                    &
                                charges, lstatic, static,             &
                                lelectrolyte, electrolyte
      USE electrostatic_base, ONLY : reference, outer
      !
      ! ... Each contribution to the energy is computed in its module
      !
      USE electrostatic, ONLY : calc_eelectrostatic, calc_ereference
      USE cavity,        ONLY : calc_ecavity
      USE pressure,      ONLY : calc_epressure
      !
      IMPLICIT NONE
      !
      ! ... Declares variables
      !
      REAL( DP ), INTENT(OUT) :: deenviron, eelectrostatic, ecavity, &
                                 epressure
      REAL( DP ) :: ereference
      !
      ! ... Initializes the variables
      !
      deenviron      = 0.D0
      eelectrostatic = 0.D0
      ecavity        = 0.D0
      epressure      = 0.D0
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
         CALL calc_ereference( setup=reference, charges=charges, potential=vreference, energy=ereference )
         !
         deenviron = deenviron - &
              & scalar_product_environ_density(electrons%density,velectrostatic)
         !
         IF ( lstatic ) THEN
            IF ( lelectrolyte ) THEN
               CALL calc_eelectrostatic( setup=outer, charges=charges, dielectric=static, &
                    & electrolyte=electrolyte, potential=velectrostatic, energy=eelectrostatic )
            ELSE
               CALL calc_eelectrostatic( setup=outer, charges=charges, dielectric=static, &
                    & potential=velectrostatic, energy=eelectrostatic )
            ENDIF
         ELSE
            IF ( lelectrolyte ) THEN
               CALL calc_eelectrostatic( setup=outer, charges=charges, electrolyte=electrolyte, &
                    & potential=velectrostatic, energy=eelectrostatic )
            ELSE
               CALL calc_eelectrostatic( setup=outer, charges=charges, potential=velectrostatic, &
                    & energy=eelectrostatic )
            ENDIF
         ENDIF
         !
         eelectrostatic = eelectrostatic - ereference
         !
      END IF
      !
      IF ( lsoftcavity ) deenviron = deenviron - &
              & scalar_product_environ_density(electrons%density,vsoftcavity)
      !
      !  if surface tension different from zero compute cavitation energy
      !
      IF ( lsurface ) CALL calc_ecavity( solvent, ecavity )
      !
      !  if pressure different from zero compute PV energy
      !
      IF ( lvolume ) CALL calc_epressure( solvent, epressure )
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
      USE kinds,        ONLY : DP
      USE environ_base, ONLY : lelectrostatic, velectrostatic,    &
                               charges, lstatic, static,          &
                               lelectrolyte, electrolyte,         &
                               lrigidcavity, lsurface, lvolume,   &
                               lsolvent, solvent
      !
      USE electrostatic_base, ONLY : outer
      !
      ! ... Each contribution to the forces is computed in its module
      !
      USE electrostatic,     ONLY : calc_felectrostatic
      USE cavity,            ONLY : calc_decavity_dboundary
      USE pressure,          ONLY : calc_depressure_dboundary
      USE dielectric,        ONLY : calc_dedielectric_dboundary
      USE generate_boundary, ONLY : calc_dboundary_dions, solvent_aware_de_dboundary
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN) :: nat
      REAL( DP ), INTENT(INOUT) :: force_environ( 3, nat )
      !
      TYPE( environ_cell ), POINTER :: cell
      !
      INTEGER :: i
      TYPE( environ_density ) :: vrigidcavity
      TYPE( environ_gradient ) :: partial
      !
      force_environ = 0.D0
      !
      IF ( lelectrostatic ) THEN
         !
         IF ( lstatic ) THEN
            IF ( lelectrolyte ) THEN
               CALL calc_felectrostatic( setup=outer, natoms=nat, charges=charges, potential=velectrostatic, &
                    & dielectric=static, electrolyte=electrolyte, forces=force_environ )
            ELSE
               CALL calc_felectrostatic( setup=outer, natoms=nat, charges=charges, potential=velectrostatic, &
                    & dielectric=static, forces=force_environ )
            END IF
         ELSE
            IF ( lelectrolyte ) THEN
               CALL calc_felectrostatic( setup=outer, natoms=nat, charges=charges, potential=velectrostatic, &
                    & electrolyte=electrolyte, forces=force_environ )
            ELSE
               CALL calc_felectrostatic( setup=outer, natoms=nat, charges=charges, potential=velectrostatic, &
                    & forces=force_environ )
            END IF
         END IF
         !
      END IF
      !
      ! ... Compute the total forces depending on the boundary
      !
      IF ( lrigidcavity ) THEN
         !
         IF ( lsolvent ) THEN
            cell => solvent % scaled % cell
         ELSE IF ( lelectrolyte ) THEN
            cell => electrolyte % density % cell
         ELSE
            RETURN
         ENDIF
         !
         CALL init_environ_density( cell, vrigidcavity )
         !
         ! ... If surface tension greater than zero, calculates cavity contribution
         !
         IF ( lsurface ) CALL calc_decavity_dboundary( solvent, vrigidcavity )
         !
         ! ... If external pressure different from zero, calculates PV contribution
         !
         IF ( lvolume ) CALL calc_depressure_dboundary( solvent, vrigidcavity )
         !
         ! ... If dielectric embedding, calcultes dielectric contribution
         !
         IF ( lstatic ) CALL calc_dedielectric_dboundary( static, velectrostatic, vrigidcavity )
         !
         ! ... If solvent-aware interface correct the potential
         !
         IF ( solvent % solvent_aware ) CALL solvent_aware_de_dboundary( solvent, vrigidcavity )
         !
         ! ... Multiply for the derivative of the boundary wrt ionic positions
         !
         CALL init_environ_gradient( cell, partial )
         !
         DO i = 1, nat
            !
            CALL calc_dboundary_dions( i, solvent, partial )
            !
            force_environ( :, i ) = force_environ( :, i ) &
                 & - scalar_product_environ_gradient_density( partial, vrigidcavity )
            !
         ENDDO
         !
         CALL destroy_environ_gradient( partial )
         !
         CALL destroy_environ_density( vrigidcavity )
         !
      END IF
      !
      RETURN
      !
!--------------------------------------------------------------------
      END SUBROUTINE calc_fenviron
!--------------------------------------------------------------------
END MODULE environ_main
