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
                                lsurface, vcavity,                    &
                                lvolume, vpressure,                   &
                                charges, lstatic, static,             &
                                lelectrolyte, electrolyte
      !
      ! ... Each contribution to the potential is computed in its module
      !
      USE cavity,        ONLY : calc_vcavity
      USE pressure,      ONLY : calc_vpressure
      USE electrostatic, ONLY : calc_velectrostatic, calc_vreference
      USE softcavity,    ONLY : calc_vsoftcavity
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
        IF ( lsurface ) vtot = vtot + vcavity % of_r
        IF ( lvolume ) vtot = vtot + vpressure % of_r
        IF ( lelectrostatic ) THEN
           vtot = vtot + velectrostatic % of_r - vreference % of_r
           IF ( lsoftcavity ) vtot = vtot + vsoftcavity % of_r
        ENDIF
        RETURN
      END IF

      ! ... If surface tension greater than zero, calculates cavity contribution

      IF ( lsurface ) THEN

        vcavity % of_r = 0.D0
        CALL calc_vcavity( solvent, vcavity )
        IF ( verbose .GE. 3 ) CALL print_environ_density( vcavity )
        vtot = vtot + vcavity % of_r

      ENDIF

      ! ... If external pressure different from zero, calculates PV contribution

      IF ( lvolume ) THEN

        vpressure % of_r = 0.D0
        CALL calc_vpressure( solvent, vpressure )
        IF ( verbose .GE. 3 ) CALL print_environ_density( vpressure )
        vtot = vtot + vpressure % of_r

      ENDIF

      ! ... If any form of electrostatic embedding is present, calculate its contribution

      IF ( lelectrostatic ) THEN

         ! ... Electrostatics is also computed inside the calling program, need to remove the reference

         CALL calc_vreference( charges, vreference )
         IF ( verbose .GE. 3 ) CALL print_environ_density( vreference )

         IF ( lstatic ) THEN

            IF ( lelectrolyte ) THEN

               CALL calc_velectrostatic( charges=charges, dielectric=static,  &
                    & electrolyte=electrolyte, potential=velectrostatic )
!               IF ( lsoftcavity ) CALL calc_vsoftcavity( charges=charges, dielectric=static, &
!                    & electrolyte=electrolyte, velectrostatic=velectrostatic, vsoftcavity=vsoftcavity )

            ELSE

               CALL calc_velectrostatic( charges=charges, dielectric=static, &
                    & potential=velectrostatic )
               IF ( lsoftcavity ) CALL calc_vsoftcavity( charges=charges, dielectric=static, &
                    & velectrostatic=velectrostatic, vsoftcavity=vsoftcavity )

            ENDIF

         ELSE

            IF ( lelectrolyte ) THEN

               CALL calc_velectrostatic( charges=charges, electrolyte=electrolyte, &
                    & potential=velectrostatic )
!               IF ( lsoftcavity ) CALL calc_vsoftcavity( charges=charges, electrolyte=electrolyte, &
!                    & velectrostatic=velectrostatic, vsoftcavity=vsoftcavity )

            ELSE

               CALL calc_velectrostatic( charges=charges, potential=velectrostatic )

            ENDIF

         ENDIF

         IF ( verbose .GE. 3 ) CALL print_environ_density( velectrostatic )
         vtot = vtot + velectrostatic % of_r - vreference % of_r

         IF ( lsoftcavity ) THEN
            IF ( verbose .GE. 3 ) CALL print_environ_density( vsoftcavity )
            vtot = vtot + vsoftcavity % of_r
         END IF

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
                                lsoftcavity, vsoftcavity,             &
                                lsurface, vcavity,                    &
                                lvolume, vpressure,                   &
                                charges, lstatic, static,             &
                                lelectrolyte, electrolyte
      !
      ! ... Each contribution to the energy is computed in its module
      !
      USE electrostatic, ONLY : calc_eelectrostatic
      USE cavity,        ONLY : calc_ecavity
      USE pressure,      ONLY : calc_epressure
      !
      IMPLICIT NONE
      !
      ! ... Declares variables
      !
      REAL( DP ), INTENT(OUT) :: deenviron, eelectrostatic, ecavity, &
                                 epressure
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
      !  if surface tension different from zero compute cavitation energy
      !
      IF ( lsurface ) THEN
         !
         deenviron = deenviron + &
              & scalar_product_environ_density(electrons%density,vcavity)
         !
         CALL calc_ecavity( solvent, ecavity )
         !
      END IF
      !
      !  if pressure different from zero compute PV energy
      !
      IF ( lvolume ) THEN
         !
         deenviron = deenviron + &
              & scalar_product_environ_density(electrons%density,vpressure)
         !
         CALL calc_epressure( solvent, epressure )
         !
      END IF
      !
      ! if electrostatic is on compute electrostatic energy
      !
      IF ( lelectrostatic ) THEN
         !
         deenviron = deenviron + &
              & scalar_product_environ_density(electrons%density,velectrostatic) - &
              & scalar_product_environ_density(electrons%density,vreference)
         !
         IF ( lsoftcavity ) deenviron = deenviron + &
              & scalar_product_environ_density(electrons%density,vsoftcavity)
         !
         IF ( lstatic ) THEN
            !
            IF ( lelectrolyte ) THEN
               !
               CALL calc_eelectrostatic( charges=charges, dielectric=static, &
                    & electrolyte=electrolyte, energy=eelectrostatic )
               !
            ELSE
               !
               CALL calc_eelectrostatic( charges=charges, dielectric=static, &
                    & energy=eelectrostatic )
               !
            ENDIF
            !
         ELSE
            !
            IF ( lelectrolyte ) THEN
               !
               CALL calc_eelectrostatic( charges=charges, electrolyte=electrolyte, &
                    & energy=eelectrostatic )
               !
            ELSE
               !
               CALL calc_eelectrostatic( charges=charges, energy=eelectrostatic )
               !
            ENDIF
            !
         ENDIF
         !
      END IF
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
      USE environ_base, ONLY : lelectrostatic, charges,           &
                               lstatic, static,                   &
                               lelectrolyte, electrolyte,         &
                               lrigidcavity, lsurface, lvolume,   &
                               solvent
      !
      ! ... Each contribution to the forces is computed in its module
      !
      USE electrostatic, ONLY : calc_felectrostatic
      USE cavity,        ONLY : calc_fcavity
      USE pressure,      ONLY : calc_fpressure
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN) :: nat
      REAL( DP ), INTENT(INOUT) :: force_environ( 3, nat )
      !
      force_environ = 0.D0
      !
      IF ( lelectrostatic ) THEN
         !
         IF ( lstatic ) THEN
            !
            IF ( lelectrolyte ) THEN
               !
               CALL calc_felectrostatic( natoms=nat, charges=charges, dielectric=static, &
                    & electrolyte=electrolyte, forces=force_environ )
               !
            ELSE
               !
               CALL calc_felectrostatic( natoms=nat, charges=charges, dielectric=static, &
                    & forces=force_environ )
               !
            END IF
            !
         ELSE
            !
            IF ( lelectrolyte ) THEN
               !
               CALL calc_felectrostatic( natoms=nat, charges=charges, electrolyte=electrolyte, &
                    & forces=force_environ )
               !
            ELSE
               !
               CALL calc_felectrostatic( natoms=nat, charges=charges, forces=force_environ )
               !
            END IF
            !
         END IF
         !
      END IF
      !
      IF ( lrigidcavity ) THEN
         !
         IF ( lsurface ) THEN
            !
            CALL calc_fcavity( nat, solvent, force_environ )
            !
         END IF
         !
         IF ( lvolume ) THEN
            !
            CALL calc_fpressure( nat, solvent, force_environ )
            !
         END IF
         !
      END IF
      !
      RETURN
!--------------------------------------------------------------------
      END SUBROUTINE calc_fenviron
!--------------------------------------------------------------------
END MODULE environ_main
