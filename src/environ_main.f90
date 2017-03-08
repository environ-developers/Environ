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
USE environ_debug,  ONLY: write_cube
!
PRIVATE
!
PUBLIC :: calc_venviron, calc_eenviron, calc_fenviron
!
CONTAINS
!--------------------------------------------------------------------
      SUBROUTINE calc_venviron( update, nnr, nspin, dr2, rho, vltot )
!--------------------------------------------------------------------
      !
      ! ... Calculates the environ contribution to the local
      !     potential. All the Environ modules need to be called here.
      !     The potentials are all computed on the dense real-space
      !     grid and added to vltot.
      !
      USE kinds,         ONLY : DP
      USE environ_ions,  ONLY : rhoions
      USE environ_base,  ONLY : use_smeared_ions, rhopol
      USE environ_base,  ONLY : verbose, eps_mode, env_static_permittivity, &
                                vsolvent, vepsilon, env_surface_tension,    &
                                vcavity, env_pressure, vpressure,           &
                                env_periodicity, vperiodic,                 &
                                env_ioncc_level, vioncc,                    &
                                env_external_charges, vextcharge,           &
                                rhoexternal, env_dielectric_regions
      !
      ! ... Each contribution to the potential is computed in its module
      !
      USE solvent,       ONLY : calc_vsolvent
      USE periodic,      ONLY : calc_vperiodic
      USE cavity,        ONLY : calc_vcavity
      USE pressure,      ONLY : calc_vpressure
      USE ioncc,         ONLY : calc_vioncc
      USE extcharge,     ONLY : calc_vextcharge
      !
      IMPLICIT NONE
      !
      ! ... Declares variables
      !
      LOGICAL, INTENT(IN)       :: update
      INTEGER, INTENT(IN)       :: nnr, nspin
      REAL( DP ), INTENT(IN)    :: dr2
      REAL( DP ), INTENT(IN)    :: rho( nnr, nspin )
      REAL( DP ), INTENT(INOUT) :: vltot( nnr )
      !
      ! ... Local variables
      !
      REAL( DP ), ALLOCATABLE   :: rhoelec(:)
      REAL( DP ), ALLOCATABLE   :: rhotot(:)
      REAL( DP ), ALLOCATABLE   :: vaux(:,:)
      REAL(DP) :: ehart, charge

      ! ... If not updating the potentials, add old potentials and exit

      IF ( .NOT. update ) THEN
        IF ( env_surface_tension .GT. 0.D0 ) vltot = vltot + vcavity
        IF ( env_pressure .NE. 0.D0 ) vltot = vltot + vpressure
        IF ( env_external_charges .GT. 0 ) vltot = vltot + vextcharge
        IF ( env_static_permittivity .GT. 1.D0 .OR. env_dielectric_regions .GT. 0 ) THEN
           vltot = vltot + vsolvent
           IF ( env_static_permittivity .GT. 1.D0 .AND. TRIM(eps_mode) .NE. 'ionic' ) &
           vltot = vltot + vepsilon
        ENDIF
        IF ( env_ioncc_level .GT. 0 ) THEN
           vltot = vltot + vioncc
        ELSE IF ( env_periodicity .NE. 3 ) THEN
           vltot = vltot + vperiodic
        END IF
        RETURN
      END IF

      ! ... Initializes variables, only spinless density is needed

      ALLOCATE( rhoelec( nnr ) )
      rhoelec( : ) = rho( : , 1 )
      IF( nspin==2 ) rhoelec( : ) = rhoelec( : ) + rho( : , 2 )
      IF ( verbose .GE. 4 ) CALL write_cube( nnr, rhoelec, 'rhoelec.cube' )

      ! ... If surface tension greater than zero, calculates cavity contribution

      IF ( env_surface_tension .GT. 0.D0 ) THEN

        vcavity = 0.D0
        CALL calc_vcavity( nnr, rhoelec, vcavity )
        IF ( verbose .GE. 3 ) CALL write_cube( nnr, vcavity, 'vcavity.cube' )
        vltot = vltot + vcavity

      ENDIF

      ! ... If external pressure different from zero, calculates PV contribution

      IF ( env_pressure .NE. 0.D0 ) THEN

        vpressure = 0.D0
        CALL calc_vpressure( nnr, rhoelec, vpressure )
        IF ( verbose .GE. 3 ) CALL write_cube( nnr, vpressure, 'vpressure.cube' )
        vltot = vltot + vpressure

      ENDIF

      ! ... If external charges are present in the system

      IF ( env_external_charges .GT. 0 ) THEN

        CALL calc_vextcharge( nnr, nspin, vextcharge )
        IF ( verbose .GE. 3 ) CALL write_cube( nnr, vextcharge, 'vextcharge.cube')
        vltot = vltot + vextcharge

      END IF

      ! ... If dielectric is not vacuum, calculates solvent contributions

      IF ( env_static_permittivity .GT. 1.D0 .OR. env_dielectric_regions .GT. 0 ) THEN

        vsolvent = 0.D0
        vepsilon = 0.D0
        CALL calc_vsolvent( nnr, nspin, dr2, rhoelec, vsolvent, vepsilon )
        IF ( verbose .GE. 3 ) CALL write_cube( nnr, vsolvent, 'vsolvent.cube' )
        IF ( verbose .GE. 3 ) CALL write_cube( nnr, vepsilon, 'vepsilon.cube' )
        vltot = vltot + vsolvent
        IF ( env_static_permittivity .GT. 1.D0 .AND. TRIM(eps_mode) .NE. 'ionic' ) &
        vltot = vltot + vepsilon

      ENDIF

      ! ... The following effects depend on the total charge distribution of the system

      ALLOCATE( rhotot( nnr ) )
      rhotot = rhoelec
      IF ( use_smeared_ions ) rhotot = rhotot + rhoions
      IF ( env_static_permittivity .GT. 1.D0 .OR. env_dielectric_regions .GT. 0 ) rhotot = rhotot + rhopol
      IF ( env_external_charges .GT. 0 ) rhotot = rhotot + rhoexternal

      ! ... If ionic counter-charge is present, add contribution
      !     note that ionic counter-charge assumes a correction for 3D periodic
      !     boundaries, so this driver is alternative to vperiodic

      IF ( env_ioncc_level .GT. 0 ) THEN

        vioncc = 0.D0
        CALL calc_vioncc( nnr, nspin, rhotot, vioncc )
        IF ( verbose .GE. 3 ) CALL write_cube( nnr, vioncc, 'vioncc.cube' )
        vltot = vltot + vioncc

      ! ... If not in 3D periodic boundary conditions, add corrective potential
      !     note that this term must be computed after the solvent one
      !     since it changes the charge distribution inside the cell

      ELSE IF ( env_periodicity .NE. 3 ) THEN

        vperiodic = 0.D0
        CALL calc_vperiodic( nnr, nspin, .TRUE., rhotot, vperiodic )
        IF ( verbose .GE. 3 ) CALL write_cube( nnr, vperiodic, 'vperiodic.cube' )
        vltot = vltot + vperiodic

      ENDIF

      DEALLOCATE( rhotot )
      DEALLOCATE( rhoelec )

      RETURN
!--------------------------------------------------------------------
      END SUBROUTINE calc_venviron
!--------------------------------------------------------------------
!--------------------------------------------------------------------
      SUBROUTINE calc_eenviron( nnr, nspin, rho, deenviron, esolvent, &
                      ecavity, epressure, eperiodic, eioncc, eextcharge )
!--------------------------------------------------------------------
      !
      ! ... Calculates the environ contribution to the Energy.
      !     We must remove \int v_environ * rhoelec that is
      !     automatically included in the energy computed as sum of
      !     Kohn-Sham eigenvalues.
      !
      USE kinds,         ONLY : DP
      USE mp,            ONLY : mp_sum
      USE mp_bands,      ONLY : intra_bgrp_comm
      USE environ_cell,  ONLY : domega
      USE environ_base,  ONLY : env_static_permittivity, eps_mode,    &
                                vsolvent, vepsilon,                   &
                                env_surface_tension, vcavity,         &
                                env_pressure, vpressure,              &
                                env_periodicity, vperiodic,           &
                                env_ioncc_level, stern_mode,          &
                                vioncc, vgamma,                       &
                                env_external_charges, vextcharge,     &
                                env_dielectric_regions
      !
      ! ... Each contribution to the energy is computed in its module
      !
      USE solvent,       ONLY : calc_esolvent
      USE periodic,      ONLY : calc_eperiodic
      USE cavity,        ONLY : calc_ecavity
      USE pressure,      ONLY : calc_epressure
      USE ioncc,         ONLY : calc_eioncc
      USE extcharge,     ONLY : calc_eextcharge
      !
      IMPLICIT NONE
      !
      ! ... Declares variables
      !
      INTEGER, INTENT(IN)     :: nnr, nspin
      REAL( DP ), INTENT(IN)  :: rho( nnr, nspin )
      REAL( DP ), INTENT(OUT) :: deenviron, esolvent, ecavity, &
                                 epressure, eperiodic, eioncc, &
                                 eextcharge
      !
      ! ... Local variables
      !
      REAL( DP ), ALLOCATABLE :: rhoelec( : )
      !
      ! ... Initializes the variables
      !
      deenviron  = 0.D0
      esolvent   = 0.D0
      ecavity    = 0.D0
      epressure  = 0.D0
      eioncc     = 0.D0
      eextcharge = 0.D0
      !
      ! ... Also in this case only spinless density is needed
      !
      ALLOCATE( rhoelec(nnr) )
      rhoelec( : ) = rho( : , 1)
      IF( nspin==2 ) rhoelec( : ) = rhoelec( : ) + rho( : , 2)
      !
      ! ... Calculates the energy corrections
      !
      !  if surface tension different from zero compute cavitation energy
      !
      IF ( env_surface_tension .GT. 0.D0 ) THEN
         !
         deenviron = deenviron -                                 &
                   SUM( vcavity( : ) * rhoelec( : ) ) * domega
         !
         CALL calc_ecavity( nnr, rhoelec, ecavity )
         !
      END IF
      !
      !  if pressure different from zero compute PV energy
      !
      IF ( env_pressure .NE. 0.D0 ) THEN
         !
         deenviron = deenviron -                                 &
                   SUM( vpressure( : ) * rhoelec( : ) ) * domega
         !
         CALL calc_epressure( nnr, rhoelec, epressure )
         !
      END IF
      !
      !  if external charges are present, compute their contribution
      !
      IF ( env_external_charges .GT. 0) THEN
         !
         deenviron = deenviron -                                 &
                   SUM( vextcharge( : ) * rhoelec( : ) ) * domega
         !
         CALL calc_eextcharge( nnr, rhoelec, eextcharge )
         !
      END IF
      !
      !  if ionic countercharge is present compute extra terms
      !  ioncc possibly includes dielectric and pbc correction terms
      !
      IF ( env_ioncc_level .GT. 0 ) THEN
         !
         deenviron = deenviron -                                 &
              SUM( vioncc( : ) * rhoelec( : ) ) * domega
         !
         IF ( env_ioncc_level .EQ. 3 .AND. TRIM(stern_mode) .NE. 'ionic' ) &
              deenviron = deenviron -                                 &
              SUM( vgamma( : ) * rhoelec( : ) ) * domega
         !
         IF ( env_static_permittivity .GT. 1.D0 .AND. TRIM(eps_mode) .NE. 'ionic' ) &
              deenviron = deenviron -                                 &
              SUM( vepsilon( : ) * rhoelec( : ) ) * domega
         !
         CALL calc_eioncc( nnr, rhoelec, eioncc )
         !
      ELSE
      !
      !  if dielectric different from vacuum compute solvent term
      !
         IF ( env_static_permittivity .GT. 1.D0 .OR. env_dielectric_regions .GT. 0 ) THEN
            !
            deenviron =  deenviron -                                &
                 SUM( vsolvent( : ) * rhoelec( : ) ) * domega
            !
            IF ( env_static_permittivity .GT. 1.D0 .AND. TRIM(eps_mode) .NE. 'ionic' ) &
                 deenviron = deenviron -                                 &
                 SUM( vepsilon( : ) * rhoelec( : ) ) * domega
            !
            CALL calc_esolvent( nnr, nspin, rhoelec, esolvent )
            !
         END IF
      !
      !  if periodic geometry compute boundary condition term
      !  note that this term must be computed after the solvent one
      !
         IF ( env_periodicity .NE. 3 ) THEN
            !
            deenviron =  deenviron -                                &
                 SUM( vperiodic( : ) * rhoelec( : ) ) * domega
            !
            CALL calc_eperiodic(nnr, rhoelec, eperiodic)
            !
         END IF
         !
      END IF
      !
      DEALLOCATE ( rhoelec )
      !
      CALL mp_sum(  deenviron, intra_bgrp_comm )
      !
      RETURN
      !
!--------------------------------------------------------------------
      END SUBROUTINE calc_eenviron
!--------------------------------------------------------------------
!--------------------------------------------------------------------
      SUBROUTINE calc_fenviron( nnr, nspin, nat, force_environ )
!--------------------------------------------------------------------
      !
      ! ... Calculates the environ contribution to the forces.
      !     Due to Hellman-Feynman only a few of Environ modules
      !     have an effect on atomic forces.
      !
      USE kinds,        ONLY : DP
      USE environ_base, ONLY : env_static_permittivity, &
                               env_periodicity, eps_mode,&
                               env_external_charges
      !
      ! ... Each contribution to the forces is computed in its module
      !
      USE solvent,      ONLY : calc_fsolvent
      USE periodic,     ONLY : calc_fperiodic
      USE extcharge,    ONLY : calc_fextcharge
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN) :: nnr, nspin, nat
      REAL( DP ), INTENT(INOUT) :: force_environ( 3, nat )
      !
      ! force_environ already contains the solvent contribution
      ! if the medium is defined on the electronic density
      !
      IF ( env_static_permittivity .GT. 1.D0 ) THEN
        !
        ! if the medium dielectric is defined on the
        ! ionic density, trow away already computed contribution
        ! and compute here the correct forces
        !
        IF ( TRIM(eps_mode) .EQ. 'ionic' ) THEN
          force_environ = 0.D0
          CALL calc_fsolvent( nnr, nat, force_environ )
        END IF
        !
      END IF
      !
      ! if not 3D periodic compute the boundary contribution to forces
      !
      IF ( env_periodicity .NE. 3 ) THEN
        !
        CALL calc_fperiodic( nnr, nat, force_environ )
        !
      END IF
      !
      ! if external charges defined wrt the center of charge,
      ! compute the extra contribution to inter-atomic forces
      !
      IF ( env_external_charges .GT. 0 ) THEN
        !
        CALL calc_fextcharge( nnr, nspin, nat, force_environ )
        !
      END IF
      !
      RETURN
!--------------------------------------------------------------------
      END SUBROUTINE calc_fenviron
!--------------------------------------------------------------------
END MODULE environ_main
