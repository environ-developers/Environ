!
! Copyright (C) 2007-2008 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! The different modules in this file contain all the Environ related variables 
! that need to be passed in and out. All module-specific variables are declared
! inside the appropriate modules. 
!
! original version by O. Andreussi, I. Dabo, and N. Marzari (MIT)
!
!----------------------------------------------------------------------------
MODULE environ_base
  !--------------------------------------------------------------------------
  !
  ! ... this module contains all the main variables needed for the 
  ! ... environ module. This include the control and debug variables,
  ! ... all physical and numerical parameters,
  ! ... the contributions to the energy and to the potential, 
  ! ... and the polarization charge density (to be used to compute forces)
  !
  USE kinds, ONLY :  DP
  SAVE
  !
  ! Global parameters
  !
  INTEGER, PARAMETER ::             &
       environ_unit = 11 
  LOGICAL ::                        &
       environ_restart
  LOGICAL ::                        &
       update_venviron
  INTEGER ::                        &
       verbose
  REAL (KIND=DP) ::                 &
       environ_thr
  INTEGER ::                        &
       environ_nskip
  INTEGER ::                        &
       ir_end
  REAL (KIND=DP) ::                 &
       e2
  !
  ! Switching function parameters
  !
  INTEGER ::                        &
       stype
  REAL (KIND=DP) ::                 &
       rhomax,                      &
       rhomin,                      &
       tbeta
  !
  ! Dielectric parameters (solvent)
  !
  REAL (KIND=DP) ::                 &
       env_static_permittivity,     &
       env_optical_permittivity
  CHARACTER (LEN=256) ::            &
       eps_mode
  REAL (KIND=DP) ::                 &
       alpha
  REAL (KIND=DP), ALLOCATABLE ::    &
       solvationrad(:)
  REAL (KIND=DP), ALLOCATABLE ::    &
       corespread(:)
  !
  ! Smeared ions parameters
  !
  LOGICAL ::                        &
       use_smeared_ions
  REAL (KIND=DP), ALLOCATABLE ::    &
       atomicspread(:)
  !
  ! Makov Payne correction variables
  !
  LOGICAL ::                        &
       add_jellium               
  !
  ! Numerical Differentiators parameters
  !
  INTEGER ::                        &
       ifdtype,                     &
       nfdpoint,                    &
       ncfd
  INTEGER, ALLOCATABLE ::           &
       icfd(:)
  !
  ! Iterative solvation parameters
  !
  CHARACTER (LEN=256) ::            &
       mixtype
  INTEGER ::                        &
       ndiis
  REAL (KIND=DP) ::                 &
       tolrhopol,                   &
       mixrhopol
  !
  ! Cavitation energy parameters
  !
  REAL (KIND=DP) ::                 &
       env_surface_tension,         &
       delta
  !
  ! PV term parameters
  !
  REAL (KIND=DP) ::                 &
       env_pressure
  !
  ! Periodicity correction parameters
  !
  INTEGER ::                        &
       env_periodicity
  INTEGER ::                        &
       slab_axis
  !
  ! Ionic countercharge parameters
  !
  INTEGER ::                        &
       env_ioncc_level,             &
       nrep
  CHARACTER (LEN=256) ::            &
       stern_mode
  REAL (KIND=DP) ::                 &
       stern_distance,              &
       stern_spread,                &
       cion,                        &
       cionmax,                     &
       rion,                        &
       zion,                        &
       rhopb,                       &
       solvent_temperature
  !
  ! External charges parameters
  !
  INTEGER ::                        &
       env_external_charges
  INTEGER, ALLOCATABLE ::           &
       extcharge_dim(:),            &
       extcharge_axis(:)
  REAL (KIND=DP), ALLOCATABLE ::    &
       extcharge_charge(:),         &
       extcharge_spread(:),         &
       extcharge_pos(:,:)      
  !
  ! Dielectric regions parameters
  !
  INTEGER ::                        &
       env_dielectric_regions
  INTEGER, ALLOCATABLE ::           &
       epsregion_dim(:),            &
       epsregion_axis(:)
  REAL (KIND=DP), ALLOCATABLE ::    &
       epsregion_eps(:,:),          &
       epsregion_width(:),          &
       epsregion_spread(:),         &
       epsregion_pos(:,:)      
  !
  ! Computed physical variables
  !
  REAL (KIND=DP) ::                 &
       deenviron,                   &
       esolvent,                    &
       ecavity,                     &
       epressure,                   &
       eperiodic,                   &
       eioncc,                      &
       eextcharge
  REAL (KIND=DP), ALLOCATABLE ::    &
       vltot_zero(:),               &
       vsolvent(:),                 &
       vepsilon(:),                 &
       vcavity(:),                  &
       vpressure(:),                &
       vperiodic(:),                &
       vioncc(:),                   &
       vgamma(:),                   &
       rhopol(:),                   &
       rhoioncc(:),                 &
       rhopolcc(:),                 &
       vextcharge(:),               &
       rhoexternal(:),              &
       epsstatic(:),                &
       epsoptical(:),               &
       gvtot0(:,:)
       ! gradient of the total electrostatic potential
       ! (elec + pol + ions) for the TDDFPT calculation
  !
  CONTAINS
     !
     ! ... the following routine copies input variables read in input
     ! ... to global variables kept in this module
     !
     SUBROUTINE environ_base_init &
                      ( assume_isolated, environ_restart_,          &
                        verbose_, environ_thr_, environ_nskip_,     &
                        environ_type_, stype_, rhomax_, rhomin_,    &
                        tbeta_, env_static_permittivity_,           &
                        env_optical_permittivity_, eps_mode_,       &
                        alpha_, solvationrad_, corespread_,         & 
                        atomicspread_, add_jellium_,                &
                        ifdtype_, nfdpoint_,                        &
                        mixtype_, ndiis_, mixrhopol_, tolrhopol_,   &
                        env_surface_tension_, delta_,               &
                        env_pressure_,                              &
                        env_ioncc_level_, nrep_,                    &
                        stern_mode_, stern_distance_, stern_spread_,&
                        cion_, cionmax_, rion_, zion_, rhopb_,      &
                        solvent_temperature_,                       &
                        env_external_charges_, extcharge_charge_,   & 
                        extcharge_dim_, extcharge_axis_,            &
                        extcharge_pos_, extcharge_spread_,          & 
                        env_dielectric_regions_, epsregion_eps_,    &
                        epsregion_dim_, epsregion_axis_,            &
                        epsregion_pos_, epsregion_spread_,          &
                        epsregion_width_ )       
        !
        USE constants,       ONLY : rydberg_si, bohr_radius_si, amu_si, fpi
        USE control_flags,   ONLY : tddfpt
        IMPLICIT NONE
        LOGICAL, INTENT(IN) :: environ_restart_, add_jellium_
        INTEGER, INTENT(IN) :: verbose_, environ_nskip_, ifdtype_, nfdpoint_,   &
                               ndiis_, stype_, env_ioncc_level_, nrep_,         &
                               env_external_charges_,                           &
                               extcharge_dim_(:), extcharge_axis_(:),           &
                               env_dielectric_regions_,                         &
                               epsregion_dim_(:), epsregion_axis_(:)
        REAL(DP), INTENT(IN) :: environ_thr_, rhomax_, rhomin_, tbeta_,         &
                               env_static_permittivity_,                        &
                               env_optical_permittivity_, mixrhopol_,           &
                               tolrhopol_, alpha_, solvationrad_(:),            &
                               corespread_(:), atomicspread_(:),                &
                               env_surface_tension_, delta_, env_pressure_,     &
                               stern_distance_, stern_spread_, cion_, cionmax_, &
                               rion_, zion_, rhopb_, solvent_temperature_,      &
                               extcharge_charge_(:), extcharge_spread_(:),      &
                               extcharge_pos_(:,:), epsregion_eps_(:,:),        &
                               epsregion_pos_(:,:), epsregion_spread_(:),       &
                               epsregion_width_(:)        
        CHARACTER( LEN = * ), INTENT(IN) :: assume_isolated, environ_type_, &
                                              eps_mode_, stern_mode_, mixtype_
        INTEGER :: i
        !
        environ_restart = environ_restart_
        verbose         = verbose_
        environ_thr     = environ_thr_
        environ_nskip   = environ_nskip_
        !
        stype    = stype_
        rhomax   = rhomax_
        rhomin   = rhomin_
        tbeta    = tbeta_
        IF ( stype .EQ. 1 ) THEN
           tbeta  = LOG( rhomax / rhomin )
        END IF
        !
        eps_mode = eps_mode_
        IF ( ALLOCATED(solvationrad) ) DEALLOCATE( solvationrad )
        ALLOCATE( solvationrad( SIZE(solvationrad_) ) )
        solvationrad( : ) = solvationrad_( : )
        IF ( ALLOCATED(corespread) ) DEALLOCATE( corespread )
        ALLOCATE( corespread( SIZE(corespread_) ) )
        corespread( : ) = corespread_( : )
        IF ( ALLOCATED(atomicspread) ) DEALLOCATE( atomicspread )
        ALLOCATE( atomicspread( SIZE(atomicspread_) ) )
        atomicspread( : ) = atomicspread_( : )
        alpha = alpha_
        add_jellium = add_jellium_
        !
        ifdtype   = ifdtype_
        nfdpoint  = nfdpoint_
        !
        mixtype   = mixtype_
        ndiis     = ndiis_
        mixrhopol = mixrhopol_
        tolrhopol = tolrhopol_
        !
        delta     = delta_
        !
        nrep      = nrep_
        stern_mode = stern_mode_
        stern_distance = stern_distance_
        stern_spread = stern_spread_
        cion      = cion_ * bohr_radius_si**3 / amu_si
        cionmax   = cionmax_ * bohr_radius_si**3 / amu_si
        rion      = rion_
        zion      = zion_
        rhopb     = rhopb_
        solvent_temperature = solvent_temperature_
        IF ( cion .LT. 0.D0 .OR. rion .LT. 0.D0) &
          & call errore ('environ_base_init','cion and rion should be positive',1)
        IF ( cion .GT. 0.D0 .AND. cionmax .LT. cion ) &
          & call errore ('environ_base_init','cionmax should be at least greater than cion',1)
        IF ( cionmax .EQ. 0.D0 .AND. rion .GT. 0.D0 ) &
          & cionmax = 0.64D0 * 3.D0 / fpi / rion**3 
        !
        SELECT CASE (TRIM(environ_type_))
        ! if a specific environ is selected use hardcoded parameters
        CASE ('vacuum')
           ! vacuum, all flags off
           env_static_permittivity = 1.D0
           env_optical_permittivity = 1.D0
           env_surface_tension = 0.D0
           env_pressure = 0.D0
           env_periodicity = 3
           env_ioncc_level = 0
        CASE ('water')
           ! water, experimental and SCCS tuned parameters
           env_static_permittivity = 78.3D0
           env_optical_permittivity = 1.776D0
           env_surface_tension = 50.D0*1.D-3*bohr_radius_si**2/rydberg_si
           env_pressure = -0.35D0*1.D9/rydberg_si*bohr_radius_si**3
           env_periodicity = 3
           env_ioncc_level = 0
           rhomax = 0.005
           rhomin = 0.0001
           tbeta = LOG( rhomax / rhomin )
        CASE ('water-cation')
           ! water, experimental and SCCS tuned parameters for cations
           env_static_permittivity = 78.3D0
           env_optical_permittivity = 1.776D0
           env_surface_tension = 5.D0*1.D-3*bohr_radius_si**2/rydberg_si
           env_pressure = 0.125D0*1.D9/rydberg_si*bohr_radius_si**3
           env_periodicity = 3
           env_ioncc_level = 0
           rhomax = 0.0035
           rhomin = 0.0002
           tbeta = LOG( rhomax / rhomin )
        CASE ('water-anion')
           ! water, experimental and SCCS tuned parameters for anions
           env_static_permittivity = 78.3D0
           env_optical_permittivity = 1.776D0
           env_surface_tension = 0.D0*1.D-3*bohr_radius_si**2/rydberg_si
           env_pressure = 0.450D0*1.D9/rydberg_si*bohr_radius_si**3
           env_periodicity = 3
           env_ioncc_level = 0
           rhomax = 0.0155
           rhomin = 0.0024
           tbeta = LOG( rhomax / rhomin )
        CASE ('input')
           ! take values from input, this is the default option
           env_static_permittivity = env_static_permittivity_
           env_optical_permittivity = env_optical_permittivity_
           env_surface_tension = &
             env_surface_tension_*1.D-3*bohr_radius_si**2/rydberg_si
           env_pressure = env_pressure_*1.D9/rydberg_si*bohr_radius_si**3
           env_periodicity = 3
           env_ioncc_level = env_ioncc_level_
        CASE DEFAULT
           call errore ('environ_base_init','unrecognized value for environ_type',1)
        END SELECT
        !
        env_external_charges = env_external_charges_
        IF ( env_external_charges .GT. 0 ) THEN
          IF ( ALLOCATED( extcharge_dim ) ) DEALLOCATE( extcharge_dim )
          ALLOCATE( extcharge_dim( env_external_charges ) )
          extcharge_dim = extcharge_dim_
          IF ( ALLOCATED( extcharge_axis ) ) DEALLOCATE( extcharge_axis )
          ALLOCATE( extcharge_axis( env_external_charges ) )
          extcharge_axis = extcharge_axis_
          IF ( ALLOCATED( extcharge_spread ) ) DEALLOCATE( extcharge_spread )
          ALLOCATE( extcharge_spread( env_external_charges ) )
          extcharge_spread = extcharge_spread_
          IF ( ALLOCATED( extcharge_charge ) ) DEALLOCATE( extcharge_charge )
          ALLOCATE( extcharge_charge( env_external_charges ) )
          extcharge_charge = extcharge_charge_
          IF ( ALLOCATED( extcharge_pos ) ) DEALLOCATE( extcharge_pos ) 
          ALLOCATE( extcharge_pos( 3, env_external_charges ) )
          DO i = 1, env_external_charges
            extcharge_pos(1:3,i) = extcharge_pos_(1:3,i)
          END DO
        END IF
        !
        env_dielectric_regions = env_dielectric_regions_
        IF ( env_dielectric_regions .GT. 0 ) THEN
          IF ( ALLOCATED( epsregion_dim ) ) DEALLOCATE( epsregion_dim )
          ALLOCATE( epsregion_dim( env_dielectric_regions ) )
          epsregion_dim = epsregion_dim_
          IF ( ALLOCATED( epsregion_axis ) ) DEALLOCATE( epsregion_axis )
          ALLOCATE( epsregion_axis( env_dielectric_regions ) )
          epsregion_axis = epsregion_axis_
          IF ( ALLOCATED( epsregion_spread ) ) DEALLOCATE( epsregion_spread )
          ALLOCATE( epsregion_spread( env_dielectric_regions ) )
          epsregion_spread = epsregion_spread_
          IF ( ALLOCATED( epsregion_width ) ) DEALLOCATE( epsregion_width )
          ALLOCATE( epsregion_width( env_dielectric_regions ) )
          epsregion_width = epsregion_width_
          IF ( ALLOCATED( epsregion_eps ) ) DEALLOCATE( epsregion_eps )
          ALLOCATE( epsregion_eps( 2, env_dielectric_regions ) )
          IF ( ALLOCATED( epsregion_pos ) ) DEALLOCATE( epsregion_pos ) 
          ALLOCATE( epsregion_pos( 3, env_dielectric_regions ) )
          DO i = 1, env_dielectric_regions
            epsregion_eps(1:2,i) = epsregion_eps_(1:2,i)
            epsregion_pos(1:3,i) = epsregion_pos_(1:3,i)
          END DO
        END IF
        ! Need to add a check on periodic corrections 
        !
        slab_axis  = 0
        SELECT CASE( trim( assume_isolated ) )
        !
        CASE( 'slabx' )
           !
           env_periodicity = 2
           slab_axis       = 1
           !
        CASE( 'slaby' )
           !
           env_periodicity = 2
           slab_axis       = 2
           !
        CASE( 'slabz' )
           !
           env_periodicity = 2
           slab_axis       = 3
           !
        CASE( 'pcc' )
           !
           env_periodicity = 0
           !
        END SELECT
        !
        ! The periodic-boundary correction methods 
        ! slabx, slaby, slabz, pcc, and esm are
        ! not implemented in TDDFPT.
        !
        IF (tddfpt) THEN
           IF ( trim( assume_isolated ) /= 'makov-payne'       .AND. &
              & trim( assume_isolated ) /= 'martyna-tuckerman' .AND. &
              & trim( assume_isolated ) /= 'none' ) &
              CALL errore ('environ_base_init', &
                       & 'The activated periodic-boundary correction method' // &
                       & ' is not implemented in TDDFPT', 1 ) 
        ENDIF 
        !
     END SUBROUTINE environ_base_init
     !
  !--------------------------------------------------------------------------
END MODULE environ_base
!----------------------------------------------------------------------------
MODULE environ_cell
  !--------------------------------------------------------------------------
  !
  ! ... this module contains the cell parameters needed for the 
  ! ... environ module
  !
  USE kinds, ONLY :  DP
  SAVE
  !
  ! number of dense real space grid, volume of single grid point, 
  ! volume of cell, lattice spacing, cell vectors in unit of alat
  !
  INTEGER ::                        &
       ntot,                        &
       ibrav
  REAL (KIND=DP) ::                 &
       domega,                      &
       omega,                       &
       alat
  REAL (KIND=DP) ::                 &
       at(3,3)
  !
  !--------------------------------------------------------------------------
END MODULE environ_cell
!----------------------------------------------------------------------------
MODULE environ_ions
  !--------------------------------------------------------------------------
  !
  ! ... this module contains the ions parameters needed for the 
  ! ... environ module
  !
  USE kinds, ONLY :  DP
  SAVE
  !
  ! number of atoms and number of different atomic types 
  !
  INTEGER ::                        &
       nat,                         &
       ntyp
  !
  ! total ionic charge and center of charge
  !
  REAL (KIND=DP) ::                 &
       zvtot
  REAL (KIND=DP) ::                 &
       avg_pos(3)
  !
  ! atom number - atom type label 
  !
  INTEGER, ALLOCATABLE ::           &
       ityp(:)
  !
  ! charge of atomic type, position of atomic number, smeared ions density
  !
  REAL (KIND=DP), ALLOCATABLE ::    &
       zv(:),                       &
       tau(:,:),                    &
       rhoions(:)
  !
  !--------------------------------------------------------------------------
END MODULE environ_ions
!----------------------------------------------------------------------------
