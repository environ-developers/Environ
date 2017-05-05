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
  ! ... all physical and numerical parameters, and
  ! ... the contributions to the energy and to the potential.
  !
  USE environ_types
  !
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
  REAL (KIND=DP) ::                 &
       e2
  !
  ! Control flags
  !
  LOGICAL ::                          &
       lsolvent
  LOGICAL ::                          &
       lelectrostatic
  LOGICAL ::                          &
       lsoftcavity
  LOGICAL ::                          &
       lrigidcavity
  LOGICAL ::                          &
       lcoredensity
  LOGICAL ::                          &
       lsmearedions
  !
  ! Internal cell parameters
  !
  TYPE( environ_cell ), TARGET ::     &
       cell
  !
  ! Internal parameters of ions
  !
  TYPE( environ_ions ), TARGET ::  &
       ions
  !
  ! Internal parameters of electrons
  !
  TYPE( environ_electrons ), TARGET ::  &
       electrons
  !
  ! Internal parameters of external charges
  !
  LOGICAL ::                            &
       lexternals
  TYPE( environ_externals ), TARGET ::  &
       externals
  !
  ! Internal paramters of total charges
  !
  TYPE( environ_charges ), TARGET :: &
       charges
  !
  ! Details of the selected system
  !
  TYPE( environ_system ), TARGET ::   &
       system
  !
  ! Details of the continuum interface
  !
  TYPE( environ_boundary ), TARGET :: &
       solvent
  !
  ! Dielectric parameters (solvent)
  !
  LOGICAL ::                        &
       lstatic
  LOGICAL ::                        &
       loptical
  REAL (KIND=DP) ::                 &
       env_static_permittivity,     &
       env_optical_permittivity
  TYPE( environ_dielectric ) ::     &
       static,                      &
       optical
  !
  ! Ionic countercharge parameters
  !
  LOGICAL ::                        &
       lelectrolyte
  INTEGER ::                        &
       env_ioncc_ntyp
  TYPE( environ_electrolyte ) ::    &
       electrolyte
  !
  ! Cavitation energy parameters
  !
  LOGICAL ::                        &
       lsurface
  REAL (KIND=DP) ::                 &
       env_surface_tension
  !
  ! PV term parameters
  !
  LOGICAL ::                        &
       lvolume
  REAL (KIND=DP) ::                 &
       env_pressure
  !
  ! Periodicity correction parameters
  !
  LOGICAL ::                        &
       lperiodic
  INTEGER ::                        &
       env_periodicity
  INTEGER ::                        &
       slab_axis
  !
  ! Temporary parameters
  !
  LOGICAL ::                        &
       add_jellium
  INTEGER ::                        &
       nrep
  !
  ! Computed physical variables
  !
  REAL (KIND=DP) ::                 &
       deenviron,                   &
       eelectrostatic,              &
       ecavity,                     &
       epressure
  TYPE ( environ_density ) ::       &
       vzero,                       &
       velectrostatic,              &
       vsoftcavity,                 &
       vcavity,                     &
       vpressure
  TYPE ( environ_gradient ) ::      &
       gelectrostatic
       ! gradient of the total electrostatic potential
       ! (elec + pol + ions) for the TDDFPT calculation
  !
  CONTAINS
     !
     ! ... the following routine copies variables read in input
     ! ... to global variables and derived data types kept in this module
     !
     SUBROUTINE set_environ_base &
                      ( nelec, nspin,                               &
                        nat, ntyp, atom_label, atomicspread,        &
                        corespread, solvationrad,                   &
                        assume_isolated, ibrav,                     &
                        environ_restart_, verbose_, environ_thr_,   &
                        environ_nskip_, environ_type,               &
                        system_ntyp, system_dim, system_axis,       &
                        stype_, rhomax_, rhomin_, tbeta,            &
                        env_static_permittivity_,                   &
                        env_optical_permittivity_, solvent_mode,    &
                        radius_mode, alpha, softness,               &
                        eps_distance, eps_spread,                   &
                        add_jellium_,                               &
                        env_surface_tension_, delta,                &
                        env_pressure_,                              &
                        env_ioncc_ntyp_, nrep_,                     &
                        stern_mode, stern_distance, stern_spread,   &
                        cion, cionmax, rion, zion, rhopb,           &
                        solvent_temperature,                        &
                        env_external_charges,                       &
                        extcharge_charge, extcharge_dim,            &
                        extcharge_axis, extcharge_pos,              &
                        extcharge_spread,                           &
                        env_dielectric_regions,                     &
                        epsregion_eps, epsregion_dim,               &
                        epsregion_axis, epsregion_pos,              &
                        epsregion_spread, epsregion_width )
        !
        USE control_flags,   ONLY : tddfpt
        IMPLICIT NONE
        CHARACTER(LEN=20)   :: sub_name = ' set_environ_base '
        LOGICAL, INTENT(IN) :: environ_restart_, add_jellium_
        INTEGER, INTENT(IN) :: nspin, nelec, nat, ntyp, ibrav,                  &
                               verbose_, environ_nskip_,                        &
                               system_ntyp, system_dim, system_axis,            &
                               stype_, env_ioncc_ntyp_, nrep_,                  &
                               env_external_charges,                            &
                               extcharge_dim(:), extcharge_axis(:),             &
                               env_dielectric_regions,                          &
                               epsregion_dim(:), epsregion_axis(:)
        REAL(DP), INTENT(IN) :: environ_thr_, rhomax_, rhomin_, tbeta,          &
                               env_static_permittivity_,                        &
                               env_optical_permittivity_,                       &
                               alpha, softness, solvationrad(:), corespread(:), &
                               atomicspread(:),                                 &
                               eps_distance, eps_spread,                        &
                               env_surface_tension_, delta, env_pressure_,      &
                               stern_distance, stern_spread, cion(:),           &
                               cionmax(:), rion(:), zion(:), rhopb,             &
                               solvent_temperature,                             &
                               extcharge_charge(:), extcharge_spread(:),        &
                               extcharge_pos(:,:), epsregion_eps(:,:),          &
                               epsregion_pos(:,:), epsregion_spread(:),         &
                               epsregion_width(:)
        CHARACTER( LEN = * ), INTENT(IN) :: assume_isolated, environ_type,      &
                               solvent_mode, radius_mode, stern_mode
        CHARACTER( LEN = 3 ), DIMENSION(:), INTENT(IN) :: atom_label
        INTEGER :: i
        INTEGER :: stype
        REAL(DP) :: rhomax, rhomin
        !
        ! Create necessary local types
        !
        ! CALL create_environ_cell(cell) THIS IS NOT NEEDED AS THERE ARE NO POINTERS OR ALLOCATABLES
        !
        CALL create_environ_electrons(electrons)
        !
        CALL create_environ_ions(ions)
        !
        CALL create_environ_system(system)
        !
        ! General flags
        !
        environ_restart = environ_restart_
        verbose         = verbose_
        environ_thr     = environ_thr_
        environ_nskip   = environ_nskip_
        !
        ! Set main environment flags
        !
        SELECT CASE (TRIM(environ_type))
        ! if a specific environ is selected use hardcoded parameters
        CASE ('vacuum')
           ! vacuum, all flags off
           env_static_permittivity = 1.D0
           env_optical_permittivity = 1.D0
           env_surface_tension = 0.D0
           env_pressure = 0.D0
           env_ioncc_ntyp = 0
           stype = 1
           rhomax = 0.005
           rhomin = 0.0001
        CASE ('water')
           ! water, experimental and SCCS tuned parameters
           env_static_permittivity = 78.3D0
           env_optical_permittivity = 1.776D0
           env_surface_tension = 50.D0*1.D-3*bohr_radius_si**2/rydberg_si
           env_pressure = -0.35D0*1.D9/rydberg_si*bohr_radius_si**3
           env_ioncc_ntyp = 0
           stype = 1
           rhomax = 0.005
           rhomin = 0.0001
        CASE ('water-cation')
           ! water, experimental and SCCS tuned parameters for cations
           env_static_permittivity = 78.3D0
           env_optical_permittivity = 1.776D0
           env_surface_tension = 5.D0*1.D-3*bohr_radius_si**2/rydberg_si
           env_pressure = 0.125D0*1.D9/rydberg_si*bohr_radius_si**3
           env_ioncc_ntyp = 0
           stype = 1
           rhomax = 0.0035
           rhomin = 0.0002
        CASE ('water-anion')
           ! water, experimental and SCCS tuned parameters for anions
           env_static_permittivity = 78.3D0
           env_optical_permittivity = 1.776D0
           env_surface_tension = 0.D0*1.D-3*bohr_radius_si**2/rydberg_si
           env_pressure = 0.450D0*1.D9/rydberg_si*bohr_radius_si**3
           env_ioncc_ntyp = 0
           stype = 1
           rhomax = 0.0155
           rhomin = 0.0024
        CASE ('input')
           ! take values from input, this is the default option
           env_static_permittivity = env_static_permittivity_
           env_optical_permittivity = env_optical_permittivity_
           env_surface_tension = &
             env_surface_tension_*1.D-3*bohr_radius_si**2/rydberg_si
           env_pressure = env_pressure_*1.D9/rydberg_si*bohr_radius_si**3
           env_ioncc_ntyp = env_ioncc_ntyp_
           stype = stype_
           rhomax = rhomax_
           rhomin = rhomin_
        CASE DEFAULT
           call errore (sub_name,'unrecognized value for environ_type',1)
        END SELECT
        !
        ! Set periodic flags according to the host code keyword
        !
        env_periodicity = 3
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
        ! Set basic logical flags
        !
        lstatic        = env_static_permittivity .GT. 1.D0
        loptical       = env_optical_permittivity .GT. 1.D0
        IF ( env_dielectric_regions .GT. 0 ) THEN
           DO i = 1, env_dielectric_regions
              lstatic  = lstatic  .OR. ( epsregion_eps(1,i) .GT. 1.D0 )
              loptical = loptical .OR. ( epsregion_eps(2,i) .GT. 1.D0 )
           ENDDO
        ENDIF
        lsurface       = env_surface_tension .GT. 0.D0
        lvolume        = env_pressure .NE. 0.D0
        lexternals     = env_external_charges .GT. 0
        lelectrolyte   = env_ioncc_ntyp .GT. 0
        lperiodic      = env_periodicity .NE. 3
        !
        ! Derived flags
        !
        lsolvent       = lstatic .OR. loptical .OR. lsurface .OR. lvolume
        lelectrostatic = lstatic .OR. loptical .OR. lelectrolyte .OR. &
                         lexternals .OR. lperiodic
        lsoftcavity    = ( lsolvent .AND. solvent_mode .NE. 'ionic' ) .OR. &
                         ( lelectrolyte .AND. stern_mode .NE. 'ionic' )
        lrigidcavity   = ( lsolvent .AND. solvent_mode .NE. 'electronic' ) .OR. &
                         ( lelectrolyte .AND. stern_mode .NE. 'electronic' )
        lcoredensity   = ( lsolvent .AND. solvent_mode .EQ. 'full' ) .OR. &
                         ( lelectrolyte .AND. stern_mode .EQ. 'full' )
        lsmearedions   = lelectrostatic
        !
        ! Create optional types
        !
        IF ( lexternals ) CALL create_environ_externals(externals)
        !
        IF ( lsolvent ) CALL create_environ_boundary(solvent)
        !
        IF ( lelectrolyte ) CALL create_environ_electrolyte(electrolyte)
        !
        IF ( lstatic ) CALL create_environ_dielectric(static)
        !
        IF ( loptical ) CALL create_environ_dielectric(optical)
        !
        IF ( lelectrostatic ) CALL create_environ_charges(charges)
        !
        ! Allocate and set basic properties of ions
        !
        CALL init_environ_ions_first( nat, ntyp, lsoftcavity, lcoredensity, lsmearedions, &
             & radius_mode, atom_label, atomicspread, corespread, solvationrad, ions )
        !
        ! Set basic properties of electrons
        !
        CALL init_environ_electrons_first( nelec, nspin, electrons )
        !
        ! Set basic properties of the selected system
        !
        CALL init_environ_system( system_ntyp, system_dim, system_axis, ions, system )
        !
        ! Collect free charges if computing electrostatics
        !
        IF ( lelectrostatic ) CALL init_environ_charges_first( ions=ions, electrons=electrons, charges=charges )
        !
        ! Allocate and set basic properties of external charges
        !
        IF ( lexternals ) THEN
           CALL init_environ_externals_first( env_external_charges, extcharge_dim, &
                & extcharge_axis, extcharge_pos, extcharge_spread, extcharge_charge, externals )
           CALL init_environ_charges_first( externals=externals, charges=charges )
        ENDIF
        !
        ! Set the parameters of the solvent boundary
        !
        IF ( lsolvent ) &
             CALL init_environ_boundary_first( solvent_mode, env_static_permittivity, stype, &
             & rhomax, rhomin, tbeta, lsurface, delta, alpha, softness, ions, solvent )
        !
        ! Set the parameters of the electrolyte and of its boundary
        !
        IF ( lelectrolyte ) CALL init_environ_electrolyte_first( env_ioncc_ntyp,      &
             & stern_mode, stype, rhomin, rhopb, tbeta, stern_distance, stern_spread, &
             & alpha, softness, ions, solvent_temperature, cion, cionmax, rion, zion, &
             & electrolyte )
        !
        ! Set the parameters of the dielectric
        !
        IF ( lstatic ) CALL init_environ_dielectric_first( env_static_permittivity,   &
             & env_dielectric_regions, epsregion_dim, epsregion_axis, epsregion_pos,  &
             & epsregion_width, epsregion_spread, epsregion_eps(1,:), solvent, static )
        !
        IF ( loptical ) CALL init_environ_dielectric_first( env_optical_permittivity,  &
             & env_dielectric_regions, epsregion_dim, epsregion_axis, epsregion_pos,   &
             & epsregion_width, epsregion_spread, epsregion_eps(2,:), solvent, optical )
        !
        ! Obsolote keywords to be moved or removed
        !
        add_jellium = add_jellium_
        !
        nrep      = nrep_
        !
        ! The periodic-boundary correction methods
        ! slabx, slaby, slabz, pcc, and esm are
        ! not implemented in TDDFPT.
        !
        IF (tddfpt) THEN
           IF ( trim( assume_isolated ) /= 'makov-payne'       .AND. &
              & trim( assume_isolated ) /= 'martyna-tuckerman' .AND. &
              & trim( assume_isolated ) /= 'none' ) &
              CALL errore (sub_name, &
                       & 'The activated periodic-boundary correction method' // &
                       & ' is not implemented in TDDFPT', 1 )
        ENDIF
        !
     END SUBROUTINE set_environ_base
     !
  !--------------------------------------------------------------------------
END MODULE environ_base
!----------------------------------------------------------------------------
