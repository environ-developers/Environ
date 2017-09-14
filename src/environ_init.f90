!
! Copyright (C) 2007-2008 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Subroutines to initialize and clean up the global ENVIRON variables
! Only variables that are required outside of the specific modules
! are initiliazed here, module-specific variables should be initialized
! inside each module by the specific initiliazation subruotine
!
! original version by O. Andreussi, I. Dabo and N. Marzari (MIT)
!
MODULE environ_init
  !
  USE environ_types
  USE environ_output
  USE environ_base
  !
  USE ions_utils
  USE boundary
  USE dielectric
  USE electrolyte_utils
  USE externals_utils
  !
  PRIVATE
  !
  PUBLIC :: set_environ_base, environ_initbase, environ_initcell,  &
       & environ_initions, environ_initelectrons, environ_initpotential, &
       & environ_clean
  !
CONTAINS
!--------------------------------------------------------------------
  SUBROUTINE set_environ_base &
!--------------------------------------------------------------------
       ( nelec, nspin,                               &
       & nat, ntyp, atom_label, atomicspread,        &
       & corespread, solvationrad,                   &
       & oldenviron_, environ_restart_, environ_thr_,&
       & environ_nskip_, environ_type,               &
       & system_ntyp, system_dim, system_axis,       &
       & stype_, rhomax_, rhomin_, tbeta,            &
       & env_static_permittivity_,                   &
       & env_optical_permittivity_, solvent_mode,    &
       & radius_mode, alpha, softness,               &
       & eps_distance, eps_spread,                   &
       & add_jellium_,                               &
       & env_surface_tension_,                       &
       & env_pressure_,                              &
       & env_ioncc_ntyp_, nrep_,                     &
       & stern_mode, stern_distance, stern_spread,   &
       & cion, cionmax, rion, zion, rhopb,           &
       & solvent_temperature,                        &
       & env_external_charges,                       &
       & extcharge_charge, extcharge_dim,            &
       & extcharge_axis, extcharge_pos,              &
       & extcharge_spread,                           &
       & env_dielectric_regions,                     &
       & epsregion_eps, epsregion_dim,               &
       & epsregion_axis, epsregion_pos,              &
       & epsregion_spread, epsregion_width )
    !
    ! ... the following routine copies variables read in input
    ! ... to global variables and derived data types kept in this module
    !
    USE electrostatic_base, ONLY : need_pbc_correction, need_gradient, &
         & need_gradlog, need_factsqrt, need_auxiliary
    !
    IMPLICIT NONE
    CHARACTER(LEN=20)   :: sub_name = ' set_environ_base '
    LOGICAL, INTENT(IN) :: oldenviron_, environ_restart_, add_jellium_
    INTEGER, INTENT(IN) :: nspin, nelec, nat, ntyp,       &
         environ_nskip_,                                  &
         system_ntyp, system_dim, system_axis,            &
         stype_, env_ioncc_ntyp_, nrep_,                  &
         env_external_charges,                            &
         extcharge_dim(:), extcharge_axis(:),             &
         env_dielectric_regions,                          &
         epsregion_dim(:), epsregion_axis(:)
    REAL(DP), INTENT(IN) :: environ_thr_, rhomax_,        &
         rhomin_, tbeta,                                  &
         env_static_permittivity_,                        &
         env_optical_permittivity_,                       &
         alpha, softness, solvationrad(:), corespread(:), &
         atomicspread(:),                                 &
         eps_distance, eps_spread,                        &
         env_surface_tension_, env_pressure_,             &
         stern_distance, stern_spread, cion(:),           &
         cionmax(:), rion(:), zion(:), rhopb,             &
         solvent_temperature,                             &
         extcharge_charge(:), extcharge_spread(:),        &
         extcharge_pos(:,:), epsregion_eps(:,:),          &
         epsregion_pos(:,:), epsregion_spread(:),         &
         epsregion_width(:)
    CHARACTER( LEN = * ), INTENT(IN) :: environ_type,     &
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
    oldenviron      = oldenviron_
    environ_restart = environ_restart_
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
       stype = 2
       rhomax = 0.005
       rhomin = 0.0001
    CASE ('water')
       ! water, experimental and SCCS tuned parameters
       env_static_permittivity = 78.3D0
       env_optical_permittivity = 1.D0 ! 1.776D0
       env_surface_tension = 50.D0*1.D-3*bohr_radius_si**2/rydberg_si
       env_pressure = -0.35D0*1.D9/rydberg_si*bohr_radius_si**3
       env_ioncc_ntyp = 0
       stype = 2
       rhomax = 0.005
       rhomin = 0.0001
    CASE ('water-cation')
       ! water, experimental and SCCS tuned parameters for cations
       env_static_permittivity = 78.3D0
       env_optical_permittivity = 1.D0 ! 1.776D0
       env_surface_tension = 5.D0*1.D-3*bohr_radius_si**2/rydberg_si
       env_pressure = 0.125D0*1.D9/rydberg_si*bohr_radius_si**3
       env_ioncc_ntyp = 0
       stype = 2
       rhomax = 0.0035
       rhomin = 0.0002
    CASE ('water-anion')
       ! water, experimental and SCCS tuned parameters for anions
       env_static_permittivity = 78.3D0
       env_optical_permittivity = 1.D0 ! 1.776D0
       env_surface_tension = 0.D0*1.D-3*bohr_radius_si**2/rydberg_si
       env_pressure = 0.450D0*1.D9/rydberg_si*bohr_radius_si**3
       env_ioncc_ntyp = 0
       stype = 2
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
    lperiodic      = need_pbc_correction
    lauxiliary     = need_auxiliary
    !
    ! Derived flags
    !
    ldielectric    = lstatic .OR. loptical
    lsolvent       = ldielectric .OR. lsurface .OR. lvolume
    lelectrostatic = ldielectric .OR. lelectrolyte .OR. &
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
    IF ( lauxiliary ) CALL create_environ_auxiliary(auxiliary)
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
    ! If needed, link auxiliary charges
    !
    IF ( lauxiliary ) CALL init_environ_charges_first( auxiliary=auxiliary, charges=charges )
    !
    ! Set the parameters of the solvent boundary
    !
    IF ( lsolvent ) THEN
       CALL init_environ_boundary_first( ldielectric, need_factsqrt, lsurface, solvent_mode, &
            & stype, rhomax, rhomin, tbeta, env_static_permittivity, alpha, softness, &
            & electrons, ions, solvent )
    ENDIF
    !
    ! Set the parameters of the electrolyte and of its boundary
    !
    IF ( lelectrolyte ) CALL init_environ_electrolyte_first( env_ioncc_ntyp, &
         & stern_mode, stype, rhomin, rhopb, tbeta, env_static_permittivity, &
         & stern_distance, stern_spread, alpha, softness, electrons, ions, &
         & solvent_temperature, cion, cionmax, rion, zion, electrolyte )
    !
    ! Set the parameters of the dielectric
    !
    IF ( lstatic ) THEN
       CALL init_environ_dielectric_first( env_static_permittivity, solvent, &
         & need_gradient, need_factsqrt, need_gradlog, static )
       IF ( env_dielectric_regions .GT. 0 ) CALL set_dielectric_regions( &
         & env_dielectric_regions, epsregion_dim, epsregion_axis, epsregion_pos, &
         & epsregion_width, epsregion_spread, epsregion_eps(1,:), static )
    END IF
    !
    IF ( loptical ) THEN
       CALL init_environ_dielectric_first( env_optical_permittivity, solvent, &
         & need_gradient, need_factsqrt, need_gradlog, optical )
       IF ( env_dielectric_regions .GT. 0 ) CALL set_dielectric_regions( &
         & env_dielectric_regions, epsregion_dim, epsregion_axis, epsregion_pos, &
         & epsregion_width, epsregion_spread, epsregion_eps(2,:), optical )
    END IF
    !
    ! Obsolote keywords to be moved or removed
    !
    add_jellium = add_jellium_
    !
    nrep      = nrep_
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE set_environ_base
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE environ_initbase( n1, n2, n3, ibrav, alat, omega, at, &
                             & nnr, ir_end, comm, me, root, e2 )
!--------------------------------------------------------------------
!
! Subroutine to initialize fundamental quantities needed by the
! environ modules. This subroutine is called by init_run.f90, thus
! only once per pw.x execution.
!
      ! ... Declares modules
      !
      ! In environ_base all the control flags plus the variables
      ! that need to be initialized
      !
      USE environ_base, ONLY : e2_ => e2,                               &
                               cell, electrons, charges,                &
                               vzero, deenviron,                        &
                               lelectrostatic, eelectrostatic,          &
                               velectrostatic, vreference,              &
                               lsoftcavity, vsoftcavity,                &
                               lelectrolyte, electrolyte,               &
                               lsolvent, solvent, lstatic, static,      &
                               loptical, optical,                       &
                               lexternals, externals,                   &
                               lsurface, ecavity, lvolume, epressure
      !
      USE electrostatic_init, ONLY : electrostatic_initbase
      !
      ! Local base initialization subroutines for the different
      ! environ contributions
      !
      IMPLICIT NONE
      !
      ! ... Input variable
      !
      INTEGER, INTENT(IN) :: nnr
      INTEGER, INTENT(IN) :: ir_end
      INTEGER, INTENT(IN) :: n1, n2, n3
      INTEGER, INTENT(IN) :: ibrav
      INTEGER, INTENT(IN) :: comm
      INTEGER, INTENT(IN) :: me
      INTEGER, INTENT(IN) :: root
      REAL(DP), INTENT(IN) :: alat
      REAL(DP), INTENT(IN) :: omega
      REAL(DP), DIMENSION(3,3), INTENT(IN) :: at
      REAL(DP), OPTIONAL, INTENT(IN) :: e2
      !
      CHARACTER( LEN = 80 ) :: label = ' '
      !
      ! ... Common initialization for simulations with Environ
      !
      e2_ = 2.D0
      IF ( PRESENT(e2) ) e2_ = e2
      !
      CALL init_environ_cell( n1, n2, n3, ibrav, alat, omega, at, nnr, ir_end, comm, me, root, cell )
      !
      ! ... Create local storage for base potential, that needs to be modified
      !
      label = 'vzero'
      CALL create_environ_density( vzero, label )
      CALL init_environ_density( cell, vzero )
      !
      deenviron = 0.0_DP
      !
      ! ... Electrostatic contribution
      !
      eelectrostatic  = 0.0_DP
      IF ( lelectrostatic ) THEN
         !
         label = 'velectrostatic'
         CALL create_environ_density( velectrostatic, label )
         CALL init_environ_density( cell, velectrostatic )
         !
         label = 'vreference'
         CALL create_environ_density( vreference, label )
         CALL init_environ_density( cell, vreference )
         !
         CALL electrostatic_initbase( cell )
         !
      END IF
      !
      ! ... Contribution to the potential due to boundary
      !
      IF ( lsoftcavity ) THEN
         !
         label = 'vsoftcavity'
         CALL create_environ_density( vsoftcavity, label )
         CALL init_environ_density( cell, vsoftcavity )
         !
      END IF
      !
      ! ... Cavity contribution
      !
      ecavity  = 0.0_DP
      !
      ! ... Pressure contribution
      !
      epressure  = 0.0_DP
      !
      ! ... Second step of initialization of some environ derived type
      !
      CALL init_environ_electrons_second( cell, electrons )
      !
      IF ( lsolvent ) CALL init_environ_boundary_second( cell, solvent )
      !
      IF ( lstatic ) CALL init_environ_dielectric_second( cell, static )
      !
      IF ( loptical ) CALL init_environ_dielectric_second( cell, optical )
      !
      IF ( lelectrolyte ) CALL init_environ_electrolyte_second( cell, electrolyte )
      !
      IF ( lexternals ) CALL init_environ_externals_second( cell, externals )
      !
      IF ( lauxiliary ) CALL init_environ_auxiliary( cell, auxiliary )
      !
      IF ( lelectrostatic ) CALL init_environ_charges_second( cell, charges )
      !
      RETURN
      !
!--------------------------------------------------------------------
   END SUBROUTINE environ_initbase
!--------------------------------------------------------------------
!--------------------------------------------------------------------
      SUBROUTINE environ_initpotential( nnr, vltot )
!--------------------------------------------------------------------
!
! Save local potential that will be overwritten by environ
!
      ! ... Declares modules
      USE environ_base,  ONLY : vzero
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN) :: nnr
      REAL( DP ), INTENT( IN ) :: vltot(nnr)
      !
      CHARACTER( LEN=80 ) :: sub_name = 'environ_initpotential'
      !
      vzero % update = .TRUE.
      !
      IF ( vzero % cell % nnr .NE. nnr ) &
           & CALL errore(sub_name,'Inconsistent size in input potential',1)
      vzero % of_r = vltot
      !
      RETURN
      !
!--------------------------------------------------------------------
    END SUBROUTINE environ_initpotential
!--------------------------------------------------------------------
!--------------------------------------------------------------------
      SUBROUTINE environ_initcell( omega, at )
!--------------------------------------------------------------------
!
! Initialize the cell-related quantities to be used in the Environ
! modules. This initialization is called by electrons.f90, thus it
! is performed at every step of electronic optimization.
!
      ! ... Declares modules
      USE environ_base,  ONLY : cell, lstatic, loptical, static, optical, &
                                lexternals, externals, lelectrostatic
      ! ... Cell-related updates
      USE dielectric,     ONLY : update_environ_dielectric
      USE electrostatic_init,  ONLY : electrostatic_initcell
!      USE extcharges,    ONLY : update_external_charges
      !
      IMPLICIT NONE
      !
      REAL( DP ), INTENT( IN ) :: omega
      REAL( DP ), INTENT( IN ) :: at(3,3)
      !
      cell%update = .TRUE.
      !
      ! ... Update cell parameters
      !
      CALL update_environ_cell( omega, at, cell )
      !
      ! ... Update fixed quantities defined inside the cell
      !
      IF ( lstatic ) CALL update_environ_dielectric( static )
      IF ( loptical ) CALL update_environ_dielectric( optical )
      !
!      IF ( lexternals ) CALL update_environ_externals( externals )
      !
      IF ( lelectrostatic ) CALL electrostatic_initcell( cell )
      !
      cell%update = .FALSE.
      !
      RETURN
      !
!--------------------------------------------------------------------
   END SUBROUTINE environ_initcell
!--------------------------------------------------------------------
!--------------------------------------------------------------------
      SUBROUTINE environ_initions( nnr, nat, ntyp, ityp, zv, tau )
!--------------------------------------------------------------------
!
! Initialize the ions-related quantities to be used in the Environ
! modules. This initialization is called by electrons.f90, thus it
! is performed at every step of electronic optimization. It may not
! be the most efficient choice, but it is a safe choice.
!
      ! ... Declares modules
      USE environ_base,      ONLY : cell, ions, electrons, system,  &
                                    lsolvent, solvent,              &
                                    lstatic, static,                &
                                    loptical, optical,              &
                                    lelectrolyte, electrolyte,      &
                                    lrigidcavity,                   &
                                    lelectrostatic
      USE boundary,          ONLY : update_environ_boundary,        &
                                    set_soft_spheres
      USE dielectric,        ONLY : update_environ_dielectric
      USE electrostatic_init,     ONLY : electrostatic_initions
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT( IN )     :: nnr, nat, ntyp
      INTEGER, INTENT( IN )     :: ityp( nat )
      REAL ( DP ), INTENT( IN ) :: zv( ntyp )
      REAL ( DP ), INTENT( IN ) :: tau( 3, nat )
      !
      INTEGER :: ia, dim, axis, icor, max_ntyp
      REAL ( DP ) :: charge, spread, dist, pos(3)
      !
      ions%update = .TRUE.
      !
      ! ... Second step of initialization, need to be moved out of here
      !
      CALL init_environ_ions_second( nat, ntyp, ityp, zv, cell, ions )
      IF ( lsolvent ) CALL set_soft_spheres( solvent )
      !
      ! ... Update ions parameters
      !
      CALL update_environ_ions( nat, tau, ions )
      CALL print_environ_ions( ions )
      !
      ! ... Update system parameters
      !
      CALL update_environ_system( system )
      CALL print_environ_system( system )
      !
      ! ... Update rigid environ properties, defined on ions
      !
      IF ( lrigidcavity ) THEN
         !
         IF ( lsolvent ) THEN
            !
            CALL update_environ_boundary( solvent )
            IF ( solvent % update_status .EQ. 2 ) CALL print_environ_boundary( solvent )
            !
            ! ... Update quantities that depend on the solvent boundary
            !
            IF ( lstatic ) THEN
               CALL update_environ_dielectric( static )
               IF ( .NOT. static % update ) CALL print_environ_dielectric( static )
            END IF
            !
            IF ( loptical ) THEN
               CALL update_environ_dielectric( optical )
               IF ( .NOT. optical % update ) CALL print_environ_dielectric( optical )
            END IF
            !
         ENDIF
         !
         IF ( lelectrolyte ) THEN
            CALL update_environ_boundary( electrolyte%boundary )
            IF ( electrolyte % boundary % update_status .EQ. 2 ) &
                 & CALL print_environ_boundary( electrolyte%boundary )
         END IF
         !
      END IF
      !
      IF ( lelectrostatic ) CALL electrostatic_initions( system )
      !
      ions%update = .FALSE.
      !
      RETURN
      !
!--------------------------------------------------------------------
   END SUBROUTINE environ_initions
!--------------------------------------------------------------------
!--------------------------------------------------------------------
   SUBROUTINE environ_initelectrons( nelec, nspin, nnr, rho )
!--------------------------------------------------------------------
!
! Initialize the electrons-related quantities to be used in the Environ
! modules. This initialization is called by electrons.f90, thus it
! is performed at every step of electronic optimization.
!
      ! ... Declares modules
      USE environ_base,      ONLY : electrons, lsolvent, solvent, &
                                    lstatic, static,              &
                                    loptical, optical,            &
                                    lelectrolyte, electrolyte,    &
                                    lsoftcavity
      USE boundary,          ONLY : update_environ_boundary
      USE dielectric,        ONLY : update_environ_dielectric
      !
      IMPLICIT NONE
      !
      REAL( DP ), INTENT( IN )  :: nelec
      INTEGER, INTENT( IN )     :: nspin, nnr
      REAL ( DP ), INTENT( IN ) :: rho( nnr, nspin )
      !
      electrons%update = .TRUE.
      !
      ! ... Update electrons parameters
      !
      CALL update_environ_electrons( nelec, nspin, nnr, rho, electrons )
      CALL print_environ_electrons( electrons )
      !
      ! ... Update soft environ properties, defined on electrons
      !
      IF ( lsoftcavity ) THEN
         !
         IF ( lsolvent ) THEN
            !
            CALL update_environ_boundary( solvent )
            IF ( solvent % update_status .EQ. 2 ) CALL print_environ_boundary( solvent )
            !
            ! ... Update quantities that depend on the solvent boundary
            !
            IF ( lstatic ) THEN
               CALL update_environ_dielectric( static )
               IF ( .NOT. static % update ) CALL print_environ_dielectric( static )
            END IF
            !
            IF ( loptical ) THEN
               CALL update_environ_dielectric( optical )
               IF ( .NOT. optical % update ) CALL print_environ_dielectric( optical )
            END IF
            !
         END IF
         !
         IF ( lelectrolyte ) THEN
            CALL update_environ_boundary( electrolyte%boundary )
            IF ( electrolyte % boundary % update_status .EQ. 2 ) &
                 & CALL print_environ_boundary( electrolyte%boundary )
         END IF
         !
      END IF
      !
      electrons%update = .FALSE.
      !
      RETURN
      !
!--------------------------------------------------------------------
   END SUBROUTINE environ_initelectrons
!--------------------------------------------------------------------
!--------------------------------------------------------------------
      SUBROUTINE environ_clean(lflag)
!--------------------------------------------------------------------
!
! Clean up all the Environ related allocated variables, and call
! clean up subroutines of specific Environ modules. The structure of
! this subroutine should mirror the one of the environ_init...
! subroutines above. This clean up is called by clean_pw.f90
!
      ! ... Declares modules
      !
      ! In environ_base all the control flags plus the variables
      ! that need to be deallocated
      !
      USE environ_base, ONLY : vzero,                                &
                               lelectrostatic, velectrostatic,       &
                               vreference, lsoftcavity, vsoftcavity, &
                               ions, electrons, system, charges,     &
                               lexternals, externals,                &
                               lstatic, static, loptical, optical,   &
                               lsolvent, solvent,                    &
                               lelectrolyte, electrolyte
      !
      USE electrostatic_init, ONLY : electrostatic_clean
      !
      ! Local clean up subroutines for the different contributions
      !
      IMPLICIT NONE
      !
      LOGICAL, INTENT(IN) :: lflag
      LOGICAL :: opnd
      !
      ! ... Deallocate environment variables
      !
      INQUIRE( unit=environ_unit, opened= opnd )
      IF ( opnd ) CLOSE(unit=environ_unit)
      !
      CALL destroy_environ_density( vzero )
      !
      ! ... environ_base variables
      !
      IF ( lelectrostatic ) THEN
         CALL destroy_environ_density( velectrostatic )
         CALL destroy_environ_density( vreference )
         CALL electrostatic_clean( lflag )
      END IF
      IF ( lsoftcavity ) &
           & CALL destroy_environ_density( vsoftcavity )
      !
      ! ... destroy derived types which were allocated in input
      !
      IF ( lelectrostatic ) CALL destroy_environ_charges( lflag, charges )
      IF ( lauxiliary ) CALL destroy_environ_auxiliary( auxiliary )
      IF ( lexternals ) CALL destroy_environ_externals( lflag, externals )
      IF ( lstatic ) CALL destroy_environ_dielectric( lflag, static )
      IF ( loptical ) CALL destroy_environ_dielectric( lflag, optical )
      IF ( lsolvent ) CALL destroy_environ_boundary( lflag, solvent )
      IF ( lelectrolyte ) CALL destroy_environ_electrolyte( lflag, electrolyte )
      !
      CALL destroy_environ_electrons( lflag, electrons )
      CALL destroy_environ_ions( lflag, ions )
      CALL destroy_environ_system( lflag, system )
      !
!!!! DOUBLE CHECK TDDFPT      IF ( .NOT. tddfpt) THEN
!!!! DOUBLE CHECK TDDFPT        IF ( ALLOCATED( rhoions ) ) DEALLOCATE( rhoions )
!!!! DOUBLE CHECK TDDFPT      END IF
      !
      RETURN
      !
!--------------------------------------------------------------------
   END SUBROUTINE environ_clean
!--------------------------------------------------------------------
END MODULE environ_init
