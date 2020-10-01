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
! Module to initilize environ-related variables
!
! Authors: Oliviero Andreussi (Department of Physics, UNT)
!          Francesco Nattino  (THEOS and NCCR-MARVEL, EPFL)
!          Ismaila Dabo       (DMSE, Penn State)
!          Nicola Marzari     (THEOS and NCCR-MARVEL, EPFL)
!
!----------------------------------------------------------------------------
MODULE environ_init
!----------------------------------------------------------------------------
  !
  USE modules_constants, ONLY : bohr_radius_si, rydberg_si
  USE core_types
  USE environ_types
  USE environ_output
  USE environ_base
  !
  USE utils_ions
  USE utils_boundary
  USE utils_dielectric
  USE utils_electrolyte
  USE utils_externals
  USE utils_charges
  USE utils_semiconductor

  !
  PRIVATE
  !
  PUBLIC :: set_environ_base, environ_initbase, environ_initcell,  &
       & environ_initions, environ_initelectrons, environ_initpotential, &
       & environ_clean, environ_clean_pw, environ_clean_tddfpt
  !
CONTAINS
  !
!--------------------------------------------------------------------
  SUBROUTINE set_environ_base &
!--------------------------------------------------------------------
! BACKWARD COMPATIBILITY
! Compatible with QE-6.0 QE-6.1.X QE-6.2.X QE-6.3
!       ( prog, nelec, nspin,                         &
! Compatible with QE-6.4.X QE-GIT
       ( prog, nelec,                                &
! END BACKWARD COMPATIBILITY
       & nat, ntyp, atom_label, atomicspread,        &
       & corespread, solvationrad,                   &
       & oldenviron_, environ_restart_, environ_thr_,&
       & environ_nskip_, environ_type,               &
       & system_ntyp, system_dim, system_axis,       &
       & stype, rhomax, rhomin, tbeta,               &
       & env_static_permittivity_,                   &
       & env_optical_permittivity_, solvent_mode,    &
       & derivatives_,                               &
       & radius_mode, alpha, softness,               &
       & solvent_distance, solvent_spread,           &
       & solvent_radius, radial_scale,               &
       & radial_spread, filling_threshold,           &
       & filling_spread,                             &
       & field_awareness, charge_asymmetry,          &
       & field_max, field_min,                       &
       & add_jellium_,                               &
       & env_surface_tension_,                       &
       & env_pressure_,                              &
       & env_confine_,                               &
       & env_electrolyte_ntyp_,                      &
       & electrolyte_linearized,                     &
       & electrolyte_entropy, electrolyte_mode,      &
       & electrolyte_distance,                       &
       & electrolyte_spread,                         &
       & cion, cionmax, rion, zion,                  &
       & electrolyte_rhomax, electrolyte_rhomin,     &
       & electrolyte_tbeta,                          &
       & electrolyte_alpha, electrolyte_softness,    &
       & ion_adsorption, ion_adsorption_energy,      &
       & temperature,                                &
       & sc_permittivity, sc_carrier_density, sc_electrode_chg,       &
       & sc_distance, sc_spread, sc_chg_thr,                     &
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
    ! ... to global variables kept in the environ_base module
    !
    USE electrostatic_base, ONLY : need_pbc_correction, need_gradient, &
         & need_factsqrt, need_auxiliary, need_electrolyte, &
         & need_semiconductor, need_outer_loop
    !
    USE core_base
    !
    IMPLICIT NONE
    CHARACTER(LEN=20)   :: sub_name = ' set_environ_base '
    LOGICAL, INTENT(IN) :: oldenviron_, environ_restart_, &
         add_jellium_, electrolyte_linearized
! BACKWARD COMPATIBILITY
! Compatible with QE-6.0 QE-6.1.X QE-6.2.X QE-6.3.X
!    INTEGER, INTENT(IN) :: nspin
! Compatible with QE-6.4.X QE-GIT
!
! END BACKWARD COMPATIBILITY
    INTEGER, INTENT(IN) :: nelec, nat, ntyp,       &
         environ_nskip_,                                  &
         system_ntyp, system_dim, system_axis,            &
         stype, env_electrolyte_ntyp_,                    &
         env_external_charges,                            &
         extcharge_dim(:), extcharge_axis(:),             &
         env_dielectric_regions,                          &
         epsregion_dim(:), epsregion_axis(:)
    REAL(DP), INTENT(IN) :: environ_thr_, rhomax,         &
         rhomin, tbeta,                                   &
         env_static_permittivity_,                        &
         env_optical_permittivity_,                       &
         alpha, softness,                                 &
         solvent_radius, radial_scale, radial_spread,     &
         filling_threshold, filling_spread,               &
         field_awareness, charge_asymmetry, field_max,    &
         field_min,                                       &
         solvationrad(:), corespread(:), atomicspread(:), &
         solvent_distance, solvent_spread,                &
         env_surface_tension_, env_pressure_,             &
         env_confine_,                                    &
         electrolyte_distance, electrolyte_spread,        &
         cion(:), cionmax, rion, zion(:),                 &
         electrolyte_rhomax, electrolyte_rhomin,          &
         electrolyte_tbeta, electrolyte_alpha,            &
         electrolyte_softness, temperature,               &
         ion_adsorption_energy,                           &
         sc_permittivity, sc_carrier_density, sc_electrode_chg,    &
         sc_distance, sc_spread, sc_chg_thr,                          &
         extcharge_charge(:), extcharge_spread(:),        &
         extcharge_pos(:,:), epsregion_eps(:,:),          &
         epsregion_pos(:,:), epsregion_spread(:),         &
         epsregion_width(:)
    CHARACTER( LEN = * ), INTENT(IN) :: prog, environ_type, &
         solvent_mode, radius_mode, electrolyte_mode,     &
         electrolyte_entropy, ion_adsorption, derivatives_
    CHARACTER( LEN = 3 ), DIMENSION(:), INTENT(IN) :: atom_label
    INTEGER :: i
    INTEGER :: stype_
    REAL(DP) :: rhomax_, rhomin_ !, at_max, at_min, tot_len
    CHARACTER( LEN = 80 ) :: label
    !
! BACKWARD COMPATIBILITY
! Compatible with QE-5.X QE-6.1.X QE-6.2.X
!    ltddfpt = tddfpt
! Compatible with QE-6.3.X and QE-GIT
    SELECT CASE ( prog(1:2) )
    CASE ( 'TD' )
       ltddfpt = .TRUE.
    CASE DEFAULT
       ltddfpt = .FALSE.
    END SELECT
! END BACKWARD COMPATIBILITY
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
    ! Set main environment flags, convert to internal units
    !
    env_static_permittivity = env_static_permittivity_
    env_optical_permittivity = env_optical_permittivity_
    env_surface_tension = &
         env_surface_tension_*1.D-3*bohr_radius_si**2/rydberg_si
    env_pressure = env_pressure_*1.D9/rydberg_si*bohr_radius_si**3
    env_confine = env_confine_
    env_electrolyte_ntyp = env_electrolyte_ntyp_
    !stype = stype_
    !rhomax = rhomax_
    !rhomin = rhomin_

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
    lconfine       = env_confine .NE. 0.D0
    lexternals     = env_external_charges .GT. 0
    lelectrolyte   = env_electrolyte_ntyp .GT. 0 .OR. need_electrolyte
    lsemiconductor = need_semiconductor
    lperiodic      = need_pbc_correction
    louterloop     = need_outer_loop
    !
    ! Derived flags
    !
    ldielectric    = lstatic .OR. loptical!.OR. lsemiconductor
    lsolvent       = ldielectric .OR. lsurface .OR. lvolume .OR. lconfine
    lelectrostatic = ldielectric .OR. lelectrolyte .OR. &
                     lexternals .OR. lperiodic
    lsoftsolvent   = lsolvent .AND. ( solvent_mode .EQ. 'electronic' .OR. &
        & solvent_mode .EQ. 'full' .OR. solvent_mode(1:2) .EQ. 'fa' )
    lsoftelectrolyte = lelectrolyte .AND. ( electrolyte_mode .EQ. 'electronic' & 
        & .OR. electrolyte_mode .EQ. 'full' .OR. electrolyte_mode(1:2) .EQ. 'fa' )
    lsoftcavity    = lsoftsolvent .OR. lsoftelectrolyte
    lrigidsolvent  = lsolvent .AND. solvent_mode .NE. 'electronic'
    lrigidelectrolyte = lelectrolyte .AND. electrolyte_mode .NE. 'electronic'
    lrigidcavity   = lrigidsolvent .OR. lrigidelectrolyte
    lcoredensity   = ( lsolvent .AND. solvent_mode .EQ. 'full' ) .OR. &
                     ( lelectrolyte .AND. electrolyte_mode .EQ. 'full' )
    lsmearedions   = lelectrostatic
    lboundary      = lsolvent .OR. lelectrolyte
    lgradient      = ldielectric .OR. ( solvent_mode(1:2) .EQ. 'fa' )
    !
    ! Create optional types
    !
    IF ( lexternals ) CALL create_environ_externals(externals)
    !
    IF ( lboundary ) CALL create_boundary_core(derivatives)
    !
    IF ( lsolvent ) THEN
       label = 'solvent'
       CALL create_environ_boundary(solvent, label)
    ENDIF
    !
    IF ( lelectrolyte ) CALL create_environ_electrolyte(electrolyte)
    !
    IF ( lsemiconductor ) CALL create_environ_semiconductor(semiconductor)
    !
    IF ( lstatic ) CALL create_environ_dielectric(static)
    !
    IF ( loptical ) CALL create_environ_dielectric(optical)
    !
    IF ( lelectrostatic .OR. lconfine ) CALL create_environ_charges(charges)
    !
    ! Allocate and set basic properties of ions
    !
    CALL init_environ_ions_first( nat, ntyp, lsoftcavity, lcoredensity, lsmearedions, &
         & radius_mode, atom_label, atomicspread, corespread, solvationrad, ions )
    !
    ! Set basic properties of electrons
    !
! BACKWARD COMPATIBILITY
! Compatible with QE-6.0 QE-6.1.X QE-6.2.X QE-6.3.X
!    CALL init_environ_electrons_first( nelec, nspin, electrons )
! Compatible with QE-6.4.X QE-GIT
    CALL init_environ_electrons_first( nelec, electrons )
! END BACKWARD COMPATIBILITY
    !
    ! Set basic properties of the selected system
    !
    CALL init_environ_system( system_ntyp, system_dim, system_axis, ions, system )
    !
    ! Collect free charges if computing electrostatics or confinement
    !
    IF ( lelectrostatic .OR. lconfine ) CALL init_environ_charges_first( electrons=electrons, charges=charges )
    IF ( lelectrostatic ) CALL init_environ_charges_first( ions=ions, charges=charges )
    !
    ! Allocate and set basic properties of external charges
    !
    IF ( lexternals ) THEN
       CALL init_environ_externals_first( env_external_charges, extcharge_dim, &
            & extcharge_axis, extcharge_pos, extcharge_spread, extcharge_charge, externals )
       CALL init_environ_charges_first( externals=externals, charges=charges )
    ENDIF
    !
    ! Setup cores needed for derivatives of boundaries
    !
    IF ( lboundary ) THEN
       lfft = .TRUE.
       IF ( derivatives_ .EQ. 'fd' ) lfd = .TRUE.
       CALL init_boundary_core( derivatives_, derivatives, fft, fd )
    ENDIF
    !
    ! Set the parameters of the solvent boundary
    !
    IF ( lsolvent ) THEN
       CALL init_environ_boundary_first( lgradient, need_factsqrt, lsurface, solvent_mode, &
            & stype, rhomax, rhomin, tbeta, env_static_permittivity, alpha, softness, &
            & solvent_distance, solvent_spread, solvent_radius, radial_scale, radial_spread, &
            & filling_threshold, filling_spread, field_awareness, charge_asymmetry, field_max, &
            & field_min, electrons, ions, system, derivatives, solvent )
    ENDIF
    !
    ! Set the parameters of the electrolyte and of its boundary
    !
    WRITE(environ_unit, *) "electrolyte_distance: ",electrolyte_distance
    IF ( lelectrolyte ) THEN
       CALL init_environ_electrolyte_first( env_electrolyte_ntyp, &
            & electrolyte_mode, stype, electrolyte_rhomax, electrolyte_rhomin, &
            & electrolyte_tbeta, env_static_permittivity, &
            & electrolyte_alpha, electrolyte_softness, electrolyte_distance, &
            & electrolyte_spread, solvent_radius, &
            & radial_scale, radial_spread, filling_threshold, filling_spread, &
            & field_awareness, charge_asymmetry, field_max, field_min, &
            & electrons, ions, system, derivatives, temperature, cion, cionmax, rion, &
            & zion, electrolyte_entropy, ion_adsorption, ion_adsorption_energy, &
            & electrolyte_linearized, electrolyte )
       CALL init_environ_charges_first( electrolyte=electrolyte, charges=charges )
    END IF
    !
    ! Set the parameters of the semiconductor
    !
    WRITE(environ_unit, *)"sc_distance: ",sc_distance
    IF ( lsemiconductor ) THEN
       CALL init_environ_semiconductor_first( temperature, sc_permittivity, &
           & sc_carrier_density , sc_electrode_chg, sc_distance, sc_spread, &
           & sc_chg_thr, system,semiconductor)
       CALL init_environ_charges_first( semiconductor=semiconductor, charges = charges )


       ! Adding a dielectric region with the permittivity of the semiconductor
       ! to the side of the slab that wil be experiencing a mott schottky correction
       ! Appears to be much more trouble than it's worth at the moment
       ! (making new arrays of n+1 in fortran is unnecisarily painful). We'll come back to this
       !
       ! env_dielectric_regions = env_dielectric_regions + 1
       ! epsregion_eps(env_dielectric_regions) = semiconductor%permittivity
       ! epsregion_dim(env_dielectric_regions) = 2
       ! epsregion_axis(env_dielectric_regions) = semiconductor%simple%axis
       !
       ! at_max = -0.4
       ! at_min = 10000.D0
       ! DO i = 1, nat
       !    IF (system%ions%tau(i,semiconductor%simple%axis) > at_max ) &
       !         & at_max = system%ions%tau(i,semiconductor%simple%axis)
       !    IF (system%ions%tau(i,semiconductor%simple%axis) < at_min ) &
       !         & at_min = system%ions%tau(i,semiconductor%simple%axis)
       ! END DO
       !
       ! tot_len = cell%alat*cell%at(semiconductor%simple%axis,semiconductor%simple%axis)
       ! epsregion_width(env_dielectric_regions) = ( tot_len- (at_max-at_min))/2.D0
       !
       ! IF ((at_max + epsregion_width(env_dielectric_regions)) >  tot_len) THEN
       !     DO i= 1, 3
       !         IF (i .EQ. semiconductor%simple%axis) THEN
       !             epsregion_pos(i,env_dielectric_regions) = at_min -epsregion_width(env_dielectric_regions)
       !         ELSE
       !             epsregion_pos(i,env_dielectric_regions) = 0.D0
       !         END IF
       !     END DO
       ! ELSE
       !     DO i= 1, 3
       !         IF (i .EQ. semiconductor%simple%axis) THEN
       !             epsregion_pos(iat,env_dielectric_regions) = at_max +epsregion_width(env_dielectric_regions)
       !         ELSE
       !             epsregion_pos(iat,env_dielectric_regions) = 0.D0
       !         END IF
       !     END DO
       ! END IF
       !
       ! epsregion_spread(env_dielectric_regions) = sc_spread
    ENDIF
    !
    ! Set the parameters of the dielectric
    !
    IF ( lstatic ) THEN
       CALL init_environ_dielectric_first( env_static_permittivity, solvent, &
         & need_gradient, need_factsqrt, need_auxiliary, static )
       IF ( env_dielectric_regions .GT. 0 ) CALL set_dielectric_regions( &
         & env_dielectric_regions, epsregion_dim, epsregion_axis, epsregion_pos, &
         & epsregion_width, epsregion_spread, epsregion_eps(1,:), static )
       CALL init_environ_charges_first( dielectric=static, charges=charges )
    END IF
    !
    IF ( loptical ) THEN
       CALL init_environ_dielectric_first( env_optical_permittivity, solvent, &
         & need_gradient, need_factsqrt, need_auxiliary, optical )
       IF ( env_dielectric_regions .GT. 0 ) CALL set_dielectric_regions( &
         & env_dielectric_regions, epsregion_dim, epsregion_axis, epsregion_pos, &
         & epsregion_width, epsregion_spread, epsregion_eps(2,:), optical )
    END IF
    !
    ! Obsolote keywords to be moved or removed
    !
    add_jellium = add_jellium_
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE set_environ_base
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE environ_initbase( alat, at, &
                             & comm, me, root, &
                             & gcutm, e2 )
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
    USE modules_constants, ONLY : e2_ => e2
    USE environ_base, ONLY :                      &
         cell, electrons, charges,                &
         vzero, deenviron,                        &
         lelectrostatic, eelectrostatic,          &
         velectrostatic, vreference,              &
         lsoftcavity, vsoftcavity,                &
         lelectrolyte, electrolyte,               &
         lsolvent, solvent, lstatic, static,      &
         loptical, optical,                       &
         lexternals, externals,                   &
         lsurface, esurface, lvolume, evolume,    &
         lconfine, vconfine, econfine,            &
         eelectrolyte
    !
    USE cell_types, ONLY : init_environ_cell
    USE core_init, ONLY : core_initbase
    USE utils_fft, ONLY : init_dfft_core
    USE core_base
    !
    ! Local base initialization subroutines for the different
    ! environ contributions
    !
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: comm
    INTEGER, INTENT(IN) :: me
    INTEGER, INTENT(IN) :: root
    REAL(DP), INTENT(IN) :: alat
    REAL(DP), DIMENSION(3,3), INTENT(IN) :: at
    REAL( DP ), INTENT(IN) :: gcutm
    REAL(DP), OPTIONAL, INTENT(IN) :: e2
    !
    CHARACTER( LEN = 80 ) :: label = ' '
    !
    ! ... Common initialization for simulations with Environ
    !
    e2_ = 2.D0
    IF ( PRESENT(e2) ) e2_ = e2
    !
    !Create the dfft type here and pass arguments to init_envrion_cell
    CALL init_dfft_core( cell, gcutm, comm, at )
    !
    CALL init_environ_cell( alat, at, &
         & me, root, cell, dfft )
    !
    ! ... Initialization of numerical cores
    !
    CALL core_initbase( cell, gcutm )
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
       WRITE( environ_unit, * )"velectrostatic started"
       CALL create_environ_density( velectrostatic, label )
       CALL init_environ_density( cell, velectrostatic )
       !
       label = 'vreference'

       CALL create_environ_density( vreference, label )
       WRITE( environ_unit, * )"vreference created"
       CALL init_environ_density( cell, vreference )
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
    esurface  = 0.0_DP
    !
    ! ... Pressure contribution
    !
    evolume   = 0.0_DP
    !
    ! ... Confinement contribution
    !
    econfine   = 0.0_DP
    IF ( lconfine ) THEN
       !
       label = 'vconfine'
       CALL create_environ_density( vconfine, label )
       CALL init_environ_density( cell, vconfine )
       !
    END IF
    !
    ! ... Non-electrostatice electrolyte contribution
    !
    eelectrolyte = 0.0_DP
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
    IF ( lsemiconductor ) CALL init_environ_semiconductor_second( cell, semiconductor)
    !
    IF ( lelectrostatic .OR. lconfine ) CALL init_environ_charges_second( cell, charges )
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
    IF ( .NOT. ASSOCIATED( vzero % cell ) ) RETURN
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
  SUBROUTINE environ_initcell( at )
!--------------------------------------------------------------------
!
! Initialize the cell-related quantities to be used in the Environ
! modules. This initialization is called by electrons.f90, thus it
! is performed at every step of electronic optimization.
!
    ! ... Declares modules
    !
    USE environ_base,        ONLY : cell, lstatic, static, &
                                    loptical, optical,  &
                                    lexternals, externals,     &
                                    lelectrolyte, electrolyte, &
                                    lelectrostatic, &
                                    lsemiconductor, semiconductor
    ! ... Cell-related updates
    !
    USE cell_types,          ONLY : update_environ_cell
    USE utils_dielectric,    ONLY : update_environ_dielectric
    USE utils_electrolyte,   ONLY : update_environ_electrolyte
    USE utils_externals,     ONLY : update_environ_externals
    USE utils_semiconductor, ONLY : update_environ_semiconductor
    !
    USE core_init,           ONLY : core_initcell
    !
    IMPLICIT NONE
    !
    REAL( DP ), INTENT( IN ) :: at(3,3)
    !
    cell%update = .TRUE.
    !
    ! ... Update cell parameters
    !
    CALL update_environ_cell( at, cell )
    !
    ! ... Update fixed quantities defined inside the cell
    !
    IF ( lstatic ) CALL update_environ_dielectric( static )
    IF ( loptical ) CALL update_environ_dielectric( optical )
    IF ( lelectrolyte ) CALL update_environ_electrolyte( electrolyte )
    IF ( lexternals ) CALL update_environ_externals( externals )
    IF ( lsemiconductor ) CALL update_environ_semiconductor( semiconductor )
    !
    ! ... Update cores
    !
    CALL core_initcell( cell )
    !
    cell%update = .FALSE.
    !
    RETURN
    !
!--------------------------------------------------------------------
   END SUBROUTINE environ_initcell
!--------------------------------------------------------------------
!--------------------------------------------------------------------
   SUBROUTINE environ_initions( nnr, nat, ntyp, ityp, zv, tau, vloc )
!--------------------------------------------------------------------
!
! Initialize the ions-related quantities to be used in the Environ
! modules. This initialization is called by electrons.f90, thus it
! is performed at every step of electronic optimization. It may not
! be the most efficient choice, but it is a safe choice.
!
     ! ... Declares modules
     USE environ_base,       ONLY : cell, ions, electrons, system,  &
                                    lsolvent, solvent,              &
                                    lstatic, static,                &
                                    loptical, optical,              &
                                    lelectrolyte, electrolyte,      &
                                    lrigidcavity,                   &
                                    lelectrostatic, charges
     USE utils_boundary,     ONLY : update_environ_boundary,        &
                                    set_soft_spheres
     USE utils_dielectric,   ONLY : update_environ_dielectric
     USE utils_electrolyte,  ONLY : update_environ_electrolyte
     USE core_init,          ONLY : core_initions
     !
     IMPLICIT NONE
     !
     INTEGER, INTENT( IN )     :: nnr, nat, ntyp
     INTEGER, INTENT( IN )     :: ityp( nat )
     REAL ( DP ), INTENT( IN ) :: zv( ntyp )
     REAL ( DP ), INTENT( IN ) :: tau( 3, nat )
     REAL ( DP ), INTENT( IN ) :: vloc( nnr, ntyp )
     !
     INTEGER :: ia, dim, axis, icor, max_ntyp
     REAL ( DP ) :: charge, spread, dist, pos(3)
     !
     ions%update = .TRUE.
     !
     ! ... Second step of initialization, need to be moved out of here
     !
     CALL init_environ_ions_second( nat, ntyp, nnr, ityp, zv, cell, vloc, ions )
     IF ( lsolvent ) CALL set_soft_spheres( solvent )
     IF ( lelectrolyte ) CALL set_soft_spheres( electrolyte%boundary )
     !
     ! ... Update ions parameters
     !
     CALL update_environ_ions( nat, tau, ions )
     CALL print_environ_ions( ions )
     !
     ! ... Update system parameters
     !
     system%update = .TRUE.
     CALL update_environ_system( system )
     CALL print_environ_system( system )
     !
     ! ... Update cores
     !
     CALL core_initions( system%pos )
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
           CALL update_environ_electrolyte( electrolyte )
           IF ( .NOT. electrolyte % update ) CALL print_environ_electrolyte( electrolyte )
        END IF
        !
     END IF
     !
     IF ( lelectrostatic .OR. lconfine ) CALL update_environ_charges( charges )
     !
     system%update = .FALSE.
     ions%update = .FALSE.
     !
     RETURN
     !
!--------------------------------------------------------------------
   END SUBROUTINE environ_initions
!--------------------------------------------------------------------
!--------------------------------------------------------------------
! BACKWARD COMPATIBILITY
! Compatible with QE-6.0 QE-6.1.X QE-6.2.X QE-6.3.X
!   SUBROUTINE environ_initelectrons( nspin, nnr, rho, nelec )
! Compatible with QE-6.4.X QE-GIT
   SUBROUTINE environ_initelectrons( nnr, rho, nelec )
! END BACKWARD COMPATIBILITY
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
                                   lsoftcavity, lsoftsolvent,    &
                                   lsoftelectrolyte,             &
                                   lelectrostatic, charges
     USE utils_boundary,    ONLY : update_environ_boundary
     USE utils_dielectric,  ONLY : update_environ_dielectric
     !
     IMPLICIT NONE
     !
     INTEGER, INTENT( IN )     :: nnr
! BACKWARD COMPATIBILITY
! Compatible with QE-6.0 QE-6.1.X QE-6.2.X QE-6.3.X
!     INTEGER, INTENT( IN )     :: nspin
!     REAL ( DP ), INTENT( IN ) :: rho( nnr, nspin )
! Compatible with QE-6.4.X QE-GIT
     REAL ( DP ), INTENT( IN ) :: rho( nnr )
! END BACKWARD COMPATIBILITY
     REAL( DP ), INTENT( IN ), OPTIONAL  :: nelec
     !
     electrons%update = .TRUE.
     !
     ! ... Update electrons parameters
     !
! BACKWARD COMPATIBILITY
! Compatible with QE-6.0 QE-6.1.X QE-6.2.X QE-6.3.X
!     CALL update_environ_electrons( nspin, nnr, rho, electrons, nelec )
! Compatible with QE-6.4.X QE-GIT
     CALL update_environ_electrons( nnr, rho, electrons, nelec )
! END BACKWARD COMPATIBILITY
     CALL print_environ_electrons( electrons )
     !
     IF ( lelectrostatic .OR. lconfine ) CALL update_environ_charges( charges )
     !
     ! ... Update soft environ properties, defined on electrons
     !
     IF ( lsoftcavity ) THEN
        !
        IF ( lsoftsolvent ) THEN
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
        IF ( lsoftelectrolyte ) THEN
           CALL update_environ_boundary( electrolyte%boundary )
           IF ( electrolyte % boundary % update_status .EQ. 2 ) &
                & CALL print_environ_boundary( electrolyte%boundary )
           CALL update_environ_electrolyte( electrolyte )
           IF ( .NOT. electrolyte % update ) CALL print_environ_electrolyte( electrolyte )
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
! clean up subroutines of specific Environ modules.
!
     !
     ! Local clean up subroutines for the different contributions
     !
     IMPLICIT NONE
     !
     LOGICAL, INTENT(IN) :: lflag
     !
     CALL environ_clean_pw(lflag)
     !
     CALL environ_clean_tddfpt(lflag)
     !
     RETURN
     !
!--------------------------------------------------------------------
   END SUBROUTINE environ_clean
!--------------------------------------------------------------------
!--------------------------------------------------------------------
   SUBROUTINE environ_clean_pw(lflag)
!--------------------------------------------------------------------
!
! Clean up all the Environ related allocated variables, and call
! clean up subroutines of specific Environ modules. The structure of
! this subroutine should mirror the one of the environ_init...
! subroutines above.
!
     ! ... Declares modules
     !
     ! In environ_base all the control flags plus the variables
     ! that need to be deallocated
     !
     USE environ_base, ONLY : vzero, lelectrostatic, vreference,      &
                              lsoftcavity, vsoftcavity, lstatic,      &
                              static, charges, lexternals, externals, &
                              lconfine, vconfine,                     &
                              lelectrolyte, electrolyte,              &
                              lsemiconductor, semiconductor,          &
                              ions, electrons, system
     !
     ! Local clean up subroutines for the different contributions
     !
     IMPLICIT NONE
     !
     LOGICAL, INTENT(IN) :: lflag
     !
     ! ... Deallocate environment variables
     !
     IF ( ASSOCIATED( vzero%cell ) ) CALL destroy_environ_density( vzero )
     !
     ! ... environ_base variables
     !
     IF ( lelectrostatic .AND. ASSOCIATED( vreference%cell ) ) &
          & CALL destroy_environ_density( vreference )
     IF ( lsoftcavity .AND. ASSOCIATED( vsoftcavity%cell ) ) &
          & CALL destroy_environ_density( vsoftcavity )
     IF ( lconfine .AND. ASSOCIATED( vconfine%cell ) ) &
          & CALL destroy_environ_density( vconfine )
     !
     ! ... destroy derived types which were allocated in input
     !
     IF ( lelectrostatic .OR. lconfine ) CALL destroy_environ_charges( lflag, charges )
     IF ( lexternals ) CALL destroy_environ_externals( lflag, externals )
     IF ( lstatic ) CALL destroy_environ_dielectric( lflag, static )
     IF ( lelectrolyte ) CALL destroy_environ_electrolyte( lflag, electrolyte )
     IF ( lsemiconductor ) CALL destroy_environ_semiconductor( lflag, semiconductor)
     !
     CALL destroy_environ_electrons( lflag, electrons )
     CALL destroy_environ_ions( lflag, ions )
     CALL destroy_environ_system( lflag, system )
     !
     RETURN
     !
!--------------------------------------------------------------------
   END SUBROUTINE environ_clean_pw
!--------------------------------------------------------------------
!--------------------------------------------------------------------
   SUBROUTINE environ_clean_tddfpt(lflag)
!--------------------------------------------------------------------
!
! Clean up all the Environ related allocated variables, and call
! clean up subroutines of specific Environ modules. These are quantities
! that may be needed by TDDFPT, thus may need to be cleaned later
!
     ! ... Declares modules
     !
     ! In environ_base all the control flags plus the variables
     ! that need to be deallocated
     !
     USE environ_base, ONLY : lelectrostatic, velectrostatic,      &
          loptical, optical, lsolvent, solvent
     !
     USE electrostatic_init, ONLY : electrostatic_clean
     !
     IMPLICIT NONE
     !
     LOGICAL, INTENT(IN) :: lflag
     !
     LOGICAL :: opnd
     !
     IF ( lflag ) THEN
        INQUIRE( unit=environ_unit, opened= opnd )
        IF ( opnd ) CLOSE(unit=environ_unit)
     END IF
     !
     ! ... environ_base variables
     !
     IF ( lelectrostatic ) THEN
        IF ( ASSOCIATED( velectrostatic%cell ) ) &
             & CALL destroy_environ_density( velectrostatic )
        CALL electrostatic_clean( lflag )
     END IF
     !
     ! ... destroy derived types which were allocated in input
     !
     IF ( loptical ) CALL destroy_environ_dielectric( lflag, optical )
     IF ( lsolvent ) CALL destroy_environ_boundary( lflag, solvent )
     !
     IF ( lboundary ) CALL destroy_boundary_core( lflag, derivatives )
     !
     RETURN
     !
!--------------------------------------------------------------------
   END SUBROUTINE environ_clean_tddfpt
!--------------------------------------------------------------------
!----------------------------------------------------------------------------
END MODULE environ_init
!----------------------------------------------------------------------------
