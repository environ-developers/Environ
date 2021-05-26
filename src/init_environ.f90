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
!----------------------------------------------------------------------------------------
!
! Authors: Oliviero Andreussi (Department of Physics, UNT)
!          Francesco Nattino  (THEOS and NCCR-MARVEL, EPFL)
!          Ismaila Dabo       (DMSE, Penn State)
!          Nicola Marzari     (THEOS and NCCR-MARVEL, EPFL)
!
!----------------------------------------------------------------------------------------
!>
!! Module to initilize environ-related variables
!!
!----------------------------------------------------------------------------------------
MODULE init_environ
    !------------------------------------------------------------------------------------
    !
    USE environ_param, ONLY: DP, e2, bohr_radius_si, rydberg_si
    !
    USE base_environ
    USE base_core, ONLY: lfd, fd, lfft_environment, environment_fft
    !
    USE base_electrostatic, ONLY: need_pbc_correction, need_gradient, &
                                  need_factsqrt, need_auxiliary, need_electrolyte, &
                                  need_semiconductor, need_outer_loop
    !
    USE utils_electrons
    USE utils_ions
    USE utils_system
    USE utils_externals
    USE utils_boundary
    USE utils_semiconductor
    USE utils_electrolyte
    USE utils_dielectric
    USE utils_mapping
    USE utils_cell
    USE utils_charges
    USE utils_core_container
    !
    USE utils_density, ONLY: create_environ_density, init_environ_density, &
                             destroy_environ_density
    !
    USE tools_mapping, ONLY: map_small_to_large
    !
    USE init_core, ONLY: core_initbase, core_initcell, core_initions, core_clean
    USE init_electrostatic, ONLY: electrostatic_clean
    !
    USE environ_output, ONLY: environ_unit, print_environ_ions, print_environ_system, &
                              print_environ_boundary, print_environ_dielectric, &
                              print_environ_electrolyte, print_environ_electrons
    !
    !------------------------------------------------------------------------------------
    !
    PRIVATE
    !
    PUBLIC :: set_environ_base, environ_initbase, environ_initpotential, &
              environ_initcell, environ_initions, environ_initelectrons, &
              environ_initresponse, environ_clean, environ_clean_pw, &
              environ_clean_tddfpt
    !
    !------------------------------------------------------------------------------------
CONTAINS
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE set_environ_base(prog, nelec, nat, ntyp, atom_label, &
                                atomicspread, corespread, solvationrad, &
                                environ_restart_, environ_thr_, environ_nskip_, &
                                environ_type, system_ntyp, system_dim, system_axis, &
                                env_nrep, stype, rhomax, rhomin, tbeta, &
                                env_static_permittivity_, env_optical_permittivity_, &
                                solvent_mode, derivatives_, radius_mode, alpha, &
                                softness, solvent_distance, solvent_spread, &
                                solvent_radius, radial_scale, radial_spread, &
                                filling_threshold, filling_spread, field_awareness, &
                                charge_asymmetry, field_max, field_min, &
                                env_surface_tension_, env_pressure_, env_confine_, &
                                env_electrolyte_ntyp_, electrolyte_linearized, &
                                electrolyte_entropy, electrolyte_mode, &
                                electrolyte_distance, electrolyte_spread, cion, &
                                cionmax, rion, zion, electrolyte_rhomax, &
                                electrolyte_rhomin, electrolyte_tbeta, &
                                electrolyte_alpha, electrolyte_softness, temperature, &
                                sc_permittivity, sc_carrier_density, sc_electrode_chg, &
                                sc_distance, sc_spread, sc_chg_thr, &
                                env_external_charges, extcharge_charge, extcharge_dim, &
                                extcharge_axis, extcharge_pos, extcharge_spread, &
                                env_dielectric_regions, epsregion_eps, epsregion_dim, &
                                epsregion_axis, epsregion_pos, epsregion_spread, &
                                epsregion_width)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        LOGICAL, INTENT(IN) :: environ_restart_, electrolyte_linearized
        !
        INTEGER, INTENT(IN) :: nelec, nat, ntyp, &
                               environ_nskip_, &
                               system_ntyp, system_dim, system_axis, &
                               env_nrep(3), &
                               stype, env_electrolyte_ntyp_, &
                               env_external_charges, &
                               extcharge_dim(:), extcharge_axis(:), &
                               env_dielectric_regions, &
                               epsregion_dim(:), epsregion_axis(:)
        !
        REAL(DP), INTENT(IN) :: environ_thr_, rhomax, &
                                rhomin, tbeta, &
                                env_static_permittivity_, &
                                env_optical_permittivity_, &
                                alpha, softness, &
                                solvent_radius, radial_scale, radial_spread, &
                                filling_threshold, filling_spread, &
                                field_awareness, charge_asymmetry, field_max, &
                                field_min, &
                                solvationrad(:), corespread(:), atomicspread(:), &
                                solvent_distance, solvent_spread, &
                                env_surface_tension_, env_pressure_, &
                                env_confine_, &
                                electrolyte_distance, electrolyte_spread, &
                                cion(:), cionmax, rion, zion(:), &
                                electrolyte_rhomax, electrolyte_rhomin, &
                                electrolyte_tbeta, electrolyte_alpha, &
                                electrolyte_softness, temperature, &
                                sc_permittivity, sc_carrier_density, sc_electrode_chg, &
                                sc_distance, sc_spread, sc_chg_thr, &
                                extcharge_charge(:), extcharge_spread(:), &
                                extcharge_pos(:, :), epsregion_eps(:, :), &
                                epsregion_pos(:, :), epsregion_spread(:), &
                                epsregion_width(:)
        !
        CHARACTER(LEN=*), INTENT(IN) :: prog, environ_type, &
                                        solvent_mode, radius_mode, &
                                        electrolyte_mode, electrolyte_entropy, &
                                        derivatives_
        !
        CHARACTER(LEN=3), INTENT(IN) :: atom_label(:)
        !
        CHARACTER(LEN=80) :: label
        INTEGER :: i
        !
        CHARACTER(LEN=20) :: sub_name = 'set_environ_base'
        !
        !--------------------------------------------------------------------------------
        ! TDDFPT flag
        !
        SELECT CASE (prog(1:2))
        CASE ('TD')
            ltddfpt = .TRUE.
        CASE DEFAULT
            ltddfpt = .FALSE.
        END SELECT
        !
        !--------------------------------------------------------------------------------
        ! Create necessary local types
        !
        CALL create_environ_electrons(system_electrons)
        !
        CALL create_environ_ions(system_ions)
        !
        CALL create_environ_system(system_system)
        !
        CALL create_environ_electrons(environment_electrons)
        !
        CALL create_environ_ions(environment_ions)
        !
        CALL create_environ_system(environment_system)
        !
        !--------------------------------------------------------------------------------
        ! General flags
        !
        environ_restart = environ_restart_
        environ_thr = environ_thr_
        environ_nskip = environ_nskip_
        !
        !--------------------------------------------------------------------------------
        ! Set main environment flags, convert to internal units
        !
        env_static_permittivity = env_static_permittivity_
        env_optical_permittivity = env_optical_permittivity_
        !
        env_surface_tension = &
            env_surface_tension_ * 1.D-3 * bohr_radius_si**2 / rydberg_si
        !
        env_pressure = env_pressure_ * 1.D9 / rydberg_si * bohr_radius_si**3
        env_confine = env_confine_
        env_electrolyte_ntyp = env_electrolyte_ntyp_
        !
        !--------------------------------------------------------------------------------
        ! Set basic logical flags
        !
        lstatic = env_static_permittivity > 1.D0
        loptical = env_optical_permittivity > 1.D0
        !
        IF (env_dielectric_regions > 0) THEN
            !
            DO i = 1, env_dielectric_regions
                lstatic = lstatic .OR. (epsregion_eps(1, i) > 1.D0)
                loptical = loptical .OR. (epsregion_eps(2, i) > 1.D0)
            END DO
            !
        END IF
        !
        lsurface = env_surface_tension > 0.D0
        lvolume = env_pressure /= 0.D0
        lconfine = env_confine /= 0.D0
        lexternals = env_external_charges > 0
        lelectrolyte = env_electrolyte_ntyp > 0 .OR. need_electrolyte
        lsemiconductor = need_semiconductor
        lperiodic = need_pbc_correction
        ldoublecell = SUM(env_nrep) > 0
        louterloop = need_outer_loop
        !
        !--------------------------------------------------------------------------------
        ! Derived flags
        !
        ldielectric = lstatic .OR. loptical
        lsolvent = ldielectric .OR. lsurface .OR. lvolume .OR. lconfine
        lelectrostatic = ldielectric .OR. lelectrolyte .OR. lexternals .OR. lperiodic
        !
        lsoftsolvent = lsolvent .AND. (solvent_mode == 'electronic' .OR. &
                                       solvent_mode == 'full' .OR. &
                                       solvent_mode(1:2) == 'fa')
        !
        lsoftelectrolyte = lelectrolyte .AND. &
                           (electrolyte_mode == 'electronic' .OR. &
                            electrolyte_mode == 'full' .OR. &
                            electrolyte_mode(1:2) == 'fa') ! field-aware
        !
        lsoftcavity = lsoftsolvent .OR. lsoftelectrolyte
        lrigidsolvent = lsolvent .AND. solvent_mode /= 'electronic'
        lrigidelectrolyte = lelectrolyte .AND. electrolyte_mode /= 'electronic'
        lrigidcavity = lrigidsolvent .OR. lrigidelectrolyte
        !
        lcoredensity = (lsolvent .AND. solvent_mode == 'full') .OR. &
                       (lelectrolyte .AND. electrolyte_mode == 'full')
        !
        lsmearedions = lelectrostatic
        lboundary = lsolvent .OR. lelectrolyte
        lgradient = ldielectric .OR. (solvent_mode(1:2) == 'fa') ! field-aware
        !
        !--------------------------------------------------------------------------------
        ! Create optional types
        !
        IF (lexternals) CALL create_environ_externals(externals)
        !
        IF (lboundary) CALL create_core_container(derivatives)
        !
        IF (lsolvent) THEN
            label = 'solvent'
            !
            CALL create_environ_boundary(solvent, label)
            !
        END IF
        !
        IF (lelectrolyte) CALL create_environ_electrolyte(electrolyte)
        !
        IF (lsemiconductor) CALL create_environ_semiconductor(semiconductor)
        !
        IF (lstatic) CALL create_environ_dielectric(static)
        !
        !--------------------------------------------------------------------------------
        ! Set response properties for TD calculations
        !
        IF (loptical) THEN
            !
            CALL create_environ_dielectric(optical)
            !
            CALL create_environ_electrons(system_response_electrons)
            !
            CALL init_environ_electrons_first(0, system_response_electrons)
            !
            CALL create_environ_charges(system_response_charges)
            !
            CALL init_environ_charges_first(electrons=system_response_electrons, &
                                            charges=system_response_charges)
            !
            CALL create_environ_electrons(environment_response_electrons)
            !
            CALL init_environ_electrons_first(0, environment_response_electrons)
            !
            CALL create_environ_charges(environment_response_charges)
            !
            CALL init_environ_charges_first(electrons=environment_response_electrons, &
                                            dielectric=optical, &
                                            charges=environment_response_charges)
            !
        END IF
        !
        IF (lelectrostatic .OR. lconfine) THEN
            !
            CALL create_environ_charges(system_charges)
            !
            CALL create_environ_charges(environment_charges)
            !
        END IF
        !
        !--------------------------------------------------------------------------------
        ! Allocate and set basic properties of ions
        !
        CALL init_environ_ions_first(nat, ntyp, lsoftcavity, lcoredensity, &
                                     lsmearedions, radius_mode, atom_label, &
                                     atomicspread, corespread, solvationrad, &
                                     system_ions)
        !
        CALL init_environ_ions_first(nat, ntyp, lsoftcavity, lcoredensity, &
                                     lsmearedions, radius_mode, atom_label, &
                                     atomicspread, corespread, solvationrad, &
                                     environment_ions)
        !
        !--------------------------------------------------------------------------------
        ! Set basic properties of electrons
        !
        CALL init_environ_electrons_first(nelec, system_electrons)
        !
        CALL init_environ_electrons_first(nelec, environment_electrons)
        !
        !--------------------------------------------------------------------------------
        ! Set basic properties of the selected system
        !
        CALL init_environ_system(system_ntyp, system_dim, system_axis, &
                                 system_ions, system_system)
        !
        CALL init_environ_system(system_ntyp, system_dim, system_axis, &
                                 environment_ions, environment_system)
        !
        !--------------------------------------------------------------------------------
        ! Collect free charges if computing electrostatics or confinement
        !
        IF (lelectrostatic .OR. lconfine) THEN
            !
            CALL init_environ_charges_first(electrons=system_electrons, &
                                            charges=system_charges)
            !
            CALL init_environ_charges_first(electrons=environment_electrons, &
                                            charges=environment_charges)
            !
        END IF
        !
        IF (lelectrostatic) THEN
            !
            CALL init_environ_charges_first(ions=system_ions, charges=system_charges)
            !
            CALL init_environ_charges_first(ions=environment_ions, &
                                            charges=environment_charges)
            !
        END IF
        !
        !--------------------------------------------------------------------------------
        ! Allocate and set basic properties of external charges
        !
        IF (lexternals) THEN
            !
            CALL init_environ_externals_first(env_external_charges, extcharge_dim, &
                                              extcharge_axis, extcharge_pos, &
                                              extcharge_spread, extcharge_charge, &
                                              externals)
            !
            CALL init_environ_charges_first(externals=externals, &
                                            charges=environment_charges)
            !
        END IF
        !
        !--------------------------------------------------------------------------------
        ! Setup cores needed for derivatives of boundaries
        !
        IF (lboundary) THEN
            lfft_environment = .TRUE.
            !
            IF (derivatives_ == 'fd') lfd = .TRUE.
            !
            CALL init_core_container(derivatives_, derivatives, fft=environment_fft, &
                                     fd=fd)
            !
        END IF
        !
        !--------------------------------------------------------------------------------
        ! Set the parameters of the solvent boundary
        !
        IF (lsolvent) THEN
            !
            CALL init_environ_boundary_first(lgradient, need_factsqrt, lsurface, &
                                             solvent_mode, stype, rhomax, rhomin, &
                                             tbeta, env_static_permittivity, alpha, &
                                             softness, solvent_distance, &
                                             solvent_spread, solvent_radius, &
                                             radial_scale, radial_spread, &
                                             filling_threshold, filling_spread, &
                                             field_awareness, charge_asymmetry, &
                                             field_max, field_min, &
                                             environment_electrons, environment_ions, &
                                             environment_system, derivatives, solvent)
            !
        END IF
        !
        !--------------------------------------------------------------------------------
        ! Set the parameters of the electrolyte and of its boundary
        !
        IF (lelectrolyte) THEN
            !
            CALL init_environ_electrolyte_first(env_electrolyte_ntyp, &
                                                electrolyte_mode, stype, &
                                                electrolyte_rhomax, &
                                                electrolyte_rhomin, electrolyte_tbeta, &
                                                env_static_permittivity, &
                                                electrolyte_alpha, &
                                                electrolyte_softness, &
                                                electrolyte_distance, &
                                                electrolyte_spread, solvent_radius, &
                                                radial_scale, radial_spread, &
                                                filling_threshold, filling_spread, &
                                                field_awareness, charge_asymmetry, &
                                                field_max, field_min, &
                                                environment_electrons, &
                                                environment_ions, environment_system, &
                                                derivatives, temperature, cion, &
                                                cionmax, rion, zion, &
                                                electrolyte_entropy, &
                                                electrolyte_linearized, electrolyte)
            !
            CALL init_environ_charges_first(electrolyte=electrolyte, &
                                            charges=environment_charges)
            !
        END IF
        !
        !--------------------------------------------------------------------------------
        ! Set the parameters of the semiconductor
        !
        IF (lsemiconductor) THEN
            !
            CALL init_environ_semiconductor_first(temperature, sc_permittivity, &
                                                  sc_carrier_density, sc_electrode_chg, &
                                                  sc_distance, sc_spread, sc_chg_thr, &
                                                  environment_system, semiconductor)
            !
            CALL init_environ_charges_first(semiconductor=semiconductor, &
                                            charges=environment_charges)
            !
        END IF
        !
        !--------------------------------------------------------------------------------
        ! Set the parameters of the dielectric
        !
        IF (lstatic) THEN
            !
            CALL init_environ_dielectric_first(env_static_permittivity, solvent, &
                                               need_gradient, need_factsqrt, &
                                               need_auxiliary, static)
            !
            IF (env_dielectric_regions > 0) &
                CALL set_dielectric_regions(env_dielectric_regions, epsregion_dim, &
                                            epsregion_axis, epsregion_pos, &
                                            epsregion_width, epsregion_spread, &
                                            epsregion_eps(1, :), static)
            !
            CALL init_environ_charges_first(dielectric=static, &
                                            charges=environment_charges)
            !
        END IF
        !
        IF (loptical) THEN
            !
            CALL init_environ_dielectric_first(env_optical_permittivity, solvent, &
                                               need_gradient, need_factsqrt, &
                                               need_auxiliary, optical)
            !
            IF (env_dielectric_regions > 0) &
                CALL set_dielectric_regions(env_dielectric_regions, epsregion_dim, &
                                            epsregion_axis, epsregion_pos, &
                                            epsregion_width, epsregion_spread, &
                                            epsregion_eps(2, :), optical)
            !
        END IF
        !
        !--------------------------------------------------------------------------------
        ! Set the parameters for double cell mapping
        !
        IF (ldoublecell) ALLOCATE (environment_cell)
        !
        CALL init_environ_mapping_first(env_nrep, mapping)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE set_environ_base
    !------------------------------------------------------------------------------------
    !>
    !! Subroutine to initialize fundamental quantities needed by the
    !! environ modules. This subroutine is called by init_run.f90, thus
    !! only once per pw.x execution.
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE environ_initbase(alat, at, comm, me, root, gcutm, e2_in)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: comm, me, root
        REAL(DP), INTENT(IN) :: alat
        REAL(DP), INTENT(IN) :: at(3, 3)
        REAL(DP), INTENT(IN) :: gcutm
        REAL(DP), OPTIONAL, INTENT(IN) :: e2_in
        !
        INTEGER :: m(3)
        !
        INTEGER :: ipol
        INTEGER :: environment_nr(3)
        REAL(DP) :: environment_at(3, 3)
        !
        CHARACTER(LEN=80) :: label = ' '
        !
        !--------------------------------------------------------------------------------
        ! Allocate buffers used by env_mp_sum
        !
        CALL env_allocate_mp_buffers()
        !
        !--------------------------------------------------------------------------------
        ! Common initialization for simulations with Environ
        !
        e2 = 2.D0
        !
        IF (PRESENT(e2_in)) e2 = e2_in
        !
        CALL init_environ_cell(gcutm, comm, alat, at, system_cell) ! create system cell
        !
        !--------------------------------------------------------------------------------
        ! Double cell and mapping
        !
        IF (ldoublecell) THEN
            !
            !----------------------------------------------------------------------------
            ! Scale environment lattice (and corresponding ffts) by 2 * nrep(i) + 1
            !
            DO ipol = 1, 3
                !
                environment_at(:, ipol) = at(:, ipol) * &
                                          (2.D0 * mapping%nrep(ipol) + 1.D0)
                !
            END DO
            !
            environment_nr(1) = system_cell%dfft%nr1 * (2 * mapping%nrep(1) + 1)
            environment_nr(2) = system_cell%dfft%nr2 * (2 * mapping%nrep(2) + 1)
            environment_nr(3) = system_cell%dfft%nr3 * (2 * mapping%nrep(3) + 1)
            !
            !----------------------------------------------------------------------------
            ! Create environment cell
            !
            CALL init_environ_cell(gcutm, comm, alat, environment_at, &
                                   environment_cell, environment_nr)
            !
        ELSE
            environment_cell => system_cell
        END IF
        !
        CALL init_environ_mapping_second(system_cell, environment_cell, mapping)
        !
        !--------------------------------------------------------------------------------
        !
        CALL core_initbase(gcutm, environment_cell, system_cell)
        ! initialize numerical cores
        !
        !--------------------------------------------------------------------------------
        ! Create local storage for base potential, that needs to be modified
        !
        label = 'vzero'
        !
        CALL create_environ_density(vzero, label)
        !
        CALL init_environ_density(system_cell, vzero)
        !
        deenviron = 0.0_DP
        !
        label = 'dvtot'
        !
        CALL create_environ_density(dvtot, label)
        !
        CALL init_environ_density(system_cell, dvtot)
        !
        !--------------------------------------------------------------------------------
        ! Electrostatic contribution
        !
        eelectrostatic = 0.0_DP
        !
        IF (lelectrostatic) THEN
            label = 'velectrostatic'
            !
            CALL create_environ_density(velectrostatic, label)
            !
            CALL init_environ_density(environment_cell, velectrostatic)
            !
            label = 'vreference'
            !
            CALL create_environ_density(vreference, label)
            !
            CALL init_environ_density(system_cell, vreference)
            !
        END IF
        !
        !--------------------------------------------------------------------------------
        ! Contribution to the potential due to boundary
        !
        IF (lsoftcavity) THEN
            !
            label = 'vsoftcavity'
            !
            CALL create_environ_density(vsoftcavity, label)
            !
            CALL init_environ_density(environment_cell, vsoftcavity)
            !
        END IF
        !
        !--------------------------------------------------------------------------------
        ! Surface & volume contributions
        !
        esurface = 0.0_DP ! cavity contribution
        evolume = 0.0_DP ! pressure contribution
        !
        !--------------------------------------------------------------------------------
        ! Confinement contribution
        !
        econfine = 0.0_DP
        !
        IF (lconfine) THEN
            label = 'vconfine'
            !
            CALL create_environ_density(vconfine, label)
            !
            CALL init_environ_density(environment_cell, vconfine)
            !
        END IF
        !
        !--------------------------------------------------------------------------------
        !
        eelectrolyte = 0.0_DP ! non-electrostatic electrolyte contribution
        !
        !--------------------------------------------------------------------------------
        ! Second step of initialization of some environ derived type
        !
        CALL init_environ_electrons_second(system_cell, system_electrons)
        !
        CALL init_environ_electrons_second(environment_cell, environment_electrons)
        !
        IF (lsolvent) CALL init_environ_boundary_second(environment_cell, solvent)
        !
        IF (lstatic) CALL init_environ_dielectric_second(environment_cell, static)
        !
        IF (loptical) THEN
            !
            CALL init_environ_electrons_second(system_cell, system_response_electrons)
            !
            CALL init_environ_charges_second(system_cell, system_response_charges)
            !
            CALL init_environ_electrons_second(environment_cell, &
                                               environment_response_electrons)
            !
            CALL init_environ_charges_second(environment_cell, &
                                             environment_response_charges)
            !
            CALL init_environ_dielectric_second(environment_cell, optical)
            !
        END IF
        !
        IF (lelectrolyte) &
            CALL init_environ_electrolyte_second(environment_cell, electrolyte)
        !
        IF (lexternals) CALL init_environ_externals_second(environment_cell, externals)
        !
        IF (lelectrostatic .OR. lconfine) THEN
            !
            CALL init_environ_charges_second(system_cell, system_charges)
            !
            CALL init_environ_charges_second(environment_cell, environment_charges)
            !
        END IF
        !
        IF (lsemiconductor) &
            CALL init_environ_semiconductor_second(environment_cell, semiconductor)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE environ_initbase
    !------------------------------------------------------------------------------------
    !>
    !! Save local potential that will be overwritten by environ
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE environ_initpotential(nnr, vltot)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: nnr
        REAL(DP), INTENT(IN) :: vltot(nnr)
        !
        CHARACTER(LEN=80) :: sub_name = 'environ_initpotential'
        !
        !--------------------------------------------------------------------------------
        !
        vzero%update = .TRUE.
        !
        IF (.NOT. ASSOCIATED(vzero%cell)) RETURN
        !
        IF (vzero%cell%nnr /= nnr) &
            CALL env_errore(sub_name, 'Inconsistent size in input potential', 1)
        !
        vzero%of_r = vltot
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE environ_initpotential
    !------------------------------------------------------------------------------------
    !>
    !! Initialize the cell-related quantities to be used in the Environ
    !! modules. This initialization is called by electrons.f90, thus it
    !! is performed at every step of electronic optimization.
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE environ_initcell(at)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP), INTENT(IN) :: at(3, 3)
        !
        INTEGER :: ipol
        REAL(DP) :: environment_at(3, 3)
        !
        !--------------------------------------------------------------------------------
        !
        system_cell%update = .TRUE.
        !
        CALL update_environ_cell(at, system_cell) ! update system cell parameters
        !
        IF (ldoublecell) THEN
            environment_cell%update = .TRUE.
            !
            DO ipol = 1, 3
                !
                environment_at(:, ipol) = at(:, ipol) * &
                                          (2.D0 * mapping%nrep(ipol) + 1.D0)
                !
            END DO
            !
            CALL update_environ_cell(environment_at, environment_cell)
            ! update environment cell parameters
            !
        END IF
        !
        CALL core_initcell(system_cell, environment_cell) ! update cores
        !
        !--------------------------------------------------------------------------------
        ! Update fixed quantities defined inside the cell
        !
        IF (lstatic) CALL update_environ_dielectric(static)
        !
        IF (loptical) CALL update_environ_dielectric(optical)
        !
        IF (lelectrolyte) CALL update_environ_electrolyte(electrolyte)
        !
        IF (lexternals) CALL update_environ_externals(externals)
        !
        IF (lsemiconductor) CALL update_environ_semiconductor(semiconductor)
        !
        system_cell%update = .FALSE.
        !
        IF (ldoublecell) environment_cell%update = .FALSE.
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE environ_initcell
    !------------------------------------------------------------------------------------
    !>
    !! Initialize the ions-related quantities to be used in the Environ
    !! modules. This initialization is called by electrons.f90, thus it
    !! is performed at every step of electronic optimization. It may not
    !! be the most efficient choice, but it is a safe choice.
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE environ_initions(nnr, nat, ntyp, ityp, zv, tau, vloc)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: nnr, nat, ntyp
        INTEGER, INTENT(IN) :: ityp(nat)
        REAL(DP), INTENT(IN) :: zv(ntyp)
        REAL(DP), INTENT(IN) :: tau(3, nat)
        REAL(DP), INTENT(IN) :: vloc(nnr, ntyp)
        !
        INTEGER :: ia, it, dim, axis, icor, max_ntyp
        REAL(DP) :: charge, spread, dist, pos(3)
        REAL(DP), ALLOCATABLE :: aux(:, :)
        !
        !--------------------------------------------------------------------------------
        !
        system_ions%update = .TRUE.
        environment_ions%update = .TRUE.
        !
        !--------------------------------------------------------------------------------
        !
        ! second step of initialization for system ions
        CALL init_environ_ions_second(nat, ntyp, nnr, ityp, zv, system_cell, &
                                      system_ions, vloc)
        !
        !--------------------------------------------------------------------------------
        ! Update system ions parameters
        !
        CALL update_environ_ions(nat, tau, system_ions)
        !
        CALL print_environ_ions(system_ions)
        !
        !--------------------------------------------------------------------------------
        ! Update system system parameters
        !
        system_system%update = .TRUE.
        !
        CALL update_environ_system(system_system)
        !
        CALL print_environ_system(system_system)
        !
        !--------------------------------------------------------------------------------
        !
        CALL update_environ_mapping(mapping, system_system%pos)
        ! update mapping with correct shift of environment cell
        !
        !--------------------------------------------------------------------------------
        ! Second step of initialization for environment ions
        !
        IF (ldoublecell) THEN
            ALLOCATE (aux(environment_cell%nnr, ntyp))
            !
            DO it = 1, ntyp
                !
                CALL map_small_to_large(mapping, nnr, environment_cell%nnr, &
                                        vloc(:, it), aux(:, it))
                !
            END DO
            !
            CALL init_environ_ions_second(nat, ntyp, environment_cell%nnr, ityp, zv, &
                                          environment_cell, environment_ions, aux)
            !
            DEALLOCATE (aux)
        ELSE
            !
            CALL init_environ_ions_second(nat, ntyp, environment_cell%nnr, ityp, zv, &
                                          environment_cell, environment_ions, vloc)
            !
        END IF
        !
        !--------------------------------------------------------------------------------
        ! Update environment ions parameters
        !
        CALL update_environ_ions(nat, tau, environment_ions)
        !
        CALL print_environ_ions(environment_ions)
        !
        !--------------------------------------------------------------------------------
        ! Update environment system parameters
        !
        environment_system%update = .TRUE.
        !
        CALL update_environ_system(environment_system)
        !
        CALL print_environ_system(environment_system)
        !
        !--------------------------------------------------------------------------------
        ! Set soft-sphere parameters
        !
        IF (lsolvent) CALL set_soft_spheres(solvent)
        !
        IF (lelectrolyte) CALL set_soft_spheres(electrolyte%boundary)
        !
        !--------------------------------------------------------------------------------
        !
        CALL core_initions(environment_system%pos) ! update cores
        !
        !--------------------------------------------------------------------------------
        ! Update rigid environ properties, defined on ions
        !
        IF (lrigidcavity) THEN
            !
            IF (lsolvent) THEN
                !
                CALL update_environ_boundary(solvent)
                !
                IF (solvent%update_status == 2) CALL print_environ_boundary(solvent)
                !
                !------------------------------------------------------------------------
                ! Update quantities that depend on the solvent boundary
                !
                IF (lstatic) THEN
                    !
                    CALL update_environ_dielectric(static)
                    !
                    IF (.NOT. static%update) CALL print_environ_dielectric(static)
                    !
                END IF
                !
                IF (loptical) THEN
                    !
                    CALL update_environ_dielectric(optical)
                    !
                    IF (.NOT. optical%update) CALL print_environ_dielectric(optical)
                    !
                END IF
                !
            END IF
            !
            IF (lelectrolyte) THEN
                !
                CALL update_environ_boundary(electrolyte%boundary)
                !
                IF (electrolyte%boundary%update_status == 2) &
                    CALL print_environ_boundary(electrolyte%boundary)
                !
                CALL update_environ_electrolyte(electrolyte)
                !
                IF (.NOT. electrolyte%update) &
                    CALL print_environ_electrolyte(electrolyte)
                !
            END IF
            !
        END IF
        !
        !--------------------------------------------------------------------------------
        ! External charges rely on the environment cell, which is defined
        ! with respect to the system origin
        !
        IF (lexternals) CALL update_environ_externals(externals)
        !
        IF (lelectrostatic .OR. lconfine) THEN
            !
            CALL update_environ_charges(system_charges)
            !
            CALL update_environ_charges(environment_charges)
            !
        END IF
        !
        system_system%update = .FALSE.
        system_ions%update = .FALSE.
        environment_system%update = .FALSE.
        environment_ions%update = .FALSE.
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE environ_initions
    !------------------------------------------------------------------------------------
    !>
    !! Initialize the electrons-related quantities to be used in the Environ
    !! modules. This initialization is called by electrons.f90, thus it
    !! is performed at every step of electronic optimization.
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE environ_initelectrons(nnr, rho, nelec)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: nnr
        REAL(DP), INTENT(IN) :: rho(nnr)
        REAL(DP), INTENT(IN), OPTIONAL :: nelec
        !
        REAL(DP), ALLOCATABLE :: aux(:)
        !
        !--------------------------------------------------------------------------------
        !
        system_electrons%update = .TRUE.
        environment_electrons%update = .TRUE.
        !
        !--------------------------------------------------------------------------------
        ! Update electrons parameters
        !
        CALL update_environ_electrons(nnr, rho, system_electrons, nelec)
        !
        system_electrons%density%label = 'small_electrons'
        !
        ! CALL write_cube(system_electrons%density) ! DEBUGGING
        !
        IF (ldoublecell) THEN
            ALLOCATE (aux(environment_cell%nnr))
            !
            CALL map_small_to_large(mapping, nnr, environment_cell%nnr, rho, aux)
            !
            CALL update_environ_electrons(environment_cell%nnr, aux, &
                                          environment_electrons, nelec)
            !
        ELSE
            CALL update_environ_electrons(nnr, rho, environment_electrons, nelec)
        END IF
        !
        environment_electrons%density%label = 'large_electrons'
        !
        ! CALL write_cube(environment_electrons%density) ! DEBUGGING
        !
        ! STOP ! DEBUGGING
        !
        CALL print_environ_electrons(system_electrons)
        !
        CALL print_environ_electrons(environment_electrons)
        !
        IF (lelectrostatic .OR. lconfine) THEN
            !
            CALL update_environ_charges(system_charges)
            !
            CALL update_environ_charges(environment_charges)
            !
        END IF
        !
        !--------------------------------------------------------------------------------
        ! Update soft environ properties, defined on electrons
        !
        IF (lsoftcavity) THEN
            !
            IF (lsoftsolvent) THEN
                !
                CALL update_environ_boundary(solvent)
                !
                IF (solvent%update_status == 2) CALL print_environ_boundary(solvent)
                !
                !------------------------------------------------------------------------
                ! Update quantities that depend on the solvent boundary
                !
                IF (lstatic) THEN
                    !
                    CALL update_environ_dielectric(static)
                    !
                    IF (.NOT. static%update) CALL print_environ_dielectric(static)
                    !
                END IF
                !
                IF (loptical) THEN
                    !
                    CALL update_environ_dielectric(optical)
                    !
                    IF (.NOT. optical%update) CALL print_environ_dielectric(optical)
                    !
                END IF
                !
            END IF
            !
            IF (lsoftelectrolyte) THEN
                !
                CALL update_environ_boundary(electrolyte%boundary)
                !
                IF (electrolyte%boundary%update_status == 2) &
                    CALL print_environ_boundary(electrolyte%boundary)
                !
                CALL update_environ_electrolyte(electrolyte)
                !
                IF (.NOT. electrolyte%update) &
                    CALL print_environ_electrolyte(electrolyte)
                !
            END IF
            !
        END IF
        !
        system_electrons%update = .FALSE.
        environment_electrons%update = .FALSE.
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE environ_initelectrons
    !------------------------------------------------------------------------------------
    !>
    !! Initialize the response charges to be used in the TDDFPT + Environ
    !! modules. This initialization is called by plugin_tddfpt_potential.f90
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE environ_initresponse(nnr, drho)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: nnr
        REAL(DP), INTENT(IN) :: drho(nnr)
        !
        REAL(DP), ALLOCATABLE :: aux(:)
        !
        !--------------------------------------------------------------------------------
        !
        system_response_electrons%update = .TRUE.
        environment_response_electrons%update = .TRUE.
        !
        !--------------------------------------------------------------------------------
        ! Update response charges in system cell
        !
        CALL update_environ_electrons(nnr, drho, system_response_electrons, 0.D0)
        !
        CALL update_environ_charges(system_response_charges)
        !
        !--------------------------------------------------------------------------------
        ! Update response charges in environment cell
        !
        IF (ldoublecell) THEN
            !
            ALLOCATE (aux(environment_cell%nnr))
            !
            CALL map_small_to_large(mapping, nnr, environment_cell%nnr, drho, aux)
            !
            CALL update_environ_electrons(environment_cell%nnr, aux, &
                                          environment_response_electrons, 0.D0)
            !
        ELSE
            !
            CALL update_environ_electrons(nnr, drho, environment_response_electrons, &
                                          0.D0)
            !
        END IF
        !
        CALL update_environ_charges(environment_response_charges)
        !
        system_response_electrons%update = .FALSE.
        environment_response_electrons%update = .FALSE.
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE environ_initresponse
    !------------------------------------------------------------------------------------
    !>
    !! Clean up all the Environ related allocated variables, and call
    !! clean up subroutines of specific Environ modules.
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE environ_clean(lflag)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        LOGICAL, INTENT(IN) :: lflag
        !
        !--------------------------------------------------------------------------------
        !
        CALL environ_clean_pw(lflag)
        !
        CALL environ_clean_tddfpt(lflag)
        !
        CALL env_deallocate_mp_buffers()
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE environ_clean
    !------------------------------------------------------------------------------------
    !>
    !! Clean up all the Environ related allocated variables, and call clean up
    !! subroutines of specific Environ modules.
    !!
    !! The structure of this subroutine mirrors the one of init_environ subroutines
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE environ_clean_pw(lflag)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        LOGICAL, INTENT(IN) :: lflag
        !
        !--------------------------------------------------------------------------------
        ! Deallocate environment variables
        !
        IF (ASSOCIATED(vzero%cell)) CALL destroy_environ_density(vzero)
        !
        IF (ASSOCIATED(dvtot%cell)) CALL destroy_environ_density(dvtot)
        !
        !--------------------------------------------------------------------------------
        ! base_environ variables
        !
        IF (lelectrostatic .AND. ASSOCIATED(vreference%cell)) &
            CALL destroy_environ_density(vreference)
        !
        IF (lsoftcavity .AND. ASSOCIATED(vsoftcavity%cell)) &
            CALL destroy_environ_density(vsoftcavity)
        !
        IF (lconfine .AND. ASSOCIATED(vconfine%cell)) &
            CALL destroy_environ_density(vconfine)
        !
        !--------------------------------------------------------------------------------
        ! Destroy derived types which were allocated in input
        !
        IF (lelectrostatic .OR. lconfine) THEN
            !
            CALL destroy_environ_charges(lflag, system_charges)
            !
            CALL destroy_environ_charges(lflag, environment_charges)
            !
        END IF
        !
        IF (lexternals) CALL destroy_environ_externals(lflag, externals)
        !
        IF (lstatic) CALL destroy_environ_dielectric(lflag, static)
        !
        IF (lelectrolyte) CALL destroy_environ_electrolyte(lflag, electrolyte)
        !
        IF (lsemiconductor) CALL destroy_environ_semiconductor(lflag, semiconductor)
        !
        CALL destroy_environ_electrons(lflag, system_electrons)
        !
        CALL destroy_environ_ions(lflag, system_ions)
        !
        CALL destroy_environ_system(lflag, system_system)
        !
        CALL destroy_environ_electrons(lflag, environment_electrons)
        !
        CALL destroy_environ_ions(lflag, environment_ions)
        !
        CALL destroy_environ_system(lflag, environment_system)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE environ_clean_pw
    !------------------------------------------------------------------------------------
    !>
    !! Clean up all the Environ-related allocated variables and call clean up
    !! subroutines of specific Environ modules. These are quantities that may
    !! be needed by TDDFPT, thus may need to be cleaned later
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE environ_clean_tddfpt(lflag)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        LOGICAL, INTENT(IN) :: lflag
        !
        LOGICAL :: opnd
        !
        !--------------------------------------------------------------------------------
        !
        IF (lflag) THEN
            INQUIRE (unit=environ_unit, opened=opnd)
            !
            IF (opnd) CLOSE (unit=environ_unit)
            !
        END IF
        !
        !--------------------------------------------------------------------------------
        ! base_environ variables
        !
        IF (lelectrostatic) THEN
            !
            IF (ASSOCIATED(velectrostatic%cell)) &
                CALL destroy_environ_density(velectrostatic)
            !
            CALL electrostatic_clean(lflag)
            !
        END IF
        !
        !--------------------------------------------------------------------------------
        ! Destroy derived types which were allocated in input
        !
        IF (loptical) THEN
            !
            CALL destroy_environ_charges(lflag, environment_response_charges)
            !
            CALL destroy_environ_electrons(lflag, environment_response_electrons)
            !
            CALL destroy_environ_charges(lflag, system_response_charges)
            !
            CALL destroy_environ_electrons(lflag, system_response_electrons)
            !
            CALL destroy_environ_dielectric(lflag, optical)
            !
        END IF
        !
        IF (lsolvent) CALL destroy_environ_boundary(lflag, solvent)
        !
        IF (lboundary) CALL destroy_core_container(lflag, derivatives)
        !
        CALL core_clean(lflag)
        !
        IF (ldoublecell) THEN
            !
            CALL destroy_environ_mapping(lflag, mapping)
            !
            CALL destroy_environ_cell(lflag, environment_cell)
            !
        END IF
        !
        CALL destroy_environ_cell(lflag, system_cell)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE environ_clean_tddfpt
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE init_environ
!----------------------------------------------------------------------------------------
