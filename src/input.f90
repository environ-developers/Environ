!----------------------------------------------------------------------------------------
!
! Copyright (C) 2018-2022 ENVIRON (www.quantum-environ.org)
!
!----------------------------------------------------------------------------------------
!
!     This file is part of Environ version 3.0
!
!     Environ 3.0 is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation, either version 2 of the License, or
!     (at your option) any later version.
!
!     Environ 3.0 is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more detail, either the file
!     `License' in the root directory of the present distribution, or
!     online at <http://www.gnu.org/licenses/>.
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
!! This module contains all routines performing initialization and broadcast
!!
!----------------------------------------------------------------------------------------
MODULE environ_input
    !------------------------------------------------------------------------------------
    !
    USE class_io, ONLY: io
    USE env_mp, ONLY: env_mp_bcast
    USE env_char_ops, ONLY: env_uppercase, env_is_substring
    !
    USE environ_param, ONLY: DP, BOHR_RADIUS_ANGS
    !
    USE env_base_input
    !
    !------------------------------------------------------------------------------------
    !
    PRIVATE
    !
    PUBLIC :: read_environ_input
    !
    !------------------------------------------------------------------------------------
    !
    INTEGER :: local_nsx = 10 ! default maximum of ion types (used to allocate arrays)
    !
    !------------------------------------------------------------------------------------
CONTAINS
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !                                   MAIN ROUTINES
    !
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !>
    !! Routine for reading Environ input files. Uses built-in Namelist functionality
    !! and derived routines for cards (external charges and dielectric regions)
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE read_environ_input(filename, nsx)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: filename
        INTEGER, OPTIONAL, INTENT(IN) :: nsx
        !
        LOGICAL :: ext
        INTEGER :: environ_unit_input
        !
        CHARACTER(LEN=80) :: local_filename = 'environ.in'
        !
        CHARACTER(LEN=80) :: routine = 'read_environ_input'
        !
        !--------------------------------------------------------------------------------
        ! Open environ input file: environ.in
        !
        IF (PRESENT(filename)) local_filename = filename
        !
        IF (PRESENT(nsx)) local_nsx = nsx
        !
        environ_unit_input = io%find_free_unit()
        INQUIRE (file=TRIM(local_filename), exist=ext)
        !
        IF (.NOT. ext) CALL io%error(routine, "Missing input file", 1)
        !
        OPEN (unit=environ_unit_input, file=TRIM(local_filename), status="old")
        !
        !--------------------------------------------------------------------------------
        ! Read values into local variables
        !
        CALL io%header("Reading input from "//TRIM(local_filename)//":")
        !
        CALL environ_read_namelist(environ_unit_input)
        !
        CALL environ_read_cards(environ_unit_input)
        !
        CALL io%writer('') ! blank line after default settings
        !
        CLOSE (environ_unit_input)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE read_environ_input
    !------------------------------------------------------------------------------------
    !>
    !! Sets default values for all variables and overwrites with provided input
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE environ_read_namelist(environ_unit_input)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: environ_unit_input
        !
        INTEGER :: ios
        !
        CHARACTER(LEN=80) :: routine = 'environ_read_namelist'
        !
        !--------------------------------------------------------------------------------
        ! Set defaults
        !
        CALL allocate_registers()
        !
        CALL environ_defaults()
        !
        CALL boundary_defaults()
        !
        CALL electrostatic_defaults()
        !
        !--------------------------------------------------------------------------------
        ! &ENVIRON namelist
        !
        ios = 0
        !
        IF (io%lnode) READ (environ_unit_input, environ, iostat=ios)
        !
        CALL env_mp_bcast(ios, io%node, io%comm)
        !
        IF (ios /= 0) &
            CALL io%error(routine, "Missing or erroneous ENVIRON namelist", ABS(ios))
        !
        CALL environ_bcast() ! broadcast &ENVIRON variables
        !
        CALL environ_checkin() ! check &ENVIRON variables
        !
        !--------------------------------------------------------------------------------
        ! &BOUNDARY namelist (only if needed)
        !
        ios = 0
        !
        IF (io%lnode .AND. need_boundary()) &
            READ (environ_unit_input, boundary, iostat=ios)
        !
        CALL env_mp_bcast(ios, io%node, io%comm)
        !
        IF (ios /= 0) &
            CALL io%error(routine, "Missing or erroneous BOUNDARY namelist", ABS(ios))
        !
        CALL boundary_bcast() ! broadcast &BOUNDARY variables
        !
        CALL boundary_checkin() ! check &BOUNDARY variables
        !
        !--------------------------------------------------------------------------------
        ! &ELECTROSTATIC namelist (only if needed)
        !
        ios = 0
        !
        IF (io%lnode .AND. need_electrostatics()) &
            READ (environ_unit_input, electrostatic, iostat=ios)
        !
        CALL env_mp_bcast(ios, io%node, io%comm)
        !
        IF (ios /= 0) &
            CALL io%error(routine, &
                          "Missing or erroneous ELECTROSTATIC namelist", ABS(ios))
        !
        CALL electrostatic_bcast() ! broadcast &ELECTROSTATIC variables
        !
        CALL electrostatic_checkin() ! check &ELECTROSTATIC variables
        !
        !--------------------------------------------------------------------------------
        ! Final setup
        !
        CALL environ_type_setup() ! set up environment according to the boundary
        !
        CALL boundary_setup() ! correct missing/erroneous boundary settings
        !
        CALL electrostatics_setup() ! correct missing/erroneous electrostatics settings
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE environ_read_namelist
    !------------------------------------------------------------------------------------
    !>
    !! Environ cards parsing routine
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE environ_read_cards(unit)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, OPTIONAL, INTENT(IN) :: unit
        !
        INTEGER :: local_unit
        LOGICAL :: tend
        !
        CHARACTER(LEN=256) :: input_line
        CHARACTER(LEN=80) :: card
        !
        CHARACTER(LEN=80) :: routine = 'environ_read_cards'
        !
        !--------------------------------------------------------------------------------
        ! Set default READ unit if none provided
        !
        IF (PRESENT(unit)) THEN
            local_unit = unit
        ELSE
            local_unit = 5
        END IF
        !
        CALL env_read_line(local_unit, input_line, end_of_file=tend)
        !
        IF (.NOT. tend) THEN
            !
            READ (input_line, *) card ! store card title only in 'card'
            !
            input_line = env_uppercase(input_line) ! force uppercase
            !
            IF (TRIM(card) == 'EXTERNAL_CHARGES') THEN
                CALL card_external_charges(local_unit, input_line)
            ELSE IF (io%lnode) THEN
                CALL io%warning("card "//TRIM(input_line)//" ignored", 1001)
            END IF
            !
            IF (TRIM(card) == 'DIELECTRIC_REGIONS') THEN
                CALL card_dielectric_regions(local_unit, input_line)
            ELSE IF (io%lnode) THEN
                CALL io%warning("card "//TRIM(input_line)//" ignored", 1001)
            END IF
            !
        END IF
        !
        !--------------------------------------------------------------------------------
        ! Final check
        !
        IF (env_external_charges > 0 .AND. .NOT. taextchg) &
            CALL io%error(routine, "Missing card external_charges", 1)
        !
        IF (env_dielectric_regions > 0 .AND. .NOT. taepsreg) &
            CALL io%error(routine, "Missing card dielectric_regions", 1)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE environ_read_cards
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !                                  DEFAULT SETTINGS
    !
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !>
    !! Variables initialization for Namelist ENVIRON
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE environ_defaults()
        !--------------------------------------------------------------------------------
        !
        environ_debug = .FALSE.
        !
        environ_restart = .FALSE.
        verbose = 0
        environ_thr = 1.D-1
        environ_nskip = 1
        !
        env_electrostatic = .FALSE.
        no_electrostatics = .FALSE.
        !
        lsurface = .FALSE.
        lvolume = .FALSE.
        !
        env_ecut = 0.D0
        !
        environ_type = 'input'
        !
        system_ntyp = 0
        system_dim = 0
        system_axis = 3
        !
        env_nrep = 0
        system_pos = 0.D0
        !
        atomicspread = -0.5D0
        !
        env_static_permittivity = 1.D0
        env_optical_permittivity = 1.D0
        !
        env_surface_tension = 0.D0
        !
        env_pressure = 0.D0
        !
        env_confine = 0.D0
        !
        env_electrolyte_ntyp = 0
        electrolyte_linearized = .FALSE.
        electrolyte_entropy = 'full'
        cion = 1.0D0
        cionmax = 0.0D0 ! if remains zero, pb or linpb
        rion = 0.D0
        zion = 0.D0
        temperature = 300.0D0
        !
        sc_permittivity = 1.D0
        sc_carrier_density = 0.D0
        !
        sc_electrode_chg = 0.D0
        sc_chg_thr = 1.D-5
        !
        env_external_charges = 0
        env_dielectric_regions = 0
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE environ_defaults
    !------------------------------------------------------------------------------------
    !>
    !! Variables initialization for Namelist BOUNDARY
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE boundary_defaults()
        !--------------------------------------------------------------------------------
        !
        solvent_mode = 'electronic'
        !
        radius_mode = 'uff'
        alpha = 1.D0
        softness = 0.5D0
        solvationrad = -3.D0
        !
        rhomax = 0.005
        rhomin = 0.0001
        !
        corespread = -0.5D0
        !
        solvent_distance = 1.D0
        solvent_spread = 0.5D0
        !
        solvent_radius = 0.D0
        radial_scale = 2.D0
        radial_spread = 0.5D0
        filling_threshold = 0.825D0
        filling_spread = 0.02D0
        !
        field_aware = .FALSE.
        field_factor = 0.08D0
        field_asymmetry = -0.32D0
        field_max = 6.D0
        field_min = 2.D0
        !
        electrolyte_mode = 'electronic'
        !
        electrolyte_distance = 0.D0
        electrolyte_spread = 0.5D0
        !
        sc_distance = 0.D0
        sc_spread = 0.5D0
        !
        electrolyte_rhomax = 0.005D0
        electrolyte_rhomin = 0.0001D0
        !
        electrolyte_alpha = 1.D0
        electrolyte_softness = 0.5D0
        !
        deriv_method = 'default'
        deriv_core = 'fft'
        !
        electrolyte_deriv_method = 'default'
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE boundary_defaults
    !------------------------------------------------------------------------------------
    !>
    !! Variables initialization for Namelist ELECTROSTATIC
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE electrostatic_defaults()
        !--------------------------------------------------------------------------------
        !
        problem = 'none'
        tol = 1.D-5
        !
        solver = 'none'
        auxiliary = 'none'
        step_type = 'optimal'
        step = 0.3D0
        maxstep = 200
        inner_solver = 'none'
        inner_tol = 1.D-10
        inner_maxstep = 200
        inner_mix = 0.5D0
        !
        mix_type = 'linear'
        ndiis = 1
        mix = 0.5D0
        !
        preconditioner = 'sqrt'
        screening_type = 'none'
        screening = 0.D0
        !
        core = 'fft'
        inner_core = 'fft'
        !
        pbc_dim = -3
        pbc_correction = 'none'
        pbc_axis = 3
        pbc_core = '1da'
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE electrostatic_defaults
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !                                    BROADCASTING
    !
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !>
    !! Broadcast variables values for Namelist ENVIRON
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE environ_bcast()
        !--------------------------------------------------------------------------------
        !
        CALL env_mp_bcast(environ_debug, io%node, io%comm)
        !
        CALL env_mp_bcast(environ_restart, io%node, io%comm)
        !
        CALL env_mp_bcast(verbose, io%node, io%comm)
        !
        CALL env_mp_bcast(environ_thr, io%node, io%comm)
        !
        CALL env_mp_bcast(environ_nskip, io%node, io%comm)
        !
        CALL env_mp_bcast(env_ecut, io%node, io%comm)
        !
        CALL env_mp_bcast(environ_type, io%node, io%comm)
        !
        CALL env_mp_bcast(system_ntyp, io%node, io%comm)
        !
        CALL env_mp_bcast(system_dim, io%node, io%comm)
        !
        CALL env_mp_bcast(system_axis, io%node, io%comm)
        !
        CALL env_mp_bcast(env_nrep, io%node, io%comm)
        !
        CALL env_mp_bcast(system_pos, io%node, io%comm)
        !
        CALL env_mp_bcast(env_electrostatic, io%node, io%comm)
        !
        CALL env_mp_bcast(no_electrostatics, io%node, io%comm)
        !
        CALL env_mp_bcast(lvolume, io%node, io%comm)
        !
        CALL env_mp_bcast(lsurface, io%node, io%comm)
        !
        CALL env_mp_bcast(atomicspread, io%node, io%comm)
        !
        CALL env_mp_bcast(env_static_permittivity, io%node, io%comm)
        !
        CALL env_mp_bcast(env_optical_permittivity, io%node, io%comm)
        !
        CALL env_mp_bcast(env_surface_tension, io%node, io%comm)
        !
        CALL env_mp_bcast(env_pressure, io%node, io%comm)
        !
        CALL env_mp_bcast(env_confine, io%node, io%comm)
        !
        CALL env_mp_bcast(env_electrolyte_ntyp, io%node, io%comm)
        !
        CALL env_mp_bcast(electrolyte_linearized, io%node, io%comm)
        !
        CALL env_mp_bcast(electrolyte_entropy, io%node, io%comm)
        !
        CALL env_mp_bcast(cion, io%node, io%comm)
        !
        CALL env_mp_bcast(cionmax, io%node, io%comm)
        !
        CALL env_mp_bcast(rion, io%node, io%comm)
        !
        CALL env_mp_bcast(zion, io%node, io%comm)
        !
        CALL env_mp_bcast(temperature, io%node, io%comm)
        !
        CALL env_mp_bcast(sc_permittivity, io%node, io%comm)
        !
        CALL env_mp_bcast(sc_carrier_density, io%node, io%comm)
        !
        CALL env_mp_bcast(sc_electrode_chg, io%node, io%comm)
        !
        CALL env_mp_bcast(sc_chg_thr, io%node, io%comm)
        !
        CALL env_mp_bcast(env_external_charges, io%node, io%comm)
        !
        CALL env_mp_bcast(env_dielectric_regions, io%node, io%comm)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE environ_bcast
    !------------------------------------------------------------------------------------
    !>
    !! Broadcast variables values for Namelist BOUNDARY
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE boundary_bcast()
        !--------------------------------------------------------------------------------
        !
        CALL env_mp_bcast(solvent_mode, io%node, io%comm)
        !
        CALL env_mp_bcast(rhomax, io%node, io%comm)
        !
        CALL env_mp_bcast(rhomin, io%node, io%comm)
        !
        CALL env_mp_bcast(radius_mode, io%node, io%comm)
        !
        CALL env_mp_bcast(alpha, io%node, io%comm)
        !
        CALL env_mp_bcast(softness, io%node, io%comm)
        !
        CALL env_mp_bcast(solvationrad, io%node, io%comm)
        !
        CALL env_mp_bcast(corespread, io%node, io%comm)
        !
        CALL env_mp_bcast(solvent_distance, io%node, io%comm)
        !
        CALL env_mp_bcast(solvent_spread, io%node, io%comm)
        !
        CALL env_mp_bcast(solvent_radius, io%node, io%comm)
        !
        CALL env_mp_bcast(radial_scale, io%node, io%comm)
        !
        CALL env_mp_bcast(radial_spread, io%node, io%comm)
        !
        CALL env_mp_bcast(filling_threshold, io%node, io%comm)
        !
        CALL env_mp_bcast(filling_spread, io%node, io%comm)
        !
        CALL env_mp_bcast(field_aware, io%node, io%comm)
        !
        CALL env_mp_bcast(field_factor, io%node, io%comm)
        !
        CALL env_mp_bcast(field_asymmetry, io%node, io%comm)
        !
        CALL env_mp_bcast(field_max, io%node, io%comm)
        !
        CALL env_mp_bcast(field_min, io%node, io%comm)
        !
        CALL env_mp_bcast(electrolyte_mode, io%node, io%comm)
        !
        CALL env_mp_bcast(electrolyte_distance, io%node, io%comm)
        !
        CALL env_mp_bcast(electrolyte_spread, io%node, io%comm)
        !
        CALL env_mp_bcast(sc_distance, io%node, io%comm)
        !
        CALL env_mp_bcast(sc_spread, io%node, io%comm)
        !
        CALL env_mp_bcast(electrolyte_rhomax, io%node, io%comm)
        !
        CALL env_mp_bcast(electrolyte_rhomin, io%node, io%comm)
        !
        CALL env_mp_bcast(electrolyte_alpha, io%node, io%comm)
        !
        CALL env_mp_bcast(electrolyte_softness, io%node, io%comm)
        !
        CALL env_mp_bcast(deriv_method, io%node, io%comm)
        !
        CALL env_mp_bcast(deriv_core, io%node, io%comm)
        !
        CALL env_mp_bcast(electrolyte_deriv_method, io%node, io%comm)
        !
        CALL env_mp_bcast(deriv_lowpass_p1, io%node, io%comm)
        !
        CALL env_mp_bcast(deriv_lowpass_p2, io%node, io%comm)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE boundary_bcast
    !------------------------------------------------------------------------------------
    !>
    !! Broadcast variables values for Namelist ELECTROSTATIC
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE electrostatic_bcast()
        !--------------------------------------------------------------------------------
        !
        CALL env_mp_bcast(problem, io%node, io%comm)
        !
        CALL env_mp_bcast(tol, io%node, io%comm)
        !
        CALL env_mp_bcast(solver, io%node, io%comm)
        !
        CALL env_mp_bcast(inner_solver, io%node, io%comm)
        !
        CALL env_mp_bcast(inner_tol, io%node, io%comm)
        !
        CALL env_mp_bcast(inner_maxstep, io%node, io%comm)
        !
        CALL env_mp_bcast(inner_mix, io%node, io%comm)
        !
        CALL env_mp_bcast(auxiliary, io%node, io%comm)
        !
        CALL env_mp_bcast(step_type, io%node, io%comm)
        !
        CALL env_mp_bcast(step, io%node, io%comm)
        !
        CALL env_mp_bcast(maxstep, io%node, io%comm)
        !
        CALL env_mp_bcast(mix_type, io%node, io%comm)
        !
        CALL env_mp_bcast(mix, io%node, io%comm)
        !
        CALL env_mp_bcast(ndiis, io%node, io%comm)
        !
        CALL env_mp_bcast(preconditioner, io%node, io%comm)
        !
        CALL env_mp_bcast(screening_type, io%node, io%comm)
        !
        CALL env_mp_bcast(screening, io%node, io%comm)
        !
        CALL env_mp_bcast(core, io%node, io%comm)
        !
        CALL env_mp_bcast(inner_core, io%node, io%comm)
        !
        CALL env_mp_bcast(pbc_dim, io%node, io%comm)
        !
        CALL env_mp_bcast(pbc_correction, io%node, io%comm)
        !
        CALL env_mp_bcast(pbc_axis, io%node, io%comm)
        !
        CALL env_mp_bcast(pbc_core, io%node, io%comm)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE electrostatic_bcast
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !                                       FLAGS
    !
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !>
    !! Check if BOUNDARY namelist needs to be read according to the ENVIRON namelist
    !!
    !------------------------------------------------------------------------------------
    LOGICAL FUNCTION need_boundary()
        !--------------------------------------------------------------------------------
        !
        need_boundary = .FALSE.
        !
        IF (environ_type /= 'input' .AND. environ_type /= 'vacuum') need_boundary = .TRUE.
        !
        IF (env_static_permittivity > 1.D0) need_boundary = .TRUE.
        !
        IF (env_optical_permittivity > 1.D0) need_boundary = .TRUE.
        !
        IF (env_surface_tension /= 0.D0) need_boundary = .TRUE.
        !
        IF (env_pressure /= 0.D0) need_boundary = .TRUE.
        !
        IF (env_confine /= 0.D0) need_boundary = .TRUE.
        !
        IF (env_electrolyte_ntyp > 0) need_boundary = .TRUE.
        !
        IF (env_dielectric_regions > 0) need_boundary = .TRUE.
        !
        IF (sc_permittivity > 1.D0) need_boundary = .TRUE.
        !
        IF (sc_carrier_density > 0) need_boundary = .TRUE.
        !
        IF (lsurface .OR. lvolume) need_boundary = .TRUE.
        !
        !--------------------------------------------------------------------------------
    END FUNCTION need_boundary
    !------------------------------------------------------------------------------------
    !>
    !! Check if ELECTROSTATIC namelist needs to be read according to the ENVIRON namelist
    !!
    !------------------------------------------------------------------------------------
    LOGICAL FUNCTION need_electrostatics()
        !--------------------------------------------------------------------------------
        !
        need_electrostatics = env_electrostatic
        !
        IF (env_static_permittivity > 1.D0) need_electrostatics = .TRUE.
        !
        IF (env_optical_permittivity > 1.D0) need_electrostatics = .TRUE.
        !
        IF (env_external_charges > 0) need_electrostatics = .TRUE.
        !
        IF (env_dielectric_regions > 0) need_electrostatics = .TRUE.
        !
        IF (env_electrolyte_ntyp > 0) need_electrostatics = .TRUE.
        !
        IF (sc_permittivity > 1.D0) need_electrostatics = .TRUE.
        !
        IF (sc_carrier_density > 0) need_electrostatics = .TRUE.
        !
        IF (no_electrostatics) need_electrostatics = .FALSE.
        !
        !--------------------------------------------------------------------------------
    END FUNCTION need_electrostatics
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !                                  INPUT VALIDATION
    !
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !>
    !! Check input values for Namelist ENVIRON
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE environ_checkin()
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER :: i
        LOGICAL :: allowed
        !
        CHARACTER(LEN=80) :: routine = 'environ_checkin'
        !
        !--------------------------------------------------------------------------------
        ! General
        !
        IF (environ_restart) CALL io%writer("* Environ restarting")
        !
        IF (verbose < 0) CALL io%error(routine, "verbose out of range", 1)
        !
        IF (environ_thr < 0.0_DP) CALL io%error(routine, "environ_thr out of range", 1)
        !
        IF (environ_nskip < 0) CALL io%error(routine, "environ_nskip out of range", 1)
        !
        IF (env_nrep(1) < 0 .OR. env_nrep(2) < 0 .OR. env_nrep(3) < 0) &
            CALL io%error(routine, "env_nrep cannot be smaller than 0", 1)
        !
        !--------------------------------------------------------------------------------
        ! Type
        !
        allowed = .FALSE.
        !
        DO i = 1, SIZE(environ_type_allowed)
            IF (TRIM(environ_type) == environ_type_allowed(i)) allowed = .TRUE.
        END DO
        !
        IF (.NOT. allowed) CALL io%invalid_opt(routine, 'environ_type', environ_type)
        !
        !--------------------------------------------------------------------------------
        ! System
        !
        IF (system_ntyp < 0) CALL io%error(routine, "system_ntype out of range", 1)
        !
        IF (system_dim < 0 .OR. system_dim > 3) &
            CALL io%error(routine, "system_dim out of range", 1)
        !
        IF (system_axis < 1 .OR. system_axis > 3) &
            CALL io%error(routine, "system_axis out of range", 1)
        !
        !--------------------------------------------------------------------------------
        ! Physical quantities
        !
        IF (env_static_permittivity < 1.0_DP) &
            CALL io%error(routine, "env_static_permittivity out of range", 1)
        !
        IF (env_optical_permittivity < 1.0_DP) &
            CALL io%error(routine, "env_optical_permittivity out of range", 1)
        !
        IF (temperature < 0.0_DP) CALL io%error(routine, "temperature out of range", 1)
        !
        !--------------------------------------------------------------------------------
        ! Electrolyte
        !
        allowed = .FALSE.
        !
        DO i = 1, SIZE(electrolyte_entropy_allowed)
            !
            IF (TRIM(electrolyte_entropy) == electrolyte_entropy_allowed(i)) &
                allowed = .TRUE.
            !
        END DO
        !
        IF (.NOT. allowed) &
            CALL io%invalid_opt(routine, 'electrolyte_entropy', electrolyte_entropy)
        !
        IF (env_electrolyte_ntyp < 0 .OR. env_electrolyte_ntyp == 1) &
            CALL io%error(routine, "env_electrolyte_ntyp out of range", 1)
        !
        DO i = 1, env_electrolyte_ntyp
            IF (cion(i) < 0.D0) CALL io%error(routine, "cion cannot be negative", 1)
        END DO
        !
        IF (cionmax < 0.D0 .OR. rion < 0.D0) &
            CALL io%error(routine, "cionmax and rion cannot be negative", 1)
        !
        IF (cionmax > 0.D0 .AND. rion > 0.D0) &
            CALL io%error(routine, "Either cionmax or rion can be set, not both", 1)
        !
        !--------------------------------------------------------------------------------
        ! Semiconductor
        !
        allowed = .FALSE.
        !
        IF (sc_permittivity < 1.D0) &
            CALL io%error(routine, "sc_permittivity out of range", 1)
        !
        IF (sc_carrier_density < 0.D0) &
            CALL io%error(routine, "sc_carrier_density cannot be negative", 1)
        !
        !--------------------------------------------------------------------------------
        ! Externals/Regions
        !
        IF (env_external_charges < 0) &
            CALL io%error(routine, "env_external_charges out of range", 1)
        !
        IF (env_dielectric_regions < 0) &
            CALL io%error(routine, "env_dielectric_regions out of range", 1)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE environ_checkin
    !------------------------------------------------------------------------------------
    !>
    !! Check input values for Namelist BOUNDARY
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE boundary_checkin()
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER :: i
        LOGICAL :: allowed
        !
        CHARACTER(LEN=80) :: routine = 'boundary_checkin'
        !
        !--------------------------------------------------------------------------------
        ! Solvent
        !
        allowed = .FALSE.
        !
        DO i = 1, SIZE(solvent_mode_allowed)
            IF (TRIM(solvent_mode) == solvent_mode_allowed(i)) allowed = .TRUE.
        END DO
        !
        IF (.NOT. allowed) CALL io%invalid_opt(routine, 'solvent_mode', solvent_mode)
        !
        IF (rhomax < 0.0_DP) CALL io%error(routine, "rhomax out of range", 1)
        !
        IF (rhomin < 0.0_DP) CALL io%error(routine, "rhomin out of range", 1)
        !
        IF (rhomax < rhomin) CALL io%error(routine, "Inconsistent rhomax and rhomin", 1)
        !
        allowed = .FALSE.
        !
        DO i = 1, SIZE(radius_mode_allowed)
            IF (TRIM(radius_mode) == radius_mode_allowed(i)) allowed = .TRUE.
        END DO
        !
        IF (.NOT. allowed) CALL io%invalid_opt(routine, 'radius_mode', radius_mode)
        !
        IF (alpha <= 0.0_DP) CALL io%error(routine, "alpha out of range", 1)
        !
        IF (softness <= 0.0_DP) CALL io%error(routine, "softness out of range", 1)
        !
        IF (solvent_distance < 0.0_DP) &
            CALL io%error(routine, "solvent_distance out of range", 1)
        !
        IF (solvent_spread <= 0.0_DP) &
            CALL io%error(routine, "solvent_spread out of range", 1)
        !
        IF (solvent_radius < 0.0_DP) &
            CALL io%error(routine, "solvent_radius out of range", 1)
        !
        IF (radial_scale < 1.0_DP) &
            CALL io%error(routine, "radial_scale out of range", 1)
        !
        IF (radial_spread <= 0.0_DP) &
            CALL io%error(routine, "radial_spread out of range", 1)
        !
        IF (filling_threshold <= 0.0_DP) &
            CALL io%error(routine, "filling_threshold out of range", 1)
        !
        IF (filling_spread <= 0.0_DP) &
            CALL io%error(routine, "filling_spread out of range", 1)
        !
        IF (field_factor < 0.0_DP) &
            CALL io%error(routine, "field_factor out of range", 1)
        !
        IF (ABS(field_asymmetry) > 1.0_DP) &
            CALL io%error(routine, "field_asymmetry out of range", 1)
        !
        IF (field_min < 0.0_DP) CALL io%error(routine, "field_min out of range", 1)
        !
        IF (field_max <= field_min) CALL io%error(routine, "field_max out of range", 1)
        !
        !--------------------------------------------------------------------------------
        ! Electrolyte
        !
        allowed = .FALSE.
        !
        DO i = 1, SIZE(electrolyte_mode_allowed)
            IF (TRIM(electrolyte_mode) == electrolyte_mode_allowed(i)) allowed = .TRUE.
        END DO
        !
        IF (.NOT. allowed) &
            CALL io%invalid_opt(routine, 'electrolyte_mode', electrolyte_mode)
        !
        IF (electrolyte_distance < 0.0_DP) &
            CALL io%error(routine, "electrolyte_distance out of range", 1)
        !
        IF (electrolyte_spread <= 0.0_DP) &
            CALL io%error(routine, "electrolyte_spread out of range", 1)
        !
        IF (electrolyte_rhomax < 0.0_DP) &
            CALL io%error(routine, "electrolyte_rhomax out of range", 1)
        !
        IF (electrolyte_rhomin < 0.0_DP) &
            CALL io%error(routine, "electrolyte_rhomin out of range", 1)
        !
        IF (electrolyte_rhomax < electrolyte_rhomin) &
            CALL io%error(routine, &
                          "Inconsistent electrolyte_rhomax and electrolyte_rhomin", 1)
        !
        IF (electrolyte_alpha <= 0.0_DP) &
            CALL io%error(routine, "electrolyte_alpha out of range", 1)
        !
        IF (electrolyte_softness <= 0.0_DP) &
            CALL io%error(routine, "electrolyte_softness out of range", 1)
        !
        !--------------------------------------------------------------------------------
        ! Semiconductor
        !
        IF (sc_distance < 0.0_DP) CALL io%error(routine, "sc_distance out of range", 1)
        !
        IF (sc_spread <= 0.0_DP) CALL io%error(routine, "sc_spread out of range", 1)
        !
        !--------------------------------------------------------------------------------
        ! Derivatives
        !
        allowed = .FALSE.
        !
        DO i = 1, SIZE(deriv_method_allowed)
            IF (TRIM(deriv_method) == deriv_method_allowed(i)) allowed = .TRUE.
        END DO
        !
        IF (.NOT. allowed) &
            CALL io%invalid_opt(routine, 'deriv_method', deriv_method)
        !
        allowed = .FALSE.
        !
        DO i = 1, SIZE(deriv_core_allowed)
            IF (TRIM(deriv_core) == deriv_core_allowed(i)) allowed = .TRUE.
        END DO
        !
        IF (.NOT. allowed) &
            CALL io%invalid_opt(routine, 'deriv_core', deriv_method)
        !
        allowed = .FALSE.
        !
        DO i = 1, SIZE(electrolyte_deriv_method_allowed)
            !
            IF (TRIM(electrolyte_deriv_method) == electrolyte_deriv_method_allowed(i)) &
                allowed = .TRUE.
            !
        END DO
        !
        IF (.NOT. allowed) &
            CALL io%invalid_opt(routine, 'electrolyte_deriv_method', &
                                electrolyte_deriv_method)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE boundary_checkin
    !------------------------------------------------------------------------------------
    !>
    !! Check input values for Namelist ELECTROSTATIC
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE electrostatic_checkin()
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER :: i
        LOGICAL :: allowed
        !
        CHARACTER(LEN=80) :: routine = 'electrostatic_checkin'
        !
        !--------------------------------------------------------------------------------
        ! Problem
        !
        allowed = .FALSE.
        !
        DO i = 1, SIZE(problem_allowed)
            IF (TRIM(problem) == problem_allowed(i)) allowed = .TRUE.
        END DO
        !
        IF (.NOT. allowed) CALL io%invalid_opt(routine, 'problem', problem)
        !
        IF (tol <= 0.0_DP) CALL io%error(routine, "tolerance out of range", 1)
        !
        !--------------------------------------------------------------------------------
        ! Solver
        !
        allowed = .FALSE.
        !
        DO i = 1, SIZE(solver_allowed)
            IF (TRIM(solver) == solver_allowed(i)) allowed = .TRUE.
        END DO
        !
        IF (.NOT. allowed) CALL io%invalid_opt(routine, 'solver', solver)
        !
        !--------------------------------------------------------------------------------
        ! Auxiliary
        !
        allowed = .FALSE.
        !
        DO i = 1, SIZE(auxiliary_allowed)
            IF (TRIM(auxiliary) == auxiliary_allowed(i)) allowed = .TRUE.
        END DO
        !
        IF (.NOT. allowed) CALL io%invalid_opt(routine, 'auxiliary', auxiliary)
        !
        !--------------------------------------------------------------------------------
        ! Step
        !
        allowed = .FALSE.
        !
        DO i = 1, SIZE(step_type_allowed)
            IF (TRIM(step_type) == step_type_allowed(i)) allowed = .TRUE.
        END DO
        !
        IF (.NOT. allowed) CALL io%invalid_opt(routine, 'step_type', step_type)
        !
        IF (step <= 0.0_DP) CALL io%error(routine, "step out of range", 1)
        !
        IF (maxstep <= 1) CALL io%error(routine, "maxstep out of range", 1)
        !
        !--------------------------------------------------------------------------------
        ! Mixing
        !
        allowed = .FALSE.
        !
        DO i = 1, SIZE(mix_type_allowed)
            IF (TRIM(mix_type) == mix_type_allowed(i)) allowed = .TRUE.
        END DO
        !
        IF (.NOT. allowed) CALL io%invalid_opt(routine, 'mix_type', mix_type)
        !
        IF (ndiis <= 0) CALL io%error(routine, "ndiis out of range", 1)
        !
        IF (mix <= 0.0_DP) CALL io%error(routine, "mix out of range", 1)
        !
        !--------------------------------------------------------------------------------
        ! Preconditioner
        !
        allowed = .FALSE.
        !
        DO i = 1, SIZE(preconditioner_allowed)
            IF (TRIM(preconditioner) == preconditioner_allowed(i)) allowed = .TRUE.
        END DO
        !
        IF (.NOT. allowed) &
            CALL io%invalid_opt(routine, 'preconditioner', preconditioner)
        !
        !--------------------------------------------------------------------------------
        ! Screening
        !
        allowed = .FALSE.
        !
        DO i = 1, SIZE(screening_type_allowed)
            IF (TRIM(screening_type) == screening_type_allowed(i)) allowed = .TRUE.
        END DO
        !
        IF (.NOT. allowed) &
            CALL io%invalid_opt(routine, 'screening_type', screening_type)
        !
        IF (screening < 0.0_DP) CALL io%error(routine, "screening out of range", 1)
        !
        allowed = .FALSE.
        !
        !--------------------------------------------------------------------------------
        ! Outer core
        !
        DO i = 1, SIZE(core_allowed)
            IF (TRIM(core) == core_allowed(i)) allowed = .TRUE.
        END DO
        !
        IF (.NOT. allowed) CALL io%invalid_opt(routine, 'core', core)
        !
        !--------------------------------------------------------------------------------
        ! PBC correction
        !
        allowed = .FALSE.
        !
        DO i = 1, SIZE(pbc_correction_allowed)
            IF (TRIM(pbc_correction) == pbc_correction_allowed(i)) allowed = .TRUE.
        END DO
        !
        IF (.NOT. allowed) &
            CALL io%invalid_opt(routine, 'pbc_correction', pbc_correction)
        !
        IF (pbc_dim < -3 .OR. pbc_dim > 3) &
            CALL io%error(routine, "pbc_dim out of range", 1)
        !
        IF (TRIM(pbc_correction) /= 'none' .AND. pbc_dim < 0) &
            CALL io%error(routine, &
                          "pbc_correction requires manually setting pbc_dim", 1)
        !
        IF (pbc_dim == 1) CALL io%error(routine, "1d pbc correction not implemented", 1)
        !
        IF (pbc_axis < 1 .OR. pbc_axis > 3) &
            CALL io%error(routine, "cell_axis out of range", 1)
        !
        !--------------------------------------------------------------------------------
        ! PBC core
        !
        allowed = .FALSE.
        !
        DO i = 1, SIZE(pbc_core_allowed)
            IF (TRIM(pbc_core) == pbc_core_allowed(i)) allowed = .TRUE.
        END DO
        !
        IF (.NOT. allowed) CALL io%invalid_opt(routine, 'pbc_core', core)
        !
        !--------------------------------------------------------------------------------
        ! Inner solver
        !
        allowed = .FALSE.
        !
        DO i = 1, SIZE(inner_solver_allowed)
            IF (TRIM(inner_solver) == inner_solver_allowed(i)) allowed = .TRUE.
        END DO
        !
        IF (.NOT. allowed) CALL io%invalid_opt(routine, 'inner_solver', inner_solver)
        !
        IF (inner_mix <= 0.0_DP) CALL io%error(routine, "inner_mix out of range", 1)
        !
        IF (inner_tol <= 0.0_DP) CALL io%error(routine, "inner_tol out of range", 1)
        !
        IF (inner_maxstep <= 1) CALL io%error(routine, "inner_maxstep out of range", 1)
        !
        !--------------------------------------------------------------------------------
        ! Inner core
        !
        allowed = .FALSE.
        !
        DO i = 1, SIZE(inner_core_allowed)
            IF (TRIM(inner_core) == inner_core_allowed(i)) allowed = .TRUE.
        END DO
        !
        IF (.NOT. allowed) CALL io%invalid_opt(routine, 'inner_core', core)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE electrostatic_checkin
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE boundary_setup()
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CHARACTER(LEN=80) :: routine = 'boundary setup'
        !
        !--------------------------------------------------------------------------------
        !
        CALL io%writer("* using "//TRIM(solvent_mode)//" solvent mode")
        !
        SELECT CASE (TRIM(solvent_mode))
            !
        CASE ('electronic', 'full', 'system')
            !
            SELECT CASE (TRIM(deriv_method))
                !
            CASE ('default')
                deriv_method = 'chain'
                !
                CALL io%writer("* setting boundary derivatives method to SCCS default")
                !
            CASE ('highmem', 'lowmem')
                !
                CALL io%error(routine, &
                              "Only 'fft' or 'chain' are allowed with electronic interfaces", 1)
                !
            END SELECT
            !
        CASE ('ionic')
            !
            SELECT CASE (TRIM(deriv_method))
                !
            CASE ('default')
                deriv_method = 'lowmem'
                !
                CALL io%writer("* setting boundary derivatives method to SSCS default")
                !
            CASE ('chain')
                !
                CALL io%error(routine, &
                              "Only 'highmem', 'lowmem', and 'fft' are allowed with ionic interfaces", 1)
                !
            END SELECT
            !
        END SELECT
        !
        IF (solvent_mode == 'system' .AND. solvent_distance == 0.0_DP) &
            CALL io%error(routine, &
                          "solvent_distance must be set (greater than zero) for system interfaces", 1)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE boundary_setup
    !------------------------------------------------------------------------------------
    !>
    !! Set values according to the environ_type keyword and boundary mode
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE environ_type_setup()
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CHARACTER(LEN=80) :: routine = 'environ_type_setup'
        !
        !--------------------------------------------------------------------------------
        !
        IF (environ_type == 'input') RETURN
        ! skip set up if read environ keywords from input
        !
        !--------------------------------------------------------------------------------
        ! Set physically meaningful global parameters
        !
        CALL io%writer("* setting up "//TRIM(environ_type)//" environment")
        !
        SELECT CASE (environ_type)
            !
        CASE ('vacuum') ! vacuum case is straightforward, all flags are off
            env_static_permittivity = 1.D0
            env_optical_permittivity = 1.D0
            env_surface_tension = 0.D0
            env_pressure = 0.D0
            !
            RETURN
            !
        CASE ('water', 'water-cation', 'water-anion') ! water experimental permittivities
            env_static_permittivity = 78.3D0
            env_optical_permittivity = 1.776D0
            !
        CASE DEFAULT
            CALL io%error(routine, "Unrecognized value for environ_type", 1)
            !
        END SELECT
        !
        !--------------------------------------------------------------------------------
        ! Depending on the boundary mode, set fitted parameters
        !
        IF (solvent_mode == 'electronic' .OR. solvent_mode == 'full') THEN
            !
            !----------------------------------------------------------------------------
            ! Self-consistent continuum solvation (SCCS)
            !
            SELECT CASE (environ_type)
                !
            CASE ('water') ! SCCS for neutrals
                env_surface_tension = 47.9D0
                env_pressure = -0.36D0
                rhomax = 0.005
                rhomin = 0.0001
                !
            CASE ('water-cation') ! SCCS for cations
                env_surface_tension = 5.D0
                env_pressure = 0.125D0
                rhomax = 0.0035
                rhomin = 0.0002
                !
            CASE ('water-anion') ! SCCS for cations
                env_surface_tension = 0.D0
                env_pressure = 0.45D0
                rhomax = 0.0155
                rhomin = 0.0024
                !
            END SELECT
            !
        ELSE IF (solvent_mode == 'ionic') THEN
            !
            !----------------------------------------------------------------------------
            ! Soft-sphere continuum solvation
            !
            radius_mode = 'uff'
            softness = 0.5D0
            env_surface_tension = 50.D0 ! NOTE THAT WE ARE USING THE
            env_pressure = -0.35D0 ! SET FOR CLUSTERS, AS IN SCCS
            !
            SELECT CASE (environ_type)
                !
            CASE ('water') ! SS for neutrals
                alpha = 1.12D0
                !
            CASE ('water-cation') ! SS for cations
                alpha = 1.1D0
                !
            CASE ('water-anion') ! SS for anions
                alpha = 0.98D0
                !
            END SELECT
            !
        END IF
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE environ_type_setup
    !------------------------------------------------------------------------------------
    !>
    !! Set problem according to the ENVIRON and ELECTROSTATIC namelists
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE electrostatics_setup()
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CHARACTER(LEN=80) :: routine = 'electrostatics_setup'
        !
        !--------------------------------------------------------------------------------
        ! Electrolyte checks
        !
        SELECT CASE (pbc_correction)
            !
        CASE ('gcs', 'ms-gcs')
            !
            IF (electrolyte_distance == 0.0_DP) &
                CALL io%error(routine, &
                              "electrolyte_distance must be set (greater than zero) for gcs correction", 1)
            !
            IF (TRIM(electrolyte_mode) /= 'system') THEN
                electrolyte_mode = 'system'
                !
                CALL io%writer("* setting electrolyte mode to "// &
                               TRIM(electrolyte_mode)//" (required for gcs correction)")
                !
            END IF
            !
        END SELECT
        !
        IF (env_electrolyte_ntyp > 0) THEN
            !
            SELECT CASE (pbc_correction)
                !
            CASE ('gcs', 'ms-gcs')
                !
            CASE DEFAULT
                !
                IF (electrolyte_linearized) THEN
                    !
                    IF (cionmax > 0.D0 .OR. rion > 0.D0) THEN
                        problem = 'linmodpb'
                    ELSE IF (problem == 'none') THEN
                        problem = 'linpb'
                    END IF
                    !
                    IF (solver == 'none') solver = 'cg'
                    !
                ELSE
                    !
                    IF (cionmax > 0.D0 .OR. rion > 0.D0) THEN
                        problem = 'modpb'
                    ELSE IF (problem == 'none') THEN
                        problem = 'pb'
                    END IF
                    !
                    IF (solver == 'none') solver = 'newton'
                    !
                    IF (inner_solver == 'none') inner_solver = 'cg'
                    !
                END IF
                !
            END SELECT
            !
        END IF
        !
        IF (TRIM(pbc_correction) == 'gcs' .OR. env_electrolyte_ntyp > 0) THEN
            !
            SELECT CASE (TRIM(electrolyte_mode))
                !
            CASE ('electronic', 'full', 'system')
                !
                SELECT CASE (TRIM(electrolyte_deriv_method))
                    !
                CASE ('default')
                    electrolyte_deriv_method = 'chain'
                    !
                    CALL io%writer( &
                        "* setting electrolyte-boundary derivatives method to SCCS default")
                    !
                CASE ('highmem', 'lowmem')
                    !
                    CALL io%error(routine, &
                                  "Only 'fft' or 'chain' are allowed with electronic interfaces", 1)
                    !
                END SELECT
                !
            CASE ('ionic')
                !
                SELECT CASE (TRIM(electrolyte_deriv_method))
                    !
                CASE ('default')
                    electrolyte_deriv_method = 'lowmem'
                    !
                    CALL io%writer( &
                        "* setting electrolyte-boundary derivatives method to SSCS default")
                    !
                CASE ('chain')
                    !
                    CALL io%error(routine, &
                                  "Only 'highmem', 'lowmem', and 'fft' are allowed with ionic interfaces", 1)
                    !
                END SELECT
                !
            END SELECT
            !
        END IF
        !
        !--------------------------------------------------------------------------------
        ! Dielectric checks
        !
        IF (env_static_permittivity > 1.D0 .OR. env_dielectric_regions > 0) THEN
            !
            IF (problem == 'none') problem = 'generalized'
            !
            SELECT CASE (pbc_correction)
                !
            CASE ('gcs', 'ms-gcs')
                !
                IF (solver /= 'fixed-point') THEN
                    solver = 'fixed-point'
                    !
                    CALL io%writer("* setting solver to fixed-point (required for gcs correction)")
                    !
                END IF
                !
            CASE DEFAULT
                IF (solver == 'none') solver = 'cg'
                !
            END SELECT
            !
        ELSE
            !
            IF (problem == 'none') problem = 'poisson'
            !
            IF (solver == 'none') solver = 'direct'
            !
        END IF
        !
        !--------------------------------------------------------------------------------
        ! Ensure correct auxiliary is used for fixed-point
        !
        IF (solver == 'fixed-point' .AND. auxiliary == 'none') auxiliary = 'full'
        !
        !--------------------------------------------------------------------------------
        ! Set inner problem
        !
        IF (inner_solver /= 'none') THEN
            !
            SELECT CASE (solver)
                !
            CASE ('fixed-point')
                !
                IF (auxiliary == 'ioncc') inner_problem = 'generalized'
                !
            CASE ('newton')
                inner_problem = 'linpb'
                !
            END SELECT
            !
        END IF
        !
        !--------------------------------------------------------------------------------
        ! Validate use of inner solver
        !
        IF (.NOT. (problem == 'pb' .OR. &
                   problem == 'modpb' .OR. &
                   problem == 'generalized') &
            .AND. (inner_solver /= 'none')) &
            CALL io%error(routine, "Only pb or modpb problems allow inner solver", 1)
        !
        !--------------------------------------------------------------------------------
        ! Validate problem/solver combination and check for PBC correction if needed
        !
        CALL validate_electrostatic_input(problem, solver)
        !
        IF (inner_solver /= 'none') &
            CALL validate_electrostatic_input(inner_problem, inner_solver)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE electrostatics_setup
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE validate_electrostatic_input(problem_in, solver_in)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CHARACTER(LEN=*), INTENT(IN) :: problem_in
        CHARACTER(LEN=*), INTENT(IN) :: solver_in
        !
        CHARACTER(LEN=80) :: routine = 'validate_electrostatic_input'
        !
        !--------------------------------------------------------------------------------
        !
        SELECT CASE (problem_in)
            !
        CASE ('poisson')
            !
        CASE ('generalized') ! generalized Poisson-Boltzmann
            !
            IF (solver_in == 'direct') &
                CALL io%error(routine, "Cannot use a direct solver for the Generalized Poisson eq.", 1)
            !
        CASE ('linpb', 'linmodpb') ! linearized Poisson-Boltzmann
            !
            SELECT CASE (solver_in)
                !
            CASE ('cg', 'sd')
                !
            CASE DEFAULT
                CALL io%error(routine, "Only gradient-based solver for the linearized Poisson-Boltzmann eq.", 1)
                !
            END SELECT
            !
            IF (pbc_correction /= 'parabolic') &
                CALL io%error(routine, "Linearized-PB problem requires parabolic pbc correction", 1)
            !
        CASE ('pb', 'modpb') ! Poisson-Boltzmann
            !
            SELECT CASE (solver_in)
                !
            CASE ('direct', 'cg', 'sd')
                CALL io%error(routine, "No direct or gradient-based solver for the full Poisson-Boltzmann eq.", 1)
                !
            END SELECT
            !
            IF (pbc_correction /= 'parabolic') &
                CALL io%error(routine, "Full-PB problem requires parabolic pbc correction", 1)
            !
        CASE DEFAULT
            CALL io%error(routine, "Unexpected keyword for electrostatic problem", 1)
            !
        END SELECT
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE validate_electrostatic_input
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !                                       CARDS
    !
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !>
    !! Description of the allowed input CARDS
    !!
    !! EXTERNAL_CHARGES (unit_option)
    !!
    !!   set external fixed charge densities and their shape
    !!
    !! Syntax:
    !!
    !!    EXTERNAL_CHARGES (unit_option)
    !!      charge(1)  x(1) y(1) z(1)  spread(1) dim(1)  axis(1)
    !!       ...       ...        ...      ...        ...
    !!      charge(n)  x(n) y(n) z(n)  spread(n) dim(n)  axis(n)
    !!
    !! Example:
    !!
    !! EXTERNAL_CHARGES (bohr)
    !!  1.0  0.0  0.0  0.0  [0.5  2  1]
    !! -1.0  0.0  0.0  5.0  [0.5  2  1]
    !!
    !! Where:
    !!
    !!   unit_option == bohr       positions are given in Bohr (DEFAULT)
    !!   unit_option == angstrom   positions are given in Angstrom
    !!
    !!      charge(i) ( real )       total charge of the density
    !!      x/y/z(i)  ( real )       cartesian position of the density
    !!      spread(i) ( real )       gaussian spread of the density (in bohr, optional, default=0.5)
    !!      dim(i)    ( integer )    0/1/2 point/line/plane of charge (optional, default=0)
    !!      axis(i)   ( integer )    1/2/3 for x/y/z direction of line/plane (optional, default=3)
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE card_external_charges(unit, input_line)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: unit
        !
        CHARACTER(LEN=*), INTENT(OUT) :: input_line
        !
        INTEGER :: i, j
        INTEGER :: nfield
        LOGICAL :: tend
        CHARACTER(LEN=256) :: field_str
        !
        CHARACTER(LEN=80) :: routine = 'card_external_charges'
        !
        !--------------------------------------------------------------------------------
        ! Validate input
        !
        IF (taextchg) CALL io%error(routine, "Two occurrences", 2)
        !
        IF (env_external_charges > local_nsx) &
            CALL io%error(routine, "nsx out of range", env_external_charges)
        !
        CALL allocate_input_extcharge(env_external_charges)
        !
        IF (env_is_substring("BOHR", input_line)) THEN
            extcharge_units = 'bohr'
        ELSE IF (env_is_substring("ANGSTROM", input_line)) THEN
            extcharge_units = 'angstrom'
        ELSE
            !
            IF (TRIM(ADJUSTL(input_line)) /= 'EXTERNAL_CHARGES') &
                CALL io%error(routine, &
                              "Invalid units for EXTERNAL_CHARGES: "//input_line, 1)
            !
            extcharge_units = 'bohr'
            !
            CALL io%writer("* setting external charge units to "//TRIM(extcharge_units))
            !
        END IF
        !
        !--------------------------------------------------------------------------------
        ! Parse card input
        !
        DO i = 1, env_external_charges
            !
            CALL env_read_line(unit, input_line, end_of_file=tend)
            !
            IF (tend) CALL io%error(routine, "End of file reading external charges", i)
            !
            CALL env_field_count(nfield, input_line)
            !
            !----------------------------------------------------------------------------
            ! Read field 1 (total charge of the external density)
            !
            CALL env_get_field(1, field_str, input_line)
            !
            READ (field_str, *) extcharge_charge(i)
            !
            !----------------------------------------------------------------------------
            ! Read fields 2-4 (x-y-z position of external density)
            !
            CALL env_get_field(2, field_str, input_line)
            !
            READ (field_str, *) extcharge_pos(1, i)
            !
            CALL env_get_field(3, field_str, input_line)
            !
            READ (field_str, *) extcharge_pos(2, i)
            !
            CALL env_get_field(4, field_str, input_line)
            !
            READ (field_str, *) extcharge_pos(3, i)
            !
            !----------------------------------------------------------------------------
            ! Optionally read field 5 (spread of the density)
            !
            IF (nfield >= 5) THEN
                !
                CALL env_get_field(5, field_str, input_line)
                !
                READ (field_str, *) extcharge_spread(i)
                !
                IF (extcharge_spread(i) < 0.D0) &
                    CALL io%error(routine, "Spread must be positive", i)
                !
            ELSE
                CALL io%writer("* setting external charge spread to 0.5 (a.u.)")
            END IF
            !
            !----------------------------------------------------------------------------
            ! Optionally read field 6 and 7 (dimensionality and direction)
            !
            IF (nfield >= 6) THEN
                !
                CALL env_get_field(6, field_str, input_line)
                !
                READ (field_str, *) extcharge_dim(i)
                !
                IF (extcharge_dim(i) < 0 .OR. extcharge_dim(i) > 2) &
                    CALL io%error(routine, "Wrong excharge dimension", i)
                !
                IF (extcharge_dim(i) > 0) THEN
                    !
                    IF (nfield == 6) THEN
                        CALL io%writer("* setting external charge axis to 3 (z-axis)")
                    ELSE
                        !
                        CALL env_get_field(7, field_str, input_line)
                        !
                        READ (field_str, *) extcharge_axis(i)
                        !
                        IF (extcharge_axis(i) < 0 .OR. extcharge_axis(i) > 3) &
                            CALL io%error(routine, "Wrong excharge axis", i)
                        !
                    END IF
                    !
                END IF
                !
            ELSE
                CALL io%writer("* setting external charge dimensions to 0D")
            END IF
            !
        END DO
        !
        !--------------------------------------------------------------------------------
        ! Convert to atomic units
        !
        taextchg = .TRUE.
        !
        DO i = 1, env_external_charges
            !
            DO j = 1, 3
                CALL convert_length(extcharge_units, extcharge_pos(j, i))
            END DO
            !
            CALL convert_length(extcharge_units, extcharge_spread(i))
            !
        END DO
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE card_external_charges
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE allocate_input_extcharge(external_charges)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: external_charges
        !
        !--------------------------------------------------------------------------------
        !
        IF (ALLOCATED(extcharge_charge)) DEALLOCATE (extcharge_charge)
        !
        IF (ALLOCATED(extcharge_pos)) DEALLOCATE (extcharge_pos)
        !
        IF (ALLOCATED(extcharge_spread)) DEALLOCATE (extcharge_spread)
        !
        IF (ALLOCATED(extcharge_dim)) DEALLOCATE (extcharge_dim)
        !
        IF (ALLOCATED(extcharge_axis)) DEALLOCATE (extcharge_axis)
        !
        ALLOCATE (extcharge_charge(external_charges))
        ALLOCATE (extcharge_pos(3, external_charges))
        ALLOCATE (extcharge_spread(external_charges))
        ALLOCATE (extcharge_dim(external_charges))
        ALLOCATE (extcharge_axis(external_charges))
        !
        extcharge_charge = 0.0_DP
        extcharge_pos = 0.0_DP
        extcharge_spread = 0.5_DP
        extcharge_dim = 0
        extcharge_axis = 3
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE allocate_input_extcharge
    !------------------------------------------------------------------------------------
    !>
    !! Description of the allowed input CARDS
    !!
    !! DIELECTRIC_REGIONS (unit_option)
    !!
    !!   set fixed dielectric regions and their shape
    !!
    !! Syntax:
    !!
    !!    DIELECTRIC_REGIONS (unit_option)
    !!      epsilon0(1) epsilonopt(1) x(1) y(1) z(1)  width(1) spread(1) dim(1)  axis(1)
    !!       ...       ...        ...      ...        ...
    !!      epsilon0(n) epsilonopt(n) x(n) y(n) z(n)  width(n) spread(n) dim(n)  axis(n)
    !!
    !! Example:
    !!
    !! DIELECTRIC_REGIONS (bohr)
    !!  80.0  2.0   0.0  0.0  10.0   5.0  1.0  2  3
    !!
    !! Where:
    !!
    !!   unit_option == bohr       positions are given in Bohr (DEFAULT)
    !!   unit_option == angstrom   positions are given in Angstrom
    !!
    !!      epsilon0(i)   ( real )    static permittivity inside the region
    !!      epsilonopt(i) ( real )    optical permittivity inside the region
    !!      x/y/z(i)      ( real )    cartesian center of the region
    !!      width(i)      ( real )    size of the region (in bohr)
    !!      spread(i)     ( real )    spread of the interface (in bohr, optional)
    !!      dim(i)     ( integer )    0/1/2 point/line/plane region (optional)
    !!      axis(i)    ( integer )    1/2/3 for x/y/z direction of line/plane (optional)
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE card_dielectric_regions(unit, input_line)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: unit
        !
        CHARACTER(LEN=*), INTENT(OUT) :: input_line
        !
        INTEGER :: i, j
        INTEGER :: nfield
        LOGICAL :: tend
        CHARACTER(LEN=256) :: field_str
        !
        CHARACTER(LEN=80) :: routine = 'card_dielectric_regions'
        !
        !--------------------------------------------------------------------------------
        !
        IF (taepsreg) CALL io%error(routine, "Two occurrences", 2)
        !
        IF (env_dielectric_regions > local_nsx) &
            CALL io%error(routine, "nsx out of range", env_dielectric_regions)
        !
        CALL allocate_input_epsregion(env_dielectric_regions)
        !
        IF (env_is_substring("BOHR", input_line)) THEN
            epsregion_units = 'bohr'
        ELSE IF (env_is_substring("ANGSTROM", input_line)) THEN
            epsregion_units = 'angstrom'
        ELSE
            !
            IF (TRIM(ADJUSTL(input_line)) /= 'DIELECTRIC_REGIONS') &
                CALL io%error(routine, &
                              "Invalid units for DIELECTRIC_REGIONS: "//input_line, 1)
            !
            epsregion_units = 'bohr'
            !
            CALL io%writer("* setting dielectric region units to "// &
                           TRIM(epsregion_units))
            !
        END IF
        !
        !--------------------------------------------------------------------------------
        ! Parse card input
        !
        DO i = 1, env_dielectric_regions
            !
            CALL env_read_line(unit, input_line, end_of_file=tend)
            !
            IF (tend) &
                CALL io%error(routine, "End of file reading dielectric regions", i)
            !
            CALL env_field_count(nfield, input_line)
            !
            !----------------------------------------------------------------------------
            ! Read field 1-2 (static and optical permettivity inside dielectric region)
            !
            CALL env_get_field(1, field_str, input_line)
            !
            READ (field_str, *) epsregion_eps(1, i)
            !
            IF (epsregion_eps(1, i) < 1.D0) &
                CALL io%error(routine, "Static permittivity must be > 1", i)
            !
            CALL env_get_field(2, field_str, input_line)
            !
            READ (field_str, *) epsregion_eps(2, i)
            !
            IF (epsregion_eps(2, i) < 1.D0) &
                CALL io%error(routine, "Optical permittivity must be > 1", i)
            !
            !----------------------------------------------------------------------------
            ! Read fields 3-5 (x-y-z position of dielectric region)
            !
            CALL env_get_field(3, field_str, input_line)
            !
            READ (field_str, *) epsregion_pos(1, i)
            !
            CALL env_get_field(4, field_str, input_line)
            !
            READ (field_str, *) epsregion_pos(2, i)
            !
            CALL env_get_field(5, field_str, input_line)
            !
            READ (field_str, *) epsregion_pos(3, i)
            !
            !----------------------------------------------------------------------------
            ! Read field 6 (size/width of the dielectric region)
            !
            CALL env_get_field(6, field_str, input_line)
            !
            READ (field_str, *) epsregion_width(i)
            !
            IF (epsregion_width(i) < 0.D0) &
                CALL io%error(routine, "Width must be positive", i)
            !
            !----------------------------------------------------------------------------
            ! Optionally read field 7 (spread of interface of the dielectric region)
            !
            IF (nfield >= 7) THEN
                !
                CALL env_get_field(7, field_str, input_line)
                !
                READ (field_str, *) epsregion_spread(i)
                !
                IF (epsregion_spread(i) < 0.D0) &
                    CALL io%error(routine, "Spread must be positive", i)
                !
            ELSE
                CALL io%writer("* setting dielectric region spread to 0.5 (a.u.)")
            END IF
            !
            !----------------------------------------------------------------------------
            ! Optionally read field 8 and 9 (dimensionality and direction)
            !
            IF (nfield >= 8) THEN
                !
                CALL env_get_field(8, field_str, input_line)
                !
                READ (field_str, *) epsregion_dim(i)
                !
                IF (epsregion_dim(i) < 0 .OR. epsregion_dim(i) > 2) &
                    CALL io%error(routine, "Wrong epsregion dimension", i)
                !
                IF (epsregion_dim(i) > 0) THEN
                    !
                    IF (nfield == 8) THEN
                        CALL io%writer("* setting dielectric region axis to 3 (z-axis)")
                    ELSE
                        !
                        CALL env_get_field(9, field_str, input_line)
                        !
                        READ (field_str, *) epsregion_axis(i)
                        !
                        IF (epsregion_axis(i) < 1 .OR. epsregion_axis(i) > 3) &
                            CALL io%error(routine, "Wrong epsregion axis", i)
                        !
                    END IF
                    !
                END IF
                !
            ELSE
                CALL io%writer("* setting dielectric region dimensions to 0D")
            END IF
            !
        END DO
        !
        !--------------------------------------------------------------------------------
        ! Convert to atomic units
        !
        taepsreg = .TRUE.
        !
        DO i = 1, env_dielectric_regions
            !
            DO j = 1, 3
                CALL convert_length(epsregion_units, epsregion_pos(j, i))
            END DO
            !
            CALL convert_length(epsregion_units, epsregion_width(i))
            !
            CALL convert_length(epsregion_units, epsregion_spread(i))
            !
        END DO
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE card_dielectric_regions
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE allocate_input_epsregion(dielectric_regions)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: dielectric_regions
        !
        !--------------------------------------------------------------------------------
        !
        IF (ALLOCATED(epsregion_eps)) DEALLOCATE (epsregion_eps)
        !
        IF (ALLOCATED(epsregion_pos)) DEALLOCATE (epsregion_pos)
        !
        IF (ALLOCATED(epsregion_width)) DEALLOCATE (epsregion_width)
        !
        IF (ALLOCATED(epsregion_spread)) DEALLOCATE (epsregion_spread)
        !
        IF (ALLOCATED(epsregion_dim)) DEALLOCATE (epsregion_dim)
        !
        IF (ALLOCATED(epsregion_axis)) DEALLOCATE (epsregion_axis)
        !
        ALLOCATE (epsregion_eps(2, dielectric_regions))
        ALLOCATE (epsregion_pos(3, dielectric_regions))
        ALLOCATE (epsregion_spread(dielectric_regions))
        ALLOCATE (epsregion_width(dielectric_regions))
        ALLOCATE (epsregion_dim(dielectric_regions))
        ALLOCATE (epsregion_axis(dielectric_regions))
        !
        epsregion_eps = 1.0_DP
        epsregion_pos = 0.0_DP
        epsregion_spread = 0.5_DP
        epsregion_width = 0.0_DP
        epsregion_dim = 0
        epsregion_axis = 3
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE allocate_input_epsregion
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !                                     UTILITIES
    !
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !>
    !! Convert input length to atomic units
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE convert_length(length_format, length)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CHARACTER(LEN=*), INTENT(IN) :: length_format
        !
        REAL(DP), INTENT(INOUT) :: length
        !
        CHARACTER(LEN=80) :: routine = 'convert_length'
        !
        !--------------------------------------------------------------------------------
        !
        SELECT CASE (length_format)
            !
        CASE ('bohr') ! input length are in a.u., do nothing
            !
        CASE ('angstrom')
            length = length / BOHR_RADIUS_ANGS ! length in A: convert to a.u.
            !
        CASE DEFAULT
            CALL io%error(routine, TRIM(length_format)//" units not implemented", 1)
            !
        END SELECT
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE convert_length
    !------------------------------------------------------------------------------------
    !>
    !! WE MAY WANT TO ADD A SECOND COMM ON IMAGES #TODO may be required for NEB
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_read_line(unit, line, nfield, field, end_of_file, error)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: unit
        CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: field
        INTEGER, OPTIONAL, INTENT(IN) :: nfield
        !
        CHARACTER(LEN=*), INTENT(OUT) :: line
        LOGICAL, OPTIONAL, INTENT(OUT) :: end_of_file, error
        !
        INTEGER :: ios
        LOGICAL :: tend, terr
        !
        CHARACTER(LEN=80) :: routine = 'env_read_line'
        !
        !--------------------------------------------------------------------------------
        !
        IF (LEN(line) < 256) &
            CALL io%error(routine, "Input line too short", MAX(LEN(line), 1))
        !
        !--------------------------------------------------------------------------------
        !
        tend = .FALSE.
        terr = .FALSE.
        !
        IF (io%lnode) THEN
            ios = 0
            !
            DO WHILE (ios == 0)
                !
                READ (unit, fmt='(A256)', iostat=ios) line
                !
                IF (ios /= 0) EXIT
                !
                line = TRIM(ADJUSTL(line))
                !
                !------------------------------------------------------------------------
                ! Skip comment lines
                !
                IF (line == ' ' .OR. (line(1:1) == '#' .OR. &
                                      line(1:1) == '!' .OR. &
                                      line(1:1) == '/')) &
                    CYCLE
                !
                !------------------------------------------------------------------------
                ! Consume unnecessary namelists
                !
                IF (line(1:1) == '&') THEN
                    !
                    CALL io%writer("* ignoring "//TRIM(line)//" namelist")
                    !
                    DO WHILE (line(1:1) /= '/')
                        !
                        READ (unit, fmt='(A256)', iostat=ios) line
                        !
                        IF (ios /= 0) EXIT
                        !
                        line = TRIM(ADJUSTL(line))
                    END DO
                    !
                    CYCLE
                END IF
                !
                IF (ios == 0) EXIT
                !
            END DO
            !
            IF (ios < 0) THEN
                tend = .TRUE.
            ELSE IF (ios > 0) THEN
                terr = .TRUE.
            END IF
            !
        END IF
        !
        CALL env_mp_bcast(tend, io%node, io%comm)
        !
        CALL env_mp_bcast(terr, io%node, io%comm)
        !
        CALL env_mp_bcast(line, io%node, io%comm)
        !
        IF (PRESENT(end_of_file)) THEN
            end_of_file = tend
        ELSE IF (tend) THEN
            CALL io%writer("End of file")
        END IF
        !
        IF (PRESENT(error)) THEN
            error = terr
        ELSE IF (terr) THEN
            CALL io%writer("Read error")
        END IF
        !
        IF (PRESENT(field) .AND. .NOT. (tend .OR. terr)) &
            CALL env_field_compare(line, nfield, field)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_read_line
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_field_count(num, line, car)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CHARACTER(LEN=*), INTENT(IN) :: line
        CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: car
        !
        INTEGER, INTENT(OUT) :: num
        !
        CHARACTER(LEN=1) :: sep1, sep2
        INTEGER :: j
        !
        !--------------------------------------------------------------------------------
        !
        num = 0
        !
        IF (.NOT. PRESENT(car)) THEN
            sep1 = CHAR(32) ! blank character
            sep2 = CHAR(9) ! tab character
            !
            DO j = 2, MAX(LEN(line), 256)
                !
                IF (line(j:j) == '!' .OR. line(j:j) == CHAR(0)) THEN
                    !
                    IF (line(j - 1:j - 1) /= sep1 .AND. line(j - 1:j - 1) /= sep2) &
                        num = num + 1
                    !
                    EXIT
                    !
                END IF
                !
                IF ((line(j:j) == sep1 .OR. line(j:j) == sep2) .AND. &
                    (line(j - 1:j - 1) /= sep1 .AND. line(j - 1:j - 1) /= sep2)) &
                    num = num + 1
                !
            END DO
            !
        ELSE
            !
            sep1 = car
            !
            DO j = 2, MAX(LEN(line), 256)
                !
                IF (line(j:j) == '!' .OR. &
                    line(j:j) == CHAR(0) .OR. line(j:j) == CHAR(32)) THEN
                    !
                    IF (line(j - 1:j - 1) /= sep1) num = num + 1
                    !
                    EXIT
                    !
                END IF
                !
                IF (line(j:j) == sep1 .AND. line(j - 1:j - 1) /= sep1) num = num + 1
                !
            END DO
            !
        END IF
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_field_count
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_field_compare(str, nf, var)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: nf
        CHARACTER(LEN=*), INTENT(IN) :: str, var
        !
        INTEGER :: nc
        !
        CHARACTER(LEN=80) :: routine = 'env_field_compare'
        !
        !--------------------------------------------------------------------------------
        !
        CALL env_field_count(nc, str)
        !
        IF (nc < nf) CALL io%error(routine, "Wrong number of fields: "//TRIM(var), 1)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_field_compare
    !------------------------------------------------------------------------------------
    !>
    !! Extract whitespace-separated nth block from string
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_get_field(n, field, str, sep)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: n
        CHARACTER(LEN=*), INTENT(IN) :: str
        CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: sep
        !
        CHARACTER(LEN=*), INTENT(OUT) :: field
        !
        INTEGER :: i, j, z ! block start and end
        INTEGER :: k ! block counter
        CHARACTER(LEN=1) :: sep1, sep2
        !
        !--------------------------------------------------------------------------------
        !
        IF (PRESENT(sep)) THEN
            sep1 = sep
            sep2 = sep ! redundant, but easy
        ELSE
            sep1 = CHAR(32) ! blank character
            sep2 = CHAR(9) ! tab char
        END IF
        !
        k = 1 ! counter for the required block
        !
        DO i = 1, LEN(str)
            !
            z = MAX(i - 1, 1) ! look for the beginning of the required block
            !
            IF (k == n) EXIT
            !
            IF ((str(i:i) == sep1 .OR. str(i:i) == sep2) .AND. &
                (str(z:z) /= sep1 .AND. str(z:z) /= sep2)) &
                k = k + 1
            !
        END DO
        !
        DO j = i, LEN(str)
            !
            z = MAX(j - 1, 1) ! look for the beginning of the next block
            !
            IF ((str(j:j) == sep1 .OR. str(j:j) == sep2) .AND. &
                (str(z:z) /= sep1 .AND. str(z:z) /= sep2)) &
                k = k + 1
            !
            IF (k > n) EXIT
            !
        END DO
        !
        IF (j <= LEN(str)) THEN
            field = TRIM(ADJUSTL(str(i:j - 1)))
            ! if we are here, the reqired block was followed by a separator
            ! and another field. We have to trash one char (a separator)
        ELSE
            field = TRIM(ADJUSTL(str(i:LEN(str))))
            ! if we are here, it was the last block in str. We have to take
            ! all the remaining chars
        END IF
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_get_field
    !------------------------------------------------------------------------------------
    !>
    !! Allocate arrays depending on nsx
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE allocate_registers()
        !--------------------------------------------------------------------------------
        !
        ALLOCATE (atomicspread(local_nsx))
        ALLOCATE (solvationrad(local_nsx))
        ALLOCATE (corespread(local_nsx))
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE allocate_registers
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE environ_input
!----------------------------------------------------------------------------------------
