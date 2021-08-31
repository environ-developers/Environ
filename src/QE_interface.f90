!----------------------------------------------------------------------------------------
!
! Copyright (C) 2021 ENVIRON (www.quantum-environ.org)
!
!----------------------------------------------------------------------------------------
!
!     This file is part of Environ version 2.0
!
!     Environ 2.0 is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation, either version 2 of the License, or
!     (at your option) any later version.
!
!     Environ 2.0 is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more detail, either the file
!     `License' in the root directory of the present distribution, or
!     online at <http://www.gnu.org/licenses/>.
!
!----------------------------------------------------------------------------------------
!
! Authors: Oliviero Andreussi (Department of Physics, UNT)
!          Edan Bainglass     (Department of Physics, UNT)
!
!----------------------------------------------------------------------------------------
!>
!! An interface module for all communications between Environ and the calling program
!!
!----------------------------------------------------------------------------------------
MODULE environ_QE_interface
    !------------------------------------------------------------------------------------
    !
    USE env_char_ops, ONLY: env_uppercase
    !
    USE env_base_io, ONLY: prog, ionode, ionode_id, comm, program_unit, environ_unit, &
                           lstdout, depth, verbose_ => verbose
    !
    USE env_io, ONLY: env_find_free_unit
    !
    USE environ_param, ONLY: DP, BOHR_RADIUS_SI, RYDBERG_SI, RYTOEV, e2
    !
    USE env_global_objects, ONLY: setup, env
    !
    !------------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    PRIVATE
    !
    PUBLIC :: init_environ_io, read_environ_input, init_environ_setup, &
              init_environ_cell, init_environ_base

    PUBLIC :: update_environ_potential, update_environ_cell, update_environ_ions, &
              update_environ_electrons, update_environ_response
    !
    PUBLIC :: calc_environ_potential, calc_environ_energy, calc_environ_force, &
              calc_environ_dpotential, calc_environ_denergy
    !
    PUBLIC :: print_environ_energies, print_environ_potential_shift, &
              print_environ_potential_warning, print_environ_summary, &
              print_environ_clocks, update_output_program_unit
    !
    PUBLIC :: environ_clean, environ_clean_first, environ_clean_second
    !
    PUBLIC :: is_tddfpt, is_environ_restart
    !
    PUBLIC :: get_environ_threshold, set_environ_restart, get_environ_nskip, &
              get_environ_dvtot_of_r, get_environ_verbose
    !
    !------------------------------------------------------------------------------------
    !
    INTERFACE calc_environ_potential
        MODULE PROCEDURE calc_environ_potential_PW, calc_environ_potential_CP
    END INTERFACE calc_environ_potential
    !
    !------------------------------------------------------------------------------------
CONTAINS
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !                               INITIALIZATION METHODS
    !
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !>
    !! Set global I/O constants
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_environ_io(prog_, ionode_, ionode_id_, comm_, program_unit_)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CHARACTER(LEN=*), INTENT(IN) :: prog_
        LOGICAL, INTENT(IN) :: ionode_
        INTEGER, INTENT(IN) :: ionode_id_
        INTEGER, INTENT(IN) :: comm_
        INTEGER, INTENT(IN) :: program_unit_
        !
        !--------------------------------------------------------------------------------
        !
        ionode = ionode_
        ionode_id = ionode_id_
        comm = comm_
        !
        program_unit = program_unit_
        environ_unit = env_find_free_unit()
        !
        prog = env_uppercase(prog_(1:2))
        !
        SELECT CASE (prog)
            !
        CASE ('PW', 'CP')
            lstdout = .TRUE.
            !
        CASE DEFAULT
            lstdout = .FALSE.
            !
        END SELECT
        !
        lstdout = lstdout .AND. ionode
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_environ_io
    !------------------------------------------------------------------------------------
    !>
    !! Initialize Environ setup object from file
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE read_environ_input()
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CHARACTER(LEN=80) :: sub_name = 'read_environ_input'
        !
        !--------------------------------------------------------------------------------
        !
        CALL setup%read_input() ! default filename = environ.in
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE read_environ_input
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_environ_setup()
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CHARACTER(LEN=80) :: sub_name = 'init_environ_setup'
        !
        !--------------------------------------------------------------------------------
        !
        CALL setup%set_flags()
        !
        CALL setup%set_numerical_base()
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_environ_setup
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_environ_cell(alat, at, comm_in, gcutm, use_internal_pbc_corr)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: comm_in
        REAL(DP), INTENT(IN) :: alat, at(3, 3), gcutm
        LOGICAL, INTENT(IN), OPTIONAL :: use_internal_pbc_corr
        !
        REAL(DP), ALLOCATABLE :: at_scaled(:, :)
        REAL(DP) :: gcutm_scaled
        !
        CHARACTER(LEN=80) :: sub_name = 'init_environ_cell'
        !
        !--------------------------------------------------------------------------------
        !
        IF (alat < 1.D-8) CALL env_errore(sub_name, 'Wrong alat', 1)
        !
        IF (alat < 1.0_DP) CALL env_warning('strange lattice parameter')
        !
        !--------------------------------------------------------------------------------
        ! Allocate buffers used by env_mp_sum
        !
        CALL env_allocate_mp_buffers()
        !
        !--------------------------------------------------------------------------------
        !
        ALLOCATE (at_scaled(3, 3))
        at_scaled = at * alat
        !
        gcutm_scaled = gcutm / alat**2
        !
        CALL setup%init_cell(gcutm_scaled, comm_in, at_scaled)
        !
        CALL setup%init_cores(gcutm_scaled, use_internal_pbc_corr)
        !
        DEALLOCATE (at_scaled)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_environ_cell
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_environ_base(nelec, nat, ntyp, atom_label, ityp, zv, e2_in)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: nelec, nat, ntyp
        INTEGER, INTENT(IN) :: ityp(nat)
        CHARACTER(LEN=3), INTENT(IN) :: atom_label(:)
        REAL(DP), INTENT(IN) :: zv(ntyp)
        REAL(DP), INTENT(IN), OPTIONAL :: e2_in
        !
        REAL(DP), ALLOCATABLE :: at_scaled(:, :)
        REAL(DP) :: gcutm_scaled
        !
        CHARACTER(LEN=80) :: sub_name = 'init_environ_base'
        !
        !--------------------------------------------------------------------------------
        ! PW uses Ryderg units (2.D0 * AU)
        ! CP uses Hartree units (e2_in = 1.D0)
        !
        IF (PRESENT(e2_in)) THEN
            e2 = e2_in
        ELSE
            e2 = 2.D0
        END IF
        !
        !--------------------------------------------------------------------------------
        !
        CALL env%init(setup, nelec, nat, ntyp, atom_label, ityp, zv)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_environ_base
    !------------------------------------------------------------------------------------
    !>
    !! Called at every ionic step
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE update_environ_potential(nnr, vltot)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: nnr
        REAL(DP), INTENT(IN) :: vltot(nnr)
        !
        !--------------------------------------------------------------------------------
        !
        CALL env%update_potential(nnr, vltot)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE update_environ_potential
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE update_environ_cell(at, alat)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP), INTENT(IN) :: at(3, 3), alat
        !
        REAL(DP), ALLOCATABLE :: at_scaled(:, :)
        !
        !--------------------------------------------------------------------------------
        !
        ALLOCATE (at_scaled(3, 3))
        at_scaled = at * alat
        !
        CALL setup%update_cell(at_scaled)
        !
        CALL env%update_cell_dependent_quantities()
        !
        CALL setup%end_cell_update()
        !
        DEALLOCATE (at_scaled)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE update_environ_cell
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE update_environ_ions(nat, tau, alat)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: nat
        REAL(DP), INTENT(IN) :: tau(3, nat), alat
        !
        REAL(DP), ALLOCATABLE :: tau_scaled(:, :)
        !
        !--------------------------------------------------------------------------------
        !
        ALLOCATE (tau_scaled(3, nat))
        tau_scaled = tau * alat
        !
        CALL env%update_ions(nat, tau_scaled)
        !
        DEALLOCATE (tau_scaled)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE update_environ_ions
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE update_environ_electrons(nnr, rho, nelec)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: nnr
        REAL(DP), INTENT(IN) :: rho(nnr)
        REAL(DP), INTENT(IN), OPTIONAL :: nelec
        !
        !--------------------------------------------------------------------------------
        !
        CALL env%update_electrons(nnr, rho, nelec)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE update_environ_electrons
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE update_environ_response(nnr, drho)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: nnr
        REAL(DP), INTENT(IN) :: drho(nnr)
        !
        !--------------------------------------------------------------------------------
        !
        CALL env%update_response(nnr, drho)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE update_environ_response
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !                               COMPUTATION METHODS
    !
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_environ_potential_PW(update, nnr, vtot, local_verbose)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        LOGICAL, INTENT(IN) :: update
        INTEGER, INTENT(IN) :: nnr
        INTEGER, INTENT(IN), OPTIONAL :: local_verbose
        !
        REAL(DP), INTENT(OUT) :: vtot(nnr)
        !
        !--------------------------------------------------------------------------------
        !
        vtot = env%vzero%of_r
        !
        CALL env%potential(update, local_verbose)
        !
        vtot = vtot + env%dvtot%of_r
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_environ_potential_PW
    !------------------------------------------------------------------------------------
    !>
    !!
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_environ_potential_CP(update, local_verbose)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        LOGICAL, INTENT(IN) :: update
        INTEGER, INTENT(IN), OPTIONAL :: local_verbose
        !
        !--------------------------------------------------------------------------------
        !
        CALL env%potential(update, local_verbose)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_environ_potential_CP
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_environ_energy(total_energy)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP), INTENT(INOUT) :: total_energy
        !
        !--------------------------------------------------------------------------------
        !
        CALL env%energy(total_energy)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_environ_energy
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_environ_force(nat, force_environ)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: nat
        !
        REAL(DP), INTENT(INOUT) :: force_environ(3, nat)
        !
        !--------------------------------------------------------------------------------
        !
        CALL env%force(nat, force_environ)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_environ_force
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_environ_dpotential(nnr, dv)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: nnr
        !
        REAL(DP), INTENT(INOUT) :: dv(nnr)
        !
        !--------------------------------------------------------------------------------
        !
        CALL env%dpotential(nnr, dv)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_environ_dpotential
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_environ_denergy(total_energy)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP), INTENT(INOUT) :: total_energy
        !
        !--------------------------------------------------------------------------------
        !
        CALL env%denergy(total_energy)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_environ_denergy
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !                                 ACCESS METHODS
    !
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE set_environ_restart(environ_restart)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        LOGICAL, INTENT(IN) :: environ_restart
        !
        !--------------------------------------------------------------------------------
        !
        setup%restart = environ_restart
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE set_environ_restart
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    REAL(DP) FUNCTION get_environ_threshold()
        !--------------------------------------------------------------------------------
        !
        get_environ_threshold = setup%threshold
        !
        !--------------------------------------------------------------------------------
    END FUNCTION get_environ_threshold
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    INTEGER FUNCTION get_environ_nskip()
        !--------------------------------------------------------------------------------
        !
        get_environ_nskip = setup%nskip
        !
        !--------------------------------------------------------------------------------
    END FUNCTION get_environ_nskip
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    FUNCTION get_environ_dvtot_of_r(nnr) RESULT(dvtot_of_r)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: nnr
        !
        REAL(DP) :: dvtot_of_r(nnr)
        !
        !--------------------------------------------------------------------------------
        !
        dvtot_of_r = env%dvtot%of_r
        !
        !--------------------------------------------------------------------------------
    END FUNCTION get_environ_dvtot_of_r
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    INTEGER FUNCTION get_environ_verbose()
        !--------------------------------------------------------------------------------
        !
        get_environ_verbose = verbose_
        !
        !--------------------------------------------------------------------------------
    END FUNCTION get_environ_verbose
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    LOGICAL FUNCTION is_environ_restart()
        !--------------------------------------------------------------------------------
        !
        is_environ_restart = setup%restart
        !
        !--------------------------------------------------------------------------------
    END FUNCTION is_environ_restart
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    LOGICAL FUNCTION is_tddfpt()
        !--------------------------------------------------------------------------------
        !
        is_tddfpt = setup%ltddfpt
        !
        !--------------------------------------------------------------------------------
    END FUNCTION is_tddfpt
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !                               CLEANING METHODS
    !
    !------------------------------------------------------------------------------------
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
        CALL environ_clean_first(lflag)
        !
        CALL environ_clean_second(lflag)
        !
        CALL env_deallocate_mp_buffers()
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
    SUBROUTINE environ_clean_first(lflag)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        LOGICAL, INTENT(IN) :: lflag
        !
        !--------------------------------------------------------------------------------
        ! Deallocate environment variables
        !
        IF (ASSOCIATED(env%vzero%cell)) CALL env%vzero%destroy()
        !
        IF (ASSOCIATED(env%dvtot%cell)) CALL env%dvtot%destroy()
        !
        !--------------------------------------------------------------------------------
        ! base_environ variables
        !
        IF (setup%lelectrostatic .AND. ASSOCIATED(env%vreference%cell)) &
            CALL env%vreference%destroy()
        !
        IF (setup%lsoftcavity .AND. ASSOCIATED(env%vsoftcavity%cell)) &
            CALL env%vsoftcavity%destroy()
        !
        IF (setup%lconfine .AND. ASSOCIATED(env%vconfine%cell)) &
            CALL env%vconfine%destroy()
        !
        !--------------------------------------------------------------------------------
        ! Destroy derived types which were allocated in input
        !
        IF (setup%lelectrostatic .OR. setup%lconfine) THEN
            !
            CALL env%system_charges%destroy(lflag)
            !
            CALL env%environment_charges%destroy(lflag)
            !
        END IF
        !
        IF (setup%lexternals) CALL env%externals%destroy(lflag)
        !
        IF (setup%lstatic) CALL env%static%destroy(lflag)
        !
        IF (setup%lelectrolyte) CALL env%electrolyte%destroy(lflag)
        !
        IF (setup%lsemiconductor) CALL env%semiconductor%destroy(lflag)
        !
        CALL env%system_electrons%destroy(lflag)
        !
        CALL env%system_ions%destroy(lflag)
        !
        CALL env%system_system%destroy(lflag)
        !
        CALL env%environment_electrons%destroy(lflag)
        !
        CALL env%environment_ions%destroy(lflag)
        !
        CALL env%environment_system%destroy(lflag)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE environ_clean_first
    !------------------------------------------------------------------------------------
    !>
    !! Clean up all the Environ-related allocated variables and call clean up
    !! subroutines of specific Environ modules. These are quantities that may
    !! be needed by TDDFPT, thus may need to be cleaned later
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE environ_clean_second(lflag)
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
        IF (setup%lelectrostatic) THEN
            !
            IF (ASSOCIATED(env%velectrostatic%cell)) CALL env%velectrostatic%destroy()
            !
            CALL electrostatic_clean(lflag)
            !
        END IF
        !
        !--------------------------------------------------------------------------------
        ! Destroy derived types which were allocated in input
        !
        IF (setup%loptical) THEN
            !
            CALL env%environment_response_charges%destroy(lflag)
            !
            CALL env%environment_response_electrons%destroy(lflag)
            !
            CALL env%system_response_charges%destroy(lflag)
            !
            CALL env%system_response_electrons%destroy(lflag)
            !
            CALL env%optical%destroy(lflag)
            !
        END IF
        !
        IF (setup%lsolvent) CALL env%solvent%destroy(lflag)
        !
        IF (setup%lboundary) CALL setup%derivatives%destroy(lflag)
        !
        IF (setup%ldoublecell) THEN
            !
            CALL setup%mapping%destroy(lflag)
            !
            CALL setup%environment_cell%destroy(lflag)
            !
        END IF
        !
        CALL setup%system_cell%destroy(lflag)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE environ_clean_second
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE electrostatic_clean(lflag)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        LOGICAL, INTENT(IN) :: lflag
        !
        !--------------------------------------------------------------------------------
        !
        CALL setup%outer%destroy(lflag)
        !
        CALL setup%reference%destroy(lflag)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE electrostatic_clean
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !                                   OUTPUT METHODS
    !
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !>
    !! Write out the different Environ contributions to the energy.
    !! Called by electrons.f90
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE print_environ_energies()
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CHARACTER(LEN=80) :: sub_name = 'print_environ_energies'
        !
        !--------------------------------------------------------------------------------
        !
        IF (ionode) THEN
            !
            SELECT CASE (prog)
                !
            CASE ('PW')
                !
                IF (setup%lelectrostatic) WRITE (program_unit, 1000) env%eelectrostatic
                !
                IF (setup%lsurface) WRITE (program_unit, 1001) env%esurface
                !
                IF (setup%lvolume) WRITE (program_unit, 1002) env%evolume
                !
                IF (setup%lconfine) WRITE (program_unit, 1003) env%econfine
                !
                IF (setup%lelectrolyte) WRITE (program_unit, 1004) env%eelectrolyte
                !
                WRITE (program_unit, 1005) env%deenviron
                !
            CASE ('CP')
                !
                IF (setup%lelectrostatic) WRITE (program_unit, 1006) env%eelectrostatic
                !
                IF (setup%lsurface) WRITE (program_unit, 1007) env%esurface
                !
                IF (setup%lvolume) WRITE (program_unit, 1008) env%evolume
                !
                IF (setup%lconfine) WRITE (program_unit, 1009) env%econfine
                !
                IF (setup%lelectrolyte) WRITE (program_unit, 1010) env%eelectrolyte
                !
            CASE DEFAULT
                CALL env_errore(sub_name, 'Wrong program calling Environ', 1)
                !
            END SELECT
            !
        END IF
        !
        !--------------------------------------------------------------------------------
        !
1000    FORMAT('     electrostatic embedding   =', F17.8, ' Ry')
1001    FORMAT('     cavitation energy         =', F17.8, ' Ry')
1002    FORMAT('     PV energy                 =', F17.8, ' Ry')
1003    FORMAT('     confinement energy        =', F17.8, ' Ry')
1004    FORMAT('     electrolyte free energy   =', F17.8, ' Ry')
1005    FORMAT('     correction to one-el term =', F17.8, ' Ry')
        !
1006    FORMAT('     electrostatic embedding = ', F14.5, ' Hartree a.u.')
1007    FORMAT('           cavitation energy = ', F14.5, ' Hartree a.u.')
1008    FORMAT('                   PV energy = ', F14.5, ' Hartree a.u.')
1009    FORMAT('     electrolyte free energy = ', F14.5, ' Hartree a.u.')
1010    FORMAT('          confinement energy = ', F14.5, ' Hartree a.u.')
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE print_environ_energies
    !------------------------------------------------------------------------------------
    !>
    !! If Gaussian nuclei are used, write out the corresponding potential shift
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE print_environ_potential_shift()
        !--------------------------------------------------------------------------------
        !
        IF (setup%lsmearedions) &
            WRITE (program_unit, 1100) env%environment_ions%potential_shift * RYTOEV
        !
1100    FORMAT(/, 5(' '), &
                'the potential shift due to the Gaussian-smeared nuclei is ', &
                F10.4, ' ev')
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE print_environ_potential_shift
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE print_environ_potential_warning()
        !--------------------------------------------------------------------------------
        !
        IF (setup%need_pbc_correction) WRITE (program_unit, 1200)
        !
1200    FORMAT(/, &
                5(' '), 'WARNING: you are using the parabolic pbc correction;', /, &
                5(' '), '         the potential shift above must be added to ', /, &
                5(' '), '         band and Fermi energies.')
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE print_environ_potential_warning
    !------------------------------------------------------------------------------------
    !!
    !> Write out the main parameters of Environ calculations, summarizing
    !! the input keywords (some info also on internal vs input units).
    !! Called by summary.f90
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE print_environ_summary()
        !--------------------------------------------------------------------------------
        !
        IF (ionode .AND. prog == 'PW') CALL setup%printout()
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE print_environ_summary
    !------------------------------------------------------------------------------------
    !>
    !! Write out the time informations of the Environ dependent calculations.
    !! Called by print_clock_pw.f90
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE print_environ_clocks(passed_unit)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN), OPTIONAL :: passed_unit
        !
        INTEGER :: actual_unit
        !
        !--------------------------------------------------------------------------------
        !
        IF (PRESENT(passed_unit)) THEN
            actual_unit = passed_unit
        ELSE
            actual_unit = program_unit
        END IF
        !
        WRITE (actual_unit, *)
        WRITE (actual_unit, '(5X,"Environ routines")')
        !
        !--------------------------------------------------------------------------------
        ! Dielectric subroutines
        !
        IF (setup%lelectrostatic) THEN
            !
            CALL env_print_clock('calc_eelect')
            !
            CALL env_print_clock('calc_velect')
            !
            CALL env_print_clock('calc_vgcs')
            !
            CALL env_print_clock('dielectric')
            !
            CALL env_print_clock('electrolyte')
            !
            CALL env_print_clock('calc_felect')
            !
        END IF
        !
        IF (setup%lsemiconductor) CALL env_print_clock('calc_vms')
        !
        !--------------------------------------------------------------------------------
        ! TDDFT
        !
        IF (setup%ltddfpt) CALL env_print_clock('calc_vsolvent_tddfpt')
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE print_environ_clocks
    !------------------------------------------------------------------------------------
    !>
    !! Sets the output file target #TODO do we need this routine?
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE update_output_program_unit(program_unit_in)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: program_unit_in
        !
        !--------------------------------------------------------------------------------
        !
        program_unit = program_unit_in
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE update_output_program_unit
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE environ_QE_interface
!----------------------------------------------------------------------------------------
