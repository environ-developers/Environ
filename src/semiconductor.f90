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
! Authors: Quinn Campbell (Sandia National Laboratories)
!          Ismaila Dabo   (DMSE, Penn State)
!          Edan Bainglass (Department of Physics, UNT)
!
!----------------------------------------------------------------------------------------
!>
!! Module containing the main routines to handle environ_semiconductor
!! derived data types. Environ_semiconductor contains all the specifications
!! and the details of the user defined semiconductor region
!!
!----------------------------------------------------------------------------------------
MODULE class_semiconductor
    !------------------------------------------------------------------------------------
    !
    USE env_base_io, ONLY: ionode, environ_unit, global_verbose
    !
    USE environ_param, ONLY: DP, BOHR_RADIUS_CM
    !
    USE class_cell
    USE class_density
    USE class_function_erfc
    USE class_functions
    !
    USE class_system
    !
    !------------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    PRIVATE
    !
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    TYPE, PUBLIC :: environ_semiconductor
        !--------------------------------------------------------------------------------
        !
        LOGICAL :: lupdate = .FALSE.
        !
        REAL(DP) :: temperature
        REAL(DP) :: permittivity
        REAL(DP) :: carrier_density
        REAL(DP) :: electrode_charge
        REAL(DP) :: charge_threshold
        !
        REAL(DP) :: slab_charge = 0.D0
        REAL(DP) :: charge = 0.D0
        REAL(DP) :: flatband_fermi = 0.D0
        REAL(DP) :: bulk_sc_fermi = 0.D0
        REAL(DP) :: surf_area_per_sq_cm = 0.D0
        !
        TYPE(environ_function_erfc) :: simple
        TYPE(environ_density) :: density
        !
        !--------------------------------------------------------------------------------
    CONTAINS
        !--------------------------------------------------------------------------------
        !
        PROCEDURE :: init => init_environ_semiconductor
        PROCEDURE :: update => update_environ_semiconductor
        PROCEDURE :: destroy => destroy_environ_semiconductor
        !
        PROCEDURE :: printout => print_environ_semiconductor
        !
        !--------------------------------------------------------------------------------
    END TYPE environ_semiconductor
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
CONTAINS
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !                                   ADMIN METHODS
    !
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_environ_semiconductor(this, temperature, sc_permittivity, &
                                          sc_carrier_density, sc_electrode_chg, &
                                          sc_distance, sc_spread, sc_chg_thr, &
                                          system, cell)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP), INTENT(IN) :: temperature, sc_permittivity, sc_electrode_chg, &
                                sc_carrier_density, sc_distance, sc_spread, sc_chg_thr
        !
        TYPE(environ_system), TARGET, INTENT(IN) :: system
        TYPE(environ_cell), INTENT(IN) :: cell
        !
        CLASS(environ_semiconductor), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: local_label = 'semiconductor'
        !
        !--------------------------------------------------------------------------------
        !
        CALL this%density%init(cell, local_label)
        !
        CALL this%simple%init(4, system%axis, system%dim, sc_distance, sc_spread, &
                              1.D0, system%pos)
        !
        this%temperature = temperature
        this%permittivity = sc_permittivity
        this%carrier_density = sc_carrier_density
        !
        this%carrier_density = this%carrier_density * 1.48D-25
        ! convert carrier density to units of (bohr)^-3
        !
        this%electrode_charge = sc_electrode_chg
        this%charge_threshold = sc_chg_thr
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_environ_semiconductor
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE update_environ_semiconductor(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_semiconductor), INTENT(INOUT) :: this
        !
        !--------------------------------------------------------------------------------
        !
        this%charge = this%density%integrate()
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE update_environ_semiconductor
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE destroy_environ_semiconductor(this, lflag)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        LOGICAL, INTENT(IN) :: lflag
        !
        CLASS(environ_semiconductor), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: sub_name = 'destroy_environ_semiconductor'
        !
        !--------------------------------------------------------------------------------
        !
        CALL this%density%destroy()
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE destroy_environ_semiconductor
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !                                   OUTPUT METHODS
    !
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !>
    !! Prints the details of the semiconductor
    !!
    !! Nested objects receive a decremented passed verbose to trigger block printing
    !!
    !! @param verbose       : (INTEGER) adds verbosity to global verbose
    !! @param debug_verbose : (INTEGER) replaces global verbose for debugging
    !! @param unit          : (INTEGER) output target (default = io%debug_unit)
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE print_environ_semiconductor(this, verbose, debug_verbose, unit)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_semiconductor), INTENT(IN) :: this
        INTEGER, INTENT(IN), OPTIONAL :: verbose, debug_verbose, unit
        !
        INTEGER :: base_verbose, local_verbose, passed_verbose, local_unit
        !
        CHARACTER(LEN=80) :: sub_name = 'print_environ_semiconductor'
        !
        !--------------------------------------------------------------------------------
        !
        IF (PRESENT(debug_verbose)) THEN
            base_verbose = debug_verbose
            !
            IF (PRESENT(verbose)) THEN
                local_verbose = verbose
            ELSE
                local_verbose = debug_verbose
            END IF
            !
            passed_verbose = verbose - 1
            !
        ELSE IF (global_verbose > 0) THEN
            base_verbose = global_verbose
            !
            IF (PRESENT(verbose)) THEN
                local_verbose = base_verbose + verbose
            ELSE
                local_verbose = base_verbose
            END IF
            !
            passed_verbose = local_verbose - base_verbose - 1
            !
        ELSE
            RETURN
        END IF
        !
        IF (PRESENT(unit)) THEN
            local_unit = unit
        ELSE
            local_unit = environ_unit
        END IF
        !
        IF (local_verbose >= 1) THEN
            !
            IF (ionode) THEN
                WRITE (local_unit, 1000)
                !
                WRITE (local_unit, 1001) &
                    this%carrier_density / BOHR_RADIUS_CM**3, this%temperature, &
                    this%permittivity
                !
                WRITE (local_unit, 1002) this%charge
            END IF
            !
            IF (local_verbose >= 3) &
                CALL this%density%printout(passed_verbose, debug_verbose, local_unit)
            !
        END IF
        !
        FLUSH (local_unit)
        !
        !--------------------------------------------------------------------------------
        !
1000    FORMAT(/, 4('%'), ' SEMICONDUCTOR ', 64('%'))
        !
1001    FORMAT(/, ' Mott-Schottky:', /, &
                ' dopant concentation   (cm^-3) = ', E15.4, /, &
                ' semiconductor temperature (K) = ', F7.1, /, &
                ' dielectric constant           = ', F7.2)
        !
1002    FORMAT(/, ' total semiconductor charge    = ', F14.7)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE print_environ_semiconductor
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE class_semiconductor
!----------------------------------------------------------------------------------------
