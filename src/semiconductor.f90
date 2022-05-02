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
    USE class_io, ONLY: io
    !
    USE environ_param, ONLY: DP, BOHR_RADIUS_CM
    !
    USE class_cell
    USE class_density
    USE class_function_erfc
    USE class_functions
    !
    USE class_semiconductor_base
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
        TYPE(environ_semiconductor_base) :: base
        !
        !--------------------------------------------------------------------------------
        !
        TYPE(environ_function_erfc) :: simple
        TYPE(environ_density) :: density
        !
        REAL(DP) :: charge = 0.D0
        !
        !--------------------------------------------------------------------------------
    CONTAINS
        !--------------------------------------------------------------------------------
        !
        PROCEDURE, PRIVATE :: create => create_environ_semiconductor
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
    SUBROUTINE create_environ_semiconductor(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_semiconductor), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: routine = 'create_environ_semiconductor'
        !
        !--------------------------------------------------------------------------------
        !
        this%lupdate = .FALSE.
        this%charge = 0.D0
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE create_environ_semiconductor
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_environ_semiconductor(this, temperature, sc_permittivity, &
                                          sc_carrier_density, sc_electrode_chg, &
                                          sc_distance, sc_spread, sc_chg_thr, &
                                          need_flatband, system, cell)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP), INTENT(IN) :: temperature, sc_permittivity, sc_electrode_chg, &
                                sc_carrier_density, sc_distance, sc_spread, sc_chg_thr
        !
        LOGICAL, INTENT(IN) :: need_flatband
        !
        TYPE(environ_system), INTENT(IN) :: system
        TYPE(environ_cell), INTENT(IN) :: cell
        !
        CLASS(environ_semiconductor), INTENT(INOUT) :: this
        !
        !--------------------------------------------------------------------------------
        !
        CALL this%create()
        !
        CALL this%base%init(temperature, sc_permittivity, sc_carrier_density, &
                            sc_electrode_chg, sc_distance, sc_spread, sc_chg_thr, &
                            need_flatband, cell%nr(3))
        !
        CALL this%density%init(cell, 'semiconductor')
        !
        CALL this%simple%init(3, system%axis, system%dim, sc_distance, sc_spread, &
                              1.D0, system%com)
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
        ! Output current state
        !
        CALL this%printout()
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE update_environ_semiconductor
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE destroy_environ_semiconductor(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_semiconductor), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: routine = 'destroy_environ_semiconductor'
        !
        !--------------------------------------------------------------------------------
        !
        CALL this%density%destroy()
        !
        CALL this%simple%destroy()
        !
        CALL this%base%destroy()
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
        CLASS(environ_semiconductor), TARGET, INTENT(IN) :: this
        INTEGER, OPTIONAL, INTENT(IN) :: verbose, debug_verbose, unit
        !
        INTEGER :: base_verbose, local_verbose, passed_verbose, local_unit
        !
        TYPE(environ_semiconductor_base), POINTER :: base
        !
        CHARACTER(LEN=80) :: routine = 'print_environ_semiconductor'
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
        ELSE IF (io%verbosity > 0) THEN
            base_verbose = io%verbosity
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
            local_unit = io%debug_unit
        END IF
        !
        base => this%base
        !
        IF (local_verbose >= 1) THEN
            !
            IF (io%lnode) THEN
                WRITE (local_unit, 1000)
                !
                WRITE (local_unit, 1001) &
                    base%carrier_density / BOHR_RADIUS_CM**3, base%temperature, &
                    base%permittivity
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
1000    FORMAT(/, 4('%'), " SEMICONDUCTOR ", 64('%'))
        !
1001    FORMAT(/, " Mott-Schottky:", /, &
                " dopant concent.    (cm^-3) = ", E18.4, /, &
                " semiconductor temp.    (K) = ", F14.1, /, &
                " dielectric constant        = ", F14.2)
        !
1002    FORMAT(/, " total semiconductor charge = ", F14.7)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE print_environ_semiconductor
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE class_semiconductor
!----------------------------------------------------------------------------------------
