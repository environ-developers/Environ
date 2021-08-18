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
    USE environ_param, ONLY: DP
    !
    USE class_cell
    USE class_density
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
        LOGICAL :: initialized = .FALSE.
        !
        REAL(DP) :: temperature
        REAL(DP) :: permittivity
        REAL(DP) :: carrier_density
        REAL(DP) :: electrode_charge
        REAL(DP) :: charge_threshold
        REAL(DP) :: slab_charge = 0.D0
        !
        TYPE(environ_function) :: simple
        TYPE(environ_density) :: density
        !
        REAL(DP) :: charge = 0.0_DP
        REAL(DP) :: flatband_fermi = 0.D0
        REAL(DP) :: bulk_sc_fermi = 0.D0
        REAL(DP) :: surf_area_per_sq_cm = 0.D0
        !
        !--------------------------------------------------------------------------------
    CONTAINS
        !--------------------------------------------------------------------------------
        !
        PROCEDURE :: create => create_environ_semiconductor
        PROCEDURE :: init_first => init_environ_semiconductor_first
        PROCEDURE :: init_second => init_environ_semiconductor_second
        PROCEDURE :: update => update_environ_semiconductor
        PROCEDURE :: destroy => destroy_environ_semiconductor
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
        CHARACTER(LEN=80) :: label = 'semiconductor'
        !
        !--------------------------------------------------------------------------------
        !
        this%lupdate = .FALSE.
        !
        CALL this%density%create(label)
        !
        this%charge = 0.D0
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE create_environ_semiconductor
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_environ_semiconductor_first(this, temperature, sc_permittivity, &
                                                sc_carrier_density, sc_electrode_chg, &
                                                sc_distance, sc_spread, sc_chg_thr, &
                                                system)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP), INTENT(IN) :: temperature, sc_permittivity, sc_electrode_chg, &
                                sc_carrier_density, sc_distance, sc_spread, sc_chg_thr
        !
        TYPE(environ_system), TARGET, INTENT(IN) :: system
        !
        CLASS(environ_semiconductor), INTENT(INOUT) :: this
        !
        !--------------------------------------------------------------------------------
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
        this%simple%f_type = 4
        this%simple%pos => system%pos
        this%simple%volume = 1.D0
        this%simple%dim = system%dim
        this%simple%axis = system%axis
        this%simple%width = sc_distance
        this%simple%spread = sc_spread
        !
        this%initialized = .FALSE.
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_environ_semiconductor_first
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_environ_semiconductor_second(this, cell)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_cell), INTENT(IN) :: cell
        !
        CLASS(environ_semiconductor), INTENT(INOUT) :: this
        !
        !--------------------------------------------------------------------------------
        !
        CALL this%density%init(cell)
        !
        this%initialized = .TRUE.
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_environ_semiconductor_second
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
        RETURN
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
        !--------------------------------------------------------------------------------
        !
        IF (this%initialized) THEN
            !
            CALL this%density%destroy()
            !
            this%initialized = .FALSE.
        END IF
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE destroy_environ_semiconductor
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE class_semiconductor
!----------------------------------------------------------------------------------------