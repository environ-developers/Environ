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
MODULE class_semiconductor_base
    !------------------------------------------------------------------------------------
    !
    USE class_io, ONLY: io
    !
    USE environ_param, ONLY: DP
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
    TYPE, PUBLIC :: environ_semiconductor_base
        !--------------------------------------------------------------------------------
        !
        REAL(DP) :: temperature = 0.D0
        REAL(DP) :: permittivity = 0.D0
        REAL(DP) :: carrier_density = 0.D0
        !
        REAL(DP) :: sc_distance = 0.D0
        REAL(DP) :: sc_spread = 0.D0
        !
        REAL(DP) :: electrode_charge = 0.D0
        REAL(DP) :: charge_threshold = 0.D0
        !
        REAL(DP) :: flatband_fermi = 0.D0
        REAL(DP) :: bulk_sc_fermi = 0.D0
        REAL(DP) :: slab_charge = 0.D0
        REAL(DP) :: surf_area_per_sq_cm = 0.D0
        REAL(DP) :: ss_chg = 0.D0
        REAL(DP) :: ss_v_cut = 0.D0
        !
        REAL(DP), ALLOCATABLE :: flatband_pot_planar_avg(:)
        !
        !--------------------------------------------------------------------------------
    CONTAINS
        !--------------------------------------------------------------------------------
        !
        PROCEDURE, PRIVATE :: create => create_environ_semiconductor_base
        PROCEDURE :: init => init_environ_semiconductor_base
        PROCEDURE :: destroy => destroy_environ_semiconductor_base
        !
        PROCEDURE :: running_average
        !
        !--------------------------------------------------------------------------------
    END TYPE environ_semiconductor_base
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
    SUBROUTINE create_environ_semiconductor_base(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_semiconductor_base), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: sub_name = 'create_environ_semiconductor_base'
        !
        !--------------------------------------------------------------------------------
        !
        IF (ALLOCATED(this%flatband_pot_planar_avg)) CALL io%create_error(sub_name)
        !
        !--------------------------------------------------------------------------------
        !
        this%temperature = 0.D0
        this%permittivity = 0.D0
        this%carrier_density = 0.D0
        this%sc_distance = 0.D0
        this%sc_spread = 0.D0
        this%electrode_charge = 0.D0
        this%charge_threshold = 0.D0
        this%flatband_fermi = 0.D0
        this%bulk_sc_fermi = 0.D0
        this%slab_charge = 0.D0
        this%surf_area_per_sq_cm = 0.D0
        this%ss_chg = 0.D0
        this%ss_v_cut = 0.D0
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE create_environ_semiconductor_base
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_environ_semiconductor_base(this, temperature, sc_permittivity, &
                                               sc_carrier_density, sc_electrode_chg, &
                                               sc_distance, sc_spread, sc_chg_thr, &
                                               need_flatband, naxis)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP), INTENT(IN) :: temperature, sc_permittivity, sc_electrode_chg, &
                                sc_carrier_density, sc_distance, sc_spread, sc_chg_thr
        !
        LOGICAL, INTENT(IN) :: need_flatband
        INTEGER, INTENT(IN) :: naxis
        !
        CLASS(environ_semiconductor_base), INTENT(INOUT) :: this
        !
        !--------------------------------------------------------------------------------
        !
        CALL this%create()
        !
        this%temperature = temperature
        this%permittivity = sc_permittivity
        this%carrier_density = sc_carrier_density
        this%sc_distance = sc_distance
        this%sc_spread = sc_spread
        !
        this%carrier_density = this%carrier_density * 1.48D-25
        ! convert carrier density to units of (bohr)^-3
        !
        this%electrode_charge = sc_electrode_chg
        this%charge_threshold = sc_chg_thr
        !
        IF (need_flatband) ALLOCATE (this%flatband_pot_planar_avg(naxis))
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_environ_semiconductor_base
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE destroy_environ_semiconductor_base(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_semiconductor_base), INTENT(INOUT) :: this
        !
        !--------------------------------------------------------------------------------
        !
        IF (ALLOCATED(this%flatband_pot_planar_avg)) &
            DEALLOCATE (this%flatband_pot_planar_avg)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE destroy_environ_semiconductor_base
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !                                  FUNCTION METHODS
    !
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE running_average(this, z_length, naxis, pot, averaged_pot)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_semiconductor_base), INTENT(IN) :: this
        REAL(DP), INTENT(IN) :: z_length
        INTEGER, INTENT(IN) :: naxis
        REAL(DP), INTENT(IN) :: pot(naxis)
        !
        REAL(DP), INTENT(OUT) :: averaged_pot(naxis)
        !
        INTEGER :: i, idx_width, start_idx, stop_idx
        !
        !--------------------------------------------------------------------------------
        ! Averaging bb
        !
        idx_width = INT(this%sc_spread / 2.0 / z_length * naxis)
        !
        DO i = 1, naxis
            start_idx = i - idx_width
            stop_idx = i + idx_width
            !
            IF (start_idx < 1) start_idx = 1
            !
            IF (stop_idx > naxis) stop_idx = naxis
            !
            averaged_pot(i) = SUM(pot(start_idx:stop_idx)) / FLOAT(stop_idx - start_idx)
        END DO
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE running_average
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE class_semiconductor_base
!----------------------------------------------------------------------------------------
