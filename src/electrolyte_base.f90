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
!          Nicola Marzari     (THEOS and NCCR-MARVEL, EPFL)
!          Edan Bainglass     (Department of Physics, UNT)
!
!----------------------------------------------------------------------------------------
!>
!!
!----------------------------------------------------------------------------------------
MODULE class_electrolyte_base
    !------------------------------------------------------------------------------------
    !
    USE class_io, ONLY: io
    !
    USE environ_param, ONLY: DP, e2, BOHR_RADIUS_SI, AMU_SI, K_BOLTZMANN_RY, fpi
    !
    USE class_cell
    !
    USE class_ioncctype
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
    TYPE, PUBLIC :: environ_electrolyte_base
        !--------------------------------------------------------------------------------
        !
        LOGICAL :: linearized = .FALSE.
        !
        CHARACTER(LEN=80) :: electrolyte_entropy
        INTEGER :: ntyp
        TYPE(environ_ioncctype), ALLOCATABLE :: ioncctype(:)
        !
        REAL(DP) :: temperature
        REAL(DP) :: k2
        REAL(DP) :: cionmax
        REAL(DP) :: permittivity
        !
        REAL(DP) :: distance
        REAL(DP) :: spread
        !
        !--------------------------------------------------------------------------------
    CONTAINS
        !--------------------------------------------------------------------------------
        !
        PROCEDURE, PRIVATE :: create => create_environ_electrolyte_base
        PROCEDURE :: init => init_environ_electrolyte_base
        PROCEDURE :: destroy => destroy_environ_electrolyte_base
        !
        !--------------------------------------------------------------------------------
    END TYPE environ_electrolyte_base
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
    SUBROUTINE create_environ_electrolyte_base(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_electrolyte_base), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: routine = 'create_environ_electrolyte_base'
        !
        !--------------------------------------------------------------------------------
        !
        IF (ALLOCATED(this%ioncctype)) CALL io%create_error(routine)
        !
        !--------------------------------------------------------------------------------
        !
        this%linearized = .FALSE.
        this%electrolyte_entropy = ''
        this%ntyp = 0
        this%temperature = 0.D0
        this%k2 = 0.D0
        this%cionmax = 0.D0
        this%permittivity = 0.D0
        this%distance = 0.D0
        this%spread = 0.D0
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE create_environ_electrolyte_base
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_environ_electrolyte_base(this, ntyp, const, distance, spread, &
                                             temperature, cbulk, cionmax, radius, z, &
                                             electrolyte_entropy, linearized, cell)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        LOGICAL, INTENT(IN) :: linearized
        INTEGER, INTENT(IN) :: ntyp
        CHARACTER(LEN=*), INTENT(IN) :: electrolyte_entropy
        REAL(DP), INTENT(IN) :: const, distance, spread, temperature, cionmax, radius
        REAL(DP), DIMENSION(ntyp), INTENT(IN) :: cbulk, z
        TYPE(environ_cell), INTENT(IN) :: cell
        !
        CLASS(environ_electrolyte_base), INTENT(INOUT) :: this
        !
        INTEGER :: i
        REAL(DP) :: neutral, sumcbulk, sum_cz2, KT
        !
        CHARACTER(LEN=80) :: routine = 'init_environ_electrolyte_base'
        !
        !--------------------------------------------------------------------------------
        !
        CALL this%create()
        !
        !--------------------------------------------------------------------------------
        ! Setup all electrolyte parameters (with checks)
        !
        this%linearized = linearized
        this%ntyp = ntyp
        this%electrolyte_entropy = TRIM(electrolyte_entropy)
        this%temperature = temperature
        this%permittivity = const
        this%distance = distance
        this%spread = spread
        !
        ALLOCATE (this%ioncctype(ntyp))
        !
        neutral = 0.D0
        sum_cz2 = 0.D0
        !
        DO i = 1, ntyp
            !
            CALL this%ioncctype(i)%init(i, cbulk(i), z(i), cell)
            !
            neutral = neutral + cbulk(i) * z(i)
            sum_cz2 = sum_cz2 + this%ioncctype(i)%cbulk * this%ioncctype(i)%z**2
        END DO
        !
        IF (neutral > 1.D-8) &
            CALL io%error(routine, "Bulk electrolyte is not neutral", 1)
        !
        kT = K_BOLTZMANN_RY * temperature
        !
        this%k2 = sum_cz2 / kT * e2 * fpi ! k^2 = eps / lambda_D^2
        !
        this%cionmax = cionmax * BOHR_RADIUS_SI**3 / AMU_SI
        !
        !--------------------------------------------------------------------------------
        ! If not given cionmax in input, but given radius, calculate cionmax
        !
        IF (cionmax == 0.D0 .AND. radius > 0.D0) &
            this%cionmax = 0.64D0 * 3.D0 / fpi / radius**3
        !
        !--------------------------------------------------------------------------------
        ! Check suitability of cionmax value
        !
        sumcbulk = SUM(this%ioncctype%cbulk)
        !
        IF (this%cionmax > 0.D0 .AND. this%cionmax <= sumcbulk) &
            CALL io%error(routine, "cionmax should be larger than the sum of cbulks", 1)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_environ_electrolyte_base
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE destroy_environ_electrolyte_base(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(environ_electrolyte_base), INTENT(INOUT) :: this
        !
        INTEGER :: i
        !
        CHARACTER(LEN=80) :: routine = 'destroy_environ_electrolyte_base'
        !
        !--------------------------------------------------------------------------------
        !
        IF (.NOT. ALLOCATED(this%ioncctype)) CALL io%destroy_error(routine)
        !
        !--------------------------------------------------------------------------------
        !
        DO i = 1, this%ntyp
            CALL this%ioncctype(i)%destroy()
        END DO
        !
        DEALLOCATE (this%ioncctype)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE destroy_environ_electrolyte_base
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE class_electrolyte_base
!----------------------------------------------------------------------------------------
