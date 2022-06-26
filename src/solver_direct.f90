!----------------------------------------------------------------------------------------
!
! Copyright (C) 2018-2021 ENVIRON (www.quantum-environ.org)
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
! Authors: Francesco Nattino  (THEOS and NCCR-MARVEL, EPFL)
!          Oliviero Andreussi (Department of Physics, UNT)
!          Nicola Marzari     (THEOS and NCCR-MARVEL, EPFL)
!          Edan Bainglass     (Department of Physics, UNT)
!
!----------------------------------------------------------------------------------------
!>
!!
!----------------------------------------------------------------------------------------
MODULE class_solver_direct
    !------------------------------------------------------------------------------------
    !
    USE class_io, ONLY: io
    !
    USE environ_param, ONLY: DP
    !
    USE class_density
    USE class_gradient
    !
    USE class_core_container
    !
    USE class_solver
    !
    USE class_charges
    USE class_electrolyte
    USE class_semiconductor
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
    TYPE, EXTENDS(electrostatic_solver), PUBLIC :: solver_direct
        !--------------------------------------------------------------------------------
        !
        CHARACTER(LEN=80) :: corrections_method = 'none'
        !
        !--------------------------------------------------------------------------------
    CONTAINS
        !--------------------------------------------------------------------------------
        !
        PROCEDURE, PRIVATE :: create_solver_direct
        PROCEDURE :: init => init_solver_direct
        !
        PROCEDURE :: poisson_charges
        PROCEDURE :: poisson_density
        !
        PROCEDURE :: grad_poisson_charges
        PROCEDURE :: grad_poisson_density
        !
        !--------------------------------------------------------------------------------
    END TYPE solver_direct
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
    SUBROUTINE create_solver_direct(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(solver_direct), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: routine = 'create_solver_direct'
        !
        !--------------------------------------------------------------------------------
        !
        this%solver_type = 'direct'
        this%corrections_method = 'none'
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE create_solver_direct
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_solver_direct(this, cores, corr_method)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(core_container), INTENT(IN) :: cores
        CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: corr_method
        !
        CLASS(solver_direct), INTENT(INOUT) :: this
        !
        !--------------------------------------------------------------------------------
        !
        CALL this%create_solver_direct()
        !
        CALL this%set_cores(cores)
        !
        IF (PRESENT(corr_method)) this%corrections_method = corr_method
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_solver_direct
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !                                   SOLVER METHODS
    !
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE poisson_charges(this, charges, v)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(solver_direct), INTENT(IN) :: this
        TYPE(environ_charges), INTENT(IN) :: charges
        !
        TYPE(environ_density), INTENT(INOUT) :: v
        !
        CHARACTER(LEN=80) :: routine = 'poisson_charges'
        !
        !--------------------------------------------------------------------------------
        !
        IF (.NOT. ASSOCIATED(charges%density%cell, v%cell)) &
            CALL io%error(routine, "Mismatch in domains of charges and potential", 1)
        !
        !--------------------------------------------------------------------------------
        !
        CALL this%cores%electrostatics%poisson(charges%density, v)
        !
        !--------------------------------------------------------------------------------
        ! PBC corrections, if needed
        !
        IF (this%cores%has_corrections) THEN
            !
            ASSOCIATE (correction => this%cores%corrections, &
                       density => charges%density, &
                       electrolyte => charges%electrolyte, &
                       semiconductor => charges%semiconductor)
                !
                SELECT CASE (this%corrections_method)
                    !
                CASE ('parabolic')
                    CALL correction%potential(density, v)
                    !
                CASE ('gcs') ! gouy-chapman-stern
                    !
                    IF (.NOT. ASSOCIATED(charges%electrolyte)) &
                        CALL io%error(routine, &
                                      "Missing electrolyte for electrochemical boundary correction", 1)
                    !
                    CALL correction%potential(electrolyte%base, density, v)
                    !
                CASE ('ms') ! mott-schottky
                    !
                    IF (.NOT. ASSOCIATED(charges%semiconductor)) &
                        CALL io%error(routine, &
                                      "Missing semiconductor for electrochemical boundary correction", 1)
                    !
                    CALL correction%potential(semiconductor%base, density, v)
                    !
                CASE ('ms-gcs') ! mott-schottky + gouy-chapman-stern
                    !
                    IF (.NOT. ASSOCIATED(charges%semiconductor)) &
                        CALL io%error(routine, &
                                      "Missing semiconductor for electrochemical boundary correction", 1)
                    !
                    CALL correction%potential(electrolyte%base, semiconductor%base, &
                                              density, v)
                    !
                CASE DEFAULT
                    CALL io%error(routine, "Unexpected corrections method", 1)
                    !
                END SELECT
                !
            END ASSOCIATE
            !
        END IF
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE poisson_charges
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE poisson_density(this, charges, v, electrolyte, semiconductor)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(solver_direct), INTENT(IN) :: this
        TYPE(environ_density), INTENT(IN) :: charges
        TYPE(environ_electrolyte), OPTIONAL, INTENT(IN) :: electrolyte
        !
        TYPE(environ_density), INTENT(INOUT) :: v
        TYPE(environ_semiconductor), OPTIONAL, INTENT(INOUT) :: semiconductor
        !
        TYPE(environ_density) :: local
        !
        CHARACTER(LEN=80) :: routine = 'poisson_density'
        !
        !--------------------------------------------------------------------------------
        !
        IF (.NOT. ASSOCIATED(charges%cell, v%cell)) &
            CALL io%error(routine, "Mismatch in domains of charges and potential", 1)
        !
        !--------------------------------------------------------------------------------
        ! Using a local variable for the potential because the routine may be
        ! called with the same argument for charges and potential
        !
        CALL local%init(charges%cell)
        !
        CALL this%cores%electrostatics%poisson(charges, local)
        !
        !--------------------------------------------------------------------------------
        ! PBC corrections, if needed
        !
        IF (this%cores%has_corrections) THEN
            !
            ASSOCIATE (correction => this%cores%corrections)
                !
                SELECT CASE (this%corrections_method)
                    !
                CASE ('parabolic')
                    CALL correction%potential(charges, local)
                    !
                CASE ('gcs') ! gouy-chapman-stern
                    !
                    IF (.NOT. PRESENT(electrolyte)) &
                        CALL io%error(routine, &
                                      "Missing electrolyte for electrochemical boundary correction", 1)
                    !
                    CALL correction%potential(electrolyte%base, charges, local)
                    !
                CASE ('ms') ! mott-schottky
                    !
                    IF (.NOT. PRESENT(semiconductor)) &
                        CALL io%error(routine, &
                                      "Missing semiconductor for electrochemical boundary correction", 1)
                    !
                    CALL correction%potential(semiconductor%base, charges, local)
                    !
                CASE ('ms-gcs') ! mott-schottky + gouy-chapman-stern
                    !
                    IF (.NOT. PRESENT(semiconductor)) &
                        CALL io%error(routine, &
                                      "Missing semiconductor for electrochemical boundary correction", 1)
                    !
                    IF (.NOT. PRESENT(electrolyte)) &
                        CALL io%error(routine, &
                                      "Missing electrolyte for electrochemical boundary correction", 1)
                    !
                    CALL correction%potential(electrolyte%base, semiconductor%base, &
                                              charges, local)
                    !
                CASE DEFAULT
                    CALL io%error(routine, "Unexpected corrections method", 1)
                    !
                END SELECT
                !
            END ASSOCIATE
            !
        END IF
        !
        v%of_r = local%of_r ! only update the potential at the end
        !
        CALL local%destroy()
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE poisson_density
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE grad_poisson_charges(this, charges, grad_v)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(solver_direct), INTENT(IN) :: this
        TYPE(environ_charges), INTENT(IN) :: charges
        !
        TYPE(environ_gradient), INTENT(INOUT) :: grad_v
        !
        CHARACTER(LEN=80) :: routine = 'grad_poisson_charges'
        !
        !--------------------------------------------------------------------------------
        !
        IF (.NOT. ASSOCIATED(charges%density%cell, grad_v%cell)) &
            CALL io%error(routine, "Mismatch in domains of charges and gradient", 1)
        !
        !--------------------------------------------------------------------------------
        !
        CALL this%cores%electrostatics%grad_poisson(charges%density, grad_v)
        !
        !--------------------------------------------------------------------------------
        ! PBC corrections, if needed
        !
        IF (this%cores%has_corrections) THEN
            !
            ASSOCIATE (correction => this%cores%corrections, &
                       density => charges%density, &
                       electrolyte => charges%electrolyte, &
                       semiconductor => charges%semiconductor)
                !
                SELECT CASE (this%corrections_method)
                    !
                CASE ('parabolic')
                    CALL correction%grad_potential(density, grad_v)
                    !
                CASE ('gcs') ! gouy-chapman-stern
                    !
                    IF (.NOT. ASSOCIATED(charges%electrolyte)) &
                        CALL io%error(routine, &
                                      "Missing electrolyte for electrochemical boundary correction", 1)
                    !
                    CALL correction%grad_potential(electrolyte%base, density, grad_v)
                    !
                CASE ('ms') ! mott-schottky
                    !
                    IF (.NOT. ASSOCIATED(charges%semiconductor)) &
                        CALL io%error(routine, &
                                      "Missing semiconductor for electrochemical boundary correction", 1)
                    !
                    CALL correction%grad_potential(semiconductor%base, density, grad_v)
                    !
                CASE ('ms-gcs') ! mott-schottky + gouy-chapman-stern
                    !
                    IF (.NOT. ASSOCIATED(charges%semiconductor)) &
                        CALL io%error(routine, &
                                      "Missing semiconductor for electrochemical boundary correction", 1)
                    !
                    CALL correction%grad_potential(electrolyte%base, semiconductor%base, &
                                                   density, grad_v)
                    !
                CASE DEFAULT
                    CALL io%error(routine, "Unexpected corrections method", 1)
                    !
                END SELECT
                !
            END ASSOCIATE
            !
        END IF
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE grad_poisson_charges
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE grad_poisson_density(this, charges, grad_v, electrolyte, semiconductor)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(solver_direct), INTENT(IN) :: this
        TYPE(environ_density), INTENT(IN) :: charges
        TYPE(environ_electrolyte), OPTIONAL, INTENT(IN) :: electrolyte
        TYPE(environ_semiconductor), OPTIONAL, INTENT(IN) :: semiconductor
        !
        TYPE(environ_gradient), INTENT(INOUT) :: grad_v
        !
        CHARACTER(LEN=80) :: routine = 'grad_poisson_density'
        !
        !--------------------------------------------------------------------------------
        !
        IF (.NOT. ASSOCIATED(charges%cell, grad_v%cell)) &
            CALL io%error(routine, "Mismatch in domains of charges and gradient", 1)
        !
        !--------------------------------------------------------------------------------
        !
        CALL this%cores%electrostatics%grad_poisson(charges, grad_v)
        !
        !--------------------------------------------------------------------------------
        ! PBC corrections, if needed
        !
        IF (this%cores%has_corrections) THEN
            !
            ASSOCIATE (correction => this%cores%corrections)
                !
                SELECT CASE (this%corrections_method)
                    !
                CASE ('parabolic')
                    CALL correction%grad_potential(charges, grad_v)
                    !
                CASE ('gcs') ! gouy-chapman-stern
                    !
                    IF (.NOT. PRESENT(electrolyte)) &
                        CALL io%error(routine, &
                                      "Missing electrolyte for electrochemical boundary correction", 1)
                    !
                    CALL correction%grad_potential(electrolyte%base, charges, grad_v)
                    !
                CASE ('ms') ! mott-schottky
                    !
                    IF (.NOT. PRESENT(semiconductor)) &
                        CALL io%error(routine, &
                                      "Missing semiconductor for electrochemical boundary correction", 1)
                    !
                    CALL correction%grad_potential(semiconductor%base, charges, grad_v)
                    !
                CASE ('ms-gcs') ! mott-schottky + gouy-chapman-stern
                    !
                    IF (.NOT. PRESENT(semiconductor)) &
                        CALL io%error(routine, &
                                      "Missing semiconductor for electrochemical boundary correction", 1)
                    !
                    CALL correction%grad_potential(electrolyte%base, semiconductor%base, &
                                                   charges, grad_v)
                    !
                CASE DEFAULT
                    CALL io%error(routine, "Unexpected corrections method", 1)
                    !
                END SELECT
                !
            END ASSOCIATE
            !
        END IF
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE grad_poisson_density
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE class_solver_direct
!----------------------------------------------------------------------------------------
