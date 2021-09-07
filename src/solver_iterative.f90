!----------------------------------------------------------------------------------------
!
! Copyright (C) 2018-2021 ENVIRON (www.quantum-environ.org)
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
! Authors: Francesco Nattino  (THEOS and NCCR-MARVEL, EPFL)
!          Oliviero Andreussi (Department of Physics, UNT)
!          Nicola Marzari     (THEOS and NCCR-MARVEL, EPFL)
!          Edan Bainglass     (Department of Physics, UNT)
!
!----------------------------------------------------------------------------------------
!>
!!
!----------------------------------------------------------------------------------------
MODULE class_solver_iterative
    !------------------------------------------------------------------------------------
    !
    USE environ_param, ONLY: DP
    !
    USE class_core_container_electrostatics
    !
    USE class_solver
    USE class_solver_direct
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
    TYPE, EXTENDS(solver_direct), PUBLIC :: solver_iterative
        !--------------------------------------------------------------------------------
        !
        CHARACTER(LEN=80) :: auxiliary
        REAL(DP) :: tol
        INTEGER :: maxiter
        !
        !--------------------------------------------------------------------------------
    CONTAINS
        !--------------------------------------------------------------------------------
        !
        PROCEDURE :: init_iterative => init_solver_iterative
        !
        !--------------------------------------------------------------------------------
    END TYPE solver_iterative
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
    SUBROUTINE init_solver_iterative(this, cores, maxiter, tol_in, auxiliary)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(container_electrostatics), TARGET, INTENT(IN) :: cores
        INTEGER, INTENT(IN) :: maxiter
        REAL(DP), INTENT(IN) :: tol_in
        CHARACTER(LEN=80), INTENT(IN), OPTIONAL :: auxiliary
        !
        CLASS(solver_iterative), INTENT(INOUT) :: this
        !
        !--------------------------------------------------------------------------------
        !
        CALL this%init_cores(cores)
        !
        this%maxiter = maxiter
        this%tol = tol_in
        !
        IF (PRESENT(auxiliary)) this%auxiliary = auxiliary
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_solver_iterative
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE class_solver_iterative
!----------------------------------------------------------------------------------------
