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
    USE class_io, ONLY: io
    !
    USE environ_param, ONLY: DP
    !
    USE class_core_container
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
    TYPE, ABSTRACT, EXTENDS(electrostatic_solver), PUBLIC :: solver_iterative
        !--------------------------------------------------------------------------------
        !
        CHARACTER(LEN=80) :: auxiliary
        REAL(DP) :: tol
        INTEGER :: maxiter
        !
        TYPE(solver_direct), POINTER :: direct
        !
        !--------------------------------------------------------------------------------
    CONTAINS
        !--------------------------------------------------------------------------------
        !
        PROCEDURE, PRIVATE :: create_iterative => create_solver_iterative
        PROCEDURE :: init_iterative => init_solver_iterative
        PROCEDURE :: destroy => destroy_solver_iterative
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
    SUBROUTINE create_solver_iterative(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(solver_iterative), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: routine = 'create_solver_iterative'
        !
        !--------------------------------------------------------------------------------
        !
        IF (ASSOCIATED(this%direct)) CALL io%create_error(routine)
        !
        !--------------------------------------------------------------------------------
        !
        this%auxiliary = ''
        this%tol = 0.D0
        this%maxiter = 0
        !
        NULLIFY (this%direct)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE create_solver_iterative
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_solver_iterative(this, cores, direct, maxiter, tol, auxiliary)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(core_container), INTENT(IN) :: cores
        TYPE(solver_direct), TARGET, INTENT(IN) :: direct
        INTEGER, INTENT(IN) :: maxiter
        REAL(DP), INTENT(IN) :: tol
        CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: auxiliary
        !
        CLASS(solver_iterative), INTENT(INOUT) :: this
        !
        !--------------------------------------------------------------------------------
        !
        CALL this%create_iterative()
        !
        CALL this%set_cores(cores)
        !
        this%direct => direct
        !
        this%maxiter = maxiter
        this%tol = tol
        !
        IF (PRESENT(auxiliary)) this%auxiliary = auxiliary
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_solver_iterative
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE destroy_solver_iterative(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(solver_iterative), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: routine = 'destroy_solver_iterative'
        !
        !--------------------------------------------------------------------------------
        !
        IF (.NOT. ASSOCIATED(this%direct)) CALL io%destroy_error(routine)
        !
        !--------------------------------------------------------------------------------
        !
        CALL this%destroy_cores()
        !
        CALL this%direct%destroy()
        !
        NULLIFY (this%direct)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE destroy_solver_iterative
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE class_solver_iterative
!----------------------------------------------------------------------------------------
