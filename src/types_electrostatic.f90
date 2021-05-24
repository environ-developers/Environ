! Copyright (C) 2018 ENVIRON (www.quantum-environment.org)
!
!    This file is part of Environ version 1.1
!
!    Environ 1.1 is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 2 of the License, or
!    (at your option) any later version.
!
!    Environ 1.1 is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more detail, either the file
!    `License' in the root directory of the present distribution, or
!    online at <http://www.gnu.org/licenses/>.
!
!----------------------------------------------------------------------------------------
!
! Authors: Oliviero Andreussi (Department of Physics, UNT)
!          Francesco Nattino  (THEOS and NCCR-MARVEL, EPFL)
!          Nicola Marzari     (THEOS and NCCR-MARVEL, EPFL)
!
!----------------------------------------------------------------------------------------
!>
!! Module containing the definition of electrostatic derived data types
!! of the basic routines to handle them
!!
!----------------------------------------------------------------------------------------
MODULE electrostatic_types
    !------------------------------------------------------------------------------------
    !
    USE modules_constants, ONLY: DP
    !
    USE core_types, ONLY: fft_core, oned_analytic_core, core_container
    !
    !------------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    TYPE gradient_solver
        !--------------------------------------------------------------------------------
        !
        LOGICAL :: lconjugate
        CHARACTER(LEN=80) :: step_type
        REAL(DP) :: step
        INTEGER :: maxstep
        CHARACTER(LEN=80) :: preconditioner
        CHARACTER(LEN=80) :: screening_type
        REAL(DP) :: screening
        REAL(DP) :: tol
        !
        !--------------------------------------------------------------------------------
    END TYPE gradient_solver
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    TYPE iterative_solver
        !--------------------------------------------------------------------------------
        !
        CHARACTER(LEN=80) :: mix_type
        REAL(DP) :: mix
        INTEGER :: maxiter
        INTEGER :: ndiis
        REAL(DP) :: tol
        !
        !--------------------------------------------------------------------------------
    END TYPE iterative_solver
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    TYPE newton_solver
        !--------------------------------------------------------------------------------
        !
        INTEGER :: maxiter
        REAL(DP) :: tol
        !
        !--------------------------------------------------------------------------------
    END TYPE newton_solver
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    TYPE electrostatic_solver
        !--------------------------------------------------------------------------------
        !
        CHARACTER(LEN=80) :: type_
        CHARACTER(LEN=80) :: auxiliary
        LOGICAL :: use_direct
        LOGICAL :: use_gradient
        !
        TYPE(gradient_solver), POINTER :: gradient => NULL()
        !
        LOGICAL :: use_iterative
        TYPE(iterative_solver), POINTER :: iterative => NULL()
        !
        LOGICAL :: use_newton
        TYPE(newton_solver), POINTER :: newton => NULL()
        !
        ! #TODO future work
        !
        ! LOGICAL :: use_lbfgs
        ! TYPE( lbfgs_solver ) :: lbfgs
        !
        !--------------------------------------------------------------------------------
    END TYPE electrostatic_solver
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    TYPE electrostatic_setup
        !--------------------------------------------------------------------------------
        !
        CHARACTER(LEN=80) :: problem
        !
        TYPE(electrostatic_solver), POINTER :: solver => NULL()
        !
        TYPE(core_container), POINTER :: core => NULL()
        !
        LOGICAL :: nested_problem
        TYPE(electrostatic_setup), POINTER :: inner => NULL()
        !
        !--------------------------------------------------------------------------------
    END TYPE electrostatic_setup
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
    ! Keeping imports private
    !
    PRIVATE :: DP, fft_core, oned_analytic_core, core_container
    !
    !------------------------------------------------------------------------------------
END MODULE electrostatic_types
!----------------------------------------------------------------------------------------
