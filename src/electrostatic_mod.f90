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
! Authors: Oliviero Andreussi (Department of Physics, UNT)
!          Francesco Nattino  (THEOS and NCCR-MARVEL, EPFL)
!          Nicola Marzari     (THEOS and NCCR-MARVEL, EPFL)
!
!----------------------------------------------------------------------------------------
!>
!! Module used to store global electrostatic-related variables and parameters.
!! Contains all the main variables needed for the electrostatic solvers.
!!
!----------------------------------------------------------------------------------------
MODULE electrostatic_base
    !------------------------------------------------------------------------------------
    !
    USE electrostatic_types
    !
    SAVE
    !
    TYPE(electrostatic_setup) :: inner, outer, reference ! Electrostatic setup(s)
    !
    TYPE(electrostatic_core) :: outer_core, inner_core, reference_core, pbc_core
    ! Electrostatic core(s)
    !
    TYPE(electrostatic_solver) :: outer_solver, inner_solver, reference_solver
    ! Electrostatic solver(s)
    !
    !------------------------------------------------------------------------------------
    ! Internal setup of numerical solvers
    !
    LOGICAL :: lgradient, lconjugate
    TYPE(gradient_solver) :: gradient, inner_gradient
    LOGICAL :: literative
    TYPE(iterative_solver) :: iterative, inner_iterative
    LOGICAL :: lnested, lnewton
    TYPE(newton_solver) :: newton
    !
    !------------------------------------------------------------------------------------
    ! General setup of periodic boundary conditions
    !
    LOGICAL :: need_pbc_correction, need_electrolyte, need_semiconductor, &
               need_outer_loop
    !
    INTEGER :: pbc_dim, pbc_axis
    !
    !------------------------------------------------------------------------------------
    !
    LOGICAL :: need_gradient, need_factsqrt, need_auxiliary
    ! Logical flags that need to be used outside
    !
    !------------------------------------------------------------------------------------
END MODULE electrostatic_base
!----------------------------------------------------------------------------------------
