! Copyright (C) 2018 ENVIRON (www.quantum-environment.org)
!
!    This file is part of Environ version 1.0
!
!    Environ 1.0 is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 2 of the License, or
!    (at your option) any later version.
!
!    Environ 1.0 is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more detail, either the file
!    `License' in the root directory of the present distribution, or
!    online at <http://www.gnu.org/licenses/>.
!
! Module used to store global electrostatic-related variables and parameters
!
! Authors: Oliviero Andreussi (Department of Physics, UNT)
!          Francesco Nattino  (THEOS and NCCR-MARVEL, EPFL)
!          Nicola Marzari     (THEOS and NCCR-MARVEL, EPFL)
!
!----------------------------------------------------------------------------
MODULE electrostatic_base
!----------------------------------------------------------------------------
  !
  ! ... this module contains all the main variables needed for the
  ! ... electrostatic solvers.
  !
  USE electrostatic_types
  !
  SAVE
  !
  ! Electrostatic setup(s)
  !
  TYPE( electrostatic_setup ) ::       &
       inner,                          &
       outer,                          &
       reference
  !
  ! Electrostatic core(s)
  !
  TYPE( electrostatic_core ) ::        &
       outer_core,                     &
       inner_core,                     &
       reference_core,                 &
       pbc_core
  !
  ! Electrostatic solver(s)
  !
  TYPE( electrostatic_solver ) ::      &
       outer_solver,                   &
       inner_solver,                   &
       reference_solver
  !
  ! Core used for numerical derivatives of dielectric
  !
  CHARACTER (LEN=80) ::                &
       boundary_core
  !
  ! Internal setup of numerical cores
  !
  LOGICAL ::                           &
       lfd
  TYPE( fd_core ) ::                   &
       fd
  LOGICAL ::                           &
       lqe_fft
  TYPE( qe_fft_core ) ::               &
       qe_fft
  LOGICAL ::                           &
       loned_analytic
  TYPE( oned_analytic_core ) ::        &
       oned_analytic
  !
  ! Internal setup of numerical solvers
  !
  LOGICAL ::                           &
       lgradient, lconjugate
  TYPE( gradient_solver ) ::           &
       gradient, inner_gradient
  LOGICAL ::                           &
       literative
  TYPE( iterative_solver ) ::          &
       iterative, inner_iterative
  LOGICAL ::                           &
       lnested, lnewton
  TYPE( newton_solver ) ::             &
       newton
  !
  ! General setup of periodic boundary conditions
  !
  LOGICAL ::                           &
       need_pbc_correction,            &
       need_electrolyte,               &
       need_semiconductor
  INTEGER ::                           &
       pbc_dim,                        &
       pbc_axis
  !
  ! Logical flags that need to be used outside
  !
  LOGICAL ::                       &
       need_gradient,              &
       need_factsqrt,              &
       need_auxiliary
  !
END MODULE electrostatic_base
