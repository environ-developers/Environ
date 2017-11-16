!
! Copyright (C) 2007-2008 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! The different modules in this file contain all the Environ related variables
! that need to be passed in and out. All module-specific variables are declared
! inside the appropriate modules.
!
! original version by O. Andreussi and N. Marzari (MIT)
!
!----------------------------------------------------------------------------
MODULE electrostatic_base
  !--------------------------------------------------------------------------
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
       lgradient
  TYPE( gradient_solver ) ::           &
       gradient
  LOGICAL ::                           &
       literative
  TYPE( iterative_solver ) ::          &
       iterative
  !
  ! General setup of periodic boundary conditions
  !
  LOGICAL ::                           &
       need_pbc_correction
  INTEGER ::                           &
       pbc_dim,                        &
       pbc_axis
  !
  ! Logical flags that need to be used outside
  !
  LOGICAL ::                       &
       need_gradient,              &
       need_factsqrt,              &
       need_auxiliary,             &
       linearized
  !
END MODULE electrostatic_base
