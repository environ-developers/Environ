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
!> Module to initialize electrostatic-related variables
!
! Authors: Oliviero Andreussi (Department of Physics, UNT)
!          Francesco Nattino  (THEOS and NCCR-MARVEL, EPFL)
!
!----------------------------------------------------------------------------
MODULE electrostatic_init
!----------------------------------------------------------------------------
  !
  USE environ_types
  USE electrostatic_types
  USE electrostatic_base
  !
  PRIVATE
  !
  PUBLIC :: set_electrostatic_base, electrostatic_initbase, electrostatic_initcell, &
       & electrostatic_initions, electrostatic_clean
  !
CONTAINS
  !
  !  Subroutine: set_electrostatic_base
  !
  !> Copies input variables read in input to global variables
  !! kept in the electrostatic_base module
!--------------------------------------------------------------------
  SUBROUTINE set_electrostatic_base                         &
       ( problem, tol, solver_type, auxiliary,              &
       step_type, step, maxstep, mix_type, ndiis, mix,      &
       preconditioner, screening_type, screening,           &
       core_type, boundary_core_, ifdtype, nfdpoint,        &
       use_internal_pbc_corr, pbc_correction, pbc_dim_,     &
       pbc_axis_, nspin, prog, inner_tol, inner_solver_type,&
       inner_maxstep, inner_mix )
!--------------------------------------------------------------------
    IMPLICIT NONE
    !
    CHARACTER(LEN=20)   :: sub_name = ' set_electrostatic_base '
    LOGICAL, INTENT(IN) :: use_internal_pbc_corr
    INTEGER, INTENT(IN) :: maxstep, ndiis, ifdtype, nfdpoint,       &
         pbc_dim_, pbc_axis_, nspin, inner_maxstep
    REAL(DP), INTENT(IN) :: tol, step, mix, screening, inner_tol,   &
         inner_mix
    CHARACTER( LEN = * ), INTENT(IN) :: problem, solver_type,       &
         auxiliary, step_type, mix_type, preconditioner,            &
         screening_type, core_type, boundary_core_,                 &
         pbc_correction, prog, inner_solver_type
    !
    INTEGER :: i
    CHARACTER( LEN = 80 ) :: local_type, local_auxiliary, local_problem, inner_problem
    !
    ! Initial setup of core flags
    !
    lfd = .FALSE.
    loned_analytic = .FALSE.
    lqe_fft = .FALSE.
    !
    ! Setup nested scheme if required
    !
    lnested = .FALSE.
    IF ( .NOT. inner_solver_type .EQ. 'none' ) lnested = .TRUE.
    !
    ! Set reference core according to calling program
    !
    CALL create_electrostatic_core( reference_core )
    SELECT CASE ( prog )
    CASE ( 'PW', 'CP', 'TD' )
       lqe_fft = .TRUE.
       local_type = "fft"
       CALL init_electrostatic_core( type = local_type, qe_fft = qe_fft, core = reference_core )
    CASE DEFAULT
       CALL errore(sub_name,'Unexpected name of host code',1)
    END SELECT
    !
    ! Numerical core for boundary derivatives NEED TO MOVE IT OUTSIDE ELECTROSTATIC
    !
    boundary_core = boundary_core_
    SELECT CASE ( TRIM( ADJUSTL( boundary_core ) ) )
    CASE ( 'fd' )
       lfd = .TRUE.
    CASE ( 'fft', 'analytic' )
       lqe_fft = .TRUE.
    END SELECT
    !
    ! Numerical core for periodic boundary correction
    !
    need_pbc_correction = .FALSE.
    need_electrolyte = .FALSE.
    CALL create_electrostatic_core( pbc_core )
    !
    ! first check keywords specfied in input
    !
    IF ( pbc_dim_ .GE. 0 ) THEN
       !
       pbc_dim = pbc_dim_
       pbc_axis = pbc_axis_
       SELECT CASE ( TRIM( ADJUSTL( pbc_correction ) ) )
       CASE ( 'none' )
       CASE ( 'parabolic' )
          need_pbc_correction = .TRUE.
          loned_analytic = .TRUE.
          local_type = '1da'
       CASE ( 'gcs', 'gouy-chapman', 'gouy-chapman-stern' )
          need_pbc_correction = .TRUE.
          need_electrolyte = .TRUE.
          loned_analytic = .TRUE.
          local_type = 'gcs'
       CASE ( 'lgcs' )
          need_pbc_correction = .TRUE.
          need_electrolyte = .TRUE.
          loned_analytic = .TRUE.
          local_type = 'lgcs'
       CASE DEFAULT
          CALL errore(sub_name,'Option not yet implemented',1)
       END SELECT
       !
    END IF
    !
    IF ( need_pbc_correction ) THEN
       IF ( loned_analytic ) &
          CALL init_electrostatic_core( type = local_type, oned_analytic = oned_analytic , core = pbc_core )
    ENDIF
    !
    ! Set up main (outer) core and inner core (if nested scheme)
    !
    CALL create_electrostatic_core( outer_core )
    IF ( lnested ) CALL create_electrostatic_core( inner_core )
    SELECT CASE ( core_type )
    CASE ( 'fft' )
       lqe_fft = .TRUE.
       CALL init_electrostatic_core( type = core_type, qe_fft = qe_fft, core = outer_core )
       IF ( lnested ) CALL init_electrostatic_core( type = core_type, qe_fft = qe_fft, core = inner_core )
    CASE( '1d-analytic', '1da' )
       loned_analytic = .TRUE.
       CALL init_electrostatic_core( type = core_type, oned_analytic = oned_analytic, core = outer_core )
       IF ( lnested ) CALL init_electrostatic_core( type = core_type, oned_analytic = oned_analytic, core = inner_core )
    CASE DEFAULT
       CALL errore(sub_name,'Unexpected value for electrostatic core_type keyword',1)
    END SELECT
    IF ( need_pbc_correction ) CALL add_correction( correction = pbc_core, core = outer_core )
    IF ( lnested .AND. need_pbc_correction ) CALL add_correction( correction = pbc_core, core = inner_core )
    !
    ! Set up active numerical cores
    !
    IF ( lfd ) CALL init_fd_core( ifdtype, nfdpoint, fd )
    IF ( lqe_fft ) CALL init_qe_fft_core( qe_fft, use_internal_pbc_corr, nspin )
    IF ( loned_analytic ) CALL init_oned_analytic_core_first( pbc_dim, pbc_axis, oned_analytic )
    !
    ! Initial setup of solver flags
    !
    lgradient = .FALSE.
    lconjugate = .FALSE.
    literative = .FALSE.
    lnewton = .FALSE.
    !
    ! Set reference solver according to calling program
    !
    CALL create_electrostatic_solver( reference_solver )
    SELECT CASE ( prog )
    CASE ( 'PW', 'CP', 'TD' )
       local_type = "direct"
       CALL init_electrostatic_solver( type = local_type, solver = reference_solver )
    CASE DEFAULT
       CALL errore(sub_name,'Unexpected name of host code',1)
    END SELECT
    !
    ! Set up main (outer) solver
    !
    CALL create_electrostatic_solver( outer_solver )
    SELECT CASE ( solver_type )
    CASE ( 'direct' )
       CALL init_electrostatic_solver( type = solver_type, solver = outer_solver )
    CASE ( 'cg', 'sd' )
       lgradient = .TRUE.
       IF ( TRIM(ADJUSTL(solver_type)) .EQ. "cg" ) lconjugate = .TRUE.
       CALL init_electrostatic_solver( type = solver_type, auxiliary = auxiliary, &
            & gradient = gradient, solver = outer_solver )
    CASE ( 'iterative' )
       literative = .TRUE.
       CALL init_electrostatic_solver( type = solver_type, auxiliary = auxiliary, &
            & iterative = iterative, solver = outer_solver )
    CASE ( 'newton' )
       lnewton = .TRUE.
       CALL init_electrostatic_solver( type = solver_type, auxiliary = auxiliary, &
            & newton = newton, solver = outer_solver )
    CASE DEFAULT
       CALL errore(sub_name,'Unexpected value for electrostatic solver keyword',1)
    END SELECT
    !
    ! Set up active outer solvers
    !
    IF ( lgradient ) CALL init_gradient_solver( lconjugate, tol, step_type, step, maxstep, &
         & preconditioner, screening_type, screening, gradient )
    IF ( literative ) CALL init_iterative_solver( tol, mix_type, mix, maxstep, ndiis, iterative )
    IF ( lnewton ) CALL init_newton_solver( tol, maxstep, newton )
    !
    ! If nested scheme, set up inner solver
    !
    IF ( lnested ) THEN
       !
       lgradient = .FALSE.
       lconjugate = .FALSE.
       literative = .FALSE.
       !
       CALL create_electrostatic_solver( inner_solver )
       !
       SELECT CASE ( solver_type )
       CASE ( 'iterative' )
          IF (auxiliary .EQ. 'ioncc') THEN
             inner_problem = 'generalized'
             SELECT CASE ( inner_solver_type )
             CASE ( 'cg', 'sd' )
                lgradient  = .TRUE.
                IF ( TRIM(ADJUSTL(inner_solver_type)) .EQ. "cg" ) lconjugate = .TRUE.
                CALL init_electrostatic_solver( type = inner_solver_type, gradient = inner_gradient, &
                     & solver = inner_solver )
             CASE ( 'iterative' )
                literative = .TRUE.
                local_auxiliary = 'full'
                CALL init_electrostatic_solver( type = inner_solver_type, auxiliary = local_auxiliary, &
                     & iterative = inner_iterative, solver = inner_solver )
             END SELECT
          ELSE
             CALL errore(sub_name,'Unexpected value for auxiliary charge in nested solver',1)
          END IF
       CASE ( 'newton' )
          inner_problem = 'linpb'
          lgradient  = .TRUE.
          lconjugate = .TRUE.
          CALL init_electrostatic_solver( type = inner_solver_type, gradient = inner_gradient, &
                     & solver = inner_solver )
       END SELECT
       IF ( lgradient ) CALL init_gradient_solver( lconjugate, inner_tol, step_type, step, &
         & inner_maxstep, preconditioner, screening_type, screening, inner_gradient )
       IF ( literative ) CALL init_iterative_solver( inner_tol, mix_type, inner_mix, inner_maxstep, &
         & ndiis, inner_iterative )
    END IF
    !
    ! Set reference setup according to calling program
    !
    CALL create_electrostatic_setup( reference )
    SELECT CASE ( prog )
    CASE ( 'PW', 'CP', 'TD' )
       local_problem = "poisson"
    CASE DEFAULT
       CALL errore(sub_name,'Unexpected name of host code',1)
    END SELECT
    CALL init_electrostatic_setup( local_problem , reference_solver, reference_core, reference )
    !
    ! Create outer setup
    !
    CALL create_electrostatic_setup( outer )
    CALL init_electrostatic_setup( problem, outer_solver, outer_core, outer )
    !
    ! If nested scheme, create inner setup
    !
    IF ( lnested ) THEN
       CALL create_electrostatic_setup( inner )
       CALL init_electrostatic_setup( inner_problem, inner_solver, inner_core, inner )
       CALL add_inner_setup( inner, outer )
    END IF
    !
    ! Set logical flags according to electrostatic set up
    !
    need_gradient = .FALSE.
    need_factsqrt = .FALSE.
    need_auxiliary = .FALSE.
    !
    CALL set_electrostatic_flags( reference, need_auxiliary, need_gradient, need_factsqrt)
    CALL set_electrostatic_flags( outer, need_auxiliary, need_gradient, need_factsqrt)
    IF (lnested ) CALL set_electrostatic_flags( inner, need_auxiliary, need_gradient, need_factsqrt)
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE set_electrostatic_base
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE electrostatic_initbase( cell )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    TYPE( environ_cell ), INTENT(IN) :: cell
    !
    IF ( loned_analytic ) CALL init_oned_analytic_core_second( cell, oned_analytic )
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE electrostatic_initbase
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE electrostatic_initcell( cell )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    TYPE( environ_cell ), INTENT(IN) :: cell
    !
    IF ( loned_analytic ) CALL update_oned_analytic_core_cell( cell%omega, cell%at, oned_analytic )
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE electrostatic_initcell
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE electrostatic_initions( system )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    TYPE( environ_system ), INTENT(IN) :: system
    !
    IF ( loned_analytic ) CALL update_oned_analytic_core_origin( system%pos, oned_analytic )
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE electrostatic_initions
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE electrostatic_clean( lflag )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    LOGICAL, INTENT(IN) :: lflag
    !
    CALL destroy_electrostatic_setup( lflag, outer )
    CALL destroy_electrostatic_setup( lflag, reference )
    CALL destroy_electrostatic_core( lflag, outer_core )
    CALL destroy_electrostatic_core( lflag, reference_core )
    IF ( need_pbc_correction ) CALL destroy_electrostatic_core( lflag, pbc_core )
    CALL destroy_electrostatic_solver( lflag, outer_solver )
    CALL destroy_electrostatic_solver( lflag, reference_solver )
    IF ( lfd ) CALL destroy_fd_core( lflag, fd )
    IF ( loned_analytic ) CALL destroy_oned_analytic_core( lflag, oned_analytic )
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE electrostatic_clean
!--------------------------------------------------------------------
!----------------------------------------------------------------------------
END MODULE electrostatic_init
!----------------------------------------------------------------------------
