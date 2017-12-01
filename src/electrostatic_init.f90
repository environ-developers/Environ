!
! Copyright (C) 2007-2008 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! original version by O. Andreussi and N. Marzari (MIT)
!
!----------------------------------------------------------------------------
MODULE electrostatic_init
!----------------------------------------------------------------------------

  USE environ_types
  USE electrostatic_types
  USE electrostatic_base

  PRIVATE

  PUBLIC :: set_electrostatic_base, electrostatic_initbase, electrostatic_initcell, &
       & electrostatic_initions, electrostatic_clean

CONTAINS
  !
  ! ... the following routine copies input variables read in input
  ! ... to global variables kept in this module
  !
  SUBROUTINE set_electrostatic_base                         &
       ( problem, tol, solver_type, auxiliary,              &
       step_type, step, maxstep, mix_type, ndiis, mix,      &
       preconditioner, screening_type, screening,           &
       core_type, boundary_core_, ifdtype, nfdpoint,        &
       assume_isolated, pbc_correction, pbc_dim_, pbc_axis_,&
       nspin, prog )
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=20)   :: sub_name = ' set_electrostatic_base '
    INTEGER, INTENT(IN) :: maxstep, ndiis, ifdtype, nfdpoint,       &
         pbc_dim_, pbc_axis_, nspin
    REAL(DP), INTENT(IN) :: tol, step, mix, screening
    CHARACTER( LEN = * ), INTENT(IN) :: problem, solver_type,       &
         auxiliary, step_type, mix_type, preconditioner,            &
         screening_type, core_type, boundary_core_,                 &
         assume_isolated, pbc_correction, prog
    !
    INTEGER :: i
    CHARACTER( LEN = 80 ) :: local_type
    !
    ! Initial setup of core flags
    !
    lfd = .FALSE.
    loned_analytic = .FALSE.
    lqe_fft = .FALSE.
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
       CASE DEFAULT
          CALL errore(sub_name,'Option not yet implemented',1)
       END SELECT
       !
    ELSE
       !
       ! or, in case not provided, check hostcode keyword
       !
       SELECT CASE( TRIM( ADJUSTL(assume_isolated) ) )
       CASE( 'slabx' )
          need_pbc_correction = .TRUE.
          loned_analytic = .TRUE.
          pbc_dim = 2
          pbc_axis = 1
       CASE( 'slaby' )
          need_pbc_correction = .TRUE.
          loned_analytic = .TRUE.
          pbc_dim = 2
          pbc_axis = 2
       CASE( 'slabz' )
          need_pbc_correction = .TRUE.
          loned_analytic = .TRUE.
          pbc_dim = 2
          pbc_axis = 3
       CASE( 'pcc' )
          need_pbc_correction = .TRUE.
          loned_analytic = .TRUE.
          pbc_dim = 0
          pbc_axis = 3
       END SELECT
       !
    END IF
    !
    IF ( need_pbc_correction ) THEN
       IF ( loned_analytic ) THEN
          local_type = "1da"
          CALL init_electrostatic_core( type = local_type, oned_analytic = oned_analytic , core = pbc_core )
       ENDIF
    ENDIF
    !
    ! Set up main (outer) core
    !
    CALL create_electrostatic_core( outer_core )
    SELECT CASE ( core_type )
    CASE ( 'fft' )
       lqe_fft = .TRUE.
       CALL init_electrostatic_core( type = core_type, qe_fft = qe_fft, core = outer_core )
    CASE( '1d-analytic', '1da' )
       loned_analytic = .TRUE.
       CALL init_electrostatic_core( type = core_type, oned_analytic = oned_analytic, core = outer_core )
    CASE DEFAULT
       CALL errore(sub_name,'Unexpected value for electrostatic core_type keyword',1)
    END SELECT
    IF ( need_pbc_correction ) CALL add_correction( correction = pbc_core, core = outer_core )
    !
    ! Set up active numerical cores
    !
    IF ( lfd ) CALL init_fd_core( ifdtype, nfdpoint, fd )
    IF ( lqe_fft ) CALL init_qe_fft_core( qe_fft, assume_isolated, nspin )
    IF ( loned_analytic ) CALL init_oned_analytic_core_first( pbc_dim, pbc_axis, oned_analytic )
    !
    ! Initial setup of solver flags
    !
    lgradient = .FALSE.
    literative = .FALSE.
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
       CALL init_electrostatic_solver( type = solver_type, auxiliary = auxiliary, &
            & gradient = gradient, solver = outer_solver )
    CASE ( 'iterative' )
       literative = .TRUE.
       CALL init_electrostatic_solver( type = solver_type, auxiliary = auxiliary, &
            & iterative = iterative, solver = outer_solver )
    CASE DEFAULT
       CALL errore(sub_name,'Unexpected value for electrostatic solver keyword',1)
    END SELECT
    !
    ! Set up active solvers
    !
    IF ( lgradient ) CALL init_gradient_solver( solver_type, tol, step_type, step, maxstep, &
         & preconditioner, screening_type, screening, gradient )
    IF ( literative ) CALL init_iterative_solver( tol, mix_type, mix, maxstep, ndiis, iterative )
    !
    ! Set reference setup according to calling program
    !
    CALL create_electrostatic_setup( reference )
    SELECT CASE ( prog )
    CASE ( 'PW', 'CP', 'TD' )
       local_type = "poisson"
    CASE DEFAULT
       CALL errore(sub_name,'Unexpected name of host code',1)
    END SELECT
    CALL init_electrostatic_setup( local_type , reference_solver, reference_core, reference )
    !
    ! Create outer setup
    !
    CALL create_electrostatic_setup( outer )
    CALL init_electrostatic_setup( problem, outer_solver, outer_core, outer )
    !
    ! Set logical flags according to electrostatic set up
    !
    need_gradient = .FALSE.
    need_factsqrt = .FALSE.
    need_auxiliary = .FALSE.
    linearized = .FALSE.
    !
    CALL set_electrostatic_flags( reference, need_auxiliary, need_gradient, need_factsqrt, linearized )
    CALL set_electrostatic_flags( outer, need_auxiliary, need_gradient, need_factsqrt, linearized )
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

    LOGICAL, INTENT(IN) :: lflag

    CALL destroy_electrostatic_setup( lflag, outer )
    CALL destroy_electrostatic_setup( lflag, reference )
    CALL destroy_electrostatic_core( lflag, outer_core )
    CALL destroy_electrostatic_core( lflag, reference_core )
    IF ( need_pbc_correction ) CALL destroy_electrostatic_core( lflag, pbc_core )
    CALL destroy_electrostatic_solver( lflag, outer_solver )
    CALL destroy_electrostatic_solver( lflag, reference_solver )
    IF ( lfd ) CALL destroy_fd_core( lflag, fd )
    IF ( loned_analytic ) CALL destroy_oned_analytic_core( lflag, oned_analytic )

    RETURN

!--------------------------------------------------------------------
  END SUBROUTINE electrostatic_clean
!--------------------------------------------------------------------
END MODULE electrostatic_init
