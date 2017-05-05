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
  USE environ_types
  !
  SAVE
  !
  ! Global parameters
  !
  CHARACTER (LEN=80) ::            &
       problem
  REAL (DP) ::                     &
       tolvelect,                  &
       tolrhoaux
  !
  ! Driver's parameters
  !
  CHARACTER (LEN=80) ::            &
       solver,                     &
       auxiliary,                  &
       step_type
  REAL (DP) ::                     &
       step
  !
  ! Iterative driver's parameters (OBSOLETE)
  !
  CHARACTER (LEN=80) ::            &
       mix_type
  REAL (DP) ::                     &
       mix
  INTEGER ::                       &
       ndiis
  !
  ! Preconditioner's parameters
  !
  CHARACTER (LEN=80) ::            &
       preconditioner,             &
       screening_type
  REAL (DP) ::                     &
       screening
  !
  ! Numerical core's parameters
  !
  CHARACTER (LEN=80) ::            &
       core
  !
  ! Core used for numerical derivatives of dielectric
  !
  CHARACTER (LEN=80) ::            &
       dielectric_core
  !
  ! Finite differences' parameters
  !
  INTEGER ::                        &
       ifdtype,                     &
       nfdpoint,                    &
       ncfd
  INTEGER, ALLOCATABLE ::           &
       icfd(:)
  !
  ! Boundary conditions
  !
  INTEGER, DIMENSION(3) ::         &
       bcindex
  REAL (DP), DIMENSION(3) ::       &
       bcplus,                     &
       bcminus
  !
  ! Logical flags
  !
  LOGICAL ::                       &
       need_gradient,              &
       need_factsqrt,              &
       need_gradlog
  !
  ! Computed physical variables
  !
  REAL (KIND=DP), ALLOCATABLE ::   &
       rhoaux(:),                  &
       rhopol(:),                  &
       rhoioncc(:)
  !
  CONTAINS
    !
    ! ... the following routine copies input variables read in input
    ! ... to global variables kept in this module
    !
    SUBROUTINE set_electrostatic_base &
         ( problem_, tolvelect_, tolrhoaux_,                  &
           solver_, auxiliary_, step_type_, step_,            &
           mix_type_, ndiis_, mix_,                           &
           preconditioner_, screening_type_, screening_,      &
           core_, dielectric_core_, ifdtype_, nfdpoint_,      &
           bcindex_, bcplus_, bcminus_ )
      !
      IMPLICIT NONE
      !
      CHARACTER(LEN=20)   :: sub_name = ' set_electrostatic_base '
      INTEGER, INTENT(IN) :: ndiis_, ifdtype_, nfdpoint_, bcindex_(3)
      REAL(DP), INTENT(IN) :: tolvelect_, tolrhoaux_, step_, mix_,        &
                             screening_, bcplus_(3), bcminus_(3)
      CHARACTER( LEN = * ), INTENT(IN) :: problem_, solver_, auxiliary_,  &
                             step_type_, mix_type_, preconditioner_,      &
                             screening_type_, core_, dielectric_core_
      !
      INTEGER :: i
      !
      problem = problem_
      tolvelect = tolvelect_
      tolrhoaux = tolrhoaux_
      !
      solver = solver_
      auxiliary = auxiliary_
      step_type = step_type_
      step = step_
      !
      mix_type = mix_type_
      ndiis = ndiis_
      mix = mix_
      !
      preconditioner = preconditioner_
      screening_type = screening_type_
      screening = screening_
      !
      core = core_
      dielectric_core = dielectric_core_
      ifdtype = ifdtype_
      nfdpoint = nfdpoint_
      !
      bcindex = bcindex_
      bcplus = bcplus_
      bcminus = bcminus_
      !
      ! Set finite differences tools
      !
      ALLOCATE( icfd(-nfdpoint:nfdpoint) )
      CALL init_fd_gradient( ifdtype, nfdpoint, ncfd, icfd )
      !
      ! Set logical flags according to electrostatic set up
      !
      need_gradient = .FALSE.
      need_factsqrt = .FALSE.
      need_gradlog  = .FALSE.
      !
      SELECT CASE ( problem )
         !
      CASE ( 'free' )
         !
      CASE ( 'generalized' )
         !
         SELECT CASE ( solver )
            !
         CASE ( 'cg','sd','iterative' )
            !
            SELECT CASE ( auxiliary )
               !
            CASE ( 'none' )
               !
               SELECT CASE ( preconditioner )
                  !
               CASE ( 'sqrt' )
                  !
                  need_factsqrt = .TRUE.
                  !
               CASE ( 'left', 'none' )
                  !
                  need_gradient = .TRUE.
                  !
               END SELECT
               !
            CASE ( 'full' )
               !
               need_gradlog = .TRUE.
               !
            END SELECT
            !
         END SELECT
         !
      END SELECT
      !
      RETURN
      !
    END SUBROUTINE set_electrostatic_base
    !
    SUBROUTINE set_electrostatic_environ()

      ! Passes informations stored in electrostatic_base to environ_base
      ! and vice-versa

      USE environ_base, ONLY : lstatic, loptical, static, optical

      IMPLICIT NONE

      IF ( lstatic ) CALL init_environ_dielectric_flags( need_gradient, need_factsqrt, need_gradlog, static )
      IF ( loptical ) CALL init_environ_dielectric_flags( need_gradient, need_factsqrt, need_gradlog, optical )

      RETURN

    END SUBROUTINE set_electrostatic_environ
    !
END MODULE electrostatic_base
