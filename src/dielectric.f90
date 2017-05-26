!
! Copyright (C) 2007-2008 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! original version by O. Andreussi and N. Marzari
!
!--------------------------------------------------------------------
MODULE dielectric
!--------------------------------------------------------------------

  USE environ_types
  USE environ_output
  USE functions
  USE environ_base, ONLY : e2
  USE electrostatic_base, ONLY : dielectric_core, nfdpoint, icfd, ncfd
  !
  IMPLICIT NONE
  !
  SAVE
  !
  PRIVATE
  !
  PUBLIC :: create_environ_dielectric, init_environ_dielectric_first, &
       & set_dielectric_regions, init_environ_dielectric_second, &
       & update_environ_dielectric, destroy_environ_dielectric
  !
CONTAINS
  !
  SUBROUTINE create_environ_dielectric(dielectric)

    IMPLICIT NONE

    TYPE( environ_dielectric ), INTENT(INOUT) :: dielectric

    CHARACTER( LEN=80 ) :: sub_name = 'create_environ_dielectric'
    CHARACTER( LEN=80 ) :: label = ' '

    dielectric%constant = 1.0_DP

    IF ( ALLOCATED( dielectric%regions ) ) CALL errore(sub_name,'Trying to create an already allocated object',1)

    label = 'background'
    CALL create_environ_density( dielectric%background, label  )
    label = 'epsilon'
    CALL create_environ_density( dielectric%epsilon, label )
    label = 'depsilon'
    CALL create_environ_density( dielectric%depsilon, label )

    NULLIFY( dielectric%boundary )

    dielectric%need_gradient = .FALSE.
    label = 'epsilon_gradient'
    CALL create_environ_gradient( dielectric%gradient, label )
    dielectric%need_factsqrt = .FALSE.
    label = 'epsilon_factsqrt'
    CALL create_environ_density( dielectric%factsqrt, label )
    dielectric%need_gradlog = .FALSE.
    label = 'epsilon_gradlog'
    CALL create_environ_gradient( dielectric%gradlog, label )
    RETURN

  END SUBROUTINE create_environ_dielectric

  SUBROUTINE init_environ_dielectric_first( constant, boundary, &
             & need_gradient, need_factsqrt, need_gradlog, dielectric )

    IMPLICIT NONE

    REAL( DP ) :: constant
    LOGICAL, INTENT(IN) :: need_gradient, need_factsqrt, need_gradlog
    TYPE( environ_boundary ), TARGET, INTENT(IN) :: boundary
    TYPE( environ_dielectric ), INTENT(INOUT) :: dielectric

    dielectric%constant = constant

    dielectric%boundary => boundary

    dielectric%need_gradient = need_gradient
    dielectric%need_factsqrt = need_factsqrt
    dielectric%need_gradlog = need_gradlog

    RETURN

  END SUBROUTINE init_environ_dielectric_first

  SUBROUTINE set_dielectric_regions( nregions, epsregion_dim, epsregion_axis, &
       & epsregion_pos, epsregion_width, epsregion_spread, epsregion_eps, dielectric )

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nregions
    INTEGER, DIMENSION(nregions), INTENT(IN) :: epsregion_dim, epsregion_axis
    REAL( DP ), DIMENSION(nregions), INTENT(IN) :: epsregion_width, epsregion_spread, epsregion_eps
    REAL( DP ), DIMENSION(3,nregions), INTENT(IN) :: epsregion_pos
    TYPE( environ_dielectric ), INTENT(INOUT) :: dielectric

    dielectric%nregions = nregions
    IF ( dielectric%nregions .GT. 0 ) &
         & CALL create_environ_functions( dielectric%nregions, 2, epsregion_dim, &
         & epsregion_axis, epsregion_pos, epsregion_spread, epsregion_width, &
         & epsregion_eps, dielectric%regions )

    RETURN

  END SUBROUTINE set_dielectric_regions

  SUBROUTINE init_environ_dielectric_second( cell, dielectric )

    IMPLICIT NONE

    TYPE( environ_cell ), INTENT(IN) :: cell
    TYPE( environ_dielectric ), INTENT(INOUT) :: dielectric

    CHARACTER( LEN=80 ) :: sub_name = 'init_environ_dielectric_second'

    CALL init_environ_density( cell, dielectric%background )
    CALL init_environ_density( cell, dielectric%epsilon )
    CALL init_environ_density( cell, dielectric%depsilon )

    IF ( dielectric%need_gradient ) CALL init_environ_gradient( cell, dielectric%gradient )
    IF ( dielectric%need_factsqrt ) CALL init_environ_density( cell, dielectric%factsqrt )
    IF ( dielectric%need_gradlog ) CALL init_environ_gradient( cell, dielectric%gradlog )

    RETURN

  END SUBROUTINE init_environ_dielectric_second

  SUBROUTINE update_environ_dielectric( dielectric )
  !
  IMPLICIT NONE
  !
  TYPE( environ_dielectric ), INTENT(INOUT) :: dielectric
  !
  ! ... Local variables
  !
  INTEGER :: i
  TYPE( environ_density ) :: local
  TYPE( environ_cell ), POINTER :: cell
  !
  CALL start_clock( 'dielectric' )
  !
  cell => dielectric % epsilon % cell
  !
  IF ( cell % update ) THEN
     !
     ! ... Cells has changed, may need to update the background
     !
     dielectric % background % of_r(:) = dielectric % constant
     !
     IF ( dielectric % nregions .GT. 0 ) THEN
        !
        ! ... Recompute background dielectric and its derivative
        !
        CALL init_environ_density( cell , local )
        DO i = 1, dielectric % nregions
           CALL density_of_functions( 1, dielectric%regions(i), local )
           dielectric % background % of_r(:) = dielectric % background % of_r(:) - &
                  ( dielectric % background % of_r(:) - dielectric % regions(i) % volume ) * local % of_r(:)
        ENDDO
        CALL destroy_environ_density( local )
        !
        IF ( verbose .GE. 3 ) CALL print_environ_density( dielectric%background )
        !
     ENDIF
     !
     ! ... Background has changed, need to update the dielectric when ready
     !
     dielectric % update = .TRUE.
     !
  ENDIF
  !
  ! ... Check if the boundary is under update (status = 1) or has been fully updated (status = 2)
  !
  IF ( dielectric % boundary % update_status .GT. 0 ) dielectric % update = .TRUE.
  !
  IF ( dielectric % update ) THEN
     !
     ! ... Update the dielectric in space and its derivatives if the boundary is ready
     !
     IF ( dielectric % boundary % update_status .EQ. 2 ) THEN
        !
        CALL epsilon_of_boundary( dielectric%boundary, dielectric%background, dielectric%epsilon, dielectric%depsilon )
        !
        IF ( dielectric%need_gradient ) CALL generate_epsilon_gradient( dielectric%epsilon, dielectric%gradient, &
             &dielectric%background, dielectric%boundary )
        IF ( dielectric%need_factsqrt ) CALL generate_epsilon_factsqrt( dielectric%epsilon, dielectric%factsqrt, &
             & dielectric%background, dielectric%boundary )
        IF ( dielectric%need_gradlog  ) CALL generate_epsilon_gradlog( dielectric%epsilon, dielectric%gradlog, &
             & dielectric%background, dielectric%boundary )
        dielectric % update = .FALSE.
        !
     ENDIF
     !
  END IF
  !
  CALL stop_clock( 'dielectric' )
  !
  RETURN
  !
  END SUBROUTINE update_environ_dielectric

  SUBROUTINE epsilon_of_boundary( boundary, background, epsilon, depsilon )

    IMPLICIT NONE

    TYPE(environ_boundary), INTENT(IN) :: boundary
    TYPE(environ_density), INTENT(IN) :: background
    TYPE(environ_density), INTENT(INOUT) :: epsilon, depsilon

    CHARACTER ( LEN=80 ) :: sub_name = 'epsilon_of_boundary'

    SELECT CASE ( boundary%type )

    CASE ( 0 )

       epsilon % of_r = 1.D0 + ( background % of_r - 1.D0 ) * boundary % scaled % of_r
       depsilon % of_r = ( background % of_r - 1.D0 ) * boundary % dscaled % of_r

    CASE ( 1 )

       epsilon % of_r = EXP( LOG( background % of_r ) * ( 1.D0 - boundary % scaled % of_r ) )
       depsilon % of_r = - epsilon % of_r * LOG( background % of_r ) * boundary % dscaled % of_r

    CASE DEFAULT

       CALL errore(sub_name,'Unkown boundary type',1)

    END SELECT

    RETURN

  END SUBROUTINE epsilon_of_boundary

  SUBROUTINE destroy_environ_dielectric(lflag,dielectric)

    IMPLICIT NONE

    LOGICAL, INTENT(IN) :: lflag
    TYPE( environ_dielectric ), INTENT(INOUT) :: dielectric
    CHARACTER (LEN=80) :: sub_name = 'destroy_environ_dielectric'

    IF ( lflag ) THEN

       ! These components were allocated first, destroy only if lflag = .TRUE.

       IF ( dielectric%nregions .GT. 0 ) THEN
          CALL destroy_environ_functions( dielectric%nregions, dielectric%regions )
       ELSE
          IF ( ALLOCATED(dielectric%regions) ) &
               & CALL errore(sub_name,'Found unexpected allocated object',1)
       END IF

       IF (.NOT.ASSOCIATED(dielectric%boundary)) &
            & CALL errore(sub_name,'Trying to destroy a non associated object',1)
       NULLIFY( dielectric%boundary )

    END IF

    CALL destroy_environ_density( dielectric%background )
    CALL destroy_environ_density( dielectric%epsilon )
    CALL destroy_environ_density( dielectric%depsilon )

    IF ( dielectric%need_gradient ) CALL destroy_environ_gradient( dielectric%gradient )
    IF ( dielectric%need_factsqrt ) CALL destroy_environ_density( dielectric%factsqrt )
    IF ( dielectric%need_gradlog ) CALL destroy_environ_gradient( dielectric%gradlog )

    RETURN

  END SUBROUTINE destroy_environ_dielectric

!--------------------------------------------------------------------
  SUBROUTINE generate_epsilon_gradient( epsilon, gradient, &
       & background, boundary )
!--------------------------------------------------------------------

    USE fd_gradient, ONLY: calc_fd_gradient
    USE functions, ONLY: gradient_of_functions

    IMPLICIT NONE

    TYPE( environ_density ), INTENT(IN) :: epsilon
    TYPE( environ_gradient ), INTENT(INOUT) :: gradient
    TYPE( environ_density ), OPTIONAL, INTENT(IN) :: background
    TYPE( environ_boundary ), OPTIONAL, INTENT(IN) :: boundary

    INTEGER, POINTER :: nnr
    TYPE( environ_cell ), POINTER :: cell

    INTEGER :: i
    CHARACTER( LEN=80 ) :: sub_name = 'generate_epsilon_gradient'

    nnr => epsilon%cell%nnr
    cell => epsilon%cell

    SELECT CASE ( dielectric_core )
       !
    CASE ( 'fft' )
       !
       CALL external_gradient(epsilon%of_r, gradient%of_r)
       !
    CASE ( 'fd' )
       !
       CALL calc_fd_gradient(nfdpoint, icfd, ncfd, nnr, epsilon%of_r, gradient%of_r)
       !
    CASE ( 'analytic' )
       !
       IF ( boundary%mode .EQ. 'ions' ) THEN
          CALL gradient_of_functions(boundary%ions%number, boundary%soft_spheres, gradient)
       ELSE
          CALL external_gradient(boundary%density%of_r,gradient%of_r)
          DO i = 1, 3
             gradient%of_r(i,:) = boundary%dscaled%of_r(:) * gradient%of_r(i,:)
          ENDDO
       ENDIF
       !
    CASE DEFAULT
       !
    END SELECT

    RETURN
!--------------------------------------------------------------------
  END SUBROUTINE generate_epsilon_gradient
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE generate_epsilon_gradlog( epsilon, gradlog, &
       & background, boundary )
!--------------------------------------------------------------------

    USE fd_gradient, ONLY : calc_fd_gradient
    USE functions, ONLY: gradient_of_functions

    IMPLICIT NONE

    TYPE( environ_density ), INTENT(IN) :: epsilon
    TYPE( environ_gradient ), INTENT(INOUT) :: gradlog
    TYPE( environ_density ), OPTIONAL, INTENT(IN) :: background
    TYPE( environ_boundary ), OPTIONAL, INTENT(IN) :: boundary

    INTEGER, POINTER :: nnr
    TYPE( environ_cell ), POINTER :: cell

    INTEGER :: i
    TYPE( environ_density ) :: log_epsilon
    CHARACTER( LEN=80 ) :: label
    CHARACTER( LEN=80 ) :: sub_name = 'generate_epsilon_gradlog'

    nnr => epsilon%cell%nnr
    cell => epsilon%cell

    SELECT CASE ( dielectric_core )
       !
    CASE ( 'fft' )
       !
       CALL init_environ_density( cell, log_epsilon )
       log_epsilon%of_r = LOG( epsilon%of_r )
       !
       CALL external_gradient(log_epsilon%of_r, gradlog%of_r)
       !
       CALL destroy_environ_density( log_epsilon )
       !
    CASE ( 'fd' )
       !
       label='logeps'
       CALL create_environ_density( log_epsilon, label )
       CALL init_environ_density( cell, log_epsilon )
       log_epsilon%of_r = LOG( epsilon%of_r )
       !
       CALL calc_fd_gradient(nfdpoint, icfd, ncfd, nnr, log_epsilon%of_r, gradlog%of_r)
       CALL update_gradient_modulus(gradlog)
       !
       CALL destroy_environ_density( log_epsilon )
       !
    CASE ( 'analytic' )
       !
       IF ( boundary%mode .EQ. 'ions' ) THEN
          CALL gradient_of_functions(boundary%ions%number, boundary%soft_spheres, gradlog)
       ELSE
          CALL external_gradient(boundary%density%of_r,gradlog%of_r)
          DO i = 1, 3
             gradlog%of_r(i,:) = boundary%dscaled%of_r(:) * gradlog%of_r(i,:)
          ENDDO
       ENDIF
       !
       DO i = 1, 3
          gradlog%of_r(i,:) = gradlog%of_r(i,:) / epsilon%of_r(:)
       ENDDO
       !
    CASE DEFAULT
       !
    END SELECT

    RETURN
!--------------------------------------------------------------------
  END SUBROUTINE generate_epsilon_gradlog
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE generate_epsilon_factsqrt( epsilon, factsqrt, &
       & background, boundary )
!--------------------------------------------------------------------

    USE functions, ONLY: gradient_of_functions, laplacian_of_functions

    IMPLICIT NONE

    TYPE( environ_density ), INTENT(IN) :: epsilon
    TYPE( environ_density ), INTENT(INOUT) :: factsqrt
    TYPE( environ_density ), OPTIONAL, INTENT(IN) :: background
    TYPE( environ_boundary ), OPTIONAL, INTENT(IN) :: boundary

    INTEGER, POINTER :: nnr
    TYPE( environ_cell ), POINTER :: cell
    TYPE( environ_density ) :: laplacian
    TYPE( environ_gradient ) :: gradient
    CHARACTER( LEN=80 ) :: sub_name = 'generate_epsilon_factsqrt'

    nnr => epsilon%cell%nnr
    cell => epsilon%cell

    SELECT CASE ( dielectric_core )
       !
    CASE ( 'fft' )
       !
       CALL errore(sub_name,'Option not yet implemented',1)
       !
    CASE ( 'fd' )
       !
       CALL errore(sub_name,'Option not yet implemented',1)
       !
    CASE ( 'analytic' )
       !
       CALL init_environ_density( cell, laplacian )
       CALL init_environ_gradient( cell, gradient )
       !
       IF ( boundary%mode .EQ. 'ions' ) THEN
          !
          CALL gradient_of_functions(boundary%ions%number, boundary%soft_spheres, gradient)
          CALL laplacian_of_functions(boundary%ions%number, boundary%soft_spheres, laplacian)
          CALL update_gradient_modulus(gradient)
          !
          factsqrt%of_r(:) = laplacian%of_r(:) - 0.5D0 * gradient%modulus%of_r(:) / epsilon%of_r(:)
          !
       ELSE
          !
          CALL external_gradient(boundary%density%of_r,gradient%of_r)
          CALL external_laplacian(boundary%density%of_r,laplacian%of_r)
          CALL update_gradient_modulus(gradient)
          !
          factsqrt%of_r(:) = boundary%dscaled%of_r(:)*laplacian%of_r(:) + &
               & gradient%modulus%of_r(:) * &
               & (boundary%d2scaled%of_r(:)-0.5D0*boundary%dscaled%of_r(:)**2/epsilon%of_r(:))
          !
       ENDIF
       factsqrt%of_r = factsqrt%of_r * 0.5D0 / e2 / fpi
       !
       CALL destroy_environ_gradient(gradient)
       CALL destroy_environ_density(laplacian)
       !
    CASE DEFAULT
       !
    END SELECT

    RETURN
!--------------------------------------------------------------------
  END SUBROUTINE generate_epsilon_factsqrt
!--------------------------------------------------------------------
!
END MODULE dielectric
