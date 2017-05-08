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
  USE environ_base,       ONLY: verbose, environ_unit, ions, e2
  USE electrostatic_base, ONLY : dielectric_core, nfdpoint, icfd, ncfd
  USE environ_debug,      ONLY: write_cube
  USE functions,          ONLY: density_of_functions
  !
  IMPLICIT NONE
  !
  SAVE
  !
  PRIVATE
  !
  PUBLIC :: update_environ_dielectric
  !
CONTAINS
  !
!--------------------------------------------------------------------
  SUBROUTINE update_environ_dielectric( dielectric )
!--------------------------------------------------------------------
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
        CALL create_environ_density( local )
        CALL init_environ_density( cell , local )
        DO i = 1, dielectric % nregions
           CALL density_of_functions( 1, dielectric%regions(i), local )
           dielectric % background % of_r(:) = dielectric % background % of_r(:) - &
                  ( dielectric % background % of_r(:) - dielectric % regions(i) % volume ) * local % of_r(:)
        ENDDO
        CALL destroy_environ_density( local )
        !
        IF ( verbose .GE. 3 ) CALL write_cube( ions, dielectric%background )
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
!        dielectric % epsilon % of_r = 1.D0 + ( dielectric % background % of_r - 1.D0 ) * &
!             & dielectric % boundary % scaled % of_r
        dielectric % epsilon % of_r = dielectric % boundary % scaled % of_r
        !
        IF ( dielectric%need_gradient ) CALL generate_epsilon_gradient( dielectric%nregions, &
             & dielectric%epsilon, dielectric%gradient, dielectric%regions, dielectric%background, dielectric%boundary )
        IF ( dielectric%need_factsqrt ) CALL generate_epsilon_factsqrt( dielectric%nregions, &
             & dielectric%epsilon, dielectric%factsqrt, dielectric%regions, dielectric%background, dielectric%boundary )
        IF ( dielectric%need_gradlog  ) CALL generate_epsilon_gradlog( dielectric%nregions, &
             & dielectric%epsilon, dielectric%gradlog, dielectric%regions, dielectric%background, dielectric%boundary )
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
!--------------------------------------------------------------------
  END SUBROUTINE update_environ_dielectric
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE generate_epsilon_gradient( nregions, epsilon, gradient, &
       & regions, background, boundary )
!--------------------------------------------------------------------

    USE electrostatic_base, ONLY : dielectric_core
    USE fd_gradient, ONLY: calc_fd_gradient
    USE functions, ONLY: gradient_of_functions

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nregions
    TYPE( environ_density ), INTENT(IN) :: epsilon
    TYPE( environ_gradient ), INTENT(OUT) :: gradient
    TYPE( environ_functions ), OPTIONAL, INTENT(IN) :: regions(nregions)
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
!       CALL init_environ_gradient(cell,gbackground)
       IF ( nregions .GT. 0 ) &
            CALL errore(sub_name,'Option not yet implemented',1)
!            & CALL gradient_of_functions(nregions,regions,gbackground)
       !
       IF ( boundary%mode .EQ. 'ions' ) THEN
          CALL gradient_of_functions(boundary%ions%number, boundary%soft_spheres, gradient)
       ELSE
          CALL external_gradient(boundary%density%of_r,gradient%of_r)
          DO i = 1, 3
             gradient%of_r(i,:) = boundary%dscaled%of_r(:) * gradient%of_r(i,:)
          ENDDO
       ENDIF
!       CALL destroy_environ_gradient(gbackground)
       !
    CASE DEFAULT
       !
    END SELECT

    RETURN
!--------------------------------------------------------------------
  END SUBROUTINE generate_epsilon_gradient
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE generate_epsilon_gradlog( nregions, epsilon, gradlog, &
       & regions, background, boundary )
!--------------------------------------------------------------------

    USE fd_gradient, ONLY : calc_fd_gradient
    USE functions, ONLY: gradient_of_functions

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nregions
    TYPE( environ_density ), INTENT(IN) :: epsilon
    TYPE( environ_gradient ), INTENT(OUT) :: gradlog
    TYPE( environ_functions ), OPTIONAL, INTENT(IN) :: regions(nregions)
    TYPE( environ_density ), OPTIONAL, INTENT(IN) :: background
    TYPE( environ_boundary ), OPTIONAL, INTENT(IN) :: boundary

    INTEGER, POINTER :: nnr
    TYPE( environ_cell ), POINTER :: cell

    INTEGER :: i
    TYPE( environ_density ) :: log_epsilon
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
       CALL init_environ_density( cell, log_epsilon )
       log_epsilon%of_r = LOG( epsilon%of_r )
       !
       CALL calc_fd_gradient(nfdpoint, icfd, ncfd, nnr, log_epsilon%of_r, gradlog%of_r)
       !
       CALL destroy_environ_density( log_epsilon )
       !
    CASE ( 'analytic' )
       !
!       CALL init_environ_gradient(cell,gbackground)
       IF ( nregions .GT. 0 ) &
            CALL errore(sub_name,'Option not yet implemented',1)
!            & CALL gradient_of_functions(nregions,regions,gbackground)
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
       !       CALL destroy_environ_gradient(gbackground)
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
  SUBROUTINE generate_epsilon_factsqrt( nregions, epsilon, factsqrt, &
       & regions, background, boundary )
!--------------------------------------------------------------------

    USE functions, ONLY: gradient_of_functions, laplacian_of_functions

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nregions
    TYPE( environ_density ), INTENT(IN) :: epsilon
    TYPE( environ_density ), INTENT(INOUT) :: factsqrt
    TYPE( environ_functions ), OPTIONAL, INTENT(IN) :: regions(nregions)
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
!       CALL init_environ_gradient(cell,gbackground)
       IF ( nregions .GT. 0 ) &
            CALL errore(sub_name,'Option not yet implemented',1)
!            & CALL gradient_of_functions(nregions,regions,gbackground)
       !
       IF ( boundary%mode .EQ. 'ions' ) THEN
          !
          CALL gradient_of_functions(boundary%ions%number, boundary%soft_spheres, gradient)
          CALL laplacian_of_functions(boundary%ions%number, boundary%soft_spheres, laplacian)
          CALL update_gradient_modulus(gradient)
          !
          factsqrt%of_r(:) = laplacian%of_r(:) - 0.5D0 * gradient%modulus(:) / epsilon%of_r(:)
          !
       ELSE
          !
          CALL external_gradient(boundary%density%of_r,gradient%of_r)
          CALL external_laplacian(boundary%density%of_r,laplacian%of_r)
          CALL update_gradient_modulus(gradient)
          !
          factsqrt%of_r(:) = boundary%dscaled%of_r(:)*laplacian%of_r(:) + &
               & gradient%modulus(:) * &
               & (boundary%d2scaled%of_r(:)-0.5D0*boundary%dscaled%of_r(:)**2/epsilon%of_r(:))
          !
       ENDIF
       factsqrt%of_r = factsqrt%of_r * 0.5D0 / e2 / fpi
       !
       CALL destroy_environ_gradient(gradient)
       CALL destroy_environ_density(laplacian)
!       CALL destroy_environ_gradient(gbackground)
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
