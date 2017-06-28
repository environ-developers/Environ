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
        CALL dielectric_of_boundary( dielectric )
        !
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

  SUBROUTINE dielectric_of_boundary( dielectric )

    USE fd_gradient, ONLY : calc_fd_gradient

    IMPLICIT NONE

    TYPE(environ_dielectric), TARGET, INTENT(INOUT) :: dielectric

    INTEGER :: ipol
    INTEGER, POINTER :: nnr
    REAL( DP ), DIMENSION( : ), POINTER :: factsqrteps, eps, deps, const
    REAL( DP ), DIMENSION( : ), POINTER :: scaled, gradscaledmod, laplscaled
    REAL( DP ), DIMENSION( :, : ), POINTER :: gradeps, gradlogeps, gradscaled

    REAL( DP ), DIMENSION( : ), ALLOCATABLE :: dlogeps
    REAL( DP ), DIMENSION( : ), ALLOCATABLE :: d2eps

    CHARACTER ( LEN=80 ) :: sub_name = 'epsilon_of_boundary'

    ! Aliases and sanity checks

    nnr => dielectric%epsilon%cell%nnr
    eps => dielectric % epsilon % of_r
    deps => dielectric % depsilon % of_r
    const => dielectric % background % of_r
    scaled => dielectric % boundary % scaled % of_r
    IF ( dielectric % need_gradient ) THEN
       IF ( .NOT. ALLOCATED( dielectric % boundary % gradient % of_r ) ) &
            & CALL errore(sub_name,'Missing required gradient of boundary',1)
       gradscaled => dielectric % boundary % gradient % of_r
       gradeps => dielectric % gradient % of_r
    ENDIF
    IF ( dielectric % need_gradlog ) THEN
       IF ( .NOT. ALLOCATED( dielectric % boundary % gradient % of_r ) ) &
            & CALL errore(sub_name,'Missing required gradient of boundary',1)
       gradscaled => dielectric % boundary % gradient % of_r
       gradlogeps => dielectric % gradlog % of_r
       ALLOCATE( dlogeps( nnr ) )
    ENDIF
    IF ( dielectric % need_factsqrt ) THEN
       IF ( .NOT. ALLOCATED( dielectric % boundary % gradient % of_r ) ) &
            & CALL errore(sub_name,'Missing required gradient of boundary',1)
       gradscaledmod => dielectric % boundary % gradient % modulus % of_r
       IF ( .NOT. ALLOCATED( dielectric % boundary % laplacian % of_r ) ) &
            & CALL errore(sub_name,'Missing required laplacian of boundary',1)
       laplscaled => dielectric % boundary % laplacian % of_r
       factsqrteps => dielectric % factsqrt % of_r
       ALLOCATE( d2eps( nnr ) )
    ENDIF

    ! Compute epsilon(r) and its derivative wrt boundary

    SELECT CASE ( dielectric%boundary%type )

    CASE ( 0, 2 )

       eps = 1.D0 + ( const - 1.D0 ) * ( 1.D0 - scaled )
       deps = ( 1.D0 - const )
       IF ( dielectric % need_gradlog ) dlogeps = deps / eps
       IF ( dielectric % need_factsqrt ) d2eps = 0.D0

    CASE ( 1 )

       eps = EXP( LOG( const ) * ( 1.D0 - scaled ) )
       deps = - eps * LOG( const )
       IF ( dielectric % need_gradlog ) dlogeps = - LOG( const )
       IF ( dielectric % need_factsqrt ) d2eps = eps * LOG( const )**2

    CASE DEFAULT

       CALL errore(sub_name,'Unkown boundary type',1)

    END SELECT

    ! If needed, compute derived quantites

    IF ( dielectric % need_gradient ) THEN
       DO ipol = 1, 3
          gradeps( ipol, : ) = deps( : ) * gradscaled( ipol, : )
       ENDDO
       CALL update_gradient_modulus( dielectric%gradient )
    END IF
    IF ( dielectric % need_gradlog ) THEN
       DO ipol = 1, 3
          gradlogeps( ipol, : ) = dlogeps( : ) * gradscaled( ipol, : )
       ENDDO
       CALL update_gradient_modulus( dielectric%gradlog )
       DEALLOCATE( dlogeps )
    ENDIF
    IF ( dielectric % need_factsqrt ) THEN
       factsqrteps(:) = 0.5D0 * ( ( d2eps( : ) - deps(:)**2 * 0.5D0 / eps(:) ) * gradscaledmod(:) + &
            & deps(:) * laplscaled(:) )
       DEALLOCATE( d2eps )
    ENDIF

    RETURN

  END SUBROUTINE dielectric_of_boundary

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

END MODULE dielectric
