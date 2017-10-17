MODULE boundary

  USE environ_types
  USE functions

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: create_environ_boundary, init_environ_boundary_first, &
       & init_environ_boundary_second, set_soft_spheres, &
       & update_environ_boundary, destroy_environ_boundary

CONTAINS

  SUBROUTINE create_environ_boundary(boundary)

    IMPLICIT NONE

    TYPE( environ_boundary ), INTENT(INOUT) :: boundary

    CHARACTER( LEN=80 ) :: sub_name = 'create_environ_boundary'
    CHARACTER( LEN=80 ) :: label = ' '

    boundary%update_status = 0

    label = 'boundary'
    CALL create_environ_density( boundary%scaled, label )
    boundary%volume = 0.D0
    boundary%surface = 0.D0

    boundary%need_electrons = .FALSE.
    NULLIFY( boundary%electrons )
    boundary%need_ions = .FALSE.
    NULLIFY( boundary%ions   )

    ! Optional components

    boundary%deriv = 0
    label = 'gradboundary'
    CALL create_environ_gradient( boundary%gradient, label )
    label = 'laplboundary'
    CALL create_environ_density( boundary%laplacian, label )
    label = 'dsurface'
    CALL create_environ_density( boundary%dsurface, label )

    ! Components required for boundary of density

    label = 'density'
    CALL create_environ_density( boundary%density, label )
    label = 'dboundary'
    CALL create_environ_density( boundary%dscaled, label )
    label = 'd2boundary'
    CALL create_environ_density( boundary%d2scaled, label )

    ! Components required for boundary of functions

    IF ( ALLOCATED( boundary%soft_spheres ) ) &
         & CALL errore(sub_name,'Trying to create an already allocated object',1)

    ! Components required for solvent-aware interface

    boundary%solvent_aware = .FALSE.
    label = 'local'
    CALL create_environ_density( boundary%local, label )
    label = 'probe'
    CALL create_environ_density( boundary%probe, label )
    label = 'emptying'
    CALL create_environ_density( boundary%emptying, label )
    label = 'demptying'
    CALL create_environ_density( boundary%emptying, label )

    RETURN

  END SUBROUTINE create_environ_boundary

  SUBROUTINE init_environ_boundary_first( need_gradient, need_laplacian, &
       & need_hessian, mode, stype, rhomax, rhomin, tbeta, const, alpha, &
       & softness, solvent_radius, radial_scale, radial_spread,          &
       & emptying_threshold, emptying_spread, electrons, ions, boundary )

    IMPLICIT NONE

    CHARACTER( LEN=80 ), INTENT(IN) :: mode
    INTEGER, INTENT(IN) :: stype
    REAL( DP ), INTENT(IN) :: rhomax, rhomin, tbeta, const
    LOGICAL, INTENT(IN) :: need_gradient, need_laplacian, need_hessian
    REAL( DP ), INTENT(IN) :: alpha
    REAL( DP ), INTENT(IN) :: softness
    REAL( DP ), INTENT(IN) :: solvent_radius
    REAL( DP ), INTENT(IN) :: radial_scale, radial_spread
    REAL( DP ), INTENT(IN) :: emptying_threshold, emptying_spread
    TYPE( environ_electrons ), TARGET, INTENT(IN) :: electrons
    TYPE( environ_ions ), TARGET, INTENT(IN) :: ions
    TYPE( environ_boundary ), INTENT(INOUT) :: boundary

    IF ( need_hessian ) THEN
       boundary%deriv = 3
    ELSE IF ( need_laplacian ) THEN
       boundary%deriv = 2
    ELSE IF ( need_gradient ) THEN
       boundary%deriv = 1
    ENDIF

    boundary%mode = mode

    boundary%need_electrons = ( mode .EQ. 'electronic' ) .OR. ( mode .EQ. 'full' )
    IF ( boundary%need_electrons ) boundary%electrons => electrons
    boundary%need_ions = ( mode .EQ. 'ionic' ) .OR. ( mode .EQ. 'full' )
    IF ( boundary%need_ions ) boundary%ions => ions

    boundary%type = stype
    boundary%rhomax = rhomax
    boundary%rhomin = rhomin
    boundary%fact = LOG( rhomax / rhomin )
    boundary%rhozero = ( rhomax + rhomin ) * 0.5_DP
    boundary%tbeta = tbeta
    boundary%deltarho = rhomax - rhomin

    boundary%const = const
    IF ( const .EQ. 1.D0 ) boundary%const = 2.D0

    boundary%alpha = alpha
    boundary%softness = softness
    IF ( boundary%need_ions .AND. .NOT. boundary%need_electrons ) &
         & ALLOCATE( boundary%soft_spheres( boundary%ions%number ) )

    boundary%solvent_aware = solvent_radius .GT. 0.D0

    boundary%solvent_probe%type = 4
    boundary%solvent_probe%pos = 0.D0
    boundary%solvent_probe%volume = 1.D0
    boundary%solvent_probe%dim = 0
    boundary%solvent_probe%axis = 1

    boundary%solvent_probe%spread = radial_spread
    boundary%solvent_probe%width = solvent_radius * radial_scale

    boundary%emptying_threshold = emptying_threshold
    boundary%emptying_spread = emptying_spread

    RETURN

  END SUBROUTINE init_environ_boundary_first

  SUBROUTINE init_environ_boundary_second( cell, boundary )

    IMPLICIT NONE

    TYPE( environ_cell ), INTENT(IN) :: cell
    TYPE( environ_boundary ), INTENT(INOUT) :: boundary

    CALL init_environ_density( cell, boundary%scaled )

    IF ( boundary%mode .NE. 'ionic' ) THEN
       CALL init_environ_density( cell, boundary%density )
       CALL init_environ_density( cell, boundary%dscaled )
       CALL init_environ_density( cell, boundary%d2scaled )
    END IF
    IF ( boundary%deriv .GE. 1 ) CALL init_environ_gradient( cell, boundary%gradient )
    IF ( boundary%deriv .GE. 2 ) CALL init_environ_density( cell, boundary%laplacian )
    IF ( boundary%deriv .GE. 3 ) CALL init_environ_density( cell, boundary%dsurface )

    IF ( boundary%solvent_aware ) THEN
       CALL init_environ_density( cell, boundary%local )
       CALL init_environ_density( cell, boundary%probe )
       CALL init_environ_density( cell, boundary%emptying )
       CALL init_environ_density( cell, boundary%demptying )
    ENDIF

    RETURN

  END SUBROUTINE init_environ_boundary_second

  SUBROUTINE set_soft_spheres( boundary )

    IMPLICIT NONE

    TYPE( environ_boundary ), INTENT(INOUT) :: boundary

    INTEGER :: i
    REAL( DP ) :: radius

    IF ( boundary % mode .NE. 'ionic' ) RETURN

    DO i = 1, boundary%ions%number
       radius = boundary%ions%iontype(boundary%ions%ityp(i))%solvationrad * boundary%alpha
       boundary%soft_spheres(i) = environ_functions(5,1,0,radius,boundary%softness,1.D0,&
            & boundary%ions%tau(:,i))
    ENDDO

    RETURN

  END SUBROUTINE set_soft_spheres

  SUBROUTINE update_environ_boundary( bound )

    USE generate_boundary, ONLY : boundary_of_density, boundary_of_functions

    IMPLICIT NONE

    TYPE( environ_boundary ), INTENT(INOUT) :: bound

    LOGICAL :: update_anything
    CHARACTER( LEN=80 ) :: sub_name = 'update_environ_boundary'

    INTEGER :: i
    TYPE( environ_cell ), POINTER :: cell
    TYPE( environ_density ) :: local

    cell => bound%density%cell

    update_anything = .FALSE.
    IF ( bound % need_ions ) update_anything = bound % ions % update
    IF ( bound % need_electrons ) update_anything = update_anything .OR. bound % electrons % update
    IF ( .NOT. update_anything ) THEN
       !
       ! ... Nothing is under update, change update_status and exit
       !
       IF ( bound % update_status .EQ. 2 ) bound % update_status = 0
       RETURN
       !
    ENDIF

    SELECT CASE ( bound % mode )

    CASE ( 'full' )

       IF ( bound % ions % update ) THEN
          !
          ! ... Compute the ionic part
          !
          CALL density_of_functions( bound%ions%number, bound%ions%core_electrons, bound%ions%core )
          !
          bound % update_status = 1 ! waiting to finish update
          !
       ENDIF

       IF ( bound % electrons % update ) THEN
          !
          ! ... Check if the ionic part has been updated
          !
          IF ( bound % update_status .EQ. 0 ) &
               & CALL errore(sub_name,'Wrong update status, possibly missing ionic update',1)
          !
          bound % density % of_r = bound % electrons % density % of_r + bound % ions % core % of_r
          !
          CALL boundary_of_density( bound % density, bound )
          !
          bound % update_status = 2 ! boundary has changed and is ready
          !
       ENDIF

    CASE ( 'electronic' )

       IF ( bound % electrons % update ) THEN
          !
          bound % density % of_r = bound % electrons % density % of_r
          !
          CALL boundary_of_density( bound % density, bound )
          !
          bound % update_status = 2 ! boundary has changes and is ready
          !
       ELSE
          !
          IF ( bound % update_status .EQ. 2 ) bound % update_status = 0 ! boundary has not changed
          RETURN
          !
       ENDIF

    CASE ( 'ionic' )

       IF ( bound % ions % update ) THEN
          !
          ! ... Only ions are needed, fully update the boundary
          !
          CALL boundary_of_functions( bound%ions%number, bound%soft_spheres, bound )
          !
          bound % update_status = 2 ! boundary has changed and is ready
          !
          RETURN
          !
       ELSE
          !
          IF ( bound % update_status .EQ. 2 ) bound % update_status = 0 ! boundary has not changed
          RETURN
          !
       ENDIF

    CASE DEFAULT

       CALL errore(sub_name,'Unrecognized boundary mode',1)

    END SELECT

    RETURN

  END SUBROUTINE update_environ_boundary

  SUBROUTINE destroy_environ_boundary(lflag, boundary)

    IMPLICIT NONE

    LOGICAL, INTENT(IN) :: lflag
    TYPE( environ_boundary ), INTENT(INOUT) :: boundary
    CHARACTER (LEN=80) :: sub_name = 'destroy_environ_boundary'

    IF ( lflag ) THEN

       ! These components were allocated first, destroy only if lflag = .TRUE.

       IF ( boundary%need_ions ) THEN
          IF ( .NOT. boundary%need_electrons ) &
               & CALL destroy_environ_functions( boundary%ions%number, boundary%soft_spheres )
          IF (.NOT.ASSOCIATED(boundary%ions)) &
               & CALL errore(sub_name,'Trying to destroy a non associated object',1)
          NULLIFY(boundary%ions)
       ELSE
          IF (ASSOCIATED(boundary%ions))&
               & CALL errore(sub_name,'Found an unexpected associated object',1)
       ENDIF

       IF ( boundary%need_electrons ) THEN
          IF (ASSOCIATED(boundary%electrons)) NULLIFY(boundary%electrons)
       ENDIF

    ENDIF

    CALL destroy_environ_density( boundary%scaled )
    IF ( boundary%mode .NE. 'ionic' ) THEN
       CALL destroy_environ_density( boundary%density )
       CALL destroy_environ_density( boundary%dscaled )
       CALL destroy_environ_density( boundary%d2scaled )
    ENDIF
    IF ( boundary%deriv .GE. 1 ) CALL destroy_environ_gradient( boundary%gradient )
    IF ( boundary%deriv .GE. 2 ) CALL destroy_environ_density( boundary%laplacian )
    IF ( boundary%deriv .GE. 3 ) CALL destroy_environ_density( boundary%dsurface )

    IF ( boundary%solvent_aware ) THEN
       CALL destroy_environ_density( boundary%local )
       CALL destroy_environ_density( boundary%probe )
       CALL destroy_environ_density( boundary%emptying )
       CALL destroy_environ_density( boundary%demptying )
    ENDIF

    RETURN

  END SUBROUTINE destroy_environ_boundary

END MODULE boundary
