MODULE boundary

  USE environ_types
  USE functions

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: create_environ_boundary, init_environ_boundary_first, &
       & init_environ_boundary_second, update_environ_boundary, destroy_environ_boundary

CONTAINS

  SUBROUTINE create_environ_boundary(boundary)

    IMPLICIT NONE

    TYPE( environ_boundary ), INTENT(INOUT) :: boundary

    CHARACTER( LEN=80 ) :: sub_name = 'create_environ_boundary'
    CHARACTER( LEN=80 ) :: label = ' '

    boundary%update_status = 0

    label = 'boundary'
    CALL create_environ_density( boundary%scaled, label )
    label = 'dboundary'
    CALL create_environ_density( boundary%dscaled, label )
    label = 'd2boundary'
    CALL create_environ_density( boundary%d2scaled, label )

    boundary%need_electrons = .FALSE.
    NULLIFY( boundary%electrons )
    boundary%need_ions = .FALSE.
    NULLIFY( boundary%ions   )
    boundary%need_theta = .FALSE.
    label = 'theta'
    CALL create_environ_density( boundary%theta, label )

    IF ( ALLOCATED( boundary%soft_spheres ) ) &
         & CALL errore(sub_name,'Trying to create an already allocated object',1)

    label = 'density'
    CALL create_environ_density( boundary%density, label )

    RETURN

  END SUBROUTINE create_environ_boundary

  SUBROUTINE init_environ_boundary_first( mode, constant, type, &
       & rhomax, rhomin, tbeta, need_theta, delta, alpha, &
       & softness, electrons, ions, boundary )

    IMPLICIT NONE

    CHARACTER( LEN=80 ), INTENT(IN) :: mode
    REAL( DP ), INTENT(IN) :: constant
    INTEGER, INTENT(IN) :: type
    REAL( DP ), INTENT(IN) :: rhomax, rhomin, tbeta
    LOGICAL, INTENT(IN) :: need_theta
    REAL( DP ), INTENT(IN) :: delta
    REAL( DP ), INTENT(IN) :: alpha
    REAL( DP ), INTENT(IN) :: softness
    TYPE( environ_electrons ), TARGET, INTENT(IN) :: electrons
    TYPE( environ_ions ), TARGET, INTENT(IN) :: ions
    TYPE( environ_boundary ), INTENT(INOUT) :: boundary

    INTEGER :: i
    REAL( DP ) :: radius

    boundary%mode = mode
    boundary%type = type
    boundary%rhomax = rhomax
    boundary%rhomin = rhomin
    boundary%fact = LOG( rhomax / rhomin )
    boundary%rhozero = ( rhomax + rhomin ) * 0.5_DP
    boundary%tbeta = tbeta
    boundary%deltarho = rhomax - rhomin

    IF ( constant .GT. 1.D0 ) THEN
       boundary%constant = constant
       boundary%scaling_factor = -1.D0/(constant-1.D0)
    ELSE
       boundary%constant = 2.D0
       boundary%scaling_factor = -1.D0
    ENDIF

    boundary%need_electrons = ( mode .EQ. 'electronic' ) .OR. ( mode .EQ. 'full' )
    IF ( boundary%need_electrons ) boundary%electrons => electrons

    boundary%need_ions = ( mode .EQ. 'ionic' ) .OR. ( mode .EQ. 'full' )
    IF ( boundary%need_ions ) boundary%ions => ions
    boundary%alpha = alpha
    boundary%softness = softness
    IF ( boundary%need_ions .AND. .NOT. boundary%need_electrons ) THEN
       ALLOCATE( boundary%soft_spheres( boundary%ions%number ) )
       DO i = 1, boundary%ions%number
          radius = boundary%ions%iontype(boundary%ions%ityp(i))%solvationrad * boundary%alpha
          boundary%soft_spheres(i) = environ_functions(2,0,1,radius,1.0_DP,boundary%softness,&
                  & boundary%ions%tau(:,i))
       ENDDO
    ENDIF

    boundary%need_theta = need_theta
    boundary%deltatheta = delta

    RETURN

  END SUBROUTINE init_environ_boundary_first

  SUBROUTINE init_environ_boundary_second( cell, boundary )

    IMPLICIT NONE

    TYPE( environ_cell ), INTENT(IN) :: cell
    TYPE( environ_boundary ), INTENT(INOUT) :: boundary

    CALL init_environ_density( cell, boundary%scaled )
    CALL init_environ_density( cell, boundary%dscaled )
    CALL init_environ_density( cell, boundary%d2scaled )

    IF ( boundary%need_electrons ) CALL init_environ_density( cell, boundary%density )
    IF ( boundary%need_theta) CALL init_environ_density( cell, boundary%theta )

    RETURN

  END SUBROUTINE init_environ_boundary_second

  SUBROUTINE update_environ_boundary( bound )

    USE generate_boundary, ONLY : boundary_of_density

    IMPLICIT NONE

    TYPE( environ_boundary ), INTENT(INOUT) :: bound

    LOGICAL :: update_anything
    CHARACTER( LEN=80 ) :: sub_name = 'update_environ_boundary'

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
          ! COMPUTE ION-DEPENDENT CAVITY
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
          IF (.NOT.ASSOCIATED(boundary%ions)) &
               & CALL errore(sub_name,'Trying to destroy a non associated object',1)
          NULLIFY(boundary%ions)
          IF ( .NOT. boundary%need_electrons ) &
               & CALL destroy_environ_functions( boundary%ions%number, boundary%soft_spheres )
       ELSE
          IF (ASSOCIATED(boundary%ions))&
               & CALL errore(sub_name,'Found an unexpected associated object',1)
       ENDIF

       IF ( boundary%need_electrons ) THEN
          IF (ASSOCIATED(boundary%electrons)) NULLIFY(boundary%electrons)
       ENDIF

    ENDIF

    CALL destroy_environ_density( boundary%scaled )
    CALL destroy_environ_density( boundary%dscaled )
    CALL destroy_environ_density( boundary%d2scaled )

    IF ( boundary%need_electrons ) CALL destroy_environ_density( boundary%density )
    IF ( boundary%need_theta ) CALL destroy_environ_density( boundary%theta )

    RETURN

  END SUBROUTINE destroy_environ_boundary

END MODULE boundary
