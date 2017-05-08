MODULE boundary

  USE environ_types

  IMPLICIT NONE

  SAVE

  PUBLIC :: update_environ_boundary

CONTAINS

  SUBROUTINE update_environ_boundary( bound )

    USE functions, ONLY : density_of_functions
    USE generate_boundary, ONLY : boundary_of_density

    IMPLICIT NONE

    TYPE( environ_boundary ), INTENT(INOUT) :: bound

    CHARACTER( LEN=80 ) :: sub_name = 'update_environ_boundary'

    IF ( .NOT. bound % ions % update .AND. .NOT. bound % electrons % update ) THEN
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

    CASE ( 'electrons' )

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

    CASE ( 'ions' )

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

END MODULE boundary
