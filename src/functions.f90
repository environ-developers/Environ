MODULE functions

  USE environ_types

  PRIVATE

  PUBLIC : density_of_functions

CONTAINS

  SUBROUTINE density_of_functions( nfunctions, functions, density )

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nfunctions
    TYPE( environ_functions ), DIMENSION(nfunctions), INTENT(IN) :: functions
    TYPE( environ_density ), INTENT(INOUT) :: density

    INTEGER :: i

    density%of_r = 0.D0

    DO i = 1, nfunctions
       ASSOCIATE( type   => functions(i)%type, &
                  pos    => functions(i)%pos, &
                  spread => functions(i)%spread, &
                  charge => functions(i)%volume, &
                  width  => functions(i)%width, &
                  dim    => functions(i)%dim, &
                  axis   => functions(i)%axis, &
                  nnr    => density%cell%nnr )
         SELECT CASE ( type )
         CASE ( 1 ) ! Gaussian
            CALL generate_gaussian(nnr, dim, axis, charge, spread, pos, density%of_r)
         CASE ( 2 ) ! Erfc
            CALL generate_erfc(nnr, dim, axis, charge, width, spread, pos, density%of_r)
         CASE ( 3 ) ! Exponential
            CALL generate_exponential(nnr, spread, pos, density%of_r)
         END SELECT
       END ASSOCIATE
    END DO

    RETURN

  END SUBROUTINE density_of_functions

  SUBROUTINE gradient_of_functions( nfunctions, functions, gradient )

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nfunctions
    TYPE( environ_functions ), DIMENSION(nfunctions), INTENT(IN) :: functions
    TYPE( environ_gradient ), INTENT(INOUT) :: gradient

    INTEGER :: i

    gradient%of_r = 0.D0

    DO i = 1, nfunctions
       ASSOCIATE( type   => functions(i)%type, &
                  pos    => functions(i)%pos, &
                  spread => functions(i)%spread, &
                  charge => functions(i)%volume, &
                  width  => functions(i)%width, &
                  dim    => functions(i)%dim, &
                  axis   => functions(i)%axis, &
                  nnr    => gradient%cell%nnr )
         SELECT CASE ( type )
         CASE ( 1 ) ! Gaussian
            CALL generate_gradgaussian(nnr, dim, axis, charge, spread, pos, gradient%of_r)
         CASE ( 2 ) ! Erfc
            CALL generate_graderfc(nnr, dim, axis, charge, width, spread, pos, gradient%of_r)
         CASE ( 3 ) ! Exponential
            CALL generate_gradexponential(nnr, spread, pos, gradient%of_r)
         END SELECT
       END ASSOCIATE
    END DO

    RETURN

  END SUBROUTINE gradient_of_functions

  SUBROUTINE laplacian_of_functions( nfunctions, functions, laplacian )

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nfunctions
    TYPE( environ_functions ), DIMENSION(nfunctions), INTENT(IN) :: functions
    TYPE( environ_density ), INTENT(INOUT) :: laplacian

    INTEGER :: i

    laplacian%of_r = 0.D0

    DO i = 1, nfunctions
       ASSOCIATE( type   => functions(i)%type, &
                  pos    => functions(i)%pos, &
                  spread => functions(i)%spread, &
                  charge => functions(i)%volume, &
                  width  => functions(i)%width, &
                  dim    => functions(i)%dim, &
                  axis   => functions(i)%axis, &
                  nnr    => laplacian%cell%nnr )
         SELECT CASE ( type )
         CASE ( 1 ) ! Gaussian
            CALL errore(sub_name,'Options not yet implemented',1)
         CASE ( 2 ) ! Erfc
            CALL errore(sub_name,'Options not yet implemented',1)
         CASE ( 3 ) ! Exponential
            CALL errore(sub_name,'Options not yet implemented',1)
         END SELECT
       END ASSOCIATE
    END DO

    RETURN

  END SUBROUTINE laplacian_of_functions

END MODULE functions
