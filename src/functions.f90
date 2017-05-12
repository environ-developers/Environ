MODULE functions

  USE environ_types
  USE generate_function

  PRIVATE

  PUBLIC :: create_environ_functions, destroy_environ_functions, &
       & density_of_functions, gradient_of_functions, laplacian_of_functions

CONTAINS

  SUBROUTINE create_environ_functions( n, type, dimm, axis, pos, width, spreadd, volume, f )

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: n
    INTEGER, INTENT(IN) :: type
    INTEGER, DIMENSION(:), INTENT(IN) :: dimm, axis
    REAL( DP ), DIMENSION(n), INTENT(IN) :: width, spreadd, volume
    REAL( DP ), DIMENSION(3,n), TARGET, INTENT(IN) :: pos
    TYPE( environ_functions ), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: f

    INTEGER :: i

    ALLOCATE( f(n) )
    DO i = 1, n
       f(i)%type = type
       f(i)%dim  = dimm(i)
       f(i)%axis = axis(i)
       f(i)%spread = spreadd(i)
       f(i)%width = width(i)
       f(i)%volume = volume(i)
       f(i)%pos => pos(:,i)
    ENDDO

    RETURN

  END SUBROUTINE create_environ_functions

  SUBROUTINE destroy_environ_functions( n, f )

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: n
    TYPE( environ_functions ), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: f

    CHARACTER( LEN=80 ) :: sub_name = 'destroy_environ_functions'

    INTEGER :: i

    IF ( .NOT. ALLOCATED( f ) ) &
         & CALL errore(sub_name,'Trying to destroy a non allocated object',1)
    IF ( SIZE(f) .NE. n ) &
         & CALL errore(sub_name,'Inconsistent size of allocated object',1)

    DO i = 1, n
       IF ( .NOT. ASSOCIATED( f(i)%pos ) ) &
            & CALL errore(sub_name,'Trying to destroy a non allocated object',1)
       NULLIFY( f(i)%pos )
    ENDDO

    DEALLOCATE( f )

    RETURN

  END SUBROUTINE destroy_environ_functions

  SUBROUTINE density_of_functions( nfunctions, functions, density )

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nfunctions
    TYPE( environ_functions ), DIMENSION(nfunctions), TARGET, INTENT(IN) :: functions
    TYPE( environ_density ), TARGET, INTENT(INOUT) :: density

    INTEGER :: i
    CHARACTER( LEN=80 ) :: sub_name = 'density_of_functions'

    INTEGER, POINTER :: type, dim, axis, nnr
    REAL( DP ), POINTER :: charge, spread, width
    REAL( DP ), DIMENSION(:), POINTER :: pos

    density%of_r = 0.D0

    DO i = 1, nfunctions
       type   => functions(i)%type
       pos    => functions(i)%pos
       spread => functions(i)%spread
       charge => functions(i)%volume
       width  => functions(i)%width
       dim    => functions(i)%dim
       axis   => functions(i)%axis
       nnr    => density%cell%nnr
       SELECT CASE ( type )
       CASE ( 1 ) ! Gaussian
          CALL generate_gaussian(nnr, dim, axis, charge, spread, pos, density%of_r)
       CASE ( 2 ) ! Erfc
          CALL generate_erfc(nnr, dim, axis, charge, width, spread, pos, density%of_r)
       CASE ( 3 ) ! Exponential
          CALL generate_exponential(nnr, spread, pos, density%of_r)
          !
          ! TO IMPLEMENT THE DISTINCTION BETWEEN NORMALIZED ERFC AND MAXED ERFC
          !
       END SELECT
    END DO

    RETURN

  END SUBROUTINE density_of_functions

  SUBROUTINE gradient_of_functions( nfunctions, functions, gradient )

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nfunctions
    TYPE( environ_functions ), DIMENSION(nfunctions), TARGET, INTENT(IN) :: functions
    TYPE( environ_gradient ), TARGET, INTENT(INOUT) :: gradient

    INTEGER :: i
    CHARACTER( LEN=80 ) :: sub_name = 'gradient_of_functions'

    INTEGER, POINTER :: type, dim, axis, nnr
    REAL( DP ), POINTER :: charge, spread, width
    REAL( DP ), DIMENSION(:), POINTER :: pos

    gradient%of_r = 0.D0

    DO i = 1, nfunctions
       type   => functions(i)%type
       pos    => functions(i)%pos
       spread => functions(i)%spread
       charge => functions(i)%volume
       width  => functions(i)%width
       dim    => functions(i)%dim
       axis   => functions(i)%axis
       nnr    => gradient%cell%nnr
       SELECT CASE ( type )
       CASE ( 1 ) ! Gaussian
          CALL generate_gradgaussian(nnr, dim, axis, charge, spread, pos, gradient%of_r)
       CASE ( 2 ) ! Erfc
          CALL generate_graderfc(nnr, dim, axis, charge, width, spread, pos, gradient%of_r)
       CASE ( 3 ) ! Exponential
          CALL generate_gradexponential(nnr, spread, pos, gradient%of_r)
       END SELECT
    END DO

    RETURN

  END SUBROUTINE gradient_of_functions

  SUBROUTINE laplacian_of_functions( nfunctions, functions, laplacian )

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nfunctions
    TYPE( environ_functions ), DIMENSION(nfunctions), TARGET, INTENT(IN) :: functions
    TYPE( environ_density ), TARGET, INTENT(INOUT) :: laplacian

    INTEGER :: i
    CHARACTER( LEN=80 ) :: sub_name = 'laplacian_of_functions'

    INTEGER, POINTER :: type, dim, axis, nnr
    REAL( DP ), POINTER :: charge, spread, width
    REAL( DP ), DIMENSION(:), POINTER :: pos

    laplacian%of_r = 0.D0

    DO i = 1, nfunctions
       type   => functions(i)%type
       pos    => functions(i)%pos
       spread => functions(i)%spread
       charge => functions(i)%volume
       width  => functions(i)%width
       dim    => functions(i)%dim
       axis   => functions(i)%axis
       nnr    => laplacian%cell%nnr
       SELECT CASE ( type )
       CASE ( 1 ) ! Gaussian
          CALL errore(sub_name,'Options not yet implemented',1)
       CASE ( 2 ) ! Erfc
          CALL errore(sub_name,'Options not yet implemented',1)
       CASE ( 3 ) ! Exponential
          CALL errore(sub_name,'Options not yet implemented',1)
       END SELECT
    END DO

    RETURN

  END SUBROUTINE laplacian_of_functions

END MODULE functions
