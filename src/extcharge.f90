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
MODULE externals_utils
!--------------------------------------------------------------------

  USE environ_types
  USE environ_output
  USE functions
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC :: create_environ_externals, init_environ_externals_first, &
       & init_environ_externals_second, update_environ_externals, destroy_environ_externals
  !
  CONTAINS
  SUBROUTINE create_environ_externals(externals)

    IMPLICIT NONE

    TYPE( environ_externals ), INTENT(INOUT) :: externals
    CHARACTER ( LEN=80 ) :: sub_name = 'create_environ_externals'
    CHARACTER ( LEN=80 ) :: label = 'externals'

    externals%update = .FALSE.
    externals%number = 0
    IF ( ALLOCATED( externals%functions ) ) CALL errore(sub_name,'Trying to create an already allocated object',1)
    CALL create_environ_density( externals%density, label )
    externals%charge = 0.D0

    RETURN

  END SUBROUTINE create_environ_externals

  SUBROUTINE init_environ_externals_first( nexternals, dims, axis, pos, &
       & spreads, charge, externals )

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nexternals
    INTEGER, DIMENSION(nexternals), INTENT(IN) :: dims, axis
    REAL( DP ), DIMENSION(3,nexternals), INTENT(IN) :: pos
    REAL( DP ), DIMENSION(nexternals), INTENT(IN) :: spreads, charge
    TYPE( environ_externals ), INTENT(INOUT) :: externals

    INTEGER :: i

    externals%number = nexternals
    ALLOCATE(externals%functions(externals%number))
    DO i = 1, externals%number
       ALLOCATE(externals%functions(i)%pos(3))
       externals%functions(i)%type   = 1
       externals%functions(i)%dim    = dims(i)
       externals%functions(i)%axis   = axis(i)
       externals%functions(i)%pos(:) = pos(:,i)
       externals%functions(i)%spread = spreads(i)
       externals%functions(i)%width  = spreads(i)
       externals%functions(i)%volume = charge(i)
    ENDDO

    RETURN

  END SUBROUTINE init_environ_externals_first

  SUBROUTINE init_environ_externals_second( cell, externals )

    IMPLICIT NONE

    TYPE( environ_cell ), INTENT(IN) :: cell
    TYPE( environ_externals ), INTENT(INOUT) :: externals

    INTEGER :: i

    IF ( externals % number .GT. 0 ) THEN
       DO i = 1, externals % number
          externals % functions(i) % pos = externals % functions(i) % pos / cell % alat
       END DO
    END IF

    CALL init_environ_density( cell, externals%density )

    RETURN

  END SUBROUTINE init_environ_externals_second

  SUBROUTINE update_environ_externals( externals )

    IMPLICIT NONE

    TYPE( environ_externals ), INTENT(INOUT) :: externals

    CALL density_of_functions( externals%number, externals%functions, externals%density, .TRUE. )

    externals % charge = integrate_environ_density( externals % density )

    RETURN

  END SUBROUTINE update_environ_externals

  SUBROUTINE destroy_environ_externals( lflag, externals )

    IMPLICIT NONE

    LOGICAL, INTENT(IN) :: lflag
    TYPE( environ_externals ), INTENT(INOUT) :: externals

    IF ( lflag ) CALL destroy_environ_functions( externals%number, externals%functions )

    CALL destroy_environ_density( externals%density )

    RETURN

  END SUBROUTINE destroy_environ_externals
!--------------------------------------------------------------------
END MODULE externals_utils
!--------------------------------------------------------------------
