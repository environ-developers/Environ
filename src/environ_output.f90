!
! Copyright (C) 2007-2008 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! This module contains all the subroutines related to the printout
! of Environ informations and results.
!
! Original version by Oliviero Andreussi and Nicola Marzari
!
!--------------------------------------------------------------------
MODULE environ_output
!--------------------------------------------------------------------

  USE environ_types

  SAVE

  LOGICAL :: ionode = .TRUE.
  INTEGER :: ionode_id

  INTEGER :: comm

  INTEGER :: program_unit
  INTEGER :: environ_unit
  INTEGER :: verbose
  INTEGER :: depth = 1

  CHARACTER( LEN = 2 ) :: prog

  INTEGER, PARAMETER :: nbibliography = 1
  CHARACTER( LEN = 80 ), DIMENSION( nbibliography ) :: bibliography

  PRIVATE

  PUBLIC :: ionode, ionode_id, comm, program_unit, environ_unit, &
       & verbose, prog, set_environ_output, environ_print_energies, &
       & environ_summary, environ_clock, write_cube, &
       & print_environ_density, print_environ_gradient, &
       & print_environ_functions, print_environ_iontype, &
       & print_environ_ions, print_environ_electrons, &
       & print_environ_externals, print_environ_system, &
       & print_environ_boundary, print_environ_dielectric

CONTAINS
!--------------------------------------------------------------------
  SUBROUTINE set_bibliography()
!--------------------------------------------------------------------
    IMPLICIT NONE

    bibliography(1) = "O. Andreussi, I. Dabo and N. Marzari, J.&
                      & Chem. Phys. 136, 064102 (2012)"

    RETURN
!--------------------------------------------------------------------
  END SUBROUTINE set_bibliography
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE set_environ_output( prog_, ionode_, ionode_id_, comm_, program_unit_ )
!--------------------------------------------------------------------

    IMPLICIT NONE

    CHARACTER( LEN=2 ), INTENT(IN) :: prog_
    LOGICAL, INTENT(IN) :: ionode_
    INTEGER, INTENT(IN) :: ionode_id_
    INTEGER, INTENT(IN) :: comm_
    INTEGER, INTENT(IN) :: program_unit_

    INTEGER, EXTERNAL :: find_free_unit

    CALL set_bibliography()

    ionode = ionode_
    ionode_id = ionode_id_
    comm = comm_

    program_unit = program_unit_
    environ_unit = find_free_unit()

    prog = prog_

    RETURN

!--------------------------------------------------------------------
  END SUBROUTINE set_environ_output
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE print_environ_cell( cell, local_verbose, local_depth )
!--------------------------------------------------------------------

    IMPLICIT NONE

    TYPE( environ_cell ), INTENT(IN) :: cell
    INTEGER, INTENT(IN), OPTIONAL :: local_verbose
    INTEGER, INTENT(IN), OPTIONAL :: local_depth

    INTEGER :: verbosity, passed_verbosity, passed_depth

    IF ( verbose .EQ. 0 ) RETURN ! environ output file has not been opened

    IF ( PRESENT(local_verbose) ) THEN
       verbosity = verbose + local_verbose
    ELSE
       verbosity = verbose
    END IF

    IF ( verbosity .EQ. 0 ) RETURN ! nothing to output

    IF ( PRESENT(local_depth) ) THEN
       passed_verbosity = verbosity - verbose - local_depth
       passed_depth = local_depth
    ELSE
       passed_verbosity = verbosity - verbose - depth
       passed_depth = depth
    END IF

    IF ( ionode .AND. verbosity .GE. 1 ) THEN
       IF ( verbosity .GE. verbose ) WRITE( UNIT = environ_unit, FMT = 1000 )
       WRITE( UNIT = environ_unit, FMT = 1001 )cell%ibrav,cell%alat,cell%omega
       IF ( verbosity .GE. 2 ) THEN
          WRITE( UNIT = environ_unit, FMT = 1002 )cell%at
          WRITE( UNIT = environ_unit, FMT = 1003 )cell%n1,cell%n2,cell%n3
          WRITE( UNIT = environ_unit, FMT = 1004 )cell%ntot,cell%nnr,cell%domega
          IF ( verbosity .GE. 3 ) THEN
             WRITE( UNIT = environ_unit, FMT = 1005 )cell%me,cell%root,cell%comm
          END IF
       END IF
    END IF

    RETURN

1000 FORMAT(/,4('%'),' CELL ',70('%'))
1001 FORMAT(1x,'bravais lattice index     = ',I3,' '&
          /,1x,'lattice spacing           = ',F12.6,' '&
          /,1x,'cell volume               = ',F12.6,' ')
1002 FORMAT(1x,'simulation cell axes      = ',3(F12.6),' '&
          /,1x,'                            ',3(F12.6),' '&
          /,1x,'                            ',3(F12.6),' ')
1003 FORMAT(1x,'real space grid dim.s     = ',3I4,' ')
1004 FORMAT(1x,'total size of grid        = ',I10,' '&
          /,1x,'size of r-space per proc. = ',I10,' '&
          /,1x,'finite element volume     = ',F12.6,' ')
1005 FORMAT(1x,'current processor index   = ',I10,' '&
          /,1x,'index of root processor   = ',I10,' '&
          /,1x,'communicator index        = ',I10,' ')
!--------------------------------------------------------------------
  END SUBROUTINE print_environ_cell
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE print_environ_density( density, local_verbose, local_depth, local_ions )
!--------------------------------------------------------------------

    USE environ_base, ONLY : ions

    IMPLICIT NONE

    TYPE( environ_density ), INTENT(IN) :: density
    INTEGER, INTENT(IN), OPTIONAL :: local_verbose
    INTEGER, INTENT(IN), OPTIONAL :: local_depth
    TYPE( environ_ions ), INTENT(IN), OPTIONAL :: local_ions

    INTEGER :: verbosity, passed_verbosity, passed_depth

    IF ( verbose .EQ. 0 ) RETURN ! environ output file has not been opened

    IF ( PRESENT(local_verbose) ) THEN
       verbosity = verbose + local_verbose
    ELSE
       verbosity = verbose
    END IF

    IF ( verbosity .EQ. 0 ) RETURN ! nothing to output

    IF ( PRESENT(local_depth) ) THEN
       passed_verbosity = verbosity - verbose - local_depth
       passed_depth = local_depth
    ELSE
       passed_verbosity = verbosity - verbose - depth
       passed_depth = depth
    END IF

    IF ( ionode .AND. verbosity .GE. 1 ) THEN
       IF ( verbosity .GE. verbose ) WRITE( UNIT = environ_unit, FMT = 1100 )
       WRITE( UNIT = environ_unit, FMT = 1101 )ADJUSTL(density%label)
       WRITE( UNIT = environ_unit, FMT = 1102 )integrate_environ_density(density)
       ! MAY ADD MAXVAL AND MINVAL
       IF ( verbosity .GE. 2 ) THEN
          CALL print_environ_cell( density%cell, passed_verbosity, passed_depth )
          IF ( PRESENT(local_ions) ) THEN
             WRITE(environ_unit,*)'present local ions'
             CALL write_cube( density, local_ions )
          ELSE IF ( ions%initialized ) THEN
             WRITE(environ_unit,*)'using stored ions'
             CALL write_cube( density, ions )
          ELSE
             WRITE(environ_unit,*)'no ions'
             CALL write_cube( density )
          END IF
       END IF
    END IF

    FLUSH( environ_unit )

    RETURN

1100 FORMAT(/,4('%'),' DENSITY ',67('%'))
1101 FORMAT(1x,'density label              = ',A80)
1102 FORMAT(1x,'integral of density        = ',G20.10)
!--------------------------------------------------------------------
  END SUBROUTINE print_environ_density
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE print_environ_gradient( gradient, local_verbose, local_depth, local_ions )
!--------------------------------------------------------------------

    USE environ_base, ONLY : ions

    IMPLICIT NONE

    TYPE( environ_gradient ), INTENT(IN) :: gradient
    INTEGER, INTENT(IN), OPTIONAL :: local_verbose
    INTEGER, INTENT(IN), OPTIONAL :: local_depth
    TYPE( environ_ions ), INTENT(IN), OPTIONAL :: local_ions

    INTEGER :: verbosity, passed_verbosity, passed_depth

    IF ( verbose .EQ. 0 ) RETURN ! environ output file has not been opened

    IF ( PRESENT(local_verbose) ) THEN
       verbosity = verbose + local_verbose
    ELSE
       verbosity = verbose
    END IF

    IF ( verbosity .EQ. 0 ) RETURN ! nothing to output

    IF ( PRESENT(local_depth) ) THEN
       passed_verbosity = verbosity - verbose - local_depth
       passed_depth = local_depth
    ELSE
       passed_verbosity = verbosity - verbose - depth
       passed_depth = depth
    END IF

    IF ( ionode .AND. verbosity .GE. 1 ) THEN
       IF ( verbosity .GE. verbose ) WRITE( UNIT = environ_unit, FMT = 1200 )
       WRITE( UNIT = environ_unit, FMT = 1201 )TRIM(ADJUSTL(gradient%label))
       WRITE( UNIT = environ_unit, FMT = 1201 )integrate_environ_density(gradient%modulus)
       ! MAY ADD MAXVAL AND MINVAL
       IF ( verbosity .GE. 2 ) THEN
          CALL print_environ_cell( gradient%cell, passed_verbosity, passed_depth )
          IF ( PRESENT(local_ions) ) THEN
             CALL write_cube( gradient%modulus, local_ions )
          ELSE IF ( ASSOCIATED(ions%tau) ) THEN
             CALL write_cube( gradient%modulus, ions )
          ELSE
             CALL write_cube( gradient%modulus )
          END IF
       END IF
    END IF

    FLUSH( environ_unit )

    RETURN

1200 FORMAT(/,4('%'),' GRADIENT ',66('%'))
1201 FORMAT(1x,'gradient label              = ',A80)
1202 FORMAT(1x,'integral of square modul.   = ',G20.10)
!--------------------------------------------------------------------
  END SUBROUTINE print_environ_gradient
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE print_environ_functions( nfunctions, functions, local_verbose, local_depth )
!--------------------------------------------------------------------

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nfunctions
    TYPE( environ_functions ), DIMENSION(nfunctions), INTENT(IN) :: functions
    INTEGER, INTENT(IN), OPTIONAL :: local_verbose
    INTEGER, INTENT(IN), OPTIONAL :: local_depth

    INTEGER :: verbosity, passed_verbosity, passed_depth
    INTEGER :: ifunctions

    CHARACTER( LEN=80 ) :: sub_name = 'print_environ_functions'

    IF ( verbose .EQ. 0 ) RETURN ! environ output file has not been opened

    IF ( PRESENT(local_verbose) ) THEN
       verbosity = verbose + local_verbose
    ELSE
       verbosity = verbose
    END IF

    IF ( verbosity .EQ. 0 ) RETURN ! nothing to output

    IF ( PRESENT(local_depth) ) THEN
       passed_verbosity = verbosity - verbose - local_depth
       passed_depth = local_depth
    ELSE
       passed_verbosity = verbosity - verbose - depth
       passed_depth = depth
    END IF

    IF ( ionode .AND. verbosity .GE. 1 ) THEN
       IF ( verbosity .GE. verbose ) WRITE( UNIT = environ_unit, FMT = 1300 )
       WRITE( UNIT = environ_unit, FMT = 1301 )nfunctions
       DO ifunctions = 1, nfunctions
          SELECT CASE (functions(ifunctions)%type )
          CASE ( 1 )
             WRITE( UNIT = environ_unit, FMT = 1302 )ifunctions,&
                  & functions(ifunctions)%dim,functions(ifunctions)%axis,&
                  & functions(ifunctions)%spread,functions(ifunctions)%volume
          CASE ( 2 )
             WRITE( UNIT = environ_unit, FMT = 1303 )ifunctions,&
                  & functions(ifunctions)%dim,functions(ifunctions)%axis,&
                  & functions(ifunctions)%width,functions(ifunctions)%spread,&
                  & functions(ifunctions)%volume
          CASE ( 3 )
             WRITE( UNIT = environ_unit, FMT = 1304 )ifunctions,&
                  & functions(ifunctions)%spread,functions(ifunctions)%volume
          CASE DEFAULT
             CALL errore(sub_name,'Unexpected function type',1)
          END SELECT
          IF ( verbosity .GE. 2 ) WRITE( UNIT = environ_unit, FMT = 1305 )&
               & functions(ifunctions)%pos
       END DO
    END IF

    FLUSH( environ_unit )

    RETURN

1300 FORMAT(/,4('%'),' FUNCTIONS ',65('%'))
1301 FORMAT(1x,'number of functions        = ',I10)
1302 FORMAT(1x,'Gaussian function, type    = 1 '&
          /,1x,'number                     = ',I3,' '&
          /,1x,'dimensionality             = ',I1,' '&
          /,1x,'axis                       = ',I1,' '&
          /,1x,'spread                     = ',F12.6,' '&
          /,1x,'integral                   = ',F12.6,' ')
1303 FORMAT(1x,'ERFC function, type        = 2 '&
          /,1x,'number                     = ',I3,' '&
          /,1x,'dimensionality             = ',I1,' '&
          /,1x,'axis                       = ',I1,' '&
          /,1x,'width                      = ',F12.6,' '&
          /,1x,'spread                     = ',F12.6,' '&
          /,1x,'integral                   = ',F12.6,' ')
1304 FORMAT(1x,'exponential function, type = 3 '&
          /,1x,'number                     = ',I3,' '&
          /,1x,'spread                     = ',F12.6,' '&
          /,1x,'integral                   = ',F12.6,' ')
1305 FORMAT(1x,'position                   = ',3(F14.7),' ')
!--------------------------------------------------------------------
  END SUBROUTINE print_environ_functions
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE print_environ_iontype( ntyp, iontype, local_verbose, local_depth )
!--------------------------------------------------------------------

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: ntyp
    TYPE( environ_iontype ), DIMENSION(ntyp), INTENT(IN) :: iontype
    INTEGER, INTENT(IN), OPTIONAL :: local_verbose
    INTEGER, INTENT(IN), OPTIONAL :: local_depth

    INTEGER :: verbosity, passed_verbosity, passed_depth
    INTEGER :: ityp

    CHARACTER( LEN=80 ) :: sub_name = 'print_environ_iontype'

    IF ( verbose .EQ. 0 ) RETURN ! environ output file has not been opened

    IF ( PRESENT(local_verbose) ) THEN
       verbosity = verbose + local_verbose
    ELSE
       verbosity = verbose
    END IF

    IF ( verbosity .EQ. 0 ) RETURN ! nothing to output

    IF ( PRESENT(local_depth) ) THEN
       passed_verbosity = verbosity - verbose - local_depth
       passed_depth = local_depth
    ELSE
       passed_verbosity = verbosity - verbose - depth
       passed_depth = depth
    END IF

    IF ( ionode .AND. verbosity .GE. 1 ) THEN
       IF ( verbosity .GE. verbose ) WRITE( UNIT = environ_unit, FMT = 1400 )
       WRITE( UNIT = environ_unit, FMT = 1401 )ntyp
       DO ityp = 1, ntyp
          WRITE( UNIT = environ_unit, FMT = 1402 )ityp,iontype(ityp)%index,&
               & iontype(ityp)%label,iontype(ityp)%atmnum,iontype(ityp)%zv
          IF ( verbosity .GE. 2 ) WRITE( UNIT = environ_unit, FMT = 1403 )&
               & iontype(ityp)%atomicspread,iontype(ityp)%corespread,&
               & iontype(ityp)%solvationrad
       END DO
    END IF

    FLUSH( environ_unit )

    RETURN

1400 FORMAT(/,4('%'),' IONTYPES ',66('%'))
1401 FORMAT(1x,'number of ionic types      = ',I10)
1402 FORMAT(1x,'iontype = ',I2,' index = ',I2,' label = ',A3,' atomic num = ',I2,' charge = ',F7.2)
1403 FORMAT(1x,'atomic spread = ',F7.2,' core spread = ',F7.2,' solvation radius = ',F7.2)
!--------------------------------------------------------------------
  END SUBROUTINE print_environ_iontype
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE print_environ_ions( ions, local_verbose, local_depth )
!--------------------------------------------------------------------

    IMPLICIT NONE

    TYPE( environ_ions ), INTENT(IN) :: ions
    INTEGER, INTENT(IN), OPTIONAL :: local_verbose
    INTEGER, INTENT(IN), OPTIONAL :: local_depth

    INTEGER :: verbosity, passed_verbosity, passed_depth

    INTEGER :: i
    CHARACTER( LEN=80 ) :: sub_name = 'print_environ_ions'

    IF ( verbose .EQ. 0 ) RETURN ! environ output file has not been opened

    IF ( PRESENT(local_verbose) ) THEN
       verbosity = verbose + local_verbose
    ELSE
       verbosity = verbose
    END IF

    IF ( verbosity .EQ. 0 ) RETURN ! nothing to output

    IF ( PRESENT(local_depth) ) THEN
       passed_verbosity = verbosity - verbose - local_depth
       passed_depth = local_depth
    ELSE
       passed_verbosity = verbosity - verbose - depth
       passed_depth = depth
    END IF

    IF ( ionode .AND. verbosity .GE. 1 ) THEN
       IF ( verbosity .GE. verbose ) WRITE( UNIT = environ_unit, FMT = 1500 )
       WRITE( UNIT = environ_unit, FMT = 1501 )ions%number
       WRITE( UNIT = environ_unit, FMT = 1502 )ions%charge
       WRITE( UNIT = environ_unit, FMT = 1503 )ions%center
       DO i = 1, ions%number
          WRITE( UNIT = environ_unit, FMT = 1504 )i,ions%ityp(i),ions%tau(:,i)
       END DO
       IF ( verbosity .GE. 2 ) THEN
          CALL print_environ_iontype(ions%ntyp,ions%iontype,passed_verbosity,passed_depth)
          IF ( ions%use_smeared_ions ) THEN
             CALL print_environ_density(ions%density,passed_verbosity,passed_depth)
             IF ( verbosity .GE. 3 ) &
                  & CALL print_environ_functions(ions%number,&
                  & ions%smeared_ions,passed_verbosity,passed_depth)
          END IF
          IF ( ions%use_core_electrons ) THEN
             IF ( verbosity .GE. 3 ) CALL print_environ_density(ions%core,passed_verbosity,passed_depth)
             IF ( verbosity .GE. 4 ) CALL print_environ_functions(ions%number,&
                  & ions%core_electrons,passed_verbosity,passed_depth)
          END IF
       END IF
    END IF

    FLUSH( environ_unit )

    RETURN

1500 FORMAT(/,4('%'),' IONS ',70('%'))
1501 FORMAT(1x,'number of ions             = ',I10)
1502 FORMAT(1x,'total ionic charge         = ',F14.7)
1503 FORMAT(1x,'ionic center of charge     = ',3(F14.7))
1504 FORMAT(1x,'ion ',I3,' type = ',I2,' coordinates = ',3(F14.7))
!--------------------------------------------------------------------
  END SUBROUTINE print_environ_ions
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE print_environ_electrons( electrons, local_verbose, local_depth )
!--------------------------------------------------------------------

    IMPLICIT NONE

    TYPE( environ_electrons ), INTENT(IN) :: electrons
    INTEGER, INTENT(IN), OPTIONAL :: local_verbose
    INTEGER, INTENT(IN), OPTIONAL :: local_depth

    INTEGER :: verbosity, passed_verbosity, passed_depth

    CHARACTER( LEN=80 ) :: sub_name = 'print_environ_electrons'

    IF ( verbose .EQ. 0 ) RETURN ! environ output file has not been opened

    IF ( PRESENT(local_verbose) ) THEN
       verbosity = verbose + local_verbose
    ELSE
       verbosity = verbose
    END IF

    IF ( verbosity .EQ. 0 ) RETURN ! nothing to output

    IF ( PRESENT(local_depth) ) THEN
       passed_verbosity = verbosity - verbose - local_depth
       passed_depth = local_depth
    ELSE
       passed_verbosity = verbosity - verbose - depth
       passed_depth = depth
    END IF

    IF ( ionode .AND. verbosity .GE. 1 ) THEN
       IF ( verbosity .GE. verbose ) WRITE( UNIT = environ_unit, FMT = 1600 )
       WRITE( UNIT = environ_unit, FMT = 1601 )electrons%number
       WRITE( UNIT = environ_unit, FMT = 1602 )electrons%charge
       WRITE( UNIT = environ_unit, FMT = 1603 )electrons%nspin
       IF ( verbosity .GE. 2 ) &
            & CALL print_environ_density(electrons%density,passed_verbosity,passed_depth)
    END IF

    FLUSH( environ_unit )

    RETURN

1600 FORMAT(/,4('%'),' ELECTRONS ',65('%'))
1601 FORMAT(1x,'number of electrons        = ',I10)
1602 FORMAT(1x,'total electronic charge    = ',F14.7)
1603 FORMAT(1x,'number of spin components  = ',I2)
!--------------------------------------------------------------------
  END SUBROUTINE print_environ_electrons
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE print_environ_externals( externals, local_verbose, local_depth )
!--------------------------------------------------------------------

    IMPLICIT NONE

    TYPE( environ_externals ), INTENT(IN) :: externals
    INTEGER, INTENT(IN), OPTIONAL :: local_verbose
    INTEGER, INTENT(IN), OPTIONAL :: local_depth

    INTEGER :: verbosity, passed_verbosity, passed_depth

    CHARACTER( LEN=80 ) :: sub_name = 'print_environ_externals'

    IF ( verbose .EQ. 0 ) RETURN ! environ output file has not been opened

    IF ( PRESENT(local_verbose) ) THEN
       verbosity = verbose + local_verbose
    ELSE
       verbosity = verbose
    END IF

    IF ( verbosity .EQ. 0 ) RETURN ! nothing to output

    IF ( PRESENT(local_depth) ) THEN
       passed_verbosity = verbosity - verbose - local_depth
       passed_depth = local_depth
    ELSE
       passed_verbosity = verbosity - verbose - depth
       passed_depth = depth
    END IF

    IF ( ionode .AND. verbosity .GE. 1 ) THEN
       IF ( verbosity .GE. verbose ) WRITE( UNIT = environ_unit, FMT = 1700 )
       WRITE( UNIT = environ_unit, FMT = 1701 )externals%number
       WRITE( UNIT = environ_unit, FMT = 1702 )externals%charge
       IF ( verbosity .GE. 2 ) THEN
          CALL print_environ_functions(externals%number,externals%functions,&
               & passed_verbosity,passed_depth)
          CALL print_environ_density(externals%density,passed_verbosity,passed_depth)
       ENDIF
    END IF

    FLUSH( environ_unit )

    RETURN

1700 FORMAT(/,4('%'),' EXTERNALS ',65('%'))
1701 FORMAT(1x,'number of external charges = ',I10)
1702 FORMAT(1x,'total external charge      = ',F14.7)
!--------------------------------------------------------------------
  END SUBROUTINE print_environ_externals
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE print_environ_charges( charges, local_verbose, local_depth )
!--------------------------------------------------------------------

    IMPLICIT NONE

    TYPE( environ_charges ), INTENT(IN) :: charges
    INTEGER, INTENT(IN), OPTIONAL :: local_verbose
    INTEGER, INTENT(IN), OPTIONAL :: local_depth

    INTEGER :: verbosity, passed_verbosity, passed_depth

    CHARACTER( LEN=80 ) :: sub_name = 'print_environ_charges'

    IF ( verbose .EQ. 0 ) RETURN ! environ output file has not been opened

    IF ( PRESENT(local_verbose) ) THEN
       verbosity = verbose + local_verbose
    ELSE
       verbosity = verbose
    END IF

    IF ( verbosity .EQ. 0 ) RETURN ! nothing to output

    IF ( PRESENT(local_depth) ) THEN
       passed_verbosity = verbosity - verbose - local_depth
       passed_depth = local_depth
    ELSE
       passed_verbosity = verbosity - verbose - depth
       passed_depth = depth
    END IF

    IF ( ionode .AND. verbosity .GE. 1 ) THEN
       IF ( verbosity .GE. verbose ) WRITE( UNIT = environ_unit, FMT = 1800 )
       WRITE( UNIT = environ_unit, FMT = 1801 )charges%number
       WRITE( UNIT = environ_unit, FMT = 1802 )charges%charge
       IF ( verbosity .GE. 2 ) THEN
          CALL print_environ_density(charges%density,passed_verbosity,passed_depth)
          IF ( charges % include_ions ) &
               & CALL print_environ_ions(charges%ions,passed_verbosity,passed_depth)
          IF ( charges % include_electrons ) &
               & CALL print_environ_electrons(charges%electrons,passed_verbosity,passed_depth)
          IF ( charges % include_externals ) &
               & CALL print_environ_externals(charges%externals,passed_verbosity,passed_depth)
          IF ( charges % include_auxiliary ) &
               & CALL print_environ_density(charges%auxiliary,passed_verbosity,passed_depth)
       ENDIF
    END IF

    FLUSH( environ_unit )

    RETURN

1800 FORMAT(/,4('%'),' CHARGES ',67('%'))
1801 FORMAT(1x,'total number of charges    = ',I10)
1802 FORMAT(1x,'total charge               = ',F14.7)
!--------------------------------------------------------------------
  END SUBROUTINE print_environ_charges
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE print_environ_system( system, local_verbose, local_depth )
!--------------------------------------------------------------------

    IMPLICIT NONE

    TYPE( environ_system ), INTENT(IN) :: system
    INTEGER, INTENT(IN), OPTIONAL :: local_verbose
    INTEGER, INTENT(IN), OPTIONAL :: local_depth

    INTEGER :: verbosity, passed_verbosity, passed_depth

    CHARACTER( LEN=80 ) :: sub_name = 'print_environ_system'

    IF ( verbose .EQ. 0 ) RETURN ! environ output file has not been opened

    IF ( PRESENT(local_verbose) ) THEN
       verbosity = verbose + local_verbose
    ELSE
       verbosity = verbose
    END IF

    IF ( verbosity .EQ. 0 ) RETURN ! nothing to output

    IF ( PRESENT(local_depth) ) THEN
       passed_verbosity = verbosity - verbose - local_depth
       passed_depth = local_depth
    ELSE
       passed_verbosity = verbosity - verbose - depth
       passed_depth = depth
    END IF

    IF ( ionode .AND. verbosity .GE. 1 ) THEN
       IF ( verbosity .GE. verbose ) WRITE( UNIT = environ_unit, FMT = 1900 )
       IF ( system%ntyp .EQ. 0 ) THEN
          WRITE( UNIT = environ_unit, FMT = 1901 )
       ELSE
          WRITE( UNIT = environ_unit, FMT = 1902 )system%ntyp
       END IF
       WRITE( UNIT = environ_unit, FMT = 1903 )system%dim,system%axis
       WRITE( UNIT = environ_unit, FMT = 1904 )system%pos,system%width
       IF ( verbosity .GE. 2 ) &
          & CALL print_environ_ions(system % ions, passed_verbosity, passed_depth )
    END IF

    FLUSH( environ_unit )

    RETURN

1900 FORMAT(/,4('%'),' SYSTEM ',68('%'))
1901 FORMAT(1x,'system is built from all present ionic types')
1902 FORMAT(1x,'system is built from the first ',I3,' ionic types')
1903 FORMAT(1x,'system defined dimension   = ',I2,' '&
          /,1x,'system defined axis        = ',I2,' ')
1904 FORMAT(1x,'system center              = ',3F14.7,' '&
          /,1x,'system width               = ',F14.7,' ')
!--------------------------------------------------------------------
  END SUBROUTINE print_environ_system
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE print_environ_boundary( boundary, local_verbose, local_depth )
!--------------------------------------------------------------------

    IMPLICIT NONE

    TYPE( environ_boundary ), INTENT(IN) :: boundary
    INTEGER, INTENT(IN), OPTIONAL :: local_verbose
    INTEGER, INTENT(IN), OPTIONAL :: local_depth

    INTEGER :: verbosity, passed_verbosity, passed_depth

    CHARACTER( LEN=80 ) :: sub_name = 'print_environ_boundary'

    IF ( verbose .EQ. 0 ) RETURN ! environ output file has not been opened

    IF ( PRESENT(local_verbose) ) THEN
       verbosity = verbose + local_verbose
    ELSE
       verbosity = verbose
    END IF

    IF ( verbosity .EQ. 0 ) RETURN ! nothing to output

    IF ( PRESENT(local_depth) ) THEN
       passed_verbosity = verbosity - verbose - local_depth
       passed_depth = local_depth
    ELSE
       passed_verbosity = verbosity - verbose - depth
       passed_depth = depth
    END IF

    IF ( ionode .AND. verbosity .GE. 1 ) THEN
       IF ( verbosity .GE. verbose ) WRITE( UNIT = environ_unit, FMT = 2000 )
       WRITE( UNIT = environ_unit, FMT = 2001 )boundary%mode
       IF ( boundary % need_electrons ) THEN
          WRITE( UNIT = environ_unit, FMT = 2002 ) boundary % type
          SELECT CASE ( boundary % type )
          CASE ( 1 )
             WRITE( UNIT = environ_unit, FMT = 2003 ) boundary % rhomax, boundary % rhomin
             IF ( verbosity .GE. 2 ) WRITE( UNIT = environ_unit, FMT = 2004 ) boundary % fact
          CASE ( 2 )
             WRITE( UNIT = environ_unit, FMT = 2005 ) boundary % rhozero, boundary % tbeta
          END SELECT
          IF ( verbosity .GE. 2 ) THEN
             CALL print_environ_density(boundary%density,passed_verbosity,passed_depth)
             IF ( boundary % need_ions ) THEN
                WRITE( UNIT = environ_unit, FMT = 2006 )
                IF ( verbosity .GE. 3 ) CALL print_environ_density(boundary%ions%core)
             END IF
             IF ( verbosity .GE. 3 ) &
                  & CALL print_environ_electrons(boundary%electrons,passed_verbosity,passed_depth)
          END IF
       ELSE
          WRITE( UNIT = environ_unit, FMT = 2007 ) boundary%alpha, boundary%softness
          IF ( verbosity .GE. 2 ) &
               & CALL print_environ_functions(boundary%ions%number,boundary%soft_spheres,&
               & passed_verbosity,passed_depth)
       END IF
       IF ( verbosity .GE. 2 ) CALL print_environ_density(boundary%scaled,passed_verbosity,passed_depth)
       IF ( verbosity .GE. 3 ) CALL print_environ_density(boundary%dscaled,passed_verbosity,passed_depth)
       IF ( verbosity .GE. 4 ) CALL print_environ_density(boundary%d2scaled,passed_verbosity,passed_depth)
       IF ( boundary%need_theta ) THEN
          WRITE( UNIT = environ_unit, FMT = 2008 ) boundary % deltatheta
          IF ( verbosity .GE. 2 ) CALL print_environ_density(boundary%theta,passed_verbosity,passed_depth)
       ENDIF
    END IF

    FLUSH( environ_unit )

    RETURN

2000 FORMAT(/,4('%'),' BOUNDARY ',66('%'))
2001 FORMAT(1x,'boundary mode              = ',A20,' ')
2002 FORMAT(1x,'boundary is built as a function of a smooth density'&
          /,1x,'function type              = ',I2,' ')
2003 FORMAT(1x,'using the optimal SCCS function with parameters '&
          /,1x,'rhomax                     = ',F14.7,' '&
          /,1x,'rhomin                     = ',F14.7,' ')
2004 FORMAT(1x,'log(rhomax/rhomin)         = ',F14.7,' ')
2005 FORMAT(1x,'using the Fattebert-Gygi function with parameters '&
          /,1x,'rhozero                    = ',F14.7,' '&
          /,1x,'2*beta                     = ',F14.7,' ')
2006 FORMAT(1x,'adding fictitious core-electrons')
2007 FORMAT(1x,'boundary is built from soft-spheres centered on ionic positions'&
          /,1x,'solvent-dependent scaling  = ',F14.7,' '&
          /,1x,'softness parameter         = ',F14.7,' ')
2008 FORMAT(1x,'also need the surface term of this boundary'&
          /,1x,'fd parameter for surface   = ',F14.7,' ')
!--------------------------------------------------------------------
  END SUBROUTINE print_environ_boundary
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE print_environ_dielectric( dielectric, local_verbose, local_depth )
!--------------------------------------------------------------------

    IMPLICIT NONE

    TYPE( environ_dielectric ), INTENT(IN) :: dielectric
    INTEGER, INTENT(IN), OPTIONAL :: local_verbose
    INTEGER, INTENT(IN), OPTIONAL :: local_depth

    INTEGER :: verbosity, passed_verbosity, passed_depth

    CHARACTER( LEN=80 ) :: sub_name = 'print_environ_dielectric'

    IF ( verbose .EQ. 0 ) RETURN ! environ output file has not been opened

    IF ( PRESENT(local_verbose) ) THEN
       verbosity = verbose + local_verbose
    ELSE
       verbosity = verbose
    END IF

    IF ( verbosity .EQ. 0 ) RETURN ! nothing to output

    IF ( PRESENT(local_depth) ) THEN
       passed_verbosity = verbosity - verbose - local_depth
       passed_depth = local_depth
    ELSE
       passed_verbosity = verbosity - verbose - depth
       passed_depth = depth
    END IF

    IF ( ionode .AND. verbosity .GE. 1 ) THEN
       IF ( verbosity .GE. verbose ) WRITE( UNIT = environ_unit, FMT = 2100 )
       IF ( dielectric % nregions .EQ. 0 ) THEN
          WRITE( UNIT = environ_unit, FMT = 2101 )dielectric%constant
       ELSE
          WRITE( UNIT = environ_unit, FMT = 2102 )dielectric%constant,dielectric%nregions
          CALL print_environ_functions(dielectric%nregions,dielectric%regions,passed_verbosity,passed_depth)
          IF ( verbosity .GE. 2 ) CALL print_environ_density(dielectric%background,passed_verbosity,passed_depth)
       END IF
       CALL print_environ_boundary(dielectric%boundary,passed_verbosity,passed_depth)
       IF ( verbosity .GE. 2 ) CALL print_environ_density(dielectric%epsilon,passed_verbosity,passed_depth)
       WRITE( UNIT = environ_unit, FMT = 2103 )dielectric%need_gradient,&
            & dielectric%need_factsqrt,dielectric%need_gradlog
       IF ( verbosity .GE. 3 ) THEN
          IF ( dielectric%need_gradient ) CALL print_environ_gradient(dielectric%gradient,passed_verbosity,passed_depth)
          IF ( dielectric%need_factsqrt ) CALL print_environ_density(dielectric%factsqrt,passed_verbosity,passed_depth)
          IF ( dielectric%need_gradlog ) CALL print_environ_gradient(dielectric%gradlog,passed_verbosity,passed_depth)
       END IF

    END IF

    FLUSH( environ_unit )

    RETURN

2100 FORMAT(/,4('%'),' DIELECTRIC ',65('%'))
2101 FORMAT(1x,'dielectric build on homogeneous background'&
          /,1x,'environment bulk permitt.  = ',F14.7,' ')
2102 FORMAT(1x,'dielectric build in the presence of dielectric regions'&
          /,1x,'environment bulk permitt.  = ',F14.7,' '&
          /,1x,'number of dielec. regions  = ',I4,' ')
2103 FORMAT(1x,'dielectric flags'&
          /,1x,'need gradient              = ',L2,' '&
          /,1x,'need factor depend. sqrt   = ',L2,' '&
          /,1x,'need gradient of logarithm = ',L2,' ')
!--------------------------------------------------------------------
  END SUBROUTINE print_environ_dielectric
!--------------------------------------------------------------------
!--------------------------------------------------------------------
      SUBROUTINE environ_print_energies( )
!--------------------------------------------------------------------
      !
      ! Write out the different Environ contributions to the energy.
      ! Called by electrons.f90
      !
      USE environ_base, ONLY : lelectrostatic, eelectrostatic, &
                               lsurface, ecavity, &
                               lvolume, epressure
      !
      IF ( ionode ) THEN
         IF ( prog .EQ. 'PW' ) THEN
            IF ( lelectrostatic ) WRITE( program_unit, 9201 ) eelectrostatic
            IF ( lsurface ) WRITE( program_unit, 9202 ) ecavity
            IF ( lvolume ) WRITE( program_unit, 9203 ) epressure
         ELSE IF ( prog .EQ. 'CP' ) THEN
            IF ( lelectrostatic ) WRITE( program_unit, 9301 ) eelectrostatic
            IF ( lsurface ) WRITE( program_unit, 9302 ) ecavity
            IF ( lvolume ) WRITE( program_unit, 9303 ) epressure
         ELSE
            WRITE(program_unit,*)'ERROR: wrong program in Environ'
            STOP
         END IF
      END IF
      !
      RETURN
      !
9201 FORMAT( '     electrostatic embedding   =',F17.8,' Ry')
9202 FORMAT( '     cavitation energy         =',F17.8,' Ry')
9203 FORMAT( '     PV energy                 =',F17.8,' Ry')
9301 FORMAT( '     electrostatic embedding = ',F14.5,' Hartree a.u.')
9302 FORMAT( '           cavitation energy = ',F14.5,' Hartree a.u.')
9303 FORMAT( '                   PV energy = ',F14.5,' Hartree a.u.')
      !
! ---------------------------------------------------------------------
      END SUBROUTINE environ_print_energies
! ---------------------------------------------------------------------
! ---------------------------------------------------------------------
      SUBROUTINE environ_summary( )
! ---------------------------------------------------------------------
      !
      ! Write out the main parameters of Environ calculations,
      ! summarizing the input keywords (some info also on internal
      ! vs input units). Called by summary.f90
      !
      USE environ_base,       ONLY : environ_thr, lsolvent, solvent,   &
                                     env_static_permittivity,          &
                                     env_optical_permittivity,         &
                                     env_surface_tension,              &
                                     env_pressure, lelectrostatic
      USE electrostatic_base, ONLY : problem, solver, auxiliary,       &
                                     tolvelect, tolrhoaux,             &
                                     core, dielectric_core,            &
                                     ifdtype, nfdpoint
      !
      IMPLICIT NONE
      !
      !
      IF( ionode ) THEN
         !
         IF ( prog .EQ. 'PW' ) THEN
            !
            WRITE( program_unit, * )
            !
            WRITE( UNIT = program_unit, FMT = 9000 )
            WRITE( UNIT = program_unit, FMT = 8000 ) bibliography(1)
            !
            ! ... Environ Summary
            !
            WRITE( UNIT = program_unit, FMT = 9001 ) environ_thr
            !
            IF ( lsolvent ) THEN
               IF ( solvent%type .EQ. 0 ) THEN
                  WRITE( UNIT = program_unit, FMT = 9002 ) 'Fatteber-Gygi'
                  WRITE( UNIT = program_unit, FMT = 9003 ) solvent%rhozero, solvent%tbeta
               ELSE IF ( solvent%type .EQ. 1 ) THEN
                  WRITE( UNIT = program_unit, FMT = 9002 ) 'SCCS'
                  WRITE( UNIT = program_unit, FMT = 9004 ) solvent%rhomax, solvent%rhomin
               ENDIF
            ENDIF
            !
            IF ( env_static_permittivity .GT. 1.D0 ) THEN
               WRITE( UNIT = program_unit, FMT = 9005 ) env_static_permittivity
               IF (tddfpt) &
                    & WRITE( UNIT = program_unit, FMT = 9006 ) env_optical_permittivity
               WRITE( UNIT = program_unit, FMT = 9007 ) TRIM( solvent%mode )
            END IF
            !
            IF ( env_surface_tension .GT. 0.D0 ) WRITE( UNIT = program_unit, FMT = 9010 )      &
                 env_surface_tension/1.D-3/bohr_radius_si**2*rydberg_si, env_surface_tension, solvent%deltatheta
            !
            IF ( env_pressure .NE. 0.D0 ) WRITE( UNIT = program_unit, FMT = 9011 )&
                 env_pressure*rydberg_si/bohr_radius_si**3*1.D-9, env_pressure
            !
            ! ... Electrostatic Summary
            !
            IF ( lelectrostatic ) THEN
               !
               WRITE( UNIT = program_unit, FMT = 9100 )
               !
               WRITE( UNIT = program_unit, FMT = 9101 )problem, solver, auxiliary
               !
               WRITE( UNIT = program_unit, FMT = 9102 )tolvelect, tolrhoaux
               !
               WRITE( UNIT = program_unit, FMT = 9103 )core, dielectric_core
               !
               IF ( core .EQ. 'fd' .OR. dielectric_core .EQ. 'fd' ) THEN
                  IF ( ifdtype .EQ. 1 ) THEN
                     WRITE( UNIT = program_unit, FMT = 9104 ) 'central diff.',nfdpoint
                  ELSE IF (ifdtype .EQ. 2 .OR. ifdtype .EQ. 3 ) THEN
                     WRITE( UNIT = program_unit, FMT = 9104 ) 'lanczos diff.',nfdpoint
                  ELSE IF (ifdtype .EQ.4 .OR. ifdtype .EQ. 5 ) THEN
                     WRITE( UNIT = program_unit, FMT = 9104 ) 'noise-robust diff.',nfdpoint
                  END IF
               END IF
               !
            END IF
            !
            WRITE( UNIT = program_unit, FMT = 8001 )
            !
         END IF
         !
      END IF
      !
8000  FORMAT(/,5x,'Plese cite',/,9x,A80,&
             /,5x,'in publications or presentations arising from this work.',/)
8001  FORMAT(/)
9000  FORMAT(/,5x,'Environ Module',/,5x,'==============')
9001  FORMAT( '     compensation onset threshold      = ',  E24.4,' ' )
9002  FORMAT( '     switching function adopted        = ',  A24,' ' )
9003  FORMAT( '     solvation density threshold       = ',  E24.4,' ' &
             /'     smoothness exponent (2 x beta)    = ',  F24.2,' ' )
9004  FORMAT( '     density limit for vacuum region   = ',  E24.4,' ' &
             /'     density limit for bulk solvent    = ',  E24.4,' ' )
9005  FORMAT( '     static permittivity               = ',  F24.2,' ' )
9006  FORMAT( '     optical permittivity              = ',  F24.4,' ' )
9007  FORMAT( '     epsilon calculation mode          = ',  A24,' ' )
9010  FORMAT( '     surface tension in input (dyn/cm) = ',  F24.2,' ' &
             /'     surface tension in internal units = ',  E24.4,' ' &
             /'     delta parameter for surface depth = ',  E24.4,' ' )
9011  FORMAT( '     external pressure in input (GPa)  = ',  F24.2,' ' &
             /'     external pressure in inter. units = ',  E24.4,' ' )
9012  FORMAT( '     correction slab geom. along axis  = ',  I24,' ' )
9100  FORMAT(/,5x,'Electrostatic Setup',/,5x,'-------------------')
9101  FORMAT( '     electrostatic problem to solve    = ',  A24,' ' &
             /'     numerical solver adopted          = ',  A24,' ' &
             /'     type of auxiliary density adopted = ',  A24,' ' )
9102  FORMAT( '     required accuracy on potential    = ',  F24.4,' ' &
             /'     required accuracy on aux. charge  = ',  F24.4,' ' )
9103  FORMAT( '     type of core tool for poisson     = ',  A24,' ' &
             /'     type of core tool for eps deriv   = ',  A24,' ' )
9104  FORMAT( '     type of numerical differentiator  = ',  A24,' ' &
             /'     number of points in num. diff.    = ',  I24,' ' )

!--------------------------------------------------------------------
      END SUBROUTINE environ_summary
!--------------------------------------------------------------------
!--------------------------------------------------------------------
      SUBROUTINE environ_clock( )
!--------------------------------------------------------------------
      !
      ! Write out the time informations of the Environ dependent
      ! calculations. Called by print_clock_pw.f90
      !
      USE environ_base,   ONLY : lelectrostatic, lsurface, lvolume
      !
      IMPLICIT NONE
      !
      WRITE( program_unit, * )
      WRITE( program_unit, '(5X,"Environ routines")' )
      ! dielectric subroutines
      IF ( lelectrostatic ) THEN
         CALL print_clock ('calc_eelect')
         CALL print_clock ('calc_velect')
         CALL print_clock ('dielectric')
         CALL print_clock ('calc_felect')
      END IF
      ! cavitation subroutines
      IF ( lsurface ) THEN
         CALL print_clock ('calc_ecav')
         CALL print_clock ('calc_vcav')
      END IF
      ! pressure subroutines
      IF ( lvolume ) THEN
         CALL print_clock ('calc_epre')
         CALL print_clock ('calc_vpre')
      END IF
      ! TDDFT
      IF (tddfpt) CALL print_clock ('calc_vsolvent_tddfpt')
      !
      RETURN
      !
!--------------------------------------------------------------------
      END SUBROUTINE environ_clock
!--------------------------------------------------------------------
!--------------------------------------------------------------------
      SUBROUTINE write_cube( f, ions )
!--------------------------------------------------------------------
      !
      USE fft_base,       ONLY : dfftp
! BACKWARD COMPATIBILITY
! Compatible with QE-5.1 QE-5.1.1 QE-5.1.2
!      USE fft_base,       ONLY : grid_gather
! Compatible with QE-5.2 QE-5.2.1
!      USE fft_base,       ONLY : gather_grid
! Compatible with QE-5.3 svn
      USE scatter_mod,    ONLY : gather_grid
! END BACKWARD COMPATIBILITY
      !
      IMPLICIT NONE
      !
      TYPE( environ_density ), TARGET, INTENT(IN) :: f
      TYPE( environ_ions ), TARGET, OPTIONAL, INTENT(IN) :: ions
      !
      INTEGER                  :: ir, ir1, ir2, ir3, num
      INTEGER                  :: ipol, iat, typ, count
      INTEGER                  :: nr1x, nr2x, nr3x
      INTEGER                  :: nr1, nr2, nr3
      !
      REAL( DP )               :: tmp, scale
      REAL( DP ), ALLOCATABLE  :: flocal( : )
      REAL( DP ), DIMENSION(3) :: origin
      !
      CHARACTER( LEN=100 ) :: filename
      REAL( DP ), POINTER :: alat
      REAL( DP ), DIMENSION(:,:), POINTER :: at
      !
      INTEGER :: nat
      INTEGER, DIMENSION(:), POINTER :: ityp
      REAL( DP ), DIMENSION(:,:), POINTER :: tau
      !
      nr1x = dfftp%nr1x
      nr2x = dfftp%nr2x
      nr3x = dfftp%nr3x
      !
      nr1 = dfftp%nr1
      nr2 = dfftp%nr2
      nr3 = dfftp%nr3
      !
      filename = TRIM(ADJUSTL(f%label))//".cube"
      !
      alat => f%cell%alat
      at => f%cell%at
      !
      IF ( PRESENT( ions ) ) THEN
         nat = ions%number
         ityp => ions%ityp
         tau => ions%tau
      ELSE
         nat = 1
      ENDIF
      !
      ALLOCATE( flocal( nr1x*nr2x*nr3x ) )
#ifdef __MPI
      flocal = 0.D0
!Compatible with QE-5.1 QE-5.1.1 QE-5.1.2
!      CALL grid_gather( f, flocal )
!Compatible with QE-svn
      CALL gather_grid( dfftp, f%of_r, flocal )
      CALL mp_sum( flocal, comm )
#else
      flocal = f%of_r
#endif
      !
      IF( ionode ) THEN
        !
        OPEN( 300, file = TRIM( filename ), status = 'unknown' )
        !
        origin=0.d0
        scale=alat!*0.52917720859d0
        WRITE(300,*)'CUBE FILE GENERATED BY PW.X'
        WRITE(300,*)'OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z'
        WRITE(300,'(i5,3f12.6)')nat,origin(1),origin(2),origin(3)
        WRITE(300,'(i5,3f12.6)')nr1,(at(ipol,1)/DBLE(nr1)*scale,ipol=1,3)
        WRITE(300,'(i5,3f12.6)')nr2,(at(ipol,2)/DBLE(nr2)*scale,ipol=1,3)
        WRITE(300,'(i5,3f12.6)')nr3,(at(ipol,3)/DBLE(nr3)*scale,ipol=1,3)
        IF ( PRESENT( ions ) ) THEN
           DO iat=1,nat
              typ=ityp(iat)
              num=ions%iontype(typ)%atmnum
              WRITE(300,'(i5,4f12.6)')num,0.d0,tau(1,iat)*scale,&
                   tau(2,iat)*scale,tau(3,iat)*scale
           ENDDO
        ELSE
           WRITE(300,'(i5,4f12.6)')1,0.d0,0.d0,0.d0,0.d0
        ENDIF
        count=0
        DO ir1=1,nr1
          DO ir2=1,nr2
            DO ir3=1,nr3
              count=count+1
              ir = ir1 + ( ir2 -1 ) * nr1 + ( ir3 - 1 ) * nr1 * nr2
              tmp = DBLE( flocal( ir ) )
              IF (ABS(tmp).LT.1.D-99) tmp = 0.D0
              IF (MOD(count,6).EQ.0) THEN
                WRITE(300,'(e12.6,1x)')tmp
              ELSE
                WRITE(300,'(e12.6,1x)',advance='no')tmp
              ENDIF
            ENDDO
          ENDDO
        ENDDO
        !
        CLOSE( 300 )
        !
      END IF
      !
      DEALLOCATE( flocal )
      !
      RETURN
      !
!--------------------------------------------------------------------
      END SUBROUTINE write_cube
!--------------------------------------------------------------------
!--------------------------------------------------------------------
END MODULE environ_output
!--------------------------------------------------------------------
