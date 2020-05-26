! Copyright (C) 2018 ENVIRON (www.quantum-environment.org)
!
!    This file is part of Environ version 1.1
!
!    Environ 1.1 is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 2 of the License, or
!    (at your option) any later version.
!
!    Environ 1.1 is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more detail, either the file
!    `License' in the root directory of the present distribution, or
!    online at <http://www.gnu.org/licenses/>.
!
!> This module provides output subroutines for Environ, including summary
!! of input parameters, timings, print of energy contributions, and print
!! subroutines for Environ derived data types
!
! Authors: Oliviero Andreussi (Department of Physics, UNT)
!          Francesco Nattino  (THEOS and NCCR-MARVEL, EPFL)
!          Nicola Marzari     (THEOS and NCCR-MARVEL, EPFL)
!
!----------------------------------------------------------------------------
MODULE environ_output
!----------------------------------------------------------------------------
  !
  USE environ_types
  !
  SAVE
  !
  LOGICAL :: ionode = .TRUE.
  INTEGER :: ionode_id
  !
  INTEGER :: comm
  !
  INTEGER :: program_unit
  INTEGER :: environ_unit
  INTEGER :: verbose
  INTEGER :: depth = 1
  !
  CHARACTER( LEN = 2 ) :: prog
  !
  INTEGER, PARAMETER :: nbibliography = 6
  CHARACTER( LEN = 80 ), DIMENSION( nbibliography ) :: bibliography
  !
  PRIVATE
  !
  PUBLIC :: ionode, ionode_id, comm, program_unit, environ_unit, &
       & verbose, prog, set_environ_output, environ_print_energies, &
       & environ_print_potential_shift, environ_print_potential_warning, &
       & environ_summary, environ_clock, write_cube, &
       & print_environ_density, print_environ_gradient, &
       & print_environ_hessian, &
       & print_environ_functions, print_environ_iontype, &
       & print_environ_ions, print_environ_electrons, &
       & print_environ_externals, print_environ_system, &
       & print_environ_boundary, print_environ_dielectric, &
       & print_environ_electrolyte, print_environ_charges, &
       & update_output_program_unit
  !
CONTAINS
!--------------------------------------------------------------------
  SUBROUTINE set_bibliography()
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    bibliography(1) = "O. Andreussi, I. Dabo and N. Marzari, J.&
                      & Chem. Phys. 136, 064102 (2012)"
    bibliography(2) = "I. Timrov, O. Andreussi, A. Biancardi, N.&
                      & Marzari, and S. Baroni, J. Chem. Phys.&
                      & 142, 034111 (2015)"
    bibliography(3) = "O. Andreussi, I. Dabo and N. Marzari, J.&
                      & Chem. Phys. 136, 064102 (2012)"
    bibliography(4) = "O. Andreussi, I. Dabo and N. Marzari, J.&
                      & Chem. Phys. 136, 064102 (2012)"
    bibliography(5) = "O. Andreussi, N.G. Hoermann, F. Nattino,&
                      & G. Fisicaro, S. Goedecker, and N. Marzari,&
                      & J. Chem. Theory Comput. 15, 1996 (2019)"
    bibliography(6) = "F. Nattino, M. Truscott, N. Marzari, and&
                      & O. Andreussi, J. Chem. Phys. 150, 041722&
                      & (2019)"
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE set_bibliography
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE set_environ_output( prog_, ionode_, ionode_id_, comm_, program_unit_ )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    CHARACTER( LEN=* ), INTENT(IN) :: prog_
    LOGICAL, INTENT(IN) :: ionode_
    INTEGER, INTENT(IN) :: ionode_id_
    INTEGER, INTENT(IN) :: comm_
    INTEGER, INTENT(IN) :: program_unit_
    !
    INTEGER, EXTERNAL :: find_free_unit
    !
    CALL set_bibliography()
    !
    ionode = ionode_
    ionode_id = ionode_id_
    comm = comm_
    !
    program_unit = program_unit_
    environ_unit = find_free_unit()
    !
    prog = prog_(1:2)
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE set_environ_output
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE update_output_program_unit( program_unit_ )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: program_unit_
    !
    program_unit = program_unit_
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE update_output_program_unit
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE print_environ_cell( cell, local_verbose, local_depth )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    TYPE( environ_cell ), INTENT(IN) :: cell
    INTEGER, INTENT(IN), OPTIONAL :: local_verbose
    INTEGER, INTENT(IN), OPTIONAL :: local_depth
    !
    INTEGER :: verbosity, passed_verbosity, passed_depth
    !
    IF ( verbose .EQ. 0 ) RETURN ! environ output file has not been opened
    !
    IF ( PRESENT(local_verbose) ) THEN
       verbosity = verbose + local_verbose
    ELSE
       verbosity = verbose
    END IF
    !
    IF ( verbosity .EQ. 0 ) RETURN ! nothing to output
    !
    IF ( PRESENT(local_depth) ) THEN
       passed_verbosity = verbosity - verbose - local_depth
       passed_depth = local_depth
    ELSE
       passed_verbosity = verbosity - verbose
       passed_depth = depth
    END IF
    !
    IF ( verbosity .GE. 1 ) THEN
       IF ( verbosity .GE. verbose .AND. ionode ) WRITE( UNIT = environ_unit, FMT = 1000 )
       IF ( ionode ) WRITE( UNIT = environ_unit, FMT = 1001 )cell%ibrav,cell%alat,cell%omega
       IF ( verbosity .GE. 3 ) THEN
          IF ( ionode ) WRITE( UNIT = environ_unit, FMT = 1002 )cell%at
          IF ( ionode ) WRITE( UNIT = environ_unit, FMT = 1003 )cell%n1,cell%n2,cell%n3
          IF ( ionode ) WRITE( UNIT = environ_unit, FMT = 1004 )cell%ntot,cell%nnr,cell%domega
          IF ( verbosity .GE. 4 ) THEN
             IF ( ionode ) WRITE( UNIT = environ_unit, FMT = 1005 )cell%me,cell%root,cell%comm
          END IF
       END IF
    END IF
    !
    RETURN
    !
1000 FORMAT(/,4('%'),' CELL ',70('%'))
1001 FORMAT(1x,'bravais lattice index      = ',I3,' '&
          /,1x,'lattice spacing            = ',F12.6,' '&
          /,1x,'cell volume                = ',F12.6,' ')
1002 FORMAT(1x,'simulation cell axes       = ',3(F12.6),' '&
          /,1x,'                             ',3(F12.6),' '&
          /,1x,'                             ',3(F12.6),' ')
1003 FORMAT(1x,'real space grid dim.s      = ',3I4,' ')
1004 FORMAT(1x,'total size of grid         = ',I10,' '&
          /,1x,'size of r-space per proc.  = ',I10,' '&
          /,1x,'finite element volume      = ',F12.6,' ')
1005 FORMAT(1x,'current processor index    = ',I10,' '&
          /,1x,'index of root processor    = ',I10,' '&
          /,1x,'communicator index         = ',I10,' ')
!--------------------------------------------------------------------
  END SUBROUTINE print_environ_cell
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE print_environ_density( density, local_verbose, local_depth, local_ions )
!--------------------------------------------------------------------
    !
    USE environ_base, ONLY : ions
    !
    IMPLICIT NONE
    !
    TYPE( environ_density ), INTENT(IN) :: density
    INTEGER, INTENT(IN), OPTIONAL :: local_verbose
    INTEGER, INTENT(IN), OPTIONAL :: local_depth
    TYPE( environ_ions ), INTENT(IN), OPTIONAL :: local_ions
    !
    INTEGER :: verbosity, passed_verbosity, passed_depth
    REAL( DP ) :: integral
    !
    IF ( verbose .EQ. 0 ) RETURN ! environ output file has not been opened
    !
    IF ( PRESENT(local_verbose) ) THEN
       verbosity = verbose + local_verbose
    ELSE
       verbosity = verbose
    END IF
    !
    IF ( verbosity .EQ. 0 ) RETURN ! nothing to output
    !
    IF ( PRESENT(local_depth) ) THEN
       passed_verbosity = verbosity - verbose - local_depth
       passed_depth = local_depth
    ELSE
       passed_verbosity = verbosity - verbose 
       passed_depth = depth
    END IF
    !
    IF ( verbosity .GE. 1 ) THEN
       IF ( verbosity .GE. verbose .AND. ionode ) WRITE( UNIT = environ_unit, FMT = 1100 )
       IF ( ionode ) WRITE( UNIT = environ_unit, FMT = 1101 )ADJUSTL(density%label)
       integral = integrate_environ_density(density)
       IF ( ionode ) WRITE( UNIT = environ_unit, FMT = 1102 )integral
       ! MAY ADD MAXVAL AND MINVAL
       IF ( verbosity .GE. 3 ) THEN
          CALL print_environ_cell( density%cell, passed_verbosity, passed_depth )
          IF ( PRESENT(local_ions) ) THEN
             IF ( ionode ) WRITE(environ_unit,*)'present local ions'
             CALL write_cube( density, local_ions )
          ELSE IF ( ions%initialized ) THEN
             IF ( ionode ) WRITE(environ_unit,*)'using stored ions'
             CALL write_cube( density, ions )
          ELSE
             IF ( ionode ) WRITE(environ_unit,*)'no ions'
             CALL write_cube( density )
          END IF
       END IF
    END IF
    !
    FLUSH( environ_unit )
    !
    RETURN
    !
1100 FORMAT(/,4('%'),' DENSITY ',67('%'))
1101 FORMAT(1x,'density label              = ',A80)
1102 FORMAT(1x,'integral of density        = ',G20.10)
!--------------------------------------------------------------------
  END SUBROUTINE print_environ_density
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE print_environ_gradient( gradient, local_verbose, local_depth, local_ions )
!--------------------------------------------------------------------
    !
    USE environ_base, ONLY : ions
    !
    IMPLICIT NONE
    !
    TYPE( environ_gradient ), INTENT(IN) :: gradient
    INTEGER, INTENT(IN), OPTIONAL :: local_verbose
    INTEGER, INTENT(IN), OPTIONAL :: local_depth
    TYPE( environ_ions ), INTENT(IN), OPTIONAL :: local_ions
    !
    TYPE( environ_cell ), POINTER :: cell
    TYPE( environ_density ) :: dens
    !
    INTEGER :: verbosity, passed_verbosity, passed_depth
    REAL( DP ) :: integral
    !
    IF ( verbose .EQ. 0 ) RETURN ! environ output file has not been opened
    !
    IF ( PRESENT(local_verbose) ) THEN
       verbosity = verbose + local_verbose
    ELSE
       verbosity = verbose
    END IF
    !
    IF ( verbosity .EQ. 0 ) RETURN ! nothing to output
    !
    IF ( PRESENT(local_depth) ) THEN
       passed_verbosity = verbosity - verbose - local_depth
       passed_depth = local_depth
    ELSE
       passed_verbosity = verbosity - verbose 
       passed_depth = depth
    END IF
    !
    IF ( verbosity .GE. 1 ) THEN
       IF ( verbosity .GE. verbose .AND. ionode ) WRITE( UNIT = environ_unit, FMT = 1200 )
       IF ( ionode ) WRITE( UNIT = environ_unit, FMT = 1201 )ADJUSTL(gradient%label)
       integral = integrate_environ_density(gradient%modulus)
       IF ( ionode ) WRITE( UNIT = environ_unit, FMT = 1202 )integral
       ! MAY ADD MAXVAL AND MINVAL
       IF ( verbosity .GE. 3 ) THEN
          CALL print_environ_cell( gradient%cell, passed_verbosity, passed_depth )
          IF ( PRESENT(local_ions) ) THEN
             CALL write_cube( gradient%modulus, local_ions )
          ELSE IF ( ASSOCIATED(ions%tau) ) THEN
             CALL write_cube( gradient%modulus, ions )
          ELSE
             CALL write_cube( gradient%modulus )
          END IF
       END IF
       IF ( verbosity .GE. 4 ) THEN
          cell => gradient % cell
          CALL init_environ_density( cell, dens )
          dens % label = TRIM(ADJUSTL(gradient % label))//'_x'
          dens % of_r( : ) = gradient % of_r( 1, : )
          CALL print_environ_density( dens, passed_verbosity, passed_depth )
          dens % label = TRIM(ADJUSTL(gradient % label))//'_y'
          dens % of_r( : ) = gradient % of_r( 2, : )
          CALL print_environ_density( dens, passed_verbosity, passed_depth )
          dens % label = TRIM(ADJUSTL(gradient % label))//'_z'
          dens % of_r( : ) = gradient % of_r( 3, : )
          CALL print_environ_density( dens, passed_verbosity, passed_depth )
          CALL destroy_environ_density( dens )
       ENDIF
    END IF
    !
    FLUSH( environ_unit )
    !
    RETURN
    !
1200 FORMAT(/,4('%'),' GRADIENT ',66('%'))
1201 FORMAT(1x,'gradient label             = ',A80)
1202 FORMAT(1x,'integral of square modul.  = ',G20.10)
!--------------------------------------------------------------------
  END SUBROUTINE print_environ_gradient
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE print_environ_hessian( hessian, local_verbose, local_depth, local_ions )
!--------------------------------------------------------------------
    !
    USE environ_base, ONLY : ions
    !
    IMPLICIT NONE
    !
    TYPE( environ_hessian ), INTENT(IN) :: hessian
    INTEGER, INTENT(IN), OPTIONAL :: local_verbose
    INTEGER, INTENT(IN), OPTIONAL :: local_depth
    TYPE( environ_ions ), INTENT(IN), OPTIONAL :: local_ions
    !
    TYPE( environ_cell ), POINTER :: cell
    TYPE( environ_density ) :: dens
    !
    INTEGER :: verbosity, passed_verbosity, passed_depth
    REAL( DP ) :: integral
    !
    IF ( verbose .EQ. 0 ) RETURN ! environ output file has not been opened
    !
    IF ( PRESENT(local_verbose) ) THEN
       verbosity = verbose + local_verbose
    ELSE
       verbosity = verbose
    END IF
    !
    IF ( verbosity .EQ. 0 ) RETURN ! nothing to output
    !
    IF ( PRESENT(local_depth) ) THEN
       passed_verbosity = verbosity - verbose - local_depth
       passed_depth = local_depth
    ELSE
       passed_verbosity = verbosity - verbose
       passed_depth = depth
    END IF
    !
    IF ( verbosity .GE. 1 ) THEN
       IF ( verbosity .GE. verbose .AND. ionode ) WRITE( UNIT = environ_unit, FMT = 1250 )
       IF ( ionode ) WRITE( UNIT = environ_unit, FMT = 1251 )ADJUSTL(hessian%label)
       integral = integrate_environ_density(hessian%laplacian)
       IF ( ionode ) WRITE( UNIT = environ_unit, FMT = 1252 )integral
       ! MAY ADD MAXVAL AND MINVAL
       IF ( verbosity .GE. 3 ) THEN
          CALL print_environ_cell( hessian%cell, passed_verbosity, passed_depth )
          IF ( PRESENT(local_ions) ) THEN
             CALL write_cube( hessian%laplacian, local_ions )
          ELSE IF ( ASSOCIATED(ions%tau) ) THEN
             CALL write_cube( hessian%laplacian, ions )
          ELSE
             CALL write_cube( hessian%laplacian )
          END IF
       END IF
       IF ( verbosity .GE. 4 ) THEN
          cell => hessian % cell
          CALL init_environ_density( cell, dens )
          dens % label = TRIM(ADJUSTL( hessian % label ))//'_xx'
          dens % of_r( : ) = hessian % of_r( 1, 1, : )
          CALL print_environ_density( dens, passed_verbosity, passed_depth )
          dens % label = TRIM(ADJUSTL( hessian % label ))//'_xy'
          dens % of_r( : ) = hessian % of_r( 1, 2, : )
          CALL print_environ_density( dens, passed_verbosity, passed_depth )
          dens % label = TRIM(ADJUSTL( hessian % label ))//'_xz'
          dens % of_r( : ) = hessian % of_r( 1, 3, : )
          CALL print_environ_density( dens, passed_verbosity, passed_depth )
          dens % label = TRIM(ADJUSTL( hessian % label ))//'_yy'
          dens % of_r( : ) = hessian % of_r( 2, 2, : )
          CALL print_environ_density( dens, passed_verbosity, passed_depth )
          dens % label = TRIM(ADJUSTL( hessian % label ))//'_yz'
          dens % of_r( : ) = hessian % of_r( 2, 3, : )
          CALL print_environ_density( dens, passed_verbosity, passed_depth )
          dens % label = TRIM(ADJUSTL( hessian % label ))//'_zz'
          dens % of_r( : ) = hessian % of_r( 3, 3, : )
          CALL print_environ_density( dens, passed_verbosity, passed_depth )
          CALL destroy_environ_density( dens )
       ENDIF
    END IF
    !
    FLUSH( environ_unit )
    !
    RETURN
    !
1250 FORMAT(/,4('%'),' HESSIAN ',67('%'))
1251 FORMAT(1x,'hessian label              = ',A80)
1252 FORMAT(1x,'integral of laplacian      = ',G20.10)
!--------------------------------------------------------------------
  END SUBROUTINE print_environ_hessian
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE print_environ_functions( nfunctions, functions, local_verbose, local_depth )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: nfunctions
    TYPE( environ_functions ), DIMENSION(nfunctions), INTENT(IN) :: functions
    INTEGER, INTENT(IN), OPTIONAL :: local_verbose
    INTEGER, INTENT(IN), OPTIONAL :: local_depth
    !
    INTEGER :: verbosity, passed_verbosity, passed_depth
    INTEGER :: ifunctions
    !
    CHARACTER( LEN=80 ) :: sub_name = 'print_environ_functions'
    !
    IF ( verbose .EQ. 0 ) RETURN ! environ output file has not been opened
    !
    IF ( PRESENT(local_verbose) ) THEN
       verbosity = verbose + local_verbose
    ELSE
       verbosity = verbose
    END IF
    !
    IF ( verbosity .EQ. 0 ) RETURN ! nothing to output
    !
    IF ( PRESENT(local_depth) ) THEN
       passed_verbosity = verbosity - verbose - local_depth
       passed_depth = local_depth
    ELSE
       passed_verbosity = verbosity - verbose
       passed_depth = depth
    END IF
    !
    IF ( verbosity .GE. 1 ) THEN
       IF ( verbosity .GE. verbose .AND. ionode ) WRITE( UNIT = environ_unit, FMT = 1300 )
       IF ( ionode ) WRITE( UNIT = environ_unit, FMT = 1301 )nfunctions
       DO ifunctions = 1, nfunctions
          SELECT CASE (functions(ifunctions)%type )
          CASE ( 1 )
             IF ( ionode ) WRITE( UNIT = environ_unit, FMT = 1302 )ifunctions,&
                  & functions(ifunctions)%dim,functions(ifunctions)%axis,&
                  & functions(ifunctions)%spread,functions(ifunctions)%volume
          CASE ( 2 )
             IF ( ionode ) WRITE( UNIT = environ_unit, FMT = 1303 )ifunctions,&
                  & functions(ifunctions)%dim,functions(ifunctions)%axis,&
                  & functions(ifunctions)%width,functions(ifunctions)%spread,&
                  & functions(ifunctions)%volume
          CASE ( 3 )
             IF ( ionode ) WRITE( UNIT = environ_unit, FMT = 1304 )ifunctions,&
                  & functions(ifunctions)%spread,functions(ifunctions)%volume
          CASE ( 4 )
             IF ( ionode ) WRITE( UNIT = environ_unit, FMT = 1305 )ifunctions,&
                  & functions(ifunctions)%dim,functions(ifunctions)%axis,&
                  & functions(ifunctions)%width,functions(ifunctions)%spread,&
                  & functions(ifunctions)%volume
          CASE ( 5 )
             IF ( ionode ) WRITE( UNIT = environ_unit, FMT = 1306 )ifunctions,&
                  & functions(ifunctions)%dim,functions(ifunctions)%axis,&
                  & functions(ifunctions)%width,functions(ifunctions)%spread,&
                  & functions(ifunctions)%volume
          CASE DEFAULT
             CALL errore(sub_name,'Unexpected function type',1)
          END SELECT
          IF ( verbosity .GE. 3 .AND. ionode ) WRITE( UNIT = environ_unit, FMT = 1307 )&
               & functions(ifunctions)%pos
       END DO
    END IF
    !
    FLUSH( environ_unit )
    !
    RETURN
    !
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
1305 FORMAT(1x,'Scaled ERFC function, type = 4 '&
          /,1x,'number                     = ',I3,' '&
          /,1x,'dimensionality             = ',I1,' '&
          /,1x,'axis                       = ',I1,' '&
          /,1x,'width                      = ',F12.6,' '&
          /,1x,'spread                     = ',F12.6,' '&
          /,1x,'max value                  = ',F12.6,' ')
1306 FORMAT(1x,'Scaled ERF function, type  = 5 '&
          /,1x,'number                     = ',I3,' '&
          /,1x,'dimensionality             = ',I1,' '&
          /,1x,'axis                       = ',I1,' '&
          /,1x,'width                      = ',F12.6,' '&
          /,1x,'spread                     = ',F12.6,' '&
          /,1x,'max value                  = ',F12.6,' ')
1307 FORMAT(1x,'position                   = ',3(F14.7),' ')
!--------------------------------------------------------------------
  END SUBROUTINE print_environ_functions
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE print_environ_iontype( ntyp, iontype, local_verbose, local_depth )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: ntyp
    TYPE( environ_iontype ), DIMENSION(ntyp), INTENT(IN) :: iontype
    INTEGER, INTENT(IN), OPTIONAL :: local_verbose
    INTEGER, INTENT(IN), OPTIONAL :: local_depth
    !
    INTEGER :: verbosity, passed_verbosity, passed_depth
    INTEGER :: ityp
    !
    CHARACTER( LEN=80 ) :: sub_name = 'print_environ_iontype'
    !
    IF ( verbose .EQ. 0 ) RETURN ! environ output file has not been opened
    !
    IF ( PRESENT(local_verbose) ) THEN
       verbosity = verbose + local_verbose
    ELSE
       verbosity = verbose
    END IF
    !
    IF ( verbosity .EQ. 0 ) RETURN ! nothing to output
    !
    IF ( PRESENT(local_depth) ) THEN
       passed_verbosity = verbosity - verbose - local_depth
       passed_depth = local_depth
    ELSE
       passed_verbosity = verbosity - verbose 
       passed_depth = depth
    END IF
    !
    IF ( verbosity .GE. 1 ) THEN
       IF ( verbosity .GE. verbose .AND. ionode ) WRITE( UNIT = environ_unit, FMT = 1400 )
       IF ( ionode ) WRITE( UNIT = environ_unit, FMT = 1401 )ntyp
       DO ityp = 1, ntyp
          IF ( ionode ) WRITE( UNIT = environ_unit, FMT = 1402 )ityp,iontype(ityp)%index,&
               & iontype(ityp)%label,iontype(ityp)%atmnum,iontype(ityp)%zv
          IF ( verbosity .GE. 3 .AND. ionode ) WRITE( UNIT = environ_unit, FMT = 1403 )&
               & iontype(ityp)%atomicspread,iontype(ityp)%corespread,&
               & iontype(ityp)%solvationrad
       END DO
    END IF
    !
    FLUSH( environ_unit )
    !
    RETURN
    !
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
    !
    IMPLICIT NONE
    !
    TYPE( environ_ions ), INTENT(IN) :: ions
    INTEGER, INTENT(IN), OPTIONAL :: local_verbose
    INTEGER, INTENT(IN), OPTIONAL :: local_depth
    !
    INTEGER :: verbosity, passed_verbosity, passed_depth
    !
    INTEGER :: i
    CHARACTER( LEN=80 ) :: sub_name = 'print_environ_ions'
    !
    IF ( verbose .EQ. 0 ) RETURN ! environ output file has not been opened
    !
    IF ( PRESENT(local_verbose) ) THEN
       verbosity = verbose + local_verbose
    ELSE
       verbosity = verbose
    END IF
    !
    IF ( verbosity .EQ. 0 ) RETURN ! nothing to output
    !
    IF ( PRESENT(local_depth) ) THEN
       passed_verbosity = verbosity - verbose - local_depth
       passed_depth = local_depth
    ELSE
       passed_verbosity = verbosity - verbose 
       passed_depth = depth
    END IF
    !
    IF ( verbosity .GE. 1 ) THEN
       IF ( verbosity .GE. verbose .AND. ionode ) WRITE( UNIT = environ_unit, FMT = 1500 )
       IF ( ionode ) WRITE( UNIT = environ_unit, FMT = 1501 )ions%number
       IF ( ionode ) WRITE( UNIT = environ_unit, FMT = 1502 )ions%center
       IF ( ionode ) WRITE( UNIT = environ_unit, FMT = 1503 )ions%charge
       IF ( ionode ) WRITE( UNIT = environ_unit, FMT = 1504 )ions%dipole
       IF ( ionode ) WRITE( UNIT = environ_unit, FMT = 1505 )ions%quadrupole_pc
       IF ( ions%use_smeared_ions .AND. ionode ) WRITE( UNIT = environ_unit, FMT = 1506 )ions%quadrupole_gauss
       DO i = 1, ions%number
          IF ( ionode ) WRITE( UNIT = environ_unit, FMT = 1507 )i,ions%ityp(i),ions%tau(:,i)*ions%alat
       END DO
       IF ( verbosity .GE. 3 ) THEN
          CALL print_environ_iontype(ions%ntyp,ions%iontype,passed_verbosity,passed_depth)
          IF ( ions%use_smeared_ions ) THEN
             CALL print_environ_density(ions%density,passed_verbosity,passed_depth)
             IF ( verbosity .GE. 4 ) &
                  & CALL print_environ_functions(ions%number,&
                  & ions%smeared_ions,passed_verbosity,passed_depth)
          END IF
          IF ( ions%use_core_electrons ) THEN
             IF ( verbosity .GE. 4 ) CALL print_environ_density(ions%core,passed_verbosity,passed_depth)
             IF ( verbosity .GE. 5 ) CALL print_environ_functions(ions%number,&
                  & ions%core_electrons,passed_verbosity,passed_depth)
          END IF
       END IF
    END IF
    !
    FLUSH( environ_unit )
    !
    RETURN
    !
1500 FORMAT(/,4('%'),' IONS ',70('%'))
1501 FORMAT(1x,'number of ions             = ',I10)
1502 FORMAT(1x,'ionic center of charge     = ',3(F14.7))
1503 FORMAT(1x,'total ionic charge         = ',F14.7)
1504 FORMAT(1x,'ionic dipole               = ',3(F14.7))
1505 FORMAT(1x,'ionic quadrupole (pc)      = ',3(F14.7))
1506 FORMAT(1x,'ionic quadrupole (gauss)   = ',3(F14.7))
1507 FORMAT(1x,'ion ',I3,' type = ',I2,' coordinates = ',3(F14.7))
!--------------------------------------------------------------------
  END SUBROUTINE print_environ_ions
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE print_environ_electrons( electrons, local_verbose, local_depth )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    TYPE( environ_electrons ), INTENT(IN) :: electrons
    INTEGER, INTENT(IN), OPTIONAL :: local_verbose
    INTEGER, INTENT(IN), OPTIONAL :: local_depth
    !
    INTEGER :: verbosity, passed_verbosity, passed_depth
    !
    CHARACTER( LEN=80 ) :: sub_name = 'print_environ_electrons'
    !
    IF ( verbose .EQ. 0 ) RETURN ! environ output file has not been opened
    !
    IF ( PRESENT(local_verbose) ) THEN
       verbosity = verbose + local_verbose
    ELSE
       verbosity = verbose
    END IF
    !
    IF ( verbosity .EQ. 0 ) RETURN ! nothing to output
    !
    IF ( PRESENT(local_depth) ) THEN
       passed_verbosity = verbosity - verbose - local_depth
       passed_depth = local_depth
    ELSE
       passed_verbosity = verbosity - verbose
       passed_depth = depth
    END IF
    !
    IF ( verbosity .GE. 1 ) THEN
       IF ( verbosity .GE. verbose .AND. ionode ) WRITE( UNIT = environ_unit, FMT = 1600 )
       IF ( ionode ) WRITE( UNIT = environ_unit, FMT = 1601 )electrons%number
       IF ( ionode ) WRITE( UNIT = environ_unit, FMT = 1602 )electrons%charge
       IF ( verbosity .GE. 3 ) &
            & CALL print_environ_density(electrons%density,passed_verbosity,passed_depth)
    END IF
    !
    FLUSH( environ_unit )
    !
    RETURN
    !
1600 FORMAT(/,4('%'),' ELECTRONS ',65('%'))
1601 FORMAT(1x,'number of electrons        = ',I10)
1602 FORMAT(1x,'total electronic charge    = ',F14.7)
!--------------------------------------------------------------------
  END SUBROUTINE print_environ_electrons
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE print_environ_externals( externals, local_verbose, local_depth )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    TYPE( environ_externals ), INTENT(IN) :: externals
    INTEGER, INTENT(IN), OPTIONAL :: local_verbose
    INTEGER, INTENT(IN), OPTIONAL :: local_depth
    !
    INTEGER :: verbosity, passed_verbosity, passed_depth
    !
    CHARACTER( LEN=80 ) :: sub_name = 'print_environ_externals'
    !
    IF ( verbose .EQ. 0 ) RETURN ! environ output file has not been opened
    !
    IF ( PRESENT(local_verbose) ) THEN
       verbosity = verbose + local_verbose
    ELSE
       verbosity = verbose
    END IF
    !
    IF ( verbosity .EQ. 0 ) RETURN ! nothing to output
    !
    IF ( PRESENT(local_depth) ) THEN
       passed_verbosity = verbosity - verbose - local_depth
       passed_depth = local_depth
    ELSE
       passed_verbosity = verbosity - verbose
       passed_depth = depth
    END IF
    !
    IF ( verbosity .GE. 1 ) THEN
       IF ( verbosity .GE. verbose .AND. ionode ) WRITE( UNIT = environ_unit, FMT = 1700 )
       IF ( ionode ) WRITE( UNIT = environ_unit, FMT = 1701 )externals%number
       IF ( ionode ) WRITE( UNIT = environ_unit, FMT = 1702 )externals%charge
       IF ( verbosity .GE. 3 ) THEN
          CALL print_environ_functions(externals%number,externals%functions,&
               & passed_verbosity,passed_depth)
          CALL print_environ_density(externals%density,passed_verbosity,passed_depth)
       ENDIF
    END IF
    !
    FLUSH( environ_unit )
    !
    RETURN
    !
1700 FORMAT(/,4('%'),' EXTERNALS ',65('%'))
1701 FORMAT(1x,'number of external charges = ',I10)
1702 FORMAT(1x,'total external charge      = ',F14.7)
!--------------------------------------------------------------------
  END SUBROUTINE print_environ_externals
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE print_environ_charges( charges, local_verbose, local_depth )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    TYPE( environ_charges ), INTENT(IN) :: charges
    INTEGER, INTENT(IN), OPTIONAL :: local_verbose
    INTEGER, INTENT(IN), OPTIONAL :: local_depth
    !
    INTEGER :: verbosity, passed_verbosity, passed_depth
    !
    CHARACTER( LEN=80 ) :: sub_name = 'print_environ_charges'
    !
    IF ( verbose .EQ. 0 ) RETURN ! environ output file has not been opened
    !
    IF ( PRESENT(local_verbose) ) THEN
       verbosity = verbose + local_verbose
    ELSE
       verbosity = verbose
    END IF
    !
    IF ( verbosity .EQ. 0 ) RETURN ! nothing to output
    !
    IF ( PRESENT(local_depth) ) THEN
       passed_verbosity = verbosity - verbose - local_depth
       passed_depth = local_depth
    ELSE
       passed_verbosity = verbosity - verbose
       passed_depth = depth
    END IF
    !
    IF ( verbosity .GE. 1 ) THEN
       IF ( verbosity .GE. verbose .AND. ionode ) WRITE( UNIT = environ_unit, FMT = 1800 )
       IF ( ionode ) WRITE( UNIT = environ_unit, FMT = 1801 )charges%number
       IF ( ionode ) WRITE( UNIT = environ_unit, FMT = 1802 )charges%charge
       IF ( verbosity .GE. 3 ) THEN
          CALL print_environ_density(charges%density,passed_verbosity,passed_depth)
          IF ( charges % include_ions ) &
               & CALL print_environ_ions(charges%ions,passed_verbosity,passed_depth)
          IF ( charges % include_electrons ) &
               & CALL print_environ_electrons(charges%electrons,passed_verbosity,passed_depth)
          IF ( charges % include_externals ) &
               & CALL print_environ_externals(charges%externals,passed_verbosity,passed_depth)
          IF ( charges % include_dielectric ) &
               & CALL print_environ_dielectric(charges%dielectric,passed_verbosity,passed_depth)
          IF ( charges % include_electrolyte ) &
               & CALL print_environ_electrolyte(charges%electrolyte,passed_verbosity,passed_depth)
       ENDIF
    END IF
    !
    FLUSH( environ_unit )
    !
    RETURN
    !
1800 FORMAT(/,4('%'),' CHARGES ',67('%'))
1801 FORMAT(1x,'total number of charges    = ',I10)
1802 FORMAT(1x,'total charge               = ',F14.7)
!--------------------------------------------------------------------
  END SUBROUTINE print_environ_charges
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE print_environ_system( system, local_verbose, local_depth )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    TYPE( environ_system ), INTENT(IN) :: system
    INTEGER, INTENT(IN), OPTIONAL :: local_verbose
    INTEGER, INTENT(IN), OPTIONAL :: local_depth
    !
    INTEGER :: verbosity, passed_verbosity, passed_depth
    !
    CHARACTER( LEN=80 ) :: sub_name = 'print_environ_system'
    !
    IF ( verbose .EQ. 0 ) RETURN ! environ output file has not been opened
    !
    IF ( PRESENT(local_verbose) ) THEN
       verbosity = verbose + local_verbose
    ELSE
       verbosity = verbose
    END IF
    !
    IF ( verbosity .EQ. 0 ) RETURN ! nothing to output
    !
    IF ( PRESENT(local_depth) ) THEN
       passed_verbosity = verbosity - verbose - local_depth
       passed_depth = local_depth
    ELSE
       passed_verbosity = verbosity - verbose 
       passed_depth = depth
    END IF
    !
    IF ( verbosity .GE. 1 ) THEN
       IF ( verbosity .GE. verbose .AND. ionode ) WRITE( UNIT = environ_unit, FMT = 1900 )
       IF ( system%ntyp .EQ. 0 ) THEN
          IF ( ionode ) WRITE( UNIT = environ_unit, FMT = 1901 )
       ELSE
          IF ( ionode ) WRITE( UNIT = environ_unit, FMT = 1902 )system%ntyp
       END IF
       IF ( ionode ) WRITE( UNIT = environ_unit, FMT = 1903 )system%dim,system%axis
       IF ( ionode ) WRITE( UNIT = environ_unit, FMT = 1904 )system%pos,system%width
       IF ( verbosity .GE. 3 ) &
          & CALL print_environ_ions(system % ions, passed_verbosity, passed_depth )
    END IF
    !
    FLUSH( environ_unit )
    !
    RETURN
    !
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
    !
    IMPLICIT NONE
    !
    TYPE( environ_boundary ), INTENT(IN) :: boundary
    INTEGER, INTENT(IN), OPTIONAL :: local_verbose
    INTEGER, INTENT(IN), OPTIONAL :: local_depth
    !
    INTEGER :: verbosity, passed_verbosity, passed_depth
    !
    CHARACTER( LEN=80 ) :: sub_name = 'print_environ_boundary'
    !
    IF ( verbose .EQ. 0 ) RETURN ! environ output file has not been opened
    !
    IF ( PRESENT(local_verbose) ) THEN
       verbosity = verbose + local_verbose
    ELSE
       verbosity = verbose
    END IF
    !
    IF ( verbosity .EQ. 0 ) RETURN ! nothing to output
    !
    IF ( PRESENT(local_depth) ) THEN
       passed_verbosity = verbosity - verbose - local_depth
       passed_depth = local_depth
    ELSE
       passed_verbosity = verbosity - verbose
       passed_depth = depth
    END IF
    !
    IF ( verbosity .GE. 1 ) THEN
       IF ( verbosity .GE. verbose .AND. ionode ) WRITE( UNIT = environ_unit, FMT = 2000 )
       IF ( ionode ) WRITE( UNIT = environ_unit, FMT = 2001 )boundary%label
       IF ( ionode ) WRITE( UNIT = environ_unit, FMT = 2002 )boundary%mode
       IF ( boundary % need_electrons ) THEN
          IF ( ionode ) WRITE( UNIT = environ_unit, FMT = 2003 ) boundary % type
          SELECT CASE ( boundary % type )
          CASE ( 0 )
             IF ( ionode ) WRITE( UNIT = environ_unit, FMT = 2004 ) boundary % rhozero, boundary % tbeta
          CASE ( 1 )
             IF ( ionode ) WRITE( UNIT = environ_unit, FMT = 2005 ) boundary % rhomax, boundary % rhomin
             IF ( verbosity .GE. 3 .AND. ionode ) WRITE( UNIT = environ_unit, FMT = 2006 ) boundary % fact
          CASE ( 2 )
             IF ( ionode ) WRITE( UNIT = environ_unit, FMT = 2007 ) boundary % rhomax, boundary % rhomin
          END SELECT
          IF ( verbosity .GE. 4 ) THEN
             CALL print_environ_density(boundary%density,passed_verbosity,passed_depth)
             IF ( boundary % need_ions ) THEN
                IF ( ionode ) WRITE( UNIT = environ_unit, FMT = 2008 )
                IF ( verbosity .GE. 4 ) CALL print_environ_density(boundary%ions%core)
             END IF
             IF ( verbosity .GE. 4 ) &
                  & CALL print_environ_electrons(boundary%electrons,passed_verbosity,passed_depth)
          END IF
          IF ( verbosity .GE. 5 ) CALL print_environ_density(boundary%dscaled,passed_verbosity,passed_depth)
          IF ( verbosity .GE. 5 ) CALL print_environ_density(boundary%d2scaled,passed_verbosity,passed_depth)
       ELSE IF ( boundary % need_ions ) THEN
          IF ( ionode ) WRITE( UNIT = environ_unit, FMT = 2009 ) boundary%alpha, boundary%softness
          IF ( verbosity .GE. 3 ) &
               & CALL print_environ_functions(boundary%ions%number,boundary%soft_spheres,&
               & passed_verbosity,passed_depth)
       ELSE IF ( boundary % need_system ) THEN
          IF ( ionode ) WRITE( UNIT = environ_unit, FMT = 2010 ) boundary%simple%pos,boundary%simple%width,&
               & boundary%simple%spread,boundary%simple%dim,boundary%simple%axis
       END IF
       IF ( ionode ) WRITE( UNIT = environ_unit, FMT = 2011 ) boundary%volume
       IF ( boundary%deriv .GE. 1 .AND. ionode ) WRITE( UNIT = environ_unit, FMT = 2012 ) boundary%surface
       IF ( verbosity .GE. 4 ) CALL print_environ_density(boundary%scaled,passed_verbosity,passed_depth)
       IF ( boundary%solvent_aware ) THEN
          IF ( ionode ) WRITE( UNIT = environ_unit, FMT = 2013 ) boundary%filling_threshold,boundary%filling_spread,&
               & boundary%solvent_probe%width,boundary%solvent_probe%spread
          IF ( verbosity .GE. 4 ) CALL print_environ_density(boundary%local,passed_verbosity,passed_depth)
          IF ( verbosity .GE. 4 ) CALL print_environ_density(boundary%filling,passed_verbosity,passed_depth)
          IF ( verbosity .GE. 5 ) CALL print_environ_density(boundary%dfilling,passed_verbosity,passed_depth)
          IF ( verbosity .GE. 5 ) CALL print_environ_density(boundary%probe,passed_verbosity,passed_depth)
       END IF
       IF ( verbosity .GE. 5 .AND. boundary%deriv .GE. 1 ) &
            & CALL print_environ_gradient(boundary%gradient,passed_verbosity,passed_depth)
       IF ( verbosity .GE. 5 .AND. boundary%deriv .GE. 2 ) &
            & CALL print_environ_density(boundary%laplacian,passed_verbosity,passed_depth)
       IF ( verbosity .GE. 5 .AND. boundary%deriv .EQ. 3 ) &
            & CALL print_environ_density(boundary%dsurface,passed_verbosity,passed_depth)
    END IF
    !
    FLUSH( environ_unit )
    !
    RETURN
    !
2000 FORMAT(/,4('%'),' BOUNDARY ',66('%'))
2001 FORMAT(1x,'boundary label             = ',A20,' ')
2002 FORMAT(1x,'boundary mode              = ',A20,' ')
2003 FORMAT(1x,'boundary is built as a function of a smooth density'&
          /,1x,'function type              = ',I2,' ')
2004 FORMAT(1x,'using the Fattebert-Gygi function with parameters '&
          /,1x,'rhozero                    = ',F14.7,' '&
          /,1x,'2*beta                     = ',F14.7,' ')
2005 FORMAT(1x,'using the optimal SCCS function with parameters '&
          /,1x,'rhomax                     = ',F14.7,' '&
          /,1x,'rhomin                     = ',F14.7,' ')
2006 FORMAT(1x,'log(rhomax/rhomin)         = ',F14.7,' ')
2007 FORMAT(1x,'using the modified SCCS function with parameters '&
          /,1x,'rhomax                     = ',F14.7,' '&
          /,1x,'rhomin                     = ',F14.7,' ')
2008 FORMAT(1x,'adding fictitious core-electrons')
2009 FORMAT(1x,'boundary is built from soft-spheres centered on ionic positions'&
          /,1x,'solvent-dependent scaling  = ',F14.7,' '&
          /,1x,'softness parameter         = ',F14.7,' ')
2010 FORMAT(1x,'boundary is built as an analytic function centered on system position'&
          /,1x,'center of the boundary     = ',3F14.7,' '&
          /,1x,'distance from the center   = ',F14.7,' '&
          /,1x,'spread of the interface    = ',F14.7,' '&
          /,1x,'dimensionality             = ',I2,' '&
          /,1x,'axis                       = ',I2,' ')
2011 FORMAT(1x,'volume of the QM region    = ',F14.7,' ')
2012 FORMAT(1x,'surface of the QM region   = ',F14.7,' ')
2013 FORMAT(1x,'using solvent-aware boundary           '&
          /,1x,'filling threshold          = ',F14.7,' '&
          /,1x,'filling spread             = ',F14.7,' '&
          /,1x,'solvent radius x rad scale = ',F14.7,' '&
          /,1x,'spread of solvent probe    = ',F14.7,' ')
!--------------------------------------------------------------------
  END SUBROUTINE print_environ_boundary
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE print_environ_dielectric( dielectric, local_verbose, local_depth )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    TYPE( environ_dielectric ), INTENT(IN) :: dielectric
    INTEGER, INTENT(IN), OPTIONAL :: local_verbose
    INTEGER, INTENT(IN), OPTIONAL :: local_depth
    !
    INTEGER :: verbosity, passed_verbosity, passed_depth
    !
    CHARACTER( LEN=80 ) :: sub_name = 'print_environ_dielectric'
    !
    IF ( verbose .EQ. 0 ) RETURN ! environ output file has not been opened
    !
    IF ( PRESENT(local_verbose) ) THEN
       verbosity = verbose + local_verbose
    ELSE
       verbosity = verbose
    END IF
    !
    IF ( verbosity .EQ. 0 ) RETURN ! nothing to output
    !
    IF ( PRESENT(local_depth) ) THEN
       passed_verbosity = verbosity - verbose - local_depth
       passed_depth = local_depth
    ELSE
       passed_verbosity = verbosity - verbose
       passed_depth = depth
    END IF
    !
    IF ( verbosity .GE. 1 ) THEN
       IF ( verbosity .GE. verbose .AND. ionode ) WRITE( UNIT = environ_unit, FMT = 2100 )
       IF ( dielectric % nregions .EQ. 0 ) THEN
          IF ( ionode ) WRITE( UNIT = environ_unit, FMT = 2101 )dielectric%constant
       ELSE
          IF ( ionode ) WRITE( UNIT = environ_unit, FMT = 2102 )dielectric%constant,dielectric%nregions
          CALL print_environ_functions(dielectric%nregions,dielectric%regions,passed_verbosity,passed_depth)
          IF ( verbosity .GE. 4 ) CALL print_environ_density(dielectric%background,passed_verbosity,passed_depth)
       END IF
!       CALL print_environ_boundary(dielectric%boundary,passed_verbosity,passed_depth)
       IF ( verbosity .GE. 3 ) CALL print_environ_density(dielectric%density,passed_verbosity,passed_depth)
       IF ( verbosity .GE. 3 ) CALL print_environ_density(dielectric%epsilon,passed_verbosity,passed_depth)
       IF ( verbosity .GE. 5 ) CALL print_environ_density(dielectric%depsilon,passed_verbosity,passed_depth)
       IF ( ionode ) WRITE( UNIT = environ_unit, FMT = 2103 )dielectric%need_gradient,&
            & dielectric%need_factsqrt
       IF ( ionode ) WRITE( UNIT = environ_unit, FMT = 2104 )dielectric%charge
       IF ( verbosity .GE. 5 ) THEN
          CALL print_environ_gradient(dielectric%gradlog,passed_verbosity,passed_depth)
          IF ( dielectric%need_gradient ) CALL print_environ_gradient(dielectric%gradient,passed_verbosity,passed_depth)
          IF ( dielectric%need_factsqrt ) CALL print_environ_density(dielectric%factsqrt,passed_verbosity,passed_depth)
       END IF
    END IF
    !
    FLUSH( environ_unit )
    !
    RETURN
    !
2100 FORMAT(/,4('%'),' DIELECTRIC ',65('%'))
2101 FORMAT(1x,'dielectric build on homogeneous background'&
          /,1x,'environment bulk permitt.  = ',F14.7,' ')
2102 FORMAT(1x,'dielectric build in the presence of dielectric regions'&
          /,1x,'environment bulk permitt.  = ',F14.7,' '&
          /,1x,'number of dielec. regions  = ',I4,' ')
2103 FORMAT(1x,'dielectric flags'&
          /,1x,'need gradient              = ',L2,' '&
          /,1x,'need factor depend. sqrt   = ',L2,' ')
2104 FORMAT(1x,'total dielectric charge    = ',F14.7)
!--------------------------------------------------------------------
  END SUBROUTINE print_environ_dielectric
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE print_environ_electrolyte( electrolyte, local_verbose, local_depth )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    TYPE( environ_electrolyte ), INTENT(IN) :: electrolyte
    INTEGER, INTENT(IN), OPTIONAL :: local_verbose
    INTEGER, INTENT(IN), OPTIONAL :: local_depth
    !
    INTEGER :: verbosity, passed_verbosity, passed_depth, ityp
    !
    CHARACTER( LEN=80 ) :: sub_name = 'print_environ_electrolyte'
    !
    IF ( verbose .EQ. 0 ) RETURN ! environ output file has not been opened
    !
    IF ( PRESENT(local_verbose) ) THEN
       verbosity = verbose + local_verbose
    ELSE
       verbosity = verbose
    END IF
    !
    IF ( verbosity .EQ. 0 ) RETURN ! nothing to output
    !
    IF ( PRESENT(local_depth) ) THEN
       passed_verbosity = verbosity - verbose - local_depth
       passed_depth = local_depth
    ELSE
       passed_verbosity = verbosity - verbose
       passed_depth = depth
    END IF
    !
    IF ( verbosity .GE. 1 ) THEN
       IF ( verbosity .GE. verbose .AND. ionode ) WRITE( UNIT = environ_unit, FMT = 3100 )
       IF ( ionode ) WRITE( UNIT = environ_unit, FMT = 3101 )electrolyte%ntyp
       IF ( ionode ) WRITE( UNIT = environ_unit, FMT = 3102 )electrolyte%temperature
       IF ( ionode ) WRITE( UNIT = environ_unit, FMT = 3103 )1.D0/SQRT(electrolyte%k2)
       IF ( electrolyte%cionmax .GT. 0.D0 ) THEN
          IF ( ionode ) WRITE( UNIT = environ_unit, FMT = 3104 )electrolyte%cionmax
       ENDIF
       DO ityp=1,electrolyte%ntyp
          IF ( ionode ) WRITE( UNIT = environ_unit, FMT = 3105 )electrolyte%ioncctype(ityp)%index
          IF ( ionode ) WRITE( UNIT = environ_unit, FMT = 3106 )electrolyte%ioncctype(ityp)%cbulk
          IF ( ionode ) WRITE( UNIT = environ_unit, FMT = 3107 )electrolyte%ioncctype(ityp)%cbulk &
                                                                 & * amu_si / bohr_radius_si**3
          IF ( ionode ) WRITE( UNIT = environ_unit, FMT = 3108 )electrolyte%ioncctype(ityp)%z
          IF ( verbosity .GE. 5 ) CALL print_environ_density(electrolyte%ioncctype(ityp)%c, &
              & passed_verbosity,passed_depth)
          IF ( verbosity .GE. 5 ) CALL print_environ_density(electrolyte%ioncctype(ityp)%cfactor, &
              & passed_verbosity,passed_depth)
       END DO
!       CALL print_environ_boundary(electrolyte%boundary,passed_verbosity,passed_depth)
       IF ( verbosity .GE. 3 ) CALL print_environ_density(electrolyte%density,passed_verbosity,passed_depth)
       IF ( verbosity .GE. 3 ) CALL print_environ_density(electrolyte%gamma,passed_verbosity,passed_depth)
       IF ( verbosity .GE. 5 ) CALL print_environ_density(electrolyte%dgamma,passed_verbosity,passed_depth)
       IF ( ionode ) WRITE( UNIT = environ_unit, FMT = 3109 )electrolyte%linearized
       IF ( ionode ) WRITE( UNIT = environ_unit, FMT = 3110 )electrolyte%charge
    END IF
    !
    FLUSH( environ_unit )
    !
    RETURN
3100 FORMAT(/,4('%'),' ELECTROLYTE ',64('%'))
3101 FORMAT(1x,'number electrol. species   = ',I4)
3102 FORMAT(1x,'solvent temperature        = ',F7.1)
3103 FORMAT(1x,'Debye length / sqrt(eps)   = ',F14.7)
3104 FORMAT(1x,'modified Poisson-Boltzmann '&
          /,1x,'maximum concentration      = ',F14.7)
3105 FORMAT(1x,'electrolyte species ',I4)
3106 FORMAT(1x,'bulk concentration  (a.u.) = ',E15.4)
3107 FORMAT(1x,'bulk concentration (mol/L) = ',F14.7)
3108 FORMAT(1x,'ionic charge               = ',F7.2)
3109 FORMAT(1x,'electrolyte flags'&
          /,1x,'linearized                 = ',L2)
3110 FORMAT(1x,'total electrolyte charge   = ',F14.7)
!--------------------------------------------------------------------
  END SUBROUTINE print_environ_electrolyte
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE environ_print_energies( )
!--------------------------------------------------------------------
    !
    ! Write out the different Environ contributions to the energy.
    ! Called by electrons.f90
    !
    USE environ_base, ONLY : lelectrostatic, eelectrostatic, &
                             lsurface, esurface, &
                             lvolume, evolume, &
                             lelectrolyte, eelectrolyte, &
                             lconfine, econfine, &
                             deenviron
    !
    CHARACTER( LEN=80 ) :: sub_name = 'environ_print_energies'
    !
    IF ( ionode ) THEN
       IF ( prog .EQ. 'PW' ) THEN
          IF ( lelectrostatic ) WRITE( program_unit, 9201 ) eelectrostatic
          IF ( lsurface ) WRITE( program_unit, 9202 ) esurface
          IF ( lvolume ) WRITE( program_unit, 9203 ) evolume
          IF ( lelectrolyte ) WRITE( program_unit, 9205 ) eelectrolyte
          IF ( lconfine ) WRITE( program_unit, 9206 ) econfine
          WRITE( program_unit, 9204 ) deenviron
       ELSE IF ( prog .EQ. 'CP' ) THEN
          IF ( lelectrostatic ) WRITE( program_unit, 9301 ) eelectrostatic
          IF ( lsurface ) WRITE( program_unit, 9302 ) esurface
          IF ( lvolume ) WRITE( program_unit, 9303 ) evolume
          IF ( lelectrolyte ) WRITE( program_unit, 9305 ) eelectrolyte
          IF ( lconfine ) WRITE( program_unit, 9306 ) econfine
          WRITE( program_unit, 9304 ) deenviron
       ELSE
          CALL errore(sub_name,'Wrong program calling Environ',1)
       END IF
    END IF
    !
    RETURN
    !
9201 FORMAT( '     electrostatic embedding   =',F17.8,' Ry')
9202 FORMAT( '     cavitation energy         =',F17.8,' Ry')
9203 FORMAT( '     PV energy                 =',F17.8,' Ry')
9206 FORMAT( '     confinement energy        =',F17.8,' Ry')
9205 FORMAT( '     electrolyte free energy   =',F17.8,' Ry')
9204 FORMAT( '     correction to one-el term =',F17.8,' Ry')
9301 FORMAT( '     electrostatic embedding = ',F14.5,' Hartree a.u.')
9302 FORMAT( '           cavitation energy = ',F14.5,' Hartree a.u.')
9303 FORMAT( '                   PV energy = ',F14.5,' Hartree a.u.')
9305 FORMAT( '     electrolyte free energy = ',F14.5,' Hartree a.u.')
9306 FORMAT( '          confinement energy = ',F14.5,' Hartree a.u.')
9304 FORMAT( '   correction to one-el term = ',F14.5,' Hartree a.u.')
    !
!----------------------------------------------------------------------
  END SUBROUTINE environ_print_energies
!----------------------------------------------------------------------
!----------------------------------------------------------------------
  SUBROUTINE environ_print_potential_shift( )
!----------------------------------------------------------------------
    !
    ! If Gaussian nuclei are used, write out the corresponding
    ! potential shift
    !
    USE environ_base, ONLY : lsmearedions, potential_shift

    IF (lsmearedions) WRITE( program_unit, 9400 ) potential_shift * rytoev

9400 FORMAT(/,5(' '),&
          'the potential shift due to the Gaussian-smeared nuclei is ',&
          F10.4,' ev')
    !
!----------------------------------------------------------------------
  END SUBROUTINE environ_print_potential_shift
!----------------------------------------------------------------------
  SUBROUTINE environ_print_potential_warning( )
!----------------------------------------------------------------------
    !
    USE electrostatic_base, ONLY : need_pbc_correction

    IF (need_pbc_correction) WRITE( program_unit, 9401 ) 

9401 FORMAT(/,&
     5(' '),'WARNING: you are using the parabolic pbc correction;',/,&
     5(' '),'         the potential shift above must be added to ',/,&
     5(' '),'         band and Fermi energies.')
    !
!----------------------------------------------------------------------
  END SUBROUTINE environ_print_potential_warning
!----------------------------------------------------------------------
!  Subroutine: environ_summary
!
!> Write out the main parameters of Environ calculations, summarizing
!! the input keywords (some info also on internal vs input units).
!! Called by summary.f90
!----------------------------------------------------------------------
  SUBROUTINE environ_summary( )
!----------------------------------------------------------------------
    USE environ_base,       ONLY : environ_thr, lsolvent, solvent,   &
                                   env_static_permittivity,          &
                                   env_optical_permittivity,         &
                                   env_surface_tension,              &
                                   env_pressure, lelectrostatic,     &
                                   ltddfpt
    USE electrostatic_base, ONLY : outer, boundary_core, lfd, fd
    !
    IMPLICIT NONE
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
             IF ( solvent%solvent_aware ) WRITE( UNIT = program_unit, FMT = 9013 )
             IF ( solvent%field_aware ) THEN
                WRITE( UNIT = program_unit, FMT = 9014 )
                WRITE( UNIT = program_unit, FMT = 9015 )solvent%field_factor,solvent%charge_asymmetry
                WRITE( UNIT = program_unit, FMT = 9016 )solvent%field_min,solvent%field_max
             ENDIF
          ENDIF
          !
          IF ( env_static_permittivity .GT. 1.D0 ) THEN
             WRITE( UNIT = program_unit, FMT = 9005 ) env_static_permittivity
             IF (ltddfpt) &
                  & WRITE( UNIT = program_unit, FMT = 9006 ) env_optical_permittivity
             WRITE( UNIT = program_unit, FMT = 9007 ) TRIM( solvent%mode )
          END IF
          !
          IF ( env_surface_tension .GT. 0.D0 ) WRITE( UNIT = program_unit, FMT = 9010 )      &
               env_surface_tension/1.D-3/bohr_radius_si**2*rydberg_si, env_surface_tension
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
             WRITE( UNIT = program_unit, FMT = 9101 )outer%problem, outer%solver%type, outer%solver%auxiliary
             !
             WRITE( UNIT = program_unit, FMT = 9102 )outer%core%type, boundary_core
             !
             IF ( lfd ) THEN
                IF ( fd%ifdtype .EQ. 1 ) THEN
                   WRITE( UNIT = program_unit, FMT = 9103 ) 'central diff.',fd%nfdpoint
                ELSE IF (fd%ifdtype .EQ. 2 .OR. fd%ifdtype .EQ. 3 ) THEN
                   WRITE( UNIT = program_unit, FMT = 9103 ) 'lanczos diff.',fd%nfdpoint
                ELSE IF (fd%ifdtype .EQ.4 .OR. fd%ifdtype .EQ. 5 ) THEN
                   WRITE( UNIT = program_unit, FMT = 9103 ) 'noise-robust diff.',fd%nfdpoint
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
8000 FORMAT(/,5x,'Plese cite',/,9x,A80,&
            /,5x,'in publications or presentations arising from this work.',/)
8001 FORMAT(/)
9000 FORMAT(/,5x,'Environ Module',/,5x,'==============')
9001 FORMAT( '     compensation onset threshold      = ',  E24.4,' ' )
9002 FORMAT( '     switching function adopted        = ',  A24,' ' )
9003 FORMAT( '     solvation density threshold       = ',  E24.4,' ' &
            /'     smoothness exponent (2 x beta)    = ',  F24.2,' ' )
9004 FORMAT( '     density limit for vacuum region   = ',  E24.4,' ' &
            /'     density limit for bulk solvent    = ',  E24.4,' ' )
9005 FORMAT( '     static permittivity               = ',  F24.2,' ' )
9006 FORMAT( '     optical permittivity              = ',  F24.4,' ' )
9007 FORMAT( '     epsilon calculation mode          = ',  A24,' ' )
9010 FORMAT( '     surface tension in input (dyn/cm) = ',  F24.2,' ' &
            /'     surface tension in internal units = ',  E24.4,' ' )
9011 FORMAT( '     external pressure in input (GPa)  = ',  F24.2,' ' &
            /'     external pressure in inter. units = ',  E24.4,' ' )
9012 FORMAT( '     correction slab geom. along axis  = ',  I24,' ' )
9013 FORMAT( '     interface is solvent aware            ' )
9014 FORMAT( '     interface is field aware            ' )
9015 FORMAT( '     field aware factor                = ', F24.2,' ' &
            /'     asymmetry of field-awareness      = ', F24.2,' ' )
9016 FORMAT( '     field limit for no correction     = ', F24.2,' ' &
            /'     field limit for full correction   = ', F24.2,' ' )
9100 FORMAT(/,5x,'Electrostatic Setup',/,5x,'-------------------')
9101 FORMAT( '     electrostatic problem to solve    = ',  A24,' ' &
            /'     numerical solver adopted          = ',  A24,' ' &
            /'     type of auxiliary density adopted = ',  A24,' ' )
9102 FORMAT( '     type of core tool for poisson     = ',  A24,' ' &
            /'     type of core tool for s(r) deriv  = ',  A24,' ' )
9103 FORMAT( '     type of numerical differentiator  = ',  A24,' ' &
            /'     number of points in num. diff.    = ',  I24,' ' )
!--------------------------------------------------------------------
  END SUBROUTINE environ_summary
!--------------------------------------------------------------------
!  Subroutine: environ_clock
!
!> Writes out the time informations of the Environ dependent
!! calculations. Called by print_clock_pw.f90
!--------------------------------------------------------------------
  SUBROUTINE environ_clock( passed_unit )
!--------------------------------------------------------------------
    USE environ_base,   ONLY : lelectrostatic, lsurface, lvolume, ltddfpt
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN), OPTIONAL :: passed_unit
    !
    INTEGER :: actual_unit
    !
    IF ( PRESENT( passed_unit ) ) THEN
       actual_unit = passed_unit
    ELSE
       actual_unit = program_unit
    ENDIF
    WRITE( actual_unit, * )
    WRITE( actual_unit, '(5X,"Environ routines")' )
    ! dielectric subroutines
    IF ( lelectrostatic ) THEN
       CALL print_clock ('calc_eelect')
       CALL print_clock ('calc_velect')
       CALL print_clock ('calc_vgcs')
       CALL print_clock ('dielectric')
       CALL print_clock ('electrolyte')
       CALL print_clock ('calc_felect')
    END IF
    ! TDDFT
    IF ( ltddfpt ) CALL print_clock ('calc_vsolvent_tddfpt')
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE environ_clock
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE write_cube( f, ions, idx, label )
!--------------------------------------------------------------------
    !
    USE fft_base,       ONLY : dfftp
! BACKWARD COMPATIBILITY
! Compatible with QE-5.1.X
!      USE fft_base,       ONLY : grid_gather
! Compatible with QE-5.2.X
!      USE fft_base,       ONLY : gather_grid
! Compatible with QE-5.3.X QE-5.4.X QE-6.X.X QE-GIT
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
    INTEGER, OPTIONAL :: idx
    CHARACTER( LEN=100 ), OPTIONAL :: label
    CHARACTER( LEN=100 ) :: filemod
    !
    nr1x = f%cell%n1x
    nr2x = f%cell%n2x
    nr3x = f%cell%n3x
    !
    nr1 = f%cell%n1
    nr2 = f%cell%n2
    nr3 = f%cell%n3
    !
    IF(PRESENT(idx)) THEN
       WRITE(filemod, '(i4.4)') idx
    ELSE
       filemod = ""
    ENDIF
    !
    IF(PRESENT(label)) THEN
       filename = TRIM(ADJUSTL(label))//TRIM(filemod)//".cube"
    ELSE
       filename = TRIM(ADJUSTL(f%label))//TRIM(filemod)//".cube"
    ENDIF
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
!Compatible with QE-5.1.X
!      CALL grid_gather( f, flocal )
!Compatible with QE-5.2.X QE-5.3.X QE-5.4.X QE-6.X.X QE-GIT
    CALL gather_grid( dfftp, f%of_r, flocal )
    CALL mp_sum( flocal, f%cell%comm )
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
!----------------------------------------------------------------------------
END MODULE environ_output
!----------------------------------------------------------------------------
