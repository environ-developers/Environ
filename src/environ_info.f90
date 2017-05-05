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
MODULE environ_info
!--------------------------------------------------------------------

PRIVATE

PUBLIC :: environ_print_energies, environ_summary, environ_clock

CONTAINS
!--------------------------------------------------------------------
      SUBROUTINE environ_print_energies( )
!--------------------------------------------------------------------
      !
      ! Write out the different Environ contributions to the energy.
      ! Called by electrons.f90
      !
      USE io_global,    ONLY : stdout
      USE environ_base, ONLY : e2
      USE environ_base, ONLY : lelectrostatic, eelectrostatic, &
                               lcavity, ecavity, &
                               lvolume, evolume
      !
      IF ( e2 .EQ. 2.D0 ) THEN
        IF ( lelectrostatic ) WRITE( stdout, 9201 ) eelectrostatic
        IF ( lcavity ) WRITE( stdout, 9202 ) ecavity
        IF ( lvolume ) WRITE( stdout, 9203 ) epressure
      ELSE IF ( e2 .EQ. 1.D0 ) THEN
        IF ( lelectrostatic ) WRITE( stdout, 9301 ) eelectrostatic
        IF ( lcavity ) WRITE( stdout, 9302 ) ecavity
        IF ( lvolume ) WRITE( stdout, 9303 ) epressure
      ELSE
        WRITE(stdout,*)'ERROR: wrong value of e2 in Environ'
        STOP
      ENDIF
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
      USE io_global,        ONLY : stdout, ionode
      USE constants,        ONLY : rydberg_si, bohr_radius_si
      USE environ_base,     ONLY : environ_thr, solvent,                &
                                   env_static_permittivity,             &
                                   env_optical_permittivity,            &
                                   env_surface_tension,                 &
                                   env_pressure
      USE control_flags,     ONLY : tddfpt
      !
      IMPLICIT NONE
      !
      !
      IF( ionode ) THEN
        !
        WRITE( stdout, * )
        !
        WRITE( UNIT = stdout,                                          &
               FMT = '(/,5x, "Environ Module",                         &
                      &/,5x, "==============")' )
        WRITE( stdout, '(/5X,"Please cite",&
         &/9X,"""O. Andreussi, I. Dabo and N. Marzari, J. Chem. Phys. 136, ",&
         &    "064102 (2012);""", &
         &/5X,"in publications or presentations arising from this work.",&
         &/)' )
        !
        WRITE( UNIT = stdout, FMT = 9001 ) environ_thr
        !
        IF ( solvent%type .EQ. 0 ) THEN
          WRITE( UNIT = stdout, FMT = 9002 ) 'Fatteber-Gygi'
          WRITE( UNIT = stdout, FMT = 9003 ) solvent%rhozero, solvent%tbeta
        ELSE IF ( solvent%type .EQ. 1 ) THEN
          WRITE( UNIT = stdout, FMT = 9002 ) 'SCCS'
          WRITE( UNIT = stdout, FMT = 9004 ) solvent%rhomax, solvent%rhomin
        ENDIF
        !
        IF ( env_static_permittivity .GT. 1.D0 ) THEN
           WRITE( UNIT = stdout, FMT = 9005 ) env_static_permittivity
           IF (tddfpt) &
         & WRITE( UNIT = stdout, FMT = 9006 ) env_optical_permittivity
           WRITE( UNIT = stdout, FMT = 9007 ) TRIM( solvent%mode )
        END IF
        !
        IF ( env_surface_tension .GT. 0.D0 ) WRITE( UNIT = stdout, FMT = 9010 )      &
           env_surface_tension/1.D-3/bohr_radius_si**2*rydberg_si, env_surface_tension, solvent%deltatheta
        !
        IF ( env_pressure .NE. 0.D0 ) WRITE( UNIT = stdout, FMT = 9011 )&
           env_pressure*rydberg_si/bohr_radius_si**3*1.D-9, env_pressure
        !
        WRITE( stdout, * )
        !
      END IF
      !
9001 FORMAT( '     compensation onset threshold      = ',  E24.4,' ' )
9002 FORMAT( '     switching function adopted        = ',  A24,' ' )
9003 FORMAT( '     solvation density threshold       = ',  E24.4,' ' &
            /'     smoothness exponent (2 x beta)    = ',  F24.2,' ')
9004 FORMAT( '     density limit for vacuum region   = ',  E24.4,' ' &
            /'     density limit for bulk solvent    = ',  E24.4,' ')
9005 FORMAT( '     static permittivity               = ',  F24.2,' ')
9006 FORMAT( '     optical permittivity              = ',  F24.4,' ')
9007 FORMAT( '     epsilon calculation mode          = ',  A24,' ' )
9008 FORMAT( '     type of numerical differentiator  = ',  A24,' ' &
            /'     number of points in num. diff.    = ',  I24,' ' )
9009 FORMAT( '     required accuracy                 = ',  E24.4,' ' &
            /'     linear mixing parameter           = ',  F24.2,' ' )
9010 FORMAT( '     surface tension in input (dyn/cm) = ',  F24.2,' ' &
            /'     surface tension in internal units = ',  E24.4,' ' &
            /'     delta parameter for surface depth = ',  E24.4,' ' )
9011 FORMAT( '     external pressure in input (GPa)  = ',  F24.2,' ' &
            /'     external pressure in inter. units = ',  E24.4,' ' )
9012 FORMAT( '     correction slab geom. along axis  = ',  I24,' ' )
9013 FORMAT( '     number of external charged items  = ',  I24,' ' )
9014 FORMAT( '     number of dielectric regions      = ',  I24,' ' )

!--------------------------------------------------------------------
      END SUBROUTINE environ_summary
!--------------------------------------------------------------------
!--------------------------------------------------------------------
      SUBROUTINE environ_clock( stdout )
!--------------------------------------------------------------------
      !
      ! Write out the time informations of the Environ dependent
      ! calculations. Called by print_clock_pw.f90
      !
      USE environ_base,   ONLY : env_static_permittivity, &
                                 env_surface_tension,     &
                                 env_pressure,            &
                                 env_periodicity
      USE control_flags,  ONLY : tddfpt
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN) :: stdout
      !
      WRITE( stdout, * )
      WRITE( stdout, '(5X,"Environ routines")' )
      ! dielectric subroutines
      IF ( env_static_permittivity .GT. 1.D0 &
           .OR. env_dielectric_regions .GT. 0 ) THEN
         CALL print_clock ('calc_eelect')
         CALL print_clock ('calc_velect')
         CALL print_clock ('dielectric')
         CALL print_clock ('calc_felect')
      END IF
      ! cavitation subroutines
      IF ( env_surface_tension .GT. 0.D0 ) THEN
         CALL print_clock ('calc_ecav')
         CALL print_clock ('calc_vcav')
      END IF
      ! pressure subroutines
      IF ( env_pressure .NE. 0.D0 ) THEN
         CALL print_clock ('calc_epre')
         CALL print_clock ('calc_vpre')
      END IF
      ! periodic subroutines
      IF ( env_periodicity .NE. 3 ) THEN
         CALL print_clock ('calc_epbc')
         CALL print_clock ('calc_vpbc')
         CALL print_clock ('calc_fpbc')
      END IF
      ! TDDFT
      IF (tddfpt) CALL print_clock ('calc_vsolvent_tddfpt')
      !
      RETURN
      !
!--------------------------------------------------------------------
      END SUBROUTINE environ_clock
!--------------------------------------------------------------------
!
!--------------------------------------------------------------------
END MODULE environ_info
!--------------------------------------------------------------------
