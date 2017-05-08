!
! Copyright (C) 2007-2008 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Subroutines to initialize and clean up the global ENVIRON variables
! Only variables that are required outside of the specific modules
! are initiliazed here, module-specific variables should be initialized
! inside each module by the specific initiliazation subruotine
!
! original version by O. Andreussii, I. Dabo and N. Marzari (MIT)
!
MODULE environ_init
  !
  USE environ_types
  !
  PRIVATE
  !
  PUBLIC :: environ_initbase, environ_initcell, environ_initions, &
     & environ_initelectrons, environ_initpotential, environ_clean
  !
CONTAINS
!
!--------------------------------------------------------------------
  SUBROUTINE environ_initbase( n1, n2, n3, ibrav, alat, omega, at, &
                             & nnr, ir_end, comm, me, root, e2 )
!--------------------------------------------------------------------
!
! Subroutine to initialize fundamental quantities needed by the
! environ modules. This subroutine is called by init_run.f90, thus
! only once per pw.x execution.
!
      ! ... Declares modules
      !
      ! In environ_base all the control flags plus the variables
      ! that need to be initialized
      !
      USE kinds,        ONLY : DP
      USE environ_base, ONLY : verbose, environ_unit, e2_ => e2,        &
                               cell, electrons, charges,                &
                               vzero, deenviron,                        &
                               lelectrostatic, eelectrostatic,          &
                               velectrostatic, vsoftcavity,             &
                               lelectrolyte, electrolyte,               &
                               lsolvent, solvent, lstatic, static,      &
                               loptical, optical,                       &
                               lexternals, externals,                   &
                               lsurface, ecavity, vcavity,              &
                               lvolume, epressure, vpressure
      !
      ! Local base initialization subroutines for the different
      ! environ contributions
      !
      IMPLICIT NONE
      !
      ! ... Input variable
      !
      INTEGER, INTENT(IN) :: nnr
      INTEGER, INTENT(IN) :: ir_end
      INTEGER, INTENT(IN) :: n1, n2, n3
      INTEGER, INTENT(IN) :: ibrav
      INTEGER, INTENT(IN) :: comm
      INTEGER, INTENT(IN) :: me
      INTEGER, INTENT(IN) :: root
      REAL(DP), INTENT(IN) :: alat
      REAL(DP), INTENT(IN) :: omega
      REAL(DP), DIMENSION(3,3), INTENT(IN) :: at
      REAL(DP), OPTIONAL, INTENT(IN) :: e2
      !
      ! ... Common initialization for simulations with Environ
      !
      IF ( verbose .GE. 1 ) &
        OPEN(unit=environ_unit,file='environ.debug',status='unknown')
      !
      e2_ = 2.D0
      IF ( PRESENT(e2) ) e2_ = e2
      !
      CALL init_environ_cell( n1, n2, n3, ibrav, alat, omega, at, nnr, ir_end, comm, me, root, cell )
      !
      ! ... Create local storage for base potential, that needs to be modified
      !
      CALL create_environ_density( vzero )
      CALL init_environ_density( cell, vzero )
      !
      deenviron = 0.0_DP
      !
      ! ... Electrostatic contribution
      !
      eelectrostatic  = 0.0_DP
      IF ( lelectrostatic ) THEN
         !
         CALL create_environ_density( velectrostatic )
         CALL init_environ_density( cell, velectrostatic )
         !
         CALL create_environ_density( vsoftcavity )
         CALL init_environ_density( cell, vsoftcavity )
         !
      END IF
      !
      ! ... Cavity contribution
      !
      ecavity  = 0.0_DP
      IF ( lsurface ) THEN
         !
         CALL create_environ_density( vcavity )
         CALL init_environ_density( cell, vcavity )
         !
      END IF
      !
      ! ... Pressure contribution
      !
      epressure  = 0.0_DP
      IF ( lvolume ) THEN
         !
         CALL create_environ_density( vpressure )
         CALL init_environ_density( cell, vpressure )
         !
      ENDIF
      !
      ! ... Second step of initialization of some environ derived type
      !
      CALL init_environ_electrons_second( cell, electrons )
      !
      IF ( lsolvent ) CALL init_environ_boundary_second( cell, solvent )
      !
      IF ( lstatic ) CALL init_environ_dielectric_second( cell, static )
      !
      IF ( loptical ) CALL init_environ_dielectric_second( cell, optical )
      !
      IF ( lelectrolyte ) CALL init_environ_electrolyte_second( cell, electrolyte )
      !
      IF ( lexternals ) CALL init_environ_externals_second( cell, externals )
      !
      IF ( lelectrostatic ) CALL init_environ_charges_second( cell, charges )
      !
      RETURN
      !
!--------------------------------------------------------------------
   END SUBROUTINE environ_initbase
!--------------------------------------------------------------------
!--------------------------------------------------------------------
      SUBROUTINE environ_initpotential( nnr, vltot )
!--------------------------------------------------------------------
!
! Save local potential that will be overwritten by environ
!
      ! ... Declares modules
      USE kinds,         ONLY : DP
      USE environ_base,  ONLY : verbose, environ_unit, vzero
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN) :: nnr
      REAL( DP ), INTENT( IN ) :: vltot(nnr)
      !
      CHARACTER( LEN=80 ) :: sub_name = 'environ_initpotential'
      !
      vzero % update = .TRUE.
      !
      IF ( vzero % cell % nnr .NE. nnr ) &
           & CALL errore(sub_name,'Inconsistent size in input potential',1)
      vzero % of_r = vltot
      !
      RETURN
      !
!--------------------------------------------------------------------
    END SUBROUTINE environ_initpotential
!--------------------------------------------------------------------
!--------------------------------------------------------------------
      SUBROUTINE environ_initcell( omega, at )
!--------------------------------------------------------------------
!
! Initialize the cell-related quantities to be used in the Environ
! modules. This initialization is called by electrons.f90, thus it
! is performed at every step of electronic optimization.
!
      ! ... Declares modules
      USE kinds,         ONLY : DP
      USE environ_base,  ONLY : verbose, environ_unit, cell,        &
                                lstatic, loptical, static, optical, &
                                lexternals, externals
      ! ... Cell-related updates
      USE dielectric,     ONLY : update_environ_dielectric
!      USE extcharges,    ONLY : update_external_charges
      !
      IMPLICIT NONE
      !
      REAL( DP ), INTENT( IN ) :: omega
      REAL( DP ), INTENT( IN ) :: at(3,3)
      !
      cell%update = .TRUE.
      !
      ! ... Update cell parameters
      !
      CALL update_environ_cell( omega, at, cell )
      !
      ! ... Update fixed quantities defined inside the cell
      !
      IF ( lstatic ) CALL update_environ_dielectric( static )
      IF ( loptical ) CALL update_environ_dielectric( optical )
      !
!      IF ( lexternals ) CALL update_environ_externals( externals )
      !
      cell%update = .FALSE.
      !
      RETURN
      !
!--------------------------------------------------------------------
   END SUBROUTINE environ_initcell
!--------------------------------------------------------------------
!--------------------------------------------------------------------
      SUBROUTINE environ_initions( nnr, nat, ntyp, ityp, zv, tau )
!--------------------------------------------------------------------
!
! Initialize the ions-related quantities to be used in the Environ
! modules. This initialization is called by electrons.f90, thus it
! is performed at every step of electronic optimization. It may not
! be the most efficient choice, but it is a safe choice.
!
      ! ... Declares modules
      USE kinds,             ONLY : DP
      USE environ_base,      ONLY : verbose, environ_unit, cell, &
                                    ions, electrons, system,     &
                                    lsolvent, solvent,           &
                                    lstatic, static,             &
                                    loptical, optical,           &
                                    lelectrolyte, electrolyte,    &
                                    lrigidcavity
      USE boundary,          ONLY : update_environ_boundary
      USE dielectric,        ONLY : update_environ_dielectric
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT( IN )     :: nnr, nat, ntyp
      INTEGER, INTENT( IN )     :: ityp( nat )
      REAL ( DP ), INTENT( IN ) :: zv( ntyp )
      REAL ( DP ), INTENT( IN ) :: tau( 3, nat )
      !
      INTEGER :: ia, dim, axis, icor, max_ntyp
      REAL ( DP ) :: charge, spread, dist, pos(3)
      !
      ions%update = .TRUE.
      !
      ! ... Second step of initialization, need to be moved out of here
      !
      CALL init_environ_ions_second( nat, ntyp, ityp, zv, cell, ions )
      !
      ! ... Update ions parameters
      !
      CALL update_environ_ions( nat, tau, ions )
      !
      ! ... Update system parameters
      !
      CALL update_environ_system( system )
      !
      ! ... Update rigid environ properties, defined on ions
      !
      IF ( lrigidcavity ) THEN
         !
         IF ( lsolvent ) THEN
            !
            CALL update_environ_boundary( solvent )
            !
            ! ... Update quantities that depend on the solvent boundary
            !
            IF ( lstatic ) CALL update_environ_dielectric( static )
            IF ( loptical ) CALL update_environ_dielectric( optical )
            !
         ENDIF
         !
         IF ( lelectrolyte ) CALL update_environ_boundary( electrolyte%boundary )
         !
      END IF
      !
      ions%update = .FALSE.
      !
      RETURN
      !
!--------------------------------------------------------------------
   END SUBROUTINE environ_initions
!--------------------------------------------------------------------
!--------------------------------------------------------------------
   SUBROUTINE environ_initelectrons( nelec, nspin, nnr, rho )
!--------------------------------------------------------------------
!
! Initialize the electrons-related quantities to be used in the Environ
! modules. This initialization is called by electrons.f90, thus it
! is performed at every step of electronic optimization.
!
      ! ... Declares modules
      USE kinds,             ONLY : DP
      USE environ_base,      ONLY : electrons, lsolvent, solvent, &
                                    lstatic, static,              &
                                    loptical, optical,            &
                                    lelectrolyte, electrolyte,    &
                                    lsoftcavity
      USE boundary,          ONLY : update_environ_boundary
      USE dielectric,        ONLY : update_environ_dielectric
      !
      IMPLICIT NONE
      !
      REAL( DP ), INTENT( IN )  :: nelec
      INTEGER, INTENT( IN )     :: nspin, nnr
      REAL ( DP ), INTENT( IN ) :: rho( nnr, nspin )
      !
      electrons%update = .TRUE.
      !
      ! ... Update electrons parameters
      !
      CALL update_environ_electrons( nelec, nspin, nnr, rho, electrons )
      !
      ! ... Update soft environ properties, defined on electrons
      !
      IF ( lsoftcavity ) THEN
         !
         IF ( lsolvent ) THEN
            !
            CALL update_environ_boundary( solvent )
            !
            ! ... Update quantities that depend on the solvent boundary
            !
            IF ( lstatic ) CALL update_environ_dielectric( static )
            IF ( loptical ) CALL update_environ_dielectric( optical )
            !
         ENDIF
         !
         IF ( lelectrolyte ) CALL update_environ_boundary( electrolyte%boundary )
         !
      END IF
      !
      electrons%update = .FALSE.
      !
      RETURN
      !
!--------------------------------------------------------------------
   END SUBROUTINE environ_initelectrons
!--------------------------------------------------------------------
!--------------------------------------------------------------------
      SUBROUTINE environ_clean(lflag)
!--------------------------------------------------------------------
!
! Clean up all the Environ related allocated variables, and call
! clean up subroutines of specific Environ modules. The structure of
! this subroutine should mirror the one of the environ_init...
! subroutines above. This clean up is called by clean_pw.f90
!
      ! ... Declares modules
      !
      ! In environ_base all the control flags plus the variables
      ! that need to be deallocated
      !
      USE environ_base, ONLY : environ_unit, vzero,                  &
                               lelectrostatic, velectrostatic,       &
                               vsoftcavity,                          &
                               lsurface, vcavity,                    &
                               lvolume, vpressure,                   &
                               ions, electrons, system, charges,     &
                               lexternals, externals,                &
                               lstatic, static, loptical, optical,   &
                               lsolvent, solvent,                    &
                               lelectrolyte, electrolyte
      !
      ! Local clean up subroutines for the different contributions
      !
      USE control_flags,  ONLY : tddfpt
      !
      IMPLICIT NONE
      !
      LOGICAL, INTENT(IN) :: lflag
      LOGICAL :: opnd
      !
      ! ... Deallocate environment variables
      !
      INQUIRE( unit=environ_unit, opened= opnd )
      IF ( opnd ) CLOSE(unit=environ_unit)
      !
      CALL destroy_environ_density( vzero )
      !
      ! ... environ_base variables
      !
      IF ( lelectrostatic ) THEN
         CALL destroy_environ_density( velectrostatic )
         CALL destroy_environ_density( vsoftcavity )
      END IF
      IF ( lsurface ) THEN
         CALL destroy_environ_density( vcavity )
      END IF
      IF ( lvolume ) THEN
         CALL destroy_environ_density( vpressure )
      END IF
      !
      ! ... destroy derived types which were allocated in input
      !
      IF ( lelectrostatic ) CALL destroy_environ_charges( lflag, charges )
      IF ( lexternals ) CALL destroy_environ_externals( lflag, externals )
      IF ( lstatic ) CALL destroy_environ_dielectric( lflag, static )
      IF ( loptical ) CALL destroy_environ_dielectric( lflag, optical )
      IF ( lsolvent ) CALL destroy_environ_boundary( lflag, solvent )
      IF ( lelectrolyte ) CALL destroy_environ_electrolyte( lflag, electrolyte )
      !
      CALL destroy_environ_electrons( lflag, electrons )
      CALL destroy_environ_ions( lflag, ions )
      CALL destroy_environ_system( lflag, system )
      !
!!!! DOUBLE CHECK TDDFPT      IF ( .NOT. tddfpt) THEN
!!!! DOUBLE CHECK TDDFPT        IF ( ALLOCATED( rhoions ) ) DEALLOCATE( rhoions )
!!!! DOUBLE CHECK TDDFPT      END IF
      !
      RETURN
      !
!--------------------------------------------------------------------
   END SUBROUTINE environ_clean
!--------------------------------------------------------------------
END MODULE environ_init
