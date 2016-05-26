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
PRIVATE
!
PUBLIC :: environ_initbase, environ_initcell, environ_initions_allocate, &
          environ_initions, environ_clean

CONTAINS
!
!--------------------------------------------------------------------
      SUBROUTINE environ_initbase( nnr, e2 )
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
                               use_smeared_ions, vltot_zero, deenviron, &
                               env_static_permittivity, esolvent,       &
                               vsolvent, vepsilon,                      &    
                               rhopol, ifdtype, nfdpoint,               &
                               env_surface_tension, ecavity, vcavity,   &
                               env_pressure, epressure, vpressure,      &
                               env_periodicity, eperiodic, vperiodic,   &
                               env_ioncc_level, eioncc, vioncc,         &
                               rhoioncc, rhopolcc,                      &
                               env_external_charges, eextcharge,        &
                               vextcharge, rhoexternal,                 &
                               env_dielectric_regions, epsstatic,       &
                               epsoptical   
      !
      ! Local base initialization subroutines for the different
      ! environ contributions
      !
      USE solvent,      ONLY : solvent_initbase
      USE periodic,     ONLY : periodic_initbase
      USE cavity,       ONLY : cavity_initbase
      USE ioncc,        ONLY : ioncc_initbase
      !
      IMPLICIT NONE      
      ! 
      ! ... Input variable
      !
      INTEGER, INTENT(IN)     :: nnr
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
      use_smeared_ions = .FALSE.
      !
      IF ( ALLOCATED( vltot_zero ) ) DEALLOCATE(vltot_zero)
      ALLOCATE( vltot_zero( nnr ) )
      vltot_zero = 0.0_DP
      !
      deenviron = 0.0_DP
      !
      ! ... Solvent contribution
      !
      esolvent  = 0.0_DP
      IF ( env_static_permittivity .GT. 1.D0 .OR. env_dielectric_regions .GT. 0 ) THEN
        !
        use_smeared_ions = .TRUE.
        !
        IF ( ALLOCATED( vsolvent ) ) DEALLOCATE(vsolvent)
        ALLOCATE( vsolvent( nnr ) )
        vsolvent = 0.0_DP
        IF ( ALLOCATED( vepsilon ) ) DEALLOCATE(vepsilon)
        ALLOCATE( vepsilon( nnr ) )
        vepsilon = 0.0_DP
        IF ( ALLOCATED( rhopol ) ) DEALLOCATE(rhopol)
        ALLOCATE( rhopol( nnr ) )
        rhopol = 0.0_DP
        !
        CALL solvent_initbase( ifdtype, nfdpoint, nnr )
        !
      END IF
      !
      ! ... Cavity contribution
      !
      ecavity  = 0.0_DP
      IF ( env_surface_tension .GT. 0.D0 ) THEN
        !
        IF ( ALLOCATED( vcavity ) ) DEALLOCATE(vcavity)
        ALLOCATE( vcavity( nnr ) )
        vcavity = 0.0_DP
        !
        CALL cavity_initbase( nnr )
        !
      END IF 
      !
      ! ... Pressure contribution
      !
      epressure  = 0.0_DP
      IF ( env_pressure .NE. 0.D0 ) THEN
        !
        IF ( ALLOCATED( vpressure ) ) DEALLOCATE(vpressure)
        ALLOCATE( vpressure( nnr ) )
        vpressure = 0.0_DP
        !
      ENDIF
      !
      ! ... Periodic correction
      !
      eperiodic = 0.0_DP
      IF ( env_periodicity .NE. 3 ) THEN
        !
        use_smeared_ions = .TRUE.
        !
        IF ( ALLOCATED( vperiodic ) ) DEALLOCATE(vperiodic)
        ALLOCATE( vperiodic( nnr ) ) 
        vperiodic = 0.0_DP
        CALL periodic_initbase( nnr )
        !
      ENDIF
      !
      ! ... Ionic Counter-Charge contribution
      !
      eioncc = 0.0_DP
      IF ( env_ioncc_level .GT. 0 ) THEN 
        !
        use_smeared_ions = .TRUE.
        !
        IF ( ALLOCATED( vioncc ) ) DEALLOCATE(vioncc)
        ALLOCATE( vioncc( nnr ) )
        vioncc = 0.0_DP
        IF ( ALLOCATED( rhoioncc ) ) DEALLOCATE(rhoioncc)
        ALLOCATE( rhoioncc( nnr ) )
        rhoioncc = 0.0_DP
        IF ( env_static_permittivity .GT. 1.D0 .OR. env_dielectric_regions .GT. 0 ) THEN
          IF ( ALLOCATED( rhopolcc ) ) DEALLOCATE(rhopolcc) 
          ALLOCATE( rhopolcc( nnr ) )
          rhopolcc = 0.0_DP
        END IF
        !
        CALL ioncc_initbase( ifdtype, nfdpoint, nnr )
        !
      END IF          
      !
      ! ... External Charges contribution
      !
      eextcharge = 0.0_DP
      IF ( env_external_charges .GT. 0 ) THEN
        !
        use_smeared_ions = .TRUE.
        !
        IF ( ALLOCATED( vextcharge ) ) DEALLOCATE(vextcharge)
        ALLOCATE( vextcharge( nnr ) )
        vextcharge = 0.0_DP
        IF ( ALLOCATED( rhoexternal ) ) DEALLOCATE(rhoexternal)
        ALLOCATE( rhoexternal( nnr ) )
        rhoexternal = 0.0_DP
        !
      END IF 
      !
      ! ... Dielectric regions, contribute to esolvent
      !
      IF ( env_dielectric_regions .GT. 0 ) THEN
        !
        IF ( ALLOCATED(epsstatic) ) DEALLOCATE(epsstatic)
        ALLOCATE( epsstatic( nnr ) )
        epsstatic = 0.0_DP
        IF ( ALLOCATED(epsoptical) ) DEALLOCATE(epsoptical) 
        ALLOCATE( epsoptical( nnr ) )
        epsoptical = 0.0_DP
        !
      END IF 
      !
      RETURN
      !
!--------------------------------------------------------------------
   END SUBROUTINE environ_initbase
!--------------------------------------------------------------------
!--------------------------------------------------------------------
      SUBROUTINE environ_initcell( nnr, n1, n2, n3, ibrav_, omega_, alat_, at_ )
!--------------------------------------------------------------------
!
! Initialize the cell-related quantities to be used in the Environ 
! modules. This initialization is called by electrons.f90, thus it
! is performed at every step of electronic optimization.
!
      ! ... Declares modules
      USE kinds,         ONLY : DP
      USE environ_base,  ONLY : env_periodicity, env_ioncc_level, &
                                env_dielectric_regions, verbose,  &
                                environ_unit
      USE environ_cell,  ONLY : ibrav, omega, alat, domega, ntot, at
      ! ... Local cell-related initializations
      USE periodic,      ONLY : periodic_initcell
      USE ioncc,         ONLY : ioncc_initcell           
      USE epsregion,     ONLY : generate_epsregion        
      !
      IMPLICIT NONE      
      !
      INTEGER, INTENT( IN )    :: nnr
      INTEGER, INTENT( IN )    :: n1, n2, n3
                                  !  Number of dense real space grid points
      INTEGER, INTENT( IN )    :: ibrav_
                                  !  Bravais lattice type
      REAL( DP ), INTENT( IN ) :: omega_
                                  !  Cell volume 
      REAL( DP ), INTENT( IN ) :: alat_
                                  !  Lattice spacing 
      REAL( DP ), INTENT( IN ) :: at_(3,3)
                                  !  Direct lattice primitive vectors 
      !
      ! ... Allocates cell parameters
      !
      ntot = n1*n2*n3
      ibrav = ibrav_
      omega = omega_
      alat = alat_
      domega = omega / DBLE(ntot)
      at = at_
      IF ( verbose .GT. 1 ) THEN
         WRITE(environ_unit,300)n1,n2,n3,ntot
         WRITE(environ_unit,301)ibrav,alat,omega
         WRITE(environ_unit,302)at
      ENDIF
      !
      ! ... Initialize periodicity-correction cell variables
      !
      IF ( env_periodicity .NE. 3 ) CALL periodic_initcell( nnr, at_ )
      !
      ! ... Initialize ionic counter-charge cell variables
      !
      IF ( env_ioncc_level .GT. 0 ) CALL ioncc_initcell( nnr, n1, n2, n3, at_ )
      !
      ! ... Initialize the dielectric regions
      !
      IF ( env_dielectric_regions .GT. 0 ) CALL generate_epsregion( nnr, alat_, omega_, at_ )
      !
      RETURN
      !
300   FORMAT(1X,' n1 = ',i5,' n2 = ',i5,' n3 = ',i5,' ntot = ',i10)
301   FORMAT(1X,' ibrav = ',i2,' alat = ',f14.8,' omega = ',f14.8)
302   FORMAT(1X,' at = ',9f10.4)
!--------------------------------------------------------------------
   END SUBROUTINE environ_initcell
!--------------------------------------------------------------------
!--------------------------------------------------------------------
      SUBROUTINE environ_initions_allocate( nat, nsp )
!--------------------------------------------------------------------
! Allocation and deallocatino of ions-related quantities are done at 
! input processing and need to be done with a separate call, otherwise 
! NEB will not work      
      ! ... Declares modules
      USE environ_ions,      ONLY : ntyp, nat_ => nat, ityp_ => ityp, &
                                    zv_ => zv, tau_ => tau
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN) :: nat, nsp
      !  
      ! ... Allocates ions paramters
      !
      nat_ = nat
      ntyp = nsp
      IF ( ALLOCATED( ityp_) ) DEALLOCATE(ityp_)
      ALLOCATE( ityp_( nat ) )
      IF ( ALLOCATED( zv_) ) DEALLOCATE(zv_) 
      ALLOCATE( zv_( ntyp ) )
      IF ( ALLOCATED( tau_ ) ) DEALLOCATE(tau_)
      ALLOCATE( tau_( 3, nat ) )
      !
      RETURN
      !
!--------------------------------------------------------------------
   END SUBROUTINE environ_initions_allocate
!--------------------------------------------------------------------
!--------------------------------------------------------------------
      SUBROUTINE environ_initions( nnr, nat, nsp, ityp, zv, tau, alat )
!--------------------------------------------------------------------
!
! Initialize the ions-related quantities to be used in the Environ 
! modules. This initialization is called by electrons.f90, thus it
! is performed at every step of electronic optimization. It may not
! be the most efficient choice, but it is a safe choice.
!
      ! ... Declares modules
      USE kinds,             ONLY : DP
      USE environ_base,      ONLY : env_periodicity,                &
                                    use_smeared_ions, atomicspread, &
                                    env_external_charges, verbose,  &
                                    environ_unit
      USE environ_ions,      ONLY : ntyp, nat_ => nat, ityp_ => ityp, &
                                    zv_ => zv, tau_ => tau, rhoions, &
                                    zvtot, avg_pos
      USE periodic,          ONLY : periodic_initions
      USE extcharge,         ONLY : generate_extcharge
      USE generate_function, ONLY : generate_gaussian
      !
      IMPLICIT NONE      
      !
      INTEGER, INTENT( IN )     :: nnr, nat, nsp
      INTEGER, INTENT( IN )     :: ityp( nat )
      REAL ( DP ), INTENT( IN ) :: zv( nsp )
      REAL ( DP ), INTENT( IN ) :: tau( 3, nat )
      REAL ( DP ), INTENT( IN ) :: alat
      !
      INTEGER :: ia, dim, axis
      REAL ( DP ) :: charge, spread, pos(3)
      !
      ! ... Allocates ions paramters
      !
      ityp_ = ityp
      zv_   = zv
      tau_  = tau
      IF ( verbose .GT. 1 ) WRITE(environ_unit,200)nat,nsp
      !
      ! ... If needed, generate a fictitious ion density using gaussians
      !
      IF ( use_smeared_ions ) THEN
        !
        IF ( .NOT. ALLOCATED(rhoions) ) ALLOCATE( rhoions( nnr ) )
        rhoions = 0.D0
        DO ia = 1, nat
          IF ( verbose .GT. 2 ) WRITE(environ_unit,201)ia,ityp(ia),zv(ityp(ia)),tau(:,ia)*alat
          pos( : ) = tau( :, ia)
          charge = - zv( ityp( ia ) )
          spread = atomicspread( ityp( ia ) )
          dim = 0 
          axis = 1
          CALL generate_gaussian( nnr, dim, axis, charge, spread, pos, rhoions ) 
        END DO
        !
      END IF
      !
      ! ... Center and amount of ionic charge used by three sub-modules, compute it here
      !
      zvtot = 0.D0
      avg_pos(:) = 0.D0
      !
      DO ia = 1, nat
        !
        zvtot = zvtot + zv(ityp(ia))
        !
        avg_pos(:) = avg_pos(:) + tau(:,ia)*zv(ityp(ia))
        !
      END DO
      !
      avg_pos(:) = avg_pos(:) / zvtot
      IF ( verbose .GT. 1 ) WRITE(environ_unit,202)avg_pos*alat
      !
      ! ... Initialize periodicity-correction ionic variables
      !
      IF ( env_periodicity .NE. 3 ) CALL periodic_initions( nnr, nat, ntyp, ityp, zv, tau, alat, rhoions )
      !
      ! ... Initialize the external charges
      !
      IF ( env_external_charges .GT. 0 ) CALL generate_extcharge( nnr, alat )
      !
      RETURN
      !
200   FORMAT(1X,'nat = ',i8,'  nsp = ',i8)
201   FORMAT(1X,'iat = ',i8,' ityp = ',i4,' zv = ',f8.2,' pos = ',3f14.8)
202   FORMAT(1X,'avg_pos = ',3f14.8)
!--------------------------------------------------------------------
   END SUBROUTINE environ_initions
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
      USE environ_base, ONLY : environ_unit, vltot_zero,             & 
                               env_static_permittivity, vsolvent,    &
                               vepsilon, rhopol,                     &
                               env_surface_tension, vcavity,         &
                               vpressure,                            &
                               env_periodicity, vperiodic,           &
                               env_ioncc_level, vioncc,              &
                               solvationrad, atomicspread,           &
                               env_external_charges, vextcharge,     &
                               rhoexternal, epsstatic, epsoptical,   &
                               env_dielectric_regions
      USE environ_ions, ONLY : ityp, zv, tau, rhoions
      !
      ! Local clean up subroutines for the different contributions
      !
      USE solvent,        ONLY : solvent_clean
      USE periodic,       ONLY : periodic_clean
      USE cavity,         ONLY : cavity_clean
      USE ioncc,          ONLY : ioncc_clean
      USE control_flags,  ONLY : tddfpt
      !
      IMPLICIT NONE      
      !
      LOGICAL, INTENT(IN) :: lflag
      LOGICAL :: opnd
      !
      ! ... Deallocates environment variables allocated in input
      !
      IF( lflag ) THEN
        !
        IF (ALLOCATED(solvationrad)) DEALLOCATE( solvationrad )
        IF (ALLOCATED(atomicspread)) DEALLOCATE( atomicspread )
        !
        ! ... environ_ions variables
        !
        IF ( ALLOCATED( ityp ) )   DEALLOCATE( ityp )
        IF ( ALLOCATED( zv ) )     DEALLOCATE( zv )
        IF ( ALLOCATED( tau ) )    DEALLOCATE( tau )
        !
      END IF 
      !
      ! ... Deallocates environment variables
      !
      INQUIRE( unit=environ_unit, opened= opnd ) 
      IF ( opnd ) CLOSE(unit=environ_unit)
      !
      IF ( ALLOCATED( vltot_zero ) )  DEALLOCATE( vltot_zero )
      !
      ! ... environ_base variables
      !
      IF ( ALLOCATED( vsolvent ) )    DEALLOCATE( vsolvent )
      IF ( ALLOCATED( vepsilon ) )    DEALLOCATE( vepsilon )
      IF ( ALLOCATED( rhopol ) )      DEALLOCATE( rhopol )
      IF ( ALLOCATED( vcavity ) )     DEALLOCATE( vcavity )
      IF ( ALLOCATED( vpressure ) )   DEALLOCATE( vpressure )
      IF ( ALLOCATED( vperiodic ) )   DEALLOCATE( vperiodic )
      IF ( ALLOCATED( vioncc ) )      DEALLOCATE( vioncc )
      IF ( ALLOCATED( vextcharge ) )  DEALLOCATE( vextcharge )
      IF ( ALLOCATED( rhoexternal ) ) DEALLOCATE( rhoexternal )
      IF ( ALLOCATED( epsstatic ) )   DEALLOCATE( epsstatic )
      IF ( ALLOCATED( epsoptical ) )  DEALLOCATE( epsoptical )
      !
      ! ... environ_ions variables
      !
      IF ( .NOT. tddfpt) THEN
        IF ( ALLOCATED( rhoions ) ) DEALLOCATE( rhoions )
      END IF
      !
      ! ... internal clean up of environ modules
      !
      IF ( env_static_permittivity .GT. 1.D0 & 
           .OR. env_dielectric_regions .GT. 0 ) CALL solvent_clean()
      IF ( env_surface_tension .GT. 0.D0 ) CALL cavity_clean()
      IF ( env_periodicity .NE. 3 ) CALL periodic_clean()
      IF ( env_ioncc_level .GT. 0 ) CALL ioncc_clean()
      !
      RETURN
      !
!--------------------------------------------------------------------
   END SUBROUTINE environ_clean
!--------------------------------------------------------------------
END MODULE environ_init
