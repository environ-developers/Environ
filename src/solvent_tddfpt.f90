!
! Copyright (C) 2001-2013 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE solvent_tddfpt
!----------------------------------------------------------------------------
      !
      ! ... Module for calculation of the response "polarization" and "dielectric" potentials.
      ! ... Inspired by Environ/src/solvent.f90
      ! ... Written by I. Timrov 09/2013
      !
      USE kinds,              ONLY : DP
      USE constants,          ONLY : fpi
      USE fd_gradient,        ONLY : init_fd_gradient, calc_fd_gradient
      USE control_flags,      ONLY : tddfpt
      USE io_global,          ONLY : stdout
      USE environ_base,       ONLY : e2, gvtot0, eps_mode
      USE generate_f_of_rho,  ONLY : generate_dielectric
      USE environ_debug,      ONLY : write_cube
      !
      IMPLICIT NONE
      LOGICAL :: initialized
      INTEGER, PARAMETER :: maxiter = 100
      REAL(DP) :: charge, ehart, polthr
      INTEGER  :: nfdpoint, ncfd
      INTEGER, ALLOCATABLE  :: icfd(:)
      REAL(DP), ALLOCATABLE :: eps(:),          & ! dielectric function
                               deps(:),         & ! derivative wrt density of the dielectric function
                               d2eps(:),        & ! second der. wrt density of the dielectric function
                               logeps(:),       & ! logarithm of the dielectric function
                               gradlogeps(:,:), & ! gradient of logarithm of the dielectric function
                               drhopol(:),      & ! response polarization charge-density
                               drhoiter(:)        ! iterative response charge-density

       SAVE

       PRIVATE

       PUBLIC :: solvent_initbase_tddfpt, solvent_clean_tddfpt, calc_vsolvent_tddfpt

CONTAINS
!--------------------------------------------------------------------
 SUBROUTINE solvent_initbase_tddfpt( ifdtype, nfdpoint_, nnr )
!--------------------------------------------------------------------
      !
      ! ... Local initializations
      !
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: ifdtype, nfdpoint_, nnr
      !
      initialized = .FALSE.
      !
      IF (ALLOCATED(drhoiter)) DEALLOCATE(drhoiter)
      ALLOCATE(drhoiter(nnr))
      drhoiter = 0.D0
      !
      IF (ALLOCATED(drhopol)) DEALLOCATE(drhopol)
      ALLOCATE(drhopol(nnr))
      drhopol = 0.D0
      !
      IF (ALLOCATED(eps)) DEALLOCATE(eps)
      ALLOCATE(eps(nnr))
      eps = 0.0_DP
      !
      IF (ALLOCATED(deps)) DEALLOCATE(deps)
      ALLOCATE(deps(nnr))
      deps = 0.0_DP
      !
      IF (ALLOCATED(d2eps)) DEALLOCATE(d2eps)
      ALLOCATE(d2eps(nnr))
      d2eps = 0.0_DP
      !
      IF (ALLOCATED(gradlogeps)) DEALLOCATE(gradlogeps)
      ALLOCATE(gradlogeps(nnr,3))
      gradlogeps = 0.0_DP
      !
      nfdpoint = nfdpoint_
      ALLOCATE( icfd(-nfdpoint:nfdpoint) )
      CALL init_fd_gradient( ifdtype, nfdpoint, ncfd, icfd )
      !
      RETURN
      !
!--------------------------------------------------------------------
 END SUBROUTINE solvent_initbase_tddfpt
!--------------------------------------------------------------------
!--------------------------------------------------------------------
 SUBROUTINE solvent_clean_tddfpt()
!--------------------------------------------------------------------
      !
      ! ... Local clean up
      !
      IMPLICIT NONE
      !
      IF(ALLOCATED(eps))        DEALLOCATE(eps)
      IF(ALLOCATED(deps))       DEALLOCATE(deps)
      IF(ALLOCATED(d2eps))       DEALLOCATE(d2eps)
      IF(ALLOCATED(gradlogeps)) DEALLOCATE(gradlogeps)
      IF(ALLOCATED(icfd))       DEALLOCATE(icfd)
      IF(ALLOCATED(drhoiter))   DEALLOCATE(drhoiter)
      IF(ALLOCATED(drhopol))    DEALLOCATE(drhopol)
      IF(ALLOCATED(gvtot0))     DEALLOCATE(gvtot0)
      !
      RETURN
      !
!---------------------------------------------------------------------------------
 END SUBROUTINE solvent_clean_tddfpt
!---------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
 SUBROUTINE calc_vsolvent_tddfpt(nnr, nspin, rho_0, drho_elec, dv_pol, dv_epsilon)
!-----------------------------------------------------------------------------------
      !
      ! ... This subroutine calculates:
      ! ... 1. the response "polarization" potential from the response polarization density
      ! ... 2. the response "dielectric" potential
      !
      USE calc_rhopol,              ONLY : iterative_rhopol
      USE environ_base,             ONLY : tolrhopol, mixrhopol
      USE environ_ions,             ONLY : rhoions
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN)     :: nnr,           & ! number of grid points in R-space
                                 nspin            ! if nspin=2 spin-polarized case (not supported)
      REAL( DP ), INTENT(IN)  :: rho_0(nnr),    & ! ground-state charge-density
                                 drho_elec(nnr)   ! response charge-density
      REAL( DP ), INTENT(OUT) :: dv_pol(nnr),   & ! response polarization potential
                                 dv_epsilon(nnr)  ! response dielectric potential
      !
      REAL( DP ), ALLOCATABLE :: gdvtot(:,:),   & ! gradient of the total response potential
                                 drhotot(:),    & ! total response charge-density
                                 drhoaux(:,:),  & ! temporary array for the density
                                 dvaux(:,:),    & ! temporary array for the potential
                                 rhodiel(:)       ! density which is used to define the cavity
      INTEGER                 :: ir
      !
      IF (.NOT. tddfpt) RETURN
      !
      CALL start_clock( 'calc_vsolvent_tddfpt' )
      !
      ! ... Generate the dielectric function and its derivative
      ! ... from the unperturbed charge-density.
      ! ... They must be generated only once contrarily to the
      ! ... ground-state SCF calculation.
      !
      IF ( .NOT. initialized ) THEN
         !
         ALLOCATE( rhodiel(nnr) )
         rhodiel = 0.D0
         !
         SELECT CASE( TRIM( eps_mode ) )
           !
           CASE( 'electronic' )
             !
             rhodiel = rho_0
             !
           CASE( 'ionic' )
             !
             CALL errore('calc_vsolvent_tddfpt', &
                   & 'eps_mode = ionic  is not supported', 1 )
             !
           CASE( 'full' )
             !
             rhodiel = rho_0 - rhoions
             !
           CASE( 'external' )
             !
             CALL errore('calc_vsolvent_tddfpt', &
                   & 'eps_mode = external is not supported', 1 )
             !
         END SELECT
         !
         eps  = 0.D0
         deps = 0.D0
         d2eps = 0.D0
         CALL generate_dielectric( nnr, rhodiel, eps, deps, d2eps, .TRUE. )
         !
         DEALLOCATE( rhoions )
         DEALLOCATE( rhodiel )
         !
         initialized = .TRUE.
         !
      ENDIF
      !
      ! ... Calculates a gradient of logarithm of dielectric permittivity.
      !
      ALLOCATE (logeps(nnr))
      logeps(:) = LOG(eps(:))
      CALL calc_fd_gradient( nfdpoint, icfd, ncfd, nnr, logeps, gradlogeps )
      DEALLOCATE( logeps )
      !
      ! ... Iterative calculation of the response polarization density.
      !
      CALL iterative_rhopol(nnr, nspin, maxiter, .FALSE., tolrhopol, &
                                 & mixrhopol, drho_elec, eps, gradlogeps, drhoiter)
      !
      ! ... Sum up the two turms and obtain the final expression for the
      ! ... response polarization density
      !
      drhopol(:) = drhoiter(:) + (1.D0-eps(:))/eps(:) * drho_elec(:)
      !
      ! ... Calculate the response "polarization" potential
      ! ... from the response polarization density
      ! ... Note: ehart and charge are computed, but they are not needed
      !
      ALLOCATE( drhoaux( nnr, nspin ) )
      ALLOCATE( dvaux( nnr, nspin ) )
      !
      drhoaux(:,1) = drhopol(:)
      IF ( nspin .EQ. 2 ) drhoaux(:,2) = 0.D0
      dvaux = 0.D0
      !
      CALL v_h_of_rho_r( drhoaux, ehart, charge, dvaux )
      !
      dv_pol(:) = dvaux(:,1)
      !
      DEALLOCATE( drhoaux )
      DEALLOCATE( dvaux )
      !
      !!!!!!!!!!!!!!!! Response dielectric potential !!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      ALLOCATE( drhotot(nnr) )
      ALLOCATE( gdvtot(3,nnr) )
      !
      ! ... Total response charge-density
      !
      drhotot = drho_elec + drhopol
      !
      ! ... Calculate a gradient of the total response potential [like in Eq.(50)]
      !
      CALL gradv_h_of_rho_r ( drhotot, gdvtot )
      !
      DO ir = 1, nnr
         dv_epsilon(ir) = - SUM( gvtot0(:,ir) * gdvtot(:,ir) ) * deps(ir) / (fpi * e2)
      END DO
      !
      DEALLOCATE( drhotot )
      DEALLOCATE( gdvtot )
      !
      CALL stop_clock( 'calc_vsolvent_tddfpt' )
      !
      RETURN
      !
!--------------------------------------------------------------------
END SUBROUTINE calc_vsolvent_tddfpt
!--------------------------------------------------------------------
!--------------------------------------------------------------------
END MODULE solvent_tddfpt
!--------------------------------------------------------------------
