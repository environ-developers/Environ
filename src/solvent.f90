!
! Copyright (C) 2007-2008 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Module for implicit solvation according to the Self-Consistent
! Continuum Solvation (SCCS) model of Andreussi et al. J. Chem. Phys. 136,
! 064102 (2012).
!
! original version by O. Andreussi, I. Dabo and N. Marzari
!
!----------------------------------------------------------------------------
MODULE solvent
!----------------------------------------------------------------------------
  !
  ! ... The variables needed to include the effects of a continuum solvent
  !     on the system.
  !
  USE kinds,             ONLY : DP
  USE io_global,         ONLY : stdout
  USE mp,                ONLY : mp_sum
  USE mp_bands,          ONLY : intra_bgrp_comm
  USE constants,         ONLY : pi, tpi, fpi
  USE environ_base,      ONLY : verbose, environ_unit, e2, env_static_permittivity, &
                                eps_mode, rhopol, alpha, solvationrad, corespread, env_periodicity, mixtype
  USE environ_cell,      ONLY : domega, ntot, omega, alat, at
  USE periodic,          ONLY : calc_gradvperiodic
  USE fd_gradient,       ONLY : init_fd_gradient, calc_fd_gradient
  USE generate_function, ONLY : generate_exponential, generate_gradexponential, &
                                generate_gaussian, generate_gradgaussian
! BACKWARD COMPATIBILITY
! Compatible with QE-5.1.X, QE-5.2.0
!
! Compatible with QE-5.2.1, QE-5.2.2, QE-5.3.0, QE-svn
  USE generate_function, ONLY : generate_erfc, generate_graderfc
  USE control_flags,     ONLY : tddfpt
! END BACKWARD COMPATIBILITY
  USE environ_debug,     ONLY : write_cube
  USE generate_f_of_rho, ONLY : generate_dielectric
  !
  IMPLICIT NONE
  !
  INTEGER, PARAMETER :: maxiter = 100
  INTEGER :: nfdpoint
  INTEGER :: ncfd
  INTEGER, ALLOCATABLE :: icfd(:)
  REAL(DP), ALLOCATABLE :: eps(:), deps(:), d2eps(:), gradlogeps(:,:)
  REAL(DP), ALLOCATABLE :: invsqrteps(:), factsqrteps(:)
  REAL(DP), ALLOCATABLE :: rhoiter(:)
  !
  SAVE

  PRIVATE

  PUBLIC :: eps, gradlogeps
  PUBLIC :: solvent_initbase, calc_vsolvent, calc_esolvent, calc_fsolvent, calc_esolvent_of_V, solvent_clean

CONTAINS
!--------------------------------------------------------------------
  SUBROUTINE solvent_initbase( ifdtype, nfdpoint_, nnr )
!--------------------------------------------------------------------

    ! ... Local initializations: iterative polarization density
    !     and real-space finite-differences gradient coefficients

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: ifdtype, nfdpoint_, nnr

    IF(ALLOCATED(rhoiter)) DEALLOCATE(rhoiter)
    ALLOCATE(rhoiter(nnr))
    rhoiter = 0.D0

    IF(ALLOCATED(eps)) DEALLOCATE(eps)
    ALLOCATE( eps( nnr ) )
    eps = 0.0_DP
    IF(ALLOCATED(deps)) DEALLOCATE(deps)
    ALLOCATE( deps( nnr ) )
    deps = 0.0_DP
    IF(ALLOCATED(d2eps)) DEALLOCATE(d2eps)
    ALLOCATE( d2eps( nnr ) )
    d2eps = 0.0_DP
    IF(ALLOCATED(gradlogeps)) DEALLOCATE(gradlogeps)
    ALLOCATE( gradlogeps( 3, nnr ) )
    gradlogeps = 0.0_DP
    IF ( mixtype .EQ. 'pcg' .OR. mixtype .EQ. 'psd' ) THEN
      IF(ALLOCATED(invsqrteps)) DEALLOCATE(invsqrteps)
      ALLOCATE( invsqrteps( nnr ) )
      invsqrteps = 0.0_DP
      IF(ALLOCATED(factsqrteps)) DEALLOCATE(factsqrteps)
      ALLOCATE( factsqrteps( nnr ) )
      factsqrteps = 0.0_DP
    ENDIF

    nfdpoint = nfdpoint_
    IF(ALLOCATED(icfd)) DEALLOCATE(icfd)
    ALLOCATE( icfd(-nfdpoint:nfdpoint) )
    CALL init_fd_gradient( ifdtype, nfdpoint, ncfd, icfd )

    RETURN

!--------------------------------------------------------------------
  END SUBROUTINE solvent_initbase
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE solvent_clean()
!--------------------------------------------------------------------

    ! ... Local clean up

    IMPLICIT NONE

    IF(ALLOCATED(eps)) DEALLOCATE(eps)
    IF(ALLOCATED(deps)) DEALLOCATE(deps)
    IF(ALLOCATED(d2eps)) DEALLOCATE(d2eps)
    IF(ALLOCATED(gradlogeps)) DEALLOCATE(gradlogeps)
    IF(ALLOCATED(rhoiter)) DEALLOCATE(rhoiter)
    IF(ALLOCATED(icfd)) DEALLOCATE(icfd)
    IF ( mixtype .EQ. 'pcg' .OR. mixtype .EQ. 'psd' ) THEN
      IF(ALLOCATED(invsqrteps)) DEALLOCATE(invsqrteps)
      IF(ALLOCATED(factsqrteps)) DEALLOCATE(factsqrteps)
    ENDIF

    RETURN

!--------------------------------------------------------------------
  END SUBROUTINE solvent_clean
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE calc_vsolvent( nnr, nspin, dr2, rhoelec, vsolvent, vepsilon )
!--------------------------------------------------------------------
    !
    ! ... Calculates the solvent contribution to the Hartree potential
    !
    USE constants,         ONLY : fpi, sqrtpi
    USE environ_ions,      ONLY : rhoions, nat, tau, zv, ityp, avg_pos
    USE environ_cell,      ONLY : at, alat, omega
    USE environ_base,      ONLY : e2, tolrhopol, mixrhopol, add_jellium,   &
                                  env_external_charges, rhoexternal,       &
                                  gvtot0
    USE calc_rhopol,       ONLY : iterative_rhopol, rhopol_of_v
    USE calc_vpol,         ONLY : iterative_vpol
    !
    IMPLICIT NONE
    !
    ! ... Declares variables
    !
    INTEGER, INTENT(IN)       :: nnr, nspin
    REAL( DP ), INTENT(IN)    :: dr2
    REAL( DP ), INTENT(IN)    :: rhoelec( nnr )
    REAL( DP ), INTENT(OUT)   :: vsolvent( nnr )
    REAL( DP ), INTENT(OUT)   :: vepsilon( nnr )
    !
    ! ... Local variables
    !
    INTEGER                   :: iter, ir, icor, na, iext, dim, axis
    !
    REAL( DP )                :: spread, charge, ehart, gvtot2
    REAL( DP )                :: pos( 3 ), pos0( 3 )
    REAL( DP )                :: polthr, rhoavg
    !
    REAL( DP ), ALLOCATABLE   :: rhodiel( : )
    REAL( DP ), ALLOCATABLE   :: logeps( : )
    !
    REAL( DP ), ALLOCATABLE   :: gradrho( :, : )
    REAL( DP ), ALLOCATABLE   :: hessrho( :, :, : )
    !
    REAL( DP ), ALLOCATABLE   :: rhosolute( : )
    REAL( DP ), ALLOCATABLE   :: rhotot( : )
    REAL( DP ), ALLOCATABLE   :: gvtot( :, : )
    !
    REAL( DP ), ALLOCATABLE   :: rhoaux( : , : )
    REAL( DP ), ALLOCATABLE   :: vaux( : , : )

    CALL start_clock( 'calc_vsolv' )

!    polthr = MIN(dr2*tolrhopol,1.D-7)
!    polthr = MAX(polthr,1.D-14)
    polthr = tolrhopol

    ! ... Calculates dielectric permittivity

    CALL start_clock( 'dielectric' )

    ALLOCATE( rhodiel( nnr ) )
    rhodiel = 0.D0

    SELECT CASE( TRIM( eps_mode ) )

    CASE( 'electronic' )

      rhodiel = rhoelec

    CASE( 'ionic' )

      DO na = 1, nat
        spread = solvationrad( ityp( na ) ) * alpha
        pos( : ) = tau( :, na )
        CALL generate_exponential( nnr, spread, pos, rhodiel )
      ENDDO

    CASE( 'full' )

      rhodiel = rhoelec
      DO na = 1, nat
        IF ( corespread( ityp( na ) ) .LE. 0.D0 ) CYCLE
        dim = 0
        axis = 1
        charge = zv( ityp( na ) )
        spread = corespread( ityp( na ) )
        pos( : ) = tau( :, na )
        CALL generate_gaussian( nnr, dim, axis, charge, spread, pos, rhodiel )
      ENDDO

    END SELECT

    eps = 0.D0
    deps = 0.D0
    d2eps = 0.D0
    CALL generate_dielectric( nnr, rhodiel, eps, deps, d2eps, .FALSE. )
    IF ( verbose .GE. 2 ) CALL write_cube( nnr, eps, 'eps.cube' )
    IF ( verbose .GE. 3 ) CALL write_cube( nnr, d2eps, 'd2eps.cube' )
    IF ( verbose .GE. 3 ) CALL write_cube( nnr, rhodiel, 'rhodiel.cube' )

    ! ... Calculates gradient of logarithm of dielectric permittivity

    ALLOCATE ( logeps( nnr ) )
    logeps(:) = LOG(eps(:))
    CALL calc_fd_gradient( nfdpoint, icfd, ncfd, nnr, logeps, gradlogeps )
    DEALLOCATE( logeps )

    CALL stop_clock( 'dielectric' )

    ! ... Calculates total molecular density

    ALLOCATE( rhosolute( nnr ) )
    rhosolute = rhoelec + rhoions
    IF ( env_external_charges .GT. 0 ) rhosolute = rhosolute + rhoexternal
    rhoavg = SUM(rhosolute)/DBLE(ntot)
    CALL mp_sum( rhoavg, intra_bgrp_comm )
    IF ( add_jellium ) rhosolute(:) = rhosolute(:) - rhoavg
    IF ( verbose .GE. 3 ) CALL write_cube( nnr, rhosolute, 'rhosolute.cube' )

    ! ... Either loop over the potential using a preconditioner
    !     pcg = preconditioned conjugate gradient
    !     psd = preconditioned steepest descent

    IF ( mixtype .EQ. 'pcg' .OR. mixtype .EQ. 'psd' ) THEN

      ! ... Calculates square root of dielectric permittivity and its laplacian

      invsqrteps(:) = 1.D0 / SQRT(eps(:))
      ALLOCATE( hessrho( 3, 3, nnr ) )
      ALLOCATE( gradrho( 3, nnr ) )
      hessrho = 0.D0
      gradrho = 0.D0
      CALL external_hessian( rhodiel, gradrho, hessrho )
      DO ir = 1, nnr
         factsqrteps(ir) = deps(ir)*(hessrho(1,1,ir)+hessrho(2,2,ir)+hessrho(3,3,ir)) + &
             & SUM(gradrho(:,ir)*gradrho(:,ir))*(d2eps(ir)-0.5D0*deps(ir)**2/eps(ir))
      ENDDO
      factsqrteps = factsqrteps * 0.5D0 / e2 / fpi
      DEALLOCATE( hessrho )
      DEALLOCATE( gradrho )
      IF ( verbose .GE. 3 ) CALL write_cube( nnr, factsqrteps, 'factsqrteps.cube' )

      ! ... Iterative calculation of polarization potential

      CALL iterative_vpol(nnr,nspin,maxiter,.FALSE.,polthr,rhosolute,invsqrteps,&
           & factsqrteps,rhopol,vsolvent)

      ! ... Calculates polarization charge corresponding to the computed potential

      CALL rhopol_of_v(nnr,nspin,eps,gradlogeps,rhosolute,vsolvent,rhopol)

    ! ... or loop over the polarization charge (original SCCS approach)

    ELSE

      ! ... Iterative calculation of polarization charges

      CALL iterative_rhopol(nnr,nspin,maxiter,.FALSE.,polthr,mixrhopol,rhosolute,eps,&
                            gradlogeps,rhoiter)

      rhopol(:) = rhoiter(:) + (1.D0-eps(:))/eps(:)*rhosolute(:)

      ! ... If jellium is present, remove homogeneous long-range polarization (constant shift)

      IF ( add_jellium ) rhopol(:) = rhopol(:) + (1.D0-env_static_permittivity)/env_static_permittivity * rhoavg
      IF ( verbose .GE. 2 ) CALL write_cube( nnr, rhopol, 'rhopol.cube' )
      IF ( verbose .GE. 3 ) CALL write_cube( nnr, rhoiter, 'rhoiter.cube' )

      ! ... Calculates corrective potential

      ALLOCATE( rhoaux( nnr, nspin ) )
      ALLOCATE( vaux( nnr, nspin ) )
      rhoaux( :, 1 ) = rhopol(:)
      IF ( nspin .EQ. 2 ) rhoaux( :, 2 ) = 0.D0
      vaux = 0.D0
      ! ... Calculate the "polarization" potential from the polarization density
      CALL v_h_of_rho_r( rhoaux, ehart, charge, vaux )
      vsolvent(:) = vaux( :, 1 )
      DEALLOCATE( rhoaux )
      DEALLOCATE( vaux )

    ENDIF

    DEALLOCATE( rhodiel )

    ! ... If verbose, print the total electrostatic potential
    ! ... from the total density (solute+pol)

    IF ( verbose .GE. 2 .AND. env_periodicity .EQ. 3 ) THEN
      ALLOCATE( rhoaux( nnr, nspin ) )
      ALLOCATE( vaux( nnr, nspin ) )
      rhoaux( :, 1 ) = rhosolute(:) + rhopol(:)
      IF ( nspin .EQ. 2 ) rhoaux( :, 2 ) = 0.D0
      vaux = 0.D0
      CALL v_h_of_rho_r( rhoaux, ehart, charge, vaux )
      CALL write_cube( nnr, vaux( :, 1 ), 'vfull.cube' )
      DEALLOCATE( rhoaux )
      DEALLOCATE( vaux )
    END IF

    CALL stop_clock( 'calc_vsolv' )

    ! ... Calculates extra term of functional derivative

    CALL start_clock( 'calc_veps' )

    ALLOCATE( rhotot ( nnr ) )
    ALLOCATE( gvtot ( 3, nnr ) )
    rhotot = rhosolute + rhopol
    !
    ! ... Calculate a gradient of the total potential [like in Eq.(50)]
    ! ... which is needed for a calculation of V_eps, see Eq.(26).
    !
    CALL gradv_h_of_rho_r ( rhotot, gvtot )
    !
    IF ( env_periodicity .NE. 3 ) CALL calc_gradvperiodic( nnr, rhotot, gvtot )

    DEALLOCATE( rhosolute )
    DEALLOCATE( rhotot )

    DO ir = 1, nnr
      gvtot2 = SUM( gvtot(:,ir) *  gvtot(:,ir) )
      vepsilon(ir) = - gvtot2 * deps(ir) / 2.D0 / fpi / e2
    END DO

! BACKWARD COMPATIBILITY
! Compatible with QE-5.1.X, QE-5.2.0
!
! Compatible with QE-5.2.1 and following releases
    IF (tddfpt) THEN
       !
       ! Save the gradient of the total electrostatic potential
       ! (elec + pol + ions) for the TDDFPT calculation.
       !
       ALLOCATE( gvtot0(3,nnr) )
       gvtot0 = gvtot
       !
    ENDIF
! END BACKWARD COMPATIBILITY

    DEALLOCATE( gvtot )

    CALL stop_clock( 'calc_veps' )

    RETURN

!--------------------------------------------------------------------
  END SUBROUTINE calc_vsolvent
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE calc_esolvent( nnr, nspin, rho, esolvent )
!--------------------------------------------------------------------
    USE environ_base,  ONLY : rhopol, env_external_charges, rhoexternal
    USE environ_ions,  ONLY : rhoions
!DEBUG
!    USE environ_ions,  ONLY : zv, nat, tau, ityp
!    USE environ_base,  ONLY : atomicspread
!DEBUG
    !
    ! ... Declares variables
    !
    INTEGER, INTENT(IN)     :: nnr, nspin
    REAL( DP ), INTENT(IN)  :: rho( nnr )
    REAL( DP ), INTENT(OUT) :: esolvent
    !
    REAL( DP ) :: ehart, charge
    REAL( DP ), ALLOCATABLE :: vaux(:,:), rhoaux(:,:)
!DEBUG
!    INTEGER :: ia
!    REAL( DP ) :: spread, pos(3), ecorrection
!DEBUG
    !
    CALL start_clock ('calc_esolv')
    !
    !  Solvent energy with gaussian-spread nuclei
    !
    ALLOCATE(vaux(nnr,nspin))
    vaux = 0.D0
    ALLOCATE(rhoaux(nnr,nspin))
    rhoaux( :, 1 ) = rhoions(:)+rho(:)
    IF ( env_external_charges .GT. 0 ) &
               rhoaux( :, 1 ) = rhoaux( :, 1 ) + rhoexternal(:)
    IF ( nspin .EQ. 2 ) rhoaux( :, 2 ) = 0.D0
    CALL v_h_of_rho_r( rhoaux, ehart, charge, vaux )
    esolvent = esolvent +                                        &
               0.5D0 * SUM( vaux( :, 1 ) * rhopol( : ) ) * domega
    DEALLOCATE( rhoaux )
    DEALLOCATE( vaux )
    !
    CALL mp_sum(  esolvent, intra_bgrp_comm )
!DEBUG
!
!    This part computes the effect of treating the nuclei as point-charges
!    instead of gaussian densities
!
!    ALLOCATE(vaux(nnr,1))
!    vaux = 0.D0
!    DO ia = 1, nat
!      charge = -zv(ityp(ia))
!      spread = atomicspread(ityp(ia))
!      pos(:) = tau(:,ia)
!      CALL calc_vsr( nnr, charge, spread, pos, vaux )
!    ENDDO
!    ecorrection = 0.5D0 * SUM( vaux( :, 1 ) * rhopol( : ) ) * domega
!    CALL mp_sum( ecorrection, intra_bgrp_comm )
!    esolvent = esolvent + ecorrection
!    DEALLOCATE( vaux )
!
!DEBUG
    !
    CALL stop_clock ('calc_esolv')
    !
    RETURN
!--------------------------------------------------------------------
  END SUBROUTINE calc_esolvent
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE calc_fsolvent( nnr, nat, f )
!--------------------------------------------------------------------
    !
    ! ... Calculates the solvent contribution to the forces when
    !     the dielectric is defined on the ionic positions
    !
    USE kinds,             ONLY : DP
    USE environ_ions,      ONLY : ntyp, ityp, zv, tau
    USE environ_base,      ONLY : atomicspread, solvationrad,         &
                                  vepsilon, vsolvent
    !
    IMPLICIT NONE
    !
    ! ... Declares variables
    !
    INTEGER, INTENT(IN)       :: nnr, nat
    REAL(DP), INTENT(INOUT)   :: f( 3, nat )
    !
    ! ... Local variables
    !
    INTEGER                 :: na, ipol, ir, dim, axis
    REAL(DP)                :: spread, charge
    REAL(DP)                :: pos(3)
    REAL(DP)                :: ftmp( 3, nat )
    REAL(DP), ALLOCATABLE   :: grhoions( : , :  )
    !
    ! ... Polarization contribution to the forces
    !
    ALLOCATE( grhoions(3,nnr) )
    !
    ftmp = 0.D0
    DO na = 1, nat
      !
      grhoions = 0.D0
      pos( : ) = tau( :, na )
      spread = atomicspread( ityp ( na ) )
      charge = zv( ityp ( na ) )
      dim = 0
      axis = 1
      CALL generate_gradgaussian( nnr, dim, axis, charge, spread, pos, grhoions )
      !
      DO ipol = 1,3
        ftmp( ipol, na ) = &
          - SUM( grhoions(ipol, : ) * vsolvent ( : ) ) &
          * domega
      END DO
      !
    END DO
    !
    ! ... Dielectric cavity contribution to the forces
    ! f(I) = int( gvtot2 * grad(rhoion(I)) * depsilon/drho ) / (8*pi*e2)
    !      = int( vepsilon * grad(rhoion(I)) )
    !
    DO na = 1, nat
      !
      grhoions = 0.D0
      pos( : ) = tau( :, na )
      spread = solvationrad( ityp( na ) )
      CALL generate_gradexponential( nnr, spread, pos, grhoions )
      !
      DO ipol = 1, 3
        ftmp( ipol, na ) = ftmp( ipol, na ) + &
          SUM( grhoions( ipol, : ) * vepsilon( : ) ) &
          * domega
      END DO
      !
    END DO
    !
    DEALLOCATE( grhoions )
    !
    CALL mp_sum( ftmp, intra_bgrp_comm )
    !
    f = f + ftmp
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE calc_fsolvent
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE calc_esolvent_of_V( nnr, nspin, rhoelec, ejellium,  &
                        & eperiodic, rhopol_of_V)
!--------------------------------------------------------------------
    !
    USE environ_ions,      ONLY : rhoions, nat, ntyp, ityp, zv, tau
    USE environ_base,      ONLY : tolrhopol, mixrhopol, rhopol, add_jellium
    USE calc_rhopol,       ONLY : iterative_rhopol
    USE periodic,          ONLY : periodic_initions, periodic_initbase, &
                                  periodic_initcell, periodic_clean
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN)   :: nnr, nspin
    REAL(DP), INTENT(IN)  :: rhoelec(nnr,2)
    REAL(DP), INTENT(OUT) :: ejellium, eperiodic
    REAL(DP), INTENT(OUT) :: rhopol_of_V(nnr)
    !
    REAL(DP), ALLOCATABLE :: rhosolute(:)
    !
    REAL(DP), ALLOCATABLE :: rhojellium(:)
    REAL(DP), ALLOCATABLE :: rhopol_jellium(:)
    !
    REAL(DP), ALLOCATABLE :: rhoperiodic(:)
    REAL(DP), ALLOCATABLE :: rhopol_periodic(:)
    !
    REAL(DP), ALLOCATABLE :: rhoaux(:,:)
    REAL(DP), ALLOCATABLE :: vaux(:,:)
    !
    REAL(DP) :: rhoavg
    REAL(DP) :: charge, ehart
    !
    LOGICAL :: reset = .FALSE.
    CHARACTER (LEN=256) :: mixtmp
    !
    ejellium = 0.D0
    eperiodic = 0.D0
    rhopol_of_V = 0.D0
    !
    IF ( mixtype .EQ. 'pcg' .OR. mixtype .EQ. 'psd' ) THEN
       reset = .TRUE.
       mixtmp = mixtype
       mixtype = 'linear'
    ENDIF
    !
    ALLOCATE(rhosolute(nnr))
    rhosolute(:) = rhoions(:) + rhoelec(:,1)
    IF (nspin.EQ.2) rhosolute(:) = rhosolute(:) + rhoelec(:,2)
    rhoavg = SUM(rhosolute)/DBLE(ntot)
    CALL mp_sum( rhoavg, intra_bgrp_comm )
    !
    ! ... If the system is charged in PBC calculates jellium polarization
    !
    IF ( ABS(rhoavg*omega) .GT. 1.D-4 .AND. add_jellium ) THEN
      ALLOCATE( rhojellium( nnr ) )
      ALLOCATE( rhopol_jellium( nnr ) )
      rhojellium = - rhoavg
      rhopol_jellium = 0.D0
      CALL iterative_rhopol(nnr,nspin,maxiter,.FALSE.,tolrhopol/1.D1,mixrhopol,rhojellium,eps,&
                            gradlogeps,rhopol_jellium)
      rhopol_jellium(:) = rhopol_jellium(:) + (1.D0/env_static_permittivity-1.D0/eps(:))*rhoavg
      IF ( verbose .GE. 3 ) CALL write_cube( nnr, rhopol_jellium, 'rhopol_jellium.cube' )
      DEALLOCATE( rhojellium )
      !
      ! Calculates energy due to jellium polarization
      !
      ALLOCATE(vaux(nnr,nspin))
      vaux = 0.D0
      ALLOCATE(rhoaux(nnr,nspin))
      rhoaux( :, 1 ) = rhopol_jellium(:)
      IF ( nspin .EQ. 2 ) rhoaux( :, 2 ) = 0.D0
      CALL v_h_of_rho_r( rhoaux, ehart, charge, vaux )
      ejellium = 0.5D0 * SUM( vaux( : , 1) * rhosolute(:) ) * domega
      CALL mp_sum( ejellium, intra_bgrp_comm )
      DEALLOCATE(vaux,rhoaux)
      !
      rhopol_of_V = rhopol_of_V + rhopol_jellium
      DEALLOCATE(rhopol_jellium)
    END IF
    !
    ! ... Computes the polarization due to periodic images of the solute
    !
    IF ( env_periodicity .NE. 2 ) THEN
      env_periodicity = 0
      CALL periodic_initbase(nnr)
      CALL periodic_initions(nnr, nat, ntyp, ityp, zv, tau, alat, rhoions)
      CALL periodic_initcell(nnr, at)
      ALLOCATE(rhoperiodic(nnr))
      ALLOCATE(rhopol_periodic(nnr))
      rhoperiodic = rhosolute + rhopol - rhopol_of_V
      rhopol_periodic = 0.D0
      CALL iterative_rhopol(nnr,nspin,maxiter,.TRUE.,tolrhopol/1.D1,mixrhopol,rhoperiodic,eps,&
                            gradlogeps,rhopol_periodic)
      IF ( verbose .GE. 3 ) CALL write_cube( nnr, rhopol_periodic, 'rhopol_periodic.cube' )
      DEALLOCATE(rhoperiodic)
      !
      ! Calculates energy due to periodic polarization
      !
      ALLOCATE(vaux(nnr,nspin))
      vaux = 0.D0
      ALLOCATE(rhoaux(nnr,nspin))
      rhoaux( :, 1 ) = rhopol_periodic(:)
      IF ( nspin .EQ. 2 ) rhoaux( :, 2 ) = 0.D0
      CALL v_h_of_rho_r( rhoaux, ehart, charge, vaux )
      eperiodic = 0.5D0 * SUM( vaux( : , 1) * rhosolute(:) ) * domega
      CALL mp_sum( eperiodic, intra_bgrp_comm )
      DEALLOCATE(vaux,rhoaux)
      !
      rhopol_of_V = rhopol_of_V + rhopol_periodic
      DEALLOCATE( rhopol_periodic )
      CALL periodic_clean()
      env_periodicity = 3
    END IF
    !
    DEALLOCATE(rhosolute)
    !
    IF ( reset ) THEN
      reset = .FALSE.
      mixtype = mixtmp
    ENDIF
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE calc_esolvent_of_V
!--------------------------------------------------------------------
!!----------------------------------------------------------------------
!  SUBROUTINE calc_vsr( nnr, charge, spread, pos, vsr )
!!----------------------------------------------------------------------
!      !
!      USE kinds,            ONLY : DP
!      USE io_global,        ONLY : stdout
!      USE cell_base,        ONLY : at, bg, alat, omega
!      USE fft_base,         ONLY : dfftp
!      USE mp,               ONLY : mp_sum
!      USE mp_bands,         ONLY : me_bgrp, intra_bgrp_comm
!      !
!      IMPLICIT NONE
!      !
!      ! ... Declares variables
!      !
!      INTEGER, INTENT(IN)       :: nnr
!      REAL( DP ), INTENT(IN)    :: charge, spread
!      REAL( DP ), INTENT(IN)    :: pos( 3 )
!      REAL( DP ), INTENT(INOUT) :: vsr( nnr )
!      !
!      ! ... Local variables
!      !
!      INTEGER                   :: i, j, k, ir, ir_end, ip
!      INTEGER                   :: index0, ntot
!      !
!      REAL( DP )                :: inv_nr1, inv_nr2, inv_nr3
!      REAL( DP )                :: scale, dist, arg
!      REAL( DP )                :: f1, f2
!      REAL( DP )                :: r( 3 ), s( 3 )
!      REAL( DP ), ALLOCATABLE   :: vlocal ( : )
!      REAL( DP ), EXTERNAL      :: qe_erfc
!      !
!      inv_nr1 = 1.D0 / DBLE( dfftp%nr1 )
!      inv_nr2 = 1.D0 / DBLE( dfftp%nr2 )
!      inv_nr3 = 1.D0 / DBLE( dfftp%nr3 )
!      !
!      index0 = 0
!      ir_end = nnr
!      !
!#if defined (__MPI)
!      DO i = 1, me_bgrp
!        index0 = index0 + dfftp%nr1x*dfftp%nr2x*dfftp%npp(i)
!      END DO
!      ir_end = MIN(nnr,dfftp%nr1x*dfftp%nr2x*dfftp%npp(me_bgrp+1))
!#endif
!      !
!      ntot = dfftp%nr1 * dfftp%nr2 * dfftp%nr3
!      !
!      scale = charge * e2
!      !
!      ALLOCATE( vlocal( nnr ) )
!      vlocal = 0.D0
!      !
!      DO ir = 1, ir_end
!         !
!         ! ... three dimensional indexes
!         !
!         i = index0 + ir - 1
!         k = i / (dfftp%nr1x*dfftp%nr2x)
!         i = i - (dfftp%nr1x*dfftp%nr2x)*k
!         j = i / dfftp%nr1x
!         i = i - dfftp%nr1x*j
!         !
!         DO ip = 1, 3
!            r(ip) = DBLE( i )*inv_nr1*at(ip,1) + &
!                    DBLE( j )*inv_nr2*at(ip,2) + &
!                    DBLE( k )*inv_nr3*at(ip,3)
!         END DO
!         !
!         r(:) = pos(:) - r(:)
!         !
!         ! ... minimum image convention
!         !
!         s(:) = MATMUL( r(:), bg(:,:) )
!         s(:) = s(:) - ANINT(s(:))
!         r(:) = MATMUL( at(:,:), s(:) )
!         !
!         dist = SQRT(SUM( r * r )) * alat
!         !
!         IF ( dist .LT. 1.D-8 ) THEN
!           vlocal( ir ) = 0.D0 ! should not be used for short distances, where it diverges
!         ELSE
!           vlocal( ir ) = qe_erfc( dist / spread ) / dist
!         END IF
!         !
!      END DO
!      !
!      ! ... multiply by charge
!      !
!      vlocal = vlocal * scale
!      !
!      vsr = vsr + vlocal
!      DEALLOCATE( vlocal )
!      !
!      RETURN
!      !
!!----------------------------------------------------------------------
!      END SUBROUTINE calc_vsr
!!----------------------------------------------------------------------
!--------------------------------------------------------------------
END MODULE solvent
!--------------------------------------------------------------------
