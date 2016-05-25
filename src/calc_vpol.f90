!
! Copyright (C) 2007-2008 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! original version by O. Andreussi and N. Marzari
! based on the precodintioning algorithm of
! G. Fisicaro, L. Genovese and S. Goedecker
!
MODULE calc_vpol

  USE kinds,        ONLY : DP
  USE constants,    ONLY : fpi, e2
  USE io_global,    ONLY : stdout
  USE mp,           ONLY : mp_sum
  USE mp_bands,     ONLY : intra_bgrp_comm
  USE environ_base, ONLY : verbose, environ_unit, env_periodicity,&
                           mixtype
  USE environ_cell, ONLY : ntot
  USE periodic,     ONLY : calc_vperiodic
  USE fd_gradient,  ONLY : init_fd_gradient, calc_fd_gradient
  USE environ_debug, ONLY : write_cube
  USE control_flags, ONLY : tddfpt

  LOGICAL :: initialized
  INTEGER :: nfdpoint
  INTEGER :: ncfd
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: residual, Ap, Apin, p, pin
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: vsolute, vperiodic
  REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: vaux

  SAVE

  PRIVATE

  PUBLIC :: iterative_vpol
  
CONTAINS

  SUBROUTINE vpol_init(nnr, nspin)

    IMPLICIT NONE

    INTEGER, INTENT(IN)     :: nnr, nspin

    initialized = .TRUE.

    ALLOCATE( residual( nnr ) ) 
    residual = 0.D0
    ALLOCATE( Ap( nnr ) )
    Ap = 0.D0
    ALLOCATE( Apin( nnr ) )
    Apin = 0.D0
    ALLOCATE( p( nnr ) ) 
    p = 0.D0
    ALLOCATE( pin( nnr ) )
    pin = 0.D0
    ALLOCATE( vsolute( nnr ) ) 
    vsolute = 0.D0
    ALLOCATE( vperiodic( nnr ) ) 
    vperiodic = 0.D0
    ALLOCATE( vaux( nnr, nspin ) ) 
    vaux = 0.D0

    RETURN

  END SUBROUTINE vpol_init


  SUBROUTINE vpol_clean

    IMPLICIT NONE

    initialized = .FALSE.

    IF ( ALLOCATED( residual  ) ) DEALLOCATE( residual )
    IF ( ALLOCATED( Ap        ) ) DEALLOCATE( Ap )
    IF ( ALLOCATED( Apin      ) ) DEALLOCATE( Apin )
    IF ( ALLOCATED( p         ) ) DEALLOCATE( p )
    IF ( ALLOCATED( pin       ) ) DEALLOCATE( pin )
    IF ( ALLOCATED( vsolute   ) ) DEALLOCATE( vsolute )
    IF ( ALLOCATED( vperiodic ) ) DEALLOCATE( vperiodic )
    IF ( ALLOCATED( vaux      ) ) DEALLOCATE( vaux )

    RETURN

  END SUBROUTINE vpol_clean

  SUBROUTINE iterative_vpol(nnr,nspin,maxiter,periodic,tol,rho,invsqrteps,factsqrteps,rhopol,vsolvent)
  
    IMPLICIT NONE

    ! ... Declare variables
    
    INTEGER, INTENT(IN)     :: nnr
    INTEGER, INTENT(IN)     :: nspin
    INTEGER, INTENT(IN)     :: maxiter
    LOGICAL, INTENT(IN)     :: periodic
    REAL(DP), INTENT(IN)    :: tol
    REAL(DP), INTENT(IN)    :: rho(nnr)
    REAL(DP), INTENT(IN)    :: invsqrteps(nnr)
    REAL(DP), INTENT(IN)    :: factsqrteps(nnr)
    REAL(DP), INTENT(INOUT) :: rhopol(nnr)
    REAL(DP), INTENT(INOUT) :: vsolvent(nnr)

    ! ... Local variables

    LOGICAL, SAVE :: first=.TRUE.
    INTEGER  :: iter, icor
    REAL(DP) :: alpha, beta, deltar, sr, srin, pAp, ehart, charge
    REAL(DP), ALLOCATABLE :: gvaux(:,:),gvtmp(:,:)
    REAL(DP), ALLOCATABLE, SAVE :: vtot(:)

    CALL start_clock( 'calc_vpol' ) 

    IF ( .NOT. initialized ) CALL vpol_init(nnr,nspin)

    vaux = 0.D0
    CALL v_h_of_rho_r( rho, ehart, charge, vaux )
    vsolute(:) = vaux(:,1)
    IF ( verbose .GE. 3 ) CALL write_cube( nnr, vsolute, 'vsolute.cube' )

    IF ( first ) THEN
       first = .FALSE.
       residual = rho
       vsolvent = 0.D0
       IF (ALLOCATED(vtot)) DEALLOCATE(vtot)
       ALLOCATE(vtot(nnr))
       vtot = 0.D0
    ELSE
       vsolvent = vtot ! this is the total potential from the previous step
       residual = rho - factsqrteps * vsolvent ! not sure why using this 
       ! apply the preconditioner to this modified residual
       Ap = residual * invsqrteps
       vaux = 0.D0
       CALL v_h_of_rho_r( Ap, ehart, charge, vaux )
       vperiodic = 0.D0
       CALL calc_vperiodic( nnr, nspin, .FALSE., Ap, vperiodic )
       p = ( vaux(:,1) + vperiodic(:) ) * invsqrteps(:)
       ! now the new residual is given by 
       residual = factsqrteps * ( vsolvent - p )  
       ! check if the guess is reasonable enough
       deltar = SUM(residual(:)*residual(:))
       CALL mp_sum( deltar , intra_bgrp_comm )
       deltar = SQRT(deltar/DBLE(ntot))
       IF ( deltar .LT. 1.D-02 ) THEN
         vsolvent = p 
       ELSE  
         ! otherwise reset to the null guess
         IF ( verbose .GE. 1 ) WRITE(environ_unit,9001)deltar
         residual = rho
         vsolvent = 0.D0
       ENDIF
    ENDIF
    
    DO iter = 0, maxiter
       
       IF ( verbose .GE. 1 ) WRITE(environ_unit,9002) iter
       Apin = Ap
       pin = p
       srin = sr

       ! ... Preconditioner steps: residual' = residual / sqrteps
       !                           laplace v' = - fpi * residual'
       !                           vprecond = v' / sqrteps

       Ap = residual * invsqrteps
       vaux = 0.D0
       CALL v_h_of_rho_r( Ap, ehart, charge, vaux )
       vperiodic = 0.D0
       CALL calc_vperiodic( nnr, nspin, .FALSE., Ap, vperiodic )
       p = ( vaux(:,1) + vperiodic(:) ) * invsqrteps(:)

       ! ... Apply generalized poisson operator to vprecond
       !     A * v = (\nabla \cdot eps \nabla) * v
       !     which for vprecond becomes
       !     A * s  = -sqrteps*laplsqrteps*s - fpi*residual

       Ap(:) = factsqrteps(:)*p(:) + residual(:)
       
       ! ... Conjugate gradient or steepest descent input

       sr = SUM( p * residual )
       CALL mp_sum( sr, intra_bgrp_comm )
       IF ( sr .LT. 1.D-30 ) THEN
          WRITE(stdout,*)'ERROR: zero energy step in polarization potential iteration'
          STOP
       ENDIF
       IF ( mixtype .EQ. 'pcg' .AND. iter .GT. 0 ) THEN 
         beta = sr / srin
         p = p + beta * pin
         Ap = Ap + beta * Apin
       END IF 
       
       ! ... Step downhill
       
       pAp = SUM( p * Ap )
       CALL mp_sum( pAp, intra_bgrp_comm ) 
       alpha = sr / pAp
       vsolvent = vsolvent + alpha * p
       residual = residual - alpha * Ap

       ! ... Some debugging print out

       IF ( verbose .GE. 1 ) WRITE(environ_unit,*)'alpha = ',alpha,' beta = ',beta
       IF ( verbose .GE. 2 ) WRITE(environ_unit,*)'sr = ',sr,' srin = ',srin,' pAp = ',pAp
       IF ( verbose .GE. 4 ) CALL write_cube( nnr, vsolvent, 'vsolvtmp.cube' )
       IF ( verbose .GE. 4 ) CALL write_cube( nnr, residual, 'residtmp.cube' )
       
       ! ... If polarization potential is converged exit

       deltar = SUM(residual(:)*residual(:))
       CALL mp_sum( deltar , intra_bgrp_comm )
       deltar = SQRT(deltar/DBLE(ntot))
       IF ( verbose .GE. 1 ) WRITE(environ_unit,9004)deltar,tol
       IF ( deltar .LT. tol .AND. iter .GT. 0 ) THEN
         IF ( verbose .GE. 1 ) WRITE(environ_unit,9005)
         EXIT
       ELSE IF ( iter .EQ. maxiter ) THEN
         WRITE(stdout,9006)
       ENDIF

    ENDDO
    IF (.not.tddfpt.AND.verbose.GE.1) WRITE(stdout, 9000) deltar, iter
    
    vtot = vsolvent
    vsolvent = vsolvent - vsolute
    IF ( verbose .GE. 3 ) CALL write_cube( nnr, vsolvent, 'vsolvent.cube' )
    CALL vpol_clean()

    CALL stop_clock( 'calc_vpol' )
    RETURN
    
9000 FORMAT('     polarization accuracy =',1PE8.1,', # of iterations = ',i3)
9001 FORMAT(' Warning: bad guess with residual norm = ',E14.6,', reset to no guess')
9002 FORMAT(' Iteration # ',i10)
9004 FORMAT(' deltar = ',E14.6,' tol = ',E14.6)
9005 FORMAT(' Potential is converged, exit!')
9006 FORMAT(' Warning: Polarization potential not converged')

  END SUBROUTINE iterative_vpol
!
END MODULE calc_vpol
