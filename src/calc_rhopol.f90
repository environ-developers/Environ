!
! Copyright (C) 2007-2008 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! original version by O. Andreussi and N. Marzari
!
MODULE calc_rhopol

  USE kinds,         ONLY : DP
  USE constants,     ONLY : fpi
  USE io_global,     ONLY : stdout
  USE mp,            ONLY : mp_sum
  USE mp_bands,      ONLY : intra_bgrp_comm
  USE environ_base,  ONLY : verbose, environ_unit, e2, env_periodicity, &
                           mixtype, ndiis, env_static_permittivity
  USE environ_cell,  ONLY : domega, ntot, omega
  USE periodic,      ONLY : calc_gradvperiodic
  USE environ_debug, ONLY : write_cube
  USE control_flags, ONLY : tddfpt

  LOGICAL :: initialized
  REAL(DP), ALLOCATABLE :: gvtot(:,:)
  REAL(DP), ALLOCATABLE :: rhotot(:), rhozero(:), vtot(:,:)
  REAL(DP), ALLOCATABLE :: rhoin(:,:), rhoout(:,:), residual(:,:), diismat(:,:)

  SAVE

  PRIVATE

  PUBLIC :: iterative_rhopol, rhopol_of_v

CONTAINS

  SUBROUTINE rhopol_init(nnr, nspin)

    IMPLICIT NONE

    INTEGER, INTENT(IN)     :: nnr, nspin

    initialized = .TRUE.

    SELECT CASE ( TRIM( mixtype ) )
    CASE ('linear')
      ndiis = 1
    CASE ('anderson')
      ndiis = 2
    CASE ('diis')
      ndiis = MAX(2,ndiis)
    CASE default
      WRITE(*,*)'ERROR: wrong mixtype in Environ namelist'
    END SELECT

    ALLOCATE( rhozero( nnr ) )
    ALLOCATE( vtot ( nnr, nspin ) )
    ALLOCATE( gvtot ( 3, nnr ) )
    ALLOCATE( rhotot ( nnr ) )
    ALLOCATE( rhoin ( ndiis, nnr ) )
    ALLOCATE( rhoout ( ndiis, nnr ) )
    ALLOCATE( residual ( ndiis, nnr ) )
    ALLOCATE( diismat ( ndiis, ndiis ) )
    rhozero = 0.D0
    vtot = 0.D0
    rhotot = 0.D0
    rhoin = 0.D0
    rhoout = 0.D0
    residual = 0.D0
    diismat = 0.D0

    RETURN

  END SUBROUTINE rhopol_init


  SUBROUTINE rhopol_clean

    IMPLICIT NONE

    initialized = .FALSE.

    IF ( ALLOCATED(rhozero) )  DEALLOCATE( rhozero )
    IF ( ALLOCATED(vtot) )     DEALLOCATE( vtot )
    IF ( ALLOCATED(gvtot) )    DEALLOCATE( gvtot )
    IF ( ALLOCATED(rhotot) )   DEALLOCATE( rhotot )
    IF ( ALLOCATED(rhoin) )    DEALLOCATE( rhoin )
    IF ( ALLOCATED(rhoout) )   DEALLOCATE( rhoout )
    IF ( ALLOCATED(residual) ) DEALLOCATE( residual )
    IF ( ALLOCATED(diismat) )  DEALLOCATE( diismat )

    RETURN

  END SUBROUTINE rhopol_clean

  SUBROUTINE iterative_rhopol(nnr,nspin,maxiter,periodic,tol,mix,rho,eps,gradlogeps,rhoiter)

    IMPLICIT NONE

    ! ... Declare variables

    INTEGER, INTENT(IN)     :: nnr
    INTEGER, INTENT(IN)     :: nspin
    INTEGER, INTENT(IN)     :: maxiter
    LOGICAL, INTENT(IN)     :: periodic
    REAL(DP), INTENT(IN)    :: tol
    REAL(DP), INTENT(IN)    :: mix
    REAL(DP), INTENT(IN)    :: rho(nnr)
    REAL(DP), INTENT(IN)    :: eps(nnr)
    REAL(DP), INTENT(IN)    :: gradlogeps(3,nnr)
    REAL(DP), INTENT(INOUT) :: rhoiter(nnr)

    ! ... Local variables

    INTEGER  :: iter, ir, ind
    REAL(DP) :: total, totpol, totzero, totiter, deltarho, ehart, charge

    CALL start_clock( 'get_rhopol' )

    IF ( .NOT. initialized ) CALL rhopol_init(nnr,nspin)

    total = SUM(rho)*domega
    CALL mp_sum( total, intra_bgrp_comm )
    totpol = (1.D0 - env_static_permittivity)/env_static_permittivity * total
    CALL mp_sum( totpol, intra_bgrp_comm )
    rhozero(:) = rho(:)/eps(:)
    IF (periodic) THEN
       gvtot = 0.D0
       CALL calc_gradvperiodic( nnr, rho, gvtot )
       DO ir = 1, nnr
         rhozero(ir) = - SUM(gvtot(:,ir)*gradlogeps(:,ir))/fpi/e2
       ENDDO
    ENDIF
    totzero = SUM((1.D0-eps(:))*rhozero(:))*domega
    CALL mp_sum( totzero, intra_bgrp_comm )
    totiter = SUM(rhoiter)*domega
    CALL mp_sum( totiter, intra_bgrp_comm )
    IF ( verbose .GE. 1 ) WRITE(environ_unit,9001) totiter

    CALL v_h_of_rho_r( rho, ehart, charge, vtot )
    IF ( verbose .GE. 3 ) CALL write_cube( nnr, vtot( :, 1), 'vsolute.cube' )

    DO iter = 1, maxiter

       IF ( verbose .GE. 1 ) WRITE(environ_unit,9002) iter
       ind = MOD(iter-1,ndiis)+1
       rhoin(ind,:) = rhoiter(:)

       ! ...Compute total charge density

       rhotot = rhozero + rhoiter

       ! ... Compute gradient of potential from total density

       CALL gradv_h_of_rho_r ( rhotot, gvtot )

       ! ... If partially periodic system add pbc correction to the gradient

       IF ( env_periodicity .NE. 3 ) CALL calc_gradvperiodic( nnr, rhotot, gvtot )

       ! ... Compute polarization charge from grad(V) and grad(eps)

       rhoiter = 0.D0
       DO ir = 1, nnr
         rhoiter(ir) = SUM(gvtot(:,ir)*gradlogeps(:,ir))/fpi/e2
       END DO
       rhoout(ind,:) = rhoiter(:)
       residual(ind,:) = rhoout(ind,:) - rhoin(ind,:)
       IF ( verbose .GE. 5 ) CALL write_cube( nnr, rhoout(ind,:), 'rhoiter.cube' )

       ! ... Mix polarization charge with previous iterations

       CALL mix_rhopol(iter, ind, nnr, ndiis, mix, mixtype, rhoin,&
                       rhoout, residual, diismat, rhoiter)

       ! ... If polarization charge is converged exit

       deltarho = SUM(residual(ind,:)*residual(ind,:))
       CALL mp_sum( deltarho , intra_bgrp_comm )
       deltarho = SQRT(deltarho)/DBLE(ntot)
       totiter = SUM(rhoiter)*domega
       CALL mp_sum( totiter, intra_bgrp_comm )
       IF ( verbose .GE. 1 ) WRITE(environ_unit,9004)deltarho,tol
       IF ( verbose .GE. 2 ) WRITE(environ_unit,9003)totiter,totzero,totpol,total
       IF (deltarho.LT.tol) THEN
         IF ( verbose .GE. 1 ) WRITE(environ_unit,9005)
         EXIT
       ELSE IF ( iter .EQ. maxiter ) THEN
         WRITE(stdout,9006)
       ENDIF

    ENDDO
    IF (.not.tddfpt.AND.verbose.GE.1) WRITE(stdout, 9000) deltarho, iter

    IF (periodic) THEN
      rhoiter = rhoiter + rhozero
    ENDIF

    CALL rhopol_clean()

    CALL stop_clock( 'get_rhopol' )

    RETURN

9000 FORMAT('     polarization accuracy =',1PE8.1,', # of iterations = ',i3)
9001 FORMAT(' Starting from polarization: rhoiter = ',F13.6)
9002 FORMAT(' Iteration # ',i10)
9003 FORMAT(' Total iterative polarization charge = ',4F13.6)
9004 FORMAT(' deltarho = ',E14.6,' tol = ',E14.6)
9005 FORMAT(' Charges are converged, exit!')
9006 FORMAT(' Warning: Polarization charge not converged')

  END SUBROUTINE iterative_rhopol
!
!--------------------------------------------------------------------
      SUBROUTINE mix_rhopol(iter, ind, nnr, ndiis, mix, mixing, rhoin, &
                            rhoout, residual, diismat, rhopol)
!--------------------------------------------------------------------

      USE  kinds,          ONLY : DP
      USE  environ_base,   ONLY : verbose

      IMPLICIT NONE

      ! Declares variables

      INTEGER, INTENT(IN)            :: iter, ind, nnr, ndiis
      REAL(DP), INTENT(IN)           :: mix
      REAL(DP), INTENT(IN)           :: rhoin(ndiis,nnr)
      REAL(DP), INTENT(IN)           :: rhoout(ndiis,nnr)
      REAL(DP), INTENT(IN)           :: residual(ndiis,nnr)
      REAL(DP), INTENT(OUT)          :: rhopol(nnr)
      REAL(DP), INTENT(INOUT)        :: diismat( ndiis , ndiis )
      CHARACTER(LEN=256), INTENT(IN) :: mixing

      ! Local variables

      INTEGER                        :: i, j, ndim
      REAL(DP)                       :: beta
      REAL(DP), ALLOCATABLE          :: cdiis( : )
!      REAL(DP), ALLOCATABLE          :: diismat( : , : )

      CALL start_clock( 'mixrhopol' )

      IF (TRIM(mixing).EQ.'linear') THEN

        IF (verbose.GE.1) WRITE(11,'(1x,a,f12.6)')'Linear mixing, mix =  ',mix

        ! ... Linear mixing of polarization charges

        rhopol(:) = mix * rhoout(1,:) + (1.D0-mix) * rhoin(1,:)

      ELSEIF (TRIM(mixing).EQ.'anderson') THEN

        IF (verbose.GE.1) WRITE(11,'(1x,a,f12.6)')'Anderson mixing, mix = ',mix

        ! ... Anderson mixing of polarization charges

        beta = SUM(residual(1,:)*(residual(1,:)-residual(2,:)))
        beta = beta/SUM(residual(1,:)*residual(2,:))
        rhopol = mix * ((1.D0-beta)*rhoout(1,:)+beta*rhoout(2,:))
        rhopol = rhopol + (1.D0-mix)*((1.D0-beta)*rhoin(1,:)         &
                 + beta*rhoin(2,:))

      ELSEIF (TRIM(mixing).EQ.'diis') THEN

        IF (verbose.GE.1) WRITE(11,'(1x,a,i12)')'DIIS mixing, ndiis = ',ndiis

        ! ... DIIS mixing of polarization charges

        ndim=MIN(iter,ndiis)
        DO i = 1,ndim
          diismat(ind,i) = SUM(residual(ind,:)*residual(i,:))
          diismat(i,ind) = diismat(ind,i)
        ENDDO
        IF (ndim.EQ.1) THEN
          rhopol(:) = mix * rhoout(1,:) + (1.D0-mix) * rhoin(1,:)
        ELSE
          ALLOCATE(cdiis(ndim))
          CALL solvediismat(ndiis,ndim,diismat,cdiis)
          rhopol = 0.D0
          DO i = 1, ndim
            rhopol = rhopol + cdiis(i) * rhoout(i,:)
          ENDDO
          DEALLOCATE(cdiis)
        ENDIF

      ELSE

        WRITE(*,*)'ERROR: unknown keyword for polarization mixing'
        STOP

      ENDIF

100   CONTINUE
      CALL stop_clock( 'mixrhopol' )
      RETURN

      END SUBROUTINE mix_rhopol

!--------------------------------------------------------------------
      SUBROUTINE solvediismat(n,l,mat,c)
!--------------------------------------------------------------------

      USE  kinds,         ONLY : DP

      IMPLICIT NONE

      INTEGER, PARAMETER :: nb = 20

      INTEGER, INTENT(IN) :: n, l
      REAL(DP), INTENT(IN) :: mat(n,n)
      REAL(DP), INTENT(OUT) :: c(l)

      INTEGER :: m, info
      INTEGER, ALLOCATABLE :: ipiv(:)
      REAL(DP), ALLOCATABLE :: a(:,:)
      REAL(DP), ALLOCATABLE :: b(:)
      REAL(DP), ALLOCATABLE :: work(:)

      LOGICAL, SAVE :: first=.true.
      INTEGER :: lwork

      lwork = nb * n
      c = 0.D0
      IF (l.LT.2) THEN
        WRITE(*,*)'ERROR: WRONG MATRIX DIMENSION IN POLARIZATION DIIS'
        STOP
      ELSE IF (l.EQ.2) THEN
        c(2) = (mat(1,2) - mat(1,1))/(2.D0*mat(1,2)-mat(1,1)-mat(2,2))
        c(1) = 1.D0 - c(2)
      ELSE
        m=l+1
        ALLOCATE(a(m,m))
        a=-1.D0
        a(1:l,1:l)=mat(1:l,1:l)
        a(m,m)=0.D0
        ALLOCATE(b(m))
        b=0.D0
        b(m)=-1.D0
        ALLOCATE(ipiv(m))
        ALLOCATE(work(lwork))
        CALL dsysv( 'U', m, 1, a, m, ipiv, b, m, work, lwork, info )
        DEALLOCATE(ipiv)
        DEALLOCATE(work)
        c(1:l)=b(1:l)
        DEALLOCATE(a)
        DEALLOCATE(b)
      ENDIF

      RETURN

!--------------------------------------------------------------------
      END SUBROUTINE solvediismat
!--------------------------------------------------------------------
!--------------------------------------------------------------------
      SUBROUTINE rhopol_of_v(nnr,nspin,eps,gradlogeps,rhosolute,vsolvent,rhopol)
!--------------------------------------------------------------------

      IMPLICIT NONE
      !
      ! Input variables
      !
      INTEGER, INTENT(IN) :: nnr, nspin
      REAL(DP), INTENT(IN) :: eps(nnr), gradlogeps(3,nnr)
      REAL(DP), INTENT(IN) :: rhosolute(nnr)
      REAL(DP), INTENT(IN) :: vsolvent(nnr)
      !
      ! Output variables
      !
      REAL(DP), INTENT(OUT) :: rhopol(nnr)
      !
      ! Local variables
      !
      INTEGER :: ir
      REAL(DP) :: charge, ehart
      REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: rhoaux, vaux, gvtot

      rhopol = 0.D0

      ALLOCATE( rhoaux( nnr, nspin ) )
      rhoaux( :, 1 ) = rhosolute(:)
      IF ( nspin .EQ. 2 ) rhoaux( :, 2 ) = 0.D0

      ALLOCATE( vaux( nnr, nspin ) )
      vaux = 0.D0
      CALL v_h_of_rho_r( rhoaux, ehart, charge, vaux )
      DEALLOCATE( rhoaux )
      vaux( :, 1 ) = vaux( :, 1 ) + vsolvent(:)

      ALLOCATE( gvtot( 3, nnr ) )
      gvtot = 0.D0
      CALL external_gradient( vaux( :, 1 ), gvtot )
      DEALLOCATE( vaux )

      DO ir = 1, nnr
         rhopol(ir) = SUM(gvtot(:,ir)*gradlogeps(:,ir)) / fpi / e2
      ENDDO
      DEALLOCATE( gvtot )

      rhopol(:) = rhopol(:) + (1.D0-eps(:))/eps(:)*rhosolute(:)
      IF ( verbose .GE. 2 ) CALL write_cube( nnr, rhopol, 'rhopol.cube' )
      !
      RETURN
      !
!--------------------------------------------------------------------
      END SUBROUTINE rhopol_of_v
!--------------------------------------------------------------------

END MODULE calc_rhopol
