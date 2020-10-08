MODULE core_fft
  !
  USE modules_constants, ONLY : DP, e2, tpi, fpi
  USE fft_interfaces, ONLY : fwfft, invfft
  USE core_types
  USE environ_types
  !
  PRIVATE
  !
  PUBLIC :: poisson_fft, gradpoisson_fft, force_fft, convolution_fft, gradient_fft, graddot_fft, laplacian_fft, hessian_fft, &
            & field_of_gradrho, hessv_h_of_rho_r
  !
  INTERFACE convolution_fft
     MODULE PROCEDURE convolution_fft_density, convolution_fft_gradient, convolution_fft_hessian
  END INTERFACE convolution_fft
  !
CONTAINS
!--------------------------------------------------------------------
  SUBROUTINE poisson_fft( fft, fin, fout )
!--------------------------------------------------------------------
    !
    ! Solves the Poisson equation \nabla \fout(r) = - 4 \pi \fin(r)
    ! Input and output functions are defined in real space
    !
    USE correction_mt, ONLY : calc_vmt
    !
    IMPLICIT NONE
    !
    TYPE( fft_core ), TARGET, INTENT(IN) :: fft
    TYPE( environ_density ), INTENT(IN) :: fin
    TYPE( environ_density ), INTENT(INOUT) :: fout
    !
    INTEGER :: ig
    COMPLEX(DP), DIMENSION(:), ALLOCATABLE :: auxr, auxg, vaux
    !
    ! ... local aliases
    !
    INTEGER, POINTER :: gstart, ngm
    REAL(DP), POINTER :: tpiba2, omega
    REAL(DP), DIMENSION(:), POINTER :: gg
    TYPE(fft_type_descriptor), POINTER :: dfft
    !
    ! ... add tests for compatilibity between input, output, and fft_core
    !
    tpiba2 => fft % cell % tpiba2
    omega => fft % cell % omega
    gstart => fft % gstart
    ngm => fft % ngm
    gg => fft % gg
    dfft => fft % dfft
    !
    ! Bring fin%of_r from R-space --> G-space
    !
    ALLOCATE( auxr( dfft%nnr ) )
    auxr = CMPLX( fin%of_r, 0.D0, kind=DP )
    CALL fwfft( 'Rho', auxr, dfft )
    ALLOCATE( auxg( ngm ) )
    auxg = auxr(dfft%nl(:))
    !
    auxr = CMPLX( 0.D0, 0.D0, kind=DP )
    !
!$omp parallel do
    DO ig = gstart, ngm
       !
       auxr(dfft%nl(ig)) = auxg(ig) / gg(ig)
       !
    ENDDO
!$omp end parallel do
    !
    auxr = auxr * e2 * fpi / tpiba2
    !
    IF ( fft%use_internal_pbc_corr ) THEN
       ALLOCATE( vaux( ngm ) )
       CALL calc_vmt( fft, auxg, vaux )
       auxr(dfft%nl(1:ngm)) = auxr(dfft%nl(1:ngm)) + vaux(1:ngm)
       DEALLOCATE( vaux )
    END IF
    !
    IF ( dfft%lgamma ) THEN
       auxr( dfft%nlm(:) ) = CMPLX( REAL(auxr(dfft%nl(:))), -AIMAG(auxr(dfft%nl(:))), kind=DP )
    END IF
    !
    ! ... transform hartree potential to real space
    !
    CALL invfft ('Rho', auxr, dfft)
    !
    fout % of_r(:) = DBLE(auxr(:))
    !
    DEALLOCATE( auxr )
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE poisson_fft
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE gradpoisson_fft( fft, fin, gout )
!--------------------------------------------------------------------
    !
    ! Solves the Poisson equation \grad \gout(r) = - 4 \pi \fin(r)
    ! where gout is the gradient of the potential
    ! Input and output functions are defined in real space
    !
    USE correction_mt, ONLY : calc_gradvmt
    !
    IMPLICIT NONE
    !
    TYPE( fft_core ), TARGET, INTENT(IN) :: fft
    TYPE( environ_density ), INTENT(IN) :: fin
    TYPE( environ_gradient ), INTENT(INOUT) :: gout
    !
    INTEGER :: ipol, ig
    REAL(DP) :: fac
    COMPLEX(DP), DIMENSION(:), ALLOCATABLE :: auxr, auxg, vaux
    !
    ! ... local aliases
    !
    INTEGER, POINTER :: ngm, gstart
    REAL(DP), POINTER :: tpiba, omega
    REAL(DP), DIMENSION(:), POINTER :: gg
    REAL(DP), DIMENSION(:,:), POINTER :: g
    TYPE(fft_type_descriptor), POINTER :: dfft
    !
    ! ... add tests for compatilibity between input, output, and fft
    !
    tpiba => fft % cell % tpiba
    omega => fft % cell % omega
    ngm => fft % ngm
    gstart => fft % gstart
    gg => fft % gg
    g => fft % g
    dfft => fft % dfft
    !
    ! ... Bring rho to G space
    !
    ALLOCATE( auxr( dfft%nnr ) )
    auxr( : ) = CMPLX( fin%of_r( : ), 0.D0, KIND=dp )
    CALL fwfft('Rho', auxr, dfft)
    ALLOCATE( auxg( ngm ) )
    auxg = auxr(dfft%nl(:))
    !
    ! ... Compute gradient of potential in G space one direction at a time
    !
    DO ipol = 1, 3
       !
       auxr(:) = CMPLX(0.0_dp,0.0_dp)
       !
!$omp parallel do private(fac)
       DO ig = gstart, ngm
          !
          auxr(dfft%nl(ig)) = CMPLX(-AIMAG(auxg(ig)),REAL(auxg(ig),kind=dp)) * g(ipol,ig) / gg(ig)
          !
       END DO
!$omp end parallel do
       !
       ! ...and add the factor e2*fpi/2\pi/a coming from the missing prefactor of
       !  V = e2 * fpi divided by the 2\pi/a factor missing in G
       !
       fac = e2 * fpi / tpiba
       auxr = auxr * fac
       !
       ! ...add martyna-tuckerman correction, if needed
       !
       IF ( fft%use_internal_pbc_corr ) THEN
          ALLOCATE( vaux( ngm ) )
          CALL calc_gradvmt( ipol, fft, auxg, vaux)
          auxr(dfft%nl(:)) = auxr(dfft%nl(:)) + vaux(:)
          DEALLOCATE( vaux )
       END IF
       !
       ! Assuming GAMMA ONLY
       !
       IF ( dfft%lgamma ) THEN
          auxr(dfft%nlm(:)) = &
               CMPLX( REAL( auxr(dfft%nl(:)) ), -AIMAG( auxr(dfft%nl(:)) ) ,kind=DP)
       END IF
       !
       ! ... bring back to R-space, (\grad_ipol a)(r) ...
       !
       CALL invfft('Rho', auxr, dfft)
       !
       gout%of_r(ipol,:) = REAL( auxr(:) )
       !
    ENDDO
    !
    DEALLOCATE(auxr)
    DEALLOCATE(auxg)
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE gradpoisson_fft
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE force_fft( fft, rho, ions, nat, force )
!--------------------------------------------------------------------
    !
    !
    USE mp, ONLY : mp_sum
    USE correction_mt, ONLY : calc_fmt
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: nat
    TYPE( fft_core ), TARGET, INTENT(IN) :: fft
    TYPE( environ_density ), INTENT(IN) :: rho
    TYPE( environ_ions ), INTENT(IN) :: ions
    REAL(DP), DIMENSION(3,nat), INTENT(OUT) :: force
    !
    INTEGER :: iat, ig, ityp
    REAL(DP) :: fact, arg
    COMPLEX(DP), DIMENSION(:), ALLOCATABLE :: auxr, auxg
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: vloc
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: ftmp
    !
    ! ... local aliases
    !
    INTEGER, POINTER :: ngm, gstart
    REAL(DP), POINTER :: tpiba, omega
    REAL(DP), DIMENSION(:,:), POINTER :: g
    TYPE(fft_type_descriptor), POINTER :: dfft
    !
    ! ... add tests for compatilibity between input, output, and fft
    !
    tpiba => fft % cell % tpiba
    omega => fft % cell % omega
    ngm => fft % ngm
    gstart => fft % gstart
    g => fft % g
    dfft => fft % dfft
    !
    ALLOCATE( auxr( dfft%nnr ) )
    !
    ! Bring vloc from R space to G space
    !
    ALLOCATE( vloc( ngm, ions%ntyp ) )
    DO ityp = 1, ions%ntyp
       !
       auxr = CMPLX( ions%vloc(ityp)%of_r, 0.D0, KIND=dp )
       CALL fwfft( 'Rho', auxr, dfft )
       vloc( :, ityp) = auxr(dfft%nl(:))
       !
    ENDDO
    !
    ! ... Bring rho to G space
    !
    auxr( : ) = CMPLX( rho%of_r( : ), 0.D0, KIND=dp )
    CALL fwfft('Rho', auxr, dfft)
    ALLOCATE( auxg( ngm ) )
    auxg = auxr(dfft%nl(:))
    DEALLOCATE( auxr )
    !
    ! aux contains now n(G)
    !
    IF (dfft%lgamma) THEN
       fact = 2.d0
    ELSE
       fact = 1.d0
    ENDIF
    !
    DO iat = 1, nat
       !
       force( :, iat ) = 0.D0
       !
       DO ig = gstart, ngm
          arg = tpi * SUM ( g(:,ig) * ions%tau(:,iat) )
          force( :, iat ) = force( :, iat ) + g( :, ig ) * vloc(ig,ions%ityp(iat)) * &
               ( SIN(arg)*DBLE(auxg(ig)) + COS(arg)*AIMAG(auxg(ig)) )
       ENDDO
       !
       force( :, iat ) = fact * force( :, iat ) * omega * tpiba
       !
    ENDDO
    WRITE(*,*)'forces lc'
    DO iat = 1, nat
       WRITE(*,'(i8,3f10.4)')iat,force(:,iat)
    ENDDO
    !
    ! ...add martyna-tuckerman correction, if needed
    !
    IF ( fft % use_internal_pbc_corr ) THEN
       ALLOCATE( ftmp( 3, nat ) )
       CALL calc_fmt( fft, auxg, ions, ftmp )
       force = force + fact * ftmp
       WRITE(*,*)'forces mt',fact
       DO iat = 1, nat
          WRITE(*,'(i8,3f10.4)')iat,fact*ftmp(:,iat)
       ENDDO
       DEALLOCATE( ftmp )
    ENDIF
    !
    CALL mp_sum( force, rho%cell%comm )
    !
    DEALLOCATE(auxg)
    DEALLOCATE(vloc)
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE force_fft
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE convolution_fft_density( fft, fa, fb, fc )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    TYPE( fft_core ), TARGET, INTENT(IN) :: fft
    TYPE( environ_density ), INTENT(IN) :: fa, fb
    TYPE( environ_density ), INTENT(INOUT) :: fc
    !
    COMPLEX( DP ), DIMENSION( : ), ALLOCATABLE :: auxr, auxg
    !
    ! ... local aliases
    !
    REAL(DP), POINTER :: omega
    TYPE(fft_type_descriptor), POINTER :: dfft
    !
    ! ... add tests for compatilibity between input, output, and fft
    !
    dfft => fft % dfft
    omega => fft % cell % omega
    !
    ! Bring fa and fb to reciprocal space
    !
    ALLOCATE( auxr( dfft%nnr ) )
    auxr(:) = CMPLX( fa%of_r(:), 0.D0, kind=DP )
    CALL fwfft('Rho', auxr, dfft)
    !
    ALLOCATE( auxg( dfft%nnr ) )
    auxg = 0.D0
    !
    auxg(dfft%nl(:)) = auxr(dfft%nl(:))
    !
    auxr(:) = CMPLX( fb%of_r(:), 0.D0, kind=DP )
    CALL fwfft('Rho', auxr, dfft)
    !
    ! Multiply fa(g)*fb(g)
    !
    auxg(dfft%nl(:)) = auxg(dfft%nl(:)) * auxr(dfft%nl(:))
    !
    DEALLOCATE( auxr )
    !
    IF ( dfft%lgamma ) auxg(dfft%nlm(:)) = &
         & CMPLX( REAL( auxg(dfft%nl(:)) ), -AIMAG( auxg(dfft%nl(:)) ) ,kind=DP)
    !
    ! Brings convolution back to real space
    !
    CALL invfft('Rho',auxg, dfft)
    !
    fc%of_r(:) = REAL( auxg(:) ) * omega
    !
    DEALLOCATE( auxg )
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE convolution_fft_density
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE convolution_fft_gradient( fft, fa, gb, gc )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    TYPE( fft_core ), TARGET, INTENT(IN) :: fft
    TYPE( environ_density ), INTENT(IN) :: fa
    TYPE( environ_gradient ), INTENT(IN) :: gb
    TYPE( environ_gradient ), INTENT(INOUT) :: gc
    !
    TYPE( environ_density ) :: local
    COMPLEX( DP ), DIMENSION( : ), ALLOCATABLE :: auxr, auxg
    !
    ! ... local aliases
    !
    REAL(DP), POINTER :: omega
    TYPE(fft_type_descriptor), POINTER :: dfft
    !
    INTEGER :: ipol
    !
    ! ... add tests for compatilibity between input, output, and fft
    !
    CALL init_environ_density( fa%cell, local )
    DO ipol = 1, 3
       local%of_r(:) = gb%of_r(ipol,:)
       CALL convolution_fft( fft, fa, local, local)
       gc%of_r(ipol,:) = local%of_r(:)
    ENDDO
    CALL destroy_environ_density( local )
!    dfft => fft % dfft
!    omega => fft % cell % omega
!    !
!    ! Bring fa and fb to reciprocal space
!    !
!    ALLOCATE( auxr( dfft%nnr ) )
!    auxr(:) = CMPLX( fa%of_r(:), 0.D0, kind=DP )
!    CALL fwfft('Rho', auxr, dfft)
!    !
!    ALLOCATE( auxg( dfft%nnr ) )
!    !
!    DO ipol = 1, 3
!       !
!       auxg(:) = CMPLX( gb%of_r(ipol,:), 0.D0, kind=DP )
!       CALL fwfft('Rho', auxg, dfft)
!       !
!       ! Multiply fa(g)*fb(g)
!       !
!       auxg(dfft%nl(:)) = auxg(dfft%nl(:)) * auxr(dfft%nl(:))
!       !
!       IF ( dfft%lgamma ) auxg(dfft%nlm(:)) = &
!         & CMPLX( REAL( auxg(dfft%nl(:)) ), -AIMAG( auxg(dfft%nl(:)) ) ,kind=DP)
!       !
!       ! Brings convolution back to real space
!       !
!       CALL invfft('Rho',auxg, dfft)
!       !
!       gc%of_r(ipol,:) = REAL( auxg(:) ) * omega
!       !
!    END DO
!    !
!    DEALLOCATE( auxr )
!    DEALLOCATE( auxg )
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE convolution_fft_gradient
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE convolution_fft_hessian( fft, fa, hb, hc )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    TYPE( fft_core ), TARGET, INTENT(IN) :: fft
    TYPE( environ_density ), INTENT(IN) :: fa
    TYPE( environ_hessian ), INTENT(IN) :: hb
    TYPE( environ_hessian ), INTENT(INOUT) :: hc
    !
    COMPLEX( DP ), DIMENSION( : ), ALLOCATABLE :: auxr, auxg
    !
    ! ... local aliases
    !
    REAL(DP), POINTER :: omega
    TYPE(fft_type_descriptor), POINTER :: dfft
    TYPE( environ_density ) :: local
    !
    INTEGER :: ipol, jpol
    !
    ! ... add tests for compatilibity between input, output, and fft
    !
    CALL init_environ_density( fa%cell, local )
    DO ipol = 1, 3
       DO jpol = 1, 3
          local%of_r(:) = hb%of_r(ipol,jpol,:)
          CALL convolution_fft( fft, fa, local, local)
          hc%of_r(ipol,jpol,:) = local%of_r(:)
       ENDDO
    ENDDO
    CALL destroy_environ_density( local )
!    dfft => fft % dfft
!    omega => fft % cell % omega
!    !
!    ! Bring fa and fb to reciprocal space
!    !
!    ALLOCATE( auxr( dfft%nnr ) )
!    auxr(:) = CMPLX( fa%of_r(:), 0.D0, kind=DP )
!    CALL fwfft('Rho', auxr, dfft)
!    !
!    ALLOCATE( auxg( dfft%nnr ) )
!    !
!    DO ipol = 1, 3
!       DO jpol = 1, 3
!          !
!          auxg(:) = CMPLX( hb%of_r(ipol,jpol,:), 0.D0, kind=DP )
!          CALL fwfft('Rho', auxg, dfft)
!          !
!          ! Multiply fa(g)*fb(g)
!          !
!          auxg(dfft%nl(:)) = auxg(dfft%nl(:)) * auxr(dfft%nl(:))
!          !
!          IF ( dfft%lgamma ) auxg(dfft%nlm(:)) = &
!               & CMPLX( REAL( auxg(dfft%nl(:)) ), -AIMAG( auxg(dfft%nl(:)) ) ,kind=DP)
!          !
!          ! Brings convolution back to real space
!          !
!          CALL invfft('Rho',auxg, dfft)
!          !
!          hc%of_r(ipol,jpol,:) = REAL( auxg(:) ) * omega
!          !
!       END DO
!    END DO
!    !
!    DEALLOCATE( auxr )
!    DEALLOCATE( auxg )
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE convolution_fft_hessian
!--------------------------------------------------------------------
!----------------------------------------------------------------------------
  SUBROUTINE gradient_fft( fft, a, ga )
!----------------------------------------------------------------------------
    !
    ! ... Calculates ga = \grad a
    ! ... input : fft     FFT descriptor and G vectors
    ! ...         a(:)     a real function on the real-space FFT grid
    ! ... output: ga(3,:)  \grad a, real, on the real-space FFT grid
    !
    !
    IMPLICIT NONE
    !
    TYPE( fft_core ), TARGET, INTENT(IN) :: fft
    TYPE( environ_density ), INTENT(IN)  :: a
    TYPE( environ_gradient ), INTENT(INOUT) :: ga
    !
    INTEGER  :: ipol
    COMPLEX(DP), ALLOCATABLE :: aux(:), gaux(:)
    !
    ! ... local aliases
    !
    REAL(DP), POINTER :: tpiba
    TYPE(fft_type_descriptor), POINTER :: dfft
    REAL(DP), DIMENSION(:,:), POINTER :: g
    !
    ! ... add tests for compatilibity between input, output, and fft
    !
    tpiba => fft % cell % tpiba
    dfft => fft % dfft
    g => fft % g
    !
    ALLOCATE(  aux( dfft%nnr ) )
    ALLOCATE( gaux( dfft%nnr ) )
    !
    aux = CMPLX( a%of_r(:), 0.0_dp, kind=DP)
    !
    ! ... bring a(r) to G-space, a(G) ...
    !
    CALL fwfft ('Rho', aux, dfft)
    !
    ! ... multiply by (iG) to get (\grad_ipol a)(G) ...
    !
    DO ipol = 1, 3
       !
       gaux(:) = (0.0_dp, 0.0_dp)
       !
       gaux(dfft%nl(:)) = g(ipol,:) * CMPLX( -AIMAG( aux(dfft%nl(:)) ), &
            REAL( aux(dfft%nl(:)) ), kind=DP)
       !
       IF ( dfft%lgamma ) THEN
          !
          gaux(dfft%nlm(:)) = CMPLX(  REAL( gaux(dfft%nl(:)) ), &
               -AIMAG( gaux(dfft%nl(:)) ), kind=DP)
          !
       END IF
       !
       ! ... bring back to R-space, (\grad_ipol a)(r) ...
       !
       CALL invfft ('Rho', gaux, dfft)
       !
       ! ...and add the factor 2\pi/a  missing in the definition of G
       !
       ga%of_r(ipol,:) = tpiba * DBLE( gaux(:) )
       !
    END DO
    !
    DEALLOCATE( gaux )
    DEALLOCATE( aux )
    !
    RETURN
    !
  END SUBROUTINE gradient_fft
  !
!----------------------------------------------------------------------------
  SUBROUTINE graddot_fft( fft, ga, da )
!----------------------------------------------------------------------------
    !
    ! ... Calculates da = \sum_i \grad_i a_i in R-space
    ! ... input : fft      FFT descriptor and G vectors
    ! ...         ga(3,:)  a real function on the real-space FFT grid
    ! ... output: da(:)    \sum_i \grad_i a_i, real, on the real-space FFT grid
    !
    IMPLICIT NONE
    !
    TYPE( fft_core ), TARGET, INTENT(IN) :: fft
    TYPE( environ_gradient ), INTENT(IN) :: ga
    TYPE( environ_density ), INTENT(INOUT) :: da
    !
    INTEGER                  :: n, ipol
    COMPLEX(DP), ALLOCATABLE :: aux(:), gaux(:)
    COMPLEX(DP) :: fp, fm, aux1, aux2
    !
    ! ... local aliases
    !
    REAL(DP), POINTER :: tpiba
    TYPE(fft_type_descriptor), POINTER :: dfft
    REAL(DP), DIMENSION(:,:), POINTER :: g
    !
    ! ... add tests for compatilibity between input, output, and fft
    !
    tpiba => fft % cell % tpiba
    dfft => fft % dfft
    g => fft % g
    !
    ALLOCATE( aux(dfft%nnr) )
    ALLOCATE( gaux(dfft%nnr) )
    !
    gaux(:) = (0.0_dp,0.0_dp)
    !
    IF ( dfft%lgamma ) THEN
       !
       ! Gamma tricks: perform 2 FFT's in a single shot
       ! x and y
       ipol = 1
       aux(:) = CMPLX( ga%of_r(ipol,:), ga%of_r(ipol+1,:), kind=DP)
       !
       ! ... bring a(ipol,r) to G-space, a(G) ...
       !
       CALL fwfft ('Rho', aux, dfft)
       !
       ! ... multiply by iG to get the gradient in G-space
       !
       DO n = 1, dfft%ngm
          !
          fp = (aux(dfft%nl(n)) + aux (dfft%nlm(n)))*0.5_dp
          fm = (aux(dfft%nl(n)) - aux (dfft%nlm(n)))*0.5_dp
          aux1 = CMPLX( REAL(fp), AIMAG(fm), kind=DP)
          aux2 = CMPLX(AIMAG(fp), -REAL(fm), kind=DP)
          gaux (dfft%nl(n)) = &
               CMPLX(0.0_dp, g(ipol  ,n),kind=DP) * aux1 + &
               CMPLX(0.0_dp, g(ipol+1,n),kind=DP) * aux2
       ENDDO
       ! z
       ipol = 3
       aux(:) = CMPLX( ga%of_r(ipol,:), 0.0_dp, kind=DP)
       !
       ! ... bring a(ipol,r) to G-space, a(G) ...
       !
       CALL fwfft ('Rho', aux, dfft)
       !
       ! ... multiply by iG to get the gradient in G-space
       ! ... fill both gaux(G) and gaux(-G) = gaux*(G)
       !
       DO n = 1, dfft%ngm
          gaux(dfft%nl(n)) = gaux(dfft%nl(n)) + g(ipol,n) * &
               CMPLX( -AIMAG( aux(dfft%nl(n)) ), &
               REAL( aux(dfft%nl(n)) ), kind=DP)
          gaux(dfft%nlm(n)) = CONJG( gaux(dfft%nl(n)) )
       END DO
       !
    ELSE
       !
       DO ipol = 1, 3
          !
          aux = CMPLX( ga%of_r(ipol,:), 0.0_dp, kind=DP)
          !
          ! ... bring a(ipol,r) to G-space, a(G) ...
          !
          CALL fwfft ('Rho', aux, dfft)
          !
          ! ... multiply by iG to get the gradient in G-space
          !
          DO n = 1, dfft%ngm
             gaux(dfft%nl(n)) = gaux(dfft%nl(n)) + g(ipol,n) * &
                  CMPLX( -AIMAG( aux(dfft%nl(n)) ), &
                  REAL( aux(dfft%nl(n)) ), kind=DP)
          END DO
          !
       END DO
       !
    END IF
    !
    ! ... bring back to R-space, (\grad_ipol a)(r) ...
    !
    CALL invfft ('Rho', gaux, dfft)
    !
    ! ... add the factor 2\pi/a  missing in the definition of G and sum
    !
    da%of_r(:) = tpiba * REAL( gaux(:) )
    !
    DEALLOCATE( aux, gaux )
    !
    RETURN
    !
  END SUBROUTINE graddot_fft

!--------------------------------------------------------------------
! Routines computing laplacian via FFT
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE laplacian_fft( fft, a, lapla )
!--------------------------------------------------------------------
    !
    ! ... Calculates lapla = laplacian(a)
    ! ... input : fft     FFT descriptor and G vectors
    ! ...         a(:)     a real function on the real-space FFT grid
    ! ... output: lapla(:) \nabla^2 a, real, on the real-space FFT grid
    !
    IMPLICIT NONE
    !
    TYPE(fft_core), TARGET, INTENT(IN) :: fft
    TYPE(environ_density), INTENT(IN) :: a
    TYPE(environ_density), INTENT(INOUT) :: lapla
    !
    INTEGER                  :: ig
    COMPLEX(DP), ALLOCATABLE :: aux(:), laux(:)
    !
    ! ... local aliases
    !
    REAL(DP), POINTER :: tpiba2
    REAL(DP), DIMENSION(:), POINTER :: gg
    TYPE(fft_type_descriptor), POINTER :: dfft
    !
    ! ... add tests for compatilibity between input, output, and fft
    !
    tpiba2 => fft % cell % tpiba2
    gg => fft % gg
    dfft => fft % dfft
    !
    ALLOCATE(  aux( dfft%nnr ) )
    ALLOCATE( laux( dfft%nnr ) )
    !
    aux = CMPLX( a%of_r(:), 0.0_dp, kind=DP)
    !
    ! ... bring a(r) to G-space, a(G) ...
    !
    CALL fwfft ('Rho', aux, dfft)
    !
    ! ... Compute the laplacian
    !
    laux(:) = (0.0_dp, 0.0_dp)
    !
    laux(dfft%nl(:)) = -gg(:)*aux(dfft%nl(:))
    !
    IF ( dfft%lgamma ) THEN
       !
       laux(dfft%nlm(:)) = CMPLX( REAL(laux(dfft%nl(:)) ), &
            -AIMAG(laux(dfft%nl(:)) ), kind=DP)
       !
    ENDIF
    !
    ! ... bring back to R-space, (\lapl a)(r) ...
    !
    CALL invfft ('Rho', laux, dfft)
    !
    ! ... add the missing factor (2\pi/a)^2 in G
    !
    lapla%of_r = tpiba2 * REAL( laux )
    !
    DEALLOCATE( laux )
    DEALLOCATE( aux )
    !
    RETURN
    !
  END SUBROUTINE laplacian_fft
  !
!--------------------------------------------------------------------
! Routines computing hessian via FFT
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE hessian_fft( fft, a, ga, ha )
!--------------------------------------------------------------------
    !
    ! ... Calculates ga = \grad a and ha = hessian(a)
    ! ... input : fft     FFT descriptor and G vectors
    ! ...         a(:)     a real function on the real-space FFT grid
    ! ... output: ga(3,:)  \grad a, real, on the real-space FFT grid
    ! ...         ha(3,3,:)  hessian(a), real, on the real-space FFT grid
    !
    IMPLICIT NONE
    !
    TYPE(fft_core), TARGET, INTENT(IN) :: fft
    TYPE(environ_density), INTENT(IN)  :: a
    TYPE(environ_gradient), INTENT(INOUT) :: ga
    TYPE(environ_hessian), INTENT(INOUT) :: ha
    !
    INTEGER                  :: ipol, jpol
    COMPLEX(DP), ALLOCATABLE :: aux(:), gaux(:), haux(:)
    !
    ! ... local aliases
    !
    REAL(DP), POINTER :: tpiba
    REAL(DP), DIMENSION(:,:), POINTER :: g
    TYPE(fft_type_descriptor), POINTER :: dfft
    !
    ! ... add tests for compatilibity between input, output, and fft
    !
    tpiba => fft % cell % tpiba
    g => fft % g
    dfft => fft % dfft
    !
    ALLOCATE(  aux( dfft%nnr ) )
    ALLOCATE( gaux( dfft%nnr ) )
    ALLOCATE( haux( dfft%nnr ) )
    !
    aux = CMPLX( a%of_r(:), 0.0_dp, kind=DP)
    !
    ! ... bring a(r) to G-space, a(G) ...
    !
    CALL fwfft ('Rho', aux, dfft)
    !
    ! ... multiply by (iG) to get (\grad_ipol a)(G) ...
    !
    DO ipol = 1, 3
       !
       gaux(:) = (0.0_dp,0.0_dp)
       !
       gaux(dfft%nl(:)) = g(ipol,:) * CMPLX( -AIMAG( aux(dfft%nl(:)) ), &
            REAL( aux(dfft%nl(:)) ), kind=DP )
       !
       IF ( dfft%lgamma ) THEN
          !
          gaux(dfft%nlm(:)) = CMPLX(  REAL( gaux(dfft%nl(:)) ), &
               -AIMAG( gaux(dfft%nl(:)) ), kind=DP)
          !
       END IF
       !
       ! ... bring back to R-space, (\grad_ipol a)(r) ...
       !
       CALL invfft ('Rho', gaux, dfft)
       !
       ! ...and add the factor 2\pi/a  missing in the definition of G
       !
       ga%of_r(ipol,:) = tpiba * REAL( gaux(:) )
       !
       ! ... compute the second derivatives
       !
       DO jpol = 1, ipol
          !
          haux(:) = (0.0_dp,0.0_dp)
          !
          haux(dfft%nl(:)) = - g(ipol,:) * g(jpol,:) * &
               CMPLX( REAL( aux(dfft%nl(:)) ), &
               AIMAG( aux(dfft%nl(:)) ), kind=DP)
          !
          IF ( dfft%lgamma ) THEN
             !
             haux(dfft%nlm(:)) = CMPLX(  REAL( haux(dfft%nl(:)) ), &
                  -AIMAG( haux(dfft%nl(:)) ), kind=DP)
             !
          END IF
          !
          ! ... bring back to R-space, (\grad_ipol a)(r) ...
          !
          CALL invfft ('Rho', haux, dfft)
          !
          ! ...and add the factor 2\pi/a  missing in the definition of G
          !
          ha%of_r(ipol, jpol, :) = tpiba * tpiba * REAL( haux(:) )
          !
          ha%of_r(jpol, ipol, :) = ha%of_r(ipol, jpol, :)
          !
       END DO
       !
    END DO
    !
    DEALLOCATE( haux )
    DEALLOCATE( gaux )
    DEALLOCATE( aux )
    !
    RETURN
    !
  END SUBROUTINE hessian_fft
!--------------------------------------------------------------------
SUBROUTINE hessv_h_of_rho_r( rho, hessv, fft )
   !-----------------------------------------------------------------
   !
   ! ... Gradient of Hartree potential in R space from a total 
   !     (spinless) density in R space n(r)
   ! wg_corr_h => calc_vmt
   USE correction_mt, ONLY : calc_vmt
   !
   IMPLICIT NONE
   !
   ! ... Declares variables
   !
   TYPE( fft_core ), INTENT(IN), TARGET :: fft
   REAL( DP ), INTENT(IN)     :: rho( fft % dfft % nnr )
   REAL( DP ), INTENT(OUT)    :: hessv( 3, 3, fft % dfft % nnr )
   !
   ! ... cell variables
   !
   REAL( DP ), POINTER :: omega, tpiba2
   !
   ! ... fft variables
   !
   INTEGER, POINTER    :: ngm, gstart
   LOGICAL, POINTER    :: do_comp_mt
   REAL( DP ), POINTER :: g(:,:), gg(:)
   !
   ! ... Local variables
   !
   COMPLEX( DP ), ALLOCATABLE :: rhoaux( : )
   COMPLEX( DP ), ALLOCATABLE :: gaux( : )
   COMPLEX( DP ), ALLOCATABLE :: rgtot(:), vaux(:)
   REAL( DP )                 :: fac, eh_corr
   INTEGER                    :: ig, ipol, jpol
   LOGICAL                    :: gamma_only
   !
   ! ... Setting pointer variables
   !
   omega => fft % cell % omega
   tpiba2 => fft % cell % tpiba2
   ngm => fft % ngm
   gstart => fft % gstart
   gg => fft % gg
   g => fft % g
   do_comp_mt => fft % use_internal_pbc_corr
   gamma_only = .TRUE.
   !
   ! ... Bring rho to G space
   !
   ALLOCATE( rhoaux( fft % dfft % nnr ) )
   rhoaux( : ) = CMPLX( rho( : ), 0.D0, KIND=dp )
   !
   CALL fwfft('Rho', rhoaux, fft%dfft)
   !
   ! ... Compute total potential in G space
   !
   ALLOCATE( gaux( fft % dfft % nnr ) )
   !
   DO ipol = 1, 3
      !
      DO jpol = 1, 3
         !
         gaux(:) = (0.0_dp,0.0_dp)
         !
         DO ig = gstart, ngm
            !
            fac = g(ipol,ig) * g(jpol,ig) / gg(ig)
            gaux(fft%dfft%nl(ig)) = CMPLX(REAL(rhoaux(fft%dfft%nl(ig))),AIMAG(rhoaux(fft%dfft%nl(ig))),kind=dp) * fac 
            !
         END DO
         !
         ! ...and add the factor e2*fpi coming from the missing prefactor of 
         !  V = e2 * fpi
         !
         fac = e2 * fpi
         gaux = gaux * fac 
         !
         ! ...add martyna-tuckerman correction, if needed
         ! 
         if (do_comp_mt) then
            ALLOCATE( vaux( ngm ), rgtot(ngm) )
            rgtot(1:ngm) = rhoaux(fft%dfft%nl(1:ngm))
            !CALL wg_corr_h (omega, ngm, rgtot, vaux, eh_corr)
            CALL calc_vmt( fft, rgtot, vaux )
            DO ig = gstart, ngm
               fac = g(ipol,ig) * g(jpol,ig) * tpiba2
               gaux(fft%dfft%nl(ig)) = gaux(fft%dfft%nl(ig)) + CMPLX(REAL(vaux(ig)),AIMAG(vaux(ig)),kind=dp)*fac 
            END DO
            DEALLOCATE( rgtot, vaux )
         end if
         !
         IF ( gamma_only ) THEN
            !
            gaux(fft%dfft%nlm(:)) = &
                 CMPLX( REAL( gaux(fft%dfft%nl(:)) ), -AIMAG( gaux(fft%dfft%nl(:)) ) ,kind=DP)
            !
         END IF
         !
         ! ... bring back to R-space
         !
         CALL invfft ('Rho', gaux, fft%dfft)
         !
         hessv(ipol,jpol,:) = REAL( gaux(:) )
         !
      ENDDO
      !
   ENDDO
   !
   DEALLOCATE(gaux)
   !
   DEALLOCATE(rhoaux)
   !
   RETURN
   !
 END SUBROUTINE hessv_h_of_rho_r
 !----------------------------------------------------------------------------
 SUBROUTINE field_of_gradrho( gradrho, e, fft )
   !----------------------------------------------------------------------------
   !
   ! ... Gradient of Hartree potential in R space from a total
   !     (spinless) density in R space n(r)
   !
   USE correction_mt, ONLY : calc_vmt
   !
   IMPLICIT NONE
   !
   ! ... Declares variables
   !
   TYPE( fft_core ), INTENT(IN), TARGET :: fft
   REAL( DP ), INTENT(IN)           :: gradrho( 3, fft % dfft % nnr )
   REAL( DP ), INTENT(OUT)          :: e( fft % dfft % nnr )
   !
   ! ... cell variables
   !
   REAL( DP ), POINTER :: omega, tpiba
   !
   ! ... fft variables
   !
   INTEGER, POINTER    :: ngm, gstart
   LOGICAL, POINTER    :: do_comp_mt
   REAL( DP ), POINTER :: g(:,:), gg(:)
   !
   ! ... Local variables
   !
   COMPLEX( DP ), ALLOCATABLE :: aux( : )
   COMPLEX( DP ), ALLOCATABLE :: eaux( : )
   COMPLEX( DP ), ALLOCATABLE :: gaux( : )
   COMPLEX( DP ), ALLOCATABLE :: rgtot(:), vaux(:)
   REAL( DP )                 :: fac, eh_corr
   INTEGER                    :: ig, ipol
   LOGICAL                    :: gamma_only
   !
   ! ... Setting pointer variables
   !
   omega => fft % cell % omega
   tpiba => fft % cell % tpiba
   ngm => fft % ngm
   gstart => fft % gstart
   gg => fft % gg
   g => fft % g
   do_comp_mt => fft % use_internal_pbc_corr
   gamma_only = .TRUE.
   !
   ! ... Bring gradrho to G space
   !
   ALLOCATE( eaux( fft % dfft % nnr ) )
   eaux( : ) = CMPLX( 0.D0, 0.D0, KIND=dp )
   !
   ALLOCATE( aux( fft % dfft % nnr ) )
   aux( : ) = CMPLX( 0.D0, 0.D0, KIND=dp )
   !
   ALLOCATE( gaux( fft % dfft % nnr ) )
   !
   IF ( do_comp_mt ) ALLOCATE( vaux( ngm ), rgtot( ngm ) )
   !
   DO ipol = 1, 3
      !
      gaux( : ) = CMPLX( gradrho( ipol, : ), 0.D0, KIND=dp )
      !
      CALL fwfft('Rho', gaux, fft % dfft)
      !
      ! ... Compute total potential in G space
      !
      DO ig = gstart, ngm
         !
         fac = g(ipol,ig) / gg(ig)
         aux(fft%dfft%nl(ig)) = CMPLX(-AIMAG(gaux(fft%dfft%nl(ig))),REAL(gaux(fft%dfft%nl(ig))),kind=dp) * fac
         !
      END DO
      !
      ! ...and add the factor e2*fpi/2\pi/a coming from the missing prefactor of
      !  V = e2 * fpi divided by the 2\pi/a factor missing in G
      !
      fac = e2 * fpi / tpiba
      aux = aux * fac
      !
      ! ...add martyna-tuckerman correction, if needed
      !
      IF ( do_comp_mt ) THEN
         !
         rgtot(1:ngm) = gaux(fft%dfft%nl(1:ngm))
         !
         !CALL wg_corr_h (omega, ngm, rgtot, vaux, eh_corr)
         CALL calc_vmt( fft, rgtot, vaux )
         !
         DO ig = gstart, ngm
            !
            fac = g(ipol,ig) * tpiba
            aux(fft%dfft%nl(ig)) = aux(fft%dfft%nl(ig)) + CMPLX(-AIMAG(vaux(ig)),REAL(vaux(ig)),kind=dp)*fac
            !
         END DO
         !
      END IF
      !
      eaux = eaux + aux
      !
   END DO
   !
   IF ( do_comp_mt ) DEALLOCATE( rgtot, vaux )
   !
   DEALLOCATE(gaux)
   !
   DEALLOCATE(aux)
   !
   IF ( gamma_only ) THEN
      !
      eaux(fft%dfft%nlm(:)) = &
           CMPLX( REAL( eaux(fft%dfft%nl(:)) ), -AIMAG( eaux(fft%dfft%nl(:)) ) ,kind=DP)
      !
   END IF
   !
   ! ... bring back to R-space, (\grad_ipol a)(r) ...
   !
   CALL invfft ('Rho', eaux, fft%dfft)
   !
   e(:) = REAL( eaux(:) )
   !
   DEALLOCATE(eaux)
   !
   RETURN
   !
 END SUBROUTINE field_of_gradrho
 
END MODULE core_fft
