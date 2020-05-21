MODULE core_fft
  !
  USE environ_types
  USE electrostatic_types
  USE fft_interfaces, ONLY : fwfft, invfft
  !
  PRIVATE
  !
  PUBLIC :: poisson_fft, gradpoisson_fft, convolution_fft, gradient_fft, graddot_fft, laplacian_fft, hessian_fft
  !
CONTAINS
!--------------------------------------------------------------------
  SUBROUTINE poisson_fft( fft, fin, fout )
!--------------------------------------------------------------------
    !
    ! Solves the Poisson equation \nabla \fout(r) = - 4 \pi \fin(r)
    ! Input and output functions are defined in real space
    !
!    USE correction_mt, ONLY : calc_vmt
    !
    IMPLICIT NONE
    !
    TYPE( fft_core ), INTENT(IN) :: fft
    TYPE( environ_density ), INTENT(IN) :: fin
    TYPE( environ_density ), INTENT(OUT) :: fout
    !
    ! ... local aliases
    !
    REAL(DP), POINTER :: tpiba2, omega, gstart
    REAL(DP), DIMENSION(:), POINTER :: gg
    TYPE(fft_dlay_descriptor), POINTER :: dfft
    !
    ! ... add tests for compatilibity between input, output, and fft_core
    !
    tpiba => fft % tpiba
    omega => fft % omega
    gg => fft % gg
    dfft => fft % dfft
    !
    ! Bring fin%of_r from R-space --> G-space
    !
    ALLOCATE( auxr( dfft%nnr ) )
    auxr = CMPLX( fin%of_r, 0.D0, kind=DP )
    CALL fwfft( 'Rho', auxr, dfft )
    ALLOCATE( auxg( dfft%ngm ) )
    auxg = auxr(dfft%nl(:))
    !
    auxr = CMPLX( 0.D0, 0.D0, kind=DP )
    auxr(dfft%nl(:)) = auxg(:)/gg(:)
!$omp parallel do
!    DO ig = gstart, ngm
!       auxr(dfft%nl(ig)) = auxg(ig) / gg(ig)
!    ENDDO
!$omp end parallel do
    auxr = auxr * e2 * fpi / tpiba2
    !
!    IF ( fft%do_comp_mt ) THEN
!       ALLOCATE( vaux( ngm ) )
!       CALL calc_vmt(omega, ngm, auxg, vaux, eh_corr)
!       auxr(dfft%nl(1:ngm)) = auxr(dfft%np(1:ngm)) + vaux(1:ngm)
!       DEALLOCATE( vaux )
!    END IF
    !
    IF ( dfft%lgamma ) THEN
       auxr( dfft%nlm(:) ) = CMPLX( REAL(auxr(dfft%nl(:)), -AIMAG(auxr(dfft%nl(:)), kind=DP )
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
    IMPLICIT NONE
    !
    TYPE( fft_core ), INTENT(IN) :: fft
    TYPE( environ_density ), INTENT(IN) :: fin
    TYPE( environ_gradient ), INTENT(OUT) :: gout
    !
    ! ... local aliases
    !
    REAL(DP), POINTER :: tpiba, omega
    REAL(DP), DIMENSION(:), POINTER :: gg
    REAL(DP), DIMENSION(:,:), POINTER :: g
    TYPE(fft_dlay_descriptor), POINTER :: dfft
    !
    ! ... add tests for compatilibity between input, output, and fft
    !
    tpiba => fft % tpiba
    omega => fft % omega
    gg => fft % gg
    g => fft % g
    dfft => fft % dfft
    !
    ! ... Bring rho to G space
    !
    ALLOCATE( aux( dfft%nnr ) )
    aux( : ) = CMPLX( fin%of_r( : ), 0.D0, KIND=dp )
    !
    CALL fwfft('Rho', aux, dfft)
    !
    ! ... Compute total potential in G space
    !
    ALLOCATE( gaux( dfft%nnr ) )
    !
    DO ipol = 1, 3
       !
       gaux(:) = (0.0_dp,0.0_dp)
       !
!$omp parallel do private(fac)
       DO ig = gstart, ngm
          !
          fac = g(ipol,ig) / gg(ig)
          gaux(dfft%nl(ig)) = CMPLX(-AIMAG(rhoaux(dfft%nl(ig))),REAL(rhoaux(dfft%nl(ig))),kind=dp) * fac
          !
       END DO
!$omp end parallel do
       !
       ! ...and add the factor e2*fpi/2\pi/a coming from the missing prefactor of
       !  V = e2 * fpi divided by the 2\pi/a factor missing in G
       !
       fac = e2 * fpi / tpiba
       gaux = gaux * fac
       !
       ! ...add martyna-tuckerman correction, if needed
       !
!       IF ( fft%do_comp_mt ) THEN
!          ALLOCATE( vaux( ngm ), rgtot(ngm) )
!          rgtot(1:ngm) = rhoaux(dfft%nl(1:ngm))
!          CALL calc_vmt(omega, ngm, rgtot, vaux, eh_corr)
!$omp parallel do private(fac)
!          DO ig = gstart, ngm
!             fac = g(ipol,ig) * tpiba
!             gaux(dfft%nl(ig)) = gaux(dfft%nl(ig)) + CMPLX(-AIMAG(vaux(ig)),REAL(vaux(ig)),kind=dp)*fac
!          END DO
!$omp end parallel do
!          DEALLOCATE( rgtot, vaux )
!       END IF
       !
       ! Assuming GAMMA ONLY
       !
       IF ( dfft%lgamma ) THEN
          gaux(dfft%nlm(:)) = &
               CMPLX( REAL( gaux(dfft%nl(:)) ), -AIMAG( gaux(dfft%nl(:)) ) ,kind=DP)
       END IF
       !
       ! ... bring back to R-space, (\grad_ipol a)(r) ...
       !
       CALL invfft ('Rho', gaux, dfft)
       !
       gout%of_r(ipol,:) = REAL( gaux(:) )
       !
    ENDDO
    !
    DEALLOCATE(gaux)
    !
    DEALLOCATE(aux)
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE gradpoisson_fft
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE convolution_fft( fft, fa, fb, fc )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    TYPE( fft_core ), INTENT(IN) :: fft
    TYPE( environ_density ), INTENT(IN) :: fa, fb
    TYPE( environ_density ), INTENT(OUT) :: fc
    !
    COMPLEX( DP ), DIMENSION( : ), ALLOCATABLE :: auxr, auxg
    !
    ! ... local aliases
    !
    REAL(DP), POINTER :: omega
    REAL(DP), DIMENSION(:,:), POINTER :: g
    TYPE(fft_dlay_descriptor), POINTER :: dfft
    !
    ! ... add tests for compatilibity between input, output, and fft
    !
    dfft => fft % dfft
    omega => fft % omega
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
  END SUBROUTINE convolution_fft
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
    TYPE( fft_core ), INTENT(IN) :: fft
    TYPE( environ_density ), INTENT(IN)  :: a
    TYPE( environ_gradient ), INTENT(OUT) :: ga
    !
    INTEGER  :: ipol
    COMPLEX(DP), ALLOCATABLE :: aux(:), gaux(:)
    !
    ! ... local aliases
    !
    REAL(DP), POINTER :: tpiba
    TYPE(fft_dlay_descriptor), POINTER :: dfft
    REAL(DP), DIMENSION(:,:), POINTER :: g
    !
    ! ... add tests for compatilibity between input, output, and fft
    !
    tpiba => fft % tpiba
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
  SUBROUTINE graddot_fft( fft, a, da )
!----------------------------------------------------------------------------
    !
    ! ... Calculates da = \sum_i \grad_i a_i in R-space
    ! ... input : fft     FFT descriptor and G vectors
    ! ...         a(3,:)   a real function on the real-space FFT grid
    ! ... output: ga(:)    \sum_i \grad_i a_i, real, on the real-space FFT grid
    !
    IMPLICIT NONE
    !
    TYPE( fft_core ), INTENT(IN) :: fft
    TYPE( environ_density ), INTENT(IN) :: a
    TYPE( environ_density ), INTENT(OUT) :: da
    !
    INTEGER                  :: n, ipol
    COMPLEX(DP), ALLOCATABLE :: aux(:), gaux(:)
    COMPLEX(DP) :: fp, fm, aux1, aux2
    !
    ! ... local aliases
    !
    REAL(DP), POINTER :: tpiba
    TYPE(fft_dlay_descriptor), POINTER :: dfft
    REAL(DP), DIMENSION(:,:), POINTER :: g
    !
    ! ... add tests for compatilibity between input, output, and fft
    !
    tpiba => fft % tpiba
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
       aux(:) = CMPLX( a%of_r(ipol,:), a%of_r(ipol+1,:), kind=DP)
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
       aux(:) = CMPLX( a%of_r(ipol,:), 0.0_dp, kind=DP)
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
          aux = CMPLX( a%of_r(ipol,:), 0.0_dp, kind=DP)
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
    TYPE(fft_core),INTENT(IN) :: fft
    TYPE(environ_density), INTENT(IN) :: a
    TYPE(environ_density), INTENT(OUT) :: lapla(dfft%nnr)
    !
    INTEGER                  :: ig
    COMPLEX(DP), ALLOCATABLE :: aux(:), laux(:)
    !
    ! ... local aliases
    !
    REAL(DP), POINTER :: tpiba2
    REAL(DP), DIMENSION(:), POINTER :: gg
    TYPE(fft_dlay_descriptor), POINTER :: dfft
    !
    ! ... add tests for compatilibity between input, output, and fft
    !
    tpiba2 => fft % tpiba2
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
    TYPE(fft_core),INTENT(IN) :: fft
    TYPE(environ_density), INTENT(IN)  :: a
    TYPE(environ_gradient), INTENT(OUT) :: ga
    TYPE(environ_hessian), INTENT(OUT) :: ha
    !
    INTEGER                  :: ipol, jpol
    COMPLEX(DP), ALLOCATABLE :: aux(:), gaux(:), haux(:)
    !
    ! ... local aliases
    !
    REAL(DP), POINTER :: tpiba
    REAL(DP), DIMENSION(:,:), POINTER :: g
    TYPE(fft_dlay_descriptor), POINTER :: dfft
    !
    ! ... add tests for compatilibity between input, output, and fft
    !
    tpiba => fft % tpiba
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
END MODULE core_fft
