MODULE core_fft
  !
  USE environ_types
  USE electrostatic_types
  !
CONTAINS
!--------------------------------------------------------------------
  SUBROUTINE poisson_fft( fft_core, fin, fout )
!--------------------------------------------------------------------
    !
    ! Solves the Poisson equation \nabla \fout(r) = - 4 \pi \fin(r)
    ! Input and output functions are defined in real space
    !
    USE libs_fft_interfaces, ONLY : libs_fwfft, libs_invfft
    USE correction_mt, ONLY : calc_vmt
    !
    IMPLICIT NONE
    !
    TYPE( fft_core_type ), INTENT(IN) :: fft_core
    !
    TYPE( environ_density ), INTENT(IN) :: fin
    !
    TYPE( environ_density ), INTENT(OUT) :: fout
    !
    ! Test compatilibity of fin, fout, fft_core
    !
    !
    ! Bring fin%of_r from R-space --> G-space
    !
    ALLOCATE( auxr( dfftp%nnr ) )
    auxr = CMPLX( fin%of_r, 0.D0, kind=DP )
    CALL libs_fwfft( 'Rho', auxr, dfftp )
    ALLOCATE( auxg( dfftp%ngm ) )
    auxg = auxr(dfftp%nl(:))
    !
    auxr = CMPLX( 0.D0, 0.D0, kind=DP )
!$omp parallel do
    DO ig = gstart, ngm
       auxr(dfftp%nl(ig)) = auxg(ig) / gg(ig)
    ENDDO
!$omp end parallel do
    auxr = auxr * e2 * fpi / tpiba2
    !
    IF ( fft_core%do_comp_mt ) THEN
       ALLOCATE( vaux( ngm ) )
       CALL calc_vmt(omega, ngm, auxg, vaux, eh_corr)
       auxr(dfftp%nl(1:ngm)) = auxr(dfftp%np(1:ngm)) + vaux(1:ngm)
       DEALLOCATE( vaux )
    END IF
    !
    IF ( fft_core % gamma_only ) THEN
       auxr( dfftp%nlm(:) ) = CMPLX( REAL(auxr(dfftp%nl(:)), -AIMAG(auxr(dfftp%nl(:)), kind=DP )
    END IF
    !
    ! ... transform hartree potential to real space
    !
    CALL libs_invfft ('Rho', auxr, dfftp)
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
  SUBROUTINE gradpoisson_fft( fft_core, fin, gout )
!--------------------------------------------------------------------
    !
    ! Solves the Poisson equation \grad \gout(r) = - 4 \pi \fin(r)
    ! where gout is the gradient of the potential
    ! Input and output functions are defined in real space
    !
    USE libs_fft_interfaces, ONLY : libs_fwfft, libs_invfft
    !
    IMPLICIT NONE
    !
    TYPE( fft_core_type ), INTENT(IN) :: fft_core
    !
    TYPE( environ_density ), INTENT(IN) :: fin
    !
    TYPE( environ_gradient ), INTENT(OUT) :: gout
    !
    ! Test compatilibity of fin, gout, fft_core
    !
    !
    ! ... Bring rho to G space
    !
    ALLOCATE( aux( dfftp%nnr ) )
    aux( : ) = CMPLX( fin%of_r( : ), 0.D0, KIND=dp )
    !
    CALL libs_fwfft('Rho', aux, dfftp)
    !
    ! ... Compute total potential in G space
    !
    ALLOCATE( gaux( dfftp%nnr ) )
    !
    DO ipol = 1, 3
       !
       gaux(:) = (0.0_dp,0.0_dp)
       !
!$omp parallel do private(fac)
       DO ig = gstart, ngm
          !
          fac = g(ipol,ig) / gg(ig)
          gaux(dfftp%nl(ig)) = CMPLX(-AIMAG(rhoaux(dfftp%nl(ig))),REAL(rhoaux(dfftp%nl(ig))),kind=dp) * fac
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
       IF ( fft_core%do_comp_mt ) THEN
          ALLOCATE( vaux( ngm ), rgtot(ngm) )
          rgtot(1:ngm) = rhoaux(dfftp%nl(1:ngm))
          CALL calc_vmt(omega, ngm, rgtot, vaux, eh_corr)
!$omp parallel do private(fac)
          DO ig = gstart, ngm
             fac = g(ipol,ig) * tpiba
             gaux(dfftp%nl(ig)) = gaux(dfftp%nl(ig)) + CMPLX(-AIMAG(vaux(ig)),REAL(vaux(ig)),kind=dp)*fac
          END DO
!$omp end parallel do
          DEALLOCATE( rgtot, vaux )
       END IF
       !
       ! Assuming GAMMA ONLY
       !
       IF ( fft_core%gamma_only ) THEN
          gaux(dfftp%nlm(:)) = &
               CMPLX( REAL( gaux(dfftp%nl(:)) ), -AIMAG( gaux(dfftp%nl(:)) ) ,kind=DP)
       END IF
       !
       ! ... bring back to R-space, (\grad_ipol a)(r) ...
       !
       CALL libs_invfft ('Rho', gaux, dfftp)
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
  SUBROUTINE convolution_fft( fft_core, fa, fb, fc )
!--------------------------------------------------------------------
    !
    USE libs_fft_interfaces, ONLY : libs_fwfft, libs_invfft
    !
    IMPLICIT NONE
    !
    TYPE( fft_core_type ), INTENT(IN) :: fft_core
    TYPE( environ_density ), INTENT(IN) :: fa, fb
    TYPE( environ_density ), INTENT(OUT) :: fc
    !
    COMPLEX( DP ), DIMENSION( : ), ALLOCATABLE :: auxr, auxg
    !
    ! Test compatilibity of fin, gout, fft_core
    !
    !
    ! Bring fa and fb to reciprocal space
    !
    ALLOCATE( auxr( nnr ) )
    auxr(:) = CMPLX( fa%of_r(:), 0.D0, kind=DP )
    CALL env_fwfft('Rho', auxr, dfftp)
    !
    ALLOCATE( auxg( nnr ) )
    auxg = 0.D0
    !
    auxg(dfftp%nl(1:ngm)) = auxr(dfftp%nl(1:ngm))
    !
    auxr(:) = CMPLX( fb%of_r(:), 0.D0, kind=DP )
    CALL libs_fwfft('Rho', auxr, dfftp)
    !
    ! Multiply fa(g)*fb(g)
    !
    auxg(dfftp%nl(1:ngm)) = auxg(dfftp%nl(1:ngm)) * auxr(dfftp%nl(1:ngm))
    !
    DEALLOCATE( auxr )
    !
    IF ( fft_core%gamma_only ) auxg(dfftp%nlm(1:ngm)) = &
         & CMPLX( REAL( auxg(dfftp%nl(1:ngm)) ), -AIMAG( auxg(dfftp%nl(1:ngm)) ) ,kind=DP)
    !
    ! Brings convolution back to real space
    !
    CALL libs_invfft('Rho',auxg, dfftp)
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
END MODULE core_fft
