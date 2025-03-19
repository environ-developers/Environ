MODULE class_core_fft_lowpass
    !------------------------------------------------------------------------------------
    !
    USE class_core_fft
    !
    USE class_io, ONLY: io
    !
    USE env_fft_interfaces, ONLY: env_fwfft, env_invfft
    !
    USE environ_param, ONLY: DP, pi, tpi, tpi2, pi_third
    !
    USE class_density
    USE class_gradient
    USE class_hessian
    !
    USE tools_math, ONLY: environ_erfc
    !
    !------------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    PRIVATE
    !
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    TYPE, EXTENDS(core_fft), PUBLIC :: core_fft_lowpass
        !--------------------------------------------------------------------------------
        !
        REAL(DP) :: filter_shift, filter_fac
        !
        !--------------------------------------------------------------------------------
    CONTAINS
        !--------------------------------------------------------------------------------
        !
        PROCEDURE :: gradient => gradient_fft_lowpass
        PROCEDURE :: divergence => divergence_fft_lowpass
        PROCEDURE :: laplacian => laplacian_fft_lowpass
        PROCEDURE :: hessian => hessian_fft_lowpass
        !
        PROCEDURE, PRIVATE :: filter_lowpass
        !
        !--------------------------------------------------------------------------------
    END TYPE core_fft_lowpass
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
CONTAINS
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !                                 DERIVATIVE METHODS
    !
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !>
    !! Calculates grad = \grad a
    !! input : f(:)      a real function on the real-space FFT grid
    !! output: grad(3,:) \grad a, real, on the real-space FFT grid
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE gradient_fft_lowpass(this, f, grad)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(core_fft_lowpass), INTENT(IN) :: this
        TYPE(environ_density), INTENT(IN) :: f
        !
        TYPE(environ_gradient), INTENT(INOUT) :: grad
        !
        INTEGER :: i
        COMPLEX(DP), DIMENSION(:), ALLOCATABLE :: aux, gaux
        !
        !--------------------------------------------------------------------------------
        !
        ASSOCIATE (dfft => this%cell%dfft, &
                   nnr => this%cell%dfft%nnr, &
                   nl => this%cell%dfft%nl)
            !
            !----------------------------------------------------------------------------
            !
            ALLOCATE (aux(nnr))
            ALLOCATE (gaux(nnr))
            !
            aux = CMPLX(f%of_r, kind=DP)
            !
            CALL env_fwfft('Rho', aux, dfft) ! bring a(r) to G-space, a(G)
            !
            !----------------------------------------------------------------------------
            ! Multiply by (iG) to get (\grad_ipol a)(G)
            !
            DO i = 1, 3
                gaux = (0.0_DP, 0.0_DP)
                !
                gaux(nl) = this%cell%g(i, :) * &
                           this%filter_lowpass(this%cell%dfft%ngm, this%cell%gg(:)) * &
                           CMPLX(-AIMAG(aux(nl)), REAL(aux(nl)), kind=DP)
                !
                IF (dfft%lgamma) &
                    gaux(dfft%nlm) = CMPLX(REAL(gaux(nl)), -AIMAG(gaux(nl)), kind=DP)
                !
                CALL env_invfft('Rho', gaux, dfft) ! bring back to R-space, (\grad_ipol a)(r)
                !
                grad%of_r(i, :) = tpi * DBLE(gaux)
                ! add the factor 2\pi/a missing in the definition of G
                !
            END DO
            !
            DEALLOCATE(aux)
            DEALLOCATE(gaux)
            !
        END ASSOCIATE
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE gradient_fft_lowpass
    !------------------------------------------------------------------------------------
    !>
    !! Calculates div = \sum_i \grad_i a_i in R-space
    !! input : grad(3,:) a real function on the real-space FFT grid
    !! output: div(:)    \sum_i \grad_i a_i, real, on the real-space FFT grid
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE divergence_fft_lowpass(this, grad, div)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(core_fft_lowpass), INTENT(IN) :: this
        TYPE(environ_gradient), INTENT(IN) :: grad
        !
        TYPE(environ_density), INTENT(INOUT) :: div
        !
        INTEGER :: n, i
        COMPLEX(DP), DIMENSION(:), ALLOCATABLE :: aux, gaux
        COMPLEX(DP) :: fp, fm, aux1, aux2
        !
        !--------------------------------------------------------------------------------
        !
        ASSOCIATE (dfft => this%cell%dfft, &
                   nnr => this%cell%dfft%nnr, &
                   nl => this%cell%dfft%nl, &
                   nlm => this%cell%dfft%nlm, &
                   ngm => this%cell%dfft%ngm, &
                   g => this%cell%g)
            !
            !----------------------------------------------------------------------------
            !
            ALLOCATE (aux(nnr))
            ALLOCATE (gaux(nnr))
            !
            gaux = (0.0_DP, 0.0_DP)
            !
            IF (dfft%lgamma) THEN
                !
                !------------------------------------------------------------------------
                ! Gamma tricks: perform 2 FFT's in a single shot
                !
                ! x and y
                !
                i = 1
                aux = CMPLX(grad%of_r(i, :), grad%of_r(i + 1, :), kind=DP)
                !
                CALL env_fwfft('Rho', aux, dfft) ! bring a(i,r) to G-space, a(G)
                !
                !------------------------------------------------------------------------
                ! Multiply by iG to get the gradient in G-space
                !
                DO n = 1, ngm
                    fp = (aux(nl(n)) + aux(nlm(n))) * 0.5_DP
                    fm = (aux(nl(n)) - aux(nlm(n))) * 0.5_DP
                    aux1 = CMPLX(REAL(fp), AIMAG(fm), kind=DP)
                    aux2 = CMPLX(AIMAG(fp), -REAL(fm), kind=DP)
                    !
                    gaux(nl(n)) = CMPLX(0.0_DP, g(i, n), kind=DP) * aux1 + &
                                  CMPLX(0.0_DP, g(i + 1, n), kind=DP) * aux2
                    !
                END DO
                !
                !------------------------------------------------------------------------
                ! z
                !
                i = 3
                aux = CMPLX(grad%of_r(i, :), kind=DP)
                !
                CALL env_fwfft('Rho', aux, dfft) ! bring a(i, r) to G-space, a(G)
                !
                !------------------------------------------------------------------------
                ! Multiply by iG to get the gradient in G-space
                ! Fill both gaux(G) and gaux(-G) = gaux*(G)
                !
                DO n = 1, ngm
                    !
                    gaux(nl(n)) = gaux(nl(n)) + &
                                  g(i, n) * &
                                  CMPLX(-AIMAG(aux(nl(n))), REAL(aux(nl(n))), kind=DP)
                    !
                    gaux(nlm(n)) = CONJG(gaux(nl(n)))
                END DO
                !
            ELSE
                !
                DO i = 1, 3
                    aux = CMPLX(grad%of_r(i, :), kind=DP)
                    !
                    CALL env_fwfft('Rho', aux, dfft) ! bring a(i,r) to G-space, a(G)
                    !
                    !--------------------------------------------------------------------
                    ! Multiply by iG to get the gradient in G-space
                    !
                    DO n = 1, ngm
                        !
                        gaux(nl(n)) = &
                            gaux(nl(n)) + &
                            g(i, n) * &
                            CMPLX(-AIMAG(aux(nl(n))), REAL(aux(nl(n))), kind=DP)
                        !
                    END DO
                    !
                END DO
                !
            END IF
            !
            CALL env_invfft('Rho', gaux, dfft) ! bring back to R-space, (\grad_ipol a)(r)
            !
            div%of_r = tpi * REAL(gaux)
            ! Add the factor 2\pi/a missing in the definition of G and sum
            !
        END ASSOCIATE
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE divergence_fft_lowpass
    !------------------------------------------------------------------------------------
    !>
    !! Calculates lapla = laplacian(a)
    !! input : f(:)       A real function on the real-space FFT grid
    !! output: lapla(:)   \nabla^2 a, real, on the real-space FFT grid
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE laplacian_fft_lowpass(this, f, lapla)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(core_fft_lowpass), INTENT(IN) :: this
        TYPE(environ_density), INTENT(IN) :: f
        !
        TYPE(environ_density), INTENT(INOUT) :: lapla
        !
        COMPLEX(DP), DIMENSION(:), ALLOCATABLE :: aux, laux
        !
        !--------------------------------------------------------------------------------
        !
        ASSOCIATE (dfft => this%cell%dfft, &
                   nnr => this%cell%dfft%nnr, &
                   nl => this%cell%dfft%nl)
            !
            !----------------------------------------------------------------------------
            !
            ALLOCATE (aux(nnr))
            ALLOCATE (laux(nnr))
            !
            aux = CMPLX(f%of_r, kind=DP)
            !
            CALL env_fwfft('Rho', aux, dfft) ! bring a(r) to G-space, a(G)
            !
            !----------------------------------------------------------------------------
            ! Compute the laplacian
            !
            laux = (0.0_DP, 0.0_DP)
            !
            laux(nl) = -this%cell%gg * aux(nl) * &
                        this%filter_lowpass(this%cell%dfft%ngm, this%cell%gg(:))
            !
            IF (dfft%lgamma) &
                laux(dfft%nlm) = CMPLX(REAL(laux(nl)), -AIMAG(laux(nl)), kind=DP)
            !
            CALL env_invfft('Rho', laux, dfft) ! bring back to R-space, (\lapl a)(r)
            !
            lapla%of_r = tpi2 * REAL(laux) ! add the missing factor (2\pi)^2 in G
            !
        END ASSOCIATE
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE laplacian_fft_lowpass
    !------------------------------------------------------------------------------------
    !>
    !! Calculates grad = \grad a and hess = hessian(a)
    !! input : f(:)        a real function on the real-space FFT grid
    !! output: grad(3,:)   \grad a, real, on the real-space FFT grid
    !!         hess(3,3,:) hessian(a), real, on the real-space FFT grid
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE hessian_fft_lowpass(this, f, grad, hess)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(core_fft_lowpass), INTENT(IN) :: this
        TYPE(environ_density), INTENT(IN) :: f
        !
        TYPE(environ_gradient), INTENT(INOUT) :: grad
        TYPE(environ_hessian), INTENT(INOUT) :: hess
        !
        INTEGER :: i, j
        COMPLEX(DP), DIMENSION(:), ALLOCATABLE :: aux, gaux, haux
        !
        !--------------------------------------------------------------------------------
        !
        ASSOCIATE (dfft => this%cell%dfft, &
                   nnr => this%cell%dfft%nnr, &
                   nl => this%cell%dfft%nl, &
                   nlm => this%cell%dfft%nlm, &
                   g => this%cell%g)
            !
            !----------------------------------------------------------------------------
            !
            ALLOCATE (aux(nnr))
            ALLOCATE (gaux(nnr))
            ALLOCATE (haux(nnr))
            !
            aux = CMPLX(f%of_r, kind=DP)
            !
            CALL env_fwfft('Rho', aux, dfft) ! bring a(r) to G-space, a(G)
            !
            !----------------------------------------------------------------------------
            ! Multiply by (iG) to get (\grad_ipol a)(G)
            !
            DO i = 1, 3
                gaux = (0.0_DP, 0.0_DP)
                !
                gaux(nl) = g(i, :) * CMPLX(-AIMAG(aux(nl)), REAL(aux(nl)), kind=DP) * &
                           this%filter_lowpass(this%cell%dfft%ngm, this%cell%gg(:))
                !
                IF (dfft%lgamma) &
                    gaux(nlm) = CMPLX(REAL(gaux(nl)), -AIMAG(gaux(nl)), kind=DP)
                !
                CALL env_invfft('Rho', gaux, dfft)
                ! bring back to R-space, (\grad_ipol a)(r)
                !
                grad%of_r(i, :) = tpi * REAL(gaux)
                ! add the factor 2\pi/a missing in the definition of G
                !
                !------------------------------------------------------------------------
                ! Compute the second derivatives
                !
                DO j = 1, i
                    haux = (0.0_DP, 0.0_DP)
                    !
                    haux(nl) = -g(i, :) * g(j, :) * &
                               CMPLX(REAL(aux(nl)), AIMAG(aux(nl)), kind=DP) * &
                               this%filter_lowpass(this%cell%dfft%ngm, this%cell%gg(:))
                    !
                    IF (dfft%lgamma) &
                        haux(nlm) = CMPLX(REAL(haux(nl)), -AIMAG(haux(nl)), kind=DP)
                    !
                    CALL env_invfft('Rho', haux, dfft)
                    ! bring back to R-space (\grad_ipol a)(r)
                    !
                    hess%of_r(i, j, :) = tpi2 * REAL(haux)
                    ! add the factor (2\pi)**2 missing in the definition of G
                    !
                    hess%of_r(j, i, :) = hess%of_r(i, j, :)
                    !
                END DO
                !
            END DO
            !
        END ASSOCIATE
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE hessian_fft_lowpass
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !                               PRIVATE HELPER METHODS
    !
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    FUNCTION filter_lowpass(this, ngm, gg)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(core_fft_lowpass), INTENT(IN) :: this
        !
        REAL(DP) :: filter_lowpass(ngm)
        INTEGER, INTENT(IN) :: ngm
        REAL(DP), INTENT(IN) :: gg(ngm)
        !
        INTEGER :: i
        !
        !--------------------------------------------------------------------------------
        !
        DO i = 1, ngm
            filter_lowpass(i) = 0.5D0 * environ_erfc( &
                ( 8.D0 * gg(i) / this%cell%gcutm - 5.D0 ) * pi_third )
        ENDDO
        !
        !--------------------------------------------------------------------------------
    END FUNCTION filter_lowpass
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE class_core_fft_lowpass
!----------------------------------------------------------------------------------------