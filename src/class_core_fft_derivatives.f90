!----------------------------------------------------------------------------------------
!
! Copyright (C) 2018-2021 ENVIRON (www.quantum-environ.org)
!
!----------------------------------------------------------------------------------------
!
!     This file is part of Environ version 2.0
!
!     Environ 2.0 is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation, either version 2 of the License, or
!     (at your option) any later version.
!
!     Environ 2.0 is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more detail, either the file
!     `License' in the root directory of the present distribution, or
!     online at <http://www.gnu.org/licenses/>.
!
!----------------------------------------------------------------------------------------
!
! Authors: Oliviero Andreussi (Department of Physics, UNT)
!          Francesco Nattino  (THEOS and NCCR-MARVEL, EPFL)
!          Nicola Marzari     (THEOS and NCCR-MARVEL, EPFL)
!          Edan Bainglass     (Department of Physics, UNT)
!
!----------------------------------------------------------------------------------------
!>
!!
!----------------------------------------------------------------------------------------
MODULE class_core_fft_derivatives
    !------------------------------------------------------------------------------------
    !
    USE env_types_fft, ONLY: env_fft_type_descriptor
    USE env_fft_main, ONLY: env_fwfft, env_invfft
    !
    USE environ_param, ONLY: DP, tpi, tpi2
    !
    USE class_cell
    USE class_density
    USE class_gradient
    USE class_hessian
    !
    USE class_core_fft
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
    TYPE, EXTENDS(core_fft), PUBLIC :: core_fft_derivatives
        !--------------------------------------------------------------------------------
        !
        !--------------------------------------------------------------------------------
    CONTAINS
        !--------------------------------------------------------------------------------
        !
        PROCEDURE :: gradient => gradient_fft
        PROCEDURE :: graddot => graddot_fft
        PROCEDURE :: laplacian => laplacian_fft
        PROCEDURE :: hessian => hessian_fft
        !
        PROCEDURE :: convolution_density, convolution_gradient, convolution_hessian
        !
        !--------------------------------------------------------------------------------
    END TYPE core_fft_derivatives
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
CONTAINS
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !                                  GENERAL METHODS
    !
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !>
    !! Calculates ga = \grad a
    !! input : a(:)     a real function on the real-space FFT grid
    !! output: ga(3,:)  \grad a, real, on the real-space FFT grid
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE gradient_fft(this, a, ga)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(core_fft_derivatives), TARGET, INTENT(IN) :: this
        TYPE(environ_density), INTENT(IN) :: a
        !
        TYPE(environ_gradient), INTENT(INOUT) :: ga
        !
        INTEGER :: ipol
        COMPLEX(DP), DIMENSION(:), ALLOCATABLE :: aux, gaux
        !
        TYPE(env_fft_type_descriptor), POINTER :: dfft
        REAL(DP), POINTER :: g(:, :)
        !
        !--------------------------------------------------------------------------------
        !
        dfft => this%cell%dfft
        g => this%g
        !
        ALLOCATE (aux(dfft%nnr))
        ALLOCATE (gaux(dfft%nnr))
        !
        aux = CMPLX(a%of_r(:), 0.0_DP, kind=DP)
        !
        CALL env_fwfft(aux, dfft) ! bring a(r) to G-space, a(G)
        !
        !--------------------------------------------------------------------------------
        ! Multiply by (iG) to get (\grad_ipol a)(G)
        !
        DO ipol = 1, 3
            gaux(:) = (0.0_DP, 0.0_DP)
            !
            gaux(dfft%nl(:)) = g(ipol, :) * CMPLX(-AIMAG(aux(dfft%nl(:))), &
                                                  REAL(aux(dfft%nl(:))), kind=DP)
            !
            IF (dfft%lgamma) THEN
                !
                gaux(dfft%nlm(:)) = CMPLX(REAL(gaux(dfft%nl(:))), &
                                          -AIMAG(gaux(dfft%nl(:))), kind=DP)
                !
            END IF
            !
            CALL env_invfft(gaux, dfft) ! bring back to R-space, (\grad_ipol a)(r)
            !
            ga%of_r(ipol, :) = tpi * DBLE(gaux(:))
            ! add the factor 2\pi/a missing in the definition of G
            !
        END DO
        !
        DEALLOCATE (gaux)
        DEALLOCATE (aux)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE gradient_fft
    !------------------------------------------------------------------------------------
    !>
    !! Calculates da = \sum_i \grad_i a_i in R-space
    !! input : ga(3,:)  a real function on the real-space FFT grid
    !! output: da(:)    \sum_i \grad_i a_i, real, on the real-space FFT grid
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE graddot_fft(this, ga, da)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(core_fft_derivatives), TARGET, INTENT(IN) :: this
        TYPE(environ_gradient), INTENT(IN) :: ga
        !
        TYPE(environ_density), INTENT(INOUT) :: da
        !
        INTEGER :: n, ipol
        COMPLEX(DP), DIMENSION(:), ALLOCATABLE :: aux, gaux
        COMPLEX(DP) :: fp, fm, aux1, aux2
        !
        TYPE(env_fft_type_descriptor), POINTER :: dfft
        REAL(DP), POINTER :: g(:, :)
        !
        !--------------------------------------------------------------------------------
        !
        dfft => this%cell%dfft
        g => this%g
        !
        ALLOCATE (aux(dfft%nnr))
        ALLOCATE (gaux(dfft%nnr))
        !
        gaux(:) = (0.0_DP, 0.0_DP)
        !
        IF (dfft%lgamma) THEN
            !
            !----------------------------------------------------------------------------
            ! Gamma tricks: perform 2 FFT's in a single shot
            !
            ! x and y
            !
            ipol = 1
            aux(:) = CMPLX(ga%of_r(ipol, :), ga%of_r(ipol + 1, :), kind=DP)
            !
            CALL env_fwfft(aux, dfft) ! bring a(ipol,r) to G-space, a(G)
            !
            !----------------------------------------------------------------------------
            ! Multiply by iG to get the gradient in G-space
            !
            DO n = 1, dfft%ngm
                fp = (aux(dfft%nl(n)) + aux(dfft%nlm(n))) * 0.5_DP
                fm = (aux(dfft%nl(n)) - aux(dfft%nlm(n))) * 0.5_DP
                aux1 = CMPLX(REAL(fp), AIMAG(fm), kind=DP)
                aux2 = CMPLX(AIMAG(fp), -REAL(fm), kind=DP)
                !
                gaux(dfft%nl(n)) = CMPLX(0.0_DP, g(ipol, n), kind=DP) * aux1 + &
                                   CMPLX(0.0_DP, g(ipol + 1, n), kind=DP) * aux2
                !
            END DO
            !
            !----------------------------------------------------------------------------
            ! z
            !
            ipol = 3
            aux(:) = CMPLX(ga%of_r(ipol, :), 0.0_DP, kind=DP)
            !
            CALL env_fwfft(aux, dfft) ! bring a(ipol, r) to G-space, a(G)
            !
            !----------------------------------------------------------------------------
            ! Multiply by iG to get the gradient in G-space
            ! Fill both gaux(G) and gaux(-G) = gaux*(G)
            !
            DO n = 1, dfft%ngm
                !
                gaux(dfft%nl(n)) = gaux(dfft%nl(n)) + &
                                   g(ipol, n) * CMPLX(-AIMAG(aux(dfft%nl(n))), &
                                                      REAL(aux(dfft%nl(n))), kind=DP)
                !
                gaux(dfft%nlm(n)) = CONJG(gaux(dfft%nl(n)))
            END DO
            !
        ELSE
            !
            DO ipol = 1, 3
                aux = CMPLX(ga%of_r(ipol, :), 0.0_DP, kind=DP)
                !
                CALL env_fwfft(aux, dfft) ! bring a(ipol,r) to G-space, a(G)
                !
                !------------------------------------------------------------------------
                ! Multiply by iG to get the gradient in G-space
                !
                DO n = 1, dfft%ngm
                    !
                    gaux(dfft%nl(n)) = gaux(dfft%nl(n)) + &
                                       g(ipol, n) * CMPLX(-AIMAG(aux(dfft%nl(n))), &
                                                          REAL(aux(dfft%nl(n))), kind=DP)
                    !
                END DO
                !
            END DO
            !
        END IF
        !
        CALL env_invfft(gaux, dfft) ! bring back to R-space, (\grad_ipol a)(r)
        !
        da%of_r(:) = tpi * REAL(gaux(:))
        ! Add the factor 2\pi/a missing in the definition of G and sum
        !
        DEALLOCATE (aux, gaux)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE graddot_fft
    !------------------------------------------------------------------------------------
    !>
    !! Calculates lapla = laplacian(a)
    !! input : a(:)     A real function on the real-space FFT grid
    !! output: lapla(:) \nabla^2 a, real, on the real-space FFT grid
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE laplacian_fft(this, a, lapla)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(core_fft_derivatives), TARGET, INTENT(IN) :: this
        TYPE(environ_density), INTENT(IN) :: a
        !
        TYPE(environ_density), INTENT(INOUT) :: lapla
        !
        INTEGER :: ig
        COMPLEX(DP), DIMENSION(:), ALLOCATABLE :: aux, laux
        !
        REAL(DP), POINTER :: gg(:)
        TYPE(env_fft_type_descriptor), POINTER :: dfft
        !
        !--------------------------------------------------------------------------------
        !
        gg => this%gg
        dfft => this%cell%dfft
        !
        ALLOCATE (aux(dfft%nnr))
        ALLOCATE (laux(dfft%nnr))
        !
        aux = CMPLX(a%of_r(:), 0.0_DP, kind=DP)
        !
        CALL env_fwfft(aux, dfft) ! bring a(r) to G-space, a(G)
        !
        !--------------------------------------------------------------------------------
        ! Compute the laplacian
        !
        laux(:) = (0.0_DP, 0.0_DP)
        !
        laux(dfft%nl(:)) = -gg(:) * aux(dfft%nl(:))
        !
        IF (dfft%lgamma) THEN
            !
            laux(dfft%nlm(:)) = CMPLX(REAL(laux(dfft%nl(:))), &
                                      -AIMAG(laux(dfft%nl(:))), kind=DP)
            !
        END IF
        !
        CALL env_invfft(laux, dfft) ! bring back to R-space, (\lapl a)(r)
        !
        lapla%of_r = tpi2 * REAL(laux) ! add the missing factor (2\pi)^2 in G
        !
        DEALLOCATE (laux)
        DEALLOCATE (aux)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE laplacian_fft
    !------------------------------------------------------------------------------------
    !>
    !! Calculates ga = \grad a and ha = hessian(a)
    !! input : a(:)     a real function on the real-space FFT grid
    !! output: ga(3,:)  \grad a, real, on the real-space FFT grid
    !!         ha(3,3,:)  hessian(a), real, on the real-space FFT grid
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE hessian_fft(this, a, ga, ha)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(core_fft_derivatives), TARGET, INTENT(IN) :: this
        TYPE(environ_density), INTENT(IN) :: a
        !
        TYPE(environ_gradient), INTENT(INOUT) :: ga
        TYPE(environ_hessian), INTENT(INOUT) :: ha
        !
        INTEGER :: ipol, jpol
        COMPLEX(DP), DIMENSION(:), ALLOCATABLE :: aux, gaux, haux
        !
        REAL(DP), POINTER :: g(:, :)
        TYPE(env_fft_type_descriptor), POINTER :: dfft
        !
        !--------------------------------------------------------------------------------
        !
        g => this%g
        dfft => this%cell%dfft
        !
        ALLOCATE (aux(dfft%nnr))
        ALLOCATE (gaux(dfft%nnr))
        ALLOCATE (haux(dfft%nnr))
        !
        aux = CMPLX(a%of_r(:), 0.0_DP, kind=DP)
        !
        CALL env_fwfft(aux, dfft) ! bring a(r) to G-space, a(G)
        !
        !--------------------------------------------------------------------------------
        ! Multiply by (iG) to get (\grad_ipol a)(G)
        !
        DO ipol = 1, 3
            gaux(:) = (0.0_DP, 0.0_DP)
            !
            gaux(dfft%nl(:)) = g(ipol, :) * CMPLX(-AIMAG(aux(dfft%nl(:))), &
                                                  REAL(aux(dfft%nl(:))), kind=DP)
            !
            IF (dfft%lgamma) THEN
                !
                gaux(dfft%nlm(:)) = CMPLX(REAL(gaux(dfft%nl(:))), &
                                          -AIMAG(gaux(dfft%nl(:))), kind=DP)
                !
            END IF
            !
            CALL env_invfft(gaux, dfft) ! bring back to R-space, (\grad_ipol a)(r)
            !
            ga%of_r(ipol, :) = tpi * REAL(gaux(:))
            ! add the factor 2\pi/a missing in the definition of G
            !
            !----------------------------------------------------------------------------
            ! Compute the second derivatives
            !
            DO jpol = 1, ipol
                haux(:) = (0.0_DP, 0.0_DP)
                !
                haux(dfft%nl(:)) = -g(ipol, :) * g(jpol, :) * &
                                   CMPLX(REAL(aux(dfft%nl(:))), &
                                         AIMAG(aux(dfft%nl(:))), kind=DP)
                !
                IF (dfft%lgamma) THEN
                    !
                    haux(dfft%nlm(:)) = CMPLX(REAL(haux(dfft%nl(:))), &
                                              -AIMAG(haux(dfft%nl(:))), kind=DP)
                    !
                END IF
                !
                CALL env_invfft(haux, dfft)
                ! bring back to R-space (\grad_ipol a)(r)
                !
                ha%of_r(ipol, jpol, :) = tpi2 * REAL(haux(:))
                ! add the factor (2\pi)**2 missing in the definition of G
                !
                ha%of_r(jpol, ipol, :) = ha%of_r(ipol, jpol, :)
                !
            END DO
            !
        END DO
        !
        DEALLOCATE (haux)
        DEALLOCATE (gaux)
        DEALLOCATE (aux)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE hessian_fft
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE convolution_density(this, fa, fb, fc)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(core_fft_derivatives), TARGET, INTENT(IN) :: this
        TYPE(environ_density), INTENT(IN) :: fa, fb
        !
        TYPE(environ_density), INTENT(INOUT) :: fc
        !
        COMPLEX(DP), DIMENSION(:), ALLOCATABLE :: auxr, auxg
        !
        REAL(DP), POINTER :: omega
        TYPE(env_fft_type_descriptor), POINTER :: dfft
        !
        !--------------------------------------------------------------------------------
        !
        dfft => this%cell%dfft
        omega => this%cell%omega
        !
        !--------------------------------------------------------------------------------
        ! Bring fa and fb to reciprocal space
        !
        ALLOCATE (auxr(dfft%nnr))
        auxr(:) = CMPLX(fa%of_r(:), 0.D0, kind=DP)
        !
        CALL env_fwfft(auxr, dfft)
        !
        ALLOCATE (auxg(dfft%nnr))
        auxg = 0.D0
        !
        auxg(dfft%nl(:)) = auxr(dfft%nl(:))
        !
        auxr(:) = CMPLX(fb%of_r(:), 0.D0, kind=DP)
        !
        CALL env_fwfft(auxr, dfft)
        !
        !--------------------------------------------------------------------------------
        ! Multiply fa(g)*fb(g)
        !
        auxg(dfft%nl(:)) = auxg(dfft%nl(:)) * auxr(dfft%nl(:))
        !
        DEALLOCATE (auxr)
        !
        IF (dfft%lgamma) auxg(dfft%nlm(:)) = &
            CMPLX(REAL(auxg(dfft%nl(:))), -AIMAG(auxg(dfft%nl(:))), kind=DP)
        !
        !--------------------------------------------------------------------------------
        ! Brings convolution back to real space
        !
        CALL env_invfft(auxg, dfft)
        !
        fc%of_r(:) = REAL(auxg(:)) * omega
        !
        DEALLOCATE (auxg)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE convolution_density
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE convolution_gradient(this, fa, gb, gc)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(core_fft_derivatives), TARGET, INTENT(IN) :: this
        TYPE(environ_density), INTENT(IN) :: fa
        TYPE(environ_gradient), INTENT(IN) :: gb
        !
        TYPE(environ_gradient), INTENT(INOUT) :: gc
        !
        TYPE(environ_density) :: local
        COMPLEX(DP), DIMENSION(:), ALLOCATABLE :: auxr, auxg
        !
        REAL(DP), POINTER :: omega
        TYPE(env_fft_type_descriptor), POINTER :: dfft
        !
        INTEGER :: ipol
        !
        !--------------------------------------------------------------------------------
        !
        CALL local%init(fa%cell)
        !
        DO ipol = 1, 3
            local%of_r(:) = gb%of_r(ipol, :)
            !
            CALL this%convolution_density(fa, local, local)
            !
            gc%of_r(ipol, :) = local%of_r(:)
        END DO
        !
        CALL local%destroy()
        !
        ! #TODO check for performance
        !
        ! dfft => this%cell%dfft
        ! omega => this%cell%omega
        ! !
        ! ! Bring fa and fb to reciprocal space
        ! !
        ! ALLOCATE (auxr(dfft%nnr))
        ! auxr(:) = CMPLX(fa%of_r(:), 0.D0, kind=DP)
        ! !
        ! CALL env_fwfft(auxr, dfft)
        ! !
        ! ALLOCATE (auxg(dfft%nnr))
        ! !
        ! DO ipol = 1, 3
        !     !
        !     auxg(:) = CMPLX(gb%of_r(ipol, :), 0.D0, kind=DP)
        !     !
        !     CALL env_fwfft(auxg, dfft)
        !     !
        !     ! Multiply fa(g)*fb(g)
        !     !
        !     auxg(dfft%nl(:)) = auxg(dfft%nl(:)) * auxr(dfft%nl(:))
        !     !
        !     IF (dfft%lgamma) auxg(dfft%nlm(:)) = &
        !         CMPLX(REAL(auxg(dfft%nl(:))), -AIMAG(auxg(dfft%nl(:))), kind=DP)
        !     !
        !     ! Brings convolution back to real space
        !     !
        !     CALL env_invfft(auxg, dfft)
        !     !
        !     gc%of_r(ipol, :) = REAL(auxg(:)) * omega
        !     !
        ! END DO
        ! !
        ! DEALLOCATE (auxr)
        ! DEALLOCATE (auxg)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE convolution_gradient
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE convolution_hessian(this, fa, hb, hc)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(core_fft_derivatives), TARGET, INTENT(IN) :: this
        TYPE(environ_density), INTENT(IN) :: fa
        TYPE(environ_hessian), INTENT(IN) :: hb
        !
        TYPE(environ_hessian), INTENT(INOUT) :: hc
        !
        COMPLEX(DP), DIMENSION(:), ALLOCATABLE :: auxr, auxg
        !
        REAL(DP), POINTER :: omega
        TYPE(env_fft_type_descriptor), POINTER :: dfft
        TYPE(environ_density) :: local
        !
        INTEGER :: ipol, jpol
        !
        !--------------------------------------------------------------------------------
        !
        CALL local%init(fa%cell)
        !
        DO ipol = 1, 3
            !
            DO jpol = 1, 3
                local%of_r(:) = hb%of_r(ipol, jpol, :)
                !
                CALL this%convolution_density(fa, local, local)
                !
                hc%of_r(ipol, jpol, :) = local%of_r(:)
            END DO
            !
        END DO
        !
        CALL local%destroy()
        !
        ! #TODO check for performance
        !
        ! dfft => this%cell%dfft
        ! omega => this%cell%omega
        ! !
        ! ! Bring fa and fb to reciprocal space
        ! !
        ! ALLOCATE (auxr(dfft%nnr))
        ! auxr(:) = CMPLX(fa%of_r(:), 0.D0, kind=DP)
        ! !
        ! CALL env_fwfft(auxr, dfft)
        ! !
        ! ALLOCATE (auxg(dfft%nnr))
        ! !
        ! DO ipol = 1, 3
        !     DO jpol = 1, 3
        !         !
        !         auxg(:) = CMPLX(hb%of_r(ipol, jpol, :), 0.D0, kind=DP)
        !         !
        !         CALL env_fwfft(auxg, dfft)
        !         !
        !         ! Multiply fa(g)*fb(g)
        !         !
        !         auxg(dfft%nl(:)) = auxg(dfft%nl(:)) * auxr(dfft%nl(:))
        !         !
        !         IF (dfft%lgamma) auxg(dfft%nlm(:)) = &
        !             CMPLX(REAL(auxg(dfft%nl(:))), -AIMAG(auxg(dfft%nl(:))), kind=DP)
        !         !
        !         ! Brings convolution back to real space
        !         !
        !         CALL env_invfft(auxg, dfft)
        !         !
        !         hc%of_r(ipol, jpol, :) = REAL(auxg(:)) * omega
        !         !
        !     END DO
        ! END DO
        ! !
        ! DEALLOCATE (auxr)
        ! DEALLOCATE (auxg)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE convolution_hessian
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE class_core_fft_derivatives
!----------------------------------------------------------------------------------------
