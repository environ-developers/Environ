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
MODULE class_core_fft
    !------------------------------------------------------------------------------------
    !
    USE class_io, ONLY: io
    USE env_mp, ONLY: env_mp_sum
    !
    USE env_fft_interfaces, ONLY: env_fwfft, env_invfft
    !
    USE environ_param, ONLY: DP, pi, tpi, tpi2, fpi, e2
    !
    USE class_cell
    USE class_density
    USE class_gradient
    USE class_hessian
    USE class_functions
    !
    USE class_core
    !
    USE tools_math, ONLY: environ_erf, environ_erfc
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
    TYPE, EXTENDS(environ_core), PUBLIC :: core_fft
        !--------------------------------------------------------------------------------
        !
        LOGICAL :: use_internal_pbc_corr = .FALSE.
        !
        REAL(DP) :: alpha = 0.D0
        REAL(DP) :: beta = 0.D0
        !
        REAL(DP), ALLOCATABLE :: correction(:)
        !
        !--------------------------------------------------------------------------------
    CONTAINS
        !--------------------------------------------------------------------------------
        !
        PROCEDURE :: create => create_core_fft
        PROCEDURE :: init => init_core_fft
        PROCEDURE :: update_cell => update_core_fft_cell
        PROCEDURE :: destroy => destroy_core_fft
        !
        PROCEDURE :: gradient => gradient_fft
        PROCEDURE :: divergence => divergence_fft
        PROCEDURE :: laplacian => laplacian_fft
        PROCEDURE :: hessian => hessian_fft
        !
        PROCEDURE :: calc_convolution_density
        PROCEDURE :: calc_convolution_gradient
        PROCEDURE :: calc_convolution_hessian
        !
        PROCEDURE :: poisson => poisson_fft
        PROCEDURE :: grad_poisson => grad_poisson_fft
        PROCEDURE :: force => force_fft
        !
        PROCEDURE :: grad_v_h_of_rho_r
        PROCEDURE :: hess_v_h_of_rho_r
        PROCEDURE :: field_of_grad_rho
        !
        PROCEDURE, PRIVATE :: update_mt_correction
        !
        !--------------------------------------------------------------------------------
    END TYPE core_fft
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
CONTAINS
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !                                   ADMIN METHODS
    !
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE create_core_fft(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(core_fft), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: routine = 'create_core_fft'
        !
        !--------------------------------------------------------------------------------
        !
        IF (ASSOCIATED(this%cell)) CALL io%create_error(routine)
        !
        !--------------------------------------------------------------------------------
        !
        this%core_type = 'fft'
        this%alpha = 0.D0
        this%beta = 0.D0
        !
        NULLIFY (this%cell)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE create_core_fft
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_core_fft(this, cell, use_internal_pbc_corr)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_cell), TARGET, INTENT(IN) :: cell
        LOGICAL, OPTIONAL, INTENT(IN) :: use_internal_pbc_corr
        !
        CLASS(core_fft), INTENT(INOUT) :: this
        !
        !--------------------------------------------------------------------------------
        !
        CALL this%create()
        !
        this%cell => cell
        !
        !--------------------------------------------------------------------------------
        !
        IF (PRESENT(use_internal_pbc_corr)) &
            this%use_internal_pbc_corr = use_internal_pbc_corr
        !
        ALLOCATE (this%correction(cell%dfft%ngm))
        !
        IF (this%use_internal_pbc_corr) THEN
            CALL this%update_mt_correction()
        ELSE
            this%correction = 0.D0
        END IF
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_core_fft
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE update_core_fft_cell(this, cell)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_cell), TARGET, INTENT(IN) :: cell
        !
        CLASS(core_fft), INTENT(INOUT) :: this
        !
        !--------------------------------------------------------------------------------
        !
        this%cell => cell
        !
        IF (this%use_internal_pbc_corr) CALL this%update_mt_correction()
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE update_core_fft_cell
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE destroy_core_fft(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(core_fft), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: routine = 'destroy_core_fft'
        !
        !--------------------------------------------------------------------------------
        !
        IF (.NOT. ASSOCIATED(this%cell)) RETURN
        !
        IF (.NOT. ALLOCATED(this%correction)) CALL io%destroy_error(routine)
        !
        !--------------------------------------------------------------------------------
        !
        DEALLOCATE (this%correction)
        NULLIFY (this%cell)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE destroy_core_fft
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
    SUBROUTINE gradient_fft(this, f, grad)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(core_fft), INTENT(IN) :: this
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
            aux = CMPLX(f%of_r, 0.0_DP, kind=DP)
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
        END ASSOCIATE
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE gradient_fft
    !------------------------------------------------------------------------------------
    !>
    !! Calculates div = \sum_i \grad_i a_i in R-space
    !! input : grad(3,:) a real function on the real-space FFT grid
    !! output: div(:)    \sum_i \grad_i a_i, real, on the real-space FFT grid
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE divergence_fft(this, grad, div)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(core_fft), INTENT(IN) :: this
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
                aux = CMPLX(grad%of_r(i, :), 0.0_DP, kind=DP)
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
                    aux = CMPLX(grad%of_r(i, :), 0.0_DP, kind=DP)
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
    END SUBROUTINE divergence_fft
    !------------------------------------------------------------------------------------
    !>
    !! Calculates lapla = laplacian(a)
    !! input : f(:)       A real function on the real-space FFT grid
    !! output: lapla(:)   \nabla^2 a, real, on the real-space FFT grid
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE laplacian_fft(this, f, lapla)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(core_fft), INTENT(IN) :: this
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
            aux = CMPLX(f%of_r, 0.0_DP, kind=DP)
            !
            CALL env_fwfft('Rho', aux, dfft) ! bring a(r) to G-space, a(G)
            !
            !----------------------------------------------------------------------------
            ! Compute the laplacian
            !
            laux = (0.0_DP, 0.0_DP)
            !
            laux(nl) = -this%cell%gg * aux(nl)
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
    END SUBROUTINE laplacian_fft
    !------------------------------------------------------------------------------------
    !>
    !! Calculates grad = \grad a and hess = hessian(a)
    !! input : f(:)        a real function on the real-space FFT grid
    !! output: grad(3,:)   \grad a, real, on the real-space FFT grid
    !!         hess(3,3,:) hessian(a), real, on the real-space FFT grid
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE hessian_fft(this, f, grad, hess)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(core_fft), INTENT(IN) :: this
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
            aux = CMPLX(f%of_r, 0.0_DP, kind=DP)
            !
            CALL env_fwfft('Rho', aux, dfft) ! bring a(r) to G-space, a(G)
            !
            !----------------------------------------------------------------------------
            ! Multiply by (iG) to get (\grad_ipol a)(G)
            !
            DO i = 1, 3
                gaux = (0.0_DP, 0.0_DP)
                !
                gaux(nl) = g(i, :) * CMPLX(-AIMAG(aux(nl)), REAL(aux(nl)), kind=DP)
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
                               CMPLX(REAL(aux(nl)), AIMAG(aux(nl)), kind=DP)
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
    END SUBROUTINE hessian_fft
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_convolution_density(this, f1, f2, f_out)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(core_fft), INTENT(IN) :: this
        TYPE(environ_density), INTENT(IN) :: f1, f2
        !
        TYPE(environ_density), INTENT(INOUT) :: f_out
        !
        COMPLEX(DP), DIMENSION(:), ALLOCATABLE :: auxr, auxg
        !
        !--------------------------------------------------------------------------------
        !
        ASSOCIATE (dfft => this%cell%dfft, &
                   nnr => this%cell%dfft%nnr, &
                   nl => this%cell%dfft%nl)
            !
            !----------------------------------------------------------------------------
            ! Bring fa and fb to reciprocal space
            !
            ALLOCATE (auxr(nnr))
            auxr = CMPLX(f1%of_r, 0.D0, kind=DP)
            !
            CALL env_fwfft('Rho', auxr, dfft)
            !
            ALLOCATE (auxg(nnr))
            auxg = 0.D0
            !
            auxg(nl) = auxr(nl)
            !
            auxr = CMPLX(f2%of_r, 0.D0, kind=DP)
            !
            CALL env_fwfft('Rho', auxr, dfft)
            !
            !----------------------------------------------------------------------------
            ! Multiply fa(g)*fb(g)
            !
            auxg(nl) = auxg(nl) * auxr(nl)
            !
            DEALLOCATE (auxr)
            !
            IF (dfft%lgamma) &
                auxg(dfft%nlm) = CMPLX(REAL(auxg(nl)), -AIMAG(auxg(nl)), kind=DP)
            !
            !----------------------------------------------------------------------------
            ! Brings convolution back to real space
            !
            CALL env_invfft('Rho', auxg, dfft)
            !
            f_out%of_r = REAL(auxg) * this%cell%omega
            !
        END ASSOCIATE
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_convolution_density
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_convolution_gradient(this, f, grad, grad_out)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(core_fft), INTENT(IN) :: this
        TYPE(environ_density), INTENT(IN) :: f
        TYPE(environ_gradient), INTENT(IN) :: grad
        !
        TYPE(environ_gradient), INTENT(INOUT) :: grad_out
        !
        INTEGER :: i
        !
        TYPE(environ_density) :: local
        COMPLEX(DP), DIMENSION(:), ALLOCATABLE :: auxr, auxg
        !
        !--------------------------------------------------------------------------------
        !
        CALL local%init(f%cell)
        !
        DO i = 1, 3
            local%of_r = grad%of_r(i, :)
            !
            CALL this%calc_convolution_density(f, local, local)
            !
            grad_out%of_r(i, :) = local%of_r
        END DO
        !
        CALL local%destroy()
        !
        ! #TODO check for performance
        !
        ! ASSOCIATE (dfft => this%cell%dfft, &
        !            nnr => this%cell%dfft%nnr, &
        !            nl => this%cell%dfft%nl)
        !     !
        !     !----------------------------------------------------------------------------
        !     ! Bring fa and fb to reciprocal space
        !     !
        !     ALLOCATE (auxr(nnr))
        !     auxr = CMPLX(f%of_r, 0.D0, kind=DP)
        !     !
        !     CALL env_fwfft('Rho', auxr, dfft)
        !     !
        !     ALLOCATE (auxg(nnr))
        !     !
        !     DO i = 1, 3
        !         auxg = CMPLX(grad%of_r(i, :), 0.D0, kind=DP)
        !         !
        !         CALL env_fwfft('Rho', auxg, dfft)
        !         !
        !         !------------------------------------------------------------------------
        !         ! Multiply fa(g)*fb(g)
        !         !
        !         auxg(nl) = auxg(nl) * auxr(nl)
        !         !
        !         IF (dfft%lgamma) &
        !             auxg(dfft%nlm) = CMPLX(REAL(auxg(nl)), -AIMAG(auxg(nl)), kind=DP)
        !         !
        !         !------------------------------------------------------------------------
        !         ! Brings convolution back to real space
        !         !
        !         CALL env_invfft('Rho', auxg, dfft)
        !         !
        !         grad_out%of_r(i, :) = REAL(auxg) * this%cell%omega
        !     END DO
        !     !
        ! END ASSOCIATE
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_convolution_gradient
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_convolution_hessian(this, f, hess, hess_out)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(core_fft), INTENT(IN) :: this
        TYPE(environ_density), INTENT(IN) :: f
        TYPE(environ_hessian), INTENT(IN) :: hess
        !
        TYPE(environ_hessian), INTENT(INOUT) :: hess_out
        !
        COMPLEX(DP), DIMENSION(:), ALLOCATABLE :: auxr, auxg
        !
        INTEGER :: i, j
        !
        TYPE(environ_density) :: local
        !
        !--------------------------------------------------------------------------------
        !
        CALL local%init(f%cell)
        !
        DO i = 1, 3
            !
            DO j = 1, 3
                local%of_r = hess%of_r(i, j, :)
                !
                CALL this%calc_convolution_density(f, local, local)
                !
                hess_out%of_r(i, j, :) = local%of_r
            END DO
            !
        END DO
        !
        CALL local%destroy()
        !
        ! #TODO check for performance
        !
        ! ASSOCIATE (dfft => this%cell%dfft, &
        !            nnr => this%cell%dfft%nnr, &
        !            nl => this%cell%dfft%nl)
        !     !
        !     !----------------------------------------------------------------------------
        !     ! Bring fa and fb to reciprocal space
        !     !
        !     ALLOCATE (auxr(nnr))
        !     auxr = CMPLX(f%of_r, 0.D0, kind=DP)
        !     !
        !     CALL env_fwfft('Rho', auxr, dfft)
        !     !
        !     ALLOCATE (auxg(nnr))
        !     !
        !     DO i = 1, 3
        !         !
        !         DO j = 1, 3
        !             auxg = CMPLX(hess%of_r(i, j, :), 0.D0, kind=DP)
        !             !
        !             CALL env_fwfft('Rho', auxg, dfft)
        !             !
        !             !--------------------------------------------------------------------
        !             ! Multiply fa(g)*fb(g)
        !             !
        !             auxg(nl) = auxg(nl) * auxr(nl)
        !             !
        !             IF (dfft%lgamma) &
        !                 auxg(dfft%nlm) = &
        !                 CMPLX(REAL(auxg(nl)), -AIMAG(auxg(nl)), kind=DP)
        !             !
        !             !--------------------------------------------------------------------
        !             ! Brings convolution back to real space
        !             !
        !             CALL env_invfft('Rho', auxg, dfft)
        !             !
        !             hess_out%of_r(i, j, :) = REAL(auxg) * this%cell%omega
        !         END DO
        !         !
        !     END DO
        !     !
        ! END ASSOCIATE
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_convolution_hessian
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !                               ELECTROSTATIC METHODS
    !
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !>
    !! Solves the Poisson equation \nabla^2 v(r) = - 4 \pi \rho(r)
    !! Input and output functions are defined in real space
    !!
    !! @ param rho : input density
    !! @ param v   : output potential => v(r)
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE poisson_fft(this, rho, v)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(core_fft), INTENT(IN) :: this
        TYPE(environ_density), INTENT(IN) :: rho
        !
        TYPE(environ_density), INTENT(INOUT) :: v
        !
        INTEGER :: i
        COMPLEX(DP), DIMENSION(:), ALLOCATABLE :: auxr, auxg
        !
        !--------------------------------------------------------------------------------
        !
        ASSOCIATE (dfft => this%cell%dfft, &
                   nl => this%cell%dfft%nl, &
                   ngm => this%cell%dfft%ngm)
            !
            !----------------------------------------------------------------------------
            ! Bring rho to G space
            !
            ALLOCATE (auxr(dfft%nnr))
            auxr = CMPLX(rho%of_r, 0.D0, kind=DP)
            !
            CALL env_fwfft('Rho', auxr, dfft)
            !
            ALLOCATE (auxg(ngm))
            auxg = auxr(nl)
            !
            auxr = CMPLX(0.D0, 0.D0, kind=DP)
            !
            auxr(nl(1)) = auxg(1) * this%correction(1) * pi ! G = 0 term
            !
!$omp parallel do
            DO i = this%cell%gstart, ngm
                auxr(nl(i)) = auxg(i) * (1.D0 / this%cell%gg(i) + this%correction(i) * pi)
            END DO
!$omp end parallel do
            !
            auxr = auxr * e2 / pi
            !
            IF (dfft%lgamma) &
                auxr(dfft%nlm) = CMPLX(REAL(auxr(nl)), -AIMAG(auxr(nl)), kind=DP)
            !
            !----------------------------------------------------------------------------
            ! Transform hartree potential to real space
            !
            CALL env_invfft('Rho', auxr, dfft)
            !
            v%of_r = DBLE(auxr)
            !
        END ASSOCIATE
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE poisson_fft
    !------------------------------------------------------------------------------------
    !>
    !! Solves the Poisson equation \nabla^2 v(r) = - 4 \pi \rho(r)
    !! Input and output functions are defined in real space
    !!
    !! @ param rho    : input density
    !! @ param grad_v : output gradient of the potential => \nabla v(r)
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE grad_poisson_fft(this, rho, grad_v)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(core_fft), INTENT(IN) :: this
        TYPE(environ_density), INTENT(IN) :: rho
        !
        TYPE(environ_gradient), INTENT(INOUT) :: grad_v
        !
        INTEGER :: i, j
        REAL(DP) :: fac
        COMPLEX(DP), DIMENSION(:), ALLOCATABLE :: auxr, auxg, vaux
        !
        !--------------------------------------------------------------------------------
        !
        ASSOCIATE (dfft => this%cell%dfft, &
                   nl => this%cell%dfft%nl, &
                   ngm => this%cell%dfft%ngm)
            !
            !----------------------------------------------------------------------------
            ! Bring rho to G space
            !
            ALLOCATE (auxr(dfft%nnr))
            auxr = CMPLX(rho%of_r, 0.D0, KIND=DP)
            !
            CALL env_fwfft('Rho', auxr, dfft)
            !
            ALLOCATE (auxg(ngm))
            auxg = auxr(nl)
            !
            !----------------------------------------------------------------------------
            ! Compute gradient of potential in G space one direction at a time
            !
            DO i = 1, 3
                auxr = CMPLX(0.0_DP, 0.0_DP)
                !
!$omp parallel do private(fac)
                DO j = this%cell%gstart, ngm
                    !
                    auxr(nl(j)) = &
                        this%cell%g(i, j) * &
                        (1.D0 / this%cell%gg(j) + this%correction(j) * pi) * &
                        CMPLX(-AIMAG(auxg(j)), REAL(auxg(j), kind=DP))
                    !
                END DO
!$omp end parallel do
                !
                auxr = auxr * 2.D0 * e2
                ! add the factor 2*e2 coming from the missing prefactor of
                ! V = e2 * 4pi divided by the 2pi factor missing in G
                !
                !------------------------------------------------------------------------
                ! Assuming GAMMA ONLY
                !
                IF (dfft%lgamma) &
                    auxr(dfft%nlm) = CMPLX(REAL(auxr(nl)), -AIMAG(auxr(nl)), kind=DP)
                !
                !------------------------------------------------------------------------
                ! Bring back to R-space, (\grad_ipol a)(r)
                !
                CALL env_invfft('Rho', auxr, dfft)
                !
                grad_v%of_r(i, :) = REAL(auxr)
            END DO
            !
        END ASSOCIATE
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE grad_poisson_fft
    !------------------------------------------------------------------------------------
    !>
    !! Compute the force in G-space using a Gaussian-smeared ion description
    !!
    !! @ param nat   : number of ions
    !! @ param rho   : input density
    !! @ param ions  : object containing the charge, spread, and position of the ions
    !! @ param force : output force
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE force_fft(this, nat, rho, ions, force)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(core_fft), INTENT(IN) :: this
        INTEGER, INTENT(IN) :: nat
        TYPE(environ_density), INTENT(IN) :: rho
        TYPE(environ_functions), INTENT(IN) :: ions
        !
        REAL(DP), INTENT(INOUT) :: force(3, nat)
        !
        INTEGER :: i, j
        !
        REAL(DP) :: fact ! symmetry factor
        REAL(DP) :: fpibg2 ! 4pi / G^2
        REAL(DP) :: t_arg ! G . R
        REAL(DP) :: e_arg ! -D^2 * G^2 / 4
        REAL(DP) :: euler_term ! sin(G . R) * Re(rho(G)) + cos(G . R) * Im(rho(G))
        REAL(DP) :: gauss_term ! Z * 4pi / G^2 * exp(-D^2 * G^2 / 4)
        REAL(DP) :: main_term ! gauss_term + any active PBC correction
        !
        REAL(DP) :: R(3) ! position vector
        !
        COMPLEX(DP), DIMENSION(:), ALLOCATABLE :: auxr, auxg
        !
        CHARACTER(LEN=80) :: routine = 'force_fft'
        !
        !--------------------------------------------------------------------------------
        !
        IF (nat /= ions%number) &
            CALL io%error(routine, &
                          'Mismatch in numbers of atoms passed in input and stored', 1)
        !
        !--------------------------------------------------------------------------------
        !
        ASSOCIATE (dfft => this%cell%dfft, &
                   ngm => this%cell%dfft%ngm, &
                   G => this%cell%g, &
                   G2 => this%cell%gg)
            !
            !----------------------------------------------------------------------------
            ! Bring rho to G space
            !
            ALLOCATE (auxr(dfft%nnr))
            !
            auxr = CMPLX(rho%of_r, 0.D0, KIND=DP)
            !
            CALL env_fwfft('Rho', auxr, dfft)
            !
            ALLOCATE (auxg(ngm))
            auxg = auxr(dfft%nl) ! aux now contains n(G)
            DEALLOCATE (auxr)
            !
            !----------------------------------------------------------------------------
            ! Set symmetry factor
            !
            IF (dfft%lgamma) THEN
                fact = 2.D0
            ELSE
                fact = 1.D0
            END IF
            !
            !----------------------------------------------------------------------------
            ! Compute force
            !
            force = 0.D0
            !
            DO i = 1, nat
                !
                ASSOCIATE (Z => ions%array(i)%volume, & ! ionic charge
                           D => ions%array(i)%spread, & ! gaussian spread of smeared ion
                           R => ions%array(i)%pos - this%cell%origin) ! ion position
                    !
                    DO j = this%cell%gstart, ngm
                        fpibg2 = fpi / (G2(j) * tpi2)
                        !
                        t_arg = tpi * SUM(G(:, j) * R)
                        !
                        euler_term = SIN(t_arg) * DBLE(auxg(j)) + &
                                     COS(t_arg) * AIMAG(auxg(j))
                        !
                        e_arg = -0.25D0 * D**2 * G2(j) * tpi2
                        gauss_term = fpibg2 * EXP(e_arg)
                        !
                        main_term = gauss_term + this%correction(j)
                        !
                        force(:, i) = force(:, i) + G(:, j) * euler_term * main_term
                    END DO
                    !
                    force(:, i) = force(:, i) * Z
                    !
                END ASSOCIATE
                !
            END DO
            !
            force = force * tpi * e2 * fact
            !
            CALL env_mp_sum(force, rho%cell%dfft%comm)
            !
        END ASSOCIATE
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE force_fft
    !------------------------------------------------------------------------------------
    !>
    !!  Gradient of Hartree potential in R space from a total
    !!  (spinless) density in R space n(r)
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE grad_v_h_of_rho_r(this, nnr, rho, grad_v)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(core_fft), INTENT(IN) :: this
        INTEGER, INTENT(IN) :: nnr
        REAL(DP), INTENT(IN) :: rho(nnr)
        !
        REAL(DP), INTENT(OUT) :: grad_v(3, nnr)
        !
        INTEGER :: i, j
        LOGICAL :: gamma_only = .TRUE.
        !
        COMPLEX(DP), DIMENSION(:), ALLOCATABLE :: auxr, auxg
        !
        CHARACTER(LEN=80) :: routine = 'grad_v_h_of_rho_r'
        !
        !--------------------------------------------------------------------------------
        !
        IF (nnr /= this%cell%nnr) CALL io%error(routine, "Mismatch in FFT domain", 1)
        !
        !--------------------------------------------------------------------------------
        !
        ASSOCIATE (dfft => this%cell%dfft, &
                   nl => this%cell%dfft%nl, &
                   g => this%cell%g, &
                   gg => this%cell%gg)
            !
            !----------------------------------------------------------------------------
            ! Bring rho to G space
            !
            ALLOCATE (auxr(nnr))
            auxr = CMPLX(rho, 0.D0, KIND=DP)
            !
            CALL env_fwfft('Rho', auxr, dfft)
            !
            !----------------------------------------------------------------------------
            ! Compute total potential in G space
            !
            ALLOCATE (auxg(nnr))
            !
            DO i = 1, 3
                auxg = (0.0_DP, 0.0_DP)
                !
                DO j = this%cell%gstart, dfft%ngm
                    !
                    auxg(nl(j)) = &
                        g(i, j) * (1.D0 / gg(j) + this%correction(j) * tpi2) * &
                        CMPLX(-AIMAG(auxr(nl(j))), REAL(auxr(nl(j))), kind=DP)
                    !
                END DO
                !
                auxg = auxg * e2 * fpi / tpi
                ! add the factor e2*fpi/2\pi coming from the missing prefactor of
                ! V = e2 * fpi divided by the 2\pi factor missing in G
                !
                IF (gamma_only) THEN
                    !
                    auxg(dfft%nlm) = &
                        CMPLX(REAL(auxg(nl)), -AIMAG(auxg(nl)), kind=DP)
                    !
                END IF
                !
                CALL env_invfft('Rho', auxg, dfft) ! bring back to R-space
                !
                grad_v(i, :) = REAL(auxg)
            END DO
            !
        END ASSOCIATE
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE grad_v_h_of_rho_r
    !------------------------------------------------------------------------------------
    !>
    !! Gradient of Hartree potential in R space from a total
    !! (spinless) density in R space n(r)
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE hess_v_h_of_rho_r(this, nnr, rho, hess_v)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(core_fft), INTENT(IN) :: this
        INTEGER, INTENT(IN) :: nnr
        REAL(DP), INTENT(IN) :: rho(nnr)
        !
        REAL(DP), INTENT(OUT) :: hess_v(3, 3, nnr)
        !
        INTEGER :: i, j, k
        LOGICAL :: gamma_only = .TRUE.
        !
        COMPLEX(DP), DIMENSION(:), ALLOCATABLE :: auxr, auxg
        !
        CHARACTER(LEN=80) :: routine = 'hess_v_h_of_rho_r'
        !
        !--------------------------------------------------------------------------------
        !
        IF (nnr /= this%cell%nnr) CALL io%error(routine, "Mismatch in FFT domain", 1)
        !
        !--------------------------------------------------------------------------------
        !
        ASSOCIATE (dfft => this%cell%dfft, &
                   nl => this%cell%dfft%nl, &
                   g => this%cell%g)
            !
            !----------------------------------------------------------------------------
            ! Bring rho to G space
            !
            ALLOCATE (auxr(nnr))
            auxr = CMPLX(rho, 0.D0, KIND=DP)
            !
            CALL env_fwfft('Rho', auxr, dfft)
            !
            !----------------------------------------------------------------------------
            ! Compute total potential in G space
            !
            ALLOCATE (auxg(nnr))
            !
            DO i = 1, 3
                !
                DO j = 1, 3
                    auxg = (0.0_DP, 0.0_DP)
                    !
                    DO k = this%cell%gstart, dfft%ngm
                        !
                        auxg(nl(k)) = &
                            g(i, k) * g(j, k) * &
                            (1.D0 / this%cell%gg(k) + this%correction(k) * tpi2) * &
                            CMPLX(DBLE(auxr(nl(k))), AIMAG(auxr(nl(k))), kind=DP)
                        !
                    END DO
                    !
                    auxg = auxg * e2 * fpi
                    ! add the factor e2*fpi coming from the missing prefactor of
                    ! V = e2 * fpi
                    !
                    IF (gamma_only) &
                        auxg(dfft%nlm) = CMPLX(DBLE(auxg(nl)), -AIMAG(auxg(nl)), kind=DP)
                    !
                    CALL env_invfft('Rho', auxg, dfft) ! bring back to R-space
                    !
                    hess_v(i, j, :) = REAL(auxg)
                END DO
                !
            END DO
            !
        END ASSOCIATE
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE hess_v_h_of_rho_r
    !------------------------------------------------------------------------------------
    !>
    !! Gradient of Hartree potential in R space from a total
    !! (spinless) density in R space n(r)
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE field_of_grad_rho(this, nnr, grad_rho, field)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(core_fft), INTENT(IN) :: this
        INTEGER, INTENT(IN) :: nnr
        REAL(DP), INTENT(IN) :: grad_rho(3, nnr)
        !
        REAL(DP), INTENT(OUT) :: field(nnr)
        !
        INTEGER :: i, j
        LOGICAL :: gamma_only = .TRUE.
        !
        COMPLEX(DP), DIMENSION(:), ALLOCATABLE :: auxr, auxg, auxe
        !
        CHARACTER(LEN=80) :: routine = 'field_of_grad_rho'
        !
        !--------------------------------------------------------------------------------
        !
        IF (nnr /= this%cell%nnr) CALL io%error(routine, "Mismatch in FFT domain", 1)
        !
        !--------------------------------------------------------------------------------
        !
        ASSOCIATE (dfft => this%cell%dfft, &
                   nl => this%cell%dfft%nl)
            !
            !----------------------------------------------------------------------------
            ! Bring gradrho to G space
            !
            ALLOCATE (auxe(nnr))
            auxe = CMPLX(0.D0, 0.D0, KIND=DP)
            !
            ALLOCATE (auxr(nnr))
            auxr = CMPLX(0.D0, 0.D0, KIND=DP)
            !
            ALLOCATE (auxg(nnr))
            !
            DO i = 1, 3
                auxg = CMPLX(grad_rho(i, :), 0.D0, KIND=DP)
                !
                CALL env_fwfft('Rho', auxg, dfft)
                !
                !------------------------------------------------------------------------
                ! Compute total potential in G space
                !
                DO j = this%cell%gstart, dfft%ngm
                    !
                    auxr(nl(j)) = &
                        this%cell%g(i, j) * &
                        (1.D0 / this%cell%gg(j) + this%correction(j) * pi) * &
                        CMPLX(-AIMAG(auxg(nl(j))), REAL(auxg(nl(j))), kind=DP)
                    !
                END DO
                !
                auxr = auxr * 2.D0 * e2
                ! add the factor 2*e2 coming from the missing prefactor of
                ! V = e2 * 4pi divided by the 2pi factor missing in G
                !
                auxe = auxe + auxr
            END DO
            !
            IF (gamma_only) &
                auxe(dfft%nlm) = CMPLX(REAL(auxe(nl)), -AIMAG(auxe(nl)), kind=DP)
            !
            CALL env_invfft('Rho', auxe, dfft)
            ! bring back to R-space (\grad_ipol a)(r)
            !
            field = REAL(auxe)
            !
        END ASSOCIATE
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE field_of_grad_rho
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
    SUBROUTINE update_mt_correction(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(core_fft), INTENT(INOUT) :: this
        !
        INTEGER :: i
        LOGICAL :: physical
        REAL(DP) :: r(3), rws, upperbound, ecutrho
        COMPLEX(DP), ALLOCATABLE :: aux(:)
        !
        CHARACTER(LEN=80) :: routine = 'update_mt_correction'
        !
        !--------------------------------------------------------------------------------
        !
        ASSOCIATE (alpha => this%alpha, &
                   beta => this%beta, &
                   mt_corr => this%correction, &
                   cell => this%cell, &
                   dfft => this%cell%dfft)
            !
            !----------------------------------------------------------------------------
            !
            ecutrho = cell%gcutm * tpi2
            !
            !----------------------------------------------------------------------------
            ! choose alpha in order to have convergence in the sum over G
            ! upperbound is a safe upper bound for the error in the sum over G
            !
            alpha = 2.9D0
            upperbound = 1._DP
            !
            DO WHILE (upperbound > 1.E-7_DP)
                alpha = alpha - 0.1_DP
                !
                IF (alpha <= 0._DP) CALL io%error(routine, "Optimal alpha not found", 1)
                !
                upperbound = e2 * SQRT(2.D0 * alpha / tpi) * &
                             environ_erfc(SQRT(ecutrho / 4.D0 / alpha))
                !
            END DO
            !
            beta = 0.5_DP / alpha ! 1._dp/alpha
            !
            ALLOCATE (aux(dfft%nnr))
            aux = (0._DP, 0._DP)
            !
            DO i = 1, cell%ir_end
                !
                CALL cell%get_min_distance(i, 0, 0, cell%origin, r, rws, physical)
                ! compute minimum distance using minimum image convention
                !
                IF (.NOT. physical) CYCLE
                !
                aux(i) = smooth_coulomb_r(alpha, SQRT(rws))
                !
            END DO
            !
            CALL env_fwfft('Rho', aux, dfft)
            !
!$omp parallel do
            DO i = 1, cell%dfft%ngm
                !
                mt_corr(i) = cell%omega * REAL(aux(dfft%nl(i))) - &
                             smooth_coulomb_g(alpha, beta, tpi2 * cell%gg(i))
                !
            END DO
!$omp end parallel do
            !
            mt_corr = mt_corr * EXP(-tpi2 * cell%gg * beta / 4._DP)**2
            !
        END ASSOCIATE
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE update_mt_correction
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    REAL(DP) FUNCTION smooth_coulomb_r(alpha, r)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP), INTENT(IN) :: alpha, r
        !
        !--------------------------------------------------------------------------------
        !
        IF (r > 1.E-6_DP) THEN
            smooth_coulomb_r = environ_erf(SQRT(alpha) * r) / r
        ELSE
            smooth_coulomb_r = 2._DP / SQRT(pi) * SQRT(alpha)
        END IF
        !
        !--------------------------------------------------------------------------------
    END FUNCTION smooth_coulomb_r
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    REAL(DP) FUNCTION smooth_coulomb_g(alpha, beta, q2)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP), INTENT(IN) :: alpha, beta, q2
        !
        !--------------------------------------------------------------------------------
        !
        IF (q2 > 1.E-6_DP) THEN
            smooth_coulomb_g = fpi * EXP(-q2 / 4._DP / alpha) / q2
        ELSE
            !
            smooth_coulomb_g = -1._DP * fpi * &
                               (1._DP / 4._DP / alpha + 2._DP * beta / 4._DP)
            !
        END IF
        !
        !--------------------------------------------------------------------------------
    END FUNCTION smooth_coulomb_g
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE class_core_fft
!----------------------------------------------------------------------------------------
