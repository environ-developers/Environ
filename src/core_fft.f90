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
    USE env_sorting, ONLY: env_hpsort_eps
    USE env_mp, ONLY: env_mp_sum
    !
    USE env_types_fft, ONLY: env_fft_type_descriptor, env_fft_stick_index
    USE env_fft_main, ONLY: env_fwfft, env_invfft
    USE env_fft_ggen, ONLY: env_fft_set_nl
    !
    USE environ_param, ONLY: DP, pi, tpi, tpi2, fpi, e2, eps8
    !
    USE class_cell
    USE class_density
    USE class_gradient
    USE class_hessian
    USE class_function
    USE class_function_gaussian
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
        INTEGER :: ngm = 0 ! local  number of G vectors (on this processor)
        ! with gamma tricks, only vectors in G>
        !
        REAL(DP) :: gcutm = 0.0_DP ! ecutrho/(2 pi/a)^2, cut-off for |G|^2
        !
        INTEGER :: gstart = 2 ! index of the first G vector whose module is > 0
        ! needed in parallel execution:
        ! gstart=2 for the proc that holds G=0
        ! gstart=1 for all others
        !
        REAL(DP), ALLOCATABLE :: gg(:)
        ! G^2 in increasing order (in units of (2pi/a)^2)
        !
        REAL(DP), ALLOCATABLE :: g(:, :)
        ! G-vectors cartesian components (in units 2pi/a)
        !
        !--------------------------------------------------------------------------------
        ! Electrostatics
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
        PROCEDURE :: graddot => graddot_fft
        PROCEDURE :: laplacian => laplacian_fft
        PROCEDURE :: hessian => hessian_fft
        !
        PROCEDURE :: convolution_density, convolution_gradient, convolution_hessian
        !
        PROCEDURE :: poisson => poisson_fft
        PROCEDURE :: gradpoisson => gradpoisson_fft
        PROCEDURE :: force => force_fft
        !
        PROCEDURE :: hessv_h_of_rho_r, field_of_gradrho
        !
        PROCEDURE, PRIVATE :: update_mt_correction
        PROCEDURE, PRIVATE :: init_gvect => env_gvect_init
        PROCEDURE, PRIVATE :: deallocate_gvect => env_deallocate_gvect
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
        CHARACTER(LEN=80) :: sub_name = 'create_core_fft'
        !
        !--------------------------------------------------------------------------------
        !
        IF (ASSOCIATED(this%cell)) CALL io%create_error(sub_name)
        !
        IF (ALLOCATED(this%g)) CALL io%create_error(sub_name)
        !
        IF (ALLOCATED(this%gg)) CALL io%create_error(sub_name)
        !
        !--------------------------------------------------------------------------------
        !
        this%core_type = 'fft'
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE create_core_fft
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_core_fft(this, gcutm, cell, use_internal_pbc_corr)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP), INTENT(IN) :: gcutm
        TYPE(environ_cell), TARGET, INTENT(IN) :: cell
        LOGICAL, INTENT(IN), OPTIONAL :: use_internal_pbc_corr
        !
        CLASS(core_fft), INTENT(INOUT) :: this
        !
        INTEGER :: ngm_g ! global number of G vectors (summed on all procs)
        ! in serial execution, ngm_g = ngm
        !
        !--------------------------------------------------------------------------------
        !
        CALL this%create()
        !
        this%gcutm = gcutm
        this%cell => cell
        this%ngm = cell%dfft%ngm
        !
        !--------------------------------------------------------------------------------
        !
        CALL this%init_gvect(ngm_g, cell%dfft%comm)
        !
        CALL env_ggen(this%cell%dfft, cell%dfft%comm, cell%at, cell%bg, this%gcutm, &
                      ngm_g, this%ngm, this%g, this%gg, this%gstart, .TRUE.)
        !
        !--------------------------------------------------------------------------------
        !
        IF (PRESENT(use_internal_pbc_corr)) &
            this%use_internal_pbc_corr = use_internal_pbc_corr
        !
        ALLOCATE (this%correction(this%ngm))
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
        CHARACTER(LEN=80) :: sub_name = 'destroy_core_fft'
        !
        !--------------------------------------------------------------------------------
        !
        IF (.NOT. ASSOCIATED(this%cell)) CALL io%destroy_error(sub_name)
        !
        !--------------------------------------------------------------------------------
        !
        IF (this%use_internal_pbc_corr) THEN
            !
            IF (.NOT. ALLOCATED(this%correction)) CALL io%destroy_error(sub_name)
            !
            DEALLOCATE (this%correction)
        END IF
        !
        NULLIFY (this%cell)
        !
        CALL this%deallocate_gvect()
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
        CLASS(core_fft), TARGET, INTENT(IN) :: this
        TYPE(environ_density), INTENT(IN) :: a
        !
        TYPE(environ_gradient), INTENT(INOUT) :: ga
        !
        INTEGER :: ipol
        COMPLEX(DP), DIMENSION(:), ALLOCATABLE :: aux, gaux
        !
        !--------------------------------------------------------------------------------
        !
        ASSOCIATE (dfft => this%cell%dfft, &
                   g => this%g)
            !
            !----------------------------------------------------------------------------
            !
            ALLOCATE (aux(dfft%nnr))
            ALLOCATE (gaux(dfft%nnr))
            !
            aux = CMPLX(a%of_r(:), 0.0_DP, kind=DP)
            !
            CALL env_fwfft(aux, dfft) ! bring a(r) to G-space, a(G)
            !
            !----------------------------------------------------------------------------
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
        END ASSOCIATE
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
        CLASS(core_fft), TARGET, INTENT(IN) :: this
        TYPE(environ_gradient), INTENT(IN) :: ga
        !
        TYPE(environ_density), INTENT(INOUT) :: da
        !
        INTEGER :: n, ipol
        COMPLEX(DP), DIMENSION(:), ALLOCATABLE :: aux, gaux
        COMPLEX(DP) :: fp, fm, aux1, aux2
        !
        !--------------------------------------------------------------------------------
        !
        ASSOCIATE (dfft => this%cell%dfft, &
                   g => this%g)
            !
            !----------------------------------------------------------------------------
            !
            ALLOCATE (aux(dfft%nnr))
            ALLOCATE (gaux(dfft%nnr))
            !
            gaux(:) = (0.0_DP, 0.0_DP)
            !
            IF (dfft%lgamma) THEN
                !
                !------------------------------------------------------------------------
                ! Gamma tricks: perform 2 FFT's in a single shot
                !
                ! x and y
                !
                ipol = 1
                aux(:) = CMPLX(ga%of_r(ipol, :), ga%of_r(ipol + 1, :), kind=DP)
                !
                CALL env_fwfft(aux, dfft) ! bring a(ipol,r) to G-space, a(G)
                !
                !------------------------------------------------------------------------
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
                !------------------------------------------------------------------------
                ! z
                !
                ipol = 3
                aux(:) = CMPLX(ga%of_r(ipol, :), 0.0_DP, kind=DP)
                !
                CALL env_fwfft(aux, dfft) ! bring a(ipol, r) to G-space, a(G)
                !
                !------------------------------------------------------------------------
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
                    !--------------------------------------------------------------------
                    ! Multiply by iG to get the gradient in G-space
                    !
                    DO n = 1, dfft%ngm
                        !
                        gaux(dfft%nl(n)) = &
                            gaux(dfft%nl(n)) + &
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
        END ASSOCIATE
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
        CLASS(core_fft), TARGET, INTENT(IN) :: this
        TYPE(environ_density), INTENT(IN) :: a
        !
        TYPE(environ_density), INTENT(INOUT) :: lapla
        !
        COMPLEX(DP), DIMENSION(:), ALLOCATABLE :: aux, laux
        !
        !--------------------------------------------------------------------------------
        !
        ASSOCIATE (dfft => this%cell%dfft, &
                   gg => this%gg)
            !
            !----------------------------------------------------------------------------
            !
            ALLOCATE (aux(dfft%nnr))
            ALLOCATE (laux(dfft%nnr))
            !
            aux = CMPLX(a%of_r(:), 0.0_DP, kind=DP)
            !
            CALL env_fwfft(aux, dfft) ! bring a(r) to G-space, a(G)
            !
            !----------------------------------------------------------------------------
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
        END ASSOCIATE
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
        CLASS(core_fft), TARGET, INTENT(IN) :: this
        TYPE(environ_density), INTENT(IN) :: a
        !
        TYPE(environ_gradient), INTENT(INOUT) :: ga
        TYPE(environ_hessian), INTENT(INOUT) :: ha
        !
        INTEGER :: ipol, jpol
        COMPLEX(DP), DIMENSION(:), ALLOCATABLE :: aux, gaux, haux
        !
        !--------------------------------------------------------------------------------
        !
        ASSOCIATE (dfft => this%cell%dfft, &
                   g => this%g)
            !
            !----------------------------------------------------------------------------
            !
            ALLOCATE (aux(dfft%nnr))
            ALLOCATE (gaux(dfft%nnr))
            ALLOCATE (haux(dfft%nnr))
            !
            aux = CMPLX(a%of_r(:), 0.0_DP, kind=DP)
            !
            CALL env_fwfft(aux, dfft) ! bring a(r) to G-space, a(G)
            !
            !----------------------------------------------------------------------------
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
                !------------------------------------------------------------------------
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
        END ASSOCIATE
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
        CLASS(core_fft), TARGET, INTENT(IN) :: this
        TYPE(environ_density), INTENT(IN) :: fa, fb
        !
        TYPE(environ_density), INTENT(INOUT) :: fc
        !
        COMPLEX(DP), DIMENSION(:), ALLOCATABLE :: auxr, auxg
        !
        !--------------------------------------------------------------------------------
        !
        ASSOCIATE (dfft => this%cell%dfft, &
                   omega => this%cell%omega)
            !
            !----------------------------------------------------------------------------
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
            !----------------------------------------------------------------------------
            ! Multiply fa(g)*fb(g)
            !
            auxg(dfft%nl(:)) = auxg(dfft%nl(:)) * auxr(dfft%nl(:))
            !
            DEALLOCATE (auxr)
            !
            IF (dfft%lgamma) auxg(dfft%nlm(:)) = &
                CMPLX(REAL(auxg(dfft%nl(:))), -AIMAG(auxg(dfft%nl(:))), kind=DP)
            !
            !----------------------------------------------------------------------------
            ! Brings convolution back to real space
            !
            CALL env_invfft(auxg, dfft)
            !
            fc%of_r(:) = REAL(auxg(:)) * omega
            !
        END ASSOCIATE
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
        CLASS(core_fft), TARGET, INTENT(IN) :: this
        TYPE(environ_density), INTENT(IN) :: fa
        TYPE(environ_gradient), INTENT(IN) :: gb
        !
        TYPE(environ_gradient), INTENT(INOUT) :: gc
        !
        INTEGER :: ipol
        !
        TYPE(environ_density) :: local
        COMPLEX(DP), DIMENSION(:), ALLOCATABLE :: auxr, auxg
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
        ! ASSOCIATE (dfft => this%cell%dfft, &
        !            omega => this%cell%omega)
        !     !
        !     !----------------------------------------------------------------------------
        !     ! Bring fa and fb to reciprocal space
        !     !
        !     ALLOCATE (auxr(dfft%nnr))
        !     auxr(:) = CMPLX(fa%of_r(:), 0.D0, kind=DP)
        !     !
        !     CALL env_fwfft(auxr, dfft)
        !     !
        !     ALLOCATE (auxg(dfft%nnr))
        !     !
        !     DO ipol = 1, 3
        !         auxg(:) = CMPLX(gb%of_r(ipol, :), 0.D0, kind=DP)
        !         !
        !         CALL env_fwfft(auxg, dfft)
        !         !
        !         !------------------------------------------------------------------------
        !         ! Multiply fa(g)*fb(g)
        !         !
        !         auxg(dfft%nl(:)) = auxg(dfft%nl(:)) * auxr(dfft%nl(:))
        !         !
        !         IF (dfft%lgamma) auxg(dfft%nlm(:)) = &
        !             CMPLX(REAL(auxg(dfft%nl(:))), -AIMAG(auxg(dfft%nl(:))), kind=DP)
        !         !
        !         !------------------------------------------------------------------------
        !         ! Brings convolution back to real space
        !         !
        !         CALL env_invfft(auxg, dfft)
        !         !
        !         gc%of_r(ipol, :) = REAL(auxg(:)) * omega
        !     END DO
        !     !
        ! END ASSOCIATE
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
        CLASS(core_fft), TARGET, INTENT(IN) :: this
        TYPE(environ_density), INTENT(IN) :: fa
        TYPE(environ_hessian), INTENT(IN) :: hb
        !
        TYPE(environ_hessian), INTENT(INOUT) :: hc
        !
        COMPLEX(DP), DIMENSION(:), ALLOCATABLE :: auxr, auxg
        !
        INTEGER :: ipol, jpol
        !
        TYPE(environ_density) :: local
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
        ! ASSOCIATE (dfft => this%cell%dfft, &
        !            omega => this%cell%omega)
        !     !
        !     !----------------------------------------------------------------------------
        !     ! Bring fa and fb to reciprocal space
        !     !
        !     ALLOCATE (auxr(dfft%nnr))
        !     auxr(:) = CMPLX(fa%of_r(:), 0.D0, kind=DP)
        !     !
        !     CALL env_fwfft(auxr, dfft)
        !     !
        !     ALLOCATE (auxg(dfft%nnr))
        !     !
        !     DO ipol = 1, 3
        !         !
        !         DO jpol = 1, 3
        !             !
        !             auxg(:) = CMPLX(hb%of_r(ipol, jpol, :), 0.D0, kind=DP)
        !             !
        !             CALL env_fwfft(auxg, dfft)
        !             !
        !             !--------------------------------------------------------------------
        !             ! Multiply fa(g)*fb(g)
        !             !
        !             auxg(dfft%nl(:)) = auxg(dfft%nl(:)) * auxr(dfft%nl(:))
        !             !
        !             IF (dfft%lgamma) auxg(dfft%nlm(:)) = &
        !                 CMPLX(REAL(auxg(dfft%nl(:))), -AIMAG(auxg(dfft%nl(:))), kind=DP)
        !             !
        !             !--------------------------------------------------------------------
        !             ! Brings convolution back to real space
        !             !
        !             CALL env_invfft(auxg, dfft)
        !             !
        !             hc%of_r(ipol, jpol, :) = REAL(auxg(:)) * omega
        !         END DO
        !         !
        !     END DO
        !     !
        ! END ASSOCIATE
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE convolution_hessian
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !                               ELECTROSTATIC METHODS
    !
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !>
    !! Solves the Poisson equation \nabla^2 \phi(r) = - 4 \pi \rho(r)
    !! Input and output functions are defined in real space
    !!
    !! @ param rho : input density
    !! @ param phi : output potential => \phi(r)
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE poisson_fft(this, rho, phi)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(core_fft), TARGET, INTENT(IN) :: this
        TYPE(environ_density), INTENT(IN) :: rho
        !
        TYPE(environ_density), INTENT(INOUT) :: phi
        !
        INTEGER :: ig
        COMPLEX(DP), DIMENSION(:), ALLOCATABLE :: auxr, auxg
        !
        !--------------------------------------------------------------------------------
        !
        ASSOCIATE (gstart => this%gstart, &
                   ngm => this%ngm, &
                   gg => this%gg, &
                   dfft => this%cell%dfft)
            !
            !----------------------------------------------------------------------------
            ! Bring rho to G space
            !
            ALLOCATE (auxr(dfft%nnr))
            auxr = CMPLX(rho%of_r, 0.D0, kind=DP)
            !
            CALL env_fwfft(auxr, dfft)
            !
            ALLOCATE (auxg(ngm))
            auxg = auxr(dfft%nl)
            !
            auxr = CMPLX(0.D0, 0.D0, kind=DP)
            !
            auxr(dfft%nl(1)) = auxg(1) * this%correction(1) * pi ! G = 0 term
            !
!$omp parallel do
            DO ig = gstart, ngm
                auxr(dfft%nl(ig)) = auxg(ig) * (1.D0 / gg(ig) + this%correction(ig) * pi)
            END DO
!$omp end parallel do
            !
            auxr = auxr * e2 / pi
            !
            IF (dfft%lgamma) &
                auxr(dfft%nlm) = CMPLX(REAL(auxr(dfft%nl)), &
                                       -AIMAG(auxr(dfft%nl)), kind=DP)
            !
            !----------------------------------------------------------------------------
            ! Transform hartree potential to real space
            !
            CALL env_invfft(auxr, dfft)
            !
            phi%of_r = DBLE(auxr)
            !
        END ASSOCIATE
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE poisson_fft
    !------------------------------------------------------------------------------------
    !>
    !! Solves the Poisson equation \nabla^2 \phi(r) = - 4 \pi \rho(r)
    !! Input and output functions are defined in real space
    !!
    !! @ param rho  : input density
    !! @ param gphi : output gradient of the potential => \nabla \phi(r)
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE gradpoisson_fft(this, rho, gphi)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(core_fft), TARGET, INTENT(IN) :: this
        TYPE(environ_density), INTENT(IN) :: rho
        !
        TYPE(environ_gradient), INTENT(INOUT) :: gphi
        !
        INTEGER :: ipol, ig
        REAL(DP) :: fac
        COMPLEX(DP), DIMENSION(:), ALLOCATABLE :: auxr, auxg, vaux
        !
        !--------------------------------------------------------------------------------
        !
        ASSOCIATE (ngm => this%ngm, &
                   gstart => this%gstart, &
                   gg => this%gg, &
                   g => this%g, &
                   dfft => this%cell%dfft)
            !
            !----------------------------------------------------------------------------
            ! Bring rho to G space
            !
            ALLOCATE (auxr(dfft%nnr))
            auxr = CMPLX(rho%of_r, 0.D0, KIND=DP)
            !
            CALL env_fwfft(auxr, dfft)
            !
            ALLOCATE (auxg(ngm))
            auxg = auxr(dfft%nl)
            !
            !----------------------------------------------------------------------------
            ! Compute gradient of potential in G space one direction at a time
            !
            DO ipol = 1, 3
                auxr = CMPLX(0.0_DP, 0.0_DP)
                !
!$omp parallel do private(fac)
                DO ig = gstart, ngm
                    !
                    auxr(dfft%nl(ig)) = g(ipol, ig) * &
                                        (1.D0 / gg(ig) + this%correction(ig) * pi) * &
                                        CMPLX(-AIMAG(auxg(ig)), REAL(auxg(ig), kind=DP))
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
                IF (dfft%lgamma) THEN
                    !
                    auxr(dfft%nlm) = &
                        CMPLX(REAL(auxr(dfft%nl)), -AIMAG(auxr(dfft%nl)), kind=DP)
                    !
                END IF
                !
                !------------------------------------------------------------------------
                ! Bring back to R-space, (\grad_ipol a)(r)
                !
                CALL env_invfft(auxr, dfft)
                !
                gphi%of_r(ipol, :) = REAL(auxr)
            END DO
            !
        END ASSOCIATE
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE gradpoisson_fft
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
        INTEGER, INTENT(IN) :: nat
        CLASS(core_fft), TARGET, INTENT(IN) :: this
        TYPE(environ_density), INTENT(IN) :: rho
        CLASS(environ_function), TARGET, INTENT(IN) :: ions(:)
        !
        REAL(DP), INTENT(INOUT) :: force(3, nat)
        !
        INTEGER :: iat, ig, ityp
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
        REAL(DP), POINTER :: Z ! ionic charge
        REAL(DP), POINTER :: D ! gaussian spread of smeared ion function
        !
        TYPE(environ_function_gaussian), POINTER :: local_ions(:)
        !
        CHARACTER(LEN=80) :: sub_name = 'force_fft'
        !
        !--------------------------------------------------------------------------------
        !
        SELECT TYPE (ions)
            !
        TYPE IS (environ_function_gaussian)
            local_ions => ions
            !
        CLASS DEFAULT
            CALL io%error(sub_name, "Unexpected function type", 1)
            !
        END SELECT
        !
        IF (nat /= SIZE(local_ions)) &
            CALL io%error(sub_name, &
                          'Mismatch in numbers of atoms passed in input and stored', 1)
        !
        !--------------------------------------------------------------------------------
        !
        ASSOCIATE (ngm => this%ngm, &
                   gstart => this%gstart, &
                   G => this%g, &
                   G2 => this%gg, &
                   dfft => this%cell%dfft)
            !
            ALLOCATE (auxr(dfft%nnr))
            !
            !----------------------------------------------------------------------------
            ! Bring rho to G space
            !
            auxr = CMPLX(rho%of_r, 0.D0, KIND=DP)
            !
            CALL env_fwfft(auxr, dfft)
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
            DO iat = 1, nat
                Z => local_ions(iat)%volume
                D => local_ions(iat)%spread
                R = local_ions(iat)%pos - this%cell%origin ! account for any origin shift
                !
                DO ig = gstart, ngm
                    fpibg2 = fpi / (G2(ig) * tpi2)
                    !
                    t_arg = tpi * SUM(G(:, ig) * R)
                    !
                    euler_term = SIN(t_arg) * DBLE(auxg(ig)) + &
                                 COS(t_arg) * AIMAG(auxg(ig))
                    !
                    e_arg = -0.25D0 * D**2 * G2(ig) * tpi2
                    gauss_term = fpibg2 * EXP(e_arg)
                    !
                    main_term = gauss_term + this%correction(ig)
                    !
                    force(:, iat) = force(:, iat) + G(:, ig) * euler_term * main_term
                END DO
                !
                force(:, iat) = force(:, iat) * Z
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
    !! Gradient of Hartree potential in R space from a total
    !! (spinless) density in R space n(r)
    !! wg_corr_h => calc_vmt
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE hessv_h_of_rho_r(this, rho, hessv)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(core_fft), INTENT(IN), TARGET :: this
        REAL(DP), INTENT(IN) :: rho(this%cell%dfft%nnr)
        !
        REAL(DP), INTENT(OUT) :: hessv(3, 3, this%cell%dfft%nnr)
        !
        COMPLEX(DP), DIMENSION(:), ALLOCATABLE :: auxr, auxg
        INTEGER :: ig, ipol, jpol
        LOGICAL :: gamma_only = .TRUE.
        !
        !--------------------------------------------------------------------------------
        !
        ASSOCIATE (ngm => this%ngm, &
                   gstart => this%gstart, &
                   gg => this%gg, &
                   g => this%g, &
                   dfft => this%cell%dfft)
            !
            !----------------------------------------------------------------------------
            ! Bring rho to G space
            !
            ALLOCATE (auxr(dfft%nnr))
            auxr = CMPLX(rho, 0.D0, KIND=DP)
            !
            CALL env_fwfft(auxr, dfft)
            !
            !----------------------------------------------------------------------------
            ! Compute total potential in G space
            !
            ALLOCATE (auxg(dfft%nnr))
            !
            DO ipol = 1, 3
                !
                DO jpol = 1, 3
                    auxg = (0.0_DP, 0.0_DP)
                    !
                    DO ig = gstart, ngm
                        !
                        auxg(dfft%nl(ig)) = &
                            g(ipol, ig) * g(jpol, ig) * &
                            (1.D0 / gg(ig) + this%correction(ig) * tpi2) * &
                            CMPLX(DBLE(auxr(dfft%nl(ig))), &
                                  AIMAG(auxr(dfft%nl(ig))), kind=DP)
                        !
                    END DO
                    !
                    !--------------------------------------------------------------------
                    ! Add the factor e2*fpi coming from the missing prefactor of
                    ! V = e2 * fpi
                    !
                    auxg = auxg * e2 * fpi
                    !
                    IF (gamma_only) THEN
                        !
                        auxg(dfft%nlm) = &
                            CMPLX(DBLE(auxg(dfft%nl)), -AIMAG(auxg(dfft%nl)), kind=DP)
                        !
                    END IF
                    !
                    CALL env_invfft(auxg, dfft) ! bring back to R-space
                    !
                    hessv(ipol, jpol, :) = REAL(auxg)
                    !
                END DO
                !
            END DO
            !
        END ASSOCIATE
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE hessv_h_of_rho_r
    !------------------------------------------------------------------------------------
    !>
    !! Gradient of Hartree potential in R space from a total
    !! (spinless) density in R space n(r)
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE field_of_gradrho(this, gradrho, e)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(core_fft), INTENT(IN), TARGET :: this
        REAL(DP), INTENT(IN) :: gradrho(3, this%cell%dfft%nnr)
        !
        REAL(DP), INTENT(OUT) :: e(this%cell%dfft%nnr)
        !
        COMPLEX(DP), DIMENSION(:), ALLOCATABLE :: auxr, auxg, auxe
        !
        INTEGER :: ig, ipol
        LOGICAL :: gamma_only = .TRUE.
        !
        !--------------------------------------------------------------------------------
        !
        ASSOCIATE (ngm => this%ngm, &
                   gstart => this%gstart, &
                   gg => this%gg, &
                   g => this%g, &
                   dfft => this%cell%dfft)
            !
            !----------------------------------------------------------------------------
            ! Bring gradrho to G space
            !
            ALLOCATE (auxe(dfft%nnr))
            auxe = CMPLX(0.D0, 0.D0, KIND=DP)
            !
            ALLOCATE (auxr(dfft%nnr))
            auxr = CMPLX(0.D0, 0.D0, KIND=DP)
            !
            ALLOCATE (auxg(dfft%nnr))
            !
            DO ipol = 1, 3
                auxg = CMPLX(gradrho(ipol, :), 0.D0, KIND=DP)
                !
                CALL env_fwfft(auxg, dfft)
                !
                !------------------------------------------------------------------------
                ! Compute total potential in G space
                !
                DO ig = gstart, ngm
                    !
                    auxr(dfft%nl(ig)) = &
                        g(ipol, ig) * (1.D0 / gg(ig) + this%correction(ig) * pi) * &
                        CMPLX(-AIMAG(auxg(dfft%nl(ig))), &
                              REAL(auxg(dfft%nl(ig))), kind=DP)
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
            DEALLOCATE (auxg)
            DEALLOCATE (auxr)
            !
            IF (gamma_only) &
                auxe(dfft%nlm) = CMPLX(REAL(auxe(dfft%nl)), &
                                       -AIMAG(auxe(dfft%nl)), kind=DP)
            !
            CALL env_invfft(auxe, dfft)
            ! bring back to R-space (\grad_ipol a)(r)
            !
            e = REAL(auxe)
            !
        END ASSOCIATE
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE field_of_gradrho
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !                               PRIVATE HELPER METHODS
    !
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !>
    !! This routine generates all the reciprocal lattice vectors
    !! contained in the sphere of radius gcutm. Furthermore it
    !! computes the indices nl which give the correspondence
    !! between the fft mesh points and the array of g vectors.
    !
    !------------------------------------------------------------------------------------
    SUBROUTINE env_ggen(dfftp, comm, at, bg, gcutm, ngm_g, ngm, g, gg, gstart, &
                        no_global_sort)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP), INTENT(IN) :: at(3, 3), bg(3, 3), gcutm
        INTEGER, INTENT(IN) :: ngm_g, comm
        LOGICAL, OPTIONAL, INTENT(IN) :: no_global_sort
        ! if no_global_sort is present (and it is true) G vectors are sorted only
        ! locally and not globally. In this case no global array needs to be
        ! allocated and sorted: saves memory and a lot of time for large systems
        !
        TYPE(env_fft_type_descriptor), INTENT(INOUT) :: dfftp
        INTEGER, INTENT(INOUT) :: ngm
        REAL(DP), INTENT(OUT) :: g(:, :), gg(:)
        INTEGER, INTENT(OUT) :: gstart
        !
        REAL(DP) :: tx(3), ty(3), t(3)
        REAL(DP), ALLOCATABLE :: tt(:)
        INTEGER :: ngm_save, n1, n2, n3, ngm_offset, ngm_max, ngm_local
        !
        REAL(DP), ALLOCATABLE :: g2sort_g(:)
        ! array containing only g vectors for the current processor
        !
        INTEGER, ALLOCATABLE :: mill_unsorted(:, :)
        ! array containing all g vectors generators, on all processors
        ! (replicated data). When no_global_sort is present and .true.,
        ! only g-vectors for the current processor are stored
        !
        INTEGER, ALLOCATABLE :: igsrt(:), g2l(:)
        !
        INTEGER :: ni, nj, nk, i, j, k, ipol, ng, igl, indsw
        INTEGER :: istart, jstart, kstart
        INTEGER :: mype, npe
        LOGICAL :: global_sort, is_local
        INTEGER, ALLOCATABLE :: ngmpe(:)
        !
        INTEGER, EXTERNAL :: env_mp_rank, env_mp_size
        !
        CHARACTER(LEN=80) :: sub_name = 'env_ggen'
        !
        !--------------------------------------------------------------------------------
        !
        global_sort = .TRUE.
        !
        IF (PRESENT(no_global_sort)) global_sort = .NOT. no_global_sort
        !
        IF (.NOT. global_sort) THEN
            ngm_max = ngm
        ELSE
            ngm_max = ngm_g
        END IF
        !
        ngm_save = ngm ! save current value of ngm
        !
        ngm = 0
        ngm_local = 0
        !
        gg(:) = gcutm + 1.D0
        ! set the total number of fft mesh points and and initial value of gg
        ! The choice of gcutm is due to the fact that we have to order the
        ! vectors after computing them.
        !
        !--------------------------------------------------------------------------------
        ! Computes all the g vectors inside a sphere
        !
        ALLOCATE (mill_unsorted(3, ngm_save))
        ALLOCATE (igsrt(ngm_max))
        ALLOCATE (g2l(ngm_max))
        ALLOCATE (g2sort_g(ngm_max))
        !
        g2sort_g(:) = 1.0D20
        !
        !--------------------------------------------------------------------------------
        !
        ALLOCATE (tt(dfftp%nr3)) ! allocate temporal array
        !
        !--------------------------------------------------------------------------------
        ! Max miller indices (same convention as in module stick_set)
        !
        ni = (dfftp%nr1 - 1) / 2
        nj = (dfftp%nr2 - 1) / 2
        nk = (dfftp%nr3 - 1) / 2
        !
        !--------------------------------------------------------------------------------
        ! Gamma-only: exclude space with x < 0
        !
        istart = 0
        !
        iloop: DO i = istart, ni
            !
            !----------------------------------------------------------------------------
            ! Gamma-only: exclude plane with x = 0, y < 0
            !
            IF (i == 0) THEN
                jstart = 0
            ELSE
                jstart = -nj
            END IF
            !
            tx(1:3) = i * bg(1:3, 1)
            !
            jloop: DO j = jstart, nj
                !
                IF (.NOT. global_sort) THEN
                    !
                    IF (env_fft_stick_index(dfftp, i, j) == 0) CYCLE jloop
                    !
                    is_local = .TRUE.
                ELSE
                    !
                    IF (dfftp%lpara .AND. env_fft_stick_index(dfftp, i, j) == 0) THEN
                        is_local = .FALSE.
                    ELSE
                        is_local = .TRUE.
                    END IF
                    !
                END IF
                !
                !------------------------------------------------------------------------
                ! Gamma-only: exclude line with x = 0, y = 0, z < 0
                !
                IF (i == 0 .AND. j == 0) THEN
                    kstart = 0
                ELSE
                    kstart = -nk
                END IF
                !
                ty(1:3) = tx(1:3) + j * bg(1:3, 2)
                !
                !------------------------------------------------------------------------
                ! Compute all the norm square
                !
                DO k = kstart, nk
                    t(1) = ty(1) + k * bg(1, 3)
                    t(2) = ty(2) + k * bg(2, 3)
                    t(3) = ty(3) + k * bg(3, 3)
                    tt(k - kstart + 1) = t(1)**2 + t(2)**2 + t(3)**2
                END DO
                !
                !------------------------------------------------------------------------
                ! Save all the norm square within cutoff
                !
                DO k = kstart, nk
                    !
                    IF (tt(k - kstart + 1) <= gcutm) THEN
                        ngm = ngm + 1
                        !
                        IF (ngm > ngm_max) &
                            CALL io%error(sub_name, 'Too many g-vectors', ngm)
                        !
                        IF (tt(k - kstart + 1) > eps8) THEN
                            g2sort_g(ngm) = tt(k - kstart + 1)
                        ELSE
                            g2sort_g(ngm) = 0.D0
                        END IF
                        !
                        IF (is_local) THEN
                            ngm_local = ngm_local + 1
                            mill_unsorted(:, ngm_local) = (/i, j, k/)
                            g2l(ngm) = ngm_local
                        ELSE
                            g2l(ngm) = 0
                        END IF
                        !
                    END IF
                    !
                END DO
                !
            END DO jloop
            !
        END DO iloop
        !
        IF (ngm /= ngm_max) &
            CALL io%error(sub_name, 'G-vectors missing!', ABS(ngm - ngm_max))
        !
        igsrt(1) = 0
        !
        IF (.NOT. global_sort) THEN
            CALL env_hpsort_eps(ngm, g2sort_g, igsrt, eps8)
        ELSE
            CALL env_hpsort_eps(ngm_g, g2sort_g, igsrt, eps8)
        END IF
        !
        DEALLOCATE (g2sort_g, tt)
        !
        IF (.NOT. global_sort) THEN
            !
            !----------------------------------------------------------------------------
            ! Compute adequate offsets in order to avoid overlap between
            ! g vectors once they are gathered on a single (global) array
            !
            mype = env_mp_rank(comm)
            npe = env_mp_size(comm)
            ALLOCATE (ngmpe(npe))
            ngmpe = 0
            ngmpe(mype + 1) = ngm
            !
            CALL env_mp_sum(ngmpe, comm)
            !
            ngm_offset = 0
            !
            DO ng = 1, mype
                ngm_offset = ngm_offset + ngmpe(ng)
            END DO
            !
            DEALLOCATE (ngmpe)
            !
        END IF
        !
        ngm = 0
        !
        ngloop: DO ng = 1, ngm_max
            !
            IF (g2l(igsrt(ng)) > 0) THEN
                !
                !------------------------------------------------------------------------
                ! Fetch the indices
                !
                i = mill_unsorted(1, g2l(igsrt(ng)))
                j = mill_unsorted(2, g2l(igsrt(ng)))
                k = mill_unsorted(3, g2l(igsrt(ng)))
                !
                ngm = ngm + 1
                !
                !------------------------------------------------------------------------
                ! Map local and global g index
                ! N.B: the global G vectors arrangement depends on the number of processors
                !
                g(1:3, ngm) = i * bg(:, 1) + j * bg(:, 2) + k * bg(:, 3)
                gg(ngm) = SUM(g(1:3, ngm)**2)
            END IF
            !
        END DO ngloop
        !
        DEALLOCATE (igsrt, g2l)
        !
        IF (ngm /= ngm_save) &
            CALL io%error(sub_name, 'G-vectors (ngm) missing!', ABS(ngm - ngm_save))
        !
        !--------------------------------------------------------------------------------
        ! Determine first nonzero g vector
        !
        IF (gg(1) <= eps8) THEN
            gstart = 2
        ELSE
            gstart = 1
        END IF
        !
        !--------------------------------------------------------------------------------
        !
        CALL env_fft_set_nl(dfftp, at, g)
        ! set nl and nls with the correct fft correspondence
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_ggen
    !------------------------------------------------------------------------------------
    !>
    !! Set local and global dimensions, allocate arrays
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_gvect_init(this, ngm_g, comm)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(core_fft), INTENT(INOUT) :: this
        INTEGER, INTENT(INOUT) :: ngm_g
        !
        INTEGER :: ngm
        !
        INTEGER, INTENT(IN) :: comm
        ! communicator of the group on which g-vecs are distributed
        !
        !--------------------------------------------------------------------------------
        ! Calculate sum over all processors
        !
        ngm = this%ngm
        ngm_g = ngm
        !
        CALL env_mp_sum(ngm_g, comm)
        !
        !--------------------------------------------------------------------------------
        ! Allocate arrays - only those that are always kept until the end
        !
        ALLOCATE (this%gg(ngm))
        ALLOCATE (this%g(3, ngm))
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_gvect_init
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_deallocate_gvect(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(core_fft), INTENT(INOUT) :: this
        !
        !--------------------------------------------------------------------------------
        !
        IF (ALLOCATED(this%gg)) DEALLOCATE (this%gg)
        !
        IF (ALLOCATED(this%g)) DEALLOCATE (this%g)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_deallocate_gvect
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE update_mt_correction(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(core_fft), TARGET, INTENT(INOUT) :: this
        !
        LOGICAL :: physical
        INTEGER :: ir, ig
        REAL(DP) :: r(3), rws, upperbound, ecutrho
        COMPLEX(DP), ALLOCATABLE :: aux(:)
        !
        CHARACTER(LEN=80) :: sub_name = 'update_mt_correction'
        !
        !--------------------------------------------------------------------------------
        !
        ASSOCIATE (alpha => this%alpha, &
                   beta => this%beta, &
                   mt_corr => this%correction, &
                   cell => this%cell, &
                   dfft => this%cell%dfft, &
                   omega => this%cell%omega)
            !
            !----------------------------------------------------------------------------
            !
            ecutrho = this%gcutm * tpi2
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
                IF (alpha <= 0._DP) CALL io%error(sub_name, 'Optimal alpha not found', 1)
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
            DO ir = 1, cell%ir_end
                !
                CALL cell%get_min_distance(ir, 0, 0, cell%origin, r, rws, physical)
                ! compute minimum distance using minimum image convention
                !
                IF (.NOT. physical) CYCLE
                !
                aux(ir) = smooth_coulomb_r(alpha, SQRT(rws))
                !
            END DO
            !
            CALL env_fwfft(aux, dfft)
            !
!$omp parallel do
            DO ig = 1, this%ngm
                !
                mt_corr(ig) = cell%omega * REAL(aux(dfft%nl(ig))) - &
                              smooth_coulomb_g(alpha, beta, tpi2 * this%gg(ig))
                !
            END DO
!$omp end parallel do
            !
            mt_corr = mt_corr * EXP(-tpi2 * this%gg * beta / 4._DP)**2
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
