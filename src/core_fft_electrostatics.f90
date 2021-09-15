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
MODULE class_core_fft_electrostatics
    !------------------------------------------------------------------------------------
    !
    USE class_io, ONLY: io
    USE env_mp, ONLY: env_mp_sum
    !
    USE env_types_fft, ONLY: env_fft_type_descriptor
    USE env_fft_main, ONLY: env_fwfft, env_invfft
    !
    USE environ_param, ONLY: DP, pi, tpi, tpi2, fpi, e2, eps8
    !
    USE class_cell
    USE class_density
    USE class_gradient
    !
    USE class_core_fft_derivatives
    !
    USE class_ions
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
    TYPE, EXTENDS(core_fft_derivatives), PUBLIC :: core_fft_electrostatics
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
        PROCEDURE :: init => init_core_fft_electrostatics
        PROCEDURE :: update_cell => update_core_fft_electrostatics_cell
        PROCEDURE :: destroy => destroy_core_fft_electrostatics
        !
        PROCEDURE, PRIVATE :: update_mt_correction
        !
        PROCEDURE :: poisson => poisson_fft
        PROCEDURE :: gradpoisson => gradpoisson_fft
        PROCEDURE :: force => force_fft
        !
        PROCEDURE :: hessv_h_of_rho_r, field_of_gradrho
        !
        !--------------------------------------------------------------------------------
    END TYPE core_fft_electrostatics
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
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_core_fft_electrostatics(this, gcutm, cell, use_internal_pbc_corr)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP), INTENT(IN) :: gcutm
        TYPE(environ_cell), TARGET, INTENT(IN) :: cell
        LOGICAL, INTENT(IN), OPTIONAL :: use_internal_pbc_corr
        !
        CLASS(core_fft_electrostatics), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: sub_name = 'init_core_fft_electrostatics'
        !
        !--------------------------------------------------------------------------------
        !
        IF (ALLOCATED(this%correction)) CALL io%create_error(sub_name)
        !
        !--------------------------------------------------------------------------------
        !
        CALL this%core_fft%init(gcutm, cell)
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
    END SUBROUTINE init_core_fft_electrostatics
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE update_core_fft_electrostatics_cell(this, cell)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(environ_cell), TARGET, INTENT(IN) :: cell
        !
        CLASS(core_fft_electrostatics), INTENT(INOUT) :: this
        !
        !--------------------------------------------------------------------------------
        !
        CALL this%core_fft%update_cell(cell)
        !
        IF (this%use_internal_pbc_corr) CALL this%update_mt_correction()
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE update_core_fft_electrostatics_cell
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE destroy_core_fft_electrostatics(this)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(core_fft_electrostatics), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: sub_name = 'destroy_core_fft_electrostatics'
        !
        !--------------------------------------------------------------------------------
        !
        CALL this%core_fft%destroy()
        !
        IF (this%use_internal_pbc_corr) THEN
            !
            IF (.NOT. ALLOCATED(this%correction)) CALL io%destroy_error(sub_name)
            !
            DEALLOCATE (this%correction)
        END IF
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE destroy_core_fft_electrostatics
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
        CLASS(core_fft_electrostatics), TARGET, INTENT(IN) :: this
        TYPE(environ_density), INTENT(IN) :: rho
        !
        TYPE(environ_density), INTENT(INOUT) :: phi
        !
        INTEGER :: ig
        COMPLEX(DP), DIMENSION(:), ALLOCATABLE :: auxr, auxg
        !
        INTEGER, POINTER :: gstart, ngm
        REAL(DP), POINTER :: gg(:)
        TYPE(env_fft_type_descriptor), POINTER :: dfft
        !
        !--------------------------------------------------------------------------------
        !
        gstart => this%gstart
        ngm => this%ngm
        gg => this%gg
        dfft => this%cell%dfft
        !
        !--------------------------------------------------------------------------------
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
            auxr(dfft%nlm) = CMPLX(REAL(auxr(dfft%nl)), -AIMAG(auxr(dfft%nl)), kind=DP)
        !
        !--------------------------------------------------------------------------------
        ! Transform hartree potential to real space
        !
        CALL env_invfft(auxr, dfft)
        !
        phi%of_r = DBLE(auxr)
        DEALLOCATE (auxr)
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
        CLASS(core_fft_electrostatics), TARGET, INTENT(IN) :: this
        TYPE(environ_density), INTENT(IN) :: rho
        !
        TYPE(environ_gradient), INTENT(INOUT) :: gphi
        !
        INTEGER :: ipol, ig
        REAL(DP) :: fac
        COMPLEX(DP), DIMENSION(:), ALLOCATABLE :: auxr, auxg, vaux
        !
        INTEGER, POINTER :: ngm, gstart
        REAL(DP), POINTER :: gg(:)
        REAL(DP), POINTER :: g(:, :)
        TYPE(env_fft_type_descriptor), POINTER :: dfft
        !
        !--------------------------------------------------------------------------------
        !
        ngm => this%ngm
        gstart => this%gstart
        gg => this%gg
        g => this%g
        dfft => this%cell%dfft
        !
        !--------------------------------------------------------------------------------
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
        !--------------------------------------------------------------------------------
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
            !----------------------------------------------------------------------------
            ! Assuming GAMMA ONLY
            !
            IF (dfft%lgamma) THEN
                !
                auxr(dfft%nlm) = &
                    CMPLX(REAL(auxr(dfft%nl)), -AIMAG(auxr(dfft%nl)), kind=DP)
                !
            END IF
            !
            !----------------------------------------------------------------------------
            ! Bring back to R-space, (\grad_ipol a)(r)
            !
            CALL env_invfft(auxr, dfft)
            !
            gphi%of_r(ipol, :) = REAL(auxr)
        END DO
        !
        DEALLOCATE (auxr)
        DEALLOCATE (auxg)
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
        CLASS(core_fft_electrostatics), TARGET, INTENT(IN) :: this
        TYPE(environ_density), INTENT(IN) :: rho
        TYPE(environ_ions), TARGET, INTENT(IN) :: ions
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
        REAL(DP), POINTER :: G(:, :), G2(:)
        !
        INTEGER, POINTER :: ngm, gstart
        TYPE(env_fft_type_descriptor), POINTER :: dfft
        !
        !--------------------------------------------------------------------------------
        !
        ngm => this%ngm
        gstart => this%gstart
        G => this%g
        G2 => this%gg
        dfft => this%cell%dfft
        !
        ALLOCATE (auxr(dfft%nnr))
        !
        !--------------------------------------------------------------------------------
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
        !--------------------------------------------------------------------------------
        ! Set symmetry factor
        !
        IF (dfft%lgamma) THEN
            fact = 2.D0
        ELSE
            fact = 1.D0
        END IF
        !
        !--------------------------------------------------------------------------------
        ! Compute force
        !
        force = 0.D0
        !
        DO iat = 1, nat
            Z => ions%iontype(ions%ityp(iat))%zv
            D => ions%iontype(ions%ityp(iat))%atomicspread
            R = ions%tau(:, iat) - this%cell%origin ! account for any origin shift
            !
            DO ig = gstart, ngm
                fpibg2 = fpi / (G2(ig) * tpi2)
                !
                t_arg = tpi * SUM(G(:, ig) * R)
                euler_term = SIN(t_arg) * DBLE(auxg(ig)) + COS(t_arg) * AIMAG(auxg(ig))
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
        DEALLOCATE (auxg)
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
        CLASS(core_fft_electrostatics), INTENT(IN), TARGET :: this
        REAL(DP), INTENT(IN) :: rho(this%cell%dfft%nnr)
        !
        REAL(DP), INTENT(OUT) :: hessv(3, 3, this%cell%dfft%nnr)
        !
        INTEGER, POINTER :: ngm, gstart
        REAL(DP), POINTER :: g(:, :), gg(:)
        TYPE(env_fft_type_descriptor), POINTER :: dfft
        !
        COMPLEX(DP), DIMENSION(:), ALLOCATABLE :: auxr, auxg
        INTEGER :: ig, ipol, jpol
        LOGICAL :: gamma_only = .TRUE.
        !
        !--------------------------------------------------------------------------------
        !
        ngm => this%ngm
        gstart => this%gstart
        gg => this%gg
        g => this%g
        dfft => this%cell%dfft
        !
        !--------------------------------------------------------------------------------
        ! Bring rho to G space
        !
        ALLOCATE (auxr(dfft%nnr))
        auxr = CMPLX(rho, 0.D0, KIND=DP)
        !
        CALL env_fwfft(auxr, dfft)
        !
        !--------------------------------------------------------------------------------
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
                        CMPLX(DBLE(auxr(dfft%nl(ig))), AIMAG(auxr(dfft%nl(ig))), kind=DP)
                    !
                END DO
                !
                !------------------------------------------------------------------------
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
        DEALLOCATE (auxg)
        DEALLOCATE (auxr)
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
        CLASS(core_fft_electrostatics), INTENT(IN), TARGET :: this
        REAL(DP), INTENT(IN) :: gradrho(3, this%cell%dfft%nnr)
        !
        REAL(DP), INTENT(OUT) :: e(this%cell%dfft%nnr)
        !
        INTEGER, POINTER :: ngm, gstart
        REAL(DP), POINTER :: g(:, :), gg(:)
        TYPE(env_fft_type_descriptor), POINTER :: dfft
        !
        COMPLEX(DP), DIMENSION(:), ALLOCATABLE :: auxr, auxg, auxe
        !
        INTEGER :: ig, ipol
        LOGICAL :: gamma_only = .TRUE.
        !
        !--------------------------------------------------------------------------------
        !
        ngm => this%ngm
        gstart => this%gstart
        gg => this%gg
        g => this%g
        dfft => this%cell%dfft
        !
        !--------------------------------------------------------------------------------
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
            !----------------------------------------------------------------------------
            ! Compute total potential in G space
            !
            DO ig = gstart, ngm
                !
                auxr(dfft%nl(ig)) = &
                    g(ipol, ig) * (1.D0 / gg(ig) + this%correction(ig) * pi) * &
                    CMPLX(-AIMAG(auxg(dfft%nl(ig))), REAL(auxg(dfft%nl(ig))), kind=DP)
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
            auxe(dfft%nlm) = CMPLX(REAL(auxe(dfft%nl)), -AIMAG(auxe(dfft%nl)), kind=DP)
        !
        CALL env_invfft(auxe, dfft)
        ! bring back to R-space (\grad_ipol a)(r)
        !
        e = REAL(auxe)
        !
        DEALLOCATE (auxe)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE field_of_gradrho
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !                                 CORRECTION METHODS
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
        CLASS(core_fft_electrostatics), TARGET, INTENT(INOUT) :: this
        !
        LOGICAL :: physical
        INTEGER :: ir, ig
        REAL(DP) :: r(3), rws, upperbound, ecutrho
        COMPLEX(DP), ALLOCATABLE :: aux(:)
        !
        REAL(DP), POINTER :: alpha, beta, omega
        REAL(DP), POINTER :: mt_corr(:)
        TYPE(env_fft_type_descriptor), POINTER :: dfft
        TYPE(environ_cell), POINTER :: cell
        !
        CHARACTER(LEN=80) :: sub_name = 'update_mt_correction'
        !
        !--------------------------------------------------------------------------------
        !
        alpha => this%alpha
        beta => this%beta
        mt_corr => this%correction
        !
        cell => this%cell
        dfft => this%cell%dfft
        omega => this%cell%omega
        !
        ecutrho = this%gcutm * tpi2
        !
        !--------------------------------------------------------------------------------
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
        DEALLOCATE (aux)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE update_mt_correction
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
END MODULE class_core_fft_electrostatics
!----------------------------------------------------------------------------------------
