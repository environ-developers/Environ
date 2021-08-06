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
    USE env_mp, ONLY: env_mp_sum
    !
    USE env_types_fft, ONLY: env_fft_type_descriptor
    USE env_fft_main, ONLY: env_fwfft, env_invfft
    !
    USE environ_param, ONLY: DP, pi, tpi, fpi, e2, eps8
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
        REAL(DP) :: alpha, beta
        REAL(DP), ALLOCATABLE :: mt_corr(:)
        !
        !--------------------------------------------------------------------------------
    CONTAINS
        !--------------------------------------------------------------------------------
        !
        PROCEDURE :: init_first => init_core_fft_electrostatics_first
        PROCEDURE :: init_second => init_core_fft_electrostatics_second
        PROCEDURE :: update_cell => update_core_fft_electrostatics_cell
        PROCEDURE :: destroy => destroy_core_fft_electrostatics
        !
        PROCEDURE :: poisson => poisson_fft
        PROCEDURE :: gradpoisson => gradpoisson_fft
        PROCEDURE :: force => force_fft
        !
        PROCEDURE :: hessv_h_of_rho_r, field_of_gradrho
        !
        PROCEDURE, PRIVATE :: update_mt_correction
        PROCEDURE, PRIVATE :: vmt => calc_vmt
        PROCEDURE, PRIVATE :: gradvmt => calc_gradvmt
        PROCEDURE, PRIVATE :: fmt => calc_fmt
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
    SUBROUTINE init_core_fft_electrostatics_first(this, use_internal_pbc_corr)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        LOGICAL, INTENT(IN), OPTIONAL :: use_internal_pbc_corr
        !
        CLASS(core_fft_electrostatics), INTENT(INOUT) :: this
        !
        !--------------------------------------------------------------------------------
        !
        CALL this%core_fft%init_first()
        !
        IF (PRESENT(use_internal_pbc_corr)) THEN
            this%use_internal_pbc_corr = use_internal_pbc_corr
        ELSE
            this%use_internal_pbc_corr = .FALSE.
        END IF
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_core_fft_electrostatics_first
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_core_fft_electrostatics_second(this, gcutm, cell)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP), INTENT(IN) :: gcutm
        TYPE(environ_cell), TARGET, INTENT(IN) :: cell
        !
        CLASS(core_fft_electrostatics), INTENT(INOUT) :: this
        !
        INTEGER :: ngm_g ! global number of G vectors (summed on all procs)
        ! in serial execution, ngm_g = ngm
        !
        !--------------------------------------------------------------------------------
        !
        CALL this%core_fft%init_second(gcutm, cell)
        !
        IF (this%use_internal_pbc_corr) THEN
            ALLOCATE (this%mt_corr(this%ngm))
            !
            CALL this%update_mt_correction()
            !
        END IF
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_core_fft_electrostatics_second
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
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE update_core_fft_electrostatics_cell
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE destroy_core_fft_electrostatics(this, lflag)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        LOGICAL, INTENT(IN) :: lflag
        !
        CLASS(core_fft_electrostatics), INTENT(INOUT) :: this
        !
        CHARACTER(LEN=80) :: sub_name = 'destroy_core_fft_electrostatics'
        !
        !--------------------------------------------------------------------------------
        !
        CALL this%core_fft%destroy(lflag)
        !
        IF (this%use_internal_pbc_corr) DEALLOCATE (this%mt_corr)
        !
        RETURN
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
    !! Solves the Poisson equation \nabla \fout(r) = - 4 \pi \fin(r)
    !! Input and output functions are defined in real space
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE poisson_fft(this, fin, fout)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(core_fft_electrostatics), TARGET, INTENT(IN) :: this
        TYPE(environ_density), INTENT(IN) :: fin
        !
        TYPE(environ_density), INTENT(INOUT) :: fout
        !
        INTEGER :: ig
        COMPLEX(DP), DIMENSION(:), ALLOCATABLE :: auxr, auxg, vaux
        !
        INTEGER, POINTER :: gstart, ngm
        REAL(DP), POINTER :: tpiba2
        REAL(DP), POINTER :: gg(:)
        TYPE(env_fft_type_descriptor), POINTER :: dfft
        !
        !--------------------------------------------------------------------------------
        !
        tpiba2 => this%cell%tpiba2
        gstart => this%gstart
        ngm => this%ngm
        gg => this%gg
        dfft => this%cell%dfft
        !
        !--------------------------------------------------------------------------------
        ! Bring fin%of_r from R-space --> G-space
        !
        ALLOCATE (auxr(dfft%nnr))
        auxr = CMPLX(fin%of_r, 0.D0, kind=DP)
        !
        CALL env_fwfft(auxr, dfft)
        !
        ALLOCATE (auxg(ngm))
        auxg = auxr(dfft%nl(:))
        !
        auxr = CMPLX(0.D0, 0.D0, kind=DP)
        !
!$omp parallel do
        DO ig = gstart, ngm
            auxr(dfft%nl(ig)) = auxg(ig) / gg(ig)
        END DO
!$omp end parallel do
        !
        auxr = auxr * e2 * fpi / tpiba2
        !
        IF (this%use_internal_pbc_corr) THEN
            ALLOCATE (vaux(ngm))
            !
            CALL this%vmt(auxg, vaux)
            !
            auxr(dfft%nl(1:ngm)) = auxr(dfft%nl(1:ngm)) + vaux(1:ngm)
            DEALLOCATE (vaux)
        END IF
        !
        IF (dfft%lgamma) THEN
            !
            auxr(dfft%nlm(:)) = CMPLX(REAL(auxr(dfft%nl(:))), &
                                      -AIMAG(auxr(dfft%nl(:))), kind=DP)
            !
        END IF
        !
        !--------------------------------------------------------------------------------
        ! Transform hartree potential to real space
        !
        CALL env_invfft(auxr, dfft)
        !
        fout%of_r(:) = DBLE(auxr(:))
        DEALLOCATE (auxr)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE poisson_fft
    !------------------------------------------------------------------------------------
    !>
    !! Solves the Poisson equation \grad \gout(r) = - 4 \pi \fin(r)
    !! where gout is the gradient of the potential
    !! Input and output functions are defined in real space
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE gradpoisson_fft(this, fin, gout)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(core_fft_electrostatics), TARGET, INTENT(IN) :: this
        TYPE(environ_density), INTENT(IN) :: fin
        !
        TYPE(environ_gradient), INTENT(INOUT) :: gout
        !
        INTEGER :: ipol, ig
        REAL(DP) :: fac
        COMPLEX(DP), DIMENSION(:), ALLOCATABLE :: auxr, auxg, vaux
        !
        INTEGER, POINTER :: ngm, gstart
        REAL(DP), POINTER :: tpiba
        REAL(DP), POINTER :: gg(:)
        REAL(DP), POINTER :: g(:, :)
        TYPE(env_fft_type_descriptor), POINTER :: dfft
        !
        !--------------------------------------------------------------------------------
        !
        tpiba => this%cell%tpiba
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
        auxr(:) = CMPLX(fin%of_r(:), 0.D0, KIND=DP)
        !
        CALL env_fwfft(auxr, dfft)
        !
        ALLOCATE (auxg(ngm))
        auxg = auxr(dfft%nl(:))
        !
        !--------------------------------------------------------------------------------
        ! Compute gradient of potential in G space one direction at a time
        !
        DO ipol = 1, 3
            auxr(:) = CMPLX(0.0_DP, 0.0_DP)
            !
!$omp parallel do private(fac)
            DO ig = gstart, ngm
                !
                auxr(dfft%nl(ig)) = &
                    CMPLX(-AIMAG(auxg(ig)), &
                          REAL(auxg(ig), kind=DP)) * g(ipol, ig) / gg(ig)
                !
            END DO
!$omp end parallel do
            !
            !----------------------------------------------------------------------------
            ! Add the factor e2*fpi/2\pi/a coming from the missing prefactor of
            ! V = e2 * fpi divided by the 2\pi/a factor missing in G
            !
            fac = e2 * fpi / tpiba
            auxr = auxr * fac
            !
            !----------------------------------------------------------------------------
            ! Add martyna-tuckerman correction, if needed
            !
            IF (this%use_internal_pbc_corr) THEN
                ALLOCATE (vaux(ngm))
                !
                CALL this%gradvmt(ipol, auxg, vaux)
                !
                auxr(dfft%nl(:)) = auxr(dfft%nl(:)) + vaux(:)
                DEALLOCATE (vaux)
            END IF
            !
            !----------------------------------------------------------------------------
            ! Assuming GAMMA ONLY
            !
            IF (dfft%lgamma) THEN
                !
                auxr(dfft%nlm(:)) = &
                    CMPLX(REAL(auxr(dfft%nl(:))), -AIMAG(auxr(dfft%nl(:))), kind=DP)
                !
            END IF
            !
            !----------------------------------------------------------------------------
            ! Bring back to R-space, (\grad_ipol a)(r) ...
            !
            CALL env_invfft(auxr, dfft)
            !
            gout%of_r(ipol, :) = REAL(auxr(:))
        END DO
        !
        DEALLOCATE (auxr)
        DEALLOCATE (auxg)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE gradpoisson_fft
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE force_fft(this, rho, ions, nat, force)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: nat
        CLASS(core_fft_electrostatics), TARGET, INTENT(IN) :: this
        TYPE(environ_density), INTENT(IN) :: rho
        TYPE(environ_ions), TARGET, INTENT(IN) :: ions
        !
        REAL(DP), INTENT(OUT) :: force(3, nat)
        !
        INTEGER :: iat, ig, ityp
        REAL(DP) :: fpibg2, e_arg, gauss_term, t_arg, euler_term, fact
        COMPLEX(DP), DIMENSION(:), ALLOCATABLE :: auxr, auxg
        REAL(DP), ALLOCATABLE :: ftmp(:, :)
        !
        INTEGER, POINTER :: ngm, gstart
        REAL(DP), POINTER :: tpiba, tpiba2, Z, D, R(:)
        REAL(DP), POINTER :: G(:, :), G2(:)
        TYPE(env_fft_type_descriptor), POINTER :: dfft
        !
        !--------------------------------------------------------------------------------
        !
        tpiba => this%cell%tpiba
        tpiba2 => this%cell%tpiba2
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
        auxg = auxr(dfft%nl(:)) ! aux now contains n(G)
        DEALLOCATE (auxr)
        !
        !--------------------------------------------------------------------------------
        ! Compute forces
        !
        IF (dfft%lgamma) THEN
            fact = 2.D0
        ELSE
            fact = 1.D0
        END IF
        !
        force = 0.D0
        !
        DO iat = 1, nat
            D => ions%iontype(ions%ityp(iat))%atomicspread
            Z => ions%iontype(ions%ityp(iat))%zv
            R => ions%tau(:, iat)
            !
            DO ig = gstart, ngm
                !
                IF (G2(ig) <= eps8) CYCLE
                !
                fpibg2 = fpi / (G2(ig) * tpiba2)
                !
                e_arg = -0.25D0 * D**2 * G2(ig) * tpiba2
                gauss_term = Z * fpibg2 * EXP(e_arg)
                !
                t_arg = tpi * SUM(G(:, ig) * R)
                euler_term = SIN(t_arg) * DBLE(auxg(ig)) + COS(t_arg) * AIMAG(auxg(ig))
                !
                force(:, iat) = force(:, iat) + G(:, ig) * gauss_term * euler_term
            END DO
            !
        END DO
        !
        force = e2 * fact * force * tpiba
        !
        !--------------------------------------------------------------------------------
        ! DEBUGGING
        !
        ! WRITE (*, *) 'forces lc' !
        ! !
        ! DO iat = 1, nat
        !     WRITE (*, '(i8,3f10.4)') iat, force(:, iat)
        ! END DO
        !
        !--------------------------------------------------------------------------------
        ! Add martyna-tuckerman correction, if needed
        !
        IF (this%use_internal_pbc_corr) THEN
            ALLOCATE (ftmp(3, nat))
            !
            CALL this%fmt(auxg, ions, ftmp)
            !
            force = force + fact * ftmp
            !
            !----------------------------------------------------------------------------
            ! DEBUGGING
            !
            ! WRITE (*, *) 'forces mt', fact
            ! !
            ! DO iat = 1, nat
            !     WRITE (*, '(i8,3f10.4)') iat, fact * ftmp(:, iat)
            ! END DO
            !
            DEALLOCATE (ftmp)
        END IF
        !
        CALL env_mp_sum(force, rho%cell%dfft%comm)
        !
        DEALLOCATE (auxg)
        !
        RETURN
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
        REAL(DP), POINTER :: tpiba2
        !
        INTEGER, POINTER :: ngm, gstart
        LOGICAL, POINTER :: do_comp_mt
        REAL(DP), POINTER :: g(:, :), gg(:)
        !
        COMPLEX(DP), DIMENSION(:), ALLOCATABLE :: rhoaux, gaux, rgtot, vaux
        REAL(DP) :: fac, eh_corr
        INTEGER :: ig, ipol, jpol
        LOGICAL :: gamma_only
        !
        !--------------------------------------------------------------------------------
        !
        tpiba2 => this%cell%tpiba2
        ngm => this%ngm
        gstart => this%gstart
        gg => this%gg
        g => this%g
        do_comp_mt => this%use_internal_pbc_corr
        gamma_only = .TRUE.
        !
        !--------------------------------------------------------------------------------
        ! Bring rho to G space
        !
        ALLOCATE (rhoaux(this%cell%dfft%nnr))
        rhoaux(:) = CMPLX(rho(:), 0.D0, KIND=DP)
        !
        CALL env_fwfft(rhoaux, this%cell%dfft)
        !
        !--------------------------------------------------------------------------------
        ! Compute total potential in G space
        !
        ALLOCATE (gaux(this%cell%dfft%nnr))
        !
        DO ipol = 1, 3
            !
            DO jpol = 1, 3
                gaux(:) = (0.0_DP, 0.0_DP)
                !
                DO ig = gstart, ngm
                    fac = g(ipol, ig) * g(jpol, ig) / gg(ig)
                    !
                    gaux(this%cell%dfft%nl(ig)) = &
                        CMPLX(REAL(rhoaux(this%cell%dfft%nl(ig))), &
                              AIMAG(rhoaux(this%cell%dfft%nl(ig))), kind=DP) * fac
                    !
                END DO
                !
                !------------------------------------------------------------------------
                ! Add the factor e2*fpi coming from the missing prefactor of
                ! V = e2 * fpi
                !
                fac = e2 * fpi
                gaux = gaux * fac
                !
                !------------------------------------------------------------------------
                ! Add martyna-tuckerman correction, if needed
                !
                IF (do_comp_mt) THEN
                    ALLOCATE (vaux(ngm), rgtot(ngm))
                    rgtot(1:ngm) = rhoaux(this%cell%dfft%nl(1:ngm))
                    !
                    CALL this%vmt(rgtot, vaux)
                    !
                    DO ig = gstart, ngm
                        fac = g(ipol, ig) * g(jpol, ig) * tpiba2
                        !
                        gaux(this%cell%dfft%nl(ig)) = &
                            gaux(this%cell%dfft%nl(ig)) + &
                            CMPLX(REAL(vaux(ig)), AIMAG(vaux(ig)), kind=DP) * fac
                        !
                    END DO
                    !
                    DEALLOCATE (rgtot, vaux)
                END IF
                !
                IF (gamma_only) THEN
                    !
                    gaux(this%cell%dfft%nlm(:)) = &
                        CMPLX(REAL(gaux(this%cell%dfft%nl(:))), &
                              -AIMAG(gaux(this%cell%dfft%nl(:))), kind=DP)
                    !
                END IF
                !
                CALL env_invfft(gaux, this%cell%dfft) ! bring back to R-space
                !
                hessv(ipol, jpol, :) = REAL(gaux(:))
                !
            END DO
            !
        END DO
        !
        DEALLOCATE (gaux)
        DEALLOCATE (rhoaux)
        !
        RETURN
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
        REAL(DP), POINTER :: tpiba
        !
        INTEGER, POINTER :: ngm, gstart
        LOGICAL, POINTER :: do_comp_mt
        REAL(DP), POINTER :: g(:, :), gg(:)
        !
        COMPLEX(DP), DIMENSION(:), ALLOCATABLE :: aux, eaux, gaux, rgtot, vaux
        !
        REAL(DP) :: fac, eh_corr
        INTEGER :: ig, ipol
        LOGICAL :: gamma_only
        !
        !--------------------------------------------------------------------------------
        !
        tpiba => this%cell%tpiba
        ngm => this%ngm
        gstart => this%gstart
        gg => this%gg
        g => this%g
        do_comp_mt => this%use_internal_pbc_corr
        gamma_only = .TRUE.
        !
        !--------------------------------------------------------------------------------
        ! Bring gradrho to G space
        !
        ALLOCATE (eaux(this%cell%dfft%nnr))
        eaux(:) = CMPLX(0.D0, 0.D0, KIND=DP)
        !
        ALLOCATE (aux(this%cell%dfft%nnr))
        aux(:) = CMPLX(0.D0, 0.D0, KIND=DP)
        !
        ALLOCATE (gaux(this%cell%dfft%nnr))
        !
        IF (do_comp_mt) ALLOCATE (vaux(ngm), rgtot(ngm))
        !
        DO ipol = 1, 3
            gaux(:) = CMPLX(gradrho(ipol, :), 0.D0, KIND=DP)
            !
            CALL env_fwfft(gaux, this%cell%dfft)
            !
            !----------------------------------------------------------------------------
            ! Compute total potential in G space
            !
            DO ig = gstart, ngm
                fac = g(ipol, ig) / gg(ig)
                !
                aux(this%cell%dfft%nl(ig)) = &
                    CMPLX(-AIMAG(gaux(this%cell%dfft%nl(ig))), &
                          REAL(gaux(this%cell%dfft%nl(ig))), kind=DP) * fac
                !
            END DO
            !
            !----------------------------------------------------------------------------
            ! Add the factor e2*fpi/2\pi/a coming from the missing prefactor of
            ! V = e2 * fpi divided by the 2\pi/a factor missing in G
            !
            fac = e2 * fpi / tpiba
            aux = aux * fac
            !
            !----------------------------------------------------------------------------
            ! Add martyna-tuckerman correction, if needed
            !
            IF (do_comp_mt) THEN
                rgtot(1:ngm) = gaux(this%cell%dfft%nl(1:ngm))
                !
                CALL this%vmt(rgtot, vaux)
                !
                DO ig = gstart, ngm
                    fac = g(ipol, ig) * tpiba
                    !
                    aux(this%cell%dfft%nl(ig)) = &
                        aux(this%cell%dfft%nl(ig)) + &
                        CMPLX(-AIMAG(vaux(ig)), REAL(vaux(ig)), kind=DP) * fac
                    !
                END DO
                !
            END IF
            !
            eaux = eaux + aux
        END DO
        !
        IF (do_comp_mt) DEALLOCATE (rgtot, vaux)
        !
        DEALLOCATE (gaux)
        DEALLOCATE (aux)
        !
        IF (gamma_only) THEN
            !
            eaux(this%cell%dfft%nlm(:)) = &
                CMPLX(REAL(eaux(this%cell%dfft%nl(:))), &
                      -AIMAG(eaux(this%cell%dfft%nl(:))), kind=DP)
            !
        END IF
        !
        CALL env_invfft(eaux, this%cell%dfft)
        ! bring back to R-space (\grad_ipol a)(r)
        !
        e(:) = REAL(eaux(:))
        !
        DEALLOCATE (eaux)
        !
        RETURN
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
        REAL(DP), POINTER :: alpha, beta, omega, tpiba2
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
        mt_corr => this%mt_corr
        !
        cell => this%cell
        dfft => this%cell%dfft
        omega => this%cell%omega
        tpiba2 => this%cell%tpiba2
        !
        ecutrho = this%gcutm * tpiba2
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
            IF (alpha <= 0._DP) CALL env_errore(sub_name, 'Optimal alpha not found', 1)
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
            aux(ir) = smooth_coulomb_r(alpha, SQRT(rws) * cell%alat)
            !
        END DO
        !
        CALL env_fwfft(aux, dfft)
        !
!$omp parallel do
        DO ig = 1, this%ngm
            !
            mt_corr(ig) = cell%omega * REAL(aux(dfft%nl(ig))) - &
                          smooth_coulomb_g(alpha, beta, cell%tpiba2 * this%gg(ig))
            !
        END DO
!$omp end parallel do
        !
        mt_corr(:) = mt_corr(:) * EXP(-cell%tpiba2 * this%gg(:) * beta / 4._DP)**2
        !
        DEALLOCATE (aux)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE update_mt_correction
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_vmt(this, rho, v)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(core_fft_electrostatics), INTENT(IN) :: this
        COMPLEX(DP), INTENT(IN) :: rho(this%ngm)
        !
        COMPLEX(DP), INTENT(OUT) :: v(this%ngm)
        !
        INTEGER :: ig
        !
        !--------------------------------------------------------------------------------
        !
        v(:) = (0._DP, 0._DP)
        !
!$omp parallel do
        DO ig = 1, this%ngm
            v(ig) = this%mt_corr(ig) * rho(ig)
        END DO
!$omp end parallel do
        !
        v = e2 * v
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_vmt
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_gradvmt(this, ipol, rho, v)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: ipol
        CLASS(core_fft_electrostatics), INTENT(IN) :: this
        COMPLEX(DP), INTENT(IN) :: rho(this%ngm)
        !
        COMPLEX(DP), INTENT(OUT) :: v(this%ngm)
        !
        INTEGER :: ig
        REAL(DP) :: fac
        !
        !--------------------------------------------------------------------------------
        !
        v(:) = (0._DP, 0._DP)
        !
!$omp parallel do private(fac)
        DO ig = this%gstart, this%ngm
            fac = this%g(ipol, ig) * this%cell%tpiba
            !
            v(ig) = this%mt_corr(ig) * &
                    CMPLX(-AIMAG(rho(ig)), REAL(rho(ig)), kind=dp) * fac
            !
        END DO
!$omp end parallel do
        !
        v = e2 * v
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_gradvmt
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_fmt(this, rho, ions, force)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CLASS(core_fft_electrostatics), INTENT(IN) :: this
        TYPE(environ_ions), INTENT(IN) :: ions
        COMPLEX(DP), INTENT(IN) :: rho(this%ngm)
        !
        REAL(DP), INTENT(OUT) :: force(3, ions%number)
        !
        INTEGER :: iat, ig
        REAL(DP) :: arg
        !
        !--------------------------------------------------------------------------------
        !
        force = 0._DP
        !
        DO iat = 1, ions%number
            !
            DO ig = 1, this%ngm
                arg = tpi * SUM(this%g(:, ig) * ions%tau(:, iat))
                !
                force(:, iat) = &
                    force(:, iat) + this%g(:, ig) * &
                    (SIN(arg) * REAL(rho(ig)) + COS(arg) * AIMAG(rho(ig))) * &
                    this%mt_corr(ig)
                !
            END DO
            !
            force(:, iat) = &
                force(:, iat) * ions%iontype(ions%ityp(iat))%zv * this%cell%tpiba
            !
        END DO
        !
        force = e2 * force
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_fmt
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
