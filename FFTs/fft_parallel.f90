!----------------------------------------------------------------------------------------
!
! Copyright (C) 2021 ENVIRON (www.quantum-environ.org)
! Copyright (C) Quantum ESPRESSO (www.quantum-espresso.org)
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
! Authors: Carlo Cavazzoni (Apr. 2009)
!          Stefano de Gironcoli (Sep-Nov 2016)
!          Modified by Edan Bainglass
!
!----------------------------------------------------------------------------------------
!>
!! Parallel 3D FFT high level Driver (Charge density and Wave Functions)
!!
!! isgn = +-1 : parallel 3d fft for rho and for the potential
!!              NOT IMPLEMENTED WITH TASK GROUPS
!!
!! isgn = +   : G-space to R-space, output = \sum_G f(G)exp(+iG*R)
!!              fft along z using pencils        (env_cft_1z)
!!              transpose across nodes           (env_fft_scatter_yz) ... and reorder for 2d
!!              fft along y using pencils        (cft_1y)
!!              transpose across nodes           (env_fft_scatter_xy) ... and reorder for 2d
!!              fft along x using pencils        (cft_1x)
!!
!! isgn = -   : R-space to G-space, output = \int_R f(R)exp(-iG*R)/Omega
!!              fft along x using pencils        (cft_1x)
!!              transpose across nodes           (env_fft_scatter_xy)
!!              fft along y using pencils        (cft_1y)
!!              transpose across nodes           (env_fft_scatter_yz)
!!              fft along z using pencils        (env_cft_1z)
!!
!! If task_group_fft_is_active the FFT acts on a number of wfcs equal to
!! dfft%nproc2, the number of Y-sections in which a plane is divided.
!! Data are reshuffled by the fft_scatter_tg routine so that each of the
!! dfft%nproc2 subgroups (made by dfft%nproc3 procs) deals with whole planes
!! of a single wavefunciton.
!!
!! The array "planes" signals whether a fft is needed along y :
!! planes(i)=0 : column f(i,*,*) empty , don't do fft along y
!! planes(i)=1 : column f(i,*,*) filled, fft along y needed
!! "empty" = no active components are present in f(i,*,*)
!!           after (isgn>0) or before (isgn<0) the fft on z direction
!!
!! Note that if isgn=+/-1 (fft on rho and pot.) all fft's are needed
!! and all planes(i) are set to 1
!!
!! This driver is based on code written by Stefano de Gironcoli for PWSCF.
!! Task Group added by Costas Bekas, Oct. 2005, adapted from the CPMD code
!! (Alessandro Curioni) and revised by Carlo Cavazzoni 2007.
!!
!----------------------------------------------------------------------------------------
MODULE env_fft_parallel
    !------------------------------------------------------------------------------------
    !
    USE env_fft_param
    !
    USE env_types_fft, ONLY: env_fft_type_descriptor
    !
    USE env_fft_scalar, ONLY: env_cft_1z, env_cft_2xy
    !
    USE env_fft_scatter
    !
    USE omp_lib
    !
    !------------------------------------------------------------------------------------
    !
    PRIVATE
    !
    PUBLIC :: env_tg_cft3s, env_many_cft3s, env_tg_cft3s_2d
    !
    !------------------------------------------------------------------------------------
CONTAINS
    !------------------------------------------------------------------------------------
    !>
    !! General purpose driver
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_tg_cft3s(f, dfft, isgn)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(env_fft_type_descriptor), INTENT(IN) :: dfft ! descriptor of fft data layout
        INTEGER, INTENT(IN) :: isgn ! fft direction
        !
        COMPLEX(DP), INTENT(INOUT) :: f(:) ! array containing data to be transformed
        !
        INTEGER :: n1, n2, n3, nx1, nx2, nx3
        INTEGER :: nnr_
        INTEGER :: nsticks_x, nsticks_y, nsticks_z
        COMPLEX(DP), ALLOCATABLE :: aux(:)
        INTEGER :: i
        !
        CHARACTER(LEN=80) :: sub_name = 'env_tg_cft3s'
        !
        !--------------------------------------------------------------------------------
        !
        n1 = dfft%nr1
        n2 = dfft%nr2
        n3 = dfft%nr3
        nx1 = dfft%nr1x
        nx2 = dfft%nr2x
        nx3 = dfft%nr3x
        !
        nnr_ = dfft%nnr
        nsticks_x = dfft%my_nr2p * dfft%my_nr3p
        nsticks_y = dfft%nr1p(dfft%mype2 + 1) * dfft%my_nr3p
        nsticks_z = dfft%nsp(dfft%mype + 1)
        !
        ALLOCATE (aux(nnr_))
        !
        IF (isgn > 0) THEN ! G -> R
            !
            !$omp parallel do
            DO i = 1, nsticks_z * nx3
                aux(i) = f(i)
            END DO
            !$omp end parallel do
            !
            CALL env_cft_1z(aux, nsticks_z, n3, nx3, isgn, f) ! Rz, Gy, Gx
            !
            CALL env_fft_scatter_yz(dfft, f, aux, nnr_, isgn) ! Gy, Gx, Rz
            !
            CALL env_cft_1z(aux, nsticks_y, n2, nx2, isgn, f) ! Ry, Gx, Rz
            !
            CALL env_fft_scatter_xy(dfft, f, aux, nnr_, isgn) ! Gx, Ry, Rz
            !
            CALL env_cft_1z(aux, nsticks_x, n1, nx1, isgn, f) ! Rx, Ry, Rz
            !
            IF (nsticks_x * nx1 < nnr_) f(nsticks_x * nx1 + 1:nnr_) = (0.0_DP, 0.0_DP)
            ! clean garbage beyond the intended dimension.
            !
        ELSE ! R -> G
            !
            CALL env_cft_1z(f, nsticks_x, n1, nx1, isgn, aux) ! Rx, Ry, Rz
            !
            CALL env_fft_scatter_xy(dfft, f, aux, nnr_, isgn) ! Gx, Ry, Rz
            !
            CALL env_cft_1z(f, nsticks_y, n2, nx2, isgn, aux) ! Ry, Gx, Rz
            !
            CALL env_fft_scatter_yz(dfft, f, aux, nnr_, isgn) ! Gy, Gx, Rz
            !
            CALL env_cft_1z(f, nsticks_z, n3, nx3, isgn, aux) ! Rz, Gy, Gx
            !
            !$omp parallel do
            DO i = 1, nsticks_z * nx3
                f(i) = aux(i)
            END DO
            !$omp end parallel do
            !
        END IF
        !
        DEALLOCATE (aux)
        !
        RETURN
        !
99      FORMAT(20('(', 2F12.9, ')'))
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_tg_cft3s
    !------------------------------------------------------------------------------------
    !>
    !! Specific driver for the new 'many' call
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_many_cft3s(f, dfft, isgn, howmany)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: howmany ! number of FFTs grouped together
        INTEGER, INTENT(IN) :: isgn ! fft direction
        !
        TYPE(env_fft_type_descriptor), INTENT(INOUT), TARGET :: dfft
        ! descriptor of fft data layout
        !
        COMPLEX(DP), INTENT(INOUT) :: f(:) ! array containing data to be transformed
        !
        INTEGER :: n1, n2, n3, nx1, nx2, nx3
        INTEGER :: nnr_
        INTEGER :: nsticks_x, nsticks_y, nsticks_z
        INTEGER :: nsticks_yx, nsticks_zx
        INTEGER :: i, j
        COMPLEX(DP), POINTER :: aux(:)
        !
        CHARACTER(LEN=80) :: sub_name = 'env_many_cft3s'
        !
        !--------------------------------------------------------------------------------
        !
        n1 = dfft%nr1
        n2 = dfft%nr2
        n3 = dfft%nr3
        nx1 = dfft%nr1x
        nx2 = dfft%nr2x
        nx3 = dfft%nr3x
        !
        nnr_ = dfft%nnr
        nsticks_x = dfft%my_nr2p * dfft%my_nr3p
        !
        nsticks_y = dfft%nr1p(dfft%mype2 + 1) * dfft%my_nr3p
        nsticks_yx = MAXVAL(dfft%nr1p) * MAXVAL(dfft%nr3p)
        nsticks_z = dfft%nsp(dfft%mype + 1)
        nsticks_zx = MAXVAL(dfft%nsp)
        !
        aux => dfft%aux
        !
#if defined(__FFT_OPENMP_TASKS)
        CALL execute_using_tasks()
#else
        CALL execute_using_threads()
#endif
        !
        RETURN
        !
99      FORMAT(20('(', 2F12.9, ')'))
        !
        !--------------------------------------------------------------------------------
    CONTAINS
        !--------------------------------------------------------------------------------
#if defined(__FFT_OPENMP_TASKS)
        !>
        !! This implementation requires thread-safe env_cft_1z and multi-threaded MPI
        !!
        !--------------------------------------------------------------------------------
        SUBROUTINE execute_using_tasks()
            !----------------------------------------------------------------------------
            !
            IMPLICIT NONE
            !
            INTEGER :: comm2s
            !
            !----------------------------------------------------------------------------
            !
#if defined(_OPENMP)&&defined(__FFT_SCALAR_THREAD_SAFE)
            CALL env_errore(sub_name//'::execute_using_tasks', &
                            'Needs thread-safe env_fft_scalar backend selected &
                            &at compile time.', 1)
#endif
            !
            IF (isgn > 0) THEN ! G -> R
                !
                !$omp parallel default(none) &
                !$omp&         private(i, j, comm2s)
                !$omp&         shared(howmany, f, nnr_, isgn, dfft)
                !$omp&         shared(nsticks_z, n3, nx3)
                !$omp&         shared(nsticks_y, n2, nx2)
                !$omp&         shared(nsticks_x, n1, nx1)
                !$omp&         shared(nsticks_zx, nsticks_yx, aux)
                !
                !$omp do
                DO i = 0, howmany - 1
                    !
                    CALL env_cft_1z(f(i * nnr_ + 1:), nsticks_z, n3, nx3, isgn, &
                                    aux(nx3 * nsticks_zx * i + 1:))
                    !
                END DO
                !$omp end do
                !
                !$omp single
                !
                CALL env_fft_scatter_many_yz(dfft, aux, f, isgn, howmany)
                !
                DO i = 0, howmany - 1
                    !$omp task depend(out:f(i*nnr_+1:(i+1)*nnr_))
                    !
                    CALL env_cft_1z(f(i * nnr_ + 1:), nsticks_y, n2, nx2, isgn, &
                                    aux(i * nnr_ + 1:), in_place=.TRUE.)
                    !
                    !$omp end task
                    !
                    !$omp task depend(in:f(i*nnr_+1:(i+1)*nnr_)) &
                    !$omp&     depend(out:aux(i*nnr_+1:(i+1)*nnr_))
                    !
                    comm2s = dfft%comm2s(i + 1)
                    !
                    CALL env_fft_scatter_xy(dfft, f(i * nnr_ + 1:), aux(i * nnr_ + 1:), &
                                            nnr_, isgn, comm2s)
                    !
                    !$omp end task
                    !
                    !$omp task depend(in:aux(i*nnr_+1:(i+1)*nnr_))
                    !
                    CALL env_cft_1z(aux(i * nnr_ + 1:), nsticks_x, n1, nx1, isgn, &
                                    f(i * nnr_ + 1:))
                    !
                    !$omp end task
                    !
                END DO
                !
                !$omp end single
                !
                IF (nsticks_x * nx1 < nnr_) THEN
                    !
                    !$omp do collapse(2)
                    DO i = 0, howmany - 1
                        !
                        DO j = nsticks_x * nx1 + 1, nnr_
                            f(j + i * nnr_) = (0.0_DP, 0.0_DP)
                        END DO
                        !
                    END DO
                    !$omp end do
                    !
                END IF
                !
                !$omp end parallel
                !
            ELSE ! R -> G
                !
                !$omp parallel default(none) &
                !$omp&         private(i, j, comm2s)
                !$omp&         shared(howmany, f, isgn, nnr_, dfft)
                !$omp&         shared(nsticks_z, n3, nx3)
                !$omp&         shared(nsticks_y, n2, nx2)
                !$omp&         shared(nsticks_x, n1, nx1)
                !$omp&         shared(nsticks_zx, nsticks_yx, aux)
                !
                !$omp single
                !
                !$omp taskgroup
                !
                DO i = 0, howmany - 1
                    !
                    !$omp task depend(out:aux(i*nnr_+1:(i+1)*nnr_))
                    !
                    CALL env_cft_1z(f(i * nnr_ + 1:), nsticks_x, n1, nx1, isgn, &
                                    aux(i * nnr_ + 1:))
                    !
                    !$omp end task
                    !
                    !$omp task depend(in:aux(i*nnr_+1:(i+1)*nnr_)) &
                    !$omp&     depend(out:f(i*nnr_+1:(i+1)*nnr_))
                    !
                    comm2s = dfft%comm2s(i + 1)
                    !
                    CALL env_fft_scatter_xy(dfft, f(i * nnr_ + 1:), aux(i * nnr_ + 1:), &
                                            nnr_, isgn, comm2s)
                    !
                    !$omp end task
                    !
                    !$omp task depend(in:f(i*nnr_+1:(i+1)*nnr_))
                    !
                    CALL env_cft_1z(f(i * nnr_ + 1:), nsticks_y, n2, nx2, isgn, &
                                    aux(i * nnr_ + 1:), in_place=.TRUE.)
                    !
                    !$omp end task
                    !
                END DO
                !
                !$omp end taskgroup
                !
                CALL env_fft_scatter_many_yz(dfft, aux, f, isgn, howmany)
                !
                !$omp end single
                !
                !$omp do
                DO i = 0, howmany - 1
                    !
                    CALL env_cft_1z(aux(nx3 * nsticks_zx * i + 1:), nsticks_z, n3, &
                                    nx3, isgn, f(i * nnr_ + 1:))
                    !
                END DO
                !$omp end do
                !
                !$omp do collapse(2)
                DO i = 0, howmany - 1
                    !
                    DO j = nsticks_z * nx3 + 1, nnr_
                        f(j + i * nnr_) = (0.0_DP, 0.0_DP)
                    END DO
                    !
                END DO
                !$omp end do
                !
                !$omp end parallel
                !
            END IF
            !
            !----------------------------------------------------------------------------
        END SUBROUTINE execute_using_tasks
        !--------------------------------------------------------------------------------
#endif
        !>
        !! This implementation requires thread-safe env_cft_1z
        !!
        !--------------------------------------------------------------------------------
        SUBROUTINE execute_using_threads()
            !----------------------------------------------------------------------------
            !
            IMPLICIT NONE
            !
            !----------------------------------------------------------------------------
            !
#if defined(_OPENMP)&&defined(__FFT_SCALAR_THREAD_SAFE)
            CALL env_errore(sub_name//'::execute_using_threads', &
                            'Needs thread-safe env_fft_scalar backend selected &
                            &at compile time.', 1)
#endif
            !
            IF (isgn > 0) THEN ! G -> R
                !
                !$omp parallel default(none) &
                !$omp&         private(i, j)
                !$omp&         shared(howmany, f, nnr_, isgn, dfft, aux)
                !$omp&         shared(nsticks_z, n3, nx3)
                !$omp&         shared(nsticks_y, n2, nx2)
                !$omp&         shared(nsticks_x, n1, nx1)
                !$omp&         shared(nsticks_zx, nsticks_yx)
                !
                !$omp do
                DO i = 0, howmany - 1
                    !
                    DO j = 1, nsticks_z * nx3
                        aux(j + i * nnr_) = f(j + i * nnr_)
                    END DO
                    !
                END DO
                !$omp end do
                !
                !$omp do
                DO i = 0, howmany - 1
                    !
                    CALL env_cft_1z(aux(i * nnr_ + 1:), nsticks_z, n3, nx3, isgn, &
                                    f(nx3 * nsticks_zx * i + 1:))
                    !
                END DO
                !$omp end do
                !
                !$omp single
                CALL env_fft_scatter_many_yz(dfft, f, dfft%aux, isgn, howmany)
                !$omp end single
                !
                !$omp do
                DO i = 0, howmany - 1
                    !
                    CALL env_cft_1z(aux(i * nnr_ + 1:), nsticks_y, n2, nx2, isgn, &
                                    f(nx2 * nsticks_yx * i + 1:))
                    !
                END DO
                !$omp end do
                !
                !$omp single
                CALL env_fft_scatter_many_xy(dfft, f, aux, isgn, howmany)
                !$omp end single
                !
                !$omp do
                DO i = 0, howmany - 1
                    !
                    CALL env_cft_1z(aux(i * nnr_ + 1:), nsticks_x, n1, nx1, isgn, &
                                    f(i * nnr_ + 1:))
                    !
                END DO
                !$omp end do
                !
                !$omp do
                DO i = 0, howmany - 1
                    !
                    IF (nsticks_x * nx1 < nnr_) THEN
                        !
                        DO j = nsticks_x * nx1 + 1, nnr_
                            f(j + i * nnr_) = (0.0_DP, 0.0_DP)
                        END DO
                        !
                    END IF
                    !
                END DO
                !$omp end do
                !
                !$omp end parallel
                !
            ELSE ! R -> G
                !
                !$omp parallel default(none) &
                !$omp&         private(i)
                !$omp&         shared(howmany, f, isgn, nnr_, dfft, aux)
                !$omp&         shared(nsticks_z, n3, nx3)
                !$omp&         shared(nsticks_y, n2, nx2)
                !$omp&         shared(nsticks_x, n1, nx1)
                !$omp&         shared(nsticks_zx, nsticks_yx)
                !
                !$omp do
                DO i = 0, howmany - 1
                    !
                    CALL env_cft_1z(f(i * nnr_ + 1:), nsticks_x, n1, nx1, isgn, &
                                    aux(i * nnr_ + 1:))
                    !
                END DO
                !$omp end do
                !
                !$omp single
                CALL env_fft_scatter_many_xy(dfft, f, aux, isgn, howmany)
                !$omp end single
                !
                !$omp do
                DO i = 0, howmany - 1
                    !
                    CALL env_cft_1z(f(nx2 * nsticks_yx * i + 1:), nsticks_y, n2, nx2, &
                                    isgn, aux(i * nnr_ + 1:))
                    !
                END DO
                !$omp end do
                !
                !$omp single
                CALL env_fft_scatter_many_yz(dfft, f, dfft%aux, isgn, howmany)
                !$omp end single
                !
                !$omp do
                DO i = 0, howmany - 1
                    !
                    CALL env_cft_1z(f(nx3 * nsticks_zx * i + 1:), nsticks_z, n3, nx3, &
                                    isgn, aux(i * nnr_ + 1:))
                    !
                END DO
                !$omp end do
                !
                !$omp do
                DO i = 0, howmany - 1
                    !
                    DO j = 0, nsticks_z - 1
                        !
                        f(i * nnr_ + j * nx3 + 1:i * nnr_ + j * nx3 + n3) = &
                            aux(i * nnr_ + j * nx3 + 1:i * nnr_ + j * nx3 + n3)
                        !
                    END DO
                    !
                END DO
                !$omp end do
                !
                !$omp end parallel
                !
            END IF
            !
            !----------------------------------------------------------------------------
        END SUBROUTINE execute_using_threads
        !--------------------------------------------------------------------------------
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_many_cft3s
    !------------------------------------------------------------------------------------
    !>
    !! General purpose driver, including Task groups parallelization
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_tg_cft3s_2d(f, dfft, isgn)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: isgn ! fft direction
        !
        TYPE(env_fft_type_descriptor), INTENT(IN) :: dfft
        ! descriptor of fft data layout
        !
        COMPLEX(DP), INTENT(INOUT) :: f(:) ! array containing data to be transformed
        !
        INTEGER :: me_p
        INTEGER :: n1, n2, n3, nx1, nx2, nx3
        COMPLEX(DP), ALLOCATABLE :: aux(:)
        INTEGER :: planes(dfft%nr1x)
        !
        CHARACTER(LEN=80) :: sub_name = 'env_tg_cft3s_2d'
        !
        !--------------------------------------------------------------------------------
        !
        IF (dfft%has_task_groups) &
            CALL env_errore(sub_name, 'Task groups on large mesh not implemented', 1)
        !
        n1 = dfft%nr1
        n2 = dfft%nr2
        n3 = dfft%nr3
        nx1 = dfft%nr1x
        nx2 = dfft%nr2x
        nx3 = dfft%nr3x
        !
        ALLOCATE (aux(dfft%nnr))
        !
        me_p = dfft%mype + 1
        !
        !--------------------------------------------------------------------------------
        ! Transpose data for the 2-D FFT on the x-y plane
        !
        ! NOGRP*dfft%nnr: The length of aux and f
        ! nr3x          : The length of each Z-stick
        ! aux           : input - output
        ! f             : working space
        ! isgn          : type of scatter
        ! dfft%nsp(me)  : holds the number of Z-sticks proc. me has.
        ! dfft%nr3p     : number of planes per processor
        !
        IF (isgn > 0) THEN
            !
            CALL env_cft_1z(f, dfft%nsp(me_p), n3, nx3, isgn, aux)
            !
            planes = dfft%iplp
            !
            ! forward scatter from stick to planes
            CALL env_fft_scatter_2d(dfft, aux, nx3, dfft%nnr, f, dfft%nsp, dfft%nr3p, &
                                    isgn)
            !
            CALL env_cft_2xy(f, dfft%my_nr3p, n1, n2, nx1, nx2, isgn, planes)
            !
        ELSE
            planes = dfft%iplp
            !
            CALL env_cft_2xy(f, dfft%my_nr3p, n1, n2, nx1, nx2, isgn, planes)
            !
            ! backward scatter from stick to planes
            CALL env_fft_scatter_2d(dfft, aux, nx3, dfft%nnr, f, dfft%nsp, dfft%nr3p, &
                                    isgn)
            !
            CALL env_cft_1z(aux, dfft%nsp(me_p), n3, nx3, isgn, f)
            !
        END IF
        !
        DEALLOCATE (aux)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_tg_cft3s_2d
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE env_fft_parallel
!----------------------------------------------------------------------------------------
