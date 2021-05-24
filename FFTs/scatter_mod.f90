!----------------------------------------------------------------------------------------
!
! Copyright (C) 2018 ENVIRON (www.quantum-environment.org)
! Copyright (C) 2011 Quantum ESPRESSO group
!
!----------------------------------------------------------------------------------------
!
! This file is part of Environ version 2.0
!
! Environ 2.0 is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 2 of the License, or
! (at your option) any later version.
!
! Environ 2.0 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more detail, either the file
! `License' in the root directory of the present distribution, or
! online at <http://www.gnu.org/licenses/>.
!
!----------------------------------------------------------------------------------------
!
! Authors: Carlo Cavazzoni, modified by Paolo Giannozzi
!          Rewritten by Stefano de Gironcoli
!
!----------------------------------------------------------------------------------------
!>
!! FFT base Module
!!
!----------------------------------------------------------------------------------------
MODULE env_scatter_mod
    !------------------------------------------------------------------------------------
    !
    USE env_fft_param
    !
    USE env_fft_types, ONLY: env_fft_type_descriptor
    !
    !------------------------------------------------------------------------------------
    INTERFACE env_gather_grid
        MODULE PROCEDURE env_gather_real_grid, env_gather_complex_grid
    END INTERFACE
    !
    INTERFACE env_scatter_grid
        MODULE PROCEDURE env_scatter_real_grid, env_scatter_complex_grid
    END INTERFACE
    !
    !------------------------------------------------------------------------------------
    !
    PRIVATE
    !
    PUBLIC :: env_gather_grid, env_scatter_grid
    !
    !------------------------------------------------------------------------------------
CONTAINS
    !------------------------------------------------------------------------------------
    !>
    !! Gathers a distributed real-space FFT grid to dfft%root, that is,
    !! the first processor of input descriptor dfft - version for real arrays
    !!
    !! REAL*8  f_in  = distributed variable (dfft%nnr)
    !! REAL*8  f_out = gathered variable (dfft%nr1x*dfft%nr2x*dfft%nr3x)
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_gather_real_grid(dfft, f_in, f_out)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP), INTENT(IN) :: f_in(:)
        TYPE(env_fft_type_descriptor), INTENT(IN) :: dfft
        !
        REAL(DP), INTENT(INOUT) :: f_out(:)
        !
#if defined(__MPI)
        INTEGER :: proc, info, offset_in, offset_aux, ir3
        INTEGER :: displs(0:dfft%nproc - 1), recvcount(0:dfft%nproc - 1)
        REAL(DP), ALLOCATABLE :: f_aux(:)
        !
        CHARACTER(LEN=80) :: sub_name = 'env_gather_real_grid'
        !
        !--------------------------------------------------------------------------------
        !
        IF (SIZE(f_in) < dfft%nnr) &
            CALL env_fft_error(sub_name, ' f_in too small ', dfft%nnr - SIZE(f_in))
        !
        CALL env_start_clock(sub_name)
        !
        ALLOCATE (f_aux(dfft%nr1x * dfft%nr2x * dfft%my_nr3p))
        !
        !--------------------------------------------------------------------------------
        ! 1) gather within the comm2 communicator
        !
        displs = 0
        !
        DO proc = 0, (dfft%nproc2 - 1)
            recvcount(proc) = dfft%nr1x * dfft%nr2p(proc + 1)
            !
            IF (proc > 0) displs(proc) = displs(proc - 1) + recvcount(proc - 1)
            !
        END DO
        !
        offset_in = 1
        offset_aux = 1
        !
        DO ir3 = 1, dfft%my_nr3p
            info = 0
            !
            CALL MPI_GATHERV(f_in(offset_in), recvcount(dfft%mype2), &
                             MPI_DOUBLE_PRECISION, f_aux(offset_aux), recvcount, &
                             displs, MPI_DOUBLE_PRECISION, dfft%root, dfft%comm2, info)
            !
            CALL env_fft_error(sub_name, 'info<>0', info)
            !
            offset_in = offset_in + dfft%nr1x * dfft%my_nr2p
            offset_aux = offset_aux + dfft%nr1x * dfft%nr2
        END DO
        !
        !--------------------------------------------------------------------------------
        ! 2) gather within the comm3 communicator
        !
        displs = 0
        !
        DO proc = 0, (dfft%nproc3 - 1)
            recvcount(proc) = dfft%nr1x * dfft%nr2x * dfft%nr3p(proc + 1)
            !
            IF (proc > 0) displs(proc) = displs(proc - 1) + recvcount(proc - 1)
            !
        END DO
        !
        info = 0
        !
        CALL MPI_GATHERV(f_aux, recvcount(dfft%mype3), MPI_DOUBLE_PRECISION, &
                         f_out, recvcount, displs, MPI_DOUBLE_PRECISION, dfft%root, &
                         dfft%comm3, info)
        !
        CALL env_fft_error(sub_name, 'info<>0', info)
        !
        !--------------------------------------------------------------------------------
        ! The following check should be performed only on processor dfft%root
        ! otherwise f_out must be allocated on all processors even if not used
        !
        info = SIZE(f_out) - displs(dfft%nproc3 - 1) - recvcount(dfft%nproc3 - 1)
        !
        IF (info < 0) CALL env_fft_error(sub_name, ' f_out too small ', -info)
        !
        DEALLOCATE (f_aux)
        !
        CALL env_stop_clock(sub_name)
#else
        !
        CALL env_fft_error(sub_name, 'do not use in serial execution', 1)
#endif
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_gather_real_grid
    !------------------------------------------------------------------------------------
    !>
    !! Gathers a distributed real-space FFT grid to dfft%root, that is,
    !! the first processor of input descriptor dfft - complex arrays
    !!
    !! COMPLEX*16  f_in  = distributed variable (dfft%nnr)
    !! COMPLEX*16  f_out = gathered variable (dfft%nr1x*dfft%nr2x*dfft%nr3x)
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_gather_complex_grid(dfft, f_in, f_out)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        COMPLEX(DP), INTENT(IN) :: f_in(:)
        TYPE(env_fft_type_descriptor), INTENT(IN) :: dfft
        !
        COMPLEX(DP), INTENT(INOUT) :: f_out(:)
        !
        COMPLEX(DP), ALLOCATABLE :: f_aux(:)
        !
#if defined(__MPI)
        INTEGER :: proc, info, offset_in, offset_aux, ir3
        INTEGER :: displs(0:dfft%nproc - 1), recvcount(0:dfft%nproc - 1)
        !
        CHARACTER(LEN=80) :: sub_name = 'env_gather_complex_grid'
        !
        !--------------------------------------------------------------------------------
        !
        CALL env_start_clock(sub_name)
        !
        IF (2 * SIZE(f_in) < dfft%nnr) &
            CALL env_fft_error(sub_name, ' f_in too small ', dfft%nnr - SIZE(f_in))
        !
        ALLOCATE (f_aux(dfft%nr1x * dfft%nr2x * dfft%my_nr3p))
        !
        !--------------------------------------------------------------------------------
        ! 1) gather within the comm2 communicator
        !
        displs = 0
        !
        DO proc = 0, (dfft%nproc2 - 1)
            recvcount(proc) = 2 * dfft%nr1x * dfft%nr2p(proc + 1)
            !
            IF (proc > 0) displs(proc) = displs(proc - 1) + recvcount(proc - 1)
            !
        END DO
        !
        offset_in = 1
        offset_aux = 1
        !
        DO ir3 = 1, dfft%my_nr3p
            info = 0
            !
            CALL MPI_GATHERV(f_in(offset_in), recvcount(dfft%mype2), &
                             MPI_DOUBLE_PRECISION, f_aux(offset_aux), recvcount, &
                             displs, MPI_DOUBLE_PRECISION, dfft%root, dfft%comm2, info)
            !
            CALL env_fft_error(sub_name, 'info<>0', info)
            !
            offset_in = offset_in + dfft%nr1x * dfft%my_nr2p
            offset_aux = offset_aux + dfft%nr1x * dfft%nr2
        END DO
        !
        !--------------------------------------------------------------------------------
        ! 2) gather within the comm3 communicator
        !
        displs = 0
        !
        DO proc = 0, (dfft%nproc3 - 1)
            recvcount(proc) = 2 * dfft%nr1x * dfft%nr2x * dfft%nr3p(proc + 1)
            !
            IF (proc > 0) displs(proc) = displs(proc - 1) + recvcount(proc - 1)
            !
        END DO
        !
        !--------------------------------------------------------------------------------
        ! The following check should be performed only on processor dfft%root
        ! otherwise f_out must be allocated on all processors even if not used
        !
        info = 2 * SIZE(f_out) - displs(dfft%nproc3 - 1) - recvcount(dfft%nproc3 - 1)
        FLUSH (stdout)
        !
        IF (info < 0) CALL env_fft_error(sub_name, ' f_out too small ', -info)
        !
        info = 0
        !
        CALL MPI_GATHERV(f_aux, recvcount(dfft%mype3), MPI_DOUBLE_PRECISION, &
                         f_out, recvcount, displs, MPI_DOUBLE_PRECISION, dfft%root, &
                         dfft%comm3, info)
        !
        CALL env_fft_error(sub_name, 'info<>0', info)
        !
        DEALLOCATE (f_aux)
        !
        CALL env_stop_clock(sub_name)
#else
        !
        CALL env_fft_error(sub_name, 'do not use in serial execution', 1)
#endif
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_gather_complex_grid
    !------------------------------------------------------------------------------------
    !>
    !! Scatters a real-space FFT grid from dfft%root, first processor of
    !! input descriptor dfft, to all others - opposite of "env_gather_grid"
    !!
    !! REAL*8  f_in  = gathered variable (dfft%nr1x*dfft%nr2x*dfft%nr3x)
    !! REAL*8  f_out = distributed variable (dfft%nnr)
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_scatter_real_grid(dfft, f_in, f_out)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP), INTENT(IN) :: f_in(:)
        TYPE(env_fft_type_descriptor), INTENT(IN) :: dfft
        !
        REAL(DP), INTENT(INOUT) :: f_out(:)
        !
        REAL(DP), ALLOCATABLE :: f_aux(:)
        !
#if defined(__MPI)
        INTEGER :: proc, info, offset_in, offset_aux, ir3
        INTEGER :: displs(0:dfft%nproc - 1), sendcount(0:dfft%nproc - 1)
        !
        CHARACTER(LEN=80) :: sub_name = 'env_scatter_real_grid'
        !
        !--------------------------------------------------------------------------------
        !
        CALL env_start_clock(sub_name)
        !
        ALLOCATE (f_aux(dfft%nr1x * dfft%nr2x * dfft%my_nr3p))
        !
        !--------------------------------------------------------------------------------
        ! 1) scatter within the comm3 communicator
        !
        displs = 0
        !
        DO proc = 0, (dfft%nproc3 - 1)
            sendcount(proc) = dfft%nr1x * dfft%nr2x * dfft%nr3p(proc + 1)
            !
            IF (proc > 0) displs(proc) = displs(proc - 1) + sendcount(proc - 1)
            !
        END DO
        !
        info = SIZE(f_in) - displs(dfft%nproc3 - 1) - sendcount(dfft%nproc3 - 1)
        !
        IF (info < 0) CALL env_fft_error(sub_name, ' f_in too small ', -info)
        !
        info = 0
        !
        CALL MPI_SCATTERV(f_in, sendcount, displs, MPI_DOUBLE_PRECISION, &
                          f_aux, sendcount(dfft%mype3), MPI_DOUBLE_PRECISION, &
                          dfft%root, dfft%comm3, info)
        !
        CALL env_fft_error('env_scatter_real_grid', 'info<>0', info)
        !
        !--------------------------------------------------------------------------------
        ! 2) scatter within the comm2 communicator
        !
        IF (SIZE(f_out) < dfft%nnr) &
            CALL env_fft_error(sub_name, ' f_out too small ', dfft%nnr - SIZE(f_out))
        !
        displs = 0
        f_out = 0.0D0
        !
        DO proc = 0, (dfft%nproc2 - 1)
            sendcount(proc) = dfft%nr1x * dfft%nr2p(proc + 1)
            !
            IF (proc > 0) displs(proc) = displs(proc - 1) + sendcount(proc - 1)
            !
        END DO
        !
        offset_in = 1
        offset_aux = 1
        !
        DO ir3 = 1, dfft%my_nr3p
            info = 0
            !
            CALL MPI_SCATTERV(f_aux(offset_aux), sendcount, displs, &
                              MPI_DOUBLE_PRECISION, f_out(offset_in), &
                              sendcount(dfft%mype2), MPI_DOUBLE_PRECISION, &
                              dfft%root, dfft%comm2, info)
            !
            CALL env_fft_error(sub_name, 'info<>0', info)
            !
            offset_in = offset_in + dfft%nr1x * dfft%my_nr2p
            offset_aux = offset_aux + dfft%nr1x * dfft%nr2
        END DO
        !
        DEALLOCATE (f_aux)
        !
        CALL env_stop_clock(sub_name)
#else
        !
        CALL env_fft_error(sub_name, 'do not use in serial execution', 1)
#endif
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_scatter_real_grid
    !------------------------------------------------------------------------------------
    !>
    !! Scatters a real-space FFT grid from dfft%root, first processor of
    !! input descriptor dfft, to all others - opposite of "env_gather_grid"
    !!
    !! COMPLEX*16  f_in  = gathered variable (dfft%nr1x*dfft%nr2x*dfft%nr3x)
    !! COMPLEX*16  f_out = distributed variable (dfft%nnr)
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_scatter_complex_grid(dfft, f_in, f_out)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        COMPLEX(DP), INTENT(IN) :: f_in(:)
        TYPE(env_fft_type_descriptor), INTENT(IN) :: dfft
        !
        COMPLEX(DP), INTENT(INOUT) :: f_out(:)
        !
        COMPLEX(DP), ALLOCATABLE :: f_aux(:)
        !
#if defined(__MPI)
        INTEGER :: proc, info, offset_in, offset_aux, ir3
        INTEGER :: displs(0:dfft%nproc - 1), sendcount(0:dfft%nproc - 1)
        !
        CHARACTER(LEN=80) :: sub_name = 'env_scatter_complex_grid'
        !
        !--------------------------------------------------------------------------------
        !
        CALL env_start_clock(sub_name)
        !
        ALLOCATE (f_aux(dfft%nr1x * dfft%nr2x * dfft%my_nr3p))
        !
        !--------------------------------------------------------------------------------
        ! 1) scatter within the comm3 communicator
        !
        displs = 0
        !
        DO proc = 0, (dfft%nproc3 - 1)
            sendcount(proc) = 2 * dfft%nr1x * dfft%nr2x * dfft%nr3p(proc + 1)
            !
            IF (proc > 0) displs(proc) = displs(proc - 1) + sendcount(proc - 1)
            !
        END DO
        !
        info = 2 * SIZE(f_in) - displs(dfft%nproc3 - 1) - sendcount(dfft%nproc3 - 1)
        !
        IF (info < 0) &
            CALL env_fft_error(sub_name, ' f_in too small ', -info)
        !
        info = 0
        !
        CALL MPI_SCATTERV(f_in, sendcount, displs, MPI_DOUBLE_PRECISION, &
                          f_aux, sendcount(dfft%mype3), MPI_DOUBLE_PRECISION, &
                          dfft%root, dfft%comm3, info)
        !
        CALL env_fft_error(sub_name, 'info<>0', info)
        !
        !--------------------------------------------------------------------------------
        ! 2) scatter within the comm2 communicator
        !
        IF (SIZE(f_out) < dfft%nnr) &
            CALL env_fft_error(sub_name, ' f_out too small ', dfft%nnr - SIZE(f_out))
        !
        displs = 0
        f_out = 0.0D0
        !
        DO proc = 0, (dfft%nproc2 - 1)
            sendcount(proc) = 2 * dfft%nr1x * dfft%nr2p(proc + 1)
            !
            IF (proc > 0) displs(proc) = displs(proc - 1) + sendcount(proc - 1)
            !
        END DO
        !
        offset_in = 1
        offset_aux = 1
        !
        DO ir3 = 1, dfft%my_nr3p
            info = 0
            !
            CALL MPI_SCATTERV(f_aux(offset_aux), sendcount, displs, &
                              MPI_DOUBLE_PRECISION, f_out(offset_in), &
                              sendcount(dfft%mype2), MPI_DOUBLE_PRECISION, &
                              dfft%root, dfft%comm2, info)
            !
            CALL env_fft_error(sub_name, 'info<>0', info)
            !
            offset_in = offset_in + dfft%nr1x * dfft%my_nr2p
            offset_aux = offset_aux + dfft%nr1x * dfft%nr2x
        END DO
        !
        !--------------------------------------------------------------------------------
        ! The following check should be performed only on processor dfft%root
        ! otherwise f_in must be allocated on all processors even if not used
        !
        DEALLOCATE (f_aux)
        !
        CALL env_stop_clock(sub_name)
#else
        !
        CALL env_fft_error(sub_name, 'do not use in serial execution', 1)
#endif
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_scatter_complex_grid
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE env_scatter_mod
!----------------------------------------------------------------------------------------
