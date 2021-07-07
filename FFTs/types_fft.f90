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
! Authors: Modified by Edan Bainglass
!
!----------------------------------------------------------------------------------------
!
#if defined(__CUDA)
#define DEV_ATTRIBUTES , DEVICE
#else
#define DEV_ATTRIBUTES
#endif
!
!----------------------------------------------------------------------------------------
!>
!!
!----------------------------------------------------------------------------------------
MODULE env_types_fft
    !------------------------------------------------------------------------------------
    !
    USE env_fft_param
    !
    USE env_fft_support, ONLY: env_good_fft_order, env_good_fft_dimension
    !
    USE env_base_stick, ONLY: env_sticks_map, env_sticks_map_allocate, &
                              env_get_sticks
    !
    USE omp_lib
    !
#if defined(__CUDA)
    USE cudafor
#endif
    !
    !------------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    !------------------------------------------------------------------------------------
    !>
    !! Data type for FFT descriptor
    !!
    !------------------------------------------------------------------------------------
    TYPE env_fft_type_descriptor
        !--------------------------------------------------------------------------------
        !
        INTEGER :: nr1 = 0 !
        INTEGER :: nr2 = 0 ! effective FFT dimensions of the 3D grid (global)
        INTEGER :: nr3 = 0 !
        INTEGER :: nr1x = 0 ! FFT grids leading dimensions
        INTEGER :: nr2x = 0 ! dimensions of the arrays for the 3D grid (global)
        INTEGER :: nr3x = 0 ! may differ from nr1 ,nr2 ,nr3 in order to boost performance
        !
        !--------------------------------------------------------------------------------
        ! Parallel layout: in reciprocal space data are organized in columns (sticks)
        !                  along the third direction and distributed across nproc
        !                  processors. In real space data are distributed in blocks
        !                  comprising sections of the Y and Z directions and complete
        !                  rows in the X direction in a matrix of  nproc2 x nproc3
        !                  processors. nproc = nproc2 x nproc3 and additional
        !                  communicators are introduced for data redistribution across
        !                  matrix columns and rows.
        !
        !--------------------------------------------------------------------------------
        ! Communicators and processor coordinates
        !
        LOGICAL :: lpara = .FALSE. ! .TRUE. if parallel FFT is active
        LOGICAL :: lgamma = .FALSE. ! .TRUE. if the grid has Gamma symmetry
        INTEGER :: root = 0 ! root processor
        !
        INTEGER :: comm = MPI_COMM_NULL ! communicator for main fft group
        INTEGER :: comm2 = MPI_COMM_NULL ! ...for fft group along the second direction
        INTEGER :: comm3 = MPI_COMM_NULL ! ...for fft group along the third direction
        INTEGER :: nproc = 1 ! number of processor in main fft group
        INTEGER :: nproc2 = 1 ! ...in fft group along the second direction
        INTEGER :: nproc3 = 1 ! ...in fft group along the third direction
        INTEGER :: mype = 0 ! my processor id (starting from 0) in fft main communicator
        INTEGER :: mype2 = 0 ! ...in fft communicator along the second direction (nproc2)
        INTEGER :: mype3 = 0 ! ...in fft communicator along the third direction (nproc3)
        !
        INTEGER, ALLOCATABLE :: iproc(:, :), iproc2(:), iproc3(:)
        ! subcommunicators proc mapping (starting from 1)
        !
        !--------------------------------------------------------------------------------
        ! FFT distributed data dimensions and indices
        !
        INTEGER :: my_nr3p = 0
        ! size of the "Z" section for this processor = nr3p( mype3 + 1 ) ~ nr3/nproc3
        !
        INTEGER :: my_nr2p = 0
        ! size of the "Y" section for this processor = nr2p( mype2 + 1 ) ~ nr2/nproc2
        !
        INTEGER :: my_i0r3p = 0 ! offset of the first "Z" element of this proc
        !                         in the nproc3 group = i0r3p( mype3 + 1 )
        !
        INTEGER :: my_i0r2p = 0 ! offset of the first "Y" element of this proc
        !                         in the nproc2 group = i0r2p( mype2 + 1 )
        !
        INTEGER, ALLOCATABLE :: nr3p(:) ! size of the "Z" section of each processor
        !                                 in the nproc3 group along Z
        !
        INTEGER, ALLOCATABLE :: nr3p_offset(:) ! offset of the "Z" section of each
        !                                        processor in the nproc3 group along Z
        !
        INTEGER, ALLOCATABLE :: nr2p(:) ! size of the "Y" section of each processor
        !                                 in the nproc2 group along Y
        !
        INTEGER, ALLOCATABLE :: nr2p_offset(:) ! offset of the "Y" section of each
        !                                        processor in the nproc2 group along Y
        !
        INTEGER, ALLOCATABLE :: nr1p(:) ! number of active "X" values ( potential )
        !                                 for a given proc in the nproc2 group
        !
        INTEGER, ALLOCATABLE :: i0r3p(:) ! offset of the first "Z" element of each
        !                                  proc in the nproc3 group (starting from 0)
        !
        INTEGER, ALLOCATABLE :: i0r2p(:) ! offset of the first "Y" element of each
        !                                  proc in the nproc2 group (starting from 0)
        !
        INTEGER, ALLOCATABLE :: ir1p(:) ! if >0 ir1p(m1) is the incremental index of
        !                                 the active ( potential ) X value of this proc
        !
        INTEGER, ALLOCATABLE :: indp(:, :) ! is the inverse of ir1p
        !
        INTEGER, POINTER DEV_ATTRIBUTES :: ir1p_d(:), indp_d(:, :), nr1p_d(:)
        ! duplicated version of the arrays declared above
        !
        INTEGER :: nst ! total number of sticks ( potential )
        !
        INTEGER, ALLOCATABLE :: nsp(:) ! number of sticks per processor ( potential )
        !                                using proc index starting from 1 that is on
        !                                proc mype -> nsp( mype + 1 )
        !
        INTEGER, ALLOCATABLE :: nsp_offset(:, :)
        ! offset of sticks per processor ( potential )
        !
        INTEGER, ALLOCATABLE :: ngl(:)
        ! per proc. no. of non zero charge density/potential components
        !
        INTEGER :: ngm ! my no. of non zero charge density/potential components
        ! ngm = dfftp%ngl(dfftp%mype + 1)
        ! with gamma sym, ngm = (dfftp%ngl(dfftp%mype + 1) + 1) / 2
        !
        INTEGER, ALLOCATABLE :: iplp(:) ! if > 0 is the iproc2 processor owning the
        !                                 active "X" value ( potential )
        !
        INTEGER :: nnp = 0 ! number of 0 and non 0 sticks in a plane ( ~nr1*nr2/nproc )
        INTEGER :: nnr = 0 ! local number of FFT grid elements  ( ~nr1*nr2*nr3/nproc )
        ! size of the arrays allocated for the FFT, local to each processor:
        ! - in parallel execution may differ from nr1x*nr2x*nr3x
        ! - not to be confused either with nr1*nr2*nr3
        !
        INTEGER :: nnr_tg = 0
        ! local number of grid elements for task group FFT ( ~nr1*nr2*nr3/proc3 )
        !
        INTEGER, ALLOCATABLE :: iss(:) ! index of the first rho stick on each proc
        !
        INTEGER, ALLOCATABLE :: isind(:)
        ! for each position in the plane indicate the stick index
        !
        INTEGER, ALLOCATABLE :: ismap(:)
        ! for each stick in the plane indicate the position
        !
        INTEGER, ALLOCATABLE :: nl(:) ! position of the G vec in the FFT grid
        !
        INTEGER, ALLOCATABLE :: nlm(:)
        ! with gamma sym. position of -G vec in the FFT grid
        !
        INTEGER, POINTER DEV_ATTRIBUTES :: ismap_d(:), nl_d(:), nlm_d(:)
        ! duplication of the variables defined above
        !
        !--------------------------------------------------------------------------------
        ! task group ALLTOALL communication layout
        !
        INTEGER, ALLOCATABLE :: tg_snd(:)
        ! number of elements to be sent in task group redistribution
        !
        INTEGER, ALLOCATABLE :: tg_rcv(:)
        ! number of elements to be received in task group redistribution
        !
        INTEGER, ALLOCATABLE :: tg_sdsp(:)
        ! send displacement for task group A2A communication
        !
        INTEGER, ALLOCATABLE :: tg_rdsp(:)
        ! receive displacement for task group A2A communicattion
        !
        LOGICAL :: has_task_groups = .FALSE.
        LOGICAL :: use_pencil_decomposition = .TRUE.
        !
        CHARACTER(LEN=12) :: rho_clock_label = ' '
        !
        INTEGER :: grid_id
        !
#if defined(__CUDA)
        INTEGER(kind=cuda_stream_kind), ALLOCATABLE :: stream_scatter_yz(:)
        INTEGER(kind=cuda_stream_kind), ALLOCATABLE :: stream_many(:)
        INTEGER :: nstream_many = 16
        !
        INTEGER(kind=cuda_stream_kind) :: a2a_comp
        INTEGER(kind=cuda_stream_kind), ALLOCATABLE :: bstreams(:)
        TYPE(cudaEvent), ALLOCATABLE :: bevents(:)
        !
        INTEGER :: batchsize = 16 ! how many ffts to batch together
        INTEGER :: subbatchsize = 4 ! size of subbatch for pipelining
        !
#if defined(__IPC)
        INTEGER :: IPC_PEER(16) ! this is used for IPC that is not imlpemented yet.
#endif
        INTEGER, ALLOCATABLE :: srh(:, :) ! Isend/recv handles by subbatch
#endif
        COMPLEX(DP), ALLOCATABLE :: aux(:)
#if defined(__FFT_OPENMP_TASKS)
        INTEGER, ALLOCATABLE :: comm2s(:) ! multiple communicator for the fft group
        !                                   along the second direction
        !
        INTEGER, ALLOCATABLE :: comm3s(:) ! multiple communicator for the fft group
        !                                   along the third direction
        !
#endif
        !
        !--------------------------------------------------------------------------------
    END TYPE
    !------------------------------------------------------------------------------------
    !
    INTEGER :: incremental_grid_identifier = 0
    !
    !------------------------------------------------------------------------------------
    !
    PRIVATE
    !
    PUBLIC :: env_fft_type_descriptor, env_fft_type_allocate, env_fft_type_deallocate, &
              env_fft_type_init, env_fft_stick_index
    !
    !------------------------------------------------------------------------------------
CONTAINS
    !------------------------------------------------------------------------------------
    !>
    !! Routine allocating arrays of env_fft_type_descriptor, called by env_fft_type_init
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_fft_type_allocate(desc, at, bg, gcutm, comm, fft_fact, nyfft)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP), DIMENSION(3, 3), INTENT(IN) :: at, bg
        REAL(DP), INTENT(IN) :: gcutm
        INTEGER, INTENT(IN) :: comm ! mype starting from 0
        INTEGER, INTENT(IN), OPTIONAL :: fft_fact(3)
        INTEGER, INTENT(IN), OPTIONAL :: nyfft
        !
        TYPE(env_fft_type_descriptor) :: desc
        INTEGER :: nx, ny, ierr, nzfft, i, nsubbatches
        INTEGER :: mype, root, nproc, iproc, iproc2, iproc3 ! mype starting from 0
        INTEGER :: color, key
        !
        CHARACTER(LEN=80) :: sub_name = 'env_fft_type_allocate'
        !
        !--------------------------------------------------------------------------------
        !
        IF (ALLOCATED(desc%nsp)) &
            CALL env_errore(sub_name, 'FFT arrays already allocated', 1)
        !
        desc%comm = comm
        !
#if defined(__MPI)
        IF (desc%comm == MPI_COMM_NULL) &
            CALL env_errore(sub_name, 'FFT communicator is null', 1)
#endif
        !
        root = 0
        mype = 0
        nproc = 1
        !
#if defined(__MPI)
        CALL MPI_COMM_RANK(comm, mype, ierr)
        !
        CALL MPI_COMM_SIZE(comm, nproc, ierr)
#endif
        !
        desc%root = root
        desc%mype = mype
        desc%nproc = nproc
        !
        IF (PRESENT(nyfft)) THEN
            !
            ! check on yfft group dimension
            CALL env_errore(sub_name, 'MOD(nproc,nyfft) /= 0', MOD(nproc, nyfft))
            !
#if defined(__MPI)
#if defined(ZCOMPACT)
            nzfft = nproc / nyfft
            color = mype / nzfft
            key = MOD(mype, nzfft)
#else
            color = MOD(mype, nyfft)
            key = mype / nyfft
#endif
            !
            !----------------------------------------------------------------------------
            ! Processes with the same key are in the same group along Y
            !
            CALL MPI_COMM_SPLIT(comm, key, color, desc%comm2, ierr)
            !
            CALL MPI_COMM_RANK(desc%comm2, desc%mype2, ierr)
            !
            CALL MPI_COMM_SIZE(desc%comm2, desc%nproc2, ierr)
            !
            !----------------------------------------------------------------------------
            ! processes with the same color are in the same group along Z
            !
            CALL MPI_COMM_SPLIT(comm, color, key, desc%comm3, ierr)
            !
            CALL MPI_COMM_RANK(desc%comm3, desc%mype3, ierr)
            !
            CALL MPI_COMM_SIZE(desc%comm3, desc%nproc3, ierr)
            !
#if defined(__FFT_OPENMP_TASKS)
            ALLOCATE (desc%comm2s(omp_get_max_threads()))
            ALLOCATE (desc%comm3s(omp_get_max_threads()))
            !
            DO i = 1, OMP_GET_MAX_THREADS()
                !
                CALL MPI_COMM_DUP(desc%comm2, desc%comm2s(i), ierr)
                !
                CALL MPI_COMM_DUP(desc%comm3, desc%comm3s(i), ierr)
                !
            END DO
            !
#endif
#else
            desc%comm2 = desc%comm
            desc%mype2 = desc%mype
            desc%nproc2 = desc%nproc
            desc%comm3 = desc%comm
            desc%mype3 = desc%mype
            desc%nproc3 = desc%nproc
#endif
            !
        END IF
        !
        ALLOCATE (desc%iproc(desc%nproc2, desc%nproc3), desc%iproc2(desc%nproc), &
                  desc%iproc3(desc%nproc))
        !
        DO iproc = 1, desc%nproc
#if defined(ZCOMPACT)
            iproc3 = MOD(iproc - 1, desc%nproc3) + 1
            iproc2 = (iproc - 1) / desc%nproc3 + 1
#else
            iproc2 = MOD(iproc - 1, desc%nproc2) + 1
            iproc3 = (iproc - 1) / desc%nproc2 + 1
#endif
            desc%iproc2(iproc) = iproc2
            desc%iproc3(iproc) = iproc3
            desc%iproc(iproc2, iproc3) = iproc
        END DO
        !
        CALL env_realspace_grid_init(desc, at, bg, gcutm, fft_fact)
        !
        ALLOCATE (desc%nr2p(desc%nproc2), desc%i0r2p(desc%nproc2))
        desc%nr2p = 0
        desc%i0r2p = 0
        ALLOCATE (desc%nr2p_offset(desc%nproc2))
        desc%nr2p_offset = 0
        ALLOCATE (desc%nr3p(desc%nproc3), desc%i0r3p(desc%nproc3))
        desc%nr3p = 0
        desc%i0r3p = 0
        ALLOCATE (desc%nr3p_offset(desc%nproc3))
        desc%nr3p_offset = 0
        !
        nx = desc%nr1x
        ny = desc%nr2x
        !
        ALLOCATE (desc%nsp(desc%nproc))
        desc%nsp = 0
        ALLOCATE (desc%nsp_offset(desc%nproc2, desc%nproc3))
        desc%nsp_offset = 0
        ALLOCATE (desc%ngl(desc%nproc))
        desc%ngl = 0
        ALLOCATE (desc%iss(desc%nproc))
        desc%iss = 0
        ALLOCATE (desc%isind(nx * ny))
        desc%isind = 0
        ALLOCATE (desc%ismap(nx * ny))
        desc%ismap = 0
        ALLOCATE (desc%nr1p(desc%nproc2))
        desc%nr1p = 0
        ALLOCATE (desc%ir1p(desc%nr1x))
        desc%ir1p = 0
        ALLOCATE (desc%indp(desc%nr1x, desc%nproc2))
        desc%indp = 0
        ALLOCATE (desc%iplp(nx))
        desc%iplp = 0
        !
        ALLOCATE (desc%tg_snd(desc%nproc2))
        desc%tg_snd = 0
        ALLOCATE (desc%tg_rcv(desc%nproc2))
        desc%tg_rcv = 0
        ALLOCATE (desc%tg_sdsp(desc%nproc2))
        desc%tg_sdsp = 0
        ALLOCATE (desc%tg_rdsp(desc%nproc2))
        desc%tg_rdsp = 0
        !
#if defined(__CUDA)
        ALLOCATE (desc%indp_d(desc%nr1x, desc%nproc2))
        desc%indp_d = 0
        !
        ALLOCATE (desc%nr1p_d(desc%nproc2))
        desc%nr1p_d = 0
        !
        ALLOCATE (desc%ir1p_d(desc%nr1x))
        desc%ir1p_d = 0
        ALLOCATE (desc%ismap_d(nx * ny))
        desc%ismap_d = 0
        !
        ALLOCATE (desc%stream_scatter_yz(desc%nproc3))
        !
        DO iproc = 1, desc%nproc3
            ierr = cudaStreamCreate(desc%stream_scatter_yz(iproc))
        END DO
        !
        ALLOCATE (desc%stream_many(desc%nstream_many))
        !
        DO i = 1, desc%nstream_many
            ierr = cudaStreamCreate(desc%stream_many(i))
            !
            IF (ierr /= 0) CALL env_errore(sub_name, 'Error creating stream', i)
            !
        END DO
        !
        ierr = cudaStreamCreate(desc%a2a_comp)
        nsubbatches = CEILING(REAL(desc%batchsize) / desc%subbatchsize)
        !
        ALLOCATE (desc%bstreams(nsubbatches))
        ALLOCATE (desc%bevents(nsubbatches))
        !
        DO i = 1, nsubbatches
            ierr = cudaStreamCreate(desc%bstreams(i))
            ierr = cudaEventCreate(desc%bevents(i))
        END DO
        !
        ALLOCATE (desc%srh(2 * nproc, nsubbatches))
#endif
        !
        incremental_grid_identifier = incremental_grid_identifier + 1
        desc%grid_id = incremental_grid_identifier
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_fft_type_allocate
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_fft_type_deallocate(desc)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(env_fft_type_descriptor) :: desc
        !
        INTEGER :: i, ierr, nsubbatches
        !
        !--------------------------------------------------------------------------------
        !
        IF (ALLOCATED(desc%nr2p)) DEALLOCATE (desc%nr2p)
        !
        IF (ALLOCATED(desc%nr2p_offset)) DEALLOCATE (desc%nr2p_offset)
        !
        IF (ALLOCATED(desc%nr3p_offset)) DEALLOCATE (desc%nr3p_offset)
        !
        IF (ALLOCATED(desc%i0r2p)) DEALLOCATE (desc%i0r2p)
        !
        IF (ALLOCATED(desc%nr3p)) DEALLOCATE (desc%nr3p)
        !
        IF (ALLOCATED(desc%i0r3p)) DEALLOCATE (desc%i0r3p)
        !
        IF (ALLOCATED(desc%nsp)) DEALLOCATE (desc%nsp)
        !
        IF (ALLOCATED(desc%nsp_offset)) DEALLOCATE (desc%nsp_offset)
        !
        IF (ALLOCATED(desc%ngl)) DEALLOCATE (desc%ngl)
        !
        IF (ALLOCATED(desc%iss)) DEALLOCATE (desc%iss)
        !
        IF (ALLOCATED(desc%isind)) DEALLOCATE (desc%isind)
        !
        IF (ALLOCATED(desc%ismap)) DEALLOCATE (desc%ismap)
        !
        IF (ALLOCATED(desc%nr1p)) DEALLOCATE (desc%nr1p)
        !
        IF (ALLOCATED(desc%ir1p)) DEALLOCATE (desc%ir1p)
        !
        IF (ALLOCATED(desc%indp)) DEALLOCATE (desc%indp)
        !
        IF (ALLOCATED(desc%iplp)) DEALLOCATE (desc%iplp)
        !
        IF (ALLOCATED(desc%iproc)) DEALLOCATE (desc%iproc)
        !
        IF (ALLOCATED(desc%iproc2)) DEALLOCATE (desc%iproc2)
        !
        IF (ALLOCATED(desc%iproc3)) DEALLOCATE (desc%iproc3)
        !
        IF (ALLOCATED(desc%tg_snd)) DEALLOCATE (desc%tg_snd)
        !
        IF (ALLOCATED(desc%tg_rcv)) DEALLOCATE (desc%tg_rcv)
        !
        IF (ALLOCATED(desc%tg_sdsp)) DEALLOCATE (desc%tg_sdsp)
        !
        IF (ALLOCATED(desc%tg_rdsp)) DEALLOCATE (desc%tg_rdsp)
        !
        IF (ALLOCATED(desc%nl)) DEALLOCATE (desc%nl)
        !
        IF (ALLOCATED(desc%nlm)) DEALLOCATE (desc%nlm)
        !
#if defined(__CUDA)
        IF (ALLOCATED(desc%ismap_d)) DEALLOCATE (desc%ismap_d)
        !
        IF (ALLOCATED(desc%ir1p_d)) DEALLOCATE (desc%ir1p_d)
        !
        IF (ALLOCATED(desc%indp_d)) DEALLOCATE (desc%indp_d)
        !
        IF (ALLOCATED(desc%nr1p_d)) DEALLOCATE (desc%nr1p_d)
        !
        IF (ALLOCATED(desc%stream_scatter_yz)) THEN
            !
            DO i = 1, desc%nproc3
                ierr = cudaStreamDestroy(desc%stream_scatter_yz(i))
            END DO
            !
            DEALLOCATE (desc%stream_scatter_yz)
        END IF
        !
        IF (ALLOCATED(desc%stream_many)) THEN
            !
            DO i = 1, desc%nstream_many
                ierr = cudaStreamDestroy(desc%stream_many(i))
            END DO
            !
            DEALLOCATE (desc%stream_many)
        END IF
        !
        IF (ALLOCATED(desc%nl_d)) DEALLOCATE (desc%nl_d)
        !
        IF (ALLOCATED(desc%nlm_d)) DEALLOCATE (desc%nlm_d)
        !
        IF (ALLOCATED(desc%srh)) DEALLOCATE (desc%srh)
        !
        ierr = cudaStreamDestroy(desc%a2a_comp)
        !
        IF (ALLOCATED(desc%bstreams)) THEN
            nsubbatches = CEILING(REAL(desc%batchsize) / desc%subbatchsize)
            !
            DO i = 1, nsubbatches
                ierr = cudaStreamDestroy(desc%bstreams(i))
                ierr = cudaEventDestroy(desc%bevents(i))
            END DO
            !
            DEALLOCATE (desc%bstreams)
            DEALLOCATE (desc%bevents)
        END IF
#endif
        !
        desc%comm = MPI_COMM_NULL
        !
#if defined(__MPI)
        IF (desc%comm2 /= MPI_COMM_NULL) CALL MPI_COMM_FREE(desc%comm2, ierr)
        !
        IF (desc%comm3 /= MPI_COMM_NULL) CALL MPI_COMM_FREE(desc%comm3, ierr)
#if defined(__FFT_OPENMP_TASKS)
        !
        DO i = 1, SIZE(desc%comm2s)
            !
            IF (desc%comm2s(i) /= MPI_COMM_NULL) CALL MPI_COMM_FREE(desc%comm2s(i), ierr)
            !
            IF (desc%comm3s(i) /= MPI_COMM_NULL) CALL MPI_COMM_FREE(desc%comm3s(i), ierr)
            !
        END DO
        !
        DEALLOCATE (desc%comm2s)
        DEALLOCATE (desc%comm3s)
#endif
#else
        desc%comm2 = MPI_COMM_NULL
        desc%comm3 = MPI_COMM_NULL
#endif
        !
        desc%nr1 = 0
        desc%nr2 = 0
        desc%nr3 = 0
        desc%nr1x = 0
        desc%nr2x = 0
        desc%nr3x = 0
        !
        desc%grid_id = 0
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_fft_type_deallocate
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_fft_type_set(desc, nst, ub, lb, idx, in1, in2, ncp, ngp, st, nmany)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(env_fft_type_descriptor) :: desc
        !
        INTEGER, INTENT(IN) :: nst ! total number of sticks
        !
        INTEGER, DIMENSION(3), INTENT(IN) :: ub, lb
        ! upper/lower bounds of real-space indices
        !
        INTEGER, INTENT(IN) :: nmany ! number of FFT bands
        INTEGER, INTENT(IN) :: idx(:) ! sorting index of the sticks
        INTEGER, INTENT(IN) :: in1(:) ! x-index of a stick
        INTEGER, INTENT(IN) :: in2(:) ! y-index of a stick
        INTEGER, INTENT(IN) :: ncp(:) ! number of rho columns per processor
        INTEGER, INTENT(IN) :: ngp(:) ! number of rho G-vectors per processor
        !
        INTEGER, INTENT(IN) :: st(lb(1):ub(1), lb(2):ub(2))
        ! stick owner of a given rho stick
        !
        INTEGER :: nsp(desc%nproc)
        INTEGER :: np, nq, i, is, iss, i1, i2, m1, m2, ip
        INTEGER :: ncpx, nr1px, nr2px, nr3px
        INTEGER :: nr1, nr2, nr3 ! size of real space grid
        INTEGER :: nr1x, nr2x, nr3x ! padded size of real space grid
        INTEGER :: ierr
        !
        CHARACTER(LEN=80) :: sub_name = 'env_fft_type_set'
        !
        !--------------------------------------------------------------------------------
        !
#if defined(__MPI)
#if defined(__FFT_OPENMP_TASKS)
        IF (nmany > OMP_GET_MAX_THREADS()) THEN
            !
            DO i = 1, SIZE(desc%comm2s)
                !
                IF (desc%comm2s(i) /= MPI_COMM_NULL) &
                    CALL MPI_COMM_FREE(desc%comm2s(i), ierr)
                !
                IF (desc%comm3s(i) /= MPI_COMM_NULL) &
                    CALL MPI_COMM_FREE(desc%comm3s(i), ierr)
                !
            END DO
            !
            DEALLOCATE (desc%comm2s)
            DEALLOCATE (desc%comm3s)
            ALLOCATE (desc%comm2s(nmany))
            ALLOCATE (desc%comm3s(nmany))
            !
            DO i = 1, nmany
                !
                CALL MPI_COMM_DUP(desc%comm2, desc%comm2s(i), ierr)
                !
                CALL MPI_COMM_DUP(desc%comm3, desc%comm3s(i), ierr)
                !
            END DO
            !
        END IF
#endif
#endif
        IF (.NOT. ALLOCATED(desc%nsp)) &
            CALL env_errore(sub_name, 'FFT arrays not yet allocated', 1)
        !
        IF (desc%nr1 == 0 .OR. desc%nr2 == 0 .OR. desc%nr3 == 0) &
            CALL env_errore(sub_name, 'FFT dimensions not yet set', 1)
        !
        !--------------------------------------------------------------------------------
        ! Set fft actual and leading dimensions to be used internally
        !
        nr1 = desc%nr1
        nr2 = desc%nr2
        nr3 = desc%nr3
        nr1x = desc%nr1x
        nr2x = desc%nr2x
        nr3x = desc%nr3x
        !
        IF ((nr1 > nr1x) .OR. (nr2 > nr2x) .OR. (nr3 > nr3x)) &
            CALL env_errore(sub_name, 'Wrong FFT dimensions', 1)
        !
        IF ((SIZE(desc%ngl) < desc%nproc) .OR. (SIZE(desc%iss) < desc%nproc) .OR. &
            (SIZE(desc%nr2p) < desc%nproc2) .OR. (SIZE(desc%i0r2p) < desc%nproc2) .OR. &
            (SIZE(desc%nr3p) < desc%nproc3) .OR. (SIZE(desc%i0r3p) < desc%nproc3)) &
            CALL env_errore(sub_name, 'Wrong descriptor dimensions', 2)
        !
        IF ((SIZE(idx) < nst) .OR. (SIZE(in1) < nst) .OR. (SIZE(in2) < nst)) &
            CALL env_errore(sub_name, 'Wrong number of stick dimensions', 3)
        !
        IF ((SIZE(ncp) < desc%nproc) .OR. (SIZE(ngp) < desc%nproc)) &
            CALL env_errore(sub_name, 'Wrong stick dimensions', 4)
        !
        !--------------------------------------------------------------------------------
        ! Set the number of "Y" values for each processor in the nproc2 group
        !
        np = nr2 / desc%nproc2
        nq = nr2 - np * desc%nproc2
        !
        desc%nr2p(1:desc%nproc2) = np
        ! assign a base value to all processors of the nproc2 group
        !
        !--------------------------------------------------------------------------------
        ! Assign an extra unit to the first nq processors of the nproc2 group
        !
        DO i = 1, nq
            desc%nr2p(i) = np + 1
        END DO
        !
        !--------------------------------------------------------------------------------
        ! Set the offset
        !
        desc%nr2p_offset(1) = 0
        !
        DO i = 1, desc%nproc2 - 1
            desc%nr2p_offset(i + 1) = desc%nr2p_offset(i) + desc%nr2p(i)
        END DO
        !
        desc%my_nr2p = desc%nr2p(desc%mype2 + 1)
        ! my_nr2p is the number of planes per processor
        ! of this processor in the Y group
        !
        !--------------------------------------------------------------------------------
        ! Find out the index of the starting plane on each proc
        !
        desc%i0r2p = 0
        !
        DO i = 2, desc%nproc2
            desc%i0r2p(i) = desc%i0r2p(i - 1) + desc%nr2p(i - 1)
        END DO
        !
        desc%my_i0r2p = desc%i0r2p(desc%mype2 + 1)
        ! my_i0r2p is the index-offset of the starting plane
        ! of this processor in the Y group
        !
        !--------------------------------------------------------------------------------
        ! Set the number of "Z" values for each processor in the nproc3 group
        !
        np = nr3 / desc%nproc3
        nq = nr3 - np * desc%nproc3
        desc%nr3p(1:desc%nproc3) = np ! assign a base value to all processors
        !
        !--------------------------------------------------------------------------------
        ! Assign an extra unit to the first nq processors of the nproc3 group
        !
        DO i = 1, nq
            desc%nr3p(i) = np + 1
        END DO
        !
        !--------------------------------------------------------------------------------
        ! Set the offset
        !
        desc%nr3p_offset(1) = 0
        !
        DO i = 1, desc%nproc3 - 1
            desc%nr3p_offset(i + 1) = desc%nr3p_offset(i) + desc%nr3p(i)
        END DO
        !
        desc%my_nr3p = desc%nr3p(desc%mype3 + 1)
        ! my_nr3p is the number of planes per processor
        ! of this processor in the Z group
        !
        !--------------------------------------------------------------------------------
        ! Find out the index of the starting plane on each proc
        !
        desc%i0r3p = 0
        !
        DO i = 2, desc%nproc3
            desc%i0r3p(i) = desc%i0r3p(i - 1) + desc%nr3p(i - 1)
        END DO
        !
        desc%my_i0r3p = desc%i0r3p(desc%mype3 + 1)
        ! my_i0r3p is the index-offset of the starting plane
        ! of this processor in the Z group
        !
        !--------------------------------------------------------------------------------
        ! Dimension of the xy plane. see ncplane
        !
        desc%nnp = nr1x * nr2x
        !
        desc%ngl(1:desc%nproc) = ngp(1:desc%nproc)
        ! local number of g vectors (rho) per processor
        !
        IF (SIZE(desc%isind) < (nr1x * nr2x)) &
            CALL env_errore(sub_name, 'Wrong descriptor dimensions, isind', 5)
        !
        IF (SIZE(desc%iplp) < (nr1x)) &
            CALL env_errore(sub_name, 'Wrong descriptor dimensions, ipl', 5)
        !
        IF (desc%my_nr3p == 0 .AND. (.NOT. desc%use_pencil_decomposition)) &
            CALL env_errore(sub_name, &
                            'There are processes with no planes. &
                            &Use pencil decomposition (-pd .true.)', 6)
        !
        !--------------------------------------------------------------------------------
        !
        ! 1. Temporarily store in the array "desc%isind" the index of the processor
        !    that own the corresponding stick (index of proc starting from 1)
        ! 2. Set the array elements of "desc%iplp" to one for that index corresponding
        !    to YZ planes containing at least one stick this are used in the FFT
        !    transform along Y
        !
        desc%isind = 0
        ! will contain the +ve or -ve of the processor number,
        ! if any, that owns the stick
        !
        desc%iplp = 0
        ! if > 0 is the nproc2 processor owning this ( potential ) X active plane
        !
        !--------------------------------------------------------------------------------
        ! Set nst to the proper number of sticks
        ! (the total number of 1d fft along z to be done)
        !
        desc%nst = 0
        !
        DO iss = 1, SIZE(idx)
            is = idx(iss)
            !
            IF (is < 1) CYCLE
            !
            i1 = in1(is)
            i2 = in2(is)
            !
            IF (st(i1, i2) > 0) THEN
                desc%nst = desc%nst + 1
                m1 = i1 + 1
                !
                IF (m1 < 1) m1 = m1 + nr1
                !
                m2 = i2 + 1
                !
                IF (m2 < 1) m2 = m2 + nr2
                !
                desc%isind(m1 + (m2 - 1) * nr1x) = st(i1, i2)
                desc%iplp(m1) = desc%iproc2(st(i1, i2))
                !
                IF (desc%lgamma) THEN
                    !
                    IF (i1 /= 0 .OR. i2 /= 0) desc%nst = desc%nst + 1
                    !
                    m1 = -i1 + 1
                    !
                    IF (m1 < 1) m1 = m1 + nr1
                    !
                    m2 = -i2 + 1
                    !
                    IF (m2 < 1) m2 = m2 + nr2
                    !
                    desc%isind(m1 + (m2 - 1) * nr1x) = st(-i1, -i2)
                    desc%iplp(m1) = desc%iproc2(st(-i1, -i2))
                END IF
                !
            END IF
            !
        END DO
        !
        !--------------------------------------------------------------------------------
        ! Count how many active X values per each nproc2 processor and set the
        ! incremental index of this one
        !
        desc%nr1p = 0
        desc%ir1p = 0
        desc%indp = 0
        !
        DO i1 = 1, nr1
            !
            IF (desc%iplp(i1) > 0) THEN
                desc%nr1p(desc%iplp(i1)) = desc%nr1p(desc%iplp(i1)) + 1
                desc%indp(desc%nr1p(desc%iplp(i1)), desc%iplp(i1)) = i1
            END IF
            !
            IF (desc%iplp(i1) == desc%mype2 + 1) desc%ir1p(i1) = desc%nr1p(desc%iplp(i1))
            !
        END DO
        !
        !--------------------------------------------------------------------------------
        ! Compute for each proc the global index ( starting from 0 ) of the first
        ! local stick ( desc%iss )
        !
        DO i = 1, desc%nproc
            !
            IF (i == 1) THEN
                desc%iss(i) = 0
            ELSE
                desc%iss(i) = desc%iss(i - 1) + ncp(i - 1)
            END IF
            !
        END DO
        !
        !--------------------------------------------------------------------------------
        ! iss(1:nproc) is the index offset of the first column of a given processor
        !
        IF (SIZE(desc%ismap) < (nst)) &
            CALL env_errore(sub_name, 'Wrong descriptor dimensions', 6)
        !
        !--------------------------------------------------------------------------------
        !
        ! 1. Set the array desc%ismap which maps stick indexes to
        !    position in the plane  ( iss )
        ! 2. Re-set the array "desc%isind", that maps position
        !    in the plane with stick indexes (it is the inverse of desc%ismap )
        !
        desc%ismap = 0
        ! will be the global xy stick index in the global list
        ! of processor-ordered sticks
        !
        nsp = 0 ! will be the number of sticks of a given processor
        !
        DO iss = 1, SIZE(desc%isind)
            ip = desc%isind(iss)
            !
            IF (ip > 0) THEN ! only operates on wave sticks
                nsp(ip) = nsp(ip) + 1
                desc%ismap(nsp(ip) + desc%iss(ip)) = iss
                !
                IF (ip == (desc%mype + 1)) THEN
                    !
                    desc%isind(iss) = nsp(ip)
                    ! isind is OVERWRITTEN as the ordered index
                    ! in this processor stick list
                    !
                ELSE
                    desc%isind(iss) = 0 ! zero otherwise...
                END IF
                !
            END IF
            !
        END DO
        !
        !--------------------------------------------------------------------------------
        ! Check number of sticks against the input value
        !
        IF (ANY(nsp(1:desc%nproc) /= ncp(1:desc%nproc))) THEN
            !
            DO ip = 1, desc%nproc
                WRITE (stdout, *) ' * ', ip, ' * ', nsp(ip), ' /= ', ncp(ip)
            END DO
            !
            CALL env_errore(sub_name, 'Inconsistent number of sticks', 7)
            !
        END IF
        !
        desc%nsp(1:desc%nproc) = nsp(1:desc%nproc) ! number of rho sticks per processor
        !
        DO ip = 1, desc%nproc3
            desc%nsp_offset(1, ip) = 0
            !
            DO i = 1, desc%nproc2 - 1
                desc%nsp_offset(i + 1, ip) = desc%nsp_offset(i, ip) + &
                                             desc%nsp(desc%iproc(i, ip))
            END DO
            !
        END DO
        !
        IF (.NOT. desc%lpara) THEN ! #TODO might not be necessary... check with Oliviero
            desc%isind = 0
            desc%iplp = 1
            !
            !----------------------------------------------------------------------------
            ! Here we are setting parameter as if we were in a serial code, sticks
            ! are along X dimension and not along Z
            !
            desc%nsp(1) = 0
            !
            DO i1 = lb(1), ub(1)
                !
                DO i2 = lb(2), ub(2)
                    m1 = i1 + 1
                    !
                    IF (m1 < 1) m1 = m1 + nr1
                    !
                    m2 = i2 + 1
                    !
                    IF (m2 < 1) m2 = m2 + nr2
                    !
                    desc%nsp(1) = desc%nsp(1) + 1
                    desc%isind(m1 + (m2 - 1) * nr1x) = 1 ! st( i1, i2 )
                END DO
                !
            END DO
            !
            !----------------------------------------------------------------------------
            ! If we are in a parallel run, but would like to use serial FFT, all
            ! tasks must have the same parameters as if serial run.
            !
            desc%nnr = nr1x * nr2x * nr3x
            desc%nnp = nr1x * nr2x
            desc%my_nr2p = nr2
            desc%nr2p = nr2
            desc%i0r2p = 0
            desc%my_nr3p = nr3
            desc%nr3p = nr3
            desc%i0r3p = 0
            desc%nsp = desc%nsp(1)
            desc%ngl = SUM(ngp)
        END IF
        !
        !--------------------------------------------------------------------------------
        ! Finally set fft local workspace dimension
        !
        nr1px = MAXVAL(desc%nr1p(1:desc%nproc2))
        ! maximum number of X values per processor in the nproc2 group
        !
        nr2px = MAXVAL(desc%nr2p(1:desc%nproc2))
        ! maximum number of planes per processor in the nproc2 group
        !
        nr3px = MAXVAL(desc%nr3p(1:desc%nproc3))
        ! maximum number of planes per processor in the nproc3 group
        !
        ncpx = MAXVAL(ncp(1:desc%nproc))
        ! maximum number of columns per processor (use potential sticks to be safe)
        !
        IF (desc%nproc == 1) THEN
            desc%nnr = nr1x * nr2x * nr3x
            desc%nnr_tg = desc%nnr * desc%nproc2
        ELSE
            !
            desc%nnr = MAX(ncpx * nr3x, nr1x * nr2px * nr3px)
            ! this is required to contain the local data in R and G space
            !
            ! this is required to use ALLTOALL instead of ALLTOALLV
            desc%nnr = MAX(desc%nnr, ncpx * nr3px * desc%nproc3, &
                           nr1px * nr2px * nr3px * desc%nproc2)
            !
            desc%nnr = MAX(1, desc%nnr)
            ! ensure that desc%nrr > 0 ( for extreme parallelism )
            !
            desc%nnr_tg = desc%nnr * desc%nproc2
        END IF
        !
        IF (desc%nr3x * desc%nsp(desc%mype + 1) > desc%nnr) &
            CALL env_errore(sub_name, 'Inconsistent desc%nnr', 1)
        !
        desc%tg_snd(1) = desc%nr3x * desc%nsp(desc%mype + 1)
        desc%tg_rcv(1) = desc%nr3x * desc%nsp(desc%iproc(1, desc%mype3 + 1))
        desc%tg_sdsp(1) = 0
        desc%tg_rdsp(1) = 0
        !
        DO i = 2, desc%nproc2
            desc%tg_snd(i) = desc%nr3x * desc%nsp(desc%mype + 1)
            desc%tg_rcv(i) = desc%nr3x * desc%nsp(desc%iproc(i, desc%mype3 + 1))
            desc%tg_sdsp(i) = desc%tg_sdsp(i - 1) + desc%nnr
            desc%tg_rdsp(i) = desc%tg_rdsp(i - 1) + desc%tg_rcv(i - 1)
        END DO
        !
#if defined(__CUDA)
        desc%ismap_d = desc%ismap
        desc%ir1p_d = desc%ir1p
        desc%indp_d = desc%indp
        desc%nr1p_d = desc%nr1p
#endif
        !
        IF (nmany > 1) ALLOCATE (desc%aux(nmany * desc%nnr))
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_fft_type_set
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_fft_type_init(dfft, smap, lgamma, lpara, comm, at, bg, gcut_in, &
                                 fft_fact, nyfft, nmany, use_pd)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        LOGICAL, INTENT(IN) :: lpara, lgamma
        REAL(DP), INTENT(IN) :: gcut_in
        REAL(DP), DIMENSION(3, 3), INTENT(IN) :: bg, at
        INTEGER, INTENT(IN) :: comm, nyfft, nmany
        INTEGER, OPTIONAL, INTENT(IN) :: fft_fact(3)
        LOGICAL, OPTIONAL, INTENT(IN) :: use_pd ! whether to use pencil decomposition
        !
        TYPE(env_fft_type_descriptor), INTENT(INOUT) :: dfft
        TYPE(env_sticks_map), INTENT(INOUT) :: smap
        !
        CHARACTER(LEN=80) :: sub_name = 'env_fft_type_init'
        !
        !--------------------------------------------------------------------------------
        ! Potential
        !
        INTEGER, ALLOCATABLE :: st(:, :)
        ! stick map, st(i, j) = number of G-vector in the stick
        ! whose x and y miller index are i and j
        !
        INTEGER, ALLOCATABLE :: nstp(:)
        ! number of sticks, nstp(ip) = number of stick for processor ip
        !
        INTEGER, ALLOCATABLE :: sstp(:) ! number of G-vectors, sstp(ip) = sum of the
        ! sticks length for processor ip = number of G-vectors owned by the processor ip
        !
        INTEGER :: nst ! local number of sticks
        !
        REAL(DP) :: gcut
        INTEGER :: ngm
        !
        !--------------------------------------------------------------------------------
        !
        gcut = gcut_in
        !
        IF (.NOT. ALLOCATED(dfft%nsp)) THEN
            !
            CALL env_fft_type_allocate(dfft, at, bg, gcut, comm, fft_fact=fft_fact, &
                                       nyfft=nyfft)
            !
        ELSE
            !
            IF (dfft%comm /= comm) &
                CALL env_errore(sub_name, &
                                'FFT already allocated with a different communicator', 1)
            !
        END IF
        !
        IF (PRESENT(use_pd)) dfft%use_pencil_decomposition = use_pd
        !
        IF ((.NOT. dfft%use_pencil_decomposition) .AND. (nyfft > 1)) &
            CALL env_errore(sub_name, &
                            'Slab decomposition and task groups not implemented', 1)
        !
        dfft%lpara = lpara ! this descriptor can be either a descriptor for a
        ! parallel FFT or a serial FFT even in parallel build
        !
        CALL env_sticks_map_allocate(smap, lgamma, dfft%lpara, dfft%nproc2, &
                                     dfft%iproc, dfft%iproc2, dfft%nr1, dfft%nr2, &
                                     dfft%nr3, bg, dfft%comm)
        !
        dfft%lgamma = smap%lgamma ! .TRUE. if the grid has Gamma symmetry
        !
        ALLOCATE (st(smap%lb(1):smap%ub(1), smap%lb(2):smap%ub(2)))
        ALLOCATE (nstp(smap%nproc))
        ALLOCATE (sstp(smap%nproc))
        !
        CALL env_get_sticks(smap, gcut, nstp, sstp, st, nst, ngm)
        !
        CALL env_fft_type_set(dfft, nst, smap%ub, smap%lb, smap%idx, smap%ist(:, 1), &
                              smap%ist(:, 2), nstp, sstp, st, nmany)
        !
        dfft%ngm = dfft%ngl(dfft%mype + 1)
        !
        IF (dfft%lgamma) dfft%ngm = (dfft%ngm + 1) / 2
        !
        IF (dfft%ngm /= ngm) CALL env_errore(sub_name, 'Wrong ngm', 1)
        !
        DEALLOCATE (st)
        DEALLOCATE (nstp)
        DEALLOCATE (sstp)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_fft_type_init
    !------------------------------------------------------------------------------------
    !>
    !! Sets optimal values for dfft%nr[123] and dfft%nr[123]x
    !! If input dfft%nr[123] are non-zero, leaves them unchanged
    !! If fft_fact is present, force nr[123] to be multiple of fft_fac([123])
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_realspace_grid_init(dfft, at, bg, gcutm, fft_fact)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP), DIMENSION(3, 3), INTENT(IN) :: at, bg
        REAL(DP), INTENT(IN) :: gcutm
        INTEGER, INTENT(IN), OPTIONAL :: fft_fact(3)
        !
        TYPE(env_fft_type_descriptor), INTENT(INOUT) :: dfft
        !
        !--------------------------------------------------------------------------------
        ! Calculate the size of the real-space dense grid for FFT
        ! first, an estimate of nr1,nr2,nr3, based on the max values
        ! of n_i indices in:   G = i*b_1 + j*b_2 + k*b_3
        ! We use G*a_i = n_i => n_i .le. |Gmax||a_i|
        !
        IF (dfft%nr1 == 0 .OR. dfft%nr2 == 0 .OR. dfft%nr3 == 0) THEN
            !
            dfft%nr1 = INT(SQRT(gcutm) * &
                           SQRT(at(1, 1)**2 + at(2, 1)**2 + at(3, 1)**2)) + 1
            !
            dfft%nr2 = INT(SQRT(gcutm) * &
                           SQRT(at(1, 2)**2 + at(2, 2)**2 + at(3, 2)**2)) + 1
            !
            dfft%nr3 = INT(SQRT(gcutm) * &
                           SQRT(at(1, 3)**2 + at(2, 3)**2 + at(3, 3)**2)) + 1
            !
#if defined(__DEBUG)
            WRITE (6, *) &
                SQRT(gcutm) * SQRT(at(1, 1)**2 + at(2, 1)**2 + at(3, 1)**2), dfft%nr1
            !
            WRITE (6, *) &
                SQRT(gcutm) * SQRT(at(1, 2)**2 + at(2, 2)**2 + at(3, 2)**2), dfft%nr2
            !
            WRITE (6, *) &
                SQRT(gcutm) * SQRT(at(1, 3)**2 + at(2, 3)**2 + at(3, 3)**2), dfft%nr3
#endif
            !
            CALL env_grid_set(dfft, bg, gcutm, dfft%nr1, dfft%nr2, dfft%nr3)
            !
            IF (PRESENT(fft_fact)) THEN
                dfft%nr1 = env_good_fft_order(dfft%nr1, fft_fact(1))
                dfft%nr2 = env_good_fft_order(dfft%nr2, fft_fact(2))
                dfft%nr3 = env_good_fft_order(dfft%nr3, fft_fact(3))
            ELSE
                dfft%nr1 = env_good_fft_order(dfft%nr1)
                dfft%nr2 = env_good_fft_order(dfft%nr2)
                dfft%nr3 = env_good_fft_order(dfft%nr3)
            END IF
            !
#if defined(__DEBUG)
        ELSE
            WRITE (stdout, '( /, 3X,"Info: using nr1, nr2, nr3 values from input" )')
#endif
        END IF
        !
        dfft%nr1x = env_good_fft_dimension(dfft%nr1)
        dfft%nr2x = dfft%nr2
        dfft%nr3x = env_good_fft_dimension(dfft%nr3)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_realspace_grid_init
    !------------------------------------------------------------------------------------
    !>
    !! This routine returns in nr1, nr2, nr3 the minimal 3D real-space FFT
    !! grid required to fit the G-vector sphere with G^2 <= gcut
    !! On input, nr1,nr2,nr3 must be set to values that match or exceed
    !! the largest i,j,k (Miller) indices in G(i,j,k) = i*b1 + j*b2 + k*b3
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE env_grid_set(dfft, bg, gcut, nr1, nr2, nr3)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(env_fft_type_descriptor), INTENT(IN) :: dfft
        REAL(DP), INTENT(IN) :: bg(3, 3), gcut
        !
        INTEGER, INTENT(INOUT) :: nr1, nr2, nr3
        !
        INTEGER :: i, j, k, nb(3)
        REAL(DP) :: gsq, g(3)
        !
        !--------------------------------------------------------------------------------
        ! Calculate moduli of G vectors and the range of indices where
        ! |G|^2 < gcut (in parallel whenever possible)
        !
        nb = 0
        !
        DO k = -nr3, nr3
            !
            !----------------------------------------------------------------------------
            ! me_image = processor number, starting from 0
            !
            IF (MOD(k + nr3, dfft%nproc) == dfft%mype) THEN
                !
                DO j = -nr2, nr2
                    !
                    DO i = -nr1, nr1
                        !
                        g(1) = DBLE(i) * bg(1, 1) + DBLE(j) * bg(1, 2) + &
                               DBLE(k) * bg(1, 3)
                        !
                        g(2) = DBLE(i) * bg(2, 1) + DBLE(j) * bg(2, 2) + &
                               DBLE(k) * bg(2, 3)
                        !
                        g(3) = DBLE(i) * bg(3, 1) + DBLE(j) * bg(3, 2) + &
                               DBLE(k) * bg(3, 3)
                        !
                        gsq = g(1)**2 + g(2)**2 + g(3)**2 ! calculate modulus
                        !
                        !----------------------------------------------------------------
                        ! Calculate maximum index
                        !
                        IF (gsq < gcut) THEN
                            nb(1) = MAX(nb(1), ABS(i))
                            nb(2) = MAX(nb(2), ABS(j))
                            nb(3) = MAX(nb(3), ABS(k))
                        END IF
                        !
                    END DO
                    !
                END DO
                !
            END IF
            !
        END DO
        !
#if defined(__MPI)
        CALL MPI_ALLREDUCE(MPI_IN_PLACE, nb, 3, MPI_INTEGER, MPI_MAX, dfft%comm, i)
#endif
        !
        !--------------------------------------------------------------------------------
        ! The size of the 3d FFT matrix depends upon the maximum indices. With
        ! the following choice, the sphere in G-space "touches" its periodic image
        !
        nr1 = 2 * nb(1) + 1
        nr2 = 2 * nb(2) + 1
        nr3 = 2 * nb(3) + 1
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE env_grid_set
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    PURE FUNCTION env_fft_stick_index(desc, i, j)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(env_fft_type_descriptor), INTENT(IN) :: desc
        INTEGER, INTENT(IN) :: i, j
        !
        INTEGER :: env_fft_stick_index
        INTEGER :: mc, m1, m2
        !
        !--------------------------------------------------------------------------------
        !
        m1 = MOD(i, desc%nr1) + 1
        !
        IF (m1 < 1) m1 = m1 + desc%nr1
        !
        m2 = MOD(j, desc%nr2) + 1
        !
        IF (m2 < 1) m2 = m2 + desc%nr2
        !
        mc = m1 + (m2 - 1) * desc%nr1x
        env_fft_stick_index = desc%isind(mc)
        !
        !--------------------------------------------------------------------------------
    END FUNCTION env_fft_stick_index
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE env_types_fft
!----------------------------------------------------------------------------------------
