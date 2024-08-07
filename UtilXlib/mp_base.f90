!
! Copyright (C) 2002-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!  Wrapper for MPI implementations that have problems with large messages
!


!  In some MPI implementation the communication subsystem
!  crashes when message exceeds a given size, so we need
!  to break down MPI communications in smaller pieces
!
#define __MSGSIZ_MAX 100000
#define __BCAST_MSGSIZ_MAX 100000

!  Some implementation of MPI (OpenMPI) if it is not well tuned for the given
!  network hardware (InfiniBand) tend to lose performance or get stuck inside
!  collective routines if processors are not well synchronized
!  A barrier fixes the problem
!
#define __USE_BARRIER


!=----------------------------------------------------------------------------=!
!
! These routines allocate buffer spaces used in env_reduce_base_real_gpu.
! These should be in env_data_buffer.f90 but need to be here becouse size is
! depends on the __MSGSIZ_MAX definition
SUBROUTINE env_allocate_buffers
    USE env_data_buffer
    IMPLICIT NONE
    INTEGER, PARAMETER :: maxb = __MSGSIZ_MAX
    !
    IF (.NOT. ALLOCATED(mp_buff_r)) ALLOCATE(mp_buff_r(maxb))
    IF (.NOT. ALLOCATED(mp_buff_i)) ALLOCATE(mp_buff_i(maxb))
    !
END SUBROUTINE env_allocate_buffers

SUBROUTINE env_deallocate_buffers
    USE env_data_buffer
    IMPLICIT NONE
    !
    IF (ALLOCATED (mp_buff_r) .AND. ALLOCATED (mp_buff_i)) DEALLOCATE(mp_buff_r, mp_buff_i)
    !
END SUBROUTINE env_deallocate_buffers

!=----------------------------------------------------------------------------=!
!

SUBROUTINE env_mp_synchronize( gid )
   USE env_parallel_include  
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: gid
#if defined __MPI && defined __USE_BARRIER
   INTEGER :: ierr
   CALL mpi_barrier( gid, ierr )
   IF( ierr /= 0 ) CALL env_errore( 'env_mp_synchronize ', ' error in mpi_barrier ', ierr )
#endif
   RETURN
END SUBROUTINE env_mp_synchronize


!=----------------------------------------------------------------------------=!
!

   SUBROUTINE env_bcast_real( array, n, root, gid )
        USE env_util_param, ONLY: DP
        USE env_parallel_include  
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: n, root, gid
        REAL(DP) :: array( n )
#if defined __MPI
        INTEGER :: msgsiz_max = __BCAST_MSGSIZ_MAX
        INTEGER :: nblk, blksiz, iblk, istart, ierr

#if defined __TRACE
        write(*,*) 'env_bcast_real IN'
#endif
        IF( n <= 0 ) GO TO 1

#if defined __USE_BARRIER
        CALL env_mp_synchronize( gid )
#endif

        IF( n <= msgsiz_max ) THEN
           CALL MPI_BCAST( array, n, MPI_DOUBLE_PRECISION, root, gid, ierr )
           IF( ierr /= 0 ) CALL env_errore( ' env_bcast_real ', ' error in mpi_bcast 1 ', ierr )
        ELSE
           nblk   = n / msgsiz_max
           blksiz = msgsiz_max
           DO iblk = 1, nblk
              istart = (iblk-1)*msgsiz_max + 1
              CALL MPI_BCAST( array( istart ), blksiz, MPI_DOUBLE_PRECISION, root, gid, ierr )
              IF( ierr /= 0 ) CALL env_errore( ' env_bcast_real ', ' error in mpi_bcast 2 ', ierr )
           END DO
           blksiz = MOD( n, msgsiz_max )
           IF( blksiz > 0 ) THEN
              istart = nblk * msgsiz_max + 1
              CALL MPI_BCAST( array( istart ), blksiz, MPI_DOUBLE_PRECISION, root, gid, ierr )
              IF( ierr /= 0 ) CALL env_errore( ' env_bcast_real ', ' error in mpi_bcast 3 ', ierr )
           END IF
        END IF

1       CONTINUE
#if defined __TRACE
        write(*,*) 'env_bcast_real OUT'
#endif

#endif

        RETURN
   END SUBROUTINE env_bcast_real 

   !------------------------------------------------------------------------------!
   SUBROUTINE env_bcast_integer(array, n, root, gid)
   !------------------------------------------------------------------------------!
   !! 
   !! Broadcast integers
   !!   
   USE env_parallel_include  
   ! 
   IMPLICIT NONE
   INTEGER, INTENT(in) :: n, root, gid
   INTEGER :: array(n)
#if defined __MPI
   INTEGER :: msgsiz_max = __MSGSIZ_MAX
   INTEGER :: nblk, blksiz, iblk, istart, ierr
#if defined __TRACE
   WRITE(*, *) 'env_bcast_integer IN'
#endif
   IF(n <= 0) GOTO 1
   ! 
#if defined __USE_BARRIER
   CALL env_mp_synchronize(gid)
#endif
   ! 
   IF(n <= msgsiz_max) THEN
     CALL MPI_BCAST(array, n, MPI_INTEGER, root, gid, ierr)
     IF(ierr /= 0) CALL env_errore(' env_bcast_integer ', ' error in mpi_bcast 1 ', ierr)
   ELSE
     nblk   = n / msgsiz_max
     blksiz = msgsiz_max
     DO iblk = 1, nblk
       istart = (iblk - 1) * msgsiz_max + 1
       CALL MPI_BCAST(array(istart), blksiz, MPI_INTEGER, root, gid, ierr )
       IF(ierr /= 0) CALL env_errore(' env_bcast_integer ', ' error in mpi_bcast 2 ', ierr)
     ENDDO
     blksiz = MOD( n, msgsiz_max )
     IF(blksiz > 0) THEN
       istart = nblk * msgsiz_max + 1
       CALL MPI_BCAST(array( istart ), blksiz, MPI_INTEGER, root, gid, ierr)
       IF(ierr /= 0) CALL env_errore(' env_bcast_integer ', ' error in mpi_bcast 3 ', ierr)
     ENDIF
   ENDIF
1  CONTINUE
#if defined __TRACE
   WRITE(*, *) 'env_bcast_integer OUT'
#endif
#endif
   RETURN
   !------------------------------------------------------------------------------!
   END SUBROUTINE env_bcast_integer
   !------------------------------------------------------------------------------!
   ! 
   !------------------------------------------------------------------------------!
   SUBROUTINE env_bcast_integer8(array, n, root, gid)
   !------------------------------------------------------------------------------!
   !! 
   !! Broadcast integers
   !!   
   USE env_util_param,     ONLY : i8b
   USE env_parallel_include
   ! 
   IMPLICIT NONE
   INTEGER, INTENT(in) :: n, root, gid
   INTEGER(KIND = i8b) :: array(n)
#if defined __MPI
   INTEGER :: msgsiz_max = __MSGSIZ_MAX
   INTEGER :: nblk, blksiz, iblk, istart, ierr
#if defined __TRACE
   WRITE(*, *) 'env_bcast_integer IN'
#endif
   IF(n <= 0) GOTO 1
   ! 
#if defined __USE_BARRIER
   CALL env_mp_synchronize(gid)
#endif
   ! 
   IF(n <= msgsiz_max) THEN
     CALL MPI_BCAST(array, n, MPI_INTEGER8, root, gid, ierr)
     IF(ierr /= 0) CALL env_errore(' env_bcast_integer8 ', ' error in mpi_bcast 1 ', ierr)
   ELSE
     nblk   = n / msgsiz_max
     blksiz = msgsiz_max
     DO iblk = 1, nblk
       istart = (iblk - 1) * msgsiz_max + 1
       CALL MPI_BCAST(array(istart), blksiz, MPI_INTEGER8, root, gid, ierr )
       IF(ierr /= 0) CALL env_errore(' env_bcast_integer8 ', ' error in mpi_bcast 2 ', ierr)
     ENDDO
     blksiz = MOD( n, msgsiz_max )
     IF(blksiz > 0) THEN
       istart = nblk * msgsiz_max + 1
       CALL MPI_BCAST(array( istart ), blksiz, MPI_INTEGER8, root, gid, ierr)
       IF(ierr /= 0) CALL env_errore(' env_bcast_integer8 ', ' error in mpi_bcast 3 ', ierr)
     ENDIF
   ENDIF
1  CONTINUE
#if defined __TRACE
   WRITE(*, *) 'env_bcast_integer OUT'
#endif
#endif
   RETURN
   !------------------------------------------------------------------------------!
   END SUBROUTINE env_bcast_integer8
   !------------------------------------------------------------------------------!
   ! 
   !------------------------------------------------------------------------------!
   SUBROUTINE env_bcast_logical( array, n, root, gid )
        USE env_parallel_include  
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: n, root, gid
        LOGICAL :: array( n )
#if defined __MPI
        INTEGER :: msgsiz_max = __MSGSIZ_MAX
        INTEGER :: nblk, blksiz, iblk, istart, ierr

#if defined __TRACE
        write(*,*) 'env_bcast_logical IN'
#endif

        IF( n <= 0 ) GO TO 1

#if defined __USE_BARRIER
        CALL env_mp_synchronize( gid )
#endif

        IF( n <= msgsiz_max ) THEN
           CALL MPI_BCAST( array, n, MPI_LOGICAL, root, gid, ierr )
           IF( ierr /= 0 ) CALL env_errore( ' env_bcast_logical ', ' error in mpi_bcast 1 ', ierr )
        ELSE
           nblk   = n / msgsiz_max
           blksiz = msgsiz_max
           DO iblk = 1, nblk
              istart = (iblk-1)*msgsiz_max + 1
              CALL MPI_BCAST( array( istart ), blksiz, MPI_LOGICAL, root, gid, ierr )
              IF( ierr /= 0 ) CALL env_errore( ' env_bcast_logical ', ' error in mpi_bcast 2 ', ierr )
           END DO
           blksiz = MOD( n, msgsiz_max )
           IF( blksiz > 0 ) THEN
              istart = nblk * msgsiz_max + 1
              CALL MPI_BCAST( array( istart ), blksiz, MPI_LOGICAL, root, gid, ierr )
              IF( ierr /= 0 ) CALL env_errore( ' env_bcast_logical ', ' error in mpi_bcast 3 ', ierr )
           END IF
        END IF

1       CONTINUE
#if defined __TRACE
        write(*,*) 'env_bcast_logical OUT'
#endif
#endif
        RETURN
   END SUBROUTINE env_bcast_logical


!
! ... "reduce"-like subroutines
!
#if defined (__USE_INPLACE_MPI)
!
!----------------------------------------------------------------------------
SUBROUTINE env_reduce_base_real( dim, ps, comm, root )
  !----------------------------------------------------------------------------
  !
  ! ... sums a distributed variable ps(dim) over the processors.
  ! ... This version uses a fixed-length buffer of appropriate (?) dim
  !
  USE env_util_param, ONLY : DP
  USE env_parallel_include  
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN)    :: dim     ! size of the array
  REAL(DP)                :: ps(dim) ! array whose elements have to be reduced
  INTEGER,  INTENT(IN)    :: comm    ! communicator
  INTEGER,  INTENT(IN)    :: root    ! if root <  0 perform a reduction to all procs
                                     ! if root >= 0 perform a reduce only to root proc.
  !
#if defined (__MPI)  
  !
  INTEGER            :: info, myid
  !
#if defined __TRACE
  write(*,*) 'env_reduce_base_real IN'
#endif
  !
  IF ( dim <= 0 ) GO TO 1  ! go to the end of the subroutine
  !
  ! ... synchronize processes
  !
#if defined __USE_BARRIER
  CALL env_mp_synchronize( comm )
#endif
  !
  IF( root >= 0 ) THEN
     CALL mpi_comm_rank( comm, myid, info )
     IF( info /= 0 ) CALL env_errore( 'env_reduce_base_real', 'error in mpi_comm_rank', info )
     IF (myid == root) THEN
        CALL MPI_REDUCE( MPI_IN_PLACE, ps, dim, MPI_DOUBLE_PRECISION, MPI_SUM, root, comm, info )
        IF( info /= 0 ) CALL env_errore( 'env_reduce_base_real', 'error in mpi_reduce 1', info )
     ELSE
        CALL MPI_REDUCE( ps, ps, dim, MPI_DOUBLE_PRECISION, MPI_SUM, root, comm, info )
        IF( info /= 0 ) CALL env_errore( 'env_reduce_base_real', 'error in mpi_reduce 1', info )
     ENDIF
  ELSE
     CALL MPI_ALLREDUCE( MPI_IN_PLACE, ps, dim, MPI_DOUBLE_PRECISION, MPI_SUM, comm, info )
     IF( info /= 0 ) CALL env_errore( 'env_reduce_base_real', 'error in mpi_allreduce 1', info )
  END IF
  !
1 CONTINUE
  !
#if defined __TRACE
  write(*,*) 'env_reduce_base_real OUT'
#endif
  !
#endif
  !
  RETURN
  !
END SUBROUTINE env_reduce_base_real
!
#else
!
!----------------------------------------------------------------------------
SUBROUTINE env_reduce_base_real( dim, ps, comm, root )
  !----------------------------------------------------------------------------
  !
  ! ... sums a distributed variable ps(dim) over the processors.
  ! ... This version uses a fixed-length buffer of appropriate (?) dim
  !
  USE env_util_param, ONLY : DP
  USE env_data_buffer, ONLY: buff => mp_buff_r
  USE env_parallel_include  
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN)    :: dim     ! size of the array
  REAL(DP)                :: ps(dim) ! array whose elements have to be reduced
  INTEGER,  INTENT(IN)    :: comm    ! communicator
  INTEGER,  INTENT(IN)    :: root    ! if root <  0 perform a reduction to all procs
                                     ! if root >= 0 perform a reduce only to root proc.
  !
#if defined (__MPI)  
  !
  INTEGER            :: info, n, nbuf, nproc, myid
  INTEGER, PARAMETER :: maxb = __MSGSIZ_MAX
  !
  !
#if defined __TRACE
  write(*,*) 'env_reduce_base_real IN'
#endif

  CALL mpi_comm_size( comm, nproc, info )
  IF( info /= 0 ) CALL env_errore( 'env_reduce_base_real', 'error in mpi_comm_size', info )

  CALL mpi_comm_rank( comm, myid, info )
  IF( info /= 0 ) CALL env_errore( 'env_reduce_base_real', 'error in mpi_comm_rank', info )
  !
  IF ( dim <= 0 .OR. nproc <= 1 ) GO TO 1  ! go to the end of the subroutine
  !
  ! ... synchronize processes
  !
#if defined __USE_BARRIER
  CALL env_mp_synchronize( comm )
#endif
  !
  nbuf = dim / maxb
  !
  DO n = 1, nbuf
     !
     IF( root >= 0 ) THEN
        CALL MPI_REDUCE( ps(1+(n-1)*maxb), buff, maxb, MPI_DOUBLE_PRECISION, MPI_SUM, root, comm, info )
        IF( info /= 0 ) CALL env_errore( 'env_reduce_base_real', 'error in mpi_reduce 1', info )
     ELSE
        CALL MPI_ALLREDUCE( ps(1+(n-1)*maxb), buff, maxb, MPI_DOUBLE_PRECISION, MPI_SUM, comm, info )
        IF( info /= 0 ) CALL env_errore( 'env_reduce_base_real', 'error in mpi_allreduce 1', info )
     END IF
     !                    
     IF( root < 0 ) THEN
        ps((1+(n-1)*maxb):(n*maxb)) = buff(1:maxb)
     ELSE IF( root == myid ) THEN
        ps((1+(n-1)*maxb):(n*maxb)) = buff(1:maxb)
     END IF
     !
  END DO
  !
  ! ... possible remaining elements < maxb
  !
  IF ( ( dim - nbuf * maxb ) > 0 ) THEN
     !
     IF( root >= 0 ) THEN
        CALL MPI_REDUCE( ps(1+nbuf*maxb), buff, (dim-nbuf*maxb), MPI_DOUBLE_PRECISION, MPI_SUM, root, comm, info )
        IF( info /= 0 ) CALL env_errore( 'env_reduce_base_real', 'error in mpi_reduce 2', info )
     ELSE
        CALL MPI_ALLREDUCE( ps(1+nbuf*maxb), buff, (dim-nbuf*maxb), MPI_DOUBLE_PRECISION, MPI_SUM, comm, info )
        IF( info /= 0 ) CALL env_errore( 'env_reduce_base_real', 'error in mpi_allreduce 2', info )
     END IF
     !
     IF( root < 0 ) THEN
        ps((1+nbuf*maxb):dim) = buff(1:(dim-nbuf*maxb))
     ELSE IF( root == myid ) THEN
        ps((1+nbuf*maxb):dim) = buff(1:(dim-nbuf*maxb))
     END IF
     !
  END IF
  !
1 CONTINUE
  !
#if defined __TRACE
  write(*,*) 'env_reduce_base_real OUT'
#endif
  !
#endif
  !
  RETURN
  !
END SUBROUTINE env_reduce_base_real
!
#endif
!
!
#if defined (__USE_INPLACE_MPI)
  !
  !----------------------------------------------------------------------------
  SUBROUTINE env_reduce_base_integer(dim, ps, comm, root)
  !----------------------------------------------------------------------------
  !!
  !! Sums a distributed variable ps(dim) over the processors.
  !! This version uses a fixed-length buffer of appropriate (?) dim
  !!
  USE env_util_param,       ONLY : DP
  USE env_parallel_include  
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(in)    :: dim
  INTEGER, INTENT(inout) :: ps(dim)
  INTEGER, INTENT(in)    :: comm    ! communicator
  INTEGER, INTENT(in)    :: root    ! if root <  0 perform a reduction to all procs
                                    ! if root >= 0 perform a reduce only to root proc.
  !
#if defined (__MPI)  
  !
  INTEGER            :: info, myid
  !
#if defined __TRACE
  WRITE(*, *) 'env_reduce_base_integer IN'
#endif
  !
  ! ... synchronize processes
  !
#if defined __USE_BARRIER
  CALL env_mp_synchronize(comm)
#endif
  !
  IF(root >= 0) THEN
    CALL mpi_comm_rank( comm, myid, info )
    IF( info /= 0 ) CALL env_errore( 'env_reduce_base_integer', 'error in mpi_comm_rank', info )
    IF (myid == root) THEN
      CALL MPI_REDUCE(MPI_IN_PLACE, ps, dim, MPI_INTEGER, MPI_SUM, root, comm, info)
      IF(info /= 0) CALL env_errore('env_reduce_base_integer', 'error in mpi_reduce 1', info)
    ELSE
      CALL MPI_REDUCE(ps, ps, dim, MPI_INTEGER, MPI_SUM, root, comm, info)
      IF(info /= 0) CALL env_errore('env_reduce_base_integer', 'error in mpi_reduce 1', info)
    ENDIF
  ELSE
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, ps, dim, MPI_INTEGER, MPI_SUM, comm, info)
    IF(info /= 0 ) CALL env_errore('env_reduce_base_integer', 'error in mpi_allreduce 1', info)
  ENDIF
  !
#if defined __TRACE
  WRITE(*, *) 'env_reduce_base_integer OUT'
#endif
  !
#endif
  !
  RETURN
  !----------------------------------------------------------------------------
  END SUBROUTINE env_reduce_base_integer
  !----------------------------------------------------------------------------
  !
  !----------------------------------------------------------------------------
  SUBROUTINE env_reduce_base_integer8(dim, ps, comm, root)
  !----------------------------------------------------------------------------
  !!
  !! Sums a distributed variable ps(dim) over the processors.
  !! This version uses a fixed-length buffer of appropriate (?) dim
  !!
  USE env_util_param,       ONLY : DP, i8b
  USE env_parallel_include
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(in)  :: dim
  INTEGER(KIND = i8b), INTENT(inout) :: ps(dim)
  INTEGER, INTENT(in)  :: comm    ! communicator
  INTEGER, INTENT(in)  :: root    ! if root <  0 perform a reduction to all procs
                                  ! if root >= 0 perform a reduce only to root proc.
  !   
#if defined (__MPI)  
  !
  INTEGER            :: info, myid
  !
#if defined __TRACE
  WRITE(*, *) 'env_reduce_base_integer IN'
#endif
  !
  ! ... synchronize processes
  !
#if defined __USE_BARRIER
  CALL env_mp_synchronize(comm)
#endif
  !
  IF(root >= 0) THEN
    CALL mpi_comm_rank( comm, myid, info )
    IF( info /= 0 ) CALL env_errore( 'env_reduce_base_integer8', 'error in mpi_comm_rank', info )
    IF (myid == root) THEN
      CALL MPI_REDUCE(MPI_IN_PLACE, ps, dim, MPI_INTEGER8, MPI_SUM, root, comm, info)
      IF(info /= 0) CALL env_errore('env_reduce_base_integer8', 'error in mpi_reduce 1', info)
    ELSE
      CALL MPI_REDUCE(ps, ps, dim, MPI_INTEGER8, MPI_SUM, root, comm, info)
      IF(info /= 0) CALL env_errore('env_reduce_base_integer8', 'error in mpi_reduce 1', info)
    ENDIF
  ELSE
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, ps, dim, MPI_INTEGER8, MPI_SUM, comm, info)
    IF(info /= 0 ) CALL env_errore('env_reduce_base_integer8', 'error in mpi_allreduce 1', info)
  ENDIF
  !
#if defined __TRACE
  WRITE(*, *) 'env_reduce_base_integer OUT'
#endif
  !
#endif
  !
  RETURN
  !----------------------------------------------------------------------------
  END SUBROUTINE env_reduce_base_integer8
  !----------------------------------------------------------------------------
#else
  !
  !----------------------------------------------------------------------------
  SUBROUTINE env_reduce_base_integer(dim, ps, comm, root)
  !----------------------------------------------------------------------------
  !!
  !! Sums a distributed variable ps(dim) over the processors.
  !! This version uses a fixed-length buffer of appropriate (?) dim
  !!
  USE env_util_param,  ONLY : DP
  USE env_data_buffer, ONLY : buff => mp_buff_i
  USE env_parallel_include  
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(in)    :: dim
  INTEGER, INTENT(inout) :: ps(dim)
  INTEGER, INTENT(in)    :: comm    ! communicator
  INTEGER, INTENT(in)    :: root    ! if root <  0 perform a reduction to all procs
                                    ! if root >= 0 perform a reduce only to root proc.
  !
#if defined (__MPI)  
  !
  INTEGER            :: info, n, nbuf, nproc, myid
  INTEGER, PARAMETER :: maxb = __MSGSIZ_MAX
  !
#if defined __TRACE
  WRITE(*, *) 'env_reduce_base_integer IN'
#endif
  !
  CALL mpi_comm_size(comm, nproc, info)
  IF(info /= 0) CALL env_errore('env_reduce_base_integer', 'error in mpi_comm_size', info)
  ! 
  CALL mpi_comm_rank(comm, myid, info)
  IF(info /= 0) CALL env_errore('env_reduce_base_integer', 'error in mpi_comm_rank', info)
  !
  IF (dim <= 0 .OR. nproc <= 1) GOTO 1  ! go to the end of the subroutine
  !
  ! ... synchronize processes
  !
#if defined __USE_BARRIER
  CALL env_mp_synchronize(comm)
#endif
  !
  nbuf = dim / maxb
  !
  DO n = 1, nbuf
    !
    IF (root >= 0) THEN
      CALL MPI_REDUCE(ps(1 + (n - 1) * maxb), buff, maxb, MPI_INTEGER, MPI_SUM, root, comm, info)
      IF (info /= 0 ) CALL env_errore('env_reduce_base_integer', 'error in mpi_reduce 1', info)
    ELSE
      CALL MPI_ALLREDUCE(ps(1 + (n - 1) * maxb), buff, maxb, MPI_INTEGER, MPI_SUM, comm, info)
      IF (info /= 0) CALL env_errore('env_reduce_base_integer', 'error in mpi_allreduce 1', info)
    ENDIF
    !
    IF (root < 0 ) THEN
      ps((1 + (n - 1) * maxb):(n * maxb)) = buff(1:maxb)
    ELSEIF( root == myid) THEN
      ps((1 + (n - 1) * maxb):(n * maxb)) = buff(1:maxb)
    ENDIF
    !
  ENDDO
  !
  ! ... possible remaining elements < maxb
  !
  IF ((dim - nbuf * maxb) > 0) THEN
    !
    IF (root >= 0) THEN
      CALL MPI_REDUCE(ps(1 + nbuf * maxb), buff, (dim - nbuf * maxb), MPI_INTEGER, MPI_SUM, root, comm, info)
      IF(info /= 0) CALL env_errore('env_reduce_base_integer', 'error in mpi_reduce 2', info)
    ELSE
      CALL MPI_ALLREDUCE(ps(1 + nbuf * maxb), buff, (dim - nbuf * maxb), MPI_INTEGER, MPI_SUM, comm, info)
      IF(info /= 0) CALL env_errore('env_reduce_base_integer', 'error in mpi_allreduce 2', info)
    ENDIF
    !
    IF(root < 0) THEN
      ps((1 + nbuf * maxb):dim) = buff(1:(dim - nbuf * maxb))
    ELSEIF(root == myid) THEN
      ps((1 + nbuf * maxb):dim) = buff(1:(dim - nbuf * maxb))
    ENDIF
    !
  ENDIF
  !
1 CONTINUE
  !
#if defined __TRACE
  WRITE(*, *) 'env_reduce_base_integer OUT'
#endif
  !
#endif
  !
  RETURN
  !----------------------------------------------------------------------------
  END SUBROUTINE env_reduce_base_integer
  !----------------------------------------------------------------------------
  !
  !----------------------------------------------------------------------------
  SUBROUTINE env_reduce_base_integer8(dim, ps, comm, root)
  !----------------------------------------------------------------------------
  !!
  !! Sums a distributed variable ps(dim) over the processors.
  !! This version uses a fixed-length buffer of appropriate (?) dim
  !!
  USE env_util_param,  ONLY : DP, i8b
  USE env_data_buffer, ONLY : buff => mp_buff_i
  USE env_parallel_include  
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(in)    :: dim
  INTEGER(KIND = i8b), INTENT(inout) :: ps(dim)
  INTEGER, INTENT(in)    :: comm    ! communicator
  INTEGER, INTENT(in)    :: root    ! if root <  0 perform a reduction to all procs
                                    ! if root >= 0 perform a reduce only to root proc.
  !
#if defined (__MPI)  
  !
  INTEGER            :: info, n, nbuf, nproc, myid
  INTEGER, PARAMETER :: maxb = __MSGSIZ_MAX
  !
#if defined __TRACE
  WRITE(*, *) 'env_reduce_base_integer8 IN'
#endif
  !
  CALL mpi_comm_size(comm, nproc, info)
  IF(info /= 0) CALL env_errore('env_reduce_base_integer8', 'error in mpi_comm_size', info)
  ! 
  CALL mpi_comm_rank(comm, myid, info)
  IF(info /= 0) CALL env_errore('env_reduce_base_integer8', 'error in mpi_comm_rank', info)
  !
  IF (dim <= 0 .OR. nproc <= 1) GOTO 1  ! go to the end of the subroutine
  !
  ! ... synchronize processes
  !
#if defined __USE_BARRIER
  CALL env_mp_synchronize(comm)
#endif
  !
  nbuf = dim / maxb
  !
  DO n = 1, nbuf
    !
    IF (root >= 0) THEN
      CALL MPI_REDUCE(ps(1 + (n - 1) * maxb), buff, maxb, MPI_INTEGER8, MPI_SUM, root, comm, info)
      IF (info /= 0 ) CALL env_errore('env_reduce_base_integer8', 'error in mpi_reduce 1', info)
    ELSE
      CALL MPI_ALLREDUCE(ps(1 + (n - 1) * maxb), buff, maxb, MPI_INTEGER8, MPI_SUM, comm, info)
      IF (info /= 0) CALL env_errore('env_reduce_base_integer8', 'error in mpi_allreduce 1', info)
    ENDIF
    !
    IF (root < 0 ) THEN
      ps((1 + (n - 1) * maxb):(n * maxb)) = buff(1:maxb)
    ELSEIF( root == myid) THEN
      ps((1 + (n - 1) * maxb):(n * maxb)) = buff(1:maxb)
    ENDIF
    !
  ENDDO
  !
  ! ... possible remaining elements < maxb
  !
  IF ((dim - nbuf * maxb) > 0) THEN
    !
    IF (root >= 0) THEN
      CALL MPI_REDUCE(ps(1 + nbuf * maxb), buff, (dim - nbuf * maxb), MPI_INTEGER8, MPI_SUM, root, comm, info)
      IF(info /= 0) CALL env_errore('env_reduce_base_integer', 'error in mpi_reduce 2', info)
    ELSE
      CALL MPI_ALLREDUCE(ps(1 + nbuf * maxb), buff, (dim - nbuf * maxb), MPI_INTEGER8, MPI_SUM, comm, info)
      IF(info /= 0) CALL env_errore('env_reduce_base_integer', 'error in mpi_allreduce 2', info)
    ENDIF
    !
    IF(root < 0) THEN
      ps((1 + nbuf * maxb):dim) = buff(1:(dim - nbuf * maxb))
    ELSEIF(root == myid) THEN
      ps((1 + nbuf * maxb):dim) = buff(1:(dim - nbuf * maxb))
    ENDIF
    !
  ENDIF
  !
1 CONTINUE
  !
#if defined __TRACE
  WRITE(*, *) 'env_reduce_base_integer8 OUT'
#endif
  !
#endif
  !
  RETURN
  !----------------------------------------------------------------------------
  END SUBROUTINE env_reduce_base_integer8
  !----------------------------------------------------------------------------
  !
#endif
  !
  ! ... "reduce"-like subroutines
  !
  !----------------------------------------------------------------------------
  SUBROUTINE env_reduce_base_real_to( dim, ps, psout, comm, root )
  !----------------------------------------------------------------------------
  !
  ! ... sums a distributed variable ps(dim) over the processors,
  ! ... and store the results in variable psout.
  ! ... This version uses a fixed-length buffer of appropriate (?) length
  !
  USE env_util_param, ONLY : DP
  USE env_parallel_include  
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN)  :: dim
  REAL(DP), INTENT(IN)  :: ps(dim)
  REAL(DP)              :: psout(dim)
  INTEGER,  INTENT(IN)  :: comm    ! communecator
  INTEGER,  INTENT(IN)  :: root    ! if root <  0 perform a reduction to all procs
                                   ! if root >= 0 perform a reduce only to root proc.
  !
#if defined (__MPI)  
  !
  INTEGER            :: info, n, nbuf, nproc, myid
  INTEGER, PARAMETER :: maxb = __MSGSIZ_MAX
  !
#if defined __TRACE
  write(*,*) 'env_reduce_base_real_to IN'
#endif

  CALL mpi_comm_size( comm, nproc, info )
  IF( info /= 0 ) CALL env_errore( 'env_reduce_base_real_to', 'error in mpi_comm_size', info )

  CALL mpi_comm_rank( comm, myid, info )
  IF( info /= 0 ) CALL env_errore( 'env_reduce_base_real_to', 'error in mpi_comm_rank', info )
  !
  IF ( dim > 0 .AND. nproc <= 1 ) THEN
     psout = ps
  END IF
  IF( dim <= 0 .OR. nproc <= 1 ) GO TO 1 ! go to the end of the subroutine
  !
  ! ... synchronize processes
  !
#if defined __USE_BARRIER
  CALL env_mp_synchronize( comm )
#endif
  !
  nbuf = dim / maxb
  !
  DO n = 1, nbuf
     !
     IF( root >= 0 ) THEN
        CALL MPI_REDUCE( ps(1+(n-1)*maxb), psout(1+(n-1)*maxb), maxb, MPI_DOUBLE_PRECISION, MPI_SUM, root, comm, info )
        IF( info /= 0 ) CALL env_errore( 'env_reduce_base_real_to', 'error in mpi_reduce 1', info )
     ELSE
        CALL MPI_ALLREDUCE( ps(1+(n-1)*maxb), psout(1+(n-1)*maxb), maxb, MPI_DOUBLE_PRECISION, MPI_SUM, comm, info )
        IF( info /= 0 ) CALL env_errore( 'env_reduce_base_real_to', 'error in mpi_allreduce 1', info )
     END IF
     !                    
  END DO
  !
  ! ... possible remaining elements < maxb
  !
  IF ( ( dim - nbuf * maxb ) > 0 ) THEN
     !
     IF( root >= 0 ) THEN
        CALL MPI_REDUCE( ps(1+nbuf*maxb), psout(1+nbuf*maxb), (dim-nbuf*maxb), MPI_DOUBLE_PRECISION, MPI_SUM, root, comm, info )
        IF( info /= 0 ) CALL env_errore( 'env_reduce_base_real_to', 'error in mpi_reduce 2', info )
     ELSE
        CALL MPI_ALLREDUCE( ps(1+nbuf*maxb), psout(1+nbuf*maxb), (dim-nbuf*maxb), MPI_DOUBLE_PRECISION, MPI_SUM, comm, info )
        IF( info /= 0 ) CALL env_errore( 'env_reduce_base_real_to', 'error in mpi_allreduce 2', info )
     END IF
     !
  END IF
  !
1 CONTINUE
  !
#if defined __TRACE
  write(*,*) 'env_reduce_base_real_to OUT'
#endif
  !
#endif
  !
  RETURN
  !
END SUBROUTINE env_reduce_base_real_to
!
!
!
!----------------------------------------------------------------------------
SUBROUTINE env_reduce_base_integer_to( dim, ps, psout, comm, root )
  !----------------------------------------------------------------------------
  !
  ! ... sums a distributed integer variable ps(dim) over the processors, and
  ! ... saves the result on the output variable psout.
  ! ... This version uses a fixed-length buffer of appropriate (?) length
  !
  USE env_util_param, ONLY : DP
  USE env_parallel_include  
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN)  :: dim
  INTEGER,  INTENT(IN)  :: ps(dim)
  INTEGER               :: psout(dim)
  INTEGER,  INTENT(IN)  :: comm    ! communecator
  INTEGER,  INTENT(IN)  :: root    ! if root <  0 perform a reduction to all procs
                                     ! if root >= 0 perform a reduce only to root proc.
  !
#if defined (__MPI)  
  !
  INTEGER            :: info, n, nbuf, nproc, myid
  INTEGER, PARAMETER :: maxb = __MSGSIZ_MAX
  !
#if defined __TRACE
  write(*,*) 'env_reduce_base_integer_to IN'
#endif

  CALL mpi_comm_size( comm, nproc, info )
  IF( info /= 0 ) CALL env_errore( 'env_reduce_base_integer_to', 'error in mpi_comm_size', info )

  CALL mpi_comm_rank( comm, myid, info )
  IF( info /= 0 ) CALL env_errore( 'env_reduce_base_integer_to', 'error in mpi_comm_rank', info )
  !
  IF ( dim > 0 .AND. nproc <= 1 ) THEN
     psout = ps
  END IF
  IF( dim <= 0 .OR. nproc <= 1 ) GO TO 1 ! go to the end of the subroutine
  !
  ! ... synchronize processes
  !
#if defined __USE_BARRIER
  CALL env_mp_synchronize( comm )
#endif
  !
  nbuf = dim / maxb
  !
  DO n = 1, nbuf
     !
     IF( root >= 0 ) THEN
        CALL MPI_REDUCE( ps(1+(n-1)*maxb), psout( 1+(n-1)*maxb ), maxb, MPI_INTEGER, MPI_SUM, root, comm, info )
        IF( info /= 0 ) CALL env_errore( 'env_reduce_base_integer_to', 'error in mpi_reduce 1', info )
     ELSE
        CALL MPI_ALLREDUCE( ps(1+(n-1)*maxb), psout( 1+(n-1)*maxb ), maxb, MPI_INTEGER, MPI_SUM, comm, info )
        IF( info /= 0 ) CALL env_errore( 'env_reduce_base_integer_to', 'error in mpi_allreduce 1', info )
     END IF
     !                    
  END DO
  !
  ! ... possible remaining elements < maxb
  !
  IF ( ( dim - nbuf * maxb ) > 0 ) THEN
     !
     IF( root >= 0 ) THEN
        CALL MPI_REDUCE( ps(1+nbuf*maxb), psout(1+nbuf*maxb), (dim-nbuf*maxb), MPI_INTEGER, MPI_SUM, root, comm, info )
        IF( info /= 0 ) CALL env_errore( 'env_reduce_base_integer_to', 'error in mpi_reduce 2', info )
     ELSE
        CALL MPI_ALLREDUCE( ps(1+nbuf*maxb), psout(1+nbuf*maxb), (dim-nbuf*maxb), MPI_INTEGER, MPI_SUM, comm, info )
        IF( info /= 0 ) CALL env_errore( 'env_reduce_base_integer_to', 'error in mpi_allreduce 2', info )
     END IF
     !
  END IF
  !
1 CONTINUE
  !
#if defined __TRACE
  write(*,*) 'env_reduce_base_integer_to OUT'
#endif
  !
#endif
  !
  RETURN
  !
END SUBROUTINE env_reduce_base_integer_to
!
!
!  Parallel MIN and MAX
!

!----------------------------------------------------------------------------
SUBROUTINE env_parallel_min_integer( dim, ps, comm, root )
  !----------------------------------------------------------------------------
  !
  ! ... compute the minimum of a distributed variable ps(dim) over the processors.
  ! ... This version uses a fixed-length buffer of appropriate (?) dim
  !
  USE env_util_param, ONLY : DP
  USE env_data_buffer, ONLY : buff => mp_buff_i
  USE env_parallel_include  
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN)    :: dim
  INTEGER                 :: ps(dim)
  INTEGER,  INTENT(IN)    :: comm    ! communecator
  INTEGER,  INTENT(IN)    :: root    ! if root <  0 perform a reduction to all procs
                                     ! if root >= 0 perform a reduce only to root proc.
  !
#if defined (__MPI)  
  !
  INTEGER            :: info, n, nbuf, nproc, myid
  INTEGER, PARAMETER :: maxb = __MSGSIZ_MAX
  !
#if defined __TRACE
  write(*,*) 'env_parallel_min_integer IN'
#endif
  !
  CALL mpi_comm_size( comm, nproc, info )
  IF( info /= 0 ) CALL env_errore( 'env_parallel_min_integer', 'error in mpi_comm_size', info )

  CALL mpi_comm_rank( comm, myid, info )
  IF( info /= 0 ) CALL env_errore( 'env_parallel_min_integer', 'error in mpi_comm_rank', info )
  !
  IF ( dim <= 0 .OR. nproc <= 1 ) GO TO 1
  !
  ! ... synchronize processes
  !
#if defined __USE_BARRIER
  CALL env_mp_synchronize( comm )
#endif
  !
  nbuf = dim / maxb
  !
  DO n = 1, nbuf
     !
     IF( root >= 0 ) THEN
        CALL MPI_REDUCE( ps(1+(n-1)*maxb), buff, maxb, MPI_INTEGER, MPI_MIN, root, comm, info )
        IF( info /= 0 ) CALL env_errore( 'env_parallel_min_integer', 'error in mpi_reduce 1', info )
     ELSE
        CALL MPI_ALLREDUCE( ps(1+(n-1)*maxb), buff, maxb, MPI_INTEGER, MPI_MIN, comm, info )
        IF( info /= 0 ) CALL env_errore( 'env_parallel_min_integer', 'error in mpi_allreduce 1', info )
     END IF
     !
     IF( root < 0 ) THEN
        ps((1+(n-1)*maxb):(n*maxb)) = buff(1:maxb)
     ELSE IF( root == myid ) THEN
        ps((1+(n-1)*maxb):(n*maxb)) = buff(1:maxb)
     END IF
     !
  END DO
  !
  ! ... possible remaining elements < maxb
  !
  IF ( ( dim - nbuf * maxb ) > 0 ) THEN
     !
     IF( root >= 0 ) THEN
        CALL MPI_REDUCE( ps(1+nbuf*maxb), buff, (dim-nbuf*maxb), MPI_INTEGER, MPI_MIN, root, comm, info )
        IF( info /= 0 ) CALL env_errore( 'env_parallel_min_integer', 'error in mpi_reduce 2', info )
     ELSE
        CALL MPI_ALLREDUCE( ps(1+nbuf*maxb), buff, (dim-nbuf*maxb), MPI_INTEGER, MPI_MIN, comm, info )
        IF( info /= 0 ) CALL env_errore( 'env_parallel_min_integer', 'error in mpi_allreduce 2', info )
     END IF
     !
     IF( root < 0 ) THEN
        ps((1+nbuf*maxb):dim) = buff(1:(dim-nbuf*maxb))
     ELSE IF( root == myid ) THEN
        ps((1+nbuf*maxb):dim) = buff(1:(dim-nbuf*maxb))
     END IF
     !
  END IF
  !
1 CONTINUE
  !
#if defined __TRACE
  write(*,*) 'env_parallel_min_integer OUT'
#endif
  !
#endif
  !
  RETURN
  !
END SUBROUTINE env_parallel_min_integer

!
!----------------------------------------------------------------------------
SUBROUTINE env_parallel_max_integer( dim, ps, comm, root )
  !----------------------------------------------------------------------------
  !
  ! ... compute the maximum of a distributed variable ps(dim) over the processors.
  ! ... This version uses a fixed-length buffer of appropriate (?) dim
  !
  USE env_util_param,  ONLY : DP
  USE env_data_buffer, ONLY : buff => mp_buff_i
  USE env_parallel_include  
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN)    :: dim
  INTEGER                 :: ps(dim)
  INTEGER,  INTENT(IN)    :: comm    ! communecator
  INTEGER,  INTENT(IN)    :: root    ! if root <  0 perform a reduction to all procs
                                     ! if root >= 0 perform a reduce only to root proc.
  !
#if defined (__MPI)  
  !
  INTEGER            :: info, n, nbuf, nproc, myid
  INTEGER, PARAMETER :: maxb = __MSGSIZ_MAX
  !
#if defined __TRACE
  write(*,*) 'env_parallel_max_integer IN'
#endif
  CALL mpi_comm_size( comm, nproc, info )
  IF( info /= 0 ) CALL env_errore( 'env_parallel_max_integer', 'error in mpi_comm_size', info )

  CALL mpi_comm_rank( comm, myid, info )
  IF( info /= 0 ) CALL env_errore( 'env_parallel_max_integer', 'error in mpi_comm_rank', info )
  !
  IF ( dim <= 0 .OR. nproc <= 1 ) GO TO 1
  !
  ! ... synchronize processes
  !
#if defined __USE_BARRIER
  CALL env_mp_synchronize( comm )
#endif
  !
  nbuf = dim / maxb
  !
  DO n = 1, nbuf
     !
     IF( root >= 0 ) THEN
        CALL MPI_REDUCE( ps(1+(n-1)*maxb), buff, maxb, MPI_INTEGER, MPI_MAX, root, comm, info )
        IF( info /= 0 ) CALL env_errore( 'env_parallel_max_integer', 'error in mpi_reduce 1', info )
     ELSE
        CALL MPI_ALLREDUCE( ps(1+(n-1)*maxb), buff, maxb, MPI_INTEGER, MPI_MAX, comm, info )
        IF( info /= 0 ) CALL env_errore( 'env_parallel_max_integer', 'error in mpi_allreduce 1', info )
     END IF
     !
     IF( root < 0 ) THEN
        ps((1+(n-1)*maxb):(n*maxb)) = buff(1:maxb)
     ELSE IF( root == myid ) THEN
        ps((1+(n-1)*maxb):(n*maxb)) = buff(1:maxb)
     END IF
     !
  END DO
  !
  ! ... possible remaining elements < maxb
  !
  IF ( ( dim - nbuf * maxb ) > 0 ) THEN
     !
     IF( root >= 0 ) THEN
        CALL MPI_REDUCE( ps(1+nbuf*maxb), buff, (dim-nbuf*maxb), MPI_INTEGER, MPI_MAX, root, comm, info )
        IF( info /= 0 ) CALL env_errore( 'env_parallel_max_integer', 'error in mpi_reduce 2', info )
     ELSE
        CALL MPI_ALLREDUCE( ps(1+nbuf*maxb), buff, (dim-nbuf*maxb), MPI_INTEGER, MPI_MAX, comm, info )
        IF( info /= 0 ) CALL env_errore( 'env_parallel_max_integer', 'error in mpi_allreduce 2', info )
     END IF
     !
     IF( root < 0 ) THEN
        ps((1+nbuf*maxb):dim) = buff(1:(dim-nbuf*maxb))
     ELSE IF( root == myid ) THEN
        ps((1+nbuf*maxb):dim) = buff(1:(dim-nbuf*maxb))
     END IF
     !
  END IF
  !
1 CONTINUE
  !
#if defined __TRACE
  write(*,*) 'env_parallel_max_integer OUT'
#endif
#endif
  !
  RETURN
  !
END SUBROUTINE env_parallel_max_integer


!----------------------------------------------------------------------------
SUBROUTINE env_parallel_min_real( dim, ps, comm, root )
  !----------------------------------------------------------------------------
  !
  ! ... compute the minimum value of a distributed variable ps(dim) over the processors.
  ! ... This version uses a fixed-length buffer of appropriate (?) dim
  !
  USE env_util_param, ONLY : DP
  USE env_data_buffer, ONLY : buff => mp_buff_r
  USE env_parallel_include  
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN)    :: dim
  REAL(DP)                :: ps(dim)
  INTEGER,  INTENT(IN)    :: comm    ! communecator
  INTEGER,  INTENT(IN)    :: root    ! if root <  0 perform a reduction to all procs
                                     ! if root >= 0 perform a reduce only to root proc.
  !
#if defined (__MPI)  
  !
  INTEGER            :: info, n, nbuf, nproc, myid
  INTEGER, PARAMETER :: maxb = __MSGSIZ_MAX
  !
#if defined __TRACE
  write(*,*) 'env_parallel_min_real IN'
#endif
  CALL mpi_comm_size( comm, nproc, info )
  IF( info /= 0 ) CALL env_errore( 'env_parallel_min_real', 'error in mpi_comm_size', info )

  CALL mpi_comm_rank( comm, myid, info )
  IF( info /= 0 ) CALL env_errore( 'env_parallel_min_real', 'error in mpi_comm_rank', info )
  !
  IF ( dim <= 0 .OR. nproc <= 1 ) GO TO 1
  !
  ! ... synchronize processes
  !
#if defined __USE_BARRIER
  CALL env_mp_synchronize( comm )
#endif
  !
  nbuf = dim / maxb
  !
  DO n = 1, nbuf
     !
     IF( root >= 0 ) THEN
        CALL MPI_REDUCE( ps(1+(n-1)*maxb), buff, maxb, MPI_DOUBLE_PRECISION, MPI_MIN, root, comm, info )
        IF( info /= 0 ) CALL env_errore( 'env_parallel_min_real', 'error in mpi_reduce 1', info )
     ELSE
        CALL MPI_ALLREDUCE( ps(1+(n-1)*maxb), buff, maxb, MPI_DOUBLE_PRECISION, MPI_MIN, comm, info )
        IF( info /= 0 ) CALL env_errore( 'env_parallel_min_real', 'error in mpi_allreduce 1', info )
     END IF
     !
     IF( root < 0 ) THEN
        ps((1+(n-1)*maxb):(n*maxb)) = buff(1:maxb)
     ELSE IF( root == myid ) THEN
        ps((1+(n-1)*maxb):(n*maxb)) = buff(1:maxb)
     END IF
     !
  END DO
  !
  ! ... possible remaining elements < maxb
  !
  IF ( ( dim - nbuf * maxb ) > 0 ) THEN
     !
     IF( root >= 0 ) THEN
        CALL MPI_REDUCE( ps(1+nbuf*maxb), buff, (dim-nbuf*maxb), MPI_DOUBLE_PRECISION, MPI_MIN, root, comm, info )
        IF( info /= 0 ) CALL env_errore( 'env_parallel_min_real', 'error in mpi_reduce 2', info )
     ELSE
        CALL MPI_ALLREDUCE( ps(1+nbuf*maxb), buff, (dim-nbuf*maxb), MPI_DOUBLE_PRECISION, MPI_MIN, comm, info )
        IF( info /= 0 ) CALL env_errore( 'env_parallel_min_real', 'error in mpi_allreduce 2', info )
     END IF
     !
     IF( root < 0 ) THEN
        ps((1+nbuf*maxb):dim) = buff(1:(dim-nbuf*maxb))
     ELSE IF( root == myid ) THEN
        ps((1+nbuf*maxb):dim) = buff(1:(dim-nbuf*maxb))
     END IF
     !
  END IF
  !
1 CONTINUE
  !
#if defined __TRACE
  write(*,*) 'env_parallel_min_real OUT'
#endif
#endif
  !
  RETURN
  !
END SUBROUTINE env_parallel_min_real

!
!----------------------------------------------------------------------------
SUBROUTINE env_parallel_max_real( dim, ps, comm, root )
  !----------------------------------------------------------------------------
  !
  ! ... compute the maximum value of a distributed variable ps(dim) over the processors.
  ! ... This version uses a fixed-length buffer of appropriate (?) dim
  !
  USE env_util_param, ONLY : DP
  USE env_data_buffer, ONLY : buff => mp_buff_r
  USE env_parallel_include  
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN)    :: dim
  REAL(DP)                :: ps(dim)
  INTEGER,  INTENT(IN)    :: comm    ! communecator
  INTEGER,  INTENT(IN)    :: root    ! if root <  0 perform a reduction to all procs
                                     ! if root >= 0 perform a reduce only to root proc.
  !
#if defined (__MPI)  
  !
  INTEGER            :: info, n, nbuf, nproc, myid
  INTEGER, PARAMETER :: maxb = __MSGSIZ_MAX
  !
#if defined __TRACE
  write(*,*) 'env_parallel_max_real IN'
#endif

  CALL mpi_comm_size( comm, nproc, info )
  IF( info /= 0 ) CALL env_errore( 'env_parallel_max_real', 'error in mpi_comm_size', info )

  CALL mpi_comm_rank( comm, myid, info )
  IF( info /= 0 ) CALL env_errore( 'env_parallel_max_real', 'error in mpi_comm_rank', info )
  !
  IF ( dim <= 0 .OR. nproc <= 1 ) GO TO 1
  !
  ! ... synchronize processes
  !
#if defined __USE_BARRIER
  CALL env_mp_synchronize( comm )
#endif
  !
  nbuf = dim / maxb
  !
  DO n = 1, nbuf
     !
     IF( root >= 0 ) THEN
        CALL MPI_REDUCE( ps(1+(n-1)*maxb), buff, maxb, MPI_DOUBLE_PRECISION, MPI_MAX, root, comm, info )
        IF( info /= 0 ) CALL env_errore( 'env_parallel_max_real', 'error in mpi_reduce 1', info )
     ELSE
        CALL MPI_ALLREDUCE( ps(1+(n-1)*maxb), buff, maxb, MPI_DOUBLE_PRECISION, MPI_MAX, comm, info )
        IF( info /= 0 ) CALL env_errore( 'env_parallel_max_real', 'error in mpi_allreduce 1', info )
     END IF
     !
     IF( root < 0 ) THEN
        ps((1+(n-1)*maxb):(n*maxb)) = buff(1:maxb)
     ELSE IF( root == myid ) THEN
        ps((1+(n-1)*maxb):(n*maxb)) = buff(1:maxb)
     END IF
     !
  END DO
  !
  ! ... possible remaining elements < maxb
  !
  IF ( ( dim - nbuf * maxb ) > 0 ) THEN
     !
     IF( root >= 0 ) THEN
        CALL MPI_REDUCE( ps(1+nbuf*maxb), buff, (dim-nbuf*maxb), MPI_DOUBLE_PRECISION, MPI_MAX, root, comm, info )
        IF( info /= 0 ) CALL env_errore( 'env_parallel_max_real', 'error in mpi_reduce 2', info )
     ELSE
        CALL MPI_ALLREDUCE( ps(1+nbuf*maxb), buff, (dim-nbuf*maxb), MPI_DOUBLE_PRECISION, MPI_MAX, comm, info )
        IF( info /= 0 ) CALL env_errore( 'env_parallel_max_real', 'error in mpi_allreduce 2', info )
     END IF
     !
     IF( root < 0 ) THEN
        ps((1+nbuf*maxb):dim) = buff(1:(dim-nbuf*maxb))
     ELSE IF( root == myid ) THEN
        ps((1+nbuf*maxb):dim) = buff(1:(dim-nbuf*maxb))
     END IF
     !
  END IF
  !
1 CONTINUE
  !
#if defined __TRACE
  write(*,*) 'env_parallel_max_real OUT'
#endif
  !
#endif
  !
  RETURN
  !
END SUBROUTINE env_parallel_max_real
