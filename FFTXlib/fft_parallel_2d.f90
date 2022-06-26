!
! Copyright (C) Quantum ESPRESSO group
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!=---------------------------------------------------------------------==!
!
!
!     Parallel 3D FFT high level Driver
!     ( Charge density and Wave Functions )
!
!     Written and maintained by Carlo Cavazzoni
!     Last update Apr. 2009
!
!!=---------------------------------------------------------------------==!
!
MODULE env_fft_parallel_2d
!
#ifdef __CUDA
   USE cudafor
#endif
   !
   USE env_fft_param
   IMPLICIT NONE
   SAVE
   !
!
CONTAINS
!
!  General purpose driver, including Task groups parallelization
!
!----------------------------------------------------------------------------
SUBROUTINE env_tg_cft3s( f, dfft, isgn )
  !----------------------------------------------------------------------------
  !
  !! ... isgn = +-1 : parallel 3d fft for rho and for the potential
  !                  NOT IMPLEMENTED WITH TASK GROUPS
  !! ... isgn = +-2 : parallel 3d fft for wavefunctions
  !
  !! ... isgn = +   : G-space to R-space, output = \sum_G f(G)exp(+iG*R)
  !! ...              fft along z using pencils        (env_cft_1z)
  !! ...              transpose across nodes           (env_fft_scatter)
  !! ...                 and reorder
  ! ...              fft along y (using planes) and x (env_cft_2xy)
  ! ... isgn = -   : R-space to G-space, output = \int_R f(R)exp(-iG*R)/Omega
  ! ...              fft along x and y(using planes)  (env_cft_2xy)
  ! ...              transpose across nodes           (env_fft_scatter)
  ! ...                 and reorder
  ! ...              fft along z using pencils        (env_cft_1z)
  !
  ! ...  The array "planes" signals whether a fft is needed along y :
  ! ...    planes(i)=0 : column f(i,*,*) empty , don't do fft along y
  ! ...    planes(i)=1 : column f(i,*,*) filled, fft along y needed
  ! ...  "empty" = no active components are present in f(i,*,*)
  ! ...            after (isgn>0) or before (isgn<0) the fft on z direction
  !
  ! ...  Note that if isgn=+/-1 (fft on rho and pot.) all fft's are needed
  ! ...  and all planes(i) are set to 1
  !
  ! This driver is based on code written by Stefano de Gironcoli for PWSCF.
  ! Task Group added by Costas Bekas, Oct. 2005, adapted from the CPMD code
  ! (Alessandro Curioni) and revised by Carlo Cavazzoni 2007.
  !
  USE env_fft_scalar, ONLY : env_cft_1z, env_cft_2xy
  USE env_fft_scatter_2d,   ONLY : env_fft_scatter
  USE env_fft_types,  ONLY : env_fft_type_descriptor
  !
  IMPLICIT NONE
  !
  COMPLEX(DP), INTENT(inout)    :: f( : )  ! array containing data to be transformed
  TYPE (env_fft_type_descriptor), INTENT(in) :: dfft
                                           ! descriptor of fft data layout
  INTEGER, INTENT(in)           :: isgn    ! fft direction
  !
  !
  INTEGER                    :: me_p
  INTEGER                    :: n1, n2, n3, nx1, nx2, nx3
  COMPLEX(DP), ALLOCATABLE   :: aux (:)
  INTEGER                    :: planes( dfft%nr1x )
  !LOGICAL                    :: use_tg
  !
  !
  IF (dfft%has_task_groups) CALL env_fftx_error__( ' env_tg_cft3s ', ' task groups on large mesh not implemented ', 1 )
  !
  n1  = dfft%nr1
  n2  = dfft%nr2
  n3  = dfft%nr3
  nx1 = dfft%nr1x
  nx2 = dfft%nr2x
  nx3 = dfft%nr3x
  !
  ALLOCATE( aux( dfft%nnr ) )
  !
  me_p = dfft%mype + 1
  !
  IF ( isgn > 0 ) THEN
     !
     IF ( isgn /= 2 ) THEN
        !
        CALL env_cft_1z( f, dfft%nsp( me_p ), n3, nx3, isgn, aux )
        !
        planes = dfft%iplp
        !
     ELSE
        !
        CALL env_cft_1z( f, dfft%nsw( me_p ), n3, nx3, isgn, aux )
        !
        planes = dfft%iplw
        !
     ENDIF
     !
     CALL env_fw_scatter( isgn ) ! forward scatter from stick to planes
     !
     CALL env_cft_2xy( f, dfft%my_nr3p, n1, n2, nx1, nx2, isgn, planes )
     !
  ELSE
     !
     IF ( isgn == -1 ) THEN
        !
        planes = dfft%iplp
        !
     ELSE IF ( isgn == -2 ) THEN
        !
        planes = dfft%iplw
        !
     ENDIF
     !
     CALL env_cft_2xy( f, dfft%my_nr3p, n1, n2, nx1, nx2, isgn, planes )
     !
     CALL env_bw_scatter( isgn )
     !
     IF ( isgn /= -2 ) THEN
        !
        CALL env_cft_1z( aux, dfft%nsp( me_p ), n3, nx3, isgn, f )
        !
     ELSE
        !
        CALL env_cft_1z( aux, dfft%nsw( me_p ), n3, nx3, isgn, f )
        !
     ENDIF
     !
  ENDIF
  !
  DEALLOCATE( aux )
  !
  RETURN
  !
CONTAINS
  !
  SUBROUTINE env_fw_scatter( iopt )

     !Transpose data for the 2-D FFT on the x-y plane
     !
     !NOGRP*dfft%nnr: The length of aux and f
     !nr3x: The length of each Z-stick
     !aux: input - output
     !f: working space
     !isgn: type of scatter
     !dfft%nsw(me) holds the number of Z-sticks proc. me has.
     !dfft%nr3p: number of planes per processor
     !
     !
     USE env_fft_scatter_2d, ONLY: env_fft_scatter
     !
     INTEGER, INTENT(in) :: iopt
     !
     IF( iopt == 2 ) THEN
        !
        CALL env_fft_scatter( dfft, aux, nx3, dfft%nnr, f, dfft%nsw, dfft%nr3p, iopt )
        !
     ELSEIF( iopt == 1 ) THEN
        !
        CALL env_fft_scatter( dfft, aux, nx3, dfft%nnr, f, dfft%nsp, dfft%nr3p, iopt )
        !
     ENDIF
     !
     RETURN
  END SUBROUTINE env_fw_scatter

  !

  SUBROUTINE env_bw_scatter( iopt )
     !
     USE env_fft_scatter_2d, ONLY: env_fft_scatter
     !
     INTEGER, INTENT(in) :: iopt
     !
     IF( iopt == -2 ) THEN
        !
        CALL env_fft_scatter( dfft, aux, nx3, dfft%nnr, f, dfft%nsw, dfft%nr3p, iopt )
        !
     ELSEIF( iopt == -1 ) THEN
        !
        CALL env_fft_scatter( dfft, aux, nx3, dfft%nnr, f, dfft%nsp, dfft%nr3p, iopt )
        !
     ENDIF
     !
     RETURN
  END SUBROUTINE env_bw_scatter
  !
END SUBROUTINE env_tg_cft3s
!
!
!
#if defined(__CUDA)
!----------------------------------------------------------------------------
SUBROUTINE env_tg_cft3s_gpu( f_d, dfft, isgn )
  !----------------------------------------------------------------------------
  !
  !! ... isgn = +-1 : parallel 3d fft for rho and for the potential
  !                  NOT IMPLEMENTED WITH TASK GROUPS
  !! ... isgn = +-2 : parallel 3d fft for wavefunctions
  !
  !! ... isgn = +   : G-space to R-space, output = \sum_G f(G)exp(+iG*R)
  !! ...              fft along z using pencils        (env_cft_1z)
  !! ...              transpose across nodes           (env_fft_scatter)
  !! ...                 and reorder
  ! ...              fft along y (using planes) and x (env_cft_2xy)
  ! ... isgn = -   : R-space to G-space, output = \int_R f(R)exp(-iG*R)/Omega
  ! ...              fft along x and y(using planes)  (env_cft_2xy)
  ! ...              transpose across nodes           (env_fft_scatter)
  ! ...                 and reorder
  ! ...              fft along z using pencils        (env_cft_1z)
  !
  ! ...  The array "planes" signals whether a fft is needed along y :
  ! ...    planes(i)=0 : column f(i,*,*) empty , don't do fft along y
  ! ...    planes(i)=1 : column f(i,*,*) filled, fft along y needed
  ! ...  "empty" = no active components are present in f(i,*,*)
  ! ...            after (isgn>0) or before (isgn<0) the fft on z direction
  !
  ! ...  Note that if isgn=+/-1 (fft on rho and pot.) all fft's are needed
  ! ...  and all planes(i) are set to 1
  !
  ! This driver is based on code written by Stefano de Gironcoli for PWSCF.
  ! Task Group added by Costas Bekas, Oct. 2005, adapted from the CPMD code
  ! (Alessandro Curioni) and revised by Carlo Cavazzoni 2007.
  !
  USE env_fft_scalar, ONLY : env_cft_1z_gpu, env_cft_2xy_gpu
  USE env_fft_scatter_2d_gpu,   ONLY : env_fft_scatter_gpu
  USE env_fft_types,  ONLY : env_fft_type_descriptor
  USE env_fft_buffers, ONLY : env_check_buffers_size, &
                            f_h => pin_space_scatter_in, &
                            aux_h => pin_space_scatter_out, &
                            aux_d => dev_space_fftparallel
  !
  IMPLICIT NONE
  !
  TYPE (env_fft_type_descriptor), INTENT(in) :: dfft
  COMPLEX(DP), DEVICE, INTENT(inout)    :: f_d( dfft%nnr ) ! array containing data to be transformed
                                           ! descriptor of fft data layout
  INTEGER, INTENT(in)           :: isgn    ! fft direction
  !
  !
  INTEGER                    :: me_p, istat
  INTEGER                    :: n1, n2, n3, nx1, nx2, nx3
  COMPLEX(DP), ALLOCATABLE   :: yf(:)
  INTEGER                    :: planes( dfft%nr1x )
  INTEGER(kind = cuda_stream_kind) :: stream  = 0
  !
  !
  n1  = dfft%nr1
  n2  = dfft%nr2
  n3  = dfft%nr3
  nx1 = dfft%nr1x
  nx2 = dfft%nr2x
  nx3 = dfft%nr3x
  !
  IF( dfft%has_task_groups ) CALL env_fftx_error__( ' env_tg_cft3s ', ' task groups in 2D + 1D decomposition not implemented ', 1 )
  !
  CALL env_check_buffers_size(dfft)
  !
  me_p = dfft%mype + 1
  !
  IF ( isgn > 0 ) THEN
     !
     IF ( isgn /= 2 ) THEN
        !
        CALL env_cft_1z_gpu( f_d, dfft%nsp( me_p ), n3, nx3, isgn, aux_d, stream )
        !
        planes = dfft%iplp
        !
     ELSE
        !
        CALL env_cft_1z_gpu( f_d, dfft%nsw( me_p ), n3, nx3, isgn, aux_d, stream )
        !
        planes = dfft%iplw
        !
     ENDIF
     !
     CALL env_fw_scatter_gpu( isgn ) ! forward scatter from stick to planes
     !
     CALL env_cft_2xy_gpu( f_d, aux_d, dfft%my_nr3p, n1, n2, nx1, nx2, isgn, stream, planes )
     !
  ELSE
     !
     IF ( isgn /= -2 ) THEN
        !
        planes = dfft%iplp
        !
     ELSE
        !
        planes = dfft%iplw
        !
     ENDIF

     CALL env_cft_2xy_gpu( f_d, aux_d, dfft%my_nr3p, n1, n2, nx1, nx2, isgn, stream, planes)
     !
     CALL env_bw_scatter_gpu( isgn )
     !
     !f_d = (0.d0, 0.d0)
     !
     IF ( isgn /= -2 ) THEN
        !
        CALL env_cft_1z_gpu( aux_d, dfft%nsp( me_p ), n3, nx3, isgn, f_d, stream )
        !
     ELSE
        !
        CALL env_cft_1z_gpu( aux_d, dfft%nsw( me_p ), n3, nx3, isgn, f_d, stream )
        !
     ENDIF
     !
  ENDIF
  !
  RETURN
  !
CONTAINS
  !
  SUBROUTINE env_fw_scatter_gpu( iopt )

     !Transpose data for the 2-D FFT on the x-y plane
     !
     !NOGRP*dfft%nnr: The length of aux and f
     !nr3x: The length of each Z-stick
     !aux: input - output
     !f: working space
     !isgn: type of scatter
     !dfft%nsw(me) holds the number of Z-sticks proc. me has.
     !dfft%nr3p: number of planes per processor
     !
     !
     USE env_fft_scatter_2d_gpu, ONLY: env_fft_scatter_gpu
     !
     INTEGER, INTENT(in) :: iopt
     !
     IF( iopt == 2 ) THEN
        !
        CALL env_fft_scatter_gpu( dfft, aux_d, aux_h, nx3, dfft%nnr, f_d, f_h, dfft%nsw, dfft%nr3p, iopt )
        !
     ELSEIF( iopt == 1 ) THEN
        !
        CALL env_fft_scatter_gpu( dfft, aux_d, aux_h, nx3, dfft%nnr, f_d, f_h, dfft%nsp, dfft%nr3p, iopt )
        !
     ENDIF
     !
     RETURN
  END SUBROUTINE env_fw_scatter_gpu

  !

  SUBROUTINE env_bw_scatter_gpu( iopt )
     !
     USE env_fft_scatter_2d_gpu, ONLY: env_fft_scatter_gpu
     !
     INTEGER, INTENT(in) :: iopt
     !
     IF( iopt == -2 ) THEN
        !
        CALL env_fft_scatter_gpu( dfft, aux_d, aux_h, nx3, dfft%nnr, f_d, f_h, dfft%nsw, dfft%nr3p, iopt )
        !
     ELSEIF( iopt == -1 ) THEN
        !
        CALL env_fft_scatter_gpu( dfft, aux_d, aux_h, nx3, dfft%nnr, f_d, f_h, dfft%nsp, dfft%nr3p, iopt )
        !
     ENDIF
     !
     RETURN
  END SUBROUTINE env_bw_scatter_gpu
  !
END SUBROUTINE env_tg_cft3s_gpu

SUBROUTINE env_many_cft3s_gpu( f_d, dfft, isgn, batchsize )
  !----------------------------------------------------------------------------
  !
  !
  !! ... isgn = +-1 : parallel 3d fft for rho and for the potential
  !                  NOT IMPLEMENTED WITH TASK GROUPS
  !! ... isgn = +-2 : parallel 3d fft for wavefunctions
  !
  !! ... isgn = +   : G-space to R-space, output = \sum_G f(G)exp(+iG*R)
  !! ...              fft along z using pencils        (env_cft_1z)
  !! ...              transpose across nodes           (env_fft_scatter)
  !! ...                 and reorder
  ! ...              fft along y (using planes) and x (env_cft_2xy)
  ! ... isgn = -   : R-space to G-space, output = \int_R f(R)exp(-iG*R)/Omega
  ! ...              fft along x and y(using planes)  (env_cft_2xy)
  ! ...              transpose across nodes           (env_fft_scatter)
  ! ...                 and reorder
  ! ...              fft along z using pencils        (env_cft_1z)
  !
  ! ...  The array "planes" signals whether a fft is needed along y :
  ! ...    planes(i)=0 : column f(i,*,*) empty , don't do fft along y
  ! ...    planes(i)=1 : column f(i,*,*) filled, fft along y needed
  ! ...  "empty" = no active components are present in f(i,*,*)
  ! ...            after (isgn>0) or before (isgn<0) the fft on z direction
  !
  ! ...  Note that if isgn=+/-1 (fft on rho and pot.) all fft's are needed
  ! ...  and all planes(i) are set to 1
  !
  ! ...  batchsize : number of 3D FFTs contained in f_d to be transformed.
  ! ...              Must be 1 < batchsize <= dfft%batchsize.
  !
  ! This driver is based on code written by Stefano de Gironcoli for PWSCF.
  ! Task Group added by Costas Bekas, Oct. 2005, adapted from the CPMD code
  ! (Alessandro Curioni) and revised by Carlo Cavazzoni 2007.
  !
  ! The GPU version is based on code written by Josh Romero, Everett Phillips
  ! and Massimiliano Fatica and revised by Pietro Bonfà.
  !
  ! The current version performs batchsize FFTs and overlaps computation
  ! with MPI communications and data transfers between host and device.
  !
  USE env_fft_scalar, ONLY : env_cft_1z_gpu, env_cft_2xy_gpu
  USE env_fft_scatter_2d_gpu,   ONLY : env_fft_scatter_many_columns_to_planes_send, &
                                   env_fft_scatter_many_columns_to_planes_store, &
                                   env_fft_scatter_many_planes_to_columns_send, &
                                   env_fft_scatter_many_planes_to_columns_store
  USE env_fft_types,  ONLY : env_fft_type_descriptor
  USE env_fft_buffers, ONLY : env_check_buffers_size, &
                            f_h => pin_space_scatter_in, &
                            aux_h => pin_space_scatter_out, &
                            aux_d => dev_space_fftparallel, &
                            aux2_h => pin_space_scatter_dblbuffer, &
                            aux2_d => dev_space_scatter_dblbuffer
  !
  IMPLICIT NONE
  !
  TYPE (env_fft_type_descriptor), INTENT(in) :: dfft
                                           ! descriptor of fft data layout
  INTEGER, INTENT(in)           :: isgn    ! fft direction
  INTEGER, INTENT(in)           :: batchsize
  COMPLEX(DP), DEVICE, INTENT(inout)    :: f_d( batchsize * dfft%nnr ) ! array containing data to be transformed
  !
  INTEGER                    :: me_p, istat, i, j, currsize
  INTEGER                    :: n1, n2, n3, nx1, nx2, nx3, ncpx, nppx, proc
  COMPLEX(DP), ALLOCATABLE   :: yf(:)
  INTEGER                    :: planes( dfft%nr1x )
  INTEGER                    :: sticks( dfft%nproc  )
  INTEGER(kind = cuda_stream_kind) :: stream  = 0
  !
  !
  n1  = dfft%nr1
  n2  = dfft%nr2
  n3  = dfft%nr3
  nx1 = dfft%nr1x
  nx2 = dfft%nr2x
  nx3 = dfft%nr3x
  !
  CALL env_check_buffers_size(dfft, batchsize)
  !
  me_p = dfft%mype + 1
  !
  ncpx = 0
  nppx = 0
  DO proc = 1, dfft%nproc
     IF ( abs(isgn) == 2 ) ncpx = max( ncpx, dfft%nsw ( proc ) )
     IF ( abs(isgn) == 1 ) ncpx = max( ncpx, dfft%nsp ( proc ) )
     nppx = max( nppx, dfft%nr3p ( proc ) )
  ENDDO
  IF ( abs(isgn) == 2 ) sticks = dfft%nsw
  IF ( abs(isgn) == 1 ) sticks = dfft%nsp
  !
  IF ( (abs(isgn) /= 2) .and. (abs(isgn) /= 1) ) &
     CALL env_fftx_error__( ' env_many_cft3s_gpu ', ' abs(isgn) /= 1 or 2 not implemented ', isgn )
  !
  IF (dfft%nproc <= 1) CALL env_fftx_error__( ' env_many_cft3s_gpu ', ' this subroutine should never be called with nproc= ', dfft%nproc )
  !
  ! FFTs are done in sub-batches of dfft%subbatchsize (default is 4)
  ! When a sub-batch has been transformed in a direction or a plane,
  ! communication between device and host is started and the next subbatch is transformed.
  ! Later, the subbatch is received on target MPI process and transformed
  ! overlapping computation with MPI communication.
  !
  IF ( isgn > 0 ) THEN
     DO j = 0, batchsize-1, dfft%subbatchsize
       ! determine whether the FFTs that are left are less than the maximum
       ! subbatchsize size.
       currsize = min(dfft%subbatchsize, batchsize - j)
       !
       IF ( isgn /= 2 ) THEN
          !
          planes = dfft%iplp
          !
       ELSE
          !
          planes = dfft%iplw
          !
       ENDIF
       !
       ! perform the FFT along one direction and, at the same time,
       ! read data spaced by dfft%nnr and store in in the output
       ! with spacing ncpx*nx3, making it easy to bach communication.
       DO i = 0, currsize - 1
         CALL env_cft_1z_gpu( f_d((j+i)*dfft%nnr + 1:), sticks(me_p), n3, nx3, isgn, aux_d(j*dfft%nnr + i*ncpx*nx3 +1:), dfft%a2a_comp )
       ENDDO
       !
       i = cudaEventRecord(dfft%bevents(j/dfft%subbatchsize + 1), dfft%a2a_comp)
       i = cudaStreamWaitEvent(dfft%bstreams(j/dfft%subbatchsize + 1), dfft%bevents(j/dfft%subbatchsize + 1), 0)

       IF (j > 0) i = cudaStreamWaitEvent(dfft%bstreams(j/dfft%subbatchsize + 1), dfft%bevents(j/dfft%subbatchsize), 0)

       CALL env_fft_scatter_many_columns_to_planes_store( dfft, aux_d(j*dfft%nnr + 1:), aux_h(j*dfft%nnr + 1:), nx3, dfft%nnr, f_d(j*dfft%nnr + 1:), &
         f_h(j*dfft%nnr + 1:), aux2_d(j*dfft%nnr + 1:), aux2_h(j*dfft%nnr + 1:), sticks, dfft%nr3p, isgn, currsize, j/dfft%subbatchsize + 1 )

     ENDDO

     DO j = 0, batchsize-1, dfft%subbatchsize
       currsize = min(dfft%subbatchsize, batchsize - j)

       CALL env_fft_scatter_many_columns_to_planes_send( dfft, aux_d(j*dfft%nnr + 1:), aux_h(j*dfft%nnr + 1:), nx3, dfft%nnr, f_d(j*dfft%nnr + 1:), &
         f_h(j*dfft%nnr + 1:), aux2_d(j*dfft%nnr + 1:), aux2_h(j*dfft%nnr + 1:), sticks, dfft%nr3p, isgn, currsize, j/dfft%subbatchsize + 1 )

       IF (currsize == dfft%subbatchsize) THEN
         CALL env_cft_2xy_gpu( f_d(j*dfft%nnr + 1:), aux_d(j*dfft%nnr + 1:), currsize * nppx, n1, n2, nx1, nx2, isgn, dfft%a2a_comp, planes )
       ELSE
         DO i = 0, currsize - 1
           CALL env_cft_2xy_gpu( f_d((j+i)*dfft%nnr + 1:), aux_d((j+i)*dfft%nnr + 1:), dfft%nr3p( me_p ), n1, n2, nx1, nx2, isgn,  &
           dfft%a2a_comp, planes )
         ENDDO
       ENDIF

     ENDDO

!     i = cudaDeviceSynchronize()

     !
  ELSE
!     i = cudaDeviceSynchronize()

     DO j = 0, batchsize-1, dfft%subbatchsize
       currsize = min(dfft%subbatchsize, batchsize - j)
       !
       IF ( isgn /= -2 ) THEN
          !
          planes = dfft%iplp
          !
       ELSE
          !
          planes = dfft%iplw
          !
       ENDIF

       IF (currsize == dfft%subbatchsize) THEN
         CALL env_cft_2xy_gpu( f_d(j*dfft%nnr + 1:), aux_d(j*dfft%nnr + 1:), currsize * nppx, n1, n2, nx1, nx2, isgn, dfft%a2a_comp, planes )
       ELSE
         DO i = 0, currsize - 1
           CALL env_cft_2xy_gpu( f_d((j+i)*dfft%nnr + 1:), aux_d((j+i)*dfft%nnr + 1:), dfft%nr3p( me_p ), n1, n2, nx1, nx2, isgn, dfft%a2a_comp, planes )
         ENDDO
       ENDIF

       IF (j > 0) i = cudaStreamWaitEvent(dfft%bstreams(j/dfft%subbatchsize + 1), dfft%bevents(j/dfft%subbatchsize), 0)

       CALL env_fft_scatter_many_planes_to_columns_store( dfft, aux_d(j*dfft%nnr + 1:), aux_h(j*dfft%nnr + 1:), nx3, dfft%nnr, f_d(j*dfft%nnr + 1:), &
         f_h(j*dfft%nnr + 1:), aux2_d(j*dfft%nnr + 1:), aux2_h(j*dfft%nnr + 1:), sticks, dfft%nr3p, isgn, currsize, j/dfft%subbatchsize + 1 )

     ENDDO

     DO j = 0, batchsize-1, dfft%subbatchsize
       currsize = min(dfft%subbatchsize, batchsize - j)

       CALL env_fft_scatter_many_planes_to_columns_send( dfft, aux_d(j*dfft%nnr + 1:), aux_h(j*dfft%nnr + 1:), nx3, dfft%nnr, f_d(j*dfft%nnr + 1:), &
         f_h(j*dfft%nnr + 1:), aux2_d(j*dfft%nnr + 1:), aux2_h(j*dfft%nnr + 1:), sticks, dfft%nr3p, isgn, currsize, j/dfft%subbatchsize + 1 )


       i = cudaEventRecord(dfft%bevents(j/dfft%subbatchsize + 1), dfft%bstreams(j/dfft%subbatchsize + 1))
       i = cudaStreamWaitEvent(dfft%a2a_comp, dfft%bevents(j/dfft%subbatchsize + 1), 0)

       DO i = 0, currsize - 1
         CALL env_cft_1z_gpu( aux_d(j*dfft%nnr + i*ncpx*nx3 + 1:), sticks( me_p ), n3, nx3, isgn, f_d((j+i)*dfft%nnr + 1:), dfft%a2a_comp )
       ENDDO

     ENDDO
!    i = cudaDeviceSynchronize()
  ENDIF
  !
  RETURN
  !
END SUBROUTINE env_many_cft3s_gpu
#endif
!
END MODULE env_fft_parallel_2d
