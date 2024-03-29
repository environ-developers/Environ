MODULE env_fft_buffers
  !
  USE env_fft_param
  IMPLICIT NONE
  SAVE
  !
  INTEGER :: current_size = 0
  COMPLEX(DP), ALLOCATABLE :: dev_space_fftparallel(:)
  COMPLEX(DP), ALLOCATABLE :: dev_space_scatter_dblbuffer(:)
  COMPLEX(DP), ALLOCATABLE :: pin_space_scatter_dblbuffer(:)
  COMPLEX(DP), ALLOCATABLE :: pin_space_scatter_in(:)
  COMPLEX(DP), ALLOCATABLE :: pin_space_scatter_out(:)
#if defined(__CUDA)
  attributes(DEVICE) :: dev_space_fftparallel
  attributes(DEVICE) :: dev_space_scatter_dblbuffer
  attributes(PINNED) :: pin_space_scatter_in
  attributes(PINNED) :: pin_space_scatter_out
  attributes(PINNED) :: pin_space_scatter_dblbuffer
#endif
  !
  PUBLIC :: env_check_buffers_size, env_deallocate_buffers
  !
CONTAINS
  !
  SUBROUTINE env_check_buffers_size( desc, howmany )
    USE env_fft_types,      ONLY : env_fft_type_descriptor
    TYPE(env_fft_type_descriptor), INTENT(in) :: desc
    INTEGER, OPTIONAL, INTENT(IN) :: howmany
    INTEGER :: howmany_
    INTEGER :: info

    IF (PRESENT(howmany)) THEN
      howmany_ = howmany
    ELSE
      howmany_ = 1
    ENDIF
    !
    IF (current_size < desc%nnr * howmany_) THEN
      !
      current_size = desc%nnr * howmany_
      !
      IF( ALLOCATED( dev_space_fftparallel ) ) DEALLOCATE( dev_space_fftparallel )
      IF( ALLOCATED( pin_space_scatter_in  ) ) DEALLOCATE( pin_space_scatter_in  )
      IF( ALLOCATED( pin_space_scatter_out ) ) DEALLOCATE( pin_space_scatter_out )
      !
      ALLOCATE(dev_space_fftparallel(current_size), STAT=info)
      IF ( info /= 0 ) CALL env_fftx_error__( ' env_fft_buffers ', ' Allocation failed ', 1 )
      ALLOCATE(pin_space_scatter_in (current_size), STAT=info)
      IF ( info /= 0 ) CALL env_fftx_error__( ' env_fft_buffers ', ' Allocation failed ', 2 )
      ALLOCATE(pin_space_scatter_out(current_size), STAT=info)
      IF ( info /= 0 ) CALL env_fftx_error__( ' env_fft_buffers ', ' Allocation failed ', 3 )
      !
      ! Slab decomposition implements double buffering
      IF ( (.not. desc%use_pencil_decomposition ) ) THEN
        IF( ALLOCATED( dev_space_scatter_dblbuffer ) ) DEALLOCATE( dev_space_scatter_dblbuffer )
        IF( ALLOCATED( pin_space_scatter_dblbuffer ) ) DEALLOCATE( pin_space_scatter_dblbuffer  )
        !
        ALLOCATE(dev_space_scatter_dblbuffer(current_size), STAT=info)
        IF ( info /= 0 ) CALL env_fftx_error__( ' env_fft_buffers ', ' Allocation failed ', 4 )
        ALLOCATE(pin_space_scatter_dblbuffer(current_size), STAT=info)
        IF ( info /= 0 ) CALL env_fftx_error__( ' env_fft_buffers ', ' Allocation failed ', 5 )
      END IF
    END IF
    !
  END SUBROUTINE env_check_buffers_size
  !
  SUBROUTINE env_deallocate_buffers()
    current_size = 0
    IF( ALLOCATED( dev_space_fftparallel ) ) DEALLOCATE( dev_space_fftparallel )
    IF( ALLOCATED( pin_space_scatter_in  ) ) DEALLOCATE( pin_space_scatter_in  )
    IF( ALLOCATED( pin_space_scatter_out ) ) DEALLOCATE( pin_space_scatter_out )
    IF( ALLOCATED( dev_space_scatter_dblbuffer ) ) DEALLOCATE( dev_space_scatter_dblbuffer  )
    IF( ALLOCATED( pin_space_scatter_dblbuffer ) ) DEALLOCATE( pin_space_scatter_dblbuffer )
  END SUBROUTINE env_deallocate_buffers

END MODULE env_fft_buffers
