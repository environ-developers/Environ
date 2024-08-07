!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
subroutine env_fft_interpolate_real (dfft_in, v_in, dfft_out, v_out )
  !
  !   This subroutine interpolates an array  v_in   defined on fft grid  dfft_in
  !                           to   an array  v_out  defined on fft grid  dfft_out
  !   v_in and v_out are assumed to be real arrays and may concide
  !
  USE env_fft_param,      ONLY : DP
  USE env_fft_types,      ONLY : env_fft_type_descriptor
  USE env_fft_interfaces, ONLY : env_fwfft, env_invfft
  ! I/O variables
  TYPE(env_fft_type_descriptor), INTENT(IN) :: dfft_in, dfft_out
  REAL(DP),INTENT(IN)  :: v_in (:) !dfft_in%nnr)
  REAL(DP),INTENT(OUT) :: v_out (:) !dfft_out%nnr)
  ! local variables
  INTEGER :: ngm
  COMPLEX(DP), ALLOCATABLE :: aux_in (:), aux_out (:)

  call env_start_clock ('interpolate')

  IF (dfft_out%grid_id == dfft_in%grid_id) THEN

     v_out (1:dfft_in%nnr) = v_in (1:dfft_in%nnr)

  ELSE

     if (dfft_in%lgamma .neqv. dfft_out%lgamma) &
        call env_fftx_error__ ('env_fft_interpolate_real','two grids with inconsistent lgamma values', 1)

     ALLOCATE (aux_in( dfft_in%nnr), aux_out(dfft_out%nnr))

     aux_in (1:dfft_in%nnr) = v_in(1:dfft_in%nnr)

     CALL env_fwfft ('Rho', aux_in, dfft_in)

     aux_out(1:dfft_out%nnr) = (0.d0, 0.d0)

     ngm = min(dfft_in%ngm, dfft_out%ngm)

     aux_out (dfft_out%nl (1:ngm) ) = aux_in (dfft_in%nl (1:ngm) )
     IF (dfft_in%lgamma) aux_out (dfft_out%nlm (1:ngm) ) = aux_in (dfft_in%nlm (1:ngm) )

     CALL env_invfft ('Rho', aux_out, dfft_out)

     v_out (1:dfft_out%nnr) = aux_out (1:dfft_out%nnr)

     DEALLOCATE (aux_in, aux_out)

  END IF

  call env_stop_clock ('interpolate')

  return

end subroutine env_fft_interpolate_real
!
subroutine env_fft_interpolate_complex (dfft_in, v_in, dfft_out, v_out )
  !
  !   This subroutine interpolates an array  v_in   defined on fft grid  dfft_in
  !                           to   an array  v_out  defined on fft grid  dfft_out
  !   v_in and v_out are assumed to be complex arrays and may concide
  !
  USE env_fft_param,      ONLY : DP
  USE env_fft_types,      ONLY : env_fft_type_descriptor
  USE env_fft_interfaces, ONLY : env_fwfft, env_invfft
  ! I/O variables
  TYPE(env_fft_type_descriptor), INTENT(IN) :: dfft_in, dfft_out
  COMPLEX(DP),INTENT(IN)  :: v_in (:) !dfft_in%nnr)
  COMPLEX(DP),INTENT(OUT) :: v_out (:) !dfft_out%nnr)
  ! local variables
  INTEGER :: ngm
  COMPLEX(DP), ALLOCATABLE :: aux_in (:)

  if (dfft_out%lgamma.OR.dfft_in%lgamma) call env_fftx_error__('env_fft_interpolate_complex','lgamma not allowed', 1)

  call env_start_clock ('interpolate')

  IF (dfft_out%grid_id == dfft_in%grid_id) THEN

     v_out (1:dfft_in%nnr) = v_in (1:dfft_in%nnr)

  ELSE

     ALLOCATE (aux_in( dfft_in%nnr))

     aux_in (1:dfft_in%nnr) = v_in(1:dfft_in%nnr)

     CALL env_fwfft ('Rho', aux_in, dfft_in)

     v_out(1:dfft_out%nnr) = (0.d0, 0.d0)

     ngm = min(dfft_in%ngm, dfft_out%ngm)

     v_out (dfft_out%nl (1:ngm) ) = aux_in (dfft_in%nl (1:ngm) )

     CALL env_invfft ('Rho', v_out, dfft_out)

     DEALLOCATE (aux_in)

  END IF

  call env_stop_clock ('interpolate')

  return

end subroutine env_fft_interpolate_complex
!
