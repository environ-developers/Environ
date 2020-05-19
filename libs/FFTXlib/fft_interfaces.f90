!
! Copyright (C) Quantum ESPRESSO group
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!=---------------------------------------------------------------------------=!
MODULE env_fft_interfaces

  IMPLICIT NONE
  PRIVATE


  PUBLIC :: env_fwfft, env_invfft, env_fft_interpolate

  
  INTERFACE env_invfft
     !! invfft is the interface to both the standard fft **invfft_x**,
     !! and to the "box-grid" version **invfft_b**, used only in CP 
     !! (the latter has an additional argument)
     
     SUBROUTINE env_invfft_y( fft_kind, f, dfft, howmany )
       USE env_fft_types,  ONLY: fft_type_descriptor
       USE env_fft_param,  ONLY :DP
       IMPLICIT NONE
       CHARACTER(LEN=*),  INTENT(IN) :: fft_kind
       TYPE(fft_type_descriptor), INTENT(IN) :: dfft
       INTEGER, OPTIONAL, INTENT(IN) :: howmany
       COMPLEX(DP) :: f(:)
     END SUBROUTINE env_invfft_y
     !
     SUBROUTINE env_invfft_b( f, dfft, ia )
       USE env_fft_smallbox_type,  ONLY: fft_box_descriptor
       USE env_fft_param,  ONLY :DP
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: ia
       TYPE(fft_box_descriptor), INTENT(IN) :: dfft
       COMPLEX(DP) :: f(:)
     END SUBROUTINE env_invfft_b
  END INTERFACE

  INTERFACE env_fwfft
     SUBROUTINE env_fwfft_y( fft_kind, f, dfft, howmany )
       USE env_fft_types,  ONLY: fft_type_descriptor
       USE env_fft_param,  ONLY :DP
       IMPLICIT NONE
       CHARACTER(LEN=*), INTENT(IN) :: fft_kind
       TYPE(fft_type_descriptor), INTENT(IN) :: dfft
       INTEGER, OPTIONAL, INTENT(IN) :: howmany
       COMPLEX(DP) :: f(:)
     END SUBROUTINE env_fwfft_y
  END INTERFACE

  INTERFACE env_fft_interpolate
     !! fft_interpolate  is the interface to utility that fourier interpolate
     !! real/complex arrays between two grids 
     
     SUBROUTINE env_fft_interpolate_real( dfft_in, v_in, dfft_out, v_out )
       USE env_fft_param,  ONLY :DP
       USE env_fft_types,  ONLY: fft_type_descriptor
       IMPLICIT NONE
       TYPE(fft_type_descriptor), INTENT(IN) :: dfft_in, dfft_out
       REAL(DP), INTENT(IN)  :: v_in(:)
       REAL(DP), INTENT(OUT) :: v_out(:)
     END SUBROUTINE env_fft_interpolate_real
     !
     SUBROUTINE env_fft_interpolate_complex( dfft_in, v_in, dfft_out, v_out )
       USE env_fft_param,  ONLY :DP
       USE env_fft_types,  ONLY: fft_type_descriptor
       IMPLICIT NONE
       TYPE(fft_type_descriptor), INTENT(IN) :: dfft_in, dfft_out
       COMPLEX(DP), INTENT(IN)  :: v_in(:)
       COMPLEX(DP), INTENT(OUT) :: v_out(:)
     END SUBROUTINE env_fft_interpolate_complex
  END INTERFACE

END MODULE env_fft_interfaces
!=---------------------------------------------------------------------------=!
