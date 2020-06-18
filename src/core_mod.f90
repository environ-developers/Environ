MODULE core_base
  !
  USE core_types
  USE fft_types
  !
  SAVE
  !
  ! Internal setup of numerical cores
  !
  LOGICAL ::                           &
       lfd
  TYPE( fd_core ) ::                   &
       fd
  LOGICAL ::                           &
       lfft
  TYPE( fft_core ) ::                  &
       sys_fft
  TYPE( fft_core ) ::          &
       environment_fft
  LOGICAL ::                           &
       loned_analytic
  TYPE( oned_analytic_core ) ::        &
       oned_analytic
  TYPE( fft_type_descriptor ) ::       &
       sys_dfft
  TYPE( fft_type_descriptor ) ::       &
       environment_dfft
  !
END MODULE core_base
