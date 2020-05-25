MODULE core_base
  !
  USE core_types
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
       fft
  LOGICAL ::                           &
       loned_analytic
  TYPE( oned_analytic_core ) ::        &
       oned_analytic
  !
END MODULE core_base
