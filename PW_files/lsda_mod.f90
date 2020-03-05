MODULE env_lsda_mod
    !
    ! ... The variables needed for the lsda calculation
    !
    USE env_kinds,      ONLY : DP
    USE env_parameters, ONLY : ntypx, npk
    !
    IMPLICIT NONE
    SAVE
    !
    LOGICAL :: &
         lsda
    REAL(DP) :: &
         magtot,                       &! total magnetization
         absmag,                       &! total absolute magnetization
         starting_magnetization(ntypx)  ! the magnetization used to start with
    INTEGER :: &
         nspin,           &! number of spin polarization: 2 if lsda, 1 other
         current_spin,    &! spin of the current kpoint
         isk(npk)          ! for each k-point: 1=spin up, 2=spin down
    !
  END MODULE env_lsda_mod