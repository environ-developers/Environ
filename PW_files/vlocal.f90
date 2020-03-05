MODULE env_vlocal
    !
    ! ... The variables needed for the local potential in reciprocal space
    !
    USE env_kinds, ONLY : DP
    USE env_parameters, ONLY : ntypx
    !
    SAVE
    !
    COMPLEX(DP), ALLOCATABLE :: &
         strf(:,:)              ! the structure factor
    REAL(DP), ALLOCATABLE :: &
         vloc(:,:)              ! the local potential for each atom type
    REAL(DP) :: &
         starting_charge(ntypx) ! the atomic charge used to start with
    !
  END MODULE env_vlocal

