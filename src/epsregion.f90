!
! Copyright (C) 2007-2008 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! original version by O. Andreussi and N. Marzari
!
!--------------------------------------------------------------------
MODULE epsregion
!--------------------------------------------------------------------

  USE kinds,          ONLY: DP
  USE constants,      ONLY: sqrtpi, fpi, pi
  USE io_global,      ONLY: stdout
  USE environ_base,   ONLY: verbose, environ_unit,                   &
                            system_pos, system_width,                &
                            env_static_permittivity,                 &
                            env_optical_permittivity,                &
                            env_dielectric_regions,                  &
                            epsregion_origin, epsregion_eps,         &
                            epsregion_dim, epsregion_axis,           &
                            epsregion_pos, epsregion_width,          &
                            epsregion_spread, epsstatic, epsoptical
  USE environ_debug,  ONLY: write_cube
! BACKWARD COMPATIBILITY
! Compatible with QE-5.1.X and QE-5.2.0
!
! Compatible with QE-5.2.1, QE-5.3.0 and svn
  USE generate_function, ONLY: generate_erfc, erfcvolume
! END BACKWARD COMPATIBILITY
  !
  IMPLICIT NONE
  !
  SAVE
  !
  PRIVATE
  !
  PUBLIC :: generate_epsregion
  !
CONTAINS
  !
!--------------------------------------------------------------------
  SUBROUTINE generate_epsregion( nnr, alat, omega, at )
!--------------------------------------------------------------------
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: nnr
  REAL(DP), INTENT(IN) :: alat, omega
  REAL(DP), DIMENSION(3,3), INTENT(IN) :: at
  !
  ! ... Local variables
  !
  LOGICAL :: shift
  INTEGER :: idr, dim, axis
  REAL(DP) :: spread, width, charge, norm
  REAL(DP), DIMENSION(3) :: pos, pos0
  !
  REAL(DP), DIMENSION(:), ALLOCATABLE :: epslocal
  !
  ALLOCATE(epslocal(nnr))
  !
  ! Use the origin of the cell as default origin, removed other options
  !
  SELECT CASE ( epsregion_origin )
  CASE ( 0 )
     shift = .false.
     pos0 = 0.D0
  CASE ( 1 )
     shift = .false.
     pos0 = system_pos
  CASE ( 2 )
     shift = .true.
     pos0 = system_pos
  END SELECT
  !
  ! ... Generate the dielectric regions (static)
  !
  epsstatic = env_static_permittivity
  !
  DO idr = 1, env_dielectric_regions
     !
     norm = SQRT(SUM(epsregion_pos(:,idr)**2))
     IF ( .NOT. shift ) THEN
       ! No bounding box
       pos(:) = epsregion_pos(:,idr)/alat + pos0(:)
       width = epsregion_width(idr)
     ELSE IF ( norm .NE. 0.D0 ) THEN
       ! The region is not centered on the system, position is defined wrt to bounding box
       ! but do not change the width of the region
       pos(:) = epsregion_pos(:,idr)/alat*(1.D0+system_width/norm) + pos0(:)
       width = epsregion_width(idr)
     ELSE
       ! The region is centered on the system, position is not changed
       ! but the width is defined in addition to the system width
       pos(:) = epsregion_pos(:,idr)/alat + pos0(:)
       width = epsregion_width(idr) + system_width
     ENDIF
     !
     dim = epsregion_dim(idr)
     axis = epsregion_axis(idr)
     spread = epsregion_spread(idr)
     epslocal = 0.D0
! BACKWARD COMPATIBILITY
! Compatible with QE-5.1.X QE-5.2.0
!     WRITE(stdout,2000)
!     STOP
! Compatible with QE-5.2.1, QE-5.3.0 and svn
     charge = erfcvolume(dim,axis,width,spread,alat,omega,at)
     CALL generate_erfc( nnr, dim, axis, charge, width, spread, pos, epslocal )
! END BACKWARD COMPATIBILITY
     epsstatic(:) = epsstatic(:) - ( epsstatic(:) - epsregion_eps(1,idr) ) * epslocal(:)
     !
  END DO
  !
  IF ( verbose .GE. 3 ) CALL write_cube( nnr, epsstatic, 'epsstatic.cube' )
  !
  ! ... Generate the dielectric regions (optical)
  !
  epsoptical = env_optical_permittivity
  !
  DO idr = 1, env_dielectric_regions
     !
     norm = SQRT(SUM(epsregion_pos(:,idr)**2))
     IF ( .NOT. shift ) THEN
       ! No bounding box
       pos(:) = epsregion_pos(:,idr)/alat + pos0(:)
       width = epsregion_width(idr)
     ELSE IF ( norm .NE. 0.D0 ) THEN
       ! The region is not centered on the system, position is defined wrt to bounding box
       ! but do not change the width of the region
       pos(:) = epsregion_pos(:,idr)/alat*(1.D0+system_width/norm) + pos0(:)
       width = epsregion_width(idr)
     ELSE
       ! The region is centered on the system, position is not changed
       ! but the width is defined in addition to the system width
       pos(:) = epsregion_pos(:,idr)/alat + pos0(:)
       width = epsregion_width(idr) + system_width
     ENDIF
     !
     dim = epsregion_dim(idr)
     axis = epsregion_axis(idr)
     spread = epsregion_spread(idr)
     epslocal = 0.D0
! BACKWARD COMPATIBILITY
! Compatible with QE-5.1.X QE-5.2.0
!     WRITE(stdout,2000)
!     STOP
! Compatible with QE-5.2.1, QE-5.3.0 and svn
     charge = erfcvolume(dim,axis,width,spread,alat,omega,at)
     CALL generate_erfc( nnr, dim, axis, charge, width, spread, pos, epslocal )
! END BACKWARD COMPATIBILITY
     epsoptical(:) = epsoptical(:) - ( epsoptical(:) - epsregion_eps(2,idr) ) * epslocal(:)
     !
  END DO
  !
  IF ( verbose .GE. 3 ) CALL write_cube( nnr, epsoptical, 'epsoptical.cube' )
  !
  DEALLOCATE(epslocal)
  !
  RETURN
  !
2000 FORMAT('ERROR: dielectric regions only available for QE-5.2.1 and later releases',i3)
  !
!--------------------------------------------------------------------
  END SUBROUTINE generate_epsregion
!--------------------------------------------------------------------
!
END MODULE epsregion
