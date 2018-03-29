! Copyright (C) 2018 ENVIRON (www.quantum-environment.org)
!
!    This file is part of Environ version 1.0
!
!    Environ 1.0 is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 2 of the License, or
!    (at your option) any later version.
!
!    Environ 1.0 is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more detail, either the file
!    `License' in the root directory of the present distribution, or
!    online at <http://www.gnu.org/licenses/>.
!
! Module to include a real-space correction of the potential in order to
! remove periodic boundary conditions on a partially periodic system.
! The correction is computed for the time being using the planar average
! approximation and is of quadratic nature: the first order proportional
! to the dipole of the system along the direction perpendicular to the slab,
! the second order proportional to the total charge of the system.
!
! Authors: Oliviero Andreussi (Department of Physics, UNT)
!          Nicola Marzari     (THEOS and NCCR-MARVEL, EPFL)
!
!----------------------------------------------------------------------------
MODULE correction_periodic
!----------------------------------------------------------------------------
  !
  ! ... The variables needed to correct periodic boundary conditions
  !     for a partially periodic system. Real space correction with planar
  !     average approximation.
  !
  USE environ_types
  USE electrostatic_types
  USE environ_output
  USE environ_base,      ONLY : e2
  !
  IMPLICIT NONE
  !
  REAL(DP), PARAMETER :: madelung(3) = (/ 2.837297479D0, 2.8883D0, 2.885D0 /)
  !
  PRIVATE
  !
  PUBLIC :: calc_v0periodic, calc_vperiodic, calc_fperiodic, calc_gradvperiodic
  !
CONTAINS
!---------------------------------------------------------------------------
  SUBROUTINE calc_v0periodic( oned_analytic, charges, v0 )
!---------------------------------------------------------------------------
    !
    ! ... Compute potential to correct for pbc along the direction
    !     orthogonal to the slab
    !
    IMPLICIT NONE
    !
    TYPE( oned_analytic_core ), TARGET, INTENT(IN) :: oned_analytic
    TYPE( environ_density ), TARGET, INTENT(IN) :: charges
    REAL( DP ), INTENT(INOUT) :: v0
    !
    INTEGER, POINTER :: nnr
    REAL( DP ), DIMENSION(:), POINTER :: rhotot
    !
    INTEGER, POINTER :: env_periodicity
    INTEGER, POINTER :: slab_axis
    REAL( DP ), POINTER :: alat, omega, axis_length
    REAL( DP ), DIMENSION(:), POINTER :: origin
    !
    INTEGER :: icor
    REAL(DP) :: fact, const
    REAL(DP) :: dipole(0:3), quadrupole(3)
    REAL(DP) :: tot_charge, tot_dipole(3), tot_quadrupole(3)
    CHARACTER( LEN = 80 ) :: sub_name = 'calc_v0periodic'
    !
    nnr => charges % cell % nnr
    rhotot => charges % of_r
    !
    alat => oned_analytic % alat
    omega => oned_analytic % omega
    env_periodicity => oned_analytic % d
    slab_axis => oned_analytic % axis
    axis_length => oned_analytic % size
    origin => oned_analytic % origin
    !
    ! ... Compute multipoles of the system wrt the chosen origin
    !
    CALL compute_dipole( nnr, 1, rhotot, origin, dipole, quadrupole )
    !
    tot_charge = dipole(0)
    tot_dipole = dipole(1:3)
    tot_quadrupole = quadrupole
    !
    ! ... Compute constant shift for PBC correction
    !
    fact = e2 * tpi / omega
    !
    SELECT CASE ( env_periodicity )
       !
    CASE ( 0 )
       !
       const = madelung(1) * tot_charge / alat * e2 - fact * SUM(tot_quadrupole(:)) / 3.D0
       !
    CASE ( 1 )
       !
       CALL errore(sub_name,'Option not yet implemented',1)
       !
    CASE ( 2 )
       !
       const = - pi / 3.D0 * tot_charge / axis_length * e2 - fact * tot_quadrupole(slab_axis)
       !
    CASE ( 3 )
       !
       const = 0.D0
       !
    CASE DEFAULT
       !
       CALL errore(sub_name,'Unexpected option in dimensionality of PBC correction',1)
       !
    END SELECT
    !
    v0 = const
    !
    RETURN
!---------------------------------------------------------------------------
  END SUBROUTINE calc_v0periodic
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
  SUBROUTINE calc_vperiodic( oned_analytic, charges, potential )
!---------------------------------------------------------------------------
    !
    ! ... Compute potential to correct for pbc along the direction
    !     orthogonal to the slab
    !
    IMPLICIT NONE
    !
    TYPE( oned_analytic_core ), TARGET, INTENT(IN) :: oned_analytic
    TYPE( environ_density ), TARGET, INTENT(IN) :: charges
    TYPE( environ_density ), INTENT(INOUT) :: potential
    !
    INTEGER, POINTER :: nnr
    REAL( DP ), DIMENSION(:), POINTER :: rhotot, vperiodic
    TYPE( environ_cell ), POINTER :: cell
    !
    INTEGER, POINTER :: env_periodicity
    INTEGER, POINTER :: slab_axis
    REAL( DP ), POINTER :: alat, omega, axis_length
    REAL( DP ), DIMENSION(:), POINTER :: origin
    REAL( DP ), DIMENSION(:,:), POINTER :: axis
    !
    INTEGER :: icor
    REAL(DP) :: fact, const
    REAL(DP) :: dipole(0:3), quadrupole(3)
    REAL(DP) :: tot_charge, tot_dipole(3), tot_quadrupole(3)
    TYPE( environ_density ), TARGET :: local
    CHARACTER( LEN = 80 ) :: sub_name = 'calc_vperiodic'
    !
    CALL start_clock ('calc_vpbc')
    !
    ! ... Aliases and sanity checks
    !
    IF ( .NOT. ASSOCIATED( potential%cell, charges%cell ) ) &
         & CALL errore(sub_name,'Missmatch in domains of potential and charges',1)
    IF ( potential % cell % nnr .NE. oned_analytic % n ) &
         & CALL errore(sub_name,'Missmatch in domains of potential and solver',1)
    cell => potential % cell
    nnr => cell % nnr
    rhotot => charges % of_r
    !
    alat => oned_analytic % alat
    omega => oned_analytic % omega
    env_periodicity => oned_analytic % d
    slab_axis => oned_analytic % axis
    axis_length => oned_analytic % size
    origin => oned_analytic % origin
    axis => oned_analytic % x
    !
    CALL init_environ_density( cell, local )
    vperiodic => local % of_r
    !
    ! ... Compute multipoles of the system wrt the chosen origin
    !
    CALL compute_dipole( nnr, 1, rhotot, origin, dipole, quadrupole )
    !
    tot_charge = dipole(0)
    tot_dipole = dipole(1:3)
    tot_quadrupole = quadrupole
    !
    ! ... Compute quadratic PBC correction
    !
    fact = e2 * tpi / omega
    !
    SELECT CASE ( env_periodicity )
       !
    CASE ( 0 )
       !
       const = madelung(1) * tot_charge / alat * e2 - fact * SUM(tot_quadrupole(:)) / 3.D0
       vperiodic = 0.D0
       DO icor = 1, 3
          vperiodic(:) = vperiodic(:) - tot_charge * axis(icor,:)**2 + &
               & 2.D0 * tot_dipole(icor) * axis(icor,:)
       ENDDO
       vperiodic = fact / 3.D0 * vperiodic + const
       !
    CASE ( 1 )
       !
       CALL errore(sub_name,'Option not yet implemented',1)
       !
    CASE ( 2 )
       !
       const = - pi / 3.D0 * tot_charge / axis_length * e2 - fact * tot_quadrupole(slab_axis)
       vperiodic(:) = - tot_charge * axis(1,:)**2 + 2.D0 * tot_dipole(slab_axis) * axis(1,:)
       vperiodic = fact * vperiodic + const
       !
    CASE ( 3 )
       !
       const = 0.D0
       vperiodic = 0.D0
       !
    CASE DEFAULT
       !
       CALL errore(sub_name,'Unexpected option in dimensionality of PBC correction',1)
       !
    END SELECT
    !
    potential % of_r = potential % of_r + vperiodic
    !
    CALL destroy_environ_density(local)
    !
    CALL stop_clock ('calc_vpbc')
    !
    RETURN
!---------------------------------------------------------------------------
  END SUBROUTINE calc_vperiodic
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
  SUBROUTINE calc_gradvperiodic( oned_analytic, charges, gvtot )
!---------------------------------------------------------------------------
    !
    ! NOTE THAT IN THIS SUBROUTINE THE IONIC DENSITY IS IMPLICIT (AND THUS
    ! SPREAD GAUSSIANS). SINCE THE GRADIENT OF THE CORRECTIVE POTENTIAL DOES
    ! NOT DEPEND ON THE QUADRUPOLE MOMENT OF RHOTOT, IT SHOULD BE INDEPENDENT
    ! ON THE SHAPE OF THE IONIC DENSITY.
    !
    IMPLICIT NONE
    !
    TYPE( oned_analytic_core ), TARGET, INTENT(IN) :: oned_analytic
    TYPE( environ_density ), TARGET, INTENT(IN) :: charges
    TYPE( environ_gradient ), INTENT(INOUT) :: gvtot
    !
    INTEGER,  POINTER  :: nnr
    REAL(DP), DIMENSION(:), POINTER  :: rhotot
    REAL(DP), DIMENSION(:,:), POINTER :: gvperiodic
    TYPE( environ_cell ), POINTER :: cell
    !
    INTEGER, POINTER :: env_periodicity
    INTEGER, POINTER :: slab_axis
    REAL( DP ), POINTER :: omega
    REAL( DP ), DIMENSION(:), POINTER :: origin
    REAL( DP ), DIMENSION(:,:), POINTER :: axis
    !
    INTEGER :: icor
    REAL(DP) :: fact
    REAL(DP) :: dipole(0:3), quadrupole(3)
    REAL(DP) :: tot_charge, tot_dipole(3)
    TYPE( environ_gradient ), TARGET :: glocal
    !
    CHARACTER( LEN = 80 ) :: sub_name = 'calc_gradvperiodic'
    !
    ! ... Aliases and sanity checks
    !
    IF ( .NOT. ASSOCIATED( gvtot%cell, charges%cell ) ) &
         & CALL errore(sub_name,'Missmatch in domains of gradient and charges',1)
    IF ( gvtot % cell % nnr .NE. oned_analytic % n ) &
         & CALL errore(sub_name,'Missmatch in domains of gradient and solver',1)
    !
    cell => gvtot % cell
    nnr => cell % nnr
    rhotot => charges % of_r
    !
    omega => oned_analytic % omega
    env_periodicity => oned_analytic % d
    slab_axis => oned_analytic % axis
    origin => oned_analytic % origin
    axis => oned_analytic % x
    !
    CALL init_environ_gradient( cell, glocal )
    gvperiodic => glocal % of_r
    !
    ! ... Compute dipole of the system with respect to the center of charge
    !
    CALL compute_dipole( nnr, 1, rhotot, origin, dipole, quadrupole )
    !
    tot_charge = dipole(0)
    tot_dipole = dipole(1:3)
    !
    ! ... Compute gradient of periodic images correction
    !
    fact = e2 * fpi / omega
    !
    SELECT CASE ( env_periodicity )
       !
    CASE ( 0 )
       !
       DO icor = 1,3
          gvperiodic(icor,:) = ( tot_dipole(icor) - tot_charge * axis(icor,:) ) / 3.D0
       ENDDO
       !
    CASE ( 1 )
       !
       CALL errore(sub_name,'Option not yet implemented',1)
       !
    CASE ( 2 )
       !
       gvperiodic(slab_axis,:) = tot_dipole(slab_axis) - tot_charge * axis(1,:)
       !
    CASE( 3 )
       !
       gvperiodic = 0.D0
       !
    CASE DEFAULT
       !
       CALL errore(sub_name,'Unexpected option',1)
       !
    END SELECT
    !
    gvperiodic = gvperiodic * fact
    !
    ! ... Sum the periodic contribution to the total gradient of the potential
    !
    gvtot % of_r = gvtot % of_r + gvperiodic
    !
    CALL destroy_environ_gradient( glocal )
    !
    RETURN
    !
!---------------------------------------------------------------------------
  END SUBROUTINE calc_gradvperiodic
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
  SUBROUTINE calc_fperiodic( oned_analytic, natoms, charges, auxiliary, f )
!---------------------------------------------------------------------------
    !
    ! ... Compute contribution to the atomic forces
    !
    IMPLICIT NONE
    !
    TYPE( oned_analytic_core ), TARGET, INTENT(IN) :: oned_analytic
    INTEGER, INTENT(IN)   :: natoms
    TYPE( environ_charges ), TARGET, INTENT(IN) :: charges
    TYPE( environ_density ), INTENT(IN) :: auxiliary
    REAL(DP), INTENT(OUT) :: f( 3, natoms )
    !
    INTEGER, POINTER :: nnr
    INTEGER, DIMENSION(:), POINTER :: ityp
    REAL( DP ), DIMENSION(:), POINTER :: rhotot
    REAL( DP ), DIMENSION(:,:), POINTER :: tau
    !
    INTEGER, POINTER :: env_periodicity
    INTEGER, POINTER :: slab_axis
    REAL( DP ), POINTER :: alat, omega
    REAL( DP ), DIMENSION(:), POINTER :: origin
    !
    INTEGER  :: i
    REAL(DP) :: fact, pos(3)
    REAL(DP) :: dipole(0:3), quadrupole(3)
    REAL(DP) :: tot_charge, tot_dipole(3)
    REAL(DP) :: ftmp( 3, natoms )
    !
    CHARACTER ( LEN = 80 ) :: sub_name = 'calc_fperiodic'
    !
    CALL start_clock ('calc_fpbc')
    !
    ! ... Aliases and sanity checks
    !
    IF ( charges % density % cell % nnr .NE. oned_analytic % n ) &
         & CALL errore(sub_name,'Missmatch in domains of charges and solver',1)
    IF ( natoms .NE. charges % ions % number ) &
         & CALL errore(sub_name,'Missmatch in numbers of atoms passed in input and stored',1)
    !
    nnr => charges % density % cell % nnr
    tau => charges % ions % tau
    ityp => charges % ions % ityp
    !
    alat => oned_analytic % alat
    omega => oned_analytic % omega
    env_periodicity => oned_analytic % d
    slab_axis => oned_analytic % axis
    origin => oned_analytic % origin
    !
    ALLOCATE( rhotot( nnr ) )
    rhotot = charges % density % of_r + auxiliary % of_r
    !
    ! ... Compute multipoles of the system with respect to the chosen origin
    !
    CALL compute_dipole( nnr, 1, rhotot, origin, dipole, quadrupole )
    DEALLOCATE( rhotot )
    !
    tot_charge = dipole(0)
    tot_dipole = dipole(1:3)
    !
    ! ... Interatomic forces, quadrupole is not needed, thus the same
    !     expression holds for point-like and gaussian nuclei
    !
    fact = e2 * fpi / omega
    ftmp = 0.D0
    !
    DO i = 1, natoms
       !
       pos(:) = ( tau( :, i ) - origin( : ) ) * alat
       !
       SELECT CASE ( env_periodicity )
          !
       CASE ( 0 )
          !
          ftmp( :, i ) = ( tot_charge * pos( : ) - tot_dipole( : ) ) / 3.D0
          !
       CASE ( 1 )
          !
          CALL errore(sub_name,'Option not yet implemented',1)
          !
       CASE ( 2 )
          !
          ftmp( slab_axis, i ) = tot_charge * pos( slab_axis ) - tot_dipole( slab_axis )
          !
       CASE ( 3 )
          !
          ftmp = 0.D0
          !
       CASE DEFAULT
          !
          CALL errore(sub_name,'Unexpected',1)
          !
       END SELECT
       !
       ftmp( :, i ) = ftmp( :, i ) * fact * charges % ions % iontype( ityp ( i ) ) % zv
       !
    END DO
    !
    f = f + ftmp
    !
    CALL stop_clock ('calc_fpbc')
    !
    RETURN
!---------------------------------------------------------------------------
  END SUBROUTINE calc_fperiodic
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
END MODULE correction_periodic
!---------------------------------------------------------------------------
