!
! Copyright (C) 2007-2008 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Module to include a real-space correction of the potential in order to
! remove periodic boundary conditions on a partially periodic system.
! The correction is computed for the time being using the planar average approximation and
! is of quadratic nature: the first order proportional to the dipole of the
! system along the direction perpendicular to the slab, the second order
! proportional to the total charge of the system.
!
! original version by O. Andreussi and N. Marzari
!
!----------------------------------------------------------------------------
MODULE periodic
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
  PUBLIC :: calc_vperiodic, calc_eperiodic, calc_fperiodic !, calc_gradvperiodic
  !
CONTAINS
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
    TYPE( environ_charges ), TARGET, INTENT(IN) :: charges
    TYPE( environ_density ), INTENT(INOUT) :: potential
    !
    INTEGER, POINTER :: nnr
    REAL( DP ), DIMENSION(:), POINTER :: avg_pos
    REAL( DP ), DIMENSION(:), POINTER :: rhotot, vperiodic
    TYPE( environ_cell ), POINTER :: cell
    !
    INTEGER, POINTER :: env_periodicity
    INTEGER, POINTER :: slab_axis
    REAL( DP ), POINTER :: alat, omega, axis_length
    REAL( DP ), DIMENSION(:), POINTER :: axis_shift
    REAL( DP ), DIMENSION(:,:), POINTER :: axis
    !
    INTEGER :: icor
    REAL(DP) :: dipole(0:3), quadrupole(3)
    REAL(DP) :: f(3), x0(3), dx(3), xd(3), fact, const
    REAL(DP) :: tot_charge, tot_dipole(3), tot_quadrupole(3)
    TYPE( environ_density ), TARGET :: local
    !
    CHARACTER( LEN = 80 ) :: sub_name = 'calc_vperiodic'
    !
    CALL start_clock ('calc_vpbc')
    !
    ! ... Aliases and sanity checks
    !
    IF ( .NOT. ASSOCIATED( potential%cell, charges%density%cell ) ) &
         & CALL errore(sub_name,'Missmatch in domains of potential and charges',1)
    IF ( potential % cell % nnr .NE. oned_analytic % n ) &
         & CALL errore(sub_name,'Missmatch in domains of potential and solver',1)
    cell => potential % cell
    nnr => cell % nnr
    avg_pos => charges % ions % center
    rhotot => charges % density % of_r
    !
    alat => oned_analytic % alat
    omega => oned_analytic % omega
    env_periodicity => oned_analytic % d
    slab_axis => oned_analytic % axis
    axis_length => oned_analytic % size
    axis_shift => oned_analytic % origin
    axis => oned_analytic % x
    !
    CALL init_environ_density( cell, local )
    vperiodic => local % of_r
    !
    ! ... Compute dipole of the system with respect to the center of charge
    !
    CALL compute_dipole( nnr, 1, rhotot, avg_pos, dipole, quadrupole )
    !
    tot_charge = dipole(0)
    tot_dipole = dipole(1:3)
    tot_quadrupole = quadrupole
    !
    x0 = avg_pos*alat ! axis_shift*alat
    xd = avg_pos*alat
    !
    ! ... Compute slab correction
    !
    dx = xd - x0
    f(1:3) = tot_dipole(1:3) + tot_charge * dx(1:3)
    fact = e2 * tpi / omega
    !
    IF ( env_periodicity .EQ. 2 ) THEN
      vperiodic(:) = (-tot_charge*axis(1,:)**2 + 2.D0*f(slab_axis)*axis(1,:))
      const = - pi / 3.D0 * tot_charge / axis_length * e2 - fact * tot_quadrupole(slab_axis)
      vperiodic = fact * vperiodic + const
    ELSE IF ( env_periodicity .EQ. 0 ) THEN
      vperiodic = 0.D0
      DO icor = 1, 3
        vperiodic(:) = vperiodic(:) - tot_charge*axis(icor,:)**2 + 2.D0*f(icor)*axis(icor,:)
      ENDDO
      const = madelung(1) * tot_charge / alat * e2 - fact * SUM(tot_quadrupole(:)) / 3.D0
      vperiodic = fact / 3.D0 * vperiodic + const
    ENDIF
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
!!---------------------------------------------------------------------------
!  SUBROUTINE calc_gradvperiodic( nnr, rhotot, gvtot )
!!---------------------------------------------------------------------------
!    !
!    ! NOTE THAT IN THIS SUBROUTINE THE IONIC DENSITY IS IMPLICIT (AND THUS
!    ! SPREAD GAUSSIANS). SINCE THE GRADIENT OF THE CORRECTIVE POTENTIAL DOES
!    ! NOT DEPEND ON THE QUADRUPOLE MOMENT OF RHOTOT, IT SHOULD BE INDEPENDENT
!    ! ON THE SHAPE OF THE IONIC DENSITY, PROVIDED THAT THE IONIC DIPOLE IS
!    ! SET TO ZERO.
!    !
!    IMPLICIT NONE
!    !
!    INTEGER,  INTENT(IN)  :: nnr
!    REAL(DP), INTENT(IN)  :: rhotot(nnr)
!    REAL(DP), INTENT(INOUT) :: gvtot(3,nnr)
!    !
!    INTEGER :: icor
!    !
!    REAL(DP) :: dipole(0:3), quadrupole(3), f(3), xd(3), x0(3), dx(3)
!    REAL(DP), ALLOCATABLE :: gvperiodic(:,:)
!    !
!    ALLOCATE( gvperiodic( 3, nnr ) )
!    gvperiodic = 0.D0
!
!    ! ... Compute dipole of the system with respect to the center of charge
!
!    CALL compute_dipole( nnr, 1, rhotot, avg_pos, dipole, quadrupole )
!
!    x0 = axis_shift*alat
!    xd = avg_pos*alat
!
!    ! ... Compute gradient of periodic images correction
!
!    dx = xd - x0
!    f(1:3) = dipole(1:3)+dipole(0)*dx(1:3)
!    IF ( env_periodicity .EQ. 2 ) THEN
!      gvperiodic(slab_axis,:) = e2*tpi*2.D0/omega*(-dipole(0)*axis(:) + f(slab_axis))
!    ELSE IF ( env_periodicity .EQ. 0 ) THEN
!      DO icor = 1,3
!        gvperiodic(icor,:) = e2*tpi*2.D0/3.D0/omega*(-dipole(0)*distance(icor,:) + f(icor))
!      ENDDO
!    END IF
!
!    ! ... Sum the periodic contribution to the total gradient of the potential
!
!    DO icor = 1,3
!      gvtot(icor,:) = gvtot(icor,:) + gvperiodic(icor,:)
!    ENDDO
!
!    DEALLOCATE(gvperiodic)
!    !
!    RETURN
!
!!---------------------------------------------------------------------------
!  END SUBROUTINE calc_gradvperiodic
!!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
  SUBROUTINE calc_eperiodic( oned_analytic, charges, potential, energy )
!---------------------------------------------------------------------------
    !
    ! ... HERE ONLY COMPUTES THE CORRECTION DUE TO GAUSSIAN NUCLEI, NEED TO MOVE
    !     OUT OF HERE AND REMOVE THIS REOUTINE
    !
    IMPLICIT NONE
    !
    TYPE( oned_analytic_core ), TARGET, INTENT(IN) :: oned_analytic
    TYPE( environ_charges ), TARGET, INTENT(INOUT) :: charges
    TYPE( environ_density ), INTENT(IN) :: potential !!! we do not need this !!!
    REAL(DP), INTENT(INOUT) :: energy
    !
    REAL(DP), POINTER :: omega
    !
    REAL(DP) :: etmp
    REAL(DP) :: tot_charge
    !
    CHARACTER( LEN = 80 ) :: sub_name = 'calc_eperiodic'
    !
    CALL start_clock ('calc_epbc')
    !
    ! ... Aliases
    !
    omega => charges % density % cell % omega
    !
    ! ... Correct for point-like ions
    !
    tot_charge = integrate_environ_density( charges % density )
    !
    etmp = charges % ions % quadrupole_correction * tot_charge * e2 * tpi / omega
    !
    energy = energy + etmp
    !
    CALL stop_clock ('calc_epbc')
    !
    RETURN
!---------------------------------------------------------------------------
  END SUBROUTINE calc_eperiodic
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
  SUBROUTINE calc_fperiodic( oned_analytic, natoms, charges, f )
!---------------------------------------------------------------------------
    !
    ! ... Compute contribution to the atomic forces
    !
    IMPLICIT NONE
    !
    TYPE( oned_analytic_core ), TARGET, INTENT(IN) :: oned_analytic
    INTEGER, INTENT(IN)   :: natoms
    TYPE( environ_charges ), TARGET, INTENT(IN) :: charges
    REAL(DP), INTENT(OUT) :: f( 3, natoms )
    !
    INTEGER, POINTER :: nnr
    INTEGER, DIMENSION(:), POINTER :: ityp
    REAL( DP ), DIMENSION(:), POINTER :: avg_pos
    REAL( DP ), DIMENSION(:), POINTER :: rhotot
    REAL( DP ), DIMENSION(:,:), POINTER :: tau
    !
    INTEGER, POINTER :: env_periodicity
    INTEGER, POINTER :: slab_axis
    REAL( DP ), POINTER :: alat, omega
    REAL( DP ), DIMENSION(:), POINTER :: axis_shift
    !
    INTEGER  :: i
    REAL(DP) :: dipole(0:3), quadrupole(3)
    REAL(DP) :: tot_charge, tot_dipole(3), tot_quadrupole(3)
    REAL(DP) :: pos(3), spread, charge, fact
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
    rhotot => charges % density % of_r
    avg_pos => charges % ions % center
    tau => charges % ions % tau
    ityp => charges % ions % ityp
    !
    alat => oned_analytic % alat
    omega => oned_analytic % omega
    env_periodicity => oned_analytic % d
    slab_axis => oned_analytic % axis
    axis_shift => oned_analytic % origin
    !
    ! ... Compute dipole of the system with respect to the center of charge
    !
    CALL compute_dipole( nnr, 1, rhotot, avg_pos, dipole, quadrupole )
    !
    tot_charge = dipole(0)
    tot_dipole = dipole(1:3)
    tot_quadrupole = quadrupole
    !
    ! ... Interatomic forces, point-like nuclei
    !
    ftmp = 0.D0
    DO i = 1, natoms
       pos(:) = ( tau( :, i ) - avg_pos( : ) ) * alat
       fact = charges % ions % iontype( ityp ( i ) ) % zv * e2 * fpi / omega
      IF ( env_periodicity .EQ. 2 ) THEN
        ftmp( slab_axis, i ) = tot_charge * pos( slab_axis ) - tot_dipole( slab_axis )
      ELSE IF ( env_periodicity .EQ. 0 ) THEN
        ftmp( :, i ) = ( tot_charge * pos( : ) - tot_dipole( : ) ) / 3.D0
      END IF
      ftmp( :, i ) = ftmp( :, i ) * fact
    END DO
!    !
!    ! ... Polarization correction for gaussian nuclei !!!! STILL NEED TO TEST IT!!!!
!    !
!    IF ( env_static_permittivity .GT. 1.D0 ) THEN
!      IF ( env_periodicity .EQ. 0 ) THEN !THE CORRECTION FOR SLAB IS ZERO (MAYBE)!
!        DO ia = 1, nat
!          ftmp( :, ia ) = ftmp( :, ia ) - fact / 2.D0 *  SUM(rhopol(:)) * domega * &
!             & ( zv( ityp ( ia ) ) * atomicspread( ityp (ia ) ) / sqrtpi ) / 3.D0
!        END DO
!      END IF
!    END IF
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
END MODULE periodic
!---------------------------------------------------------------------------
