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
  USE kinds,             ONLY : dp
  USE io_global,         ONLY : stdout
  USE mp,                ONLY : mp_sum
  USE mp_bands,          ONLY : intra_bgrp_comm
  USE constants,         ONLY : tpi, pi, sqrtpi
  USE environ_base,      ONLY : verbose, environ_unit, e2, env_periodicity, slab_axis
  USE environ_cell,      ONLY : omega, ntot, domega, alat
  USE environ_ions,      ONLY : zvtot, avg_pos
  USE generate_function, ONLY : generate_axis, generate_gradgaussian, generate_distance
  USE environ_debug,     ONLY : write_cube
  !
  IMPLICIT NONE
  !
  REAL (DP) :: axis_shift(3), axis_lenght
  REAL(DP) :: const
  REAL(DP) :: quadrupole_ions_pc(3), quadrupole_ions_gauss(3)
  REAL(DP) :: tot_charge, tot_dipole(3), tot_quadrupole(3)
  REAL(DP), ALLOCATABLE :: axis(:), distance(:,:)
  REAL(DP), PARAMETER :: madelung(3) = (/ 2.837297479D0, 2.8883D0, 2.885D0 /)
  REAL(DP) :: v0
  !
  SAVE

  PRIVATE

  PUBLIC :: periodic_initbase, periodic_initcell, periodic_initions, periodic_clean, &
            calc_vperiodic, calc_gradvperiodic, calc_eperiodic, calc_fperiodic

CONTAINS
!----------------------------------------------------------------------------
  SUBROUTINE periodic_initbase( nnr )
!----------------------------------------------------------------------------
    !
    ! ... Allocate internal variables
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: nnr
    !
    IF ( env_periodicity .EQ. 2 ) THEN
      IF (ALLOCATED(axis))   DEALLOCATE(axis)
      ALLOCATE( axis(nnr) )
      axis = 0.D0
    ELSE IF ( env_periodicity .EQ. 0 ) THEN
      IF (ALLOCATED(distance))   DEALLOCATE(distance)
      ALLOCATE( distance(3,nnr) )
      distance = 0.D0
    END IF
    !
    RETURN
    !
!---------------------------------------------------------------------------
  END SUBROUTINE periodic_initbase
!---------------------------------------------------------------------------
!----------------------------------------------------------------------------
  SUBROUTINE periodic_initions( nnr, nat, ntyp, ityp, zv, tau, alat_, rhoions )
!----------------------------------------------------------------------------
    !
    ! ... Set up the center of charge of the system, to be used to compute
    !     the dipole of the cell
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN)  :: nnr
    INTEGER, INTENT(IN)  :: nat
    INTEGER, INTENT(IN)  :: ntyp
    INTEGER, INTENT(IN)  :: ityp(nat)
    REAL(DP), INTENT(IN) :: tau(3,nat), zv(ntyp)
    REAL(DP), INTENT(IN) :: alat_
    REAL(DP), INTENT(IN) :: rhoions(nnr)
    !
    INTEGER  :: ia, ip
    REAL(DP) :: pos(3), charge, spread
    REAL(DP) :: dipole(0:3), quadrupole(3)
    !
    quadrupole_ions_pc = 0.D0
    !
    DO ia = 1, nat
       !
       DO ip = 1, 3
          !
          ! NOTE THAT zv IS POSITIVE, THUS HAS THE OPPOSITE SIGN
          !
          quadrupole_ions_pc(ip) = quadrupole_ions_pc(ip) &
                & - zv(ityp(ia))*( ( tau(ip,ia) - avg_pos(ip) )*alat_ )**2
          !
       END DO
       !
    END DO
    !
    CALL compute_dipole( nnr, 1, rhoions, avg_pos, dipole, quadrupole )
    quadrupole_ions_gauss = quadrupole
    !
    RETURN
    !
!----------------------------------------------------------------------------
  END SUBROUTINE periodic_initions
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
  SUBROUTINE periodic_initcell( nnr, at )
!----------------------------------------------------------------------------
    !
    ! ... Generate a linear function of the axis perpendicular to the slab
    !     centered on a specific position
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: nnr
    REAL(DP), INTENT(IN) :: at(3,3)
    !
    axis_shift=avg_pos
    IF ( env_periodicity .EQ. 2 ) THEN
      axis_lenght = at( slab_axis, slab_axis ) * alat
      CALL generate_axis( nnr, slab_axis, axis_shift, axis )
    ELSE IF ( env_periodicity .EQ. 0 ) THEN
      CALL generate_distance( nnr, axis_shift, distance )
    END IF
    !
    RETURN
    !
!---------------------------------------------------------------------------
  END SUBROUTINE periodic_initcell
!---------------------------------------------------------------------------
!----------------------------------------------------------------------------
  SUBROUTINE periodic_clean()
!----------------------------------------------------------------------------
    !
    ! ... Clean up of local variables
    !
    IMPLICIT NONE
    !
    IF (ALLOCATED(axis)) DEALLOCATE(axis)
    IF (ALLOCATED(distance)) DEALLOCATE(distance)
    !
    RETURN
    !
!---------------------------------------------------------------------------
  END SUBROUTINE periodic_clean
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
  SUBROUTINE calc_vperiodic( nnr, nspin, correct_ions, rhotot, vperiodic )
!---------------------------------------------------------------------------
    !
    ! ... Compute potential to correct for pbc along the direction
    !     orthogonal to the slab
    !
    IMPLICIT NONE
    !
    INTEGER,  INTENT(IN)  :: nnr, nspin
    LOGICAL,  INTENT(IN)  :: correct_ions
    REAL(DP), INTENT(IN)  :: rhotot(nnr)
    REAL(DP), INTENT(OUT) :: vperiodic(nnr)
    !
    INTEGER :: icor
    REAL(DP) :: dipole(0:3), quadrupole(3)
    REAL(DP) :: f(3), x0(3), dx(3), xd(3), fact
    !
    REAL(DP), ALLOCATABLE :: vzero(:,:), rhoaux(:,:)
    REAL( DP ) :: ehart, charge
    !
    CALL start_clock ('calc_vpbc')
    !
    vperiodic = 0.D0
    !
    IF ( verbose .GE. 3 ) CALL write_cube( nnr, rhotot,  'rhotot.cube'  )
    IF ( verbose .GE. 2 ) THEN
      ! ... Print the total electrostatic potential before the correction
      ALLOCATE(vzero(nnr,nspin))
      vzero = 0.D0
      ALLOCATE(rhoaux(nnr,nspin))
      rhoaux(:,1) = rhotot
      IF (nspin .EQ. 2) rhoaux(:,2) = 0.D0
      CALL v_h_of_rho_r( rhoaux, ehart, charge, vzero )
      CALL write_cube( nnr, vzero(:,1), 'vfull_pbc.cube')
      DEALLOCATE(rhoaux)
    END IF
    !
    ! ... Compute dipole of the system with respect to the center of charge
    !
    CALL compute_dipole( nnr, 1, rhotot, avg_pos, dipole, quadrupole )
    !
    tot_charge = dipole(0)
    tot_dipole(1:3) = dipole(1:3)
    tot_quadrupole(:) = quadrupole(:)
    !
    ! ... If needed, use the correct (point-like charges) quadrupole moment for the ions
    !
    IF ( correct_ions ) tot_quadrupole = tot_quadrupole - quadrupole_ions_gauss + quadrupole_ions_pc
    !
    x0 = axis_shift*alat
    xd = avg_pos*alat
    !
    ! ... Compute slab correction
    !
    dx = xd - x0
    f(1:3) = tot_dipole(1:3)+tot_charge*dx(1:3)
    fact = e2 * tpi / omega
    IF ( env_periodicity .EQ. 2 ) THEN
      vperiodic(:) = (-tot_charge*axis(:)**2 + 2.D0*f(slab_axis)*axis(:))
      const = - pi / 3.D0 * tot_charge / axis_lenght * e2 - fact * tot_quadrupole(slab_axis)
      vperiodic = fact * vperiodic + const
    ELSE IF ( env_periodicity .EQ. 0 ) THEN
      vperiodic = 0.D0
      DO icor = 1, 3
        vperiodic(:) = vperiodic(:) - tot_charge*distance(icor,:)**2 + 2.D0*f(icor)*distance(icor,:)
      ENDDO
      const = madelung(1) * tot_charge / alat * e2 - fact * SUM(tot_quadrupole(:)) / 3.D0
      vperiodic = fact / 3.D0 * vperiodic + const
    ENDIF
    !
    IF ( verbose .GE. 3 ) CALL write_cube( nnr, vperiodic, 'vperiodic.cube')
    IF ( verbose .GE. 2 ) THEN
      ! ... Print the total electrostatic potential after the correction
      vzero(:,1) = vzero(:,1) + vperiodic
      IF ( env_periodicity .EQ. 2 ) vzero = vzero - fact * ( quadrupole_ions_gauss(slab_axis) - quadrupole_ions_pc(slab_axis) )
      CALL write_cube( nnr, vzero(:,1), 'vfull_free.cube')
      DEALLOCATE(vzero)
    END IF
    !
    CALL stop_clock ('calc_vpbc')
    !
    RETURN
!---------------------------------------------------------------------------
  END SUBROUTINE calc_vperiodic
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
  SUBROUTINE calc_gradvperiodic( nnr, rhotot, gvtot )
!---------------------------------------------------------------------------
    !
    ! NOTE THAT IN THIS SUBROUTINE THE IONIC DENSITY IS IMPLICIT (AND THUS
    ! SPREAD GAUSSIANS). SINCE THE GRADIENT OF THE CORRECTIVE POTENTIAL DOES
    ! NOT DEPEND ON THE QUADRUPOLE MOMENT OF RHOTOT, IT SHOULD BE INDEPENDENT
    ! ON THE SHAPE OF THE IONIC DENSITY, PROVIDED THAT THE IONIC DIPOLE IS
    ! SET TO ZERO.
    !
    IMPLICIT NONE
    !
    INTEGER,  INTENT(IN)  :: nnr
    REAL(DP), INTENT(IN)  :: rhotot(nnr)
    REAL(DP), INTENT(INOUT) :: gvtot(3,nnr)
    !
    INTEGER :: icor
    !
    REAL(DP) :: dipole(0:3), quadrupole(3), f(3), xd(3), x0(3), dx(3)
    REAL(DP), ALLOCATABLE :: gvperiodic(:,:)
    !
    ALLOCATE( gvperiodic( 3, nnr ) )
    gvperiodic = 0.D0

    ! ... Compute dipole of the system with respect to the center of charge

    CALL compute_dipole( nnr, 1, rhotot, avg_pos, dipole, quadrupole )

    x0 = axis_shift*alat
    xd = avg_pos*alat

    ! ... Compute gradient of periodic images correction

    dx = xd - x0
    f(1:3) = dipole(1:3)+dipole(0)*dx(1:3)
    IF ( env_periodicity .EQ. 2 ) THEN
      gvperiodic(slab_axis,:) = e2*tpi*2.D0/omega*(-dipole(0)*axis(:) + f(slab_axis))
    ELSE IF ( env_periodicity .EQ. 0 ) THEN
      DO icor = 1,3
        gvperiodic(icor,:) = e2*tpi*2.D0/3.D0/omega*(-dipole(0)*distance(icor,:) + f(icor))
      ENDDO
    END IF

    ! ... Sum the periodic contribution to the total gradient of the potential

    DO icor = 1,3
      gvtot(icor,:) = gvtot(icor,:) + gvperiodic(icor,:)
    ENDDO

    DEALLOCATE(gvperiodic)
    !
    RETURN

!---------------------------------------------------------------------------
  END SUBROUTINE calc_gradvperiodic
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
  SUBROUTINE calc_eperiodic( nnr, rho, eperiodic )
!---------------------------------------------------------------------------
    !
    ! ... Compute contribution to the total energy of the system
    !
    USE environ_base,  ONLY : vperiodic
    USE environ_ions,  ONLY : nat, ityp, zv, tau
    USE environ_base,  ONLY : env_static_permittivity, rhopol, atomicspread
    USE environ_base,  ONLY : env_external_charges, rhoexternal
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: nnr
    REAL(DP), INTENT(IN)  :: rho(nnr)
    REAL(DP), INTENT(OUT) :: eperiodic
    !
    INTEGER :: ia
    REAL(DP) :: fact, etmp, qtmp
    REAL(DP) :: pos(3), pos2, pos_dip
    !
    CALL start_clock ('calc_epbc')
    !
    eperiodic = 0.5D0 * SUM( vperiodic( : ) * rho( : ) ) * domega
    !
    CALL mp_sum( eperiodic, intra_bgrp_comm )
    !
    ! ... Ionic contribution to the energy, point-like nuclei
    !
    fact = e2*pi/omega
    etmp = 0.D0
    !
    ! NOTE THAT zv HAVE THE OPPOSITE SIGN WRT THE NUCLEAR CHARGES
    !
    DO ia = 1, nat
      pos( : )  = ( tau( :, ia ) - avg_pos(:) ) * alat
      IF ( env_periodicity .EQ. 2 ) THEN
        pos2 = pos( slab_axis )**2
!NO_ION_DIPOLE        pos_dip = tot_dipole( slab_axis ) * pos( slab_axis )
!NO_ION_DIPOLE        etmp = etmp - ( - tot_charge * pos2 + 2 * pos_dip ) * fact * zv ( ityp ( ia ) )
        etmp = etmp + tot_charge * pos2 * fact * zv ( ityp ( ia ) )
      ELSE IF ( env_periodicity .EQ. 0 ) THEN
        pos2 = SUM(pos(:)*pos(:))
!NO_ION_DIPOLE        pos_dip = SUM(pos(:)*tot_dipole(:))
!NO_ION_DIPOLE        etmp = etmp - ( - tot_charge * pos2 + 2 * pos_dip ) * fact * zv ( ityp ( ia ) ) / 3.D0
        etmp = etmp + tot_charge * pos2 * fact * zv ( ityp ( ia ) ) / 3.D0
      END IF
    END DO
    etmp = etmp - const * zvtot / 2.D0
    !
    eperiodic = eperiodic + etmp
    !
    ! ... Polarization correction for gaussian nuclei
    !
    IF ( env_static_permittivity .GT. 1.D0 ) THEN
      IF ( env_periodicity .EQ. 2 ) THEN
        qtmp = (quadrupole_ions_gauss(slab_axis) - quadrupole_ions_pc(slab_axis))
      ELSE IF ( env_periodicity .EQ. 0 ) THEN
        qtmp = SUM(quadrupole_ions_gauss(:) - quadrupole_ions_pc(:))/3.D0
      END IF
      etmp = - fact *  SUM(rhopol(:)) * domega * qtmp
      CALL mp_sum( etmp, intra_bgrp_comm )
      eperiodic = eperiodic + etmp
    END IF
    !
    ! ... If external charges are present, add correction for pbc
    !
    IF ( env_external_charges .GT. 0 ) THEN
      etmp = 0.5D0 * SUM( vperiodic( : ) * rhoexternal( : ) ) * domega
      CALL mp_sum( etmp, intra_bgrp_comm )
      eperiodic = eperiodic + etmp
    END IF
    !
    CALL stop_clock ('calc_epbc')
    !
    RETURN
!---------------------------------------------------------------------------
  END SUBROUTINE calc_eperiodic
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
  SUBROUTINE calc_fperiodic( nnr, nat, f )
!---------------------------------------------------------------------------
    !
    ! ... Compute contribution to the atomic forces
    !
    USE environ_base,  ONLY : atomicspread, vperiodic
    USE environ_ions,  ONLY : ntyp, ityp, zv, tau
    USE environ_base,  ONLY : env_static_permittivity, rhopol, atomicspread
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN)   :: nnr
    INTEGER, INTENT(IN)   :: nat
    REAL(DP), INTENT(OUT) :: f( 3, nat )
    !
    INTEGER  :: ipol, ia
    REAL(DP) :: pos(3), spread, charge, fact
    REAL(DP) :: ftmp( 3, nat )
    REAL(DP), ALLOCATABLE :: grhoion( :, : )
    !
    CALL start_clock ('calc_fpbc')
    !
    ! ... Interatomic forces, point-like nuclei
    !
    fact = e2*tpi*2.D0/omega
    ftmp = 0.D0
    DO ia = 1, nat
      IF ( env_periodicity .EQ. 2 ) THEN
        pos( slab_axis )  = ( tau( slab_axis, ia ) - avg_pos( slab_axis ) ) * alat
        ftmp( slab_axis, ia ) = ( -tot_charge * pos( slab_axis ) + tot_dipole( slab_axis ) ) * fact
      ELSE IF ( env_periodicity .EQ. 0 ) THEN
        pos( : )  = ( tau( :, ia ) - avg_pos(:) ) * alat
        ftmp( :, ia ) = ( - tot_charge * pos( : ) + tot_dipole( : ) ) * fact / 3.D0
      END IF
      ftmp( :, ia ) = ftmp( :, ia ) * zv( ityp ( ia ) )
    END DO
    !
    ! ... Polarization correction for gaussian nuclei !!!! STILL NEED TO TEST IT!!!!
    !
    IF ( env_static_permittivity .GT. 1.D0 ) THEN
      IF ( env_periodicity .EQ. 0 ) THEN !THE CORRECTION FOR SLAB IS ZERO (MAYBE)!
        DO ia = 1, nat
          ftmp( :, ia ) = ftmp( :, ia ) - fact / 2.D0 *  SUM(rhopol(:)) * domega * &
             & ( zv( ityp ( ia ) ) * atomicspread( ityp (ia ) ) / sqrtpi ) / 3.D0
        END DO
      END IF
    END IF
    !
    f = f + ftmp ! plus or minus? plus should be correct
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
