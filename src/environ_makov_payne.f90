!
! Copyright (C) 2007-2008 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! original version by O. Andreussi and N. Marzari
!
MODULE environ_mp
  !
  USE kinds,        ONLY : DP
  USE io_global,    ONLY : stdout
  USE constants,    ONLY : pi, rytoev, au_debye
  USE environ_cell, ONLY : alat, ibrav
  USE environ_ions, ONLY : rhoions, nat, ntyp, ityp, zv, tau
  USE environ_base, ONLY : e2, env_static_permittivity, env_dielectric_regions, rhopol
  !
  IMPLICIT NONE
  !
  REAL(DP) :: el_charge, el_dipole(0:3), el_quadrupole(3)
  REAL(DP) :: sol_charge, sol_dipole(0:3), sol_quadrupole(3) 
  REAL(DP) :: pol_charge, pol_dipole(0:3), pol_quadrupole(3)
  REAL(DP) :: ion_charge, ion_dipole(0:3), ion_quadrupole_pc(3), ion_quadrupole_gauss(3)
  REAL(DP) :: tot_charge, tot_dipole(0:3), tot_quadrupole(3)
  REAL(DP) :: ejellium, eperiodic
  !
  REAL(DP), PARAMETER :: madelung(3) = (/ 2.8373D0, 2.8883D0, 2.885D0 /)
  !
  SAVE
  !
  PRIVATE
  !
  PUBLIC :: environ_makov_payne
  !
CONTAINS
!---------------------------------------------------------------------------
SUBROUTINE environ_makov_payne( nnr, nspin, rhoelec, x0, etot )
!---------------------------------------------------------------------------
  !
  USE solvent,      ONLY : calc_esolvent_of_V
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN)  :: nnr
  INTEGER, INTENT(IN)  :: nspin
  REAL(DP), INTENT(IN) :: rhoelec(nnr,nspin)
  REAL(DP), INTENT(IN) :: x0(3)
  REAL(DP), INTENT(IN) :: etot
  !
  INTEGER :: ia, ip
  REAL(DP) :: aa, bb
  REAL(DP) :: corr_quad
  REAL(DP) :: corr1, corr2
  REAL(DP) :: corr1_pol, corr2_pol
  REAL(DP) :: corr1_tot, corr2_tot
  REAL(DP) :: corr1_int, corr2_int
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: rhopol_of_V
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: rhoaux
  !
  IF ( env_static_permittivity .EQ. 1.D0 &
       .AND. env_dielectric_regions .EQ. 0 ) RETURN
  !
  ! ... Compute the properties of the electronic density
  !
  CALL compute_dipole( nnr, nspin, rhoelec, x0, el_dipole, el_quadrupole )
  !
  ! ... Compute here the properties of the ion density, treating it 
  !     as gaussians. Compute the correction needed to go from gaussian to 
  !     pointlike nuclei (only affect quadrupole moment)
  !
  CALL compute_dipole( nnr, 1, rhoions, x0, ion_dipole, ion_quadrupole_gauss )
  ion_charge = ion_dipole(0)
  !
  ion_quadrupole_pc = 0.D0
  !
  DO ia = 1, nat
     !
     DO ip = 1, 3
        !
        ion_quadrupole_pc(ip) = ion_quadrupole_pc(ip) - &
          zv(ityp(ia))*( ( tau(ip,ia) - x0(ip) )*alat )**2
        !
     END DO
     !
  END DO
  !
  ! ... Need to redo some of the calculations in write_dipole
  !
  sol_dipole  = -el_dipole - ion_dipole
  sol_charge = sol_dipole(0)
  sol_quadrupole = -el_quadrupole  - ion_quadrupole_pc
  !
  corr1 = - madelung(ibrav) / alat * sol_charge**2 / 2.0D0 * e2
  !
  aa = SUM(sol_quadrupole(:))
  bb = SUM(sol_dipole(1:3)*sol_dipole(1:3)) 
  !
  corr2 = ( 2.D0 / 3.D0 * pi )*( sol_charge*aa - bb ) / alat**3 * e2
  !
  ! ... Compute the volume dependent polarization density and its contributions
  !     to the energy : ejellium and eperiodic
  !
  ALLOCATE( rhopol_of_V (nnr) )
  !
  CALL calc_esolvent_of_V( nnr, nspin, rhoelec, ejellium, eperiodic, &
                             & rhopol_of_V )
  !
  ! ... Remove the volume dependent polarization density from the total
  !
  ALLOCATE( rhoaux (nnr) )
  rhoaux = rhopol - rhopol_of_V
  !
  ! ... Compute properties of cleaned polarization density
  !
  CALL compute_dipole( nnr, 1, rhoaux, x0, pol_dipole, pol_quadrupole )
  pol_charge = pol_dipole(0) 
  !
  rhoaux = rhopol - rhopol_of_V + rhoelec(:,1)
  IF ( nspin .EQ. 2 ) rhoaux = rhoaux + rhoelec(:,2)
  !
  ! ... Compute properties of full density (adding the ion contribution 
  !     at a second time as pointlike nuclei)
  !
  CALL compute_dipole( nnr, 1, rhoaux, x0, tot_dipole, tot_quadrupole )
  tot_dipole = tot_dipole + ion_dipole
  tot_charge = tot_dipole(0)
  tot_quadrupole = tot_quadrupole + ion_quadrupole_pc
  !
  DEALLOCATE( rhoaux )
  DEALLOCATE( rhopol_of_V )
  !
  ! ... Now compute Environ corrections
  !
  corr_quad = ( pi * e2 / alat**3 ) * pol_charge * & 
              & SUM( ion_quadrupole_gauss(:) - ion_quadrupole_pc(:) ) / 3.D0  
  !
  WRITE( stdout, '(/,5X,"*******  ENVIRON-MAKOV-PAYNE CORRECTION  *******")' )
  !
  WRITE( stdout, '(/5X,"Dipole moments (with respect to x0):")' )
  WRITE( stdout, '( 5X,"Envir",3F9.4," au (Ha),", 3F9.4," Debye")' ) &
    (-pol_dipole(ip),ip = 1, 3), (-pol_dipole(ip)*au_debye,ip = 1, 3 )
  WRITE( stdout, '( 5X,"Full ",3F9.4," au (Ha),", 3F9.4," Debye")' ) &
    (-tot_dipole(ip),ip = 1, 3), (-tot_dipole(ip)*au_debye,ip = 1, 3 )
  !
  WRITE( stdout, '(/5X,"Environ quadrupole moment",F20.8," a.u. (Ha)")' )  &
    -SUM(pol_quadrupole(:))
  WRITE( stdout, '( 5X,"  FullB quadrupole moment",F20.8," a.u. (Ha)")' )  &
    -SUM(tot_quadrupole(:))
  !
  corr1_pol = -corr1 + corr1 / env_static_permittivity
  aa = SUM( - pol_quadrupole(:) - sol_quadrupole(:) ) / 2.D0 +  &
       & SUM(sol_quadrupole(:)) / 2.D0 / env_static_permittivity
  bb = - sol_dipole(1)*pol_dipole(1) - sol_dipole(2)*pol_dipole(2) - & 
       & sol_dipole(3)*pol_dipole(3) 
  corr2_pol = ( 2.D0 / 3.D0 * pi )*( sol_charge*aa - bb ) / alat**3 * e2 + corr_quad
  !
  WRITE( stdout,'(/5X,"Environ correction     ",F14.8," Ry = ",F6.3, &
     &              " eV (1st order diel, 1/a0)")' ) -corr1_pol, -corr1_pol*rytoev
  WRITE( stdout,'( 5X,"                       ",F14.8," Ry = ",F6.3, &
     &              " eV (2nd order diel, 1/a0^3)")' ) -corr2_pol, -corr2_pol*rytoev
  WRITE( stdout,'( 5X,"                       ",F14.8," Ry = ",F6.3, &
     &              " eV (jellium, 1/a0^3)")' ) -ejellium, -ejellium*rytoev
  WRITE( stdout,'( 5X,"                       ",F14.8," Ry = ",F6.3, &
     &              " eV (periodic solutes, 1/a0^3)")' ) -eperiodic, -eperiodic*rytoev
  WRITE( stdout,'( 5X,"                       ",F14.8," Ry = ",F6.3, &
     &              " eV (total diel)")' ) -corr1_pol-corr2_pol-ejellium-eperiodic, &
                                 & (-corr1_pol-corr2_pol - ejellium - eperiodic)*rytoev
  !
  WRITE( stdout,'(/"!    Envir+Makov-Payne energy  = ",F16.8," Ry")' ) &
      etot - corr1 - corr2 - corr1_pol - corr2_pol - ejellium - eperiodic
  !
  RETURN
  !
!---------------------------------------------------------------------------
END SUBROUTINE environ_makov_payne
!---------------------------------------------------------------------------
END MODULE environ_mp
