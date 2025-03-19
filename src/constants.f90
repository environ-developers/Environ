!----------------------------------------------------------------------------------------
!
! Copyright (C) 2018-2022 ENVIRON (www.quantum-environ.org)
! Copyright (C) 2002-2006 Quantum ESPRESSO (www.quantum-espresso.org)
!
!----------------------------------------------------------------------------------------
!
!     This file is part of Environ version 3.0
!
!     Environ 3.0 is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation, either version 2 of the License, or
!     (at your option) any later version.
!
!     Environ 3.0 is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more detail, either the file
!     `License' in the root directory of the present distribution, or
!     online at <http://www.gnu.org/licenses/>.
!
!----------------------------------------------------------------------------------------
!
! Authors: Oliviero Andreussi (Department of Physics, UNT)
!          Edan Bainglass     (Department of Physics, UNT)
!
!----------------------------------------------------------------------------------------
!>
!!
!----------------------------------------------------------------------------------------
MODULE environ_param
    !------------------------------------------------------------------------------------
    !
    USE env_util_param, ONLY: DP
    !
    !------------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    !------------------------------------------------------------------------------------
    !
    INTEGER, PARAMETER :: nccx = 10
    ! maximum number of countercharge species in the electrolyte
    !
    !------------------------------------------------------------------------------------
    ! Mathematical constants
    !
    REAL(DP), PARAMETER :: pi = 3.14159265358979323846_DP
    REAL(DP), PARAMETER :: tpi = 2.0_DP * pi
    REAL(DP), PARAMETER :: tpi2 = tpi * tpi
    REAL(DP), PARAMETER :: fpi = 4.0_DP * pi
    REAL(DP), PARAMETER :: sqrtpi = 1.77245385090551602729_DP
    REAL(DP), PARAMETER :: sqrt2 = 1.41421356237309504880_DP
    REAL(DP), PARAMETER :: pi_third = pi / 3.D0
    !
    !------------------------------------------------------------------------------------
    ! Physical constants, SI (NIST CODATA 2006), Web Version 5.1
    ! http://physics.nist.gov/constants
    !
    REAL(DP), PARAMETER :: K_BOLTZMANN_SI = 1.3806504E-23_DP ! J K^-1
    REAL(DP), PARAMETER :: ELECTRON_SI = 1.602176487E-19_DP ! C
    REAL(DP), PARAMETER :: HARTREE_SI = 4.35974394E-18_DP ! J
    REAL(DP), PARAMETER :: RYDBERG_SI = HARTREE_SI / 2.0_DP ! J
    REAL(DP), PARAMETER :: BOHR_RADIUS_SI = 0.52917720859E-10_DP ! m
    REAL(DP), PARAMETER :: AMU_SI = 1.660538782E-27_DP ! Kg
    REAL(DP), PARAMETER :: BOHR_RADIUS_CM = BOHR_RADIUS_SI * 100.0_DP
    REAL(DP), PARAMETER :: BOHR_RADIUS_ANGS = BOHR_RADIUS_CM * 1.0E8_DP
    REAL(DP), PARAMETER :: ELECTRONVOLT_SI = 1.602176487E-19_DP
    REAL(DP), PARAMETER :: ELECTRONMASS_SI = 9.10938215E-31_DP
    !
    !------------------------------------------------------------------------------------
    ! Physical constants, atomic units:
    ! AU for "Hartree" atomic units (e = m = hbar = 1)
    ! RY for "Rydberg" atomic units (e^2=2, m=1/2, hbar=1)
    !
    REAL(DP), PARAMETER :: K_BOLTZMANN_AU = K_BOLTZMANN_SI / HARTREE_SI
    REAL(DP), PARAMETER :: K_BOLTZMANN_RY = K_BOLTZMANN_SI / RYDBERG_SI
    !
    !------------------------------------------------------------------------------------
    ! Unit conversion factors: energy and masses
    !
    REAL(DP), PARAMETER :: AUTOEV = HARTREE_SI / ELECTRONVOLT_SI
    REAL(DP), PARAMETER :: RYTOEV = AUTOEV / 2.0_DP
    REAL(DP), PARAMETER :: AMU_AU = AMU_SI / ELECTRONMASS_SI
    REAL(DP), PARAMETER :: AMU_RY = AMU_AU / 2.0_DP
    !
    !------------------------------------------------------------------------------------
    ! Unit conversion factors: 1 debye = 10^-18 esu*cm
    !                                  = 3.3356409519*10^-30 C*m
    !                                  = 0.208194346 e*A
    !
    ! ( 1 esu = (0.1/c) Am, c=299792458 m/s)
    !
    REAL(DP), PARAMETER :: DEBYE_SI = 3.3356409519_DP * 1.0E-30_DP ! C*m
    REAL(DP), PARAMETER :: AU_DEBYE = ELECTRON_SI * BOHR_RADIUS_SI / DEBYE_SI
    !
    REAL(DP), PARAMETER :: eV_to_kelvin = ELECTRONVOLT_SI / K_BOLTZMANN_SI
    REAL(DP), PARAMETER :: ry_to_kelvin = RYDBERG_SI / K_BOLTZMANN_SI
    !
    !------------------------------------------------------------------------------------
    ! Zero up to a given accuracy
    !
    REAL(DP), PARAMETER :: eps8 = 1.0E-8_DP
    !
    !------------------------------------------------------------------------------------
    !
    REAL(DP), PARAMETER :: AVOGADRO = 6.02214129D+23
    !
    !------------------------------------------------------------------------------------
    ! Madelung constants
    !
    REAL(DP), PARAMETER :: madelung(3) = (/2.837297479D0, 2.8883D0, 2.885D0/)
    !
    !------------------------------------------------------------------------------------
    ! Rydberg vs Hartree units (e2 = 2 in Ry)
    !
    REAL(DP) :: e2 = 2.D0
    !
    !------------------------------------------------------------------------------------
END MODULE environ_param
!----------------------------------------------------------------------------------------
