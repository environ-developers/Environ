!
! Copyright (C) 2007-2008 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! The different modules in this file contain all the Environ related variables
! that need to be passed in and out. All module-specific variables are declared
! inside the appropriate modules.
!
! original version by O. Andreussi, I. Dabo, and N. Marzari (MIT)
!
!----------------------------------------------------------------------------
MODULE environ_base
  !--------------------------------------------------------------------------
  !
  ! ... this module contains all the main variables needed for the
  ! ... environ module. This include the control and debug variables,
  ! ... all physical and numerical parameters, and
  ! ... the contributions to the energy and to the potential.
  !
  USE environ_types
  !
  SAVE
  !
  ! Global parameters
  !
  LOGICAL ::                        &
       oldenviron
  LOGICAL ::                        &
       environ_restart
  LOGICAL ::                        &
       update_venviron
  REAL (KIND=DP) ::                 &
       environ_thr
  INTEGER ::                        &
       environ_nskip
  REAL (KIND=DP) ::                 &
       e2
  !
  ! Control flags
  !
  LOGICAL ::                          &
       lsolvent
  LOGICAL ::                          &
       lelectrostatic
  LOGICAL ::                          &
       lsoftcavity
  LOGICAL ::                          &
       lrigidcavity
  LOGICAL ::                          &
       lcoredensity
  LOGICAL ::                          &
       lsmearedions
  !
  ! Internal cell parameters
  !
  TYPE( environ_cell ), TARGET ::     &
       cell
  !
  ! Internal parameters of ions
  !
  TYPE( environ_ions ), TARGET ::  &
       ions
  !
  ! Internal parameters of electrons
  !
  TYPE( environ_electrons ), TARGET ::  &
       electrons
  !
  ! Internal parameters of external charges
  !
  LOGICAL ::                            &
       lexternals
  TYPE( environ_externals ), TARGET ::  &
       externals
  !
  ! Internal paramters of total charges
  !
  TYPE( environ_charges ), TARGET :: &
       charges
  !
  ! Details of the selected system
  !
  TYPE( environ_system ), TARGET ::   &
       system
  !
  ! Details of the continuum interface
  !
  TYPE( environ_boundary ), TARGET :: &
       solvent
  !
  ! Dielectric parameters (solvent)
  !
  LOGICAL ::                        &
       lstatic
  LOGICAL ::                        &
       loptical
  LOGICAL ::                        &
       ldielectric
  REAL (KIND=DP) ::                 &
       env_static_permittivity,     &
       env_optical_permittivity
  TYPE( environ_dielectric ) ::     &
       static,                      &
       optical
  !
  ! Ionic countercharge parameters
  !
  LOGICAL ::                        &
       lelectrolyte
  INTEGER ::                        &
       env_ioncc_ntyp
  TYPE( environ_electrolyte ) ::    &
       electrolyte
  !
  ! Cavitation energy parameters
  !
  LOGICAL ::                        &
       lsurface
  REAL (KIND=DP) ::                 &
       env_surface_tension
  !
  ! PV term parameters
  !
  LOGICAL ::                        &
       lvolume
  REAL (KIND=DP) ::                 &
       env_pressure
  !
  ! Periodicity correction parameters
  !
  LOGICAL ::                        &
       lperiodic
  INTEGER ::                        &
       env_periodicity
  INTEGER ::                        &
       slab_axis
  !
  ! Temporary parameters
  !
  LOGICAL ::                        &
       add_jellium
  INTEGER ::                        &
       nrep
  !
  ! Computed physical variables
  !
  REAL (KIND=DP) ::                 &
       deenviron,                   &
       eelectrostatic,              &
       ecavity,                     &
       epressure,                   &
       eelectrolyte
  TYPE ( environ_density ) ::       &
       vzero,                       &
       velectrostatic,              &
       vreference,                  &
       vsoftcavity,                 &
       vsoftelectrolyte
  TYPE ( environ_gradient ) ::      &
       gelectrostatic
  ! gradient of the total electrostatic potential
  ! (elec + pol + ions) for the TDDFPT calculation
  !
  !--------------------------------------------------------------------------
END MODULE environ_base
!----------------------------------------------------------------------------
