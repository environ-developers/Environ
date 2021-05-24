! Copyright (C) 2018 ENVIRON (www.quantum-environment.org)
!
!    This file is part of Environ version 1.1
!
!    Environ 1.1 is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 2 of the License, or
!    (at your option) any later version.
!
!    Environ 1.1 is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more detail, either the file
!    `License' in the root directory of the present distribution, or
!    online at <http://www.gnu.org/licenses/>.
!
!----------------------------------------------------------------------------------------
!
! Authors: Oliviero Andreussi (Department of Physics, UNT)
!          Francesco Nattino  (THEOS and NCCR-MARVEL, EPFL)
!          Ismaila Dabo       (DMSE, Penn State)
!          Nicola Marzari     (THEOS and NCCR-MARVEL, EPFL)
!
!----------------------------------------------------------------------------------------
!>
!! This module contains all the main variables needed for the environ module. This
!! includes the control and debug variables, all physical and numerical parameters,
!! and the contributions to the energy and to the potential.
!!
!----------------------------------------------------------------------------------------
MODULE base_environ
    !------------------------------------------------------------------------------------
    !
    USE environ_param, ONLY: DP
    !
    USE types_cell
    USE types_physical
    USE types_core, ONLY: core_container
    USE types_representation, ONLY: environ_density
    !
    !------------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    SAVE
    !
    !------------------------------------------------------------------------------------
    ! Environ execution parameters
    !
    LOGICAL :: environ_restart, update_venviron
    REAL(DP) :: environ_thr
    INTEGER :: environ_nskip
    !
    !------------------------------------------------------------------------------------
    ! Control flags
    !
    LOGICAL :: ltddfpt, lsolvent, lelectrostatic, lsoftsolvent, lsoftelectrolyte, &
               lsoftcavity, lrigidsolvent, lrigidelectrolyte, lrigidcavity, &
               lcoredensity, lsmearedions, lgradient
    !
    !------------------------------------------------------------------------------------
    ! Internal parameters for mapping between environment and system cell
    !
    LOGICAL :: ldoublecell
    TYPE(environ_mapping) :: mapping
    TYPE(environ_cell), TARGET :: system_cell
    TYPE(environ_cell), POINTER :: environment_cell
    TYPE(environ_ions), TARGET :: system_ions, environment_ions
    TYPE(environ_electrons), TARGET :: system_electrons, environment_electrons
    !
    !------------------------------------------------------------------------------------
    !
    LOGICAL :: lexternals
    TYPE(environ_externals), TARGET :: externals ! external charges
    TYPE(environ_charges), TARGET :: system_charges, environment_charges
    TYPE(environ_system), TARGET :: system_system, environment_system
    !
    !------------------------------------------------------------------------------------
    ! Response properties
    !
    TYPE(environ_electrons), TARGET :: system_response_electrons, &
                                       environment_response_electrons
    TYPE(environ_charges), TARGET :: system_response_charges, &
                                     environment_response_charges
    !
    !------------------------------------------------------------------------------------
    ! Details of the continuum interface
    !
    LOGICAL :: lboundary
    TYPE(core_container), TARGET :: derivatives
    TYPE(environ_boundary), TARGET :: solvent
    !
    !------------------------------------------------------------------------------------
    ! Dielectric parameters (solvent)
    !
    LOGICAL :: lstatic, loptical, ldielectric
    REAL(DP) :: env_static_permittivity, env_optical_permittivity
    TYPE(environ_dielectric) :: static, optical
    !
    !------------------------------------------------------------------------------------
    ! Ionic countercharge parameters
    !
    LOGICAL :: lelectrolyte
    INTEGER :: env_electrolyte_ntyp
    TYPE(environ_electrolyte) :: electrolyte
    !
    !------------------------------------------------------------------------------------
    ! Semiconductor parameters
    !
    LOGICAL :: lsemiconductor, louterloop
    TYPE(environ_semiconductor) :: semiconductor
    !
    !------------------------------------------------------------------------------------
    ! Cavitation energy parameters
    !
    LOGICAL :: lsurface
    REAL(DP) :: env_surface_tension
    !
    !------------------------------------------------------------------------------------
    ! PV term parameters
    !
    LOGICAL :: lvolume
    REAL(DP) :: env_pressure
    !
    !------------------------------------------------------------------------------------
    ! Confinement potential parameters
    !
    LOGICAL :: lconfine
    REAL(DP) :: env_confine
    !
    !------------------------------------------------------------------------------------
    ! Periodicity correction parameters
    !
    LOGICAL :: lperiodic
    !
    !------------------------------------------------------------------------------------
    ! Temporary parameters
    !
    INTEGER :: nrep
    INTEGER :: niter ! stores the iteration of environ for debugging purposes
    !
    !------------------------------------------------------------------------------------
    ! Computed physical variables
    !
    REAL(DP) :: deenviron, eelectrostatic, esurface, evolume, eelectrolyte, &
                econfine, potential_shift
    !
    TYPE(environ_density) :: vzero, velectrostatic, vreference, dvtot, vconfine, &
                             vsoftcavity
    !
    !------------------------------------------------------------------------------------
    ! Keeping imports private
    !
    PRIVATE :: DP, environ_cell, environ_mapping, environ_iontype, environ_ions, &
               environ_electrons, environ_externals, environ_charges, environ_system, &
               environ_boundary, environ_dielectric, environ_ioncctype, &
               environ_electrolyte, environ_semiconductor, core_container, &
               environ_density
    !
    !------------------------------------------------------------------------------------
END MODULE base_environ
!----------------------------------------------------------------------------------------
