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
! Authors: Oliviero Andreussi (Department of Physics, UNT)
!          Francesco Nattino  (THEOS and NCCR-MARVEL, EPFL)
!          Nicola Marzari     (THEOS and NCCR-MARVEL, EPFL)
!
!----------------------------------------------------------------------------------------
!>
!! This module contains the main control and parameter variables from QE Modules,
!! the definitions of Environ derived data types and the routines to handle the
!! basic derived data types (cell, density, gradient, hessian, electrons, system)
!!
!----------------------------------------------------------------------------------------
MODULE physical_types
    !------------------------------------------------------------------------------------
    !
    USE modules_constants, ONLY: DP
    !
    USE cell_types, ONLY: environ_cell
    USE core_types, ONLY: boundary_core
    USE representation_types
    !
    ! BACKWARD COMPATIBILITY
    ! Compatible with QE-5.X QE-6.1.X QE-6.2.X
    ! USE control_flags, ONLY: tddfpt
    ! Compatible with QE-6.3.X and QE-GIT \
    ! END BACKWARD COMPATIBILITY
    !
    !------------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    TYPE environ_iontype
        !--------------------------------------------------------------------------------
        !
        INTEGER :: index
        INTEGER :: atmnum
        CHARACTER(LEN=3) :: label
        REAL(DP) :: zv
        REAL(DP) :: atomicspread
        REAL(DP) :: corespread
        REAL(DP) :: solvationrad
        !
        !--------------------------------------------------------------------------------
    END TYPE environ_iontype
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    TYPE environ_ions
        !--------------------------------------------------------------------------------
        !
        LOGICAL :: initialized = .FALSE.
        LOGICAL :: update = .FALSE.
        INTEGER :: number = 0
        REAL(DP) :: center(3)
        REAL(DP) :: alat
        !
        !--------------------------------------------------------------------------------
        ! Specifications of point-like ions
        !
        INTEGER :: ntyp = 0
        INTEGER, ALLOCATABLE :: ityp(:)
        REAL(DP), POINTER :: tau(:, :)
        TYPE(environ_iontype), ALLOCATABLE :: iontype(:)
        !
        !--------------------------------------------------------------------------------
        ! Parameters of the fictitious gaussian ionic density
        ! needed by electrostatic calculations
        !
        LOGICAL :: use_smeared_ions = .FALSE.
        TYPE(environ_functions), ALLOCATABLE :: smeared_ions(:)
        TYPE(environ_density) :: density
        !
        !--------------------------------------------------------------------------------
        ! Parameters of the density of core electrons
        !
        LOGICAL :: use_core_electrons = .FALSE.
        TYPE(environ_functions), ALLOCATABLE :: core_electrons(:)
        TYPE(environ_density) :: core
        TYPE(environ_density), ALLOCATABLE :: vloc(:)
        !
        REAL(DP) :: charge = 0.0_DP
        REAL(DP) :: quadrupole_correction
        REAL(DP) :: selfenergy_correction
        REAL(DP) :: dipole(3)
        REAL(DP) :: quadrupole_pc(3)
        REAL(DP) :: quadrupole_gauss(3)
        !
        !--------------------------------------------------------------------------------
    END TYPE environ_ions
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    TYPE environ_electrons
        !--------------------------------------------------------------------------------
        !
        LOGICAL :: update = .FALSE.
        LOGICAL :: initialized = .FALSE.
        INTEGER :: number = 0
        ! BACKWARD COMPATIBILITY
        ! Compatible with QE-6.0 QE-6.1.X QE-6.2.X QE-6.3.X
        ! INTEGER :: nspin = 1
        ! Compatible with QE-6.4.X QE-GIT
        !
        ! END BACKWARD COMPATIBILITY
        !
        TYPE(environ_density) :: density
        REAL(DP) :: charge = 0.0_DP
        !
        !--------------------------------------------------------------------------------
    END TYPE environ_electrons
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    TYPE environ_externals
        !--------------------------------------------------------------------------------
        !
        LOGICAL :: update = .FALSE.
        LOGICAL :: initialized = .FALSE.
        INTEGER :: number = 0
        !
        TYPE(environ_functions), ALLOCATABLE :: functions(:)
        TYPE(environ_density) :: density
        REAL(DP) :: charge = 0.0_DP
        !
        !--------------------------------------------------------------------------------
    END TYPE environ_externals
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    TYPE environ_charges
        !--------------------------------------------------------------------------------
        !
        LOGICAL :: include_ions = .FALSE.
        TYPE(environ_ions), POINTER :: ions => NULL()
        !
        LOGICAL :: include_electrons = .FALSE.
        TYPE(environ_electrons), POINTER :: electrons => NULL()
        !
        LOGICAL :: include_externals = .FALSE.
        TYPE(environ_externals), POINTER :: externals => NULL()
        !
        LOGICAL :: include_dielectric = .FALSE.
        TYPE(environ_dielectric), POINTER :: dielectric => NULL()
        !
        LOGICAL :: include_electrolyte = .FALSE.
        TYPE(environ_electrolyte), POINTER :: electrolyte => NULL()
        !
        LOGICAL :: include_semiconductor = .FALSE.
        TYPE(environ_semiconductor), POINTER :: semiconductor => NULL()
        !
        !--------------------------------------------------------------------------------
        ! Total smooth free charge
        !
        INTEGER :: number = 0
        REAL(DP) :: charge = 0.0_DP
        TYPE(environ_density) :: density
        LOGICAL :: initialized = .FALSE.
        !
        !--------------------------------------------------------------------------------
    END TYPE environ_charges
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    TYPE environ_system
        !--------------------------------------------------------------------------------
        !
        LOGICAL :: update = .FALSE.
        INTEGER :: ntyp
        INTEGER :: dim
        INTEGER :: axis
        REAL(DP) :: pos(3)
        REAL(DP) :: width
        !
        TYPE(environ_ions), POINTER :: ions
        !
        !--------------------------------------------------------------------------------
    END TYPE environ_system
    !------------------------------------------------------------------------------------
    !>
    !! #TODO rename scaled variable
    !!
    !------------------------------------------------------------------------------------
    TYPE environ_boundary
        !--------------------------------------------------------------------------------
        !
        CHARACTER(LEN=80) :: label ! Boundary label
        CHARACTER(LEN=80) :: mode ! Choice of the interface
        INTEGER :: update_status = 0
        LOGICAL :: initialized = .FALSE.
        !
        !--------------------------------------------------------------------------------
        ! Parameters for the electrons-dependent interface
        !
        LOGICAL :: need_electrons
        TYPE(environ_electrons), POINTER :: electrons
        !
        !--------------------------------------------------------------------------------
        ! Parameters for the ions-dependent interface
        !
        LOGICAL :: need_ions
        TYPE(environ_ions), POINTER :: ions
        !
        !--------------------------------------------------------------------------------
        ! Parameters for the system-dependent interface
        !
        LOGICAL :: need_system
        TYPE(environ_system), POINTER :: system
        !
        !--------------------------------------------------------------------------------
        ! #TODO change scaled to interface
        TYPE(environ_density) :: scaled ! scaled switching function of interface
        ! varying from 1 (QM region) to 0 (environment region)
        !
        INTEGER :: deriv = 0
        TYPE(environ_gradient) :: gradient
        TYPE(environ_density) :: laplacian
        TYPE(environ_density) :: dsurface
        TYPE(environ_hessian) :: hessian
        !
        TYPE(boundary_core), POINTER :: core
        !
        !--------------------------------------------------------------------------------
        ! Global properties of the boundary
        !
        REAL(DP) :: volume
        REAL(DP) :: surface
        !
        !--------------------------------------------------------------------------------
        ! Components needed for boundary of density
        !
        INTEGER :: type_
        REAL(DP) :: rhomax, rhomin, fact
        REAL(DP) :: rhozero, deltarho, tbeta
        REAL(DP) :: const
        TYPE(environ_density) :: density
        !
        TYPE(environ_density) :: dscaled ! #TODO dinterface
        TYPE(environ_density) :: d2scaled
        !
        !--------------------------------------------------------------------------------
        ! Components needed for boundary of functions
        !
        REAL(DP) :: alpha ! solvent-dependent scaling factor
        REAL(DP) :: softness ! sharpness of the interface
        TYPE(environ_functions), ALLOCATABLE :: soft_spheres(:)
        !
        !--------------------------------------------------------------------------------
        !
        TYPE(environ_functions) :: simple ! Components needed for boundary of system
        !
        !--------------------------------------------------------------------------------
        ! Components needed for solvent-aware boundary
        !
        LOGICAL :: solvent_aware
        TYPE(environ_functions) :: solvent_probe
        REAL(DP) :: filling_threshold, filling_spread
        !
        TYPE(environ_density) :: local
        TYPE(environ_density) :: probe
        TYPE(environ_density) :: filling
        TYPE(environ_density) :: dfilling
        !
        !--------------------------------------------------------------------------------
        ! Components needed for field-aware boundary
        !
        LOGICAL :: field_aware
        REAL(DP) :: field_factor, charge_asymmetry, field_max, field_min
        !
        TYPE(environ_density) :: normal_field
        REAL(DP), ALLOCATABLE :: ion_field(:)
        TYPE(environ_functions), ALLOCATABLE :: local_spheres(:)
        TYPE(environ_density), ALLOCATABLE :: dion_field_drho(:)
        REAL(DP), ALLOCATABLE :: partial_of_ion_field(:, :, :)
        !
        !--------------------------------------------------------------------------------
    END TYPE environ_boundary
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    TYPE environ_dielectric
        !--------------------------------------------------------------------------------
        !
        LOGICAL :: update = .FALSE.
        LOGICAL :: initialized = .FALSE.
        !
        !--------------------------------------------------------------------------------
        ! Basic properties of the dielectric space from input
        !
        INTEGER :: nregions
        TYPE(environ_functions), ALLOCATABLE :: regions(:)
        !
        REAL(DP) :: constant
        TYPE(environ_density) :: background
        TYPE(environ_gradient) :: gradbackground
        TYPE(environ_density) :: laplbackground
        !
        !--------------------------------------------------------------------------------
        !
        TYPE(environ_boundary), POINTER :: boundary
        ! boundary is the pointer to the object controlling the interface
        ! between the QM and the continuum region
        !
        !--------------------------------------------------------------------------------
        ! The dielectric function over space is built from the boundary of the
        ! continuum environment and the basic dielectric properties of space
        !
        TYPE(environ_density) :: epsilon
        !
        TYPE(environ_density) :: depsilon
        ! this is needed in the extra term of kohn-sham/forces
        !
        !--------------------------------------------------------------------------------
        ! Quantities related to the dielectric permittivity and
        ! they may be needed by the different solvers
        !
        LOGICAL :: need_gradient = .FALSE.
        TYPE(environ_gradient) :: gradient
        !
        LOGICAL :: need_factsqrt = .FALSE.
        TYPE(environ_density) :: factsqrt
        !
        LOGICAL :: need_gradlog = .FALSE.
        TYPE(environ_gradient) :: gradlog
        !
        !--------------------------------------------------------------------------------
        ! Dielectric polarization charges and individual components
        !
        TYPE(environ_density) :: density
        LOGICAL :: need_auxiliary = .FALSE.
        TYPE(environ_density) :: iterative
        REAL(DP) :: charge = 0.0_DP
        !
        !--------------------------------------------------------------------------------
    END TYPE environ_dielectric
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    TYPE environ_ioncctype
        !--------------------------------------------------------------------------------
        !
        INTEGER :: index
        REAL(DP) :: cbulk ! bulk concentration
        REAL(DP) :: z ! charge
        !
        TYPE(environ_density) :: c ! local concentration
        TYPE(environ_density) :: cfactor ! exp(-z\phi\beta) or 1 - z\phi\beta
        TYPE(environ_density) :: potential
        !
        !--------------------------------------------------------------------------------
    END TYPE environ_ioncctype
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    TYPE environ_electrolyte
        !--------------------------------------------------------------------------------
        !
        LOGICAL :: update = .FALSE.
        LOGICAL :: initialized = .FALSE.
        CHARACTER(LEN=80) :: electrolyte_entropy
        CHARACTER(LEN=80) :: ion_adsorption
        LOGICAL :: linearized = .FALSE.
        INTEGER :: ntyp
        TYPE(environ_ioncctype), ALLOCATABLE :: ioncctype(:)
        !
        REAL(DP) :: temperature
        REAL(DP) :: k2
        REAL(DP) :: cionmax
        REAL(DP) :: permittivity
        REAL(DP) :: adsorption_energy
        !
        TYPE(environ_boundary) :: boundary
        TYPE(environ_density) :: density
        !
        !--------------------------------------------------------------------------------
        ! The electrolyte switch function and related quantities
        !
        TYPE(environ_density) :: gamma
        TYPE(environ_density) :: dgamma
        !
        TYPE(environ_functions) :: function_
        !
        TYPE(environ_density) :: de_dboundary_second_order
        REAL(DP) :: energy_second_order
        REAL(DP) :: charge = 0.0_DP
        !
        !--------------------------------------------------------------------------------
    END TYPE environ_electrolyte
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    TYPE environ_semiconductor
        !--------------------------------------------------------------------------------
        !
        LOGICAL :: update = .FALSE.
        LOGICAL :: initialized = .FALSE.
        !
        REAL(DP) :: temperature
        REAL(DP) :: permittivity
        REAL(DP) :: carrier_density
        REAL(DP) :: electrode_charge
        REAL(DP) :: charge_threshold
        REAL(DP) :: slab_charge = 0.D0
        !
        TYPE(environ_functions) :: simple
        TYPE(environ_density) :: density
        !
        REAL(DP) :: charge = 0.0_DP
        REAL(DP) :: flatband_fermi = 0.D0
        REAL(DP) :: bulk_sc_fermi = 0.D0
        REAL(DP) :: surf_area_per_sq_cm = 0.D0
        !
        !--------------------------------------------------------------------------------
    END TYPE environ_semiconductor
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE physical_types
!----------------------------------------------------------------------------------------
