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
! This module contains all variables in the environ.in input file
! together with the routines performing initialization and broadcast
!
! Authors: Oliviero Andreussi (Department of Physics, UNT)
!          Francesco Nattino  (THEOS and NCCR-MARVEL, EPFL)
!          Ismaila Dabo       (DMSE, Penn State)
!          Nicola Marzari     (THEOS and NCCR-MARVEL, EPFL)
!
!----------------------------------------------------------------------------
MODULE environ_input
!----------------------------------------------------------------------------
  !
  USE modules_constants,  ONLY : DP, bohr_radius_angs, nsx
  USE modules_parser,     ONLY : env_field_count, env_read_line, env_get_field, parse_unit
  USE mp,                 ONLY : mp_bcast
  !
  USE environ_output, ONLY : ionode, ionode_id, comm, program_unit, &
       & verbose_ => verbose, environ_unit
  !
  IMPLICIT NONE
  !
  SAVE
  !
!=----------------------------------------------------------------------------=!
!  ENVIRON Cards Parameters
!=----------------------------------------------------------------------------=!
!
! Local parameters of external charges
!
  LOGICAL :: taextchg = .false.
  INTEGER, ALLOCATABLE  :: extcharge_dim(:)
  INTEGER, ALLOCATABLE  :: extcharge_axis(:)
  REAL(DP), ALLOCATABLE :: extcharge_charge(:)
  REAL(DP), ALLOCATABLE :: extcharge_spread(:)
  REAL(DP), ALLOCATABLE :: extcharge_pos(:,:)
  CHARACTER(len=80) :: external_charges = 'bohr'
    ! atomic_positions = 'bohr' | 'angstrom'
    ! select the units for the position of external charges being read from stdin
!
! Local parameters of dielectric regions
!
  LOGICAL :: taepsreg = .false.
  INTEGER, ALLOCATABLE  :: epsregion_dim(:)
  INTEGER, ALLOCATABLE  :: epsregion_axis(:)
  REAL(DP), ALLOCATABLE :: epsregion_eps(:,:)
  REAL(DP), ALLOCATABLE :: epsregion_width(:)
  REAL(DP), ALLOCATABLE :: epsregion_spread(:)
  REAL(DP), ALLOCATABLE :: epsregion_pos(:,:)
  CHARACTER(len=80) :: dielectric_regions = 'bohr'
    ! atomic_positions = 'bohr' | 'angstrom'
    ! select the units for the position of dielectric regions being read from stdin
!
!=----------------------------------------------------------------------------=!
!  ENVIRON Namelist Input Parameters
!=----------------------------------------------------------------------------=!
!
! Global parameters
!
        LOGICAL  :: oldenviron = .FALSE.
        ! use legacy code to compare exact numbers with Environ_0.2
        LOGICAL  :: environ_restart = .FALSE.
        ! restart a previous calculation: environ contributions are computed during
        ! initialization
        INTEGER  :: verbose = 0
        ! verbosity  0: only prints summary of polarization charge calculation;
        !    1: prints an extra file with details of iterative convergence;
        !    2: prints 3D cube files of physical properties
        REAL(DP) :: environ_thr = 1.d-1
        ! how early in scf should the corrective pot start being calculated
        INTEGER  :: environ_nskip = 1
        ! how many steps should environ skip before starting to compute the
        ! additional potentials
!
! Predefined environ types
!
        CHARACTER( LEN = 80 ) :: environ_type = 'input'
        CHARACTER( LEN = 80 ) :: environ_type_allowed(5)
        DATA environ_type_allowed / 'vacuum', 'water', 'water-cation', &
          'water-anion', 'input' /
        ! keyword to set up all the environment parameters at once to a specific set
        ! vacuum = all the flags are off (perm=1.d0, surf=0.0, pres=0.0)
        ! water = parameters optimized for water solutions in Andreussi et al.
        !         J. Chem. Phys. 136, 064102 (perm=78, surf=50, pres=-0.35)
        ! water-cation = parameters optimized for aqueous solvation of cations
        !         Dupont et al. J. Chem. Phys. 139, 214110 (perm=78, surf=, pres=)
        ! water-anion = parameters optimized for aqueous solvation of anions
        !         Dupont et al. J. Chem. Phys. 139, 214110 (perm=78, surf=, pres=)
        ! input = do not use any predefined set, use parameters from input
!
! System specification
!
        INTEGER :: system_ntyp = 0
        ! specify the atom types that are used to determine the origin and
        ! size of the system (types up to system_ntyp are used, all atoms are
        ! used by default or if system_ntyp == 0)
        INTEGER :: system_dim = 0
        ! dimensionality of the system, used to determine size (only ortogonally to
        ! periodic dimensions) and position (0 = 0D, 1 = 1D, 2 = 2D)
        INTEGER :: system_axis = 3
        ! main axis of 1D or 2D systems (1 = x, 2 = y, 3 = z)
!
! Generic keyword to specify modification of electrostatic embedding (e.g. PBC correction)
!
        LOGICAL :: env_electrostatic = .false.
        ! generic keyword that flags the need to read the electrostatic namelist
        REAL(DP) :: atomicspread(nsx) = -0.5D0
        ! gaussian spreads of the atomic density of charge, in internal units (a.u.)
        LOGICAL :: add_jellium = .false.
        ! depending on periodic boundary corrections, one may need to explicitly
        ! polarize the compensatinig jellium background
!
! Dielectric solvent parameters
!
        REAL(DP) :: env_static_permittivity = 1.D0
        ! static dielectric permittivity of the solvation model. If set equal
        ! to one (=vacuum) no dielectric effects
        REAL(DP) :: env_optical_permittivity = 1.D0
        ! optical dielectric permittivity of the solvation model. If set equal
        ! to one (=vacuum) no dielectric effects. Needed only for the TDDFTPT.
!
! Cavitation energy parameters
!
        REAL(DP) :: env_surface_tension = 0.D0
        ! solvent surface tension, if equal to zero no cavitation term
!
! PV energy parameters
!
        REAL(DP) :: env_pressure = 0.D0
        ! external pressure for PV energy, if equal to zero no pressure term
!
! Confine energy parameters
!
        REAL(DP) :: env_confine = 0.D0
        ! confinement potential
!
! Ionic countercharge parameters
!
        LOGICAL :: electrolyte_linearized = .false.
        ! solve linear-regime poisson-boltzmann problem
        INTEGER :: env_electrolyte_ntyp = 0
        ! number of counter-charge species in the electrolyte ( if != 0 must be >= 2 )
        CHARACTER( LEN = 80 ) :: electrolyte_entropy = 'full'
        CHARACTER( LEN = 80 ) :: electrolyte_entropy_allowed(2)
        DATA electrolyte_entropy_allowed / 'ions', 'full' /
        ! keyword to set the electrolyte entropy terms that are affected by the
        ! Stern-layer correction.
        ! ions = only ionic terms ( Ringe et al. J. Chem. Theory Comput. 12, 4052 )
        ! full = all terms ( Dabo et al. arXiv 0901.0096 )
        CHARACTER( LEN = 80 ) :: ion_adsorption = 'none'
        CHARACTER( LEN = 80 ) :: ion_adsorption_allowed(4)
        DATA ion_adsorption_allowed / 'none', 'anion', 'cation', 'repulsion' /
        ! include asymmetric adsorption of electrolyte.
        ! ( Baskin and Prendergast J. Electrochem. Soc. 164, E3438 )
        REAL(DP) :: cion(nsx) = 1.D0
        ! molar concentration of ionic countercharge (M=mol/L)
        REAL(DP) :: cionmax = 1.D3
        ! maximum molar concentration of ionic countercharge (M=mol/L)
        REAL(DP) :: rion = 0.D0
        ! mean atomic radius of ionic countercharge (a.u.)
        REAL(DP) :: zion(nsx) = 1.D0
        ! valence of ionic countercharge
        REAL(DP) :: temperature = 300.D0
        ! temperature of the solution
        REAL(DP) :: ion_adsorption_energy = 0.D0
        ! adsorption energy of electrolyte (Ry)
!
! External charges parameters, the remaining parameters are read from
! card EXTERNAL_CHARGES
!
        INTEGER :: env_external_charges = 0
        ! number of fixed external gaussian points/lines/planes of charges to be used
        ! in the calculation
!
! Dielectric regions parameters, the remaining parameters are read from
! card DIELECTRIC_REGIONS
!
        INTEGER :: env_dielectric_regions = 0
        ! number of fixed dielectric regions in the calculation
!
        NAMELIST /environ/                                             &
             oldenviron, environ_restart, verbose, environ_thr,        &
             environ_nskip, environ_type,                              &
             system_ntyp, system_dim, system_axis,                     &
             env_electrostatic, atomicspread, add_jellium,             &
             env_static_permittivity, env_optical_permittivity,        &
             env_surface_tension,                                      &
             env_pressure,                                             &
             env_confine,                                              &
             env_electrolyte_ntyp, cion, cionmax, rion, zion,          &
             temperature, electrolyte_linearized, electrolyte_entropy, &
             ion_adsorption, ion_adsorption_energy,                    &
             env_external_charges, env_dielectric_regions
!
!=----------------------------------------------------------------------------=!
!  BOUNDARY Namelist Input Parameters
!=----------------------------------------------------------------------------=!
!
! Global parameters
!
!
! Soft boundary (electronic) parameters
!
        INTEGER :: stype = 2
        ! type of switching functions used in the solvation models
        !    0: original Fattebert-Gygi
        !    1: ultrasoft switching function (only exponential part used for non-electrostatic)
        !    2: ultrasoft switching function as defined in Andreussi et al. JCP 2012
!
! Rigid boundary (ionic) parameters
!
        CHARACTER( LEN = 80 ) :: radius_mode = 'uff'
        CHARACTER( LEN = 80 ) :: radius_mode_allowed(3)
        DATA radius_mode_allowed / 'pauling', 'bondi', 'uff' /
        ! type of hardcoded solvation radii to be used when solvent_mode = 'ionic'
        ! pauling = R.C. Weast, ed., Handbook of chemistry and physics (CRC Press, Cleveland, 1981)
        ! bondi   = A. Bondi, J. Phys. Chem. 68, 441 (1964)
        ! uff     = A.K. Rapp/'{e} et al. J. Am. Chem. Soc. 114(25) pp.10024-10035 (1992)
        REAL(DP) :: solvationrad(nsx) = -3.D0
        ! solvationrad radius of the solvation shell for each species when the
        ! ionic dielectric function is adopted, in internal units (a.u.)
!
! Full boundary parameters
!
        REAL(DP) :: corespread(nsx) = -0.5D0
        ! gaussian spreads of the core electrons, in internal units (a.u.), to
        ! be used when solvent_mode = 'full'
!
! Solvent-aware boundary parameters
!
        REAL(DP) :: solvent_radius = 0.D0
        ! size of the solvent, used to decide whether to fill a continuum
        ! void or not. If set equal to 0.D0, use the standard algorithm
        REAL(DP) :: radial_scale = 2.D0
        ! compute the filled fraction on a spherical volume scaled wrt solvent
        ! size
        REAL(DP) :: radial_spread = 0.5D0
        ! spread of the step function used to evaluate occupied volume
        REAL(DP) :: filling_threshold = 0.825D0
        ! threshold to decide whether to fill a continuum void or not, to be
        ! compared with the filled fraction: if filled fraction .GT. threshold
        ! THEN fill gridpoint
        REAL(DP) :: filling_spread = 0.02D0
        ! spread of the switching function used to decide whether the continuum
        ! void should be filled or not
!
! Numerical core's parameters
!
        CHARACTER( LEN = 80 ) :: boundary_core = 'analytic'
        CHARACTER( LEN = 80 ) :: boundary_core_allowed(4)
        DATA boundary_core_allowed / 'fft', 'fd', 'analytic', 'highmem' /
        ! choice of the core numerical methods to be exploited for the quantities derived from the dielectric
        ! fft       = fast Fourier transforms
        ! fd        = finite differences in real space
        ! analytic  = analytic derivatives for as much as possible (and FFTs for the rest)
        ! highmem   = analytic derivatives for soft-sphere computed by storing all spherical functions and derivatives
!
! Finite differences' parameters
!
        INTEGER  :: ifdtype = 1
        ! type of numerical differentiator: 1=central differences,
        ! 2=low-noise lanczos (m=2), 3=low-noise lanczos (m=4),
        ! 4=smooth noise-robust (n=2), 5=smooth noise-robust (n=4)
        INTEGER  :: nfdpoint = 2
        ! number of points used in the numerical differentiator
        ! N = 2*nfdpoint+1
!
! Solvent boundary parameters
!
        CHARACTER( LEN = 80 ) :: solvent_mode = 'electronic'
        CHARACTER( LEN = 80 ) :: solvent_mode_allowed(8)
        DATA solvent_mode_allowed / 'electronic', 'ionic', 'full', 'external', &
                              & 'system', 'elec-sys', 'ionic-sys', 'full-sys' /
        ! solvent_mode method for calculating the density that sets
        ! the dielectric constant
        ! electronic = dielectric depends self-consist. on electronic density
        ! ionic = dielectric defined on a fictitious ionic density, generated
        !         as the sum of spherical error functions centered on atomic
        !         positions of width specified in input by solvationrad(ityp)
        ! full  = similar to electronic, but an extra density is added to
        !         represent the core electrons and the nuclei. This extra
        !         density is defined as the sum of gaussian functions centered
        !         on atomic positions of width equal to corespread(ityp)
        ! system = simplified regular dielectric defined to be outside a distance
        !         solvent_distance from the specified system
        ! elec-sys = similar to electronic, but on top of the system dielectric
        ! ionic-sys = similar to ionic, but on top of the system dielectric
        ! full-sys = similar to full, but on top of the system dielectric
!
! Soft solvent boundary (electronic) parameters
!
        REAL(DP) :: rhomax = 0.005
        ! first parameter of the sw function, roughly corresponding
        ! to the density threshold of the solvation model
        REAL(DP) :: rhomin = 0.0001
        ! second parameter of the sw function when stype=1 or 2
        REAL(DP) :: tbeta = 4.8
        ! second parameter of the sw function when stype=0
!
! Rigid solvent boundary (ionic) parameters
!
        REAL(DP) :: alpha = 1.D0
        ! scaling factor for ionic radii when solvent_mode = 'ionic'
        REAL(DP) :: softness = 0.5D0
        ! spread of the rigid interfaces
!
! Simplified solvent boundary (system) parameters
!
        REAL(DP) :: solvent_distance = 1.D0
        ! distance from the system where the boundary starts if required from solvent_mode
        REAL(DP) :: solvent_spread = 0.5D0
        ! spread of the boundary interface if defined on system position and width
!
! Stern boundary parameters
!
        CHARACTER( LEN = 80 ) :: electrolyte_mode = 'electronic'
        CHARACTER( LEN = 80 ) :: electrolyte_mode_allowed(8)
        DATA electrolyte_mode_allowed / 'electronic', 'ionic', 'full', 'external', &
                                & 'system', 'elec-sys', 'ionic-sys', 'full-sys' /
        ! electrolyte_mode method for calculating the density that sets
        ! the onset of ionic countercharge ( see solvent_mode above )
!
! Soft Stern boundary (electronic) parameters
!
        REAL(DP) :: electrolyte_rhomax = 0.005D0
        ! first parameter of the Stern sw function, roughly corresponding
        ! to the density threshold of the ionic countercharge.
        REAL(DP) :: electrolyte_rhomin = 0.0001D0
        ! second parameter of the Stern sw function when stype=1 or 2
        REAL(DP) :: electrolyte_tbeta = 4.8D0
        ! second parameter of the Stern sw function when stype=0
!
! Rigid Stern boundary (ionic) parameters
!
        REAL(DP) :: electrolyte_alpha = 1.D0
        ! scaling factor for ionic radii when electrolyte_mode = 'ionic'
        REAL(DP) :: electrolyte_softness = 0.5D0
        ! spread of the rigid Stern interfaces
!
! Simplified Stern boundary (system) parameters
!
        REAL(DP) :: electrolyte_distance = 0.D0
        ! distance from the system where the electrolyte boundary starts
        REAL(DP) :: electrolyte_spread = 0.5D0
        ! spread of the interfaces for the electrolyte boundary
!
        NAMELIST /boundary/                      &
             solvent_mode,                       &
             radius_mode, alpha, softness,       &
             solvationrad,                       &
             stype, rhomax, rhomin, tbeta,       &
             corespread,                         &
             solvent_distance, solvent_spread,   &
             solvent_radius, radial_scale,       &
             radial_spread, filling_threshold,   &
             filling_spread,                     &
             electrolyte_mode, electrolyte_distance,         &
             electrolyte_spread, electrolyte_rhomax,         &
             electrolyte_rhomin, electrolyte_tbeta,          &
             electrolyte_alpha, electrolyte_softness,        &
             boundary_core,                      &
             ifdtype, nfdpoint
!
!=----------------------------------------------------------------------------=!
!  ELECTROSTATIC Namelist Input Parameters
!=----------------------------------------------------------------------------=!
!
! Global parameters
!
        CHARACTER( LEN = 80 ) :: problem = 'none'
        CHARACTER( LEN = 80 ) :: problem_allowed(6)
        DATA problem_allowed / 'poisson', 'generalized', 'pb', 'modpb', 'linpb', 'linmodpb' /
        ! type of electrostatic problem:
        ! poisson     = standard poisson equation, with or without boundary conditions (default)
        ! generalized = generalized poisson equation
        ! pb          = poisson-boltzmann equation (non-linear)
        ! modpb       = modified poisson-boltzmann equation (non-linear)
        ! linpb       = linearized poisson-boltzmann equation (debye-huckel)
        ! linmodpb    = linearized modified poisson-boltzmann equation
        REAL(DP) :: tol = 1.D-5
        ! convergence threshold for electrostatic potential or auxiliary charge
        REAL(DP) :: inner_tol = 1.D-5
        ! same as tol for inner loop in nested algorithms
!
! Driver's parameters
!
        CHARACTER( LEN = 80 ) :: solver = 'none'
        CHARACTER( LEN = 80 ) :: solver_allowed(7)
        DATA solver_allowed / 'cg', 'sd', 'iterative', 'lbfgs', 'newton', 'nested', 'direct' /
        ! type of numerical solver
        ! direct    = for simple problems with analytic or direct solution
        ! cg        = conjugate gradient (default)
        ! sd        = steepest descent
        ! iterative = fixed-point search
        ! lbfgs     = low-memory bfgs
        ! newton    = newton's method (only for non-linear problem)
        ! nested    = double iterations (only for non-linear problem)
        CHARACTER( LEN = 80 ) :: auxiliary = 'none'
        CHARACTER( LEN = 80 ) :: auxiliary_allowed(4)
        DATA auxiliary_allowed / 'none', 'full', 'pol', 'ioncc' /
        ! solve with respect to the potential or with respect to an auxiliary charge density
        ! none  = solve for the potential (default)
        ! full  = solve for the auxiliary charge density
        ! pol   = in a nested scheme, solve the inner (pol) cycle in terms of the auxiliary charge
        ! ioncc = in a nested scheme, solve the outer (ioncc) cycle in terms of the auxiliary charge
        CHARACTER( LEN = 80 ) :: step_type = 'optimal'
        CHARACTER( LEN = 80 ) :: step_type_allowed(3)
        DATA step_type_allowed / 'optimal', 'input', 'random' /
        ! how to choose the step size in gradient descent algorithms or iterative mixing
        ! optimal = step size that minimize the cost function on the descent direction
        ! input   = fixed step size as defined in input (step keyword)
        ! random  = random step size within zero and twice the optima value
        REAL(DP) :: step = 0.3
        ! step size to be used if step_type = 'input' (inherits the tasks of the old mixrhopol)
        INTEGER :: maxstep = 200
        ! maximum number of steps to be performed by gradient or iterative solvers
        CHARACTER( LEN = 80 ) :: inner_solver = 'none'
        CHARACTER( LEN = 80 ) :: inner_solver_allowed(5)
        DATA inner_solver_allowed / 'none', 'cg', 'sd', 'iterative', 'direct' /
        ! type of numerical solver for inner loop in nested algorithms
        INTEGER :: inner_maxstep = 200
        ! same as maxstep for inner loop in nested algorithms
!
! Iterative driver's parameters (OBSOLETE)
!
        CHARACTER( LEN = 80 ) :: mix_type = 'linear'
        CHARACTER( LEN = 80 ) :: mix_type_allowed(4)
        DATA mix_type_allowed / 'linear', 'anderson', 'diis', 'broyden' /
        ! mixing method for iterative calculations
        ! 'linear', 'anderson', 'diis', 'broyden'
        INTEGER :: ndiis=1
        ! order of DIIS interpolation of iterative calculation
        REAL(DP) :: mix = 0.5
        ! mixing parameter to be used in the iterative driver
        REAL(DP) :: inner_mix = 0.5
        ! same as mix but for inner loop in nested algorithm
!
! Preconditioner's parameters
!
        CHARACTER( LEN = 80 ) :: preconditioner = 'sqrt'
        CHARACTER( LEN = 80 ) :: preconditioner_allowed(3)
        DATA preconditioner_allowed / 'none', 'sqrt', 'left' /
        ! type of preconditioner
        ! none      = no preconditioner
        ! left      = left linear preconditioner eps nabla v = r
        ! sqrt      = sqrt preconditioner sqrt(eps) nabla ( sqrt(eps) * v ) = r
        CHARACTER( LEN = 80 ) :: screening_type = 'none'
        CHARACTER( LEN = 80 ) :: screening_type_allowed(4)
        DATA screening_type_allowed / 'none', 'input', 'linear', 'optimal' /
        ! use the screened coulomb Green's function instead of the vacuum one
        ! none      = unscreened coulomb
        ! input     = screened coulomb with screening lenght provided in input
        ! linear    = screened coulomb with screening lenght from linear component of the problem
        ! optimal   = screened coulomb with screening lenght optimized (to be defined)
        REAL(DP) :: screening = 0.D0
        ! screening lenght to be used if screening_type = 'input'
!
! Numerical core's parameters
!
        CHARACTER( LEN = 80 ) :: core = 'fft'
        CHARACTER( LEN = 80 ) :: core_allowed(1)
        DATA core_allowed / 'fft' /
        ! choice of the core numerical methods to be exploited for the different operations
        ! fft = fast Fourier transforms (default)
        ! to be implemented : wavelets (from big-DFT) and multigrid
!
! Periodic correction keywords
!
        INTEGER :: pbc_dim = -3
        ! dimensionality of the simulation cell
        ! periodic boundary conditions on 3/2/1/0 sides of the cell
        CHARACTER( LEN = 80 ) :: pbc_correction = 'none'
        CHARACTER( LEN = 80 ) :: pbc_correction_allowed(3)
        DATA pbc_correction_allowed / 'none', 'parabolic', 'gcs' /
        ! type of periodic boundary condition correction to be used
        ! parabolic = point-counter-charge type of correction
        INTEGER :: pbc_axis = 3
        ! choice of the sides with periodic boundary conditions
        ! 1 = x, 2 = y, 3 = z, where
        ! if pbc_dim = 2, cell_axis is orthogonal to 2D plane
        ! if pbc_dim = 1, cell_axis is along the 1D direction
!
        NAMELIST /electrostatic/                 &
             problem, tol, solver, auxiliary,    &
             step_type, step, maxstep,           &
             mix_type, mix, ndiis,               &
             preconditioner,                     &
             screening_type, screening,          &
             core,                               &
             pbc_dim, pbc_correction, pbc_axis,  &
             inner_tol, inner_solver,            &
             inner_maxstep, inner_mix
!
CONTAINS
!--------------------------------------------------------------------
! BACKWARD COMPATIBILITY
! Compatible with QE-6.0 QE-6.1.X QE-6.2.X QE-6.3.X
!  SUBROUTINE read_environ(prog,nelec,nspin,nat,ntyp,atom_label,use_internal_pbc_corr,ion_radius)
! Compatible with QE-6.4.X QE-GIT
  SUBROUTINE read_environ(prog,nelec,nat,ntyp,atom_label,use_internal_pbc_corr,ion_radius)
! END BACKWARD COMPATIBILITY
!--------------------------------------------------------------------
    !
    USE environ_init, ONLY : set_environ_base
    USE electrostatic_init, ONLY : set_electrostatic_base
    !
    CHARACTER(len=*), INTENT(IN) :: prog
    LOGICAL, INTENT(IN) :: use_internal_pbc_corr
    INTEGER, INTENT(IN) :: nelec, nat, ntyp
! BACKWARD COMPATIBILITY
! Compatible with QE-6.0 QE-6.1.X QE-6.2.X QE-6.3.X
!    INTEGER, INTENT(IN) :: nspin
! Compatible with QE-6.4.X QE-GIT
!
! END BACKWARD COMPATIBILITY
    CHARACTER(len=3), DIMENSION(:), INTENT(IN) :: atom_label
    REAL( DP ), DIMENSION(:), INTENT(IN), OPTIONAL :: ion_radius
    !
    INTEGER, EXTERNAL :: find_free_unit
    !
    LOGICAL :: ext
    INTEGER :: environ_unit_input
    INTEGER :: is
    !
    ! ... Open environ input file: environ.in
    !
    environ_unit_input = find_free_unit()
    INQUIRE(file="environ.in",exist=ext)
    IF (.NOT.ext) CALL errore( 'read_environ',&
         & ' missing environ.in file ', 1 )
    OPEN(unit=environ_unit_input,file="environ.in",status="old")
    !
    ! ... Read environ namelists
    !
    CALL environ_read_namelist( environ_unit_input )
    !
    ! ... Read environ cards
    !
    CALL environ_read_cards( environ_unit_input )
    !
    ! ... Close environ input file
    !
    CLOSE( environ_unit_input )
    !
    ! ... If passed from input, overwrites atomic spread
    ! (USED IN CP TO HAVE CONSISTENT RADII FOR ELECTROSTATICS)
    !
    IF ( PRESENT(ion_radius) ) THEN
       !
       DO is = 1, ntyp
          atomicspread(is) = ion_radius(is)
       ENDDO
       !
    ENDIF
    !
    ! ... Set verbosity and open debug file
    !
    verbose_ = verbose
    !
    IF ( verbose_ .GE. 1 ) &
         OPEN(unit=environ_unit,file='environ.debug',status='unknown')
    !
    ! ... Set module variables according to input
    !
    !
    ! ... Set electrostatic first as it does not depend on anything else
    !
    CALL set_electrostatic_base ( problem, tol, solver, auxiliary,       &
                                  step_type, step, maxstep, mix_type,    &
                                  ndiis, mix, preconditioner,            &
                                  screening_type, screening, core,       &
                                  boundary_core, ifdtype, nfdpoint,      &
                                  use_internal_pbc_corr, pbc_correction, &
! BACKWARD COMPATIBILITY
! Compatible with QE-6.0 QE-6.1.X QE-6.2.X QE-6.3.X
!                                  pbc_dim, pbc_axis, nspin, prog,        &
! Compatible with QE-6.4.X QE-GIT
                                  pbc_dim, pbc_axis, prog,        &
! END BACKWARD COMPATIBILITY
                                  inner_tol, inner_solver, inner_maxstep,&
                                  inner_mix )
    !
    ! ... Then set environ base
    !
! BACKWARD COMPATIBILITY
! Compatible with QE-6.0 QE-6.1.X QE-6.2.X QE-6.3.X
!    CALL set_environ_base  ( prog, nelec, nspin,                         &
! Compatible with QE-6.4.X QE-GIT
    CALL set_environ_base  ( prog, nelec,                                &
! END BACKWARD COMPATIBILITY
                             nat, ntyp, atom_label, atomicspread,        &
                             corespread, solvationrad,                   &
                             oldenviron, environ_restart, environ_thr,   &
                             environ_nskip, environ_type,                &
                             system_ntyp, system_dim, system_axis,       &
                             stype, rhomax, rhomin, tbeta,               &
                             env_static_permittivity,                    &
                             env_optical_permittivity,                   &
                             solvent_mode,                               &
                             radius_mode, alpha, softness,               &
                             solvent_distance, solvent_spread,           &
                             solvent_radius, radial_scale,               &
                             radial_spread, filling_threshold,           &
                             filling_spread,                             &
                             add_jellium,                                &
                             env_surface_tension,                        &
                             env_pressure,                               &
                             env_confine,                                &
                             env_electrolyte_ntyp,                       &
                             electrolyte_linearized, electrolyte_entropy,&
                             electrolyte_mode, electrolyte_distance,     &
                             electrolyte_spread, cion, cionmax, rion,    &
                             zion, electrolyte_rhomax,                   &
                             electrolyte_rhomin, electrolyte_tbeta,      &
                             electrolyte_alpha, electrolyte_softness,    &
                             ion_adsorption, ion_adsorption_energy,      &
                             temperature,                                &
                             env_external_charges,                       &
                             extcharge_charge, extcharge_dim,            &
                             extcharge_axis, extcharge_pos,              &
                             extcharge_spread,                           &
                             env_dielectric_regions,                     &
                             epsregion_eps, epsregion_dim,               &
                             epsregion_axis, epsregion_pos,              &
                             epsregion_spread, epsregion_width )
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE read_environ
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE environ_read_namelist( environ_unit_input )
!--------------------------------------------------------------------
    !
    !  Environ namelist parsing routine
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: environ_unit_input
    !
    LOGICAL :: lboundary, lelectrostatic
    INTEGER :: ios
    !
    ! ... Set the defauls
    !
    CALL environ_defaults()
    !
    CALL boundary_defaults()
    !
    CALL electrostatic_defaults()
    !
    ! ... Read the &ENVIRON namelist
    !
    ios = 0
    IF( ionode ) READ( environ_unit_input, environ, iostat = ios )
    CALL mp_bcast( ios, ionode_id, comm )
    IF( ios /= 0 ) CALL errore( ' read_environ ', &
         & ' reading namelist environ ', ABS(ios) )
    !
    ! ... Broadcast &ENVIRON variables
    !
    CALL environ_bcast()
    !
    ! ... Check &ENVIRON variables
    !
    CALL environ_checkin()
    !
    ! ... Fix some &BOUNDARY defaults depending on &ENVIRON
    !
    CALL fix_boundary(lboundary)
    !
    ! ... Read the &BOUNDARY namelist only if needed
    !
    ios = 0
    IF( ionode ) THEN
       IF ( lboundary ) READ( environ_unit_input, boundary, iostat = ios )
    END IF
    CALL mp_bcast( ios, ionode_id, comm )
    IF( ios /= 0 ) CALL errore( ' read_environ ', &
         & ' reading namelist boundary ', ABS(ios) )       !
    !
    ! ... Broadcast &BOUNDARY variables
    !
    CALL boundary_bcast()
    !
    ! ... Check &BOUNDARY variables
    !
    CALL boundary_checkin()
    !
    ! ... Set predefined envinron_types, also according to the boundary
    !
    CALL set_environ_type()
    !
    ! ... Fix some &ELECTROSTATIC defaults depending on &ENVIRON and &BOUNDARY
    !
    CALL fix_electrostatic(lelectrostatic)
    !
    ! ... Read the &ELECTROSTATIC namelist only if needed
    !
    ios = 0
    IF( ionode ) THEN
       IF ( lelectrostatic ) READ( environ_unit_input, electrostatic, iostat = ios )
    END IF
    CALL mp_bcast( ios, ionode_id, comm )
    IF( ios /= 0 ) CALL errore( ' read_environ ', &
         & ' reading namelist electrostatic ', ABS(ios) )
    !
    ! ... Broadcast &ELECTROSTATIC variables
    !
    CALL electrostatic_bcast()
    !
    ! ... Set electrostatic problem
    !
    CALL set_electrostatic_problem( )
    !
    ! ... Check &ELECTROSTATIC variables
    !
    CALL electrostatic_checkin()
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE environ_read_namelist
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE environ_defaults( )
!--------------------------------------------------------------------
    !
    !  Variables initialization for Namelist ENVIRON
    !
    IMPLICIT NONE
    !
    oldenviron = .false.
    environ_restart = .false.
    verbose       = 0
    environ_thr   = 1.D-1
    environ_nskip = 1
    environ_type  = 'input'
    !
    system_ntyp = 0
    system_dim = 0
    system_axis = 3
    !
    env_electrostatic = .false.
    atomicspread(:) = -0.5D0
    add_jellium = .false.
    !
    env_static_permittivity = 1.D0
    env_optical_permittivity = 1.D0
    !
    env_surface_tension = 0.D0
    !
    env_pressure = 0.D0
    !
    env_confine = 0.D0
    !
    env_electrolyte_ntyp = 0
    electrolyte_linearized = .false.
    electrolyte_entropy = 'full'
    cion(:) = 1.0D0
    cionmax = 0.0D0 ! if remains zero, pb or linpb
    rion = 0.D0
    zion(:) = 0.D0
    temperature = 300.0D0
    !
    ion_adsorption = 'none'
    ion_adsorption_energy = 0.D0
    !
    env_external_charges = 0
    env_dielectric_regions = 0
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE environ_defaults
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE boundary_defaults( )
!--------------------------------------------------------------------
    !
    !  Variables initialization for Namelist BOUNDARY
    !
    IMPLICIT NONE
    !
    solvent_mode = 'electronic'
    !
    radius_mode     = 'uff'
    alpha           = 1.D0
    softness        = 0.5D0
    solvationrad(:) = -3.D0
    !
    stype   = 2
    rhomax  = 0.005
    rhomin  = 0.0001
    tbeta   = 4.8
    !
    corespread(:)   = -0.5D0
    !
    solvent_distance = 1.D0
    solvent_spread   = 0.5D0
    !
    solvent_radius     = 0.D0
    radial_scale       = 2.D0
    radial_spread      = 0.5D0
    filling_threshold  = 0.825D0
    filling_spread     = 0.02D0
    !
    electrolyte_mode = 'electronic'
    !
    electrolyte_distance = 0.D0
    electrolyte_spread = 0.5D0
    !
    electrolyte_rhomax = 0.005D0
    electrolyte_rhomin = 0.0001D0
    electrolyte_tbeta = 4.8D0
    !
    electrolyte_alpha = 1.D0
    electrolyte_softness = 0.5D0
    !
    boundary_core = 'analytic'
    ifdtype  = 1
    nfdpoint = 2
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE boundary_defaults
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE electrostatic_defaults( )
!--------------------------------------------------------------------
    !
    !  Variables initialization for Namelist ELECTROSTATIC
    !
    IMPLICIT NONE
    !
    problem = 'none'
    tol = 1.D-5
    !
    solver = 'none'
    auxiliary = 'none'
    step_type = 'optimal'
    step = 0.3D0
    maxstep = 200
    inner_solver = 'none'
    inner_tol = 1.D-10
    inner_maxstep = 200
    inner_mix = 0.5D0
    !
    mix_type = 'linear'
    ndiis = 1
    mix = 0.5D0
    !
    preconditioner = 'sqrt'
    screening_type = 'none'
    screening = 0.D0
    !
    core = 'fft'
    !
    pbc_dim = -3
    pbc_correction = 'none'
    pbc_axis = 3
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE electrostatic_defaults
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE environ_bcast()
!--------------------------------------------------------------------
    !
    !  Broadcast variables values for Namelist ENVIRON
    !
    IMPLICIT NONE
    !
    CALL mp_bcast( oldenviron,                 ionode_id, comm )
    CALL mp_bcast( environ_restart,            ionode_id, comm )
    CALL mp_bcast( verbose,                    ionode_id, comm )
    CALL mp_bcast( environ_thr,                ionode_id, comm )
    CALL mp_bcast( environ_nskip,              ionode_id, comm )
    CALL mp_bcast( environ_type,               ionode_id, comm )
    !
    CALL mp_bcast( system_ntyp,                ionode_id, comm )
    CALL mp_bcast( system_dim,                 ionode_id, comm )
    CALL mp_bcast( system_axis,                ionode_id, comm )
    !
    CALL mp_bcast( env_electrostatic,          ionode_id, comm )
    CALL mp_bcast( atomicspread,               ionode_id, comm )
    CALL mp_bcast( add_jellium,                ionode_id, comm )
    !
    CALL mp_bcast( env_static_permittivity,    ionode_id, comm )
    CALL mp_bcast( env_optical_permittivity,   ionode_id, comm )
    !
    CALL mp_bcast( env_surface_tension,        ionode_id, comm )
    !
    CALL mp_bcast( env_pressure,               ionode_id, comm )
    !
    CALL mp_bcast( env_confine,                ionode_id, comm )
    !
    CALL mp_bcast( env_electrolyte_ntyp,       ionode_id, comm )
    CALL mp_bcast( electrolyte_linearized,           ionode_id, comm )
    CALL mp_bcast( electrolyte_entropy,              ionode_id, comm )
    CALL mp_bcast( cion,                       ionode_id, comm )
    CALL mp_bcast( cionmax,                    ionode_id, comm )
    CALL mp_bcast( rion,                       ionode_id, comm )
    CALL mp_bcast( zion,                       ionode_id, comm )
    CALL mp_bcast( temperature,        ionode_id, comm )
    !
    CALL mp_bcast( ion_adsorption,             ionode_id, comm )
    CALL mp_bcast( ion_adsorption_energy,      ionode_id, comm )
    !
    CALL mp_bcast( env_external_charges,       ionode_id, comm )
    CALL mp_bcast( env_dielectric_regions,     ionode_id, comm )
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE environ_bcast
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE boundary_bcast()
!--------------------------------------------------------------------
    !
    !  Broadcast variables values for Namelist BOUNDARY
    !
    IMPLICIT NONE
    !
    CALL mp_bcast( solvent_mode,               ionode_id, comm )
    !
    CALL mp_bcast( stype,                      ionode_id, comm )
    CALL mp_bcast( rhomax,                     ionode_id, comm )
    CALL mp_bcast( rhomin,                     ionode_id, comm )
    CALL mp_bcast( tbeta,                      ionode_id, comm )
    !
    CALL mp_bcast( radius_mode,                ionode_id, comm )
    CALL mp_bcast( alpha,                      ionode_id, comm )
    CALL mp_bcast( softness,                   ionode_id, comm )
    CALL mp_bcast( solvationrad,               ionode_id, comm )
    !
    CALL mp_bcast( corespread,                 ionode_id, comm )
    !
    CALL mp_bcast( solvent_distance,           ionode_id, comm )
    CALL mp_bcast( solvent_spread,             ionode_id, comm )
    !
    CALL mp_bcast( solvent_radius,             ionode_id, comm )
    CALL mp_bcast( radial_scale,               ionode_id, comm )
    CALL mp_bcast( radial_spread,              ionode_id, comm )
    CALL mp_bcast( filling_threshold,          ionode_id, comm )
    CALL mp_bcast( filling_spread,             ionode_id, comm )
    !
    CALL mp_bcast( electrolyte_mode,                 ionode_id, comm )
    !
    CALL mp_bcast( electrolyte_distance,             ionode_id, comm )
    CALL mp_bcast( electrolyte_spread,               ionode_id, comm )
    !
    CALL mp_bcast( electrolyte_rhomax,               ionode_id, comm )
    CALL mp_bcast( electrolyte_rhomin,               ionode_id, comm )
    CALL mp_bcast( electrolyte_tbeta,                ionode_id, comm )
    !
    CALL mp_bcast( electrolyte_alpha,                ionode_id, comm )
    CALL mp_bcast( electrolyte_softness,             ionode_id, comm )
    !
    CALL mp_bcast( boundary_core,              ionode_id, comm )
    CALL mp_bcast( ifdtype,                    ionode_id, comm )
    CALL mp_bcast( nfdpoint,                   ionode_id, comm )
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE boundary_bcast
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE electrostatic_bcast()
!--------------------------------------------------------------------
    !
    !  Broadcast variables values for Namelist ELECTROSTATIC
    !
    IMPLICIT NONE
    !
    CALL mp_bcast( problem,                    ionode_id, comm )
    CALL mp_bcast( tol,                        ionode_id, comm )
    !
    CALL mp_bcast( solver,                     ionode_id, comm )
    CALL mp_bcast( inner_solver,               ionode_id, comm )
    CALL mp_bcast( inner_tol,                  ionode_id, comm )
    CALL mp_bcast( inner_maxstep,              ionode_id, comm )
    CALL mp_bcast( inner_mix,                  ionode_id, comm )

    CALL mp_bcast( auxiliary,                  ionode_id, comm )
    CALL mp_bcast( step_type,                  ionode_id, comm )
    CALL mp_bcast( step,                       ionode_id, comm )
    CALL mp_bcast( maxstep,                    ionode_id, comm )
    !
    CALL mp_bcast( mix_type,                   ionode_id, comm )
    CALL mp_bcast( mix,                        ionode_id, comm )
    CALL mp_bcast( ndiis,                      ionode_id, comm )
    !
    CALL mp_bcast( preconditioner,             ionode_id, comm )
    CALL mp_bcast( screening_type,             ionode_id, comm )
    CALL mp_bcast( screening,                  ionode_id, comm )
    !
    CALL mp_bcast( core,                       ionode_id, comm )
    !
    CALL mp_bcast( pbc_dim,                    ionode_id, comm )
    CALL mp_bcast( pbc_correction,             ionode_id, comm )
    CALL mp_bcast( pbc_axis,                   ionode_id, comm )
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE electrostatic_bcast
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE environ_checkin()
!--------------------------------------------------------------------
    !
    !  Check input values for Namelist ENVIRON
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=20) :: sub_name = ' environ_checkin '
    INTEGER           :: i
    LOGICAL           :: allowed = .FALSE.
    !
    IF ( oldenviron ) &
         CALL infomsg( sub_name,' use some old legacy code for environ' )
    IF ( environ_restart ) &
         CALL infomsg( sub_name,' environ restarting' )
    IF( verbose < 0 ) &
         CALL errore( sub_name,' verbose out of range ', 1 )
    IF( environ_thr < 0.0_DP ) &
         CALL errore( sub_name,' environ_thr out of range ', 1 )
    IF( environ_nskip < 0 ) &
         CALL errore( sub_name,' environ_nskip out of range ', 1 )
    allowed = .FALSE.
    DO i = 1, SIZE( environ_type_allowed )
       IF( TRIM(environ_type) == environ_type_allowed(i) ) allowed = .TRUE.
    END DO
    IF( .NOT. allowed ) &
         CALL errore( sub_name, ' environ_type '''// &
         & TRIM(environ_type)//''' not allowed ', 1 )
    !
    IF( system_ntyp < 0 ) &
         CALL errore( sub_name,' system_ntype out of range ', 1 )
    IF( system_dim < 0 .OR. system_dim > 3 ) &
         CALL errore( sub_name,' system_dim out of range ', 1 )
    IF( system_axis < 1 .OR. system_axis > 3 ) &
         CALL errore( sub_name,' system_axis out of range ', 1 )
    !
    IF( env_static_permittivity < 1.0_DP ) &
         CALL errore( sub_name,' env_static_permittivity out of range ', 1 )
    IF( env_optical_permittivity < 1.0_DP ) &
         CALL errore( sub_name,' env_optical_permittivity out of range ', 1 )
    !
    IF( env_surface_tension < 0.0_DP ) &
         CALL errore( sub_name,' env_surface_tension out of range ', 1 )
    !
    IF( env_electrolyte_ntyp < 0 .OR. env_electrolyte_ntyp .EQ. 1 ) &
         CALL errore( sub_name,' env_electrolyte_ntyp out of range ', 1 )
    allowed = .FALSE.
    DO i = 1, SIZE( electrolyte_entropy_allowed )
       IF( TRIM(electrolyte_entropy) == electrolyte_entropy_allowed(i) ) allowed = .TRUE.
    END DO
    IF( .NOT. allowed ) &
         CALL errore( sub_name, ' electrolyte_entropy '''// &
         & TRIM(electrolyte_entropy)//''' not allowed ', 1 )
    IF( temperature < 0.0_DP ) &
         CALL errore( sub_name,' temperature out of range ', 1 )
    DO i = 1, env_electrolyte_ntyp
       IF ( cion(i) .LT. 0.D0 ) THEN
          CALL errore( sub_name, ' cion cannot be negative ', 1 )
       END IF
    END DO
    IF ( cionmax .LT. 0.D0 .OR. rion .LT. 0.D0 ) &
         CALL errore( sub_name,'cionmax and rion cannot be negative ', 1 )
    IF ( cionmax .GT. 0.D0 .AND. rion .GT. 0.D0 ) &
         CALL errore( sub_name,'either cionmax or rion can be set ', 1 )
    allowed = .FALSE.
    DO i = 1, SIZE( ion_adsorption_allowed )
       IF( TRIM(ion_adsorption) == ion_adsorption_allowed(i) ) allowed = .TRUE.
    END DO
    IF( .NOT. allowed ) &
         CALL errore( sub_name, ' ion_adsorption '''// &
         & TRIM(ion_adsorption)//''' not allowed ', 1 )
    IF ( ion_adsorption_energy .LT. 0D0 ) &
         CALL errore( sub_name,'ion_adsorption_energy must be positive', 1 )
    IF( .NOT. TRIM(ion_adsorption) .EQ. 'none' ) &
         & CALL errore( sub_name,'ion_adsorption not implemented', 1 )
    !
    IF ( env_external_charges < 0 ) &
         CALL errore( sub_name,' env_external_charges out of range ', 1 )
    !
    IF ( env_dielectric_regions < 0 ) &
         CALL errore( sub_name,' env_dielectric_regions out of range ', 1 )
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE environ_checkin
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE boundary_checkin()
!--------------------------------------------------------------------
    !
    !  Check input values for Namelist BOUNDARY
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=20) :: sub_name = ' boundary_checkin '
    INTEGER           :: i
    LOGICAL           :: allowed = .FALSE.
    !
    allowed = .FALSE.
    DO i = 1, SIZE( solvent_mode_allowed )
       IF( TRIM(solvent_mode) == solvent_mode_allowed(i) ) allowed = .TRUE.
    END DO
    IF( .NOT. allowed ) &
         CALL errore( sub_name, ' solvent_mode '''// &
         & TRIM(solvent_mode)//''' not allowed ', 1 )
    !
    IF( stype > 2 ) &
         CALL errore( sub_name,' stype out of range ', 1 )
    IF( rhomax < 0.0_DP ) &
         CALL errore( sub_name,' rhomax out of range ', 1 )
    IF( rhomin < 0.0_DP ) &
         CALL errore( sub_name,' rhomin out of range ', 1 )
    IF( rhomax < rhomin ) &
         CALL errore( sub_name,' inconsistent rhomax and rhomin', 1 )
    IF( tbeta < 0.0_DP ) &
         CALL errore( sub_name,' tbeta out of range ', 1 )
    !
    allowed = .FALSE.
    DO i = 1, SIZE( radius_mode_allowed )
       IF( TRIM(radius_mode) == radius_mode_allowed(i) ) allowed = .TRUE.
    END DO
    IF( .NOT. allowed ) &
         CALL errore( sub_name, ' radius_mode '''// &
         & TRIM(radius_mode)//''' not allowed ', 1 )
    IF( alpha <= 0.0_DP ) &
         CALL errore( sub_name,' alpha out of range ', 1 )
    IF( softness <= 0.0_DP ) &
         CALL errore( sub_name,' softness out of range ', 1 )
    !
    IF( solvent_spread <= 0.0_DP ) &
         CALL errore( sub_name,' solvent_spread out of range ', 1 )
    !
    IF ( solvent_radius < 0.0_DP ) &
         CALL errore( sub_name, 'solvent_radius out of range ', 1 )
    IF ( radial_scale < 1.0_DP ) &
         CALL errore( sub_name, 'radial_scale out of range ', 1 )
    IF ( radial_spread <= 0.0_DP ) &
         CALL errore( sub_name, 'radial_spread out of range ', 1 )
    IF ( filling_threshold <= 0.0_DP ) &
         CALL errore( sub_name, 'filling_threshold out of range ', 1 )
    IF ( filling_spread <= 0.0_DP ) &
         CALL errore( sub_name, 'filling_spread out of range ', 1 )
    !
    allowed = .FALSE.
    DO i = 1, SIZE( electrolyte_mode_allowed )
       IF( TRIM(electrolyte_mode) == electrolyte_mode_allowed(i) ) allowed = .TRUE.
    END DO
    IF( .NOT. allowed ) &
         CALL errore( sub_name, ' electrolyte_mode '''// &
         & TRIM(electrolyte_mode)//''' not allowed ', 1 )
    IF( electrolyte_distance < 0.0_DP ) &
         CALL errore( sub_name,' electrolyte_distance out of range ', 1 )
    IF( electrolyte_spread <= 0.0_DP ) &
         CALL errore( sub_name,' electrolyte_spread out of range ', 1 )
    IF( electrolyte_rhomax < 0.0_DP ) &
         CALL errore( sub_name,' electrolyte_rhomax out of range ', 1 )
    IF( electrolyte_rhomin < 0.0_DP ) &
         CALL errore( sub_name,' electrolyte_rhomin out of range ', 1 )
    IF( electrolyte_rhomax < electrolyte_rhomin ) &
         CALL errore( sub_name,' inconsistent electrolyte_rhomax and electrolyte_rhomin', 1 )
    IF( electrolyte_tbeta < 0.0_DP ) &
         CALL errore( sub_name,' electrolyte_tbeta out of range ', 1 )
    IF( electrolyte_alpha <= 0.0_DP ) &
         CALL errore( sub_name,' electrolyte_alpha out of range ', 1 )
    IF( electrolyte_softness <= 0.0_DP ) &
         CALL errore( sub_name,' electrolyte_softness out of range ', 1 )
    !
    allowed = .FALSE.
    DO i = 1, SIZE( boundary_core_allowed )
       IF( TRIM(boundary_core) == boundary_core_allowed(i) ) allowed = .TRUE.
    END DO
    IF( .NOT. allowed ) &
         CALL errore( sub_name, ' boundary_core '''// &
         & TRIM(core)//''' not allowed ', 1 )
    !
    IF( ifdtype < 1 ) &
         CALL errore( sub_name,' ifdtype out of range ', 1 )
    IF( nfdpoint < 1 ) &
         CALL errore( sub_name,' nfdpoint out of range ', 1 )
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE boundary_checkin
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE electrostatic_checkin()
!--------------------------------------------------------------------
    !
    !  Check input values for Namelist ELECTROSTATIC
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=20) :: sub_name = ' electrostatic_checkin '
    INTEGER           :: i
    LOGICAL           :: allowed = .FALSE.
    !
    allowed = .FALSE.
    DO i = 1, SIZE( problem_allowed )
       IF( TRIM(problem) == problem_allowed(i) ) allowed = .TRUE.
    END DO
    IF( .NOT. allowed ) &
         CALL errore( sub_name, ' problem '''// &
         & TRIM(problem)//''' not allowed ', 1 )
    IF( tol <= 0.0_DP ) &
         CALL errore( sub_name,' tolerance out of range ', 1 )
    !
    allowed = .FALSE.
    DO i = 1, SIZE( solver_allowed )
       IF( TRIM(solver) == solver_allowed(i) ) allowed = .TRUE.
    END DO
    IF( .NOT. allowed ) &
         CALL errore( sub_name, ' solver '''// &
         & TRIM(solver)//''' not allowed ',1)
    allowed = .FALSE.
    DO i = 1, SIZE( auxiliary_allowed )
       IF( TRIM(auxiliary) == auxiliary_allowed(i) ) allowed = .TRUE.
    END DO
    IF( .NOT. allowed ) &
         CALL errore( sub_name, ' auxiliary '''// &
         & TRIM(auxiliary)//''' not allowed ', 1 )
    allowed = .FALSE.
    DO i = 1, SIZE( step_type_allowed )
       IF( TRIM(step_type) == step_type_allowed(i) ) allowed = .TRUE.
    END DO
    IF( .NOT. allowed ) &
         CALL errore( sub_name, ' step_type '''// &
         & TRIM(step_type)//''' not allowed ', 1 )
    IF( step <= 0.0_DP ) &
         CALL errore( sub_name,' step out of range ', 1 )
    IF( maxstep <= 1 ) &
         CALL errore( sub_name,' maxstep out of range ', 1 )
    !
    allowed = .FALSE.
    DO i = 1, SIZE( mix_type_allowed )
       IF( TRIM(mix_type) == mix_type_allowed(i) ) allowed = .TRUE.
    END DO
    IF( .NOT. allowed ) &
         CALL errore( sub_name, ' mix_type '''// &
         & TRIM(mix_type)//''' not allowed ', 1 )
    IF( ndiis <= 0 ) &
         CALL errore( sub_name,' ndiis out of range ', 1 )
    IF( mix <= 0.0_DP ) &
         CALL errore( sub_name,' mix out of range ', 1 )
    !
    allowed = .FALSE.
    DO i = 1, SIZE( preconditioner_allowed )
       IF( TRIM(preconditioner) == preconditioner_allowed(i) ) allowed = .TRUE.
    END DO
    IF( .NOT. allowed ) &
         CALL errore( sub_name, ' preconditioner '''// &
         & TRIM(preconditioner)//''' not allowed ', 1 )
    allowed = .FALSE.
    DO i = 1, SIZE( screening_type_allowed )
       IF( TRIM(screening_type) == screening_type_allowed(i) ) allowed = .TRUE.
    END DO
    IF( .NOT. allowed ) &
         CALL errore( sub_name, ' screening_type '''// &
         & TRIM(screening_type)//''' not allowed ', 1 )
    IF( screening < 0.0_DP ) &
         CALL errore( sub_name,' screening out of range ', 1 )
    !
    allowed = .FALSE.
    DO i = 1, SIZE( core_allowed )
       IF( TRIM(core) == core_allowed(i) ) allowed = .TRUE.
    END DO
    IF( .NOT. allowed ) &
         CALL errore( sub_name, ' core '''// &
         & TRIM(core)//''' not allowed ', 1 )
    !
    IF( pbc_dim < -3 .OR. pbc_dim > 3 ) &
         CALL errore( sub_name,' pbc_dim out of range ', 1 )
    IF( pbc_axis < 1 .OR. pbc_axis > 3 ) &
         CALL errore( sub_name,' cell_axis out of range ', 1 )
    allowed = .FALSE.
    DO i = 1, SIZE( pbc_correction_allowed )
       IF( TRIM(pbc_correction) == pbc_correction_allowed(i) ) allowed = .TRUE.
    END DO
    IF( .NOT. allowed ) &
         CALL errore( sub_name, ' pbc_correction '''// &
         & TRIM(pbc_correction)//''' not allowed ', 1 )
    IF( TRIM(pbc_correction) .EQ. 'gcs' .AND. TRIM(electrolyte_mode) .NE. 'system' ) &
       & CALL errore( sub_name, 'Only system boundary for gcs correction', 1)
    allowed = .FALSE.
    DO i = 1, SIZE( inner_solver_allowed )
       IF( TRIM(inner_solver) == inner_solver_allowed(i) ) allowed = .TRUE.
    END DO
    IF( .NOT. allowed ) &
         CALL errore( sub_name, ' inner solver '''// &
         & TRIM(inner_solver)//''' not allowed ',1)
    IF( inner_mix <= 0.0_DP ) &
         CALL errore( sub_name,' inner_mix out of range ', 1 )
    IF( inner_tol <= 0.0_DP ) &
         CALL errore( sub_name,' inner_tol out of range ', 1 )
    IF( inner_maxstep <= 1 ) &
         CALL errore( sub_name,' inner_maxstep out of range ', 1 )
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE electrostatic_checkin
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE fix_boundary( lboundary )
!--------------------------------------------------------------------
    !
    !  Check if BOUNDARY needs to be read and reset defaults
    !  according to the ENVIRON namelist
    !
    IMPLICIT NONE
    !
    LOGICAL, INTENT(OUT) :: lboundary
    !
    CHARACTER(LEN=20) :: sub_name = ' fix_boundary '
    !
    lboundary = .FALSE.
    !
    IF ( environ_type .NE. 'input' .AND. environ_type .NE. 'vacuum' ) &
         & lboundary = .TRUE.
    IF ( env_static_permittivity .GT. 1.D0 .OR. env_optical_permittivity .GT. 1.D0 ) &
         & lboundary = .TRUE.
    IF ( env_surface_tension .GT. 0.D0 ) lboundary = .TRUE.
    IF ( env_pressure .NE. 0.D0 ) lboundary = .TRUE.
    IF ( env_confine .NE. 0.D0 ) lboundary = .TRUE.
    IF ( env_electrolyte_ntyp .GT. 0 ) lboundary = .TRUE.
    IF ( env_dielectric_regions .GT. 0 ) lboundary = .TRUE.
    !
    IF ( solvent_mode .EQ. 'ionic' .AND. boundary_core .NE. 'analytic' ) THEN
       IF ( ionode ) WRITE(program_unit,*)'Only analytic boundary_core for ionic solvent_mode'
       boundary_core = 'analytic'
    ENDIF
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE fix_boundary
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE set_environ_type( )
!--------------------------------------------------------------------
    !
    !  Set values according to the environ_type keyword and boundary mode
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=20) :: sub_name = ' set_environ_type '
    !
    ! Skip set up if read environ keywords from input
    !
    IF ( TRIM(ADJUSTL(environ_type)) .EQ. 'input' ) RETURN
    !
    ! Vacuum case is straightforward, all flags are off
    !
    IF ( TRIM(ADJUSTL(environ_type)) .EQ. 'vacuum' ) THEN
       env_static_permittivity = 1.D0
       env_optical_permittivity = 1.D0
       env_surface_tension = 0.D0
       env_pressure = 0.D0
       RETURN
    ENDIF
    !
    ! First set global physically meaningful parameters
    !
    SELECT CASE ( TRIM(ADJUSTL(environ_type)) )
       !
    CASE ('water', 'water-cation', 'water-anion' )
       ! water experimental permittivities
       env_static_permittivity = 78.3D0
       env_optical_permittivity = 1.D0 ! 1.776D0
    CASE DEFAULT
       call errore (sub_name,'unrecognized value for environ_type',1)
    END SELECT
    !
    ! Depending on the boundary mode, set fitted parameters
    !
    IF ( TRIM(ADJUSTL(solvent_mode)) .EQ. 'electronic' .OR. &
         & TRIM(ADJUSTL(solvent_mode)) .EQ. 'full' ) THEN
       !
       ! Self-consistent continuum solvation (SCCS)
       !
       SELECT CASE ( TRIM(ADJUSTL(environ_type)) )
          !
       CASE ( 'water' )
          ! SCCS for neutrals
          env_surface_tension = 50.D0
          env_pressure = -0.35D0
          rhomax = 0.005
          rhomin = 0.0001
       CASE ( 'water-cation' )
          ! SCCS for cations
          env_surface_tension = 5.D0
          env_pressure = 0.125D0
          rhomax = 0.0035
          rhomin = 0.0002
       CASE( 'water-anion' )
          ! SCCS for cations
          env_surface_tension = 0.D0
          env_pressure = 0.450D0
          rhomax = 0.0155
          rhomin = 0.0024
       END SELECT
       !
    ELSE IF ( solvent_mode .EQ. 'ionic' ) THEN
       !
       ! Soft-sphere continuum solvation
       !
       radius_mode = 'uff'
       softness = 0.5D0
       env_surface_tension = 50.D0 !! NOTE THAT WE ARE USING THE
       env_pressure = -0.35D0      !! SET FOR CLUSTERS, AS IN SCCS
       !
       SELECT CASE ( TRIM(ADJUSTL(environ_type)) )
          !
       CASE ( 'water' )
          ! SS for neutrals
          alpha = 1.12D0
       CASE ( 'water-cation' )
          ! SS for cations
          alpha = 1.10D0
       CASE( 'water-anion' )
          ! SS for anions
          alpha = 0.98D0
       END SELECT
       !
    END IF
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE set_environ_type
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE fix_electrostatic(lelectrostatic)
!--------------------------------------------------------------------
    !
    !  Set values according to the ENVIRON namelist
    !
    IMPLICIT NONE
    !
    LOGICAL, INTENT(OUT) :: lelectrostatic
    !
    CHARACTER(LEN=20) :: sub_name = ' fix_electrostatic '
    !
    lelectrostatic = env_electrostatic
    IF ( env_static_permittivity .GT. 1.D0 .OR. env_optical_permittivity .GT. 1.D0 ) &
         & lelectrostatic = .TRUE.
    IF ( env_external_charges .GT. 0 ) lelectrostatic = .TRUE.
    IF ( env_dielectric_regions .GT. 0 ) lelectrostatic = .TRUE.
    IF ( env_electrolyte_ntyp .GT. 0 ) lelectrostatic = .TRUE.
    !
!--------------------------------------------------------------------
  END SUBROUTINE fix_electrostatic
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE set_electrostatic_problem( )
!--------------------------------------------------------------------
    !
    !  Set problem according to the ENVIRON and ELECTROSTATIC namelists
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=80) :: sub_name = ' set_electrostatic_problem '
    !
    IF ( env_electrolyte_ntyp .GT. 0 ) THEN
       IF ( .NOT. TRIM(pbc_correction) == 'gcs' ) THEN
          IF ( electrolyte_linearized ) THEN
             IF (problem == 'none') problem = 'linpb'
             IF (solver == 'none' ) solver  = 'cg'
             IF ( cionmax .GT. 0.D0 .OR. rion .GT. 0.D0 ) problem = 'linmodpb'
          ELSE
             IF (problem == 'none') problem = 'pb'
             IF (solver == 'none' ) solver  = 'newton'
             IF (inner_solver == 'none') inner_solver = 'cg'
             IF ( cionmax .GT. 0.D0 .OR. rion .GT. 0.D0 ) problem = 'modpb'
          END IF
       END IF
    END IF
    !
    IF ( env_static_permittivity > 1.D0 &
         .OR. env_dielectric_regions > 0 ) THEN
         IF (problem == 'none') problem = 'generalized'
       IF ( .NOT. TRIM(pbc_correction) == 'gcs' ) THEN
         IF (solver == 'none' ) solver = 'cg'
       ELSE
         IF (solver == 'none' ) solver = 'iterative'
         IF (solver == 'iterative' &
            & .AND. auxiliary == 'none' ) auxiliary = 'full'
         IF (solver .NE. 'iterative') &
            & CALL errore( sub_name, 'GCS correction requires iterative solver', 1)
       END IF
    ELSE
       IF (problem == 'none') problem = 'poisson'
       IF (solver == 'none' ) solver = 'direct'
    ENDIF
    !
    IF (.NOT. (problem == 'pb' .OR. problem == 'modpb') &
        .AND. (inner_solver .NE. 'none')) &
        CALL errore( sub_name, 'Only pb or modpb problems allow inner solver', 1)
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE set_electrostatic_problem
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE environ_read_cards( unit )
!--------------------------------------------------------------------
    !
    !  Environ cards parsing routine
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN), optional  :: unit
    !
    CHARACTER(len=256)         :: input_line
    CHARACTER(len=80)          :: card
    CHARACTER(len=1), EXTERNAL :: capital
    LOGICAL                    :: tend
    INTEGER                    :: i
    !
    parse_unit = unit
    !
    ! CALL environ_card_default_values( )
    !
100 CALL env_read_line( input_line, end_of_file=tend )
    !
    IF( tend ) GOTO 120
    IF( input_line == ' ' .OR. input_line(1:1) == '#' .OR. &
         input_line(1:1) == '!' ) GOTO 100
    !
    READ (input_line, *) card
    !
    DO i = 1, len_trim( input_line )
       input_line( i : i ) = capital( input_line( i : i ) )
    ENDDO
    !
    IF ( trim(card) == 'EXTERNAL_CHARGES' ) THEN
       !
       CALL card_external_charges( input_line )
       !
    ELSE IF ( trim(card) == 'DIELECTRIC_REGIONS' ) THEN
       !
       CALL card_dielectric_regions( input_line )
       !
    ELSE
       !
       IF ( ionode ) WRITE( program_unit,'(A)') 'Warning: card '//trim(input_line)//' ignored'
       !
    ENDIF
    !
    ! ... END OF LOOP ... !
    !
    GOTO 100
    !
120 CONTINUE
    !
    ! ... Check
    !
    IF ( env_external_charges .GT. 0 .AND. .NOT. taextchg ) &
         CALL errore( ' environ_read_cards  ', ' missing card external_charges', 0 )
    IF ( env_dielectric_regions .GT. 0 .AND. .NOT. taepsreg ) &
         CALL errore( ' environ_read_cards  ', ' missing card dielectric_regions', 0 )
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE environ_read_cards
!--------------------------------------------------------------------
!--------------------------------------------------------------------
   SUBROUTINE card_external_charges( input_line )
!--------------------------------------------------------------------
     !
     ! ... Description of the allowed input CARDS
     !
     ! EXTERNAL_CHARGES (unit_option)
     !
     !   set external fixed charge densities and their shape
     !
     ! Syntax:
     !
     !    EXTERNAL_CHARGES (unit_option)
     !      charge(1)  x(1) y(1) z(1)  spread(1) dim(1)  axis(1)
     !       ...       ...        ...      ...        ...
     !      charge(n)  x(n) y(n) z(n)  spread(n) dim(n)  axis(n)
     !
     ! Example:
     !
     ! EXTERNAL_CHARGES (bohr)
     !  1.0  0.0  0.0  0.0  [0.5  2  1]
     ! -1.0  0.0  0.0  5.0  [0.5  2  1]
     !
     ! Where:
     !
     !   unit_option == bohr       positions are given in Bohr (DEFAULT)
     !   unit_option == angstrom   positions are given in Angstrom
     !
     !      charge(i) ( real )       total charge of the density
     !      x/y/z(i)  ( real )       cartesian position of the density
     !      spread(i) ( real )       gaussian spread of the density (in bohr, optional, default=0.5)
     !      dim(i)    ( integer )    0/1/2 point/line/plane of charge (optional, default=0)
     !      axis(i)   ( integer )    1/2/3 for x/y/z direction of line/plane (optional, default=3)
     !
     !USE wrappers, ONLY: feval_infix
     !
     IMPLICIT NONE
     !
     CHARACTER(len=256) :: input_line
     INTEGER            :: ie, ix, ierr, nfield
     LOGICAL            :: tend
     LOGICAL, EXTERNAL  :: matches
     CHARACTER(len=4)   :: lb_pos
     CHARACTER(len=256) :: field_str
     !
     IF ( taextchg ) THEN
        CALL errore( ' card_external_charges  ', ' two occurrences', 2 )
     ENDIF
     IF ( env_external_charges > nsx ) THEN
        CALL errore( ' card_external_charges ', ' nsx out of range ', env_external_charges )
     ENDIF
     !
     CALL allocate_input_extcharge(env_external_charges)
     !
     IF ( matches( "BOHR", input_line ) ) THEN
        external_charges = 'bohr'
     ELSEIF ( matches( "ANGSTROM", input_line ) ) THEN
        external_charges = 'angstrom'
     ELSE
        IF ( trim( adjustl( input_line ) ) /= 'EXTERNAL_CHARGES' ) THEN
           CALL errore( 'read_cards ', &
                & 'unknown option for EXTERNAL_CHARGES: '&
                & // input_line, 1 )
        ENDIF
        CALL infomsg( 'read_cards ', &
             & 'No units specified in EXTERNAL_CHARGES card' )
        external_charges = 'bohr'
        CALL infomsg( 'read_cards ', &
             & 'EXTERNAL_CHARGES: units set to '//TRIM(external_charges) )
     ENDIF
     !
     DO ie = 1, env_external_charges
        !
        CALL env_read_line( input_line, end_of_file = tend )
        IF ( tend ) CALL errore( 'environ_cards', &
             'end of file reading external charges', ie )
        !
        CALL env_field_count( nfield, input_line )
        !
        ! ... read field 1 (total charge of the external density)
        !
        CALL get_field(1, field_str, input_line)
        !extcharge_charge(ie) = feval_infix(ierr, field_str )
        read(field_str,*) extcharge_charge(ie)
        !
        ! ... read fields 2-4 (x-y-z position of external density)
        !
        CALL get_field(2, field_str, input_line)
        !extcharge_pos(1,ie) = feval_infix(ierr, field_str )
        read(field_str,*) extcharge_pos(1,ie)
        CALL get_field(3, field_str, input_line)
        !extcharge_pos(2,ie) = feval_infix(ierr, field_str )
        read(field_str,*) extcharge_pos(2,ie)
        CALL get_field(4, field_str, input_line)
        !extcharge_pos(3,ie) = feval_infix(ierr, field_str )
        read(field_str,*) extcharge_pos(3,ie)
        !
        ! ... optionally read field 5 (spread of the density)
        !
        IF ( nfield >= 5 ) THEN
           CALL get_field(5, field_str, input_line)
           !extcharge_spread(ie) = feval_infix(ierr, field_str )
           read(field_str,*) extcharge_spread(ie)
           IF ( extcharge_spread(ie) .LT. 0.D0 ) &
                CALL errore( ' card_external_charges  ', ' spread must be positive', ie )
        ENDIF
        !
        ! ... optionally read field 6 and 7 (dimensionality and direction)
        !
        IF ( nfield >= 6 ) THEN
           CALL env_get_field(6, field_str, input_line)
           READ(field_str, *) extcharge_dim(ie)
           IF ( extcharge_dim(ie) .LT. 0 .OR. extcharge_dim(ie) .GT. 2 ) &
                CALL errore( ' card_external_charges  ', ' wrong excharge dimension ', ie )
           IF ( extcharge_dim(ie) .GT. 0 ) THEN
              IF ( nfield == 6 ) &
                   CALL errore('environ_cards',&
                   'missing axis direction of partially periodic external charge', ie)
              CALL env_get_field(7, field_str, input_line)
              READ(field_str, *) extcharge_axis(ie)
              IF ( extcharge_axis(ie) .LT. 0 .OR. extcharge_axis(ie) .GT. 3 ) &
                   CALL errore( ' card_external_charges  ', ' wrong excharge axis ', ie )
           ENDIF
        ENDIF
        !
     ENDDO
     taextchg = .true.
     !
     DO ie = 1, env_external_charges
        DO ix = 1, 3
           CALL convert_length( external_charges, extcharge_pos(ix, ie))
        ENDDO
        CALL convert_length( external_charges, extcharge_spread(ie))
     ENDDO
     !
     RETURN
     !
!--------------------------------------------------------------------
   END SUBROUTINE card_external_charges
!--------------------------------------------------------------------
!--------------------------------------------------------------------
   SUBROUTINE allocate_input_extcharge(env_external_charges)
!--------------------------------------------------------------------
     !
     IMPLICIT NONE
     !
     INTEGER, INTENT(in) :: env_external_charges
     !
     IF ( allocated( extcharge_dim    ) ) DEALLOCATE( extcharge_dim    )
     IF ( allocated( extcharge_axis   ) ) DEALLOCATE( extcharge_axis   )
     IF ( allocated( extcharge_charge ) ) DEALLOCATE( extcharge_charge )
     IF ( allocated( extcharge_spread ) ) DEALLOCATE( extcharge_spread )
     IF ( allocated( extcharge_pos    ) ) DEALLOCATE( extcharge_pos    )
     !
     ALLOCATE( extcharge_dim    ( env_external_charges ) )
     ALLOCATE( extcharge_axis   ( env_external_charges ) )
     ALLOCATE( extcharge_charge ( env_external_charges ) )
     ALLOCATE( extcharge_spread ( env_external_charges ) )
     ALLOCATE( extcharge_pos ( 3, env_external_charges ) )
     !
     extcharge_dim    = 0
     extcharge_axis   = 3
     extcharge_charge = 0.0_DP
     extcharge_spread = 0.5_DP
     extcharge_pos    = 0.0_DP
     !
     RETURN
     !
!--------------------------------------------------------------------
   END SUBROUTINE allocate_input_extcharge
!--------------------------------------------------------------------
   !
   !----------------------------------------------------------------------
   !
   !----------------------------------------------------------------------
   !
!--------------------------------------------------------------------
   SUBROUTINE card_dielectric_regions( input_line )
!--------------------------------------------------------------------
     !
     !
     ! ... Description of the allowed input CARDS
     !
     ! DIELECTRIC_REGIONS (unit_option)
     !
     !   set fixed dielectric regions and their shape
     !
     ! Syntax:
     !
     !    DIELECTRIC_REGIONS (unit_option)
     !      epsilon0(1) epsilonopt(1) x(1) y(1) z(1)  width(1) spread(1) dim(1)  axis(1)
     !       ...       ...        ...      ...        ...
     !      epsilon0(n) epsilonopt(n) x(n) y(n) z(n)  width(n) spread(n) dim(n)  axis(n)
     !
     ! Example:
     !
     ! DIELECTRIC_REGIONS (bohr)
     !  80.0  2.0   0.0  0.0  10.0   5.0  1.0  2  3
     !
     ! Where:
     !
     !   unit_option == bohr       positions are given in Bohr (DEFAULT)
     !   unit_option == angstrom   positions are given in Angstrom
     !
     !      epsilon0(i)   ( real )    static permittivity inside the region
     !      epsilonopt(i) ( real )    optical permittivity inside the region
     !      x/y/z(i)      ( real )    cartesian center of the region
     !      width(i)      ( real )    size of the region (in bohr)
     !      spread(i)     ( real )    spread of the interface (in bohr, optional)
     !      dim(i)     ( integer )    0/1/2 point/line/plane region (optional)
     !      axis(i)    ( integer )    1/2/3 for x/y/z direction of line/plane (optional)
     !
     !USE wrappers, ONLY: feval_infix
     !
     IMPLICIT NONE
     !
     CHARACTER(len=256) :: input_line
     INTEGER            :: ie, ix, ierr, nfield
     LOGICAL            :: tend
     LOGICAL, EXTERNAL  :: matches
     CHARACTER(len=4)   :: lb_pos
     CHARACTER(len=256) :: field_str
     !
     IF ( taepsreg ) THEN
        CALL errore( ' card_dielectric_regions  ', ' two occurrences', 2 )
     ENDIF
     IF ( env_dielectric_regions > nsx ) THEN
        CALL errore( ' card_dielectric_regions ', ' nsx out of range ', env_dielectric_regions )
     ENDIF
     !
     CALL allocate_input_epsregion(env_dielectric_regions)
     !
     IF ( matches( "BOHR", input_line ) ) THEN
        dielectric_regions = 'bohr'
     ELSEIF ( matches( "ANGSTROM", input_line ) ) THEN
        dielectric_regions = 'angstrom'
     ELSE
        IF ( trim( adjustl( input_line ) ) /= 'DIELECTRIC_REGIONS' ) THEN
           CALL errore( 'read_cards ', &
                & 'unknown option for DIELECTRIC_REGIONS: '&
                & // input_line, 1 )
        ENDIF
        CALL infomsg( 'read_cards ', &
             & 'No units specified in DIELECTRIC_REGIONS card' )
        dielectric_regions = 'bohr'
        CALL infomsg( 'read_cards ', &
             & 'DIELECTRIC_REGIONS: units set to '//TRIM(dielectric_regions) )
     ENDIF
     !
     DO ie = 1, env_dielectric_regions
        !
        CALL env_read_line( input_line, end_of_file = tend )
        IF ( tend ) CALL errore( 'environ_cards', &
             'end of file reading dielectric regions', ie )
        !
        CALL env_field_count( nfield, input_line )
        !
        ! ... read field 1-2 (static and optical permettivity inside dielectric region)
        !
        CALL get_field(1, field_str, input_line)
        !epsregion_eps(1,ie) = feval_infix(ierr, field_str )
        READ(field_str, *) epsregion_eps(1,ie)
        IF ( epsregion_eps(1,ie) .LT. 1.D0 ) &
             CALL errore( ' card_dielectric_regions  ', ' static permittivity must be .gt. 1', ie )
        CALL get_field(2, field_str, input_line)
        !epsregion_eps(2,ie) = feval_infix(ierr, field_str )
        READ(field_str, *) epsregion_eps(2,ie)
        IF ( epsregion_eps(2,ie) .LT. 1.D0 ) &
             CALL errore( ' card_dielectric_regions  ', ' optical permittivity must be .gt. 1', ie )
        !
        ! ... read fields 3-5 (x-y-z position of dielectric region)
        !
        CALL get_field(3, field_str, input_line)
        !epsregion_pos(1,ie) = feval_infix(ierr, field_str )
        READ(field_str, *) epsregion_pos(1,ie)
        CALL get_field(4, field_str, input_line)
        !epsregion_pos(2,ie) = feval_infix(ierr, field_str )
        READ(field_str, *) epsregion_pos(2,ie)
        CALL get_field(5, field_str, input_line)
        !epsregion_pos(3,ie) = feval_infix(ierr, field_str )
        READ(field_str, *) epsregion_pos(3,ie)
        !
        ! ... read field 6 (size/width of the dielectric region)
        !
        CALL get_field(6, field_str, input_line)
        !epsregion_width(ie) = feval_infix(ierr, field_str )
        READ(field_str, *) epsregion_width(ie)
        IF ( epsregion_width(ie) .LT. 0.D0 ) &
             CALL errore( ' card_dielectric_regions  ', ' width must be positive', ie )
        !
        ! ... optionally read field 7 (spread of interface of the dielectric region)
        !
        IF ( nfield >= 7 ) THEN
           CALL get_field(7, field_str, input_line)
           !epsregion_spread(ie) = feval_infix(ierr, field_str )
           READ(field_str, *) epsregion_spread(ie)
           IF ( epsregion_spread(ie) .LT. 0.D0 ) &
                CALL errore( ' card_dielectric_regions ', ' spread must be positive', ie )
        ENDIF
        !
        ! ... optionally read field 7 and 8 (dimensionality and direction)
        !
        IF ( nfield >= 8 ) THEN
           CALL env_get_field(8, field_str, input_line)
           READ(field_str, *) epsregion_dim(ie)
           IF ( epsregion_dim(ie) .LT. 0 .OR. epsregion_dim(ie) .GT. 2 ) &
                CALL errore( ' card_dielectric_regions ', ' wrong epsregion dimension ', ie )
           IF ( epsregion_dim(ie) .GT. 0 ) THEN
              IF ( nfield == 8 ) &
                   CALL errore('environ_cards',&
                   'missing axis direction of partially periodic dielectric region', ie)
              CALL env_get_field(9, field_str, input_line)
              READ(field_str, *) epsregion_axis(ie)
              IF ( epsregion_axis(ie) .LT. 1 .OR. epsregion_axis(ie) .GT. 3 ) &
                   CALL errore( ' card_dielectric_regions ', ' wrong epsregion axis ', ie )
           ENDIF
        ENDIF
        !
     ENDDO
     taepsreg = .true.
     !
     DO ie = 1, env_dielectric_regions
        DO ix = 1, 3
           CALL convert_length( dielectric_regions, epsregion_pos(ix, ie))
        ENDDO
        CALL convert_length( dielectric_regions, epsregion_width(ie))
        CALL convert_length( dielectric_regions, epsregion_spread(ie))
     ENDDO
     !
     RETURN
     !
!--------------------------------------------------------------------
   END SUBROUTINE card_dielectric_regions
!--------------------------------------------------------------------
!--------------------------------------------------------------------
   SUBROUTINE allocate_input_epsregion(env_dielectric_regions)
!--------------------------------------------------------------------
     !
     IMPLICIT NONE
     !
     INTEGER, INTENT(in) :: env_dielectric_regions
     !
     IF ( allocated( epsregion_dim    ) ) DEALLOCATE( epsregion_dim    )
     IF ( allocated( epsregion_axis   ) ) DEALLOCATE( epsregion_axis   )
     IF ( allocated( epsregion_eps    ) ) DEALLOCATE( epsregion_eps    )
     IF ( allocated( epsregion_width  ) ) DEALLOCATE( epsregion_width  )
     IF ( allocated( epsregion_spread ) ) DEALLOCATE( epsregion_spread )
     IF ( allocated( epsregion_pos    ) ) DEALLOCATE( epsregion_pos    )
     !
     ALLOCATE( epsregion_dim    ( env_dielectric_regions ) )
     ALLOCATE( epsregion_axis   ( env_dielectric_regions ) )
     ALLOCATE( epsregion_eps ( 2, env_dielectric_regions ) )
     ALLOCATE( epsregion_width  ( env_dielectric_regions ) )
     ALLOCATE( epsregion_spread ( env_dielectric_regions ) )
     ALLOCATE( epsregion_pos ( 3, env_dielectric_regions ) )
     !
     epsregion_dim    = 0
     epsregion_axis   = 3
     epsregion_eps    = 1.0_DP
     epsregion_width  = 0.0_DP
     epsregion_spread = 0.5_DP
     epsregion_pos    = 0.0_DP
     !
     RETURN
     !
!--------------------------------------------------------------------
   END SUBROUTINE allocate_input_epsregion
!--------------------------------------------------------------------
!--------------------------------------------------------------------
   SUBROUTINE convert_length(length_format, length)
!--------------------------------------------------------------------
     !
     ! ... convert input length to atomic units
     !
     IMPLICIT NONE
     CHARACTER (len=*), INTENT(in)  :: length_format
     REAL (DP), INTENT(inout) :: length
     !
     SELECT CASE( length_format )
     CASE( 'bohr' )
        !
        ! ... input length are in a.u., do nothing
        !
        length = length
        !
     CASE( 'angstrom' )
        !
        ! ... length in A: convert to a.u.
        !
        length = length / bohr_radius_angs
        !
     CASE DEFAULT
        !
        CALL errore( 'iosys','length_format=' // &
             & trim( length_format ) // ' not implemented', 1 )
        !
     END SELECT
     !
!--------------------------------------------------------------------
   END SUBROUTINE convert_length
!--------------------------------------------------------------------
!  !
!  ! Two more wrappers for eval_infix (simple algebric expression parser)
!  !
!  FUNCTION feval_infix(fierr, fstr)
!     USE ISO_C_BINDING
!     IMPLICIT NONE
!     REAL(DP) :: feval_infix
!     INTEGER :: fierr
!     CHARACTER(len=*) :: fstr
!     INTEGER :: filen
!     !
!     INTERFACE
!     FUNCTION ceval_infix(cierr, cstr, cilen) BIND(C, name="eval_infix")
!     !REAL(kind=c_double) FUNCTION ceval_infix(cierr, cstr, cilen) BIND(C, name="eval_infix")
!     !  double eval_infix( int *ierr, const char *strExpression, int len )
!       USE ISO_C_BINDING
!       REAL(kind=c_double) :: ceval_infix
!       INTEGER(kind=c_int)    :: cierr
!       CHARACTER(kind=c_char) :: cstr(*)
!       INTEGER(kind=c_int),VALUE :: cilen
!     END FUNCTION ceval_infix
!     END INTERFACE
!     !
!     INTEGER(kind=c_int) :: cierr
!     INTEGER(kind=c_int) :: cilen
!     CHARACTER(len=len_trim(fstr)+1,kind=c_char) :: cstr
!     !
!     INTEGER :: i
!     !
!     filen = len_trim(fstr)
!     cilen = INT(filen, kind=c_int)
!     DO i = 1,filen
!       cstr(i:i) = fstr(i:i)
!     ENDDO
!     cstr(filen+1:filen+1)=C_NULL_CHAR
!     !
!     feval_infix = REAL( ceval_infix(cierr, cstr, cilen), kind=DP)
!     fierr = INT(cierr)
!     RETURN
!   END FUNCTION feval_infix
!!----------------------------------------------------------------------------
END MODULE environ_input
!----------------------------------------------------------------------------
