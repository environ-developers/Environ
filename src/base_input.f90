!----------------------------------------------------------------------------------------
!
! Copyright (C) 2018 ENVIRON (www.quantum-environ.org)
!
!----------------------------------------------------------------------------------------
!
!     This file is part of Environ version 2.0
!
!     Environ 2.0 is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation, either version 2 of the License, or
!     (at your option) any later version.
!
!     Environ 2.0 is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more detail, either the file
!     `License' in the root directory of the present distribution, or
!     online at <http://www.gnu.org/licenses/>.
!
!----------------------------------------------------------------------------------------
!
! Authors: Oliviero Andreussi (Department of Physics, UNT)
!          Francesco Nattino  (THEOS and NCCR-MARVEL, EPFL)
!          Ismaila Dabo       (DMSE, Penn State)
!          Nicola Marzari     (THEOS and NCCR-MARVEL, EPFL)
!          Edan Bainglass     (Department of Physics, UNT)
!
!----------------------------------------------------------------------------------------
!>
!! This module contains all variables in the environ.in input file
!!
!----------------------------------------------------------------------------------------
MODULE env_base_input
    !------------------------------------------------------------------------------------
    !
    USE environ_param, ONLY: DP
    !
    !------------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    !=---------------------------------------------------------------------------------=!
!     ENVIRON Cards Parameters
    !=---------------------------------------------------------------------------------=!
    !
    ! Local parameters of external charges
    !
    LOGICAL :: taextchg = .FALSE.
    !
    CHARACTER(LEN=80) :: extcharge_units = 'bohr' ! atomic positions (bohr|angstrom)
    !
    REAL(DP), ALLOCATABLE :: extcharge_charge(:) ! total charge of density
    REAL(DP), ALLOCATABLE :: extcharge_pos(:, :) ! cartesian position of density
    REAL(DP), ALLOCATABLE :: extcharge_spread(:) ! gaussian density spread (bohr)
    INTEGER, ALLOCATABLE :: extcharge_dim(:) ! point/line/plane of charge
    INTEGER, ALLOCATABLE :: extcharge_axis(:) ! x/y/z direction of line/plane
    !
    !------------------------------------------------------------------------------------
    ! Local parameters of dielectric regions
    !
    LOGICAL :: taepsreg = .FALSE.
    !
    CHARACTER(LEN=80) :: epsregion_units = 'bohr' ! atomic positions (bohr|angstrom)
    !
    REAL(DP), ALLOCATABLE :: epsregion_eps(:, :) ! permittivity inside region
    REAL(DP), ALLOCATABLE :: epsregion_pos(:, :) ! cartesian center of region
    REAL(DP), ALLOCATABLE :: epsregion_width(:) ! region size (bohr)
    REAL(DP), ALLOCATABLE :: epsregion_spread(:) ! interface spread (bohr)
    INTEGER, ALLOCATABLE :: epsregion_dim(:) ! point/line/plane region
    INTEGER, ALLOCATABLE :: epsregion_axis(:) ! x/y/z direction of line/plane
    !
    !=---------------------------------------------------------------------------------=!
!     ENVIRON Namelist Input Parameters
    !=---------------------------------------------------------------------------------=!
    !
    LOGICAL :: environ_debug = .FALSE.
    !
    LOGICAL :: environ_restart = .FALSE.
    ! restart a previous calculation: environ contributions are computed during
    ! initialization
    !
    INTEGER :: verbose = 0 ! verbosity
    ! 0: only prints summary of polarization charge calculation;
    ! 1: prints an extra file with details of iterative convergence;
    ! 2: prints 3D cube files of physical properties
    !
    REAL(DP) :: environ_thr = 1.D-1 ! when in scf to start calculating corr. pot.
    !
    INTEGER :: environ_nskip = 1 ! # steps to skip before starting add. pot. computation
    !
    !------------------------------------------------------------------------------------
    ! Energy cutoff used for internal FFT-grid generation
    !
    REAL(DP) :: env_ecut = 0.D0
    ! may be used when coupled with non-FFT calling programs or via the tester program
    !
    !------------------------------------------------------------------------------------
    ! Predefined environ types
    !
    CHARACTER(LEN=80) :: environ_type = 'input'
    CHARACTER(LEN=80) :: environ_type_allowed(5)
    !
    DATA environ_type_allowed/'vacuum', 'water', 'water-cation', 'water-anion', 'input'/
    !
    ! sets all the environment parameters at once to a specific set
    !
    ! vacuum = all the flags are off (perm=1.d0, surf=0.0, pres=0.0)
    !
    ! water = parameters optimized for water solutions in Andreussi et al.
    !         J. Chem. Phys. 136, 064102 (perm=78, surf=50, pres=-0.35)
    !
    ! water-cation = parameters optimized for aqueous solvation of cations
    !         Dupont et al. J. Chem. Phys. 139, 214110 (perm=78, surf=, pres=)
    !
    ! water-anion = parameters optimized for aqueous solvation of anions
    !         Dupont et al. J. Chem. Phys. 139, 214110 (perm=78, surf=, pres=)
    !
    ! input = do not use any predefined set, use parameters from input
    !
    !------------------------------------------------------------------------------------
    ! System specification
    !
    INTEGER :: system_ntyp = 0
    ! specify the atom types that are used to determine the origin and
    ! size of the system (types up to system_ntyp are used, all atoms are
    ! used by default or if system_ntyp == 0)
    !
    INTEGER :: system_dim = 0
    ! dimensionality of the system, used to determine size (only ortogonally to
    ! periodic dimensions) and position (0 = 0D, 1 = 1D, 2 = 2D)
    !
    INTEGER :: system_axis = 3 ! main axis of 1D or 2D systems (1 = x, 2 = y, 3 = z)
    !
    !------------------------------------------------------------------------------------
    ! Environment cell specifications
    !
    INTEGER :: env_nrep(3) = 0
    ! number of additional replicas of the system cell on each side along the three axis
    ! nrep = 1 means there is one more cell on the left and on the right of the cell
    ! the environment cell is (2*nrep+1) times the system cell along the three axis
    !
    REAL(DP) :: system_pos(3) = 0.D0 ! specified for finite-difference debugging
    ! if not specified, system is fixed on center of charge
    !
    !------------------------------------------------------------------------------------
    ! Modification of electrostatic embedding (e.g. PBC correction)
    !
    LOGICAL :: env_electrostatic = .FALSE. ! flag electrostatic namelist reading
    !
    REAL(DP), ALLOCATABLE :: atomicspread(:)
    ! atomic charge density gaussian spread (a.u.)
    !
    !------------------------------------------------------------------------------------
    ! Dielectric solvent parameters
    !
    REAL(DP) :: env_static_permittivity = 1.D0
    ! static dielectric permittivity of the solvation model
    ! if set equal to one (=vacuum), no dielectric effects
    !
    REAL(DP) :: env_optical_permittivity = 1.D0
    ! optical dielectric permittivity of the solvation model
    ! if set equal to one (=vacuum), no dielectric effects
    ! needed only for the TDDFTPT
    !
    !------------------------------------------------------------------------------------
    ! Cavitation energy parameters
    !
    REAL(DP) :: env_surface_tension = 0.D0 ! solvent surface tension
    ! if equal to zero, no cavitation term
    !
    !------------------------------------------------------------------------------------
    ! PV energy parameters
    !
    REAL(DP) :: env_pressure = 0.D0 ! external pressure for PV energy
    ! if equal to zero no pressure term
    !
    !------------------------------------------------------------------------------------
    ! Confine energy parameters
    !
    REAL(DP) :: env_confine = 0.D0 ! confinement potential
    !
    !------------------------------------------------------------------------------------
    ! Ionic countercharge parameters
    !
    LOGICAL :: electrolyte_linearized = .FALSE.
    ! solve linear-regime poisson-boltzmann problem
    !
    INTEGER :: env_electrolyte_ntyp = 0
    ! number of countercharge species in the electrolyte
    ! if != 0, must be >= 2
    !
    CHARACTER(LEN=80) :: electrolyte_entropy = 'full'
    CHARACTER(LEN=80) :: electrolyte_entropy_allowed(2)
    !
    DATA electrolyte_entropy_allowed/'ions', 'full'/
    !
    ! sets the electrolyte entropy terms that are affected by the Stern-layer correction
    !
    ! ions = only ionic terms ( Ringe et al. J. Chem. Theory Comput. 12, 4052 )
    !
    ! full = all terms ( Dabo et al. arXiv 0901.0096 )
    !
    REAL(DP), ALLOCATABLE :: cion(:) ! molar concentration of ionic countercharge (M=mol/L)
    REAL(DP) :: cionmax = 1.D3 ! maximum molar concentration of ionic countercharge (M=mol/L)
    REAL(DP) :: rion = 0.D0 ! mean atomic radius of ionic countercharge (a.u.)
    REAL(DP), ALLOCATABLE :: zion(:) ! valence of ionic countercharge
    REAL(DP) :: temperature = 300.D0 ! temperature of the solution
    !
    !------------------------------------------------------------------------------------
    ! Semiconductor parameters
    !
    REAL(DP) :: sc_permittivity = 1.D0 ! dielectric permittivity of the semiconductor
    !
    REAL(DP) :: sc_carrier_density = 0.D0
    ! concentration of charge carriers within the semiconductor (cm^-3)
    !
    REAL(DP) :: sc_electrode_chg = 0.D0 ! the total charge on the electrode (e)
    !
    REAL(DP) :: sc_chg_thr = 1.D-4
    ! threshold for an outer loop of chg optimization in qe
    !
    !------------------------------------------------------------------------------------
    ! External charges parameters not read from EXTERNAL_CHARGES card
    !
    INTEGER :: env_external_charges = 0
    ! number of fixed external gaussian points/lines/planes of charges to be used
    ! in the calculation
    !
    !------------------------------------------------------------------------------------
    ! Dielectric regions parameters not read from DIELECTRIC_REGIONS card
    !
    INTEGER :: env_dielectric_regions = 0
    ! number of fixed dielectric regions in the calculation
    !
    !------------------------------------------------------------------------------------
    !
    NAMELIST /environ/ &
        environ_debug, environ_restart, verbose, environ_thr, environ_nskip, env_ecut, &
        environ_type, system_ntyp, system_dim, system_axis, env_nrep, system_pos, &
        env_electrostatic, atomicspread, env_static_permittivity, &
        env_optical_permittivity, env_surface_tension, env_pressure, env_confine, &
        env_electrolyte_ntyp, cion, cionmax, rion, zion, temperature, &
        electrolyte_linearized, electrolyte_entropy, sc_permittivity, &
        sc_carrier_density, sc_electrode_chg, sc_chg_thr, env_external_charges, &
        env_dielectric_regions
    !
    !=---------------------------------------------------------------------------------=!
!     BOUNDARY Namelist Input Parameters
    !=---------------------------------------------------------------------------------=!
    !
    ! Soft boundary (electronic) parameters
    !
    INTEGER :: stype = 2 ! type of switching functions used in the solvation models
    ! 0: original Fattebert-Gygi
    ! 1: ultrasoft switching function (only exponential part used for non-electrostatic)
    ! 2: ultrasoft switching function as defined in Andreussi et al. JCP 2012
    !
    !------------------------------------------------------------------------------------
    ! Rigid boundary (ionic) parameters
    !
    CHARACTER(LEN=80) :: radius_mode = 'uff'
    CHARACTER(LEN=80) :: radius_mode_allowed(4)
    !
    DATA radius_mode_allowed/'pauling', 'bondi', 'uff', 'muff'/
    !
    ! type of hardcoded solvation radii to be used when solvent_mode = 'ionic'
    !
    ! pauling = R.C. Weast, ed., Handbook of chemistry and physics
    !           (CRC Press, Cleveland, 1981)
    !
    ! bondi   = A. Bondi, J. Phys. Chem. 68, 441 (1964)
    !
    ! uff     = A.K. Rapp/'{e} et al. J. Am. Chem. Soc. 114(25) pp.10024-10035 (1992)
    !
    ! muff    = uff with local modifications (Nitrogen, see Fisicaro JCTC (2017)
    !
    REAL(DP), ALLOCATABLE :: solvationrad(:)
    ! solvationrad radius of the solvation shell for each species when the
    ! ionic dielectric function is adopted, in internal units (a.u.)
    !
    !------------------------------------------------------------------------------------
    ! Full boundary parameters
    !
    REAL(DP), ALLOCATABLE :: corespread(:)
    ! gaussian spreads of the core electrons, in internal units (a.u.), to
    ! be used when solvent_mode = 'full'
    !
    !------------------------------------------------------------------------------------
    ! Solvent-aware boundary parameters
    !
    REAL(DP) :: solvent_radius = 0.D0
    ! size of the solvent, used to decide whether to fill a continuum
    ! void or not. If set equal to 0.D0, use the standard algorithm
    !
    REAL(DP) :: radial_scale = 2.D0
    ! compute the filled fraction on a spherical volume scaled w.r.t solvent size
    !
    REAL(DP) :: radial_spread = 0.5D0
    ! spread of the step function used to evaluate occupied volume
    !
    REAL(DP) :: filling_threshold = 0.825D0
    ! threshold to decide whether to fill a continuum void or not, to be
    ! compared with the filled fraction: if filled fraction > threshold
    ! THEN fill gridpoint
    !
    REAL(DP) :: filling_spread = 0.02D0
    ! spread of the switching function used to decide whether the continuum
    ! void should be filled or not
    !
    !------------------------------------------------------------------------------------
    ! Field-aware boundary parameters
    !
    LOGICAL :: field_aware = .FALSE.
    ! switch to turn on the field-awareness scaling factor
    ! NOTE: only works with ionic boundary mode
    !
    REAL(DP) :: field_factor = 0.08D0
    ! maximum scaling factor possible by the field-aware model
    !
    REAL(DP) :: field_asymmetry = -0.32D0
    ! charge asymmetry factor. Positive values result in more field-awareness
    ! for positive charges, and vice versa
    !
    REAL(DP) :: field_max = 6.D0
    ! maximum flux value for switching function
    !
    REAL(DP) :: field_min = 2.D0
    ! minimum flux value for switching function
    !
    !------------------------------------------------------------------------------------
    ! Derivative core's parameters
    !
    CHARACTER(LEN=80) :: deriv_method = 'default'
    CHARACTER(LEN=80) :: deriv_method_allowed(5)
    !
    DATA deriv_method_allowed/'default', 'fft', 'chain', 'highmem', 'lowmem'/
    !
    ! algorithms for computing derivatives
    !
    ! fft       = fast Fourier transforms
    !
    ! chain     = chain-rule derivatives for as much as possible (FFTs for the rest)
    !
    ! highmem   = analytic derivatives for soft-sphere computed by storing all spherical
    !             functions and derivatives
    !
    ! lowmem    = more efficient analytic derivatives
    !
    CHARACTER(LEN=80) :: deriv_core = 'fft'
    CHARACTER(LEN=80) :: deriv_core_allowed(1)
    !
    DATA deriv_core_allowed/'fft'/
    !
    ! choice of the core numerical methods to be exploited for derivatives
    !
    ! fft = fast Fourier transforms (default)
    !
    ! to be implemented : wavelets (from big-DFT) and multigrid #TODO future work
    !
    !------------------------------------------------------------------------------------
    ! Solvent boundary parameters
    !
    CHARACTER(LEN=80) :: solvent_mode = 'electronic'
    CHARACTER(LEN=80) :: solvent_mode_allowed(5)
    !
    DATA solvent_mode_allowed/'electronic', 'ionic', 'full', 'external', 'system'/
    !
    ! solvent_mode method for calculating the density that sets the dielectric constant
    !
    ! electronic = dielectric depends self-consist. on electronic density
    !
    ! ionic = dielectric defined on a fictitious ionic density, generated
    !         as the sum of spherical error functions centered on atomic
    !         positions of width specified in input by solvationrad(ityp)
    !
    ! full  = similar to electronic, but an extra density is added to
    !         represent the core electrons and the nuclei. This extra
    !         density is defined as the sum of gaussian functions centered
    !         on atomic positions of width equal to corespread(ityp)
    !
    ! system = simplified regular dielectric defined to be outside a distance
    !          solvent_distance from the specified system
    !
    !------------------------------------------------------------------------------------
    ! Soft solvent boundary (electronic) parameters
    !
    REAL(DP) :: rhomax = 0.005D0
    ! first parameter of the sw function, roughly corresponding to the density
    ! threshold of the solvation model
    !
    REAL(DP) :: rhomin = 0.0001D0 ! second parameter of the sw function when stype=1 or 2
    !
    REAL(DP) :: tbeta = 4.8D0 ! second parameter of the sw function when stype=0
    !
    !------------------------------------------------------------------------------------
    ! Rigid solvent boundary (ionic) parameters
    !
    REAL(DP) :: alpha = 1.D0 ! scaling factor for ionic radii when solvent_mode = 'ionic'
    REAL(DP) :: softness = 0.5D0 ! spread of the rigid interfaces
    !
    !------------------------------------------------------------------------------------
    ! Simplified solvent boundary (system) parameters
    !
    REAL(DP) :: solvent_distance = 1.D0
    ! distance from the system where the boundary starts if required from solvent_mode
    !
    REAL(DP) :: solvent_spread = 0.5D0
    ! spread of the boundary interface if defined on system position and width
    !
    !------------------------------------------------------------------------------------
    ! Stern boundary parameters
    !
    CHARACTER(LEN=80) :: electrolyte_mode = 'electronic'
    CHARACTER(LEN=80) :: electrolyte_mode_allowed(5)
    !
    DATA electrolyte_mode_allowed/'electronic', 'ionic', 'full', 'external', &
        'system'/
    !
    ! electrolyte_mode method for calculating the density that sets the onset of
    ! ionic countercharge ( see solvent_mode above )
    !
    !------------------------------------------------------------------------------------
    ! Soft Stern boundary (electronic) parameters
    !
    REAL(DP) :: electrolyte_rhomax = 0.005D0
    ! first parameter of the Stern sw function, roughly corresponding
    ! to the density threshold of the ionic countercharge.
    !
    REAL(DP) :: electrolyte_rhomin = 0.0001D0
    ! second parameter of the Stern sw function when stype=1 or 2
    !
    REAL(DP) :: electrolyte_tbeta = 4.8D0
    ! second parameter of the Stern sw function when stype=0
    !
    !------------------------------------------------------------------------------------
    ! Rigid Stern boundary (ionic) parameters
    !
    REAL(DP) :: electrolyte_alpha = 1.D0
    ! scaling factor for ionic radii when electrolyte_mode = 'ionic'
    !
    REAL(DP) :: electrolyte_softness = 0.5D0 ! spread of the rigid Stern interfaces
    !
    !------------------------------------------------------------------------------------
    ! Simplified Stern boundary (system) parameters
    !
    REAL(DP) :: electrolyte_distance = 0.D0
    ! distance from the system where the electrolyte boundary starts
    !
    REAL(DP) :: electrolyte_spread = 0.5D0
    ! spread of the interfaces for the electrolyte boundary
    !
    !------------------------------------------------------------------------------------
    ! Electrolyte boundary derivatives
    !
    CHARACTER(LEN=80) :: electrolyte_deriv_method = 'default'
    CHARACTER(LEN=80) :: electrolyte_deriv_method_allowed(5)
    !
    DATA electrolyte_deriv_method_allowed/'default', 'fft', 'chain', 'highmem', 'lowmem'/
    !
    ! algorithms for computing derivatives on the electrolyte boundary
    !
    ! fft       = fast Fourier transforms
    !
    ! chain     = chain-rule derivatives for as much as possible (FFTs for the rest)
    !
    ! highmem   = analytic derivatives for soft-sphere computed by storing all spherical
    !             functions and derivatives
    !
    ! lowmem    = more efficient analytic derivatives
    !
    !------------------------------------------------------------------------------------
    ! Mott Schottky boundary (system parameters
    !
    REAL(DP) :: sc_distance = 0.D0
    ! distance from the system where the mott schottky boundary starts
    !
    REAL(DP) :: sc_spread = 0.5D0
    ! spread of the interfaces for the mott schottky boundary
    !
    !------------------------------------------------------------------------------------
    !
    NAMELIST /boundary/ &
        solvent_mode, radius_mode, alpha, softness, solvationrad, stype, rhomax, &
        rhomin, tbeta, corespread, solvent_distance, solvent_spread, solvent_radius, &
        radial_scale, radial_spread, filling_threshold, filling_spread, field_aware, &
        field_factor, field_asymmetry, field_max, field_min, electrolyte_mode, &
        electrolyte_distance, electrolyte_spread, electrolyte_rhomax, &
        electrolyte_rhomin, electrolyte_tbeta, electrolyte_alpha, &
        electrolyte_softness, deriv_method, deriv_core, sc_distance, sc_spread
    !
    !=---------------------------------------------------------------------------------=!
!     ELECTROSTATIC Namelist Input Parameters
    !=---------------------------------------------------------------------------------=!
    !
    CHARACTER(LEN=80) :: problem = 'none'
    CHARACTER(LEN=80) :: problem_allowed(7)
    !
    DATA problem_allowed/ &
        'none', 'poisson', 'generalized', 'pb', 'modpb', 'linpb', 'linmodpb'/
    !
    ! type of electrostatic problem
    !
    ! poisson     = standard poisson equation, with or without
    !               boundary conditions (default)
    !
    ! generalized = generalized poisson equation
    !
    ! pb          = poisson-boltzmann equation (non-linear)
    !
    ! modpb       = modified poisson-boltzmann equation (non-linear)
    !
    ! linpb       = linearized poisson-boltzmann equation (debye-huckel)
    !
    ! linmodpb    = linearized modified poisson-boltzmann equation
    !
    CHARACTER(LEN=80) :: inner_problem = 'none'
    !
    ! type of electrostatic problem for inner loop in nested algorithms
    !
    ! generalized = generalized poisson equation
    !
    ! linpb       = linearized poisson-boltzmann equation (debye-huckel)
    !
    !------------------------------------------------------------------------------------
    !
    REAL(DP) :: tol = 1.D-5
    ! convergence threshold for electrostatic potential or auxiliary charge
    !
    REAL(DP) :: inner_tol = 1.D-10 ! same as tol for inner loop in nested algorithms
    !
    !------------------------------------------------------------------------------------
    ! Driver's parameters
    !
    CHARACTER(LEN=80) :: solver = 'none'
    CHARACTER(LEN=80) :: solver_allowed(8)
    !
    DATA solver_allowed/ &
        'none', 'cg', 'sd', 'fixed-point', 'lbfgs', 'newton', 'nested', 'direct'/
    !
    ! type of numerical solver
    !
    ! direct    = for simple problems with analytic or direct solution
    !
    ! cg        = conjugate gradient (default)
    !
    ! sd        = steepest descent
    !
    ! fp        = fixed-point search
    !
    ! lbfgs     = low-memory bfgs
    !
    ! newton    = newton's method (only for non-linear problem)
    !
    ! nested    = double iterations (only for non-linear problem)
    !
    CHARACTER(LEN=80) :: auxiliary = 'none'
    CHARACTER(LEN=80) :: auxiliary_allowed(4)
    !
    DATA auxiliary_allowed/'none', 'full', 'pol', 'ioncc'/
    !
    ! solve with respect to the potential or with respect to an auxiliary charge density
    !
    ! none  = solve for the potential (default)
    !
    ! full  = solve for the auxiliary charge density
    !
    ! pol   = in a nested scheme, solve the inner (pol) cycle in terms of
    !         the auxiliary charge
    !
    ! ioncc = in a nested scheme, solve the outer (ioncc) cycle in terms of
    !         the auxiliary charge
    !
    CHARACTER(LEN=80) :: step_type = 'optimal'
    CHARACTER(LEN=80) :: step_type_allowed(3)
    !
    DATA step_type_allowed/'optimal', 'input', 'random'/
    !
    ! how to choose the step size in gradient descent algorithms or iterative mixing
    !
    ! optimal = step size that minimize the cost function on the descent direction
    !
    ! input   = fixed step size as defined in input (step keyword)
    !
    ! random  = random step size within zero and twice the optima value
    !
    REAL(DP) :: step = 0.3D0
    ! step size to be used if step_type = 'input'
    ! (inherits the tasks of the old mixrhopol)
    !
    INTEGER :: maxstep = 200
    ! maximum number of steps to be performed by gradient or iterative solvers
    !
    CHARACTER(LEN=80) :: inner_solver = 'none'
    CHARACTER(LEN=80) :: inner_solver_allowed(5)
    !
    DATA inner_solver_allowed/'none', 'cg', 'sd', 'fixed-point', 'direct'/
    ! type of numerical solver for inner loop in nested algorithms
    !
    INTEGER :: inner_maxstep = 200
    ! same as maxstep for inner loop in nested algorithms
    !
    !------------------------------------------------------------------------------------
    ! Iterative driver's parameters (OBSOLETE)
    !
    CHARACTER(LEN=80) :: mix_type = 'linear'
    CHARACTER(LEN=80) :: mix_type_allowed(4)
    !
    DATA mix_type_allowed/'linear', 'anderson', 'diis', 'broyden'/
    ! mixing method for iterative calculations: linear | anderson | diis | broyden
    !
    INTEGER :: ndiis = 1 ! order of DIIS interpolation of iterative calculation
    REAL(DP) :: mix = 0.5D0 ! mixing parameter to be used in the iterative driver
    REAL(DP) :: inner_mix = 0.5D0 ! same as mix but for inner loop in nested algorithm
    !
    !------------------------------------------------------------------------------------
    ! Preconditioner's parameters
    !
    CHARACTER(LEN=80) :: preconditioner = 'sqrt'
    CHARACTER(LEN=80) :: preconditioner_allowed(3)
    !
    DATA preconditioner_allowed/'none', 'sqrt', 'left'/
    !
    ! type of preconditioner
    !
    ! none      = no preconditioner
    !
    ! left      = left linear preconditioner eps nabla v = r
    !
    ! sqrt      = sqrt preconditioner sqrt(eps) nabla ( sqrt(eps) * v ) = r
    !
    CHARACTER(LEN=80) :: screening_type = 'none'
    CHARACTER(LEN=80) :: screening_type_allowed(4)
    !
    DATA screening_type_allowed/'none', 'input', 'linear', 'optimal'/
    !
    ! use the screened coulomb Green's function instead of the vacuum one
    !
    ! none      = unscreened coulomb
    !
    ! input     = screened coulomb with screening length provided in input
    !
    ! linear    = screened coulomb with screening length from linear component
    !             of the problem
    !
    ! optimal   = screened coulomb with screening length optimized (to be defined)
    !
    REAL(DP) :: screening = 0.D0
    ! screening length to be used if screening_type = 'input'
    !
    !------------------------------------------------------------------------------------
    ! Numerical core's parameters
    !
    CHARACTER(LEN=80) :: core = 'fft'
    CHARACTER(LEN=80) :: core_allowed(1)
    !
    DATA core_allowed/'fft'/
    !
    ! choice of the core numerical methods to be exploited for electrostatics
    !
    ! fft = fast Fourier transforms (default)
    !
    ! to be implemented : wavelets (from big-DFT) and multigrid #TODO future work
    !
    !------------------------------------------------------------------------------------
    ! Inner numerical core's parameters
    !
    CHARACTER(LEN=80) :: inner_core = 'fft'
    CHARACTER(LEN=80) :: inner_core_allowed(1)
    !
    DATA inner_core_allowed/'fft'/
    !
    ! choice of the core numerical methods to be exploited for nested electrostatics
    !
    ! fft = fast Fourier transforms (default)
    !
    ! to be implemented : wavelets (from big-DFT) and multigrid #TODO future work
    !
    !------------------------------------------------------------------------------------
    ! Periodic correction keywords
    !
    CHARACTER(LEN=80) :: pbc_correction = 'none'
    CHARACTER(LEN=80) :: pbc_correction_allowed(4)
    !
    DATA pbc_correction_allowed/'none', 'parabolic', 'gcs', 'ms'/
    !
    ! type of periodic boundary condition correction to be used
    !
    ! parabolic = point-counter-charge type of correction
    !
    ! gcs       = Gouy-Chapman-Stern correction for electrolyte
    !
    ! ms        = mott-schottky correction for semiconductor
    !
    INTEGER :: pbc_dim = -3 ! dimensionality of the simulation cell
    ! periodic boundary conditions on 3/2/1/0 sides of the cell
    !
    INTEGER :: pbc_axis = 3 ! choice of the sides with periodic boundary conditions
    ! 1 = x, 2 = y, 3 = z, where
    ! if pbc_dim = 2, cell_axis is orthogonal to 2D plane
    ! if pbc_dim = 1, cell_axis is along the 1D direction
    !
    CHARACTER(LEN=80) :: pbc_core = '1da'
    CHARACTER(LEN=80) :: pbc_core_allowed(1)
    !
    DATA pbc_core_allowed/'1da'/
    !
    ! choice of the core numerical methods to be exploited for pbc corrections
    !
    ! 1da = 1d-analytic
    !
    !------------------------------------------------------------------------------------
    !
    NAMELIST /electrostatic/ &
        problem, tol, solver, auxiliary, step_type, step, maxstep, mix_type, mix, &
        ndiis, preconditioner, screening_type, screening, core, pbc_dim, &
        pbc_correction, pbc_axis, pbc_core, inner_tol, inner_solver, inner_maxstep, &
        inner_mix, inner_core
    !
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE env_base_input
!----------------------------------------------------------------------------------------
