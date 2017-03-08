!
! Copyright (C) 2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!=----------------------------------------------------------------------------=!
!
MODULE environ_input
!
!=----------------------------------------------------------------------------=!
!  this module contains all variables in namelist &ENVIRON and related routines
!  performing initialization and broadcast - Written by Oliviero Andreussi
!=----------------------------------------------------------------------------=!
  !
  USE kinds,      ONLY : DP
  USE parameters, ONLY : nsx
  !
  USE io_global,  ONLY : ionode, ionode_id, stdout
  USE parser,     ONLY : field_count, read_line, get_field, parse_unit
  USE mp,         ONLY : mp_bcast
  USE mp_images,  ONLY : intra_image_comm
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
  INTEGER :: extcharge_origin = 0
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
  INTEGER :: epsregion_origin = 0
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
        CHARACTER( LEN = 80 ) :: environ_type = 'input'
        CHARACTER( LEN = 80 ) :: environ_type_allowed(5)
        DATA enviton_type_allowed / 'vacuum', 'water', 'water-cation', &
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
! Switching function parameters
!
        INTEGER :: stype = 1
        ! type of switching functions used in the solvation models
        !    0: original Fattebert-Gygi
        !    1: ultrasoft dielectric function as defined in Andreussi et al.
        REAL(DP) :: rhomax = 0.005
        ! first parameter of the sw function, roughly corresponding
        ! to the density threshold of the solvation model
        REAL(DP) :: rhomin = 0.0001
        ! second parameter of the sw function when stype=1
        REAL(DP) :: tbeta = 4.8
        ! second parameter of the sw function when stype=0
!
! Dielectric solvent parameters
!
        REAL(DP) :: env_static_permittivity = 1.D0
        ! static dielectric permittivity of the solvation model. If set equal
        ! to one (=vacuum) no dielectric effects
        REAL(DP) :: env_optical_permittivity = 1.D0
        ! optical dielectric permittivity of the solvation model. If set equal
        ! to one (=vacuum) no dielectric effects. Needed only for the TDDFTPT.
        CHARACTER( LEN = 80 ) :: eps_mode = 'electronic'
        CHARACTER( LEN = 80 ) :: eps_mode_allowed(4)
        DATA eps_mode_allowed / 'electronic', 'ionic', 'full', 'external' /
        ! eps_mode method for calculating the density that sets
        ! the dielectric constant
        ! electronic = dielectric depends self-consist. on electronic density
        ! ionic = dielectric defined on a fictitious ionic density, generated
        !         as the sum of exponential functions centered on atomic
        !         positions of width specified in input by solvationrad(ityp)
        ! full  = similar to electronic, but an extra density is added to
        !         represent the core electrons and the nuclei. This extra
        !         density is defined as the sum of gaussian functions centered
        !         on atomic positions of width equal to corespread(ityp)
        REAL(DP) :: alpha = 1.D0
        ! scaling factor for ionic radii when eps_mode = 'ionic'
        REAL(DP) :: solvationrad(nsx) = 3.D0
        ! solvationrad radius of the solvation shell for each species when the
        ! ionic dielectric function is adopted, in internal units (a.u.)
        REAL(DP) :: corespread(nsx) = 0.5D0
        ! gaussian spreads of the core electrons, in internal units (a.u.), to
        ! be used when eps_mode = 'full'
        REAL(DP) :: atomicspread(nsx) = 0.5D0
        ! gaussian spreads of the atomic density of charge, in internal units (a.u.)
        LOGICAL :: add_jellium = .false.
        ! depending on periodic boundary corrections, one may need to explicitly
        ! polarize the compensatinig jellium background
!
! Cavitation energy parameters
!
        REAL(DP) :: env_surface_tension = 0.D0
        ! solvent surface tension, if equal to zero no cavitation term
        REAL(DP) :: delta = 0.00001D0
        ! finite difference parameter to compute the molecular surface
!
! PV energy parameters
!
        REAL(DP) :: env_pressure = 0.D0
        ! external pressure for PV energy, if equal to zero no pressure term
!
! Ionic countercharge parameters
!
        INTEGER :: env_ioncc_level = 0
        ! level of accuracy in ioncc
        INTEGER :: nrep = 0
        ! number of replicas of unit cell along slab_axis
        CHARACTER( LEN = 80 ) :: stern_mode = 'electronic'
        CHARACTER( LEN = 80 ) :: stern_mode_allowed(4)
        DATA stern_mode_allowed / 'electronic', 'ionic', 'full', 'external' /
        ! stern_mode method for calculating the density that sets
        ! the onset of ionic countercharge
        ! electronic = onset depends self-consist. on electronic density
        ! ionic = onset defined on a fictitious ionic density, generated
        !         as the sum of exponential functions centered on atomic
        !         positions of width specified in input by solvationrad(ityp)
        ! full  = similar to electronic, but an extra density is added to
        !         represent the core electrons and the nuclei. This extra
        !         density is defined as the sum of gaussian functions centered
        !         on atomic positions of width equal to corespread(ityp)
        REAL(DP) :: stern_distance = 0.D0
        ! onset distance of countercharge, if ioncc_level = 1|2
        REAL(DP) :: stern_spread = 0.5D0
        ! spread of countercharge onset, if ioncc_level = 2
        REAL(DP) :: cion = 1.D0
        ! molar concentration of ionic countercharge (M=mol/L)
        REAL(DP) :: cionmax = 1.D3
        ! maximum molar concentration of ionic countercharge (M=mol/L)
        REAL(DP) :: rion = 0.D0
        ! mean atomic radius of ionic countercharge (a.u.)
        REAL(DP) :: zion = 1.D0
        ! valence of ionic countercharge
        REAL(DP) :: rhopb = 0.0001D0
        ! density threshold for the onset of ionic countercharge
        REAL(DP) :: solvent_temperature = 300.D0
        ! temperature of the solution
!
! Periodic correction keywords
!
        INTEGER :: env_periodicity = 3
        ! dimensionality of the simulation cell
        ! periodic boundary conditions on 3/2/1/0 sides of the cell
        CHARACTER( LEN = 80 ) :: pbc_correction = 'parabolic'
        CHARACTER( LEN = 80 ) :: pbc_correction_allowed(1)
        DATA pbc_correction_allowed / 'parabolic' /
        ! type of periodic boundary condition correction to be used
        ! parabolic = point-counter-charge type of correction
        INTEGER :: cell_axis = 3
        ! choice of the sides with periodic boundary conditions
        ! 1 = x, 2 = y, 3 = z, where
        ! if env_periodicity = 2, cell_axis is orthogonal to 2D plane
        ! if env_periodicity = 1, cell_axis is along the 1D direction
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

        NAMELIST /environ/                                             &
             environ_restart, verbose, environ_thr, environ_nskip,     &
             environ_type,                                             &
             system_ntyp, system_dim, system_axis,                     &
             stype, rhomax, rhomin, tbeta,                             &
             env_static_permittivity, env_optical_permittivity,        &
             eps_mode,                                                 &
             alpha, solvationrad, corespread, atomicspread,            &
             add_jellium,                                              &
             env_surface_tension, delta,                               &
             env_pressure,                                             &
             env_ioncc_level, nrep,                                    &
             stern_mode, stern_distance, stern_spread,                 &
             cion, cionmax, rion, zion, rhopb, solvent_temperature,    &
             env_periodicity, pbc_correction, cell_axis,               &
             env_external_charges, env_dielectric_regions
!
!=----------------------------------------------------------------------------=!
!  ELECTROSTATIC Namelist Input Parameters
!=----------------------------------------------------------------------------=!
!
! Global parameters
!
        CHARACTER( LEN = 80 ) :: problem = 'poisson'
        CHARACTER( LEN = 80 ) :: problem_allowed(7)
        DATA problem_allowed / 'poisson', 'free', 'generalized', 'pb', 'modpb', &
          'linpb', 'linmodpb' /
        ! type of electrostatic problem:
        ! poisson     = standard poisson equation with pbc (default, no electrostatic correction)
        ! free        = standard poitton equation with free bc
        ! generalized = generalized poisson equation
        ! pb          = poisson-boltzmann equation (non-linear)
        ! modpb       = modified poisson-boltzmann equation (non-linear)
        ! linpb       = linearized poisson-boltzmann equation (debye-huckel)
        ! linmodpb    = linearized modified poisson-boltzmann equation
        REAL(DP) :: tolvelect = 1.D-6
        ! convergence threshold for electrostatic potential
        REAL(DP) :: tolrhoaux = 1.D-10
        ! convergence threshold for auxiliary charge
!
! Driver's parameters
!
        CHARACTER( LEN = 80 ) :: solver = 'cg'
        CHARACTER( LEN = 80 ) :: solver_allowed(6)
        DATA solver_allowed / 'cg', 'sd', 'iterative', 'lbfgs', 'newton', 'nested' /
        ! type of numerical solver
        ! cg        = conjugate gradient (default)
        ! sd        = steepest descent
        ! iterative = fixed-point search (OBSOLETE)
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
        DATA step_allowed / 'optimal', 'input', 'random' /
        ! how to choose the step size in gradient descent algorithms or iterative mixing
        ! optimal = step size that minimize the cost function on the descent direction
        ! input   = fixed step size as defined in input (step keyword)
        ! random  = random step size within zero and twice the optima value
        REAL(DP) :: step = 0.3
        ! step size to be used if step_type = 'input' (inherits the tasks of the old mixrhopol)
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
        REAL(DP) :: mix = 0.3
        ! mixing parameter to be used in the iterative driver
!
! Preconditioner's parameters
!
        CHARACTER( LEN = 80 ) :: preconditioner = 'right'
        CHARACTER( LEN = 80 ) :: preconditioner_allowed(4)
        DATA preconditioner_allowed / 'none', 'right', 'sqrt', 'left' /
        ! type of preconditioner
        ! none      = no preconditioner
        ! right     = right linear preconditioner eps nabla v = r
        ! sqrt      = sqrt preconditioner sqrt(eps) nabla ( sqrt(eps) * v ) = r
        ! left      = left linear preconditioner nabla ( eps * v ) = r
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
        CHARACTER( LEN = 80 ) :: core_allowed(2)
        DATA core_allowed / 'fft', 'fd' /
        ! choice of the core numerical methods to be exploited for the different operations
        ! fft = fast Fourier transforms (default)
        ! fd  = finite differences in real space
        ! to be implemented : wavelets (from big-DFT) and multigrid
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
! Boundary conditions
!
        INTEGER, DIMENSION(3) :: bc_index = (/ 0, 0, 0 /)
        ! boundary conditions on each main axis
        ! 0 = periodic
        ! 1 = free ! SOME CASES ONLY
        ! 2 = dirichelet  ! TO IMPLEMENT
        ! 3 = von neumann ! TO IMPLEMENT
        ! 4 = electrochem ! TO IMPLEMENT
        REAL(DP), DIMENSION(3) :: bcplus  = (/ 0.D0, 0.D0, 0.D0 /)
        REAL(DP), DIMENSION(3) :: bcminus = (/ 0.D0, 0.D0, 0.D0 /)
        ! values of the fixed boundary conditions at the positive and negative boundaries
        ! in internal units

        NAMELIST /electrostatic/                 &
             problem, tolvelec, tolrhoaux,       &
             solver, auxiliary, step_type, step, &
             mix_type, mix, ndiis,               &
             preconditioner,                     &
             screening_type, screening,          &
             core,                               &
             ifdtype, nfdpoint,                  &
             bc_index, bcplus, bcminus

  CONTAINS
     !
     SUBROUTINE read_environ(nat,ntyp,assume_isolated,ibrav)
       !
       USE environ_base, ONLY : set_environ_base
       USE environ_init, ONLY : environ_initions_allocate
       !
       CHARACTER(len=80), INTENT(IN) :: assume_isolated
       INTEGER, INTENT(IN) :: nat, ntyp, ibrav
       !
       INTEGER, EXTERNAL :: find_free_unit
       !
       LOGICAL :: ext
       INTEGER :: environ_unit_input
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
       CALL environ_read_namelist( environ_unit_input, assume_isolated, ibrav )
       !
       ! ... Read environ cards
       !
       CALL environ_read_cards( environ_unit_input )
       !
       ! ... Close environ input file
       !
       CLOSE( environ_unit_input )
       !
       ! ... Set module variables according to input
       !
       CALL set_environ_base  ( assume_isolated, environ_restart,           &
                                verbose, environ_thr, environ_nskip,        &
                                environ_type,                               &
                                system_ntyp, system_dim, system_axis,       &
                                stype, rhomax, rhomin, tbeta,               &
                                env_static_permittivity,                    &
                                env_optical_permittivity, eps_mode,         &
                                alpha, solvationrad(1:ntyp),                &
                                corespread(1:ntyp), atomicspread(1:ntyp),   &
                                add_jellium,                                &
                                env_surface_tension, delta,                 &
                                env_pressure,                               &
                                env_ioncc_level, nrep,                      &
                                stern_mode, stern_distance, stern_spread,   &
                                cion, cionmax, rion, zion, rhopb,           &
                                solvent_temperature,                        &
                                env_periodicity, pbc_correction, cell_axis, &
                                env_external_charges, extcharge_origin,     &
                                extcharge_charge, extcharge_dim,            &
                                extcharge_axis, extcharge_pos,              &
                                extcharge_spread,                           &
                                env_dielectric_regions, epsregion_origin,   &
                                epsregion_eps, epsregion_dim,               &
                                epsregion_axis, epsregion_pos,              &
                                epsregion_spread, epsregion_width )
       !
       CALL environ_initions_allocate( nat, ntyp )
       !
       RETURN
       !
     END SUBROUTINE read_environ
     !
     !=----------------------------------------------------------------------=!
     !
     !  Environ namelist parsing routine
     !
     !=----------------------------------------------------------------------=!
     !
     !-----------------------------------------------------------------------
     SUBROUTINE environ_read_namelist( environ_unit_input, assume_isolated, ibrav )
       !-----------------------------------------------------------------------
       !
       IMPLICIT NONE
       !
       CHARACTER(len=80), INTENT(IN) :: assume_isolated
       INTEGER, INTENT(IN) :: environ_unit_input, ibrav
       !
       INTEGER :: ios
       !
       ! ... Set the defauls
       !
       CALL environ_defaults()
       !
       CALL electrostatic_defaults()
       !
       ! ... Read the &ENVIRON namelist
       !
       ios = 0
       IF( ionode ) READ( environ_unit_input, environ, iostat = ios )
       CALL mp_bcast( ios, ionode_id, intra_image_comm )
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
       ! ... Fix some &ELECTROSTATIC defaults depending on the &ENVIRON input
       !
       CALL fixelect()
       !
       ! ... Read the &ELECTROSTATIC namelist
       !
       ios = 0
       IF( ionode ) READ( environ_unit_input, electrostatic, iostat = ios )
       CALL mp_bcast( ios, ionode_id, intra_image_comm )
       IF( ios /= 0 ) CALL errore( ' read_environ ', &
                                 & ' reading namelist electrostatic ', ABS(ios) )
       !
       ! ... Broadcast &ELECTROSTATIC variables
       !
       CALL electrostatic_bcast()
       !
       ! ... Check &ELECTROSTATIC variables
       !
       CALL electrostatic_checkin()
       !
       RETURN
       !
     END SUBROUTINE environ_read_namelist
     !
     !=----------------------------------------------------------------------=!
     !
     !  Variables initialization for Namelist ENVIRON
     !
     !=----------------------------------------------------------------------=!
     !
     !-----------------------------------------------------------------------
     SUBROUTINE environ_defaults( )
       !-----------------------------------------------------------------------
       !
       IMPLICIT NONE
       !
       !
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
       stype   = 1
       rhomax  = 0.005
       rhomin  = 0.0001
       tbeta   = 4.8
       !
       env_static_permittivity = 1.D0
       env_optical_permittivity = 1.D0
       eps_mode        = 'electronic'
       alpha           = 1.D0
       solvationrad(:) = 3.D0
       corespread(:)   = 0.5D0
       atomicspread(:) = 0.5D0
       add_jellium = .false.
       !
       env_surface_tension = 0.D0
       delta = 0.00001D0
       !
       env_pressure = 0.D0
       !
       env_ioncc_level = 0
       nrep = 0
       stern_mode = 'electronic'
       stern_distance = 0.D0
       stern_spread = 0.5D0
       cion = 1.0D0
       cionmax = 1.0D3
       rion = 0.D0
       zion = 1.0D0
       rhopb = 0.0001D0
       solvent_temperature = 300.0D0
       !
       env_periodicity = 3          !!! using assume_isolated instead
       pbc_correction = 'parabolic' !!! using assume_isolated instead
       cell_axis = 3                !!! using assume_isolated instead
       !
       env_external_charges = 0
       env_dielectric_regions = 0
       !
       RETURN
       !
     END SUBROUTINE
     !
     !=----------------------------------------------------------------------=!
     !
     !  Variables initialization for Namelist ELECTROSTATIC
     !
     !=----------------------------------------------------------------------=!
     !
     !-----------------------------------------------------------------------
     SUBROUTINE electrostatic_defaults( )
       !-----------------------------------------------------------------------
       !
       IMPLICIT NONE
       !
       !
       problem = 'poisson'
       tolvelect = 1.D-6
       tolrhoaux = 1.D-10
       !
       solver = 'cg'
       auxiliary = 'none'
       step_type = 'optimal'
       step = 0.3D0
       !
       mix_type = 'linear'
       ndiis = 1
       mix = 0.3D0
       !
       preconditioner = 'right'
       screening_type = 'none'
       screening = 0.D0
       !
       core = 'fft'
       ifdtype  = 1
       nfdpoint = 2
       !
       bcindex = ( / 0, 0, 0 /)
       bcplus  = ( / 0.D0, 0.D0, 0.D0 / )
       bcminus = ( / 0.D0, 0.D0, 0.D0 / )
       !
       RETURN
       !
     END SUBROUTINE
     !
     !=----------------------------------------------------------------------=!
     !
     !  Broadcast variables values for Namelist ENVIRON
     !
     !=----------------------------------------------------------------------=!
     !
     !-----------------------------------------------------------------------
     SUBROUTINE environ_bcast()
       !-----------------------------------------------------------------------
       !
       IMPLICIT NONE
       !
       CALL mp_bcast( environ_restart,            ionode_id, intra_image_comm )
       CALL mp_bcast( verbose,                    ionode_id, intra_image_comm )
       CALL mp_bcast( environ_thr,                ionode_id, intra_image_comm )
       CALL mp_bcast( environ_nskip,              ionode_id, intra_image_comm )
       CALL mp_bcast( environ_type,               ionode_id, intra_image_comm )
       !
       CALL mp_bcast( system_ntyp,                ionode_id, intra_image_comm )
       CALL mp_bcast( system_dim,                 ionode_id, intra_image_comm )
       CALL mp_bcast( system_axis,                ionode_id, intra_image_comm )
       !
       CALL mp_bcast( stype,                      ionode_id, intra_image_comm )
       CALL mp_bcast( rhomax,                     ionode_id, intra_image_comm )
       CALL mp_bcast( rhomin,                     ionode_id, intra_image_comm )
       CALL mp_bcast( tbeta,                      ionode_id, intra_image_comm )
       !
       CALL mp_bcast( env_static_permittivity,    ionode_id, intra_image_comm )
       CALL mp_bcast( env_optical_permittivity,   ionode_id, intra_image_comm )
       CALL mp_bcast( eps_mode,                   ionode_id, intra_image_comm )
       CALL mp_bcast( alpha,                      ionode_id, intra_image_comm )
       CALL mp_bcast( solvationrad,               ionode_id, intra_image_comm )
       CALL mp_bcast( corespread,                 ionode_id, intra_image_comm )
       CALL mp_bcast( atomicspread,               ionode_id, intra_image_comm )
       CALL mp_bcast( add_jellium,                ionode_id, intra_image_comm )
       !
       CALL mp_bcast( env_surface_tension,        ionode_id, intra_image_comm )
       CALL mp_bcast( delta,                      ionode_id, intra_image_comm )
       !
       CALL mp_bcast( env_pressure,               ionode_id, intra_image_comm )
       !
       CALL mp_bcast( env_ioncc_level,            ionode_id, intra_image_comm )
       CALL mp_bcast( nrep,                       ionode_id, intra_image_comm )
       CALL mp_bcast( stern_mode,                 ionode_id, intra_image_comm )
       CALL mp_bcast( stern_distance,             ionode_id, intra_image_comm )
       CALL mp_bcast( stern_spread,               ionode_id, intra_image_comm )
       CALL mp_bcast( cion,                       ionode_id, intra_image_comm )
       CALL mp_bcast( cionmax,                    ionode_id, intra_image_comm )
       CALL mp_bcast( rion,                       ionode_id, intra_image_comm )
       CALL mp_bcast( zion,                       ionode_id, intra_image_comm )
       CALL mp_bcast( rhopb,                      ionode_id, intra_image_comm )
       CALL mp_bcast( solvent_temperature,        ionode_id, intra_image_comm )
       !
       CALL mp_bcast( env_periodicity,            ionode_id, intra_image_comm )
       CALL mp_bcast( pbc_correction,             ionode_id, intra_image_comm )
       CALL mp_bcast( cell_axis,                  ionode_id, intra_image_comm )
       !
       CALL mp_bcast( env_external_charges,       ionode_id, intra_image_comm )
       CALL mp_bcast( env_dielectric_regions,     ionode_id, intra_image_comm )
       !
       RETURN
       !
     END SUBROUTINE environ_bcast
     !
     !=----------------------------------------------------------------------=!
     !
     !  Broadcast variables values for Namelist ELECTROSTATIC
     !
     !=----------------------------------------------------------------------=!
     !
     !-----------------------------------------------------------------------
     SUBROUTINE electrostatic_bcast()
       !-----------------------------------------------------------------------
       !
       IMPLICIT NONE
       !
       CALL mp_bcast( problem,                    ionode_id, intra_image_comm )
       CALL mp_bcast( tolvelect,                  ionode_id, intra_image_comm )
       CALL mp_bcast( tolrhoaux,                  ionode_id, intra_image_comm )
       !
       CALL mp_bcast( solver,                     ionode_id, intra_image_comm )
       CALL mp_bcast( auxiliary,                  ionode_id, intra_image_comm )
       CALL mp_bcast( step_type,                  ionode_id, intra_image_comm )
       CALL mp_bcast( step,                       ionode_id, intra_image_comm )
       CALL mp_bcast( mix_type,                   ionode_id, intra_image_comm )
       CALL mp_bcast( mix,                        ionode_id, intra_image_comm )
       CALL mp_bcast( ndiis,                      ionode_id, intra_image_comm )
       !
       CALL mp_bcast( preconditioner,             ionode_id, intra_image_comm )
       CALL mp_bcast( screening_type,             ionode_id, intra_image_comm )
       CALL mp_bcast( screening,                  ionode_id, intra_image_comm )
       !
       CALL mp_bcast( core,                       ionode_id, intra_image_comm )
       CALL mp_bcast( ifdtype,                    ionode_id, intra_image_comm )
       CALL mp_bcast( nfdpoint,                   ionode_id, intra_image_comm )
       !
       CALL mp_bcast( bcindex,                    ionode_id, intra_image_comm )
       CALL mp_bcast( bcplus,                     ionode_id, intra_image_comm )
       CALL mp_bcast( bcminus,                    ionode_id, intra_image_comm )
       !
       RETURN
       !
     END SUBROUTINE environ_bcast
     !
     !=----------------------------------------------------------------------=!
     !
     !  Check input values for Namelist ENVIRON
     !
     !=----------------------------------------------------------------------=!
     !
     !-----------------------------------------------------------------------
     SUBROUTINE environ_checkin()
       !-----------------------------------------------------------------------
       !
       IMPLICIT NONE
       !
       CHARACTER(LEN=20) :: sub_name = ' environ_checkin '
       INTEGER           :: i
       LOGICAL           :: allowed = .FALSE.
       !
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
       IF( stype > 1 ) &
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
       IF( env_static_permittivity < 1.0_DP ) &
          CALL errore( sub_name,' env_static_permittivity out of range ', 1 )
       IF( env_optical_permittivity < 1.0_DP ) &
          CALL errore( sub_name,' env_optical_permittivity out of range ', 1 )
       allowed = .FALSE.
       DO i = 1, SIZE( eps_mode_allowed )
          IF( TRIM(eps_mode) == eps_mode_allowed(i) ) allowed = .TRUE.
       END DO
       IF( .NOT. allowed ) &
          CALL errore( sub_name, ' eps_mode '''// &
                       & TRIM(eps_mode)//''' not allowed ', 1 )
       IF( alpha <= 0.0_DP ) &
          CALL errore( sub_name,' alpha out of range ', 1 )
       !
       IF( env_surface_tension < 0.0_DP ) &
          CALL errore( sub_name,' env_surface_tension out of range ', 1 )
       IF( delta <= 0.0_DP ) &
          CALL errore( sub_name,' delta out of range ', 1 )
       !
       IF( env_pressure < 0.0_DP ) &
          CALL errore( sub_name,' env_pressure out of range ', 1 )
       !
       IF( env_ioncc_level < 0 .OR. env_ioncc_level > 3 ) &
          CALL errore( sub_name,' env_ioncc_level out of range ', 1 )
       IF( nrep < 0 ) &
          CALL errore( sub_name,' nrep out of range ', 1 )
       allowed = .FALSE.
       DO i = 1, SIZE( stern_mode_allowed )
          IF( TRIM(stern_mode) == stern_mode_allowed(i) ) allowed = .TRUE.
       END DO
       IF( .NOT. allowed ) &
          CALL errore( sub_name, ' stern_mode '''// &
                       & TRIM(stern_mode)//''' not allowed ', 1 )
       IF( stern_distance < 0.0_DP ) &
          CALL errore( sub_name,' stern_distance out of range ', 1 )
       IF( stern_spread <= 0.0_DP ) &
          CALL errore( sub_name,' stern_spread out of range ', 1 )
       IF( cion < 0.0_DP ) &
          CALL errore( sub_name,' cion out of range ', 1 )
       IF( cionmax < 0.0_DP ) &
          CALL errore( sub_name,' cionmax out of range ', 1 )
       IF( rion < 0.0_DP ) &
          CALL errore( sub_name,' rion out of range ', 1 )
       IF( zion < 0.0_DP ) &
          CALL errore( sub_name,' zion out of range ', 1 )
       IF( rhopb <= 0.0_DP ) &
          CALL errore( sub_name,' rhopb out of range ', 1 )
       IF( solvent_temperature < 0.0_DP ) &
          CALL errore( sub_name,' solvent_temperature out of range ', 1 )
       !
       IF( env_periodicity < 0 .OR. env_periodicity > 3 ) &
          CALL errore( sub_name,' env_periodicity out of range ', 1 )
       IF( cell_axis < 1 .OR. cell_axis > 3 ) &
          CALL errore( sub_name,' cell_axis out of range ', 1 )
       allowed = .FALSE.
       DO i = 1, SIZE( pbc_correction_allowed )
          IF( TRIM(pbc_correction) == pbc_correction_allowed(i) ) allowed = .TRUE.
       END DO
       IF( .NOT. allowed ) &
          CALL errore( sub_name, ' pbc_correction '''// &
                       & TRIM(pbc_correction)//''' not allowed ', 1 )
       !
       IF ( env_external_charges < 0 ) &
          CALL errore( sub_name,' env_external_charges out of range ', 1 )
       !
       IF ( env_dielectric_regions < 0 ) &
          CALL errore( sub_name,' env_dielectric_regions out of range ', 1 )
       !
       RETURN
       !
     END SUBROUTINE
     !
     !=----------------------------------------------------------------------=!
     !
     !  Check input values for Namelist ELECTROSTATIC
     !
     !=----------------------------------------------------------------------=!
     !
     !-----------------------------------------------------------------------
     SUBROUTINE electrostatic_checkin()
       !-----------------------------------------------------------------------
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
       IF( tolvelect <= 0.0_DP ) &
          CALL errore( sub_name,' tolvelect out of range ', 1 )
       IF( tolrhoaux <= 0.0_DP ) &
          CALL errore( sub_name,' tolrhoaux out of range ', 1 )
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
       IF( ifdtype < 1 ) &
          CALL errore( sub_name,' ifdtype out of range ', 1 )
       IF( nfdpoint < 1 ) &
          CALL errore( sub_name,' nfdpoint out of range ', 1 )
       !
       RETURN
       !
     END SUBROUTINE electrostatic_checkin
     !
     !=----------------------------------------------------------------------=!
     !
     !  Set values according to the ENVIRON namelist
     !
     !=----------------------------------------------------------------------=!
     !
     !-----------------------------------------------------------------------
     SUBROUTINE fixelect( )
       !-----------------------------------------------------------------------
       !
       IMPLICIT NONE
       !
       CHARACTER(LEN=20) :: sub_name = ' fixelect '
       !
       IF ( env_static_permittivity > 1.D0 &
            .OR. env_dielectric_regions > 0 ) &
          problem = 'generalized'
       !
       IF ( env_ioncc_level > 0 ) THEN
          problem = 'modpb'
          solver = 'lbfgs'
       END IF
       !
       IF ( env_periodicity == 0 ) THEN
          bcindex = 1
       ELSE IF ( env_periodicity == 1 ) THEN
          bcindex = 1
          bcindex( cell_axis ) = 0
       ELSE IF ( env_periodicity == 2 ) THEN
          bcindex = 0
          bcindex( cell_axis ) = 1
       ENDIF
       !
       RETURN
       !
     END SUBROUTINE
     !
     !=----------------------------------------------------------------------=!
     !
     !  Environ cards parsing routine
     !
     !=----------------------------------------------------------------------=!
     !
     !-----------------------------------------------------------------------
     SUBROUTINE environ_read_cards( unit )
       !-----------------------------------------------------------------------
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
100    CALL read_line( input_line, end_of_file=tend )
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
          IF ( ionode ) &
             WRITE( stdout,'(A)') 'Warning: card '//trim(input_line)//' ignored'
          !
       ENDIF
       !
       ! ... END OF LOOP ... !
       !
       GOTO 100
       !
120       CONTINUE
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
     END SUBROUTINE environ_read_cards
     !
   !----------------------------------------------------------------------
   !
   ! ... Description of the allowed input CARDS
   !
   ! EXTERNAL_CHARGES (unit_option origin_option)
   !
   !   set external fixed charge densities and their shape
   !
   ! Syntax:
   !
   !    EXTERNAL_CHARGES (unit_option origin_option)
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
   !   origin_option == cell     positions are given wrt cell's origin (DEFAULT)
   !   origin_option == bound    positions are given wrt system's bounding box
   !   origin_option == system   positions are given wrt system's origin
   !
   !      charge(i) ( real )       total charge of the density
   !      x/y/z(i)  ( real )       cartesian position of the density
   !      spread(i) ( real )       gaussian spread of the density (in bohr, optional, default=0.5)
   !      dim(i)    ( integer )    0/1/2 point/line/plane of charge (optional, default=0)
   !      axis(i)   ( integer )    1/2/3 for x/y/z direction of line/plane (optional, default=3)
   !
   !----------------------------------------------------------------------
   !
   SUBROUTINE card_external_charges( input_line )
      !
      USE wrappers, ONLY: feval_infix
      !
      IMPLICIT NONE
      !
      CHARACTER(len=256) :: input_line
      INTEGER            :: ie, ierr, nfield
      LOGICAL            :: tend
      LOGICAL, EXTERNAL  :: matches
      CHARACTER(len=4)   :: lb_pos
      CHARACTER(len=256) :: field_str
      !
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
      IF ( matches( "CEL", input_line ) ) THEN
         extcharge_origin = 0
      ELSEIF ( matches( "SYS", input_line ) ) THEN
         extcharge_origin = 1
      ELSEIF ( matches( "BOUND", input_line ) ) THEN
         extcharge_origin = 2
      ELSE
         CALL infomsg( 'read_cards ', &
            & 'No origin specified in EXTERNAL_CHARGES card' )
            extcharge_origin = 0
         CALL infomsg( 'read_cards ', &
            & 'EXTERNAL_CHARGES: origin set to cell origin')
      ENDIF
      !
      DO ie = 1, env_external_charges
         !
         CALL read_line( input_line, end_of_file = tend )
         IF ( tend ) CALL errore( 'environ_cards', &
              'end of file reading external charges', ie )
         !
         CALL field_count( nfield, input_line )
         !
         ! ... read field 1 (total charge of the external density)
         !
         CALL get_field(1, field_str, input_line)
         extcharge_charge(ie) = feval_infix(ierr, field_str )
         !
         ! ... read fields 2-4 (x-y-z position of external density)
         !
         CALL get_field(2, field_str, input_line)
         extcharge_pos(1,ie) = feval_infix(ierr, field_str )
         CALL get_field(3, field_str, input_line)
         extcharge_pos(2,ie) = feval_infix(ierr, field_str )
         CALL get_field(4, field_str, input_line)
         extcharge_pos(3,ie) = feval_infix(ierr, field_str )
         !
         ! ... optionally read field 5 (spread of the density)
         !
         IF ( nfield >= 5 ) THEN
           CALL get_field(5, field_str, input_line)
           extcharge_spread(ie) = feval_infix(ierr, field_str )
           IF ( extcharge_spread(ie) .LT. 0.D0 ) &
             CALL errore( ' card_external_charges  ', ' spread must be positive', ie )
         ENDIF
         !
         ! ... optionally read field 6 and 7 (dimensionality and direction)
         !
         IF ( nfield >= 6 ) THEN
           CALL get_field(6, field_str, input_line)
           READ(field_str, *) extcharge_dim(ie)
           IF ( extcharge_dim(ie) .LT. 0 .OR. extcharge_dim(ie) .GT. 2 ) &
             CALL errore( ' card_external_charges  ', ' wrong excharge dimension ', ie )
           IF ( extcharge_dim(ie) .GT. 0 ) THEN
             IF ( nfield == 6 ) &
             CALL errore('environ_cards',&
             'missing axis direction of partially periodic external charge', ie)
             CALL get_field(7, field_str, input_line)
             READ(field_str, *) extcharge_axis(ie)
             IF ( extcharge_axis(ie) .LT. 0 .OR. extcharge_axis(ie) .GT. 3 ) &
               CALL errore( ' card_external_charges  ', ' wrong excharge axis ', ie )
           ENDIF
         ENDIF
         !
      ENDDO
      taextchg = .true.
      !
      CALL convert_pos( external_charges, env_external_charges, extcharge_pos )
      !
      RETURN
      !
   END SUBROUTINE card_external_charges
   !
   !-----------------------------------------------------------------------------
   SUBROUTINE allocate_input_extcharge(env_external_charges)
   !-----------------------------------------------------------------------------
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
   END SUBROUTINE allocate_input_extcharge
   !
   !----------------------------------------------------------------------
   !
   ! ... Description of the allowed input CARDS
   !
   ! DIELECTRIC_REGIONS (unit_option origin_option)
   !
   !   set fixed dielectric regions and their shape
   !
   ! Syntax:
   !
   !    DIELECTRIC_REGIONS (unit_option origin_option)
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
   !   origin_option == cell     positions are given wrt cell's origin (DEFAULT)
   !   origin_option == bound    positions are given wrt system's bounding box
   !   origin_option == system   positions are given wrt system's origin
   !
   !      epsilon0(i)   ( real )    static permittivity inside the region
   !      epsilonopt(i) ( real )    optical permittivity inside the region
   !      x/y/z(i)      ( real )    cartesian center of the region
   !      width(i)      ( real )    size of the region (in bohr)
   !      spread(i)     ( real )    spread of the interface (in bohr, optional)
   !      dim(i)     ( integer )    0/1/2 point/line/plane region (optional)
   !      axis(i)    ( integer )    1/2/3 for x/y/z direction of line/plane (optional)
   !
   !----------------------------------------------------------------------
   !
   SUBROUTINE card_dielectric_regions( input_line )
      !
      USE wrappers, ONLY: feval_infix
      !
      IMPLICIT NONE
      !
      CHARACTER(len=256) :: input_line
      INTEGER            :: ie, ierr, nfield
      LOGICAL            :: tend
      LOGICAL, EXTERNAL  :: matches
      CHARACTER(len=4)   :: lb_pos
      CHARACTER(len=256) :: field_str
      !
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
      IF ( matches( "CEL", input_line ) ) THEN
         epsregion_origin = 0
      ELSEIF ( matches( "SYS", input_line ) ) THEN
         epsregion_origin = 1
      ELSEIF ( matches( "BOUND", input_line ) ) THEN
         epsregion_origin = 2
      ELSE
         CALL infomsg( 'read_cards ', &
            & 'No origin specified in DIELECTRIC_REGIONS card' )
         epsregion_origin = 0
         CALL infomsg( 'read_cards ', &
            & 'DIELECTRIC_REGIONS: origin set to cell origin' )
      ENDIF
      !
      DO ie = 1, env_dielectric_regions
         !
         CALL read_line( input_line, end_of_file = tend )
         IF ( tend ) CALL errore( 'environ_cards', &
              'end of file reading dielectric regions', ie )
         !
         CALL field_count( nfield, input_line )
         !
         ! ... read field 1-2 (static and optical permettivity inside dielectric region)
         !
         CALL get_field(1, field_str, input_line)
         epsregion_eps(1,ie) = feval_infix(ierr, field_str )
         IF ( epsregion_eps(1,ie) .LT. 1.D0 ) &
           CALL errore( ' card_dielectric_regions  ', ' static permittivity must be .gt. 1', ie )
         CALL get_field(2, field_str, input_line)
         epsregion_eps(2,ie) = feval_infix(ierr, field_str )
         IF ( epsregion_eps(2,ie) .LT. 1.D0 ) &
           CALL errore( ' card_dielectric_regions  ', ' optical permittivity must be .gt. 1', ie )
         !
         ! ... read fields 3-5 (x-y-z position of dielectric region)
         !
         CALL get_field(3, field_str, input_line)
         epsregion_pos(1,ie) = feval_infix(ierr, field_str )
         CALL get_field(4, field_str, input_line)
         epsregion_pos(2,ie) = feval_infix(ierr, field_str )
         CALL get_field(5, field_str, input_line)
         epsregion_pos(3,ie) = feval_infix(ierr, field_str )
         !
         ! ... read field 6 (size/width of the dielectric region)
         !
         CALL get_field(6, field_str, input_line)
         epsregion_width(ie) = feval_infix(ierr, field_str )
         IF ( epsregion_width(ie) .LT. 0.D0 ) &
           CALL errore( ' card_dielectric_regions  ', ' width must be positive', ie )
         !
         ! ... optionally read field 7 (spread of interface of the dielectric region)
         !
         IF ( nfield >= 7 ) THEN
           CALL get_field(7, field_str, input_line)
           epsregion_spread(ie) = feval_infix(ierr, field_str )
           IF ( epsregion_spread(ie) .LT. 0.D0 ) &
             CALL errore( ' card_dielectric_regions ', ' spread must be positive', ie )
         ENDIF
         !
         ! ... optionally read field 7 and 8 (dimensionality and direction)
         !
         IF ( nfield >= 8 ) THEN
           CALL get_field(8, field_str, input_line)
           READ(field_str, *) epsregion_dim(ie)
           IF ( epsregion_dim(ie) .LT. 0 .OR. epsregion_dim(ie) .GT. 2 ) &
             CALL errore( ' card_dielectric_regions ', ' wrong epsregion dimension ', ie )
           IF ( epsregion_dim(ie) .GT. 0 ) THEN
             IF ( nfield == 8 ) &
             CALL errore('environ_cards',&
             'missing axis direction of partially periodic dielectric region', ie)
             CALL get_field(9, field_str, input_line)
             READ(field_str, *) epsregion_axis(ie)
             IF ( epsregion_axis(ie) .LT. 1 .OR. epsregion_axis(ie) .GT. 3 ) &
               CALL errore( ' card_dielectric_regions ', ' wrong epsregion axis ', ie )
           ENDIF
         ENDIF
         !
      ENDDO
      taepsreg = .true.
      !
      CALL convert_pos( dielectric_regions, env_dielectric_regions, epsregion_pos )
      !
      RETURN
      !
   END SUBROUTINE card_dielectric_regions
   !
   !-----------------------------------------------------------------------------
   SUBROUTINE allocate_input_epsregion(env_dielectric_regions)
   !-----------------------------------------------------------------------------
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
   END SUBROUTINE allocate_input_epsregion
   !
   !-----------------------------------------------------------------------
   SUBROUTINE convert_pos (pos_format, n, pos)
   !-----------------------------------------------------------------------
     !
     ! ... convert input positions to atomic units
     !
     USE kinds,         ONLY : DP
     USE constants,     ONLY : bohr_radius_angs
     IMPLICIT NONE
     CHARACTER (len=*), INTENT(in)  :: pos_format
     INTEGER, INTENT(in)  :: n
     REAL (DP), INTENT(inout) :: pos(3,n)
     !
     SELECT CASE( pos_format )
     CASE( 'bohr' )
        !
        ! ... input positions are in a.u., do nothing
        !
        pos = pos
        !
     CASE( 'angstrom' )
        !
        ! ... positions in A: convert to a.u.
        !
        pos = pos / bohr_radius_angs
        !
     CASE DEFAULT
        !
        CALL errore( 'iosys','pos_format=' // &
                   & trim( pos_format ) // ' not implemented', 1 )
        !
     END SELECT
     !
   END SUBROUTINE convert_pos
!=----------------------------------------------------------------------------=!
!
END MODULE environ_input
!
!=----------------------------------------------------------------------------=!
