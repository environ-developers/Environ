MODULE environ_types

  USE kinds,             ONLY : DP
  USE constants,         ONLY : rydberg_si, bohr_radius_si, amu_si, fpi, tpi, sqrtpi
  USE mp,                ONLY : mp_sum
  USE control_flags,     ONLY : tddfpt

  TYPE environ_cell

     ! Global properties of the simulation cell

     LOGICAL :: update = .FALSE.
     INTEGER :: ibrav
     INTEGER :: ntot, n1, n2, n3
     REAL( DP ) :: alat
     REAL( DP ) :: omega
     REAL( DP ) :: domega
     REAL( DP ), DIMENSION( 3, 3 ) :: at

     ! Properties of the processor-specific partition

     INTEGER :: nnr    ! size of processor-specific allocated fields
     INTEGER :: ir_end ! actual physical size of processor-specific allocated field
     INTEGER :: comm   ! parallel communicator
     INTEGER :: me     ! index of processor
     INTEGER :: root   ! index of root

  END TYPE environ_cell

  TYPE environ_density

     ! Optionally have an associated logical status

     LOGICAL :: update = .FALSE.

     ! Optionally have an associated label, used for printout and debugs

     CHARACTER( LEN=80 ) :: label = ' '

     ! Each quantity in real-space is associated with its definition domain

     TYPE( environ_cell ), POINTER :: cell => NULL()

     ! The quantity in real-space, local to each processor

     REAL( DP ), DIMENSION(:), ALLOCATABLE :: of_r

  END TYPE environ_density

  TYPE environ_gradient

     ! Optionally have an associated logical status

     LOGICAL :: update = .FALSE.

     ! Optionally have an associated label, used for printout and debugs

     CHARACTER( LEN=80 ) :: label = ' '

     ! Each quantity in real-space is associated with its definition domain

     TYPE( environ_cell ), POINTER :: cell => NULL()

     ! The quantity in real-space, local to each processor

     REAL( DP ), DIMENSION(:,:), ALLOCATABLE :: of_r

     TYPE( environ_density ) :: modulus

  END TYPE environ_gradient

  TYPE environ_hessian

     ! Optionally have an associated logical status

     LOGICAL :: update = .FALSE.

     ! Optionally have an associated label, used for printout and debugs

     CHARACTER( LEN=80 ) :: label = ' '

     ! Each quantity in real-space is associated with its definition domain

     TYPE( environ_cell ), POINTER :: cell => NULL()

     ! The quantity in real-space, local to each processor

     REAL( DP ), DIMENSION(:,:,:), ALLOCATABLE :: of_r

     REAL( DP ), DIMENSION(:), ALLOCATABLE :: laplacian

  END TYPE environ_hessian

  TYPE environ_functions

     INTEGER :: type
     INTEGER :: axis, dim
     REAL( DP ) :: width, spread, volume
     REAL( DP ), DIMENSION(:), POINTER :: pos
     ! environ_functions are not designed to be mobile,
     ! thus position can be included in the definition
     ! of the type

  END TYPE environ_functions

  TYPE environ_iontype

     INTEGER :: index
     INTEGER :: atmnum
     CHARACTER( LEN=3 ) :: label
     REAL( DP ) :: zv
     REAL( DP ) :: atomicspread
     REAL( DP ) :: corespread
     REAL( DP ) :: solvationrad

  END TYPE environ_iontype

  TYPE environ_ions

     LOGICAL :: initialized = .FALSE.
     LOGICAL :: update = .FALSE.
     INTEGER :: number = 0
     REAL( DP ), DIMENSION(3) :: center

     ! Specifications of point-like ions

     INTEGER :: ntyp = 0
     INTEGER, DIMENSION(:), ALLOCATABLE :: ityp
     REAL( DP ), DIMENSION(:,:), POINTER :: tau
     TYPE( environ_iontype ), DIMENSION(:), ALLOCATABLE :: iontype

     ! Parameters of the fictitious gaussian ionic density
     ! needed by electrostatic calculations

     LOGICAL :: use_smeared_ions = .FALSE.
     TYPE( environ_functions ), DIMENSION(:), ALLOCATABLE :: smeared_ions
     TYPE( environ_density ) :: density

     ! Parameters of the density of core electrons

     LOGICAL :: use_core_electrons = .FALSE.
     TYPE( environ_functions ), DIMENSION(:), ALLOCATABLE :: core_electrons
     TYPE( environ_density ) :: core

     REAL( DP ) :: charge = 0.0_DP

  END TYPE environ_ions

  TYPE environ_electrons

     LOGICAL :: update = .FALSE.
     INTEGER :: number = 0
     INTEGER :: nspin = 1

     REAL( DP ) :: nelec

     TYPE( environ_density ) :: density

     REAL( DP ) :: charge = 0.0_DP

  END TYPE environ_electrons

  TYPE environ_externals

     LOGICAL :: update = .FALSE.
     INTEGER :: number = 0

     TYPE( environ_functions ), DIMENSION(:), ALLOCATABLE :: functions

     TYPE( environ_density ) :: density

     REAL( DP ) :: charge = 0.0_DP

  END TYPE environ_externals

  TYPE environ_auxiliary

     LOGICAL :: update = .FALSE.
     INTEGER :: number = 0

     TYPE( environ_density ) :: density

     TYPE( environ_density ) :: fixed
     TYPE( environ_density ) :: iterative

     REAL( DP ) :: charge = 0.0_DP

  END TYPE environ_auxiliary

  TYPE environ_charges

     ! Ionic charges

     LOGICAL :: include_ions = .FALSE.
     TYPE( environ_ions ), POINTER :: ions => NULL()

     ! Electrons

     LOGICAL :: include_electrons = .FALSE.
     TYPE( environ_electrons ), POINTER :: electrons => NULL()

     ! External charges

     LOGICAL :: include_externals = .FALSE.
     TYPE( environ_externals ), POINTER :: externals => NULL()

     ! Auxiliary charges

     LOGICAL :: include_auxiliary = .FALSE.
     TYPE( environ_auxiliary ), POINTER :: auxiliary => NULL()

     ! Total smooth free charge

     INTEGER :: number = 0
     REAL( DP ) :: charge = 0.0_DP
     TYPE( environ_density ) :: density

  END TYPE environ_charges

  TYPE environ_system

     INTEGER :: ntyp
     INTEGER :: dim
     INTEGER :: axis

     REAL( DP ) :: pos(3)
     REAL( DP ) :: width

     TYPE( environ_ions ), POINTER :: ions

  END TYPE environ_system

  TYPE environ_boundary

     ! Choice of the interface

     CHARACTER (LEN=80) :: mode

     ! Update status

     INTEGER :: update_status = 0

     ! Parameters for the electrons-dependent interface

     LOGICAL :: need_electrons
     TYPE( environ_electrons ), POINTER :: electrons

     ! Parameters for the ions-dependent interface

     LOGICAL :: need_ions
     TYPE( environ_ions ), POINTER :: ions

     ! scaled switching function of interface
     ! varying from 1 (QM region) to 0 (environment region)

     TYPE( environ_density ) :: scaled

     INTEGER :: deriv = 0
     TYPE( environ_gradient ) :: gradient
     TYPE( environ_density ) :: laplacian
     TYPE( environ_density ) :: dsurface

     ! Components needed for boundary of density

     INTEGER :: type
     REAL( DP ) :: rhomax, rhomin, fact
     REAL( DP ) :: rhozero, deltarho, tbeta
     REAL( DP ) :: const
     TYPE( environ_density ) :: density

     TYPE( environ_density ) :: dscaled
     TYPE( environ_density ) :: d2scaled

     ! Components needed for boundary of functions

     REAL( DP ) :: alpha ! solvent-dependent scaling factor
     REAL( DP ) :: softness ! sharpness of the interface
     TYPE( environ_functions ), DIMENSION(:), ALLOCATABLE :: soft_spheres

  END TYPE environ_boundary

  TYPE environ_dielectric

     ! Update status

     LOGICAL :: update = .FALSE.

     ! Basic properties of the dielectric space from input

     INTEGER :: nregions
     TYPE( environ_functions ), DIMENSION(:), ALLOCATABLE :: regions

     REAL( DP ) :: constant
     TYPE( environ_density ) :: background

     ! Boundary is the pointer to the object controlling
     ! the interface between the QM and the continuum region

     TYPE( environ_boundary ), POINTER :: boundary

     ! The dielectric function over space is built from the
     ! boundary of the continuum environment and the basic dielectric
     ! properties of space

     TYPE( environ_density ) :: epsilon
     TYPE( environ_density ) :: depsilon

     ! Quantities related to the dielectric permittivity and
     ! thay may be needed by the different solvers

     LOGICAL :: need_gradient = .FALSE.
     TYPE( environ_gradient ) :: gradient

     LOGICAL :: need_factsqrt = .FALSE.
     TYPE( environ_density ) :: factsqrt

     LOGICAL :: need_gradlog = .FALSE.
     TYPE( environ_gradient ) :: gradlog

  END TYPE environ_dielectric

  TYPE environ_ioncctype

     INTEGER :: index
     REAL( DP ) :: cbulk   ! bulk concentration
     REAL( DP ) :: radius  ! radius
     REAL( DP ) :: cmax    ! maximum allowed concentration
     REAL( DP ) :: z       ! charge
     REAL( DP ) :: mu      ! chemical potential
     REAL( DP ) :: epsilon ! dielectric constant

  END TYPE environ_ioncctype

  TYPE environ_electrolyte

     INTEGER :: ntyp
     TYPE( environ_ioncctype ), DIMENSION(:), ALLOCATABLE :: ioncctype

     REAL( DP ) :: temperature

     TYPE( environ_boundary ) :: boundary

     TYPE( environ_density ) :: density

  END TYPE environ_electrolyte

CONTAINS
!----------------------------------------------------------------------------------------------------------------------------------------
!- CELL ---------------------------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------------------------------
  SUBROUTINE init_environ_cell( n1, n2, n3, ibrav, alat, omega, at, nnr, ir_end, comm, me, root, cell )

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: n1, n2, n3, ibrav
    INTEGER, INTENT(IN) :: nnr, ir_end, comm, me, root
    REAL( DP ), INTENT(IN) :: alat, omega, at(3,3)
    TYPE( environ_cell ), INTENT(INOUT) :: cell
    CHARACTER( LEN=80 ) :: sub_name = 'init_environ_cell'

    cell % n1 = n1
    cell % n2 = n2
    cell % n3 = n3
    cell % ibrav = ibrav
    cell % alat = alat
    cell % omega = omega
    cell % at = at

    cell % nnr = nnr
    cell % ir_end = ir_end
    cell % comm = comm
    cell % me   = me
    cell % root = root

    cell % ntot = cell % n1 * cell % n2 * cell % n3
    cell % domega = cell % omega / cell % ntot

    RETURN

  END SUBROUTINE init_environ_cell

  SUBROUTINE update_environ_cell( omega, at, cell )

    IMPLICIT NONE

    REAL( DP ), INTENT(IN) :: omega, at(3,3)
    TYPE( environ_cell ), INTENT(INOUT) :: cell
    CHARACTER( LEN=80 ) :: sub_name = 'update_environ_cell'

    cell % omega = omega
    cell % at = at

    cell % domega = cell % omega / cell % ntot

    RETURN

  END SUBROUTINE update_environ_cell
!----------------------------------------------------------------------------------------------------------------------------------------
!- DENSITY ------------------------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------------------------------
  SUBROUTINE create_environ_density(density,local_label)

    IMPLICIT NONE

    TYPE( environ_density ), INTENT(INOUT) :: density
    CHARACTER( LEN=80 ), INTENT(IN), OPTIONAL :: local_label
    CHARACTER( LEN=80 ) :: sub_name = 'create_environ_density'

    CHARACTER ( LEN=80 ) :: label = 'density'
    IF ( PRESENT(local_label) ) THEN
       density%label = local_label
    ELSE
       density%label = label
    END IF

    NULLIFY(density%cell)
    IF ( ALLOCATED( density%of_r ) ) CALL errore(sub_name,'Trying to create an already allocated object',1)

    RETURN

  END SUBROUTINE create_environ_density

  SUBROUTINE init_environ_density( cell, density )

    IMPLICIT NONE

    TYPE( environ_cell ), TARGET, INTENT(IN) :: cell
    TYPE( environ_density ), INTENT(INOUT) :: density
    CHARACTER( LEN=80 ) :: sub_name = 'init_environ_density'

    density%update = .FALSE.

    IF ( ASSOCIATED( density%cell ) ) CALL errore(sub_name,'Trying to associate an associated object',1)
    density%cell => cell

    IF ( ALLOCATED( density%of_r ) ) CALL errore(sub_name,'Trying to allocate an allocated object',1)
    ALLOCATE(density%of_r(density%cell%nnr))
    density%of_r = 0.D0

    RETURN

  END SUBROUTINE init_environ_density

  FUNCTION integrate_environ_density(density) RESULT(integral)

    IMPLICIT NONE

    TYPE( environ_density ), INTENT(IN) :: density

    REAL( DP ) :: integral

    integral = SUM(density%of_r(1:density%cell%ir_end))
    CALL mp_sum( integral, density%cell%comm )
    integral = integral * density%cell%domega

    RETURN

  END FUNCTION integrate_environ_density

  FUNCTION scalar_product_environ_density(density1, density2) RESULT(scalar_product)

    IMPLICIT NONE

    TYPE( environ_density ), INTENT(IN) :: density1, density2

    INTEGER, POINTER :: ir_end
    REAL( DP ) :: scalar_product
    CHARACTER( LEN=80 ) :: fun_name = 'scalar_product_environ_density'

    IF ( .NOT.ASSOCIATED(density1%cell,density2%cell) ) &
         & CALL errore(fun_name,'operation on fields with inconsistent domains',1)
    ir_end => density1 % cell % ir_end
    scalar_product = DOT_PRODUCT(density1%of_r(1:ir_end),density2%of_r(1:ir_end))
    CALL mp_sum( scalar_product, density1%cell%comm )
    scalar_product = scalar_product * density1%cell%domega

    RETURN

  END FUNCTION scalar_product_environ_density

  FUNCTION quadratic_mean_environ_density(density) RESULT(quadratic_mean)

    IMPLICIT NONE

    TYPE( environ_density ), INTENT(IN) :: density

    REAL( DP ) :: quadratic_mean

    quadratic_mean = DOT_PRODUCT(density%of_r,density%of_r)
    CALL mp_sum( quadratic_mean, density%cell%comm )
    quadratic_mean = SQRT( quadratic_mean ) / density % cell % ntot

    RETURN

  END FUNCTION quadratic_mean_environ_density

  SUBROUTINE destroy_environ_density(density)

    IMPLICIT NONE

    TYPE( environ_density ), INTENT(INOUT) :: density
    CHARACTER (LEN=80) :: sub_name = 'destroy_environ_density'

    density%update = .FALSE.

    IF (.NOT.ASSOCIATED(density%cell)) &
         & CALL errore(sub_name,'Trying to destroy a non associated object',1)
    NULLIFY(density%cell)

    IF (.NOT.ALLOCATED(density%of_r)) &
         & CALL errore(sub_name,'Trying to destroy a non allocated object',1)
    DEALLOCATE( density%of_r )

    RETURN

  END SUBROUTINE destroy_environ_density
!----------------------------------------------------------------------------------------------------------------------------------------
!- GRADIENT -----------------------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------------------------------
  SUBROUTINE create_environ_gradient(gradient,label)

    IMPLICIT NONE

    TYPE( environ_gradient ), INTENT(INOUT) :: gradient
    CHARACTER (LEN=80), INTENT(IN), OPTIONAL :: label

    CHARACTER (LEN=80) :: sub_name = 'destroy_environ_density'
    CHARACTER (LEN=80) :: modulus_label

    IF ( PRESENT(label) ) THEN
       gradient%label = label
       modulus_label = TRIM(ADJUSTL(label))//'_modulus'
    ELSE
       gradient%label = 'gradient'
       modulus_label = 'gradient_modulus'
    END IF

    NULLIFY( gradient%cell )
    IF ( ALLOCATED( gradient%of_r ) ) CALL errore(sub_name,'Trying to create an already allocated object',1)
    CALL create_environ_density( gradient%modulus, modulus_label )

    RETURN

  END SUBROUTINE create_environ_gradient

  SUBROUTINE init_environ_gradient( cell, gradient )

    IMPLICIT NONE

    TYPE( environ_cell ), TARGET, INTENT(IN) :: cell
    TYPE( environ_gradient ), INTENT(INOUT) :: gradient
    CHARACTER( LEN=80 ) :: sub_name = 'init_environ_gradient'

    gradient%update = .FALSE.

    IF ( ASSOCIATED( gradient%cell ) ) CALL errore(sub_name,'Trying to associate an associated object',1)
    gradient%cell => cell

    IF ( ALLOCATED( gradient%of_r ) ) CALL errore(sub_name,'Trying to allocate an allocated object',1)
    ALLOCATE(gradient%of_r(3,gradient%cell%nnr))
    gradient%of_r = 0.D0

    CALL init_environ_density( cell, gradient%modulus )

    RETURN

  END SUBROUTINE init_environ_gradient

  SUBROUTINE update_gradient_modulus(gradient)

    IMPLICIT NONE

    TYPE( environ_gradient ), INTENT(INOUT) :: gradient
    INTEGER, POINTER :: ir_end

    ir_end => gradient % cell % ir_end
    gradient%modulus%of_r(1:ir_end) = gradient%of_r(1,1:ir_end)**2 + &
         & gradient%of_r(2,1:ir_end)**2 + gradient%of_r(3,1:ir_end)**2

    RETURN

  END SUBROUTINE update_gradient_modulus

  SUBROUTINE destroy_environ_gradient(gradient)

    IMPLICIT NONE

    TYPE( environ_gradient ), INTENT(INOUT) :: gradient
    CHARACTER (LEN=80) :: sub_name = 'destroy_environ_gradient'

    gradient%update = .FALSE.

    IF (.NOT.ASSOCIATED(gradient%cell)) &
         & CALL errore(sub_name,'Trying to destroy a non associated object',1)
    NULLIFY(gradient%cell)

    IF (.NOT.ALLOCATED(gradient%of_r)) &
         & CALL errore(sub_name,'Trying to destroy a non allocated object',1)
    DEALLOCATE( gradient%of_r )

    CALL destroy_environ_density( gradient%modulus )

    RETURN

  END SUBROUTINE destroy_environ_gradient

  SUBROUTINE scalar_product_environ_gradient( gradA, gradB, dens )

    IMPLICIT NONE

    TYPE( environ_gradient ), INTENT(IN) :: gradA, gradB
    TYPE( environ_density ), INTENT(INOUT) :: dens

    INTEGER :: ir
    CHARACTER( LEN=80 ) :: sub_name = 'scalar_product_environ_gradient'

    dens%of_r = 0.D0
    IF ( .NOT. ASSOCIATED(gradA%cell,gradB%cell) ) &
         & CALL errore(sub_name,'Missmatch in domain of input gradients',1)
    IF ( .NOT. ASSOCIATED(gradA%cell,dens%cell) ) &
         & CALL errore(sub_name,'Missmatch in domain of input and output',1)

    DO ir = 1, dens%cell%ir_end
       dens%of_r(ir) = SUM(gradA%of_r(:,ir)*gradB%of_r(:,ir))
    END DO

    RETURN

  END SUBROUTINE scalar_product_environ_gradient
!----------------------------------------------------------------------------------------------------------------------------------------
!- HESSIAN ------------------------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------------------------------
  SUBROUTINE create_environ_hessian(hessian,label)

    IMPLICIT NONE

    TYPE( environ_hessian ), INTENT(INOUT) :: hessian
    CHARACTER (LEN=80), OPTIONAL, INTENT(IN) :: label
    CHARACTER (LEN=80) :: sub_name = 'destroy_environ_hessian'

    IF ( PRESENT(label) ) THEN
       hessian%label = label
    ELSE
       hessian%label = 'hessian'
    END IF

    NULLIFY(hessian%cell)
    IF ( ALLOCATED( hessian%of_r ) ) CALL errore(sub_name,'Trying to create an already allocated object',1)
    IF ( ALLOCATED( hessian%laplacian ) ) CALL errore(sub_name,'Trying to create an already allocated object',1)

    RETURN

  END SUBROUTINE create_environ_hessian

  SUBROUTINE init_environ_hessian( cell, hessian )

    IMPLICIT NONE

    TYPE( environ_cell ), TARGET, INTENT(IN) :: cell
    TYPE( environ_hessian ), INTENT(INOUT) :: hessian
    CHARACTER( LEN=80 ) :: sub_name = 'init_environ_hessian'

    hessian%update = .FALSE.

    IF ( ASSOCIATED( hessian%cell ) ) CALL errore(sub_name,'Trying to associate an associated object',1)
    hessian%cell => cell

    IF ( ALLOCATED( hessian%of_r ) ) CALL errore(sub_name,'Trying to allocate an allocated object',1)
    ALLOCATE(hessian%of_r(3,3,hessian%cell%nnr))
    hessian%of_r = 0.D0

    IF ( ALLOCATED( hessian%laplacian ) ) CALL errore(sub_name,'Trying to allocate an allocated object',1)
    ALLOCATE(hessian%laplacian(hessian%cell%nnr))
    hessian%laplacian = 0.D0

    RETURN

  END SUBROUTINE init_environ_hessian

  SUBROUTINE update_hessian_laplacian(hessian)

    IMPLICIT NONE

    TYPE( environ_hessian ), INTENT(INOUT) :: hessian
    INTEGER, POINTER :: ir_end

    ir_end => hessian % cell % ir_end
    hessian%laplacian(1:ir_end) = hessian%of_r(1,1,1:ir_end)**2 + &
         & hessian%of_r(2,2,1:ir_end)**2 + hessian%of_r(3,3,1:ir_end)**2 

    RETURN

  END SUBROUTINE update_hessian_laplacian

  SUBROUTINE destroy_environ_hessian(hessian)

    IMPLICIT NONE

    TYPE( environ_hessian ), INTENT(INOUT) :: hessian
    CHARACTER (LEN=80) :: sub_name = 'destroy_environ_hessian'

    hessian%update = .FALSE.

    IF (.NOT.ASSOCIATED(hessian%cell)) &
         & CALL errore(sub_name,'Trying to destroy a non associated object',1)
    NULLIFY(hessian%cell)

    IF (.NOT.ALLOCATED(hessian%of_r)) &
         & CALL errore(sub_name,'Trying to destroy a non allocated object',1)
    DEALLOCATE( hessian%of_r )

    IF (.NOT.ALLOCATED(hessian%laplacian)) &
         & CALL errore(sub_name,'Trying to destroy a non allocated object',1)
    DEALLOCATE( hessian%laplacian )

    RETURN

  END SUBROUTINE destroy_environ_hessian
!----------------------------------------------------------------------------------------------------------------------------------------
!- ELECTRONS ----------------------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------------------------------
  SUBROUTINE create_environ_electrons(electrons)

    IMPLICIT NONE

    TYPE( environ_electrons ), INTENT(INOUT) :: electrons
    CHARACTER( LEN = 80 ) :: label = 'electrons'

    electrons%update = .FALSE.
    electrons%number = 0
    electrons%nspin  = 1
    electrons%charge = 0.D0
    CALL create_environ_density( electrons%density,label )

    RETURN

  END SUBROUTINE create_environ_electrons

  SUBROUTINE init_environ_electrons_first( nelec, nspin, electrons )

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nelec, nspin
    TYPE( environ_electrons ), INTENT(INOUT) :: electrons

    electrons%number = nelec
    electrons%nspin = nspin

    RETURN

  END SUBROUTINE init_environ_electrons_first

  SUBROUTINE init_environ_electrons_second( cell, electrons )

    IMPLICIT NONE

    TYPE( environ_cell ), INTENT(IN) :: cell
    TYPE( environ_electrons ), INTENT(INOUT) :: electrons

    CALL init_environ_density( cell, electrons%density )

    RETURN

  END SUBROUTINE init_environ_electrons_second

  SUBROUTINE update_environ_electrons( nelec, nspin, nnr, rho, electrons )

    IMPLICIT NONE

    REAL( DP ), INTENT(IN) :: nelec
    INTEGER, INTENT(IN) :: nspin, nnr
    REAL( DP ), DIMENSION(nnr,nspin), INTENT(IN) :: rho
    TYPE( environ_electrons ), INTENT(INOUT) :: electrons

    REAL( DP ), PARAMETER :: tol = 1.D-8
    REAL( DP ) :: charge
    CHARACTER( LEN= 80 ) :: sub_name = 'update_environ_electrons'

    ! Check on dimensions

    IF ( nspin .NE. electrons%nspin ) CALL errore(sub_name,'Missmatch in spin size',1)

    IF ( nnr .NE. electrons%density%cell%nnr ) CALL errore(sub_name,'Missmatch in grid size',1)

    ! Assign input density to electrons%density%of_r

    electrons%density%of_r(:) = rho(:,1)
    IF ( electrons%nspin .EQ. 2 ) electrons%density%of_r(:) = electrons%density%of_r(:) + rho(:,2)

    ! Update and check integral of electronic density

    electrons%charge = integrate_environ_density( electrons%density )
!    electrons%nelec = nelec
    electrons%nelec = electrons%charge
    electrons%number = INT(nelec)
    IF ( ABS(electrons%charge-electrons%nelec) .GT. tol ) &
             & CALL errore(sub_name,'Missmatch in integrated electronic charge',1)

    RETURN

  END SUBROUTINE update_environ_electrons

  SUBROUTINE destroy_environ_electrons( lflag, electrons )

    IMPLICIT NONE

    LOGICAL, INTENT(IN) :: lflag
    TYPE( environ_electrons ), INTENT(INOUT) :: electrons

    CALL destroy_environ_density( electrons%density )

    RETURN

  END SUBROUTINE destroy_environ_electrons
!----------------------------------------------------------------------------------------------------------------------------------------
!- SYSTEM -------------------------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------------------------------
  SUBROUTINE create_environ_system(system)

    IMPLICIT NONE

    TYPE( environ_system ), INTENT(INOUT) :: system

    NULLIFY( system%ions )

    RETURN

  END SUBROUTINE create_environ_system

  SUBROUTINE init_environ_system( ntyp, dim, axis, ions, system)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: ntyp, dim, axis
    TYPE( environ_ions ), TARGET, INTENT(IN) :: ions
    TYPE( environ_system ), INTENT(INOUT) :: system
    CHARACTER (LEN=80) :: sub_name = 'init_environ_system'

    system%ntyp = ntyp
    system%dim = dim
    system%axis = axis

    system%ions => ions

    RETURN

  END SUBROUTINE init_environ_system

  SUBROUTINE update_environ_system( system )

    ! Given the system definition compute position (center of charge)
    ! and width (maximum distance from center) of the system

    IMPLICIT NONE

    TYPE( environ_system ), INTENT(INOUT) :: system
    CHARACTER (LEN=80) :: sub_name = 'update_environ_system'

    INTEGER :: i, icor, max_ntyp
    REAL( DP ) :: charge, dist
    INTEGER, POINTER :: ityp
    REAL( DP ), POINTER :: zv

    IF (.NOT.ASSOCIATED(system%ions)) CALL errore(sub_name,'Trying to use a non associated object',1)

    system%pos = 0.D0
    system%width = 0.D0

    max_ntyp = system%ntyp
    IF ( system%ntyp .EQ. 0 ) max_ntyp = system%ions%ntyp

    charge = 0.D0
    DO i = 1, system%ions%number
       ityp => system%ions%ityp(i)
       IF ( ityp .GT. max_ntyp ) CYCLE
       zv => system%ions%iontype(ityp)%zv
       charge = charge + zv
       system%pos(:) = system%pos(:) + system%ions%tau(:,i) * zv
    ENDDO
    system%pos(:) = system%pos(:) / charge

    system%width = 0.D0
    DO i = 1, system%ions%number
       ityp => system%ions%ityp(i)
       IF ( ityp .GT. max_ntyp ) CYCLE
       DO icor = 1, 3
          IF ( ( system%dim .EQ. 1 .AND. icor .EQ. system%axis ) &
               .OR. ( system%dim .EQ. 2 .AND. icor .NE. system%axis ) ) CYCLE
          dist = dist + (system%ions%tau(icor,i)-system%pos(icor))**2
       ENDDO
       system%width = MAX(system%width,dist) ! need to modify it into a smooth maximum to compute derivatives
    ENDDO
    system%width = SQRT(system%width)

    RETURN

  END SUBROUTINE update_environ_system

  SUBROUTINE destroy_environ_system(lflag,system)

    IMPLICIT NONE

    LOGICAL, INTENT(IN) :: lflag
    TYPE( environ_system ), INTENT(INOUT) :: system
    CHARACTER (LEN=80) :: sub_name = 'destroy_environ_system'

    IF ( lflag ) THEN
       IF (.NOT.ASSOCIATED(system%ions)) &
            & CALL errore(sub_name,'Trying to destroy a non associated object',1)
       NULLIFY( system%ions )
    END IF

    RETURN

  END SUBROUTINE destroy_environ_system
!----------------------------------------------------------------------------------------------------------------------------------------
!- AUXILIARY ----------------------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------------------------------
    SUBROUTINE create_environ_auxiliary(auxiliary)

    IMPLICIT NONE

    TYPE( environ_auxiliary ), INTENT(INOUT) :: auxiliary
    CHARACTER (LEN=80) :: sub_name = 'create_environ_auxiliary'
    CHARACTER( LEN=80 ) :: label

    auxiliary%update = .FALSE.
    auxiliary%number = 0
    label = 'auxiliary'
    CALL create_environ_density( auxiliary%density, label )
    label = 'aux-fixed'
    CALL create_environ_density( auxiliary%fixed, label )
    label = 'aux-iterative'
    CALL create_environ_density( auxiliary%iterative, label )
    auxiliary%charge = 0.D0

    RETURN

  END SUBROUTINE create_environ_auxiliary

  SUBROUTINE init_environ_auxiliary( cell, auxiliary )

    IMPLICIT NONE

    TYPE( environ_cell ), INTENT(IN) :: cell
    TYPE( environ_auxiliary ), INTENT(INOUT) :: auxiliary

    INTEGER :: i

    auxiliary%number = 1
    CALL init_environ_density( cell, auxiliary%density )
    CALL init_environ_density( cell, auxiliary%fixed )
    CALL init_environ_density( cell, auxiliary%iterative )

    RETURN

  END SUBROUTINE init_environ_auxiliary

  SUBROUTINE update_environ_auxiliary( auxiliary )

    IMPLICIT NONE

    TYPE( environ_auxiliary ), INTENT(INOUT) :: auxiliary

    auxiliary%density%of_r = auxiliary%fixed%of_r + auxiliary%iterative%of_r
    auxiliary%charge = integrate_environ_density(auxiliary%density)

    RETURN

  END SUBROUTINE update_environ_auxiliary

  SUBROUTINE destroy_environ_auxiliary( auxiliary )

    IMPLICIT NONE

    TYPE( environ_auxiliary ), INTENT(INOUT) :: auxiliary

    CALL destroy_environ_density( auxiliary%density )
    CALL destroy_environ_density( auxiliary%fixed )
    CALL destroy_environ_density( auxiliary%iterative )

    RETURN

  END SUBROUTINE destroy_environ_auxiliary
!----------------------------------------------------------------------------------------------------------------------------------------
!- CHARGES ------------------------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------------------------------
  SUBROUTINE create_environ_charges(charges)

    IMPLICIT NONE

    TYPE( environ_charges ) :: charges
    CHARACTER( LEN = 80 ) :: label = 'charges'

    charges%include_ions = .FALSE.
    NULLIFY( charges%ions )

    charges%include_electrons = .FALSE.
    NULLIFY( charges%electrons )

    charges%include_externals = .FALSE.
    NULLIFY( charges%externals )

    charges%include_auxiliary = .FALSE.
    NULLIFY( charges%auxiliary )

    charges%number = 0
    charges%charge = 0.D0
    CALL create_environ_density( charges%density, label )

    RETURN

  END SUBROUTINE create_environ_charges

  SUBROUTINE init_environ_charges_first( charges, electrons, ions, externals, auxiliary )

    IMPLICIT NONE

    TYPE( environ_charges ), INTENT(INOUT) :: charges
    TYPE( environ_electrons ), OPTIONAL, TARGET, INTENT(IN) :: electrons
    TYPE( environ_ions ),      OPTIONAL, TARGET, INTENT(IN) :: ions
    TYPE( environ_externals ), OPTIONAL, TARGET, INTENT(IN) :: externals
    TYPE( environ_auxiliary ), OPTIONAL, TARGET, INTENT(IN) :: auxiliary

    IF ( PRESENT(ions) ) THEN
       charges%include_ions = .TRUE.
       charges%ions => ions
    END IF

    IF ( PRESENT(electrons) ) THEN
       charges%include_electrons = .TRUE.
       charges%electrons => electrons
    ENDIF

    IF ( PRESENT(externals) ) THEN
       charges%include_externals = .TRUE.
       charges%externals => externals
    ENDIF

    IF ( PRESENT(auxiliary) ) THEN
       charges%include_auxiliary = .TRUE.
       charges%auxiliary => auxiliary
    ENDIF

  END SUBROUTINE init_environ_charges_first

  SUBROUTINE init_environ_charges_second( cell, charges )

    IMPLICIT NONE

    TYPE( environ_cell ), INTENT( IN ) :: cell
    TYPE( environ_charges ), INTENT( INOUT ) :: charges


    CALL init_environ_density( cell, charges%density )

    RETURN

  END SUBROUTINE init_environ_charges_second

  SUBROUTINE update_environ_charges( charges )

    IMPLICIT NONE

    TYPE( environ_charges ), INTENT( INOUT ) :: charges

    REAL( DP ) :: local_charge
    CHARACTER( LEN = 80 ) :: sub_name = 'update_environ_charges'

    charges % number = 0
    charges % charge = 0.D0
    charges % density % of_r = 0.D0

    IF ( charges % include_electrons ) THEN
       IF ( .NOT. ASSOCIATED( charges % electrons ) ) &
            & CALL errore(sub_name,'Missing expected charge component',1)
       charges % number = charges % number + charges % electrons % number
       charges % charge = charges % charge + charges % electrons % charge
       charges % density % of_r = charges % density % of_r + charges % electrons % density % of_r
    ENDIF

    IF ( charges % include_ions ) THEN
       IF ( .NOT. ASSOCIATED( charges % ions ) ) &
            & CALL errore(sub_name,'Missing expected charge component',1)
       charges % number = charges % number + charges % ions % number
       charges % charge = charges % charge + charges % ions % charge
       charges % density % of_r = charges % density % of_r + charges % ions % density % of_r
    ENDIF

    IF ( charges % include_externals ) THEN
       IF ( .NOT. ASSOCIATED( charges % externals ) ) &
            & CALL errore(sub_name,'Missing expected charge component',1)
       charges % number = charges % number + charges % externals % number
       charges % charge = charges % charge + charges % externals  % charge
       charges % density % of_r = charges % density % of_r + charges % externals % density % of_r
    ENDIF

    IF ( charges % include_auxiliary ) THEN
       IF ( .NOT. ASSOCIATED( charges % auxiliary ) ) &
            & CALL errore(sub_name,'Missing expected charge component',1)
       charges % number = charges % number + charges % auxiliary % number
       charges % charge = charges % charge + charges % auxiliary % charge
       charges % density % of_r = charges % density % of_r + charges % auxiliary % density % of_r
    ENDIF

    local_charge = integrate_environ_density(charges%density)
    IF ( ABS(local_charge-charges%charge) .GT. 1.D-8 ) CALL errore(sub_name,'Inconsistent integral of total charge',1)

    RETURN

  END SUBROUTINE update_environ_charges

  SUBROUTINE destroy_environ_charges( lflag, charges )

    IMPLICIT NONE

    LOGICAL, INTENT(IN) :: lflag
    TYPE( environ_charges ) :: charges
    CHARACTER( LEN=80 ) :: sub_name = 'destroy_environ_charges'

    IF ( lflag ) THEN

       ! These components were allocated first, destroy only if lflag = .TRUE.

       IF (ASSOCIATED(charges%ions)) NULLIFY(charges%ions)

       IF (ASSOCIATED(charges%electrons)) NULLIFY( charges%electrons )

       IF (ASSOCIATED(charges%externals)) NULLIFY( charges%externals )

       IF (ASSOCIATED(charges%auxiliary)) NULLIFY( charges%auxiliary )

    END IF

    CALL destroy_environ_density( charges%density )

    RETURN

  END SUBROUTINE destroy_environ_charges

END MODULE environ_types
