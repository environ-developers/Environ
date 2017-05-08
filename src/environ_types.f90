MODULE environ_types

  USE kinds,       ONLY : DP
  USE constants,   ONLY : rydberg_si, bohr_radius_si, amu_si, fpi

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

     CHARACTER( LEN=80 ) :: label = " "

     ! Each quantity in real-space is associated with its definition domain

     TYPE( environ_cell ), POINTER :: cell => NULL()

     ! The quantity in real-space, local to each processor

     REAL( DP ), DIMENSION(:), ALLOCATABLE :: of_r

  END TYPE environ_density

  TYPE environ_gradient

     ! Optionally have an associated logical status

     LOGICAL :: update = .FALSE.

     ! Optionally have an associated label, used for printout and debugs

     CHARACTER( LEN=80 ) :: label = " "

     ! Each quantity in real-space is associated with its definition domain

     TYPE( environ_cell ), POINTER :: cell => NULL()

     ! The quantity in real-space, local to each processor

     REAL( DP ), DIMENSION(:,:), ALLOCATABLE :: of_r

     REAL( DP ), DIMENSION(:), ALLOCATABLE :: modulus

  END TYPE environ_gradient

  TYPE environ_hessian

     ! Optionally have an associated logical status

     LOGICAL :: update = .FALSE.

     ! Optionally have an associated label, used for printout and debugs

     CHARACTER( LEN=80 ) :: label = " "

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

  TYPE environ_charges

     ! Ionic charges

     LOGICAL :: include_ions = .FALSE.
     INTEGER :: nions = 0
     REAL( DP ) :: zions = 0.0_DP
     TYPE( environ_ions ), POINTER :: ions => NULL()

     ! Electrons

     LOGICAL :: include_electrons = .FALSE.
     INTEGER :: nelectrons = 0
     REAL( DP ) :: zelectrons = 0.0_DP
     TYPE( environ_electrons ), POINTER :: electrons => NULL()

     ! External charges

     LOGICAL :: include_externals = .FALSE.
     INTEGER :: nexternals = 0
     REAL( DP ) :: zexternals = 0.0_DP
     TYPE( environ_externals ), POINTER :: externals => NULL()

     ! Auxiliary charges

     LOGICAL :: include_auxiliary = .FALSE.
     INTEGER :: nauxiliary = 0
     REAL( DP ) :: zauxiliary = 0.0_DP
     TYPE( environ_density ), POINTER :: auxiliary => NULL()

     ! Total smooth free charge

     INTEGER :: ntot = 0
     REAL( DP ) :: ztot = 0.0_DP
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

     ! Parameters for the density-dependent interface

     INTEGER :: type
     REAL( DP ) :: rhomax, rhomin, fact
     REAL( DP ) :: rhozero, deltarho, tbeta
     TYPE( environ_density ) :: density

     ! Parameters for the electrons-dependent interface

     LOGICAL :: need_electrons
     TYPE( environ_electrons ), POINTER :: electrons

     ! Parameters for the ions-dependent interface

     LOGICAL :: need_ions
     REAL( DP ) :: alpha ! solvent-dependent scaling factor
     REAL( DP ) :: softness ! sharpness of the interface
     TYPE( environ_functions ), DIMENSION(:), ALLOCATABLE :: soft_spheres
     TYPE( environ_ions ), POINTER :: ions

     ! scaled switching function of interface
     ! varying from 1 (QM region) to 0 (environment region)
     ! WARNING::: for consistency with previous version
     ! for the time being this is instead the dielectric function
     ! and we store here the bulk dielectric constant and the
     ! scaling factor

     REAL( DP ) :: constant
     REAL( DP ) :: scaling_factor

     TYPE( environ_density ) :: scaled
     TYPE( environ_density ) :: dscaled
     TYPE( environ_density ) :: d2scaled

     ! surface function of interface, computed from the finite
     ! difference of two scaled swtiching functions

     LOGICAL :: need_theta
     REAL( DP ) :: deltatheta
     TYPE( environ_density ) :: theta

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
  SUBROUTINE init_environ_cell( n1, n2, n3, ibrav, omega, at, nnr, ir_end, comm, me, root, cell )

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: n1, n2, n3, ibrav
    INTEGER, INTENT(IN) :: nnr, ir_end, comm, me, root
    REAL( DP ), INTENT(IN) :: omega, at(3,3)
    TYPE( environ_cell ), INTENT(INOUT) :: cell
    CHARACTER( LEN=80 ) :: sub_name = 'init_environ_cell'

    cell % n1 = n1
    cell % n2 = n2
    cell % n3 = n3
    cell % ibrav = ibrav
    cell % omega = omega
    cell % at = at

    cell % nnr = nnr
    cell % ir_end = ir_end
    cell % comm = comm
    cell % me   = me
    cell % root = root

    RETURN

  END SUBROUTINE init_environ_cell

  SUBROUTINE update_environ_cell( omega, at, cell )

    IMPLICIT NONE

    REAL( DP ), INTENT(IN) :: omega, at(3,3)
    TYPE( environ_cell ), INTENT(INOUT) :: cell
    CHARACTER( LEN=80 ) :: sub_name = 'update_environ_cell'

    cell % omega = omega
    cell % at = at

    RETURN

  END SUBROUTINE update_environ_cell
!----------------------------------------------------------------------------------------------------------------------------------------
!- DENSITY ------------------------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------------------------------
  SUBROUTINE create_environ_density(density)

    IMPLICIT NONE

    TYPE( environ_density ), INTENT(INOUT) :: density
    CHARACTER( LEN=80 ) :: sub_name = 'create_environ_density'

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

    REAL( DP ) :: scalar_product
    CHARACTER( LEN=80 ) :: fun_name = 'scalar_product_environ_density'

    IF ( .NOT.ASSOCIATED(density1%cell,density2%cell) ) &
        & CALL errore(fun_name,'operation on fields with inconsistent domains',1)
    scalar_product = DOT_PRODUCT(density1%of_r,density2%of_r)
    CALL mp_sum( scalar_product, density1%cell%comm )

    RETURN

  END FUNCTION scalar_product_environ_density

  FUNCTION quadratic_mean_environ_density(density) RESULT(quadratic_mean)

    IMPLICIT NONE

    TYPE( environ_density ), INTENT(IN) :: density

    REAL( DP ) :: quadratic_mean

    quadratic_mean = DOT_PRODUCT(density%of_r,density%of_r)
    CALL mp_sum( quadratic_mean, density%cell%comm )
    quadratic_mean = SQRT( quadratic_mean / density % cell % ntot )

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
  SUBROUTINE create_environ_gradient(gradient)

    IMPLICIT NONE

    TYPE( environ_gradient ), INTENT(INOUT) :: gradient
    CHARACTER (LEN=80) :: sub_name = 'destroy_environ_density'

    NULLIFY(gradient%cell)
    IF ( ALLOCATED( gradient%of_r ) ) CALL errore(sub_name,'Trying to create an already allocated object',1)
    IF ( ALLOCATED( gradient%modulus ) ) CALL errore(sub_name,'Trying to create an already allocated object',1)

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

    IF ( ALLOCATED( gradient%modulus ) ) CALL errore(sub_name,'Trying to allocate an allocated object',1)
    ALLOCATE(gradient%modulus(gradient%cell%nnr))
    gradient%modulus = 0.D0

    RETURN

  END SUBROUTINE init_environ_gradient

  SUBROUTINE update_gradient_modulus(gradient)

    IMPLICIT NONE

    TYPE( environ_gradient ), INTENT(INOUT) :: gradient
    INTEGER, POINTER :: ir_end

    ir_end => gradient % cell % ir_end
    gradient%modulus(1:ir_end) = gradient%of_r(1,1:ir_end)**2 + &
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

    IF (.NOT.ALLOCATED(gradient%modulus)) &
         & CALL errore(sub_name,'Trying to destroy a non allocated object',1)
    DEALLOCATE( gradient%modulus )

    RETURN

  END SUBROUTINE destroy_environ_gradient
!----------------------------------------------------------------------------------------------------------------------------------------
!- HESSIAN ------------------------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------------------------------
  SUBROUTINE create_environ_hessian(hessian)

    IMPLICIT NONE

    TYPE( environ_hessian ), INTENT(INOUT) :: hessian
    CHARACTER (LEN=80) :: sub_name = 'destroy_environ_hessian'

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
!- IONS ---------------------------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------------------------------
  SUBROUTINE create_environ_ions(ions)

    IMPLICIT NONE

    TYPE( environ_ions ), INTENT(INOUT) :: ions
    CHARACTER (LEN=80) :: sub_name = 'create_environ_ions'

    ions%update = .FALSE.

    IF ( ALLOCATED( ions%ityp ) ) CALL errore(sub_name,'Trying to create an already allocated object',1)
    IF ( ALLOCATED( ions%iontype ) ) CALL errore(sub_name,'Trying to create an already allocated object',1)

    NULLIFY( ions%tau )

    ions%use_smeared_ions = .FALSE.
    IF ( ALLOCATED( ions%smeared_ions ) ) CALL errore(sub_name,'Trying to create an already allocated object',1)
    CALL create_environ_density( ions%density )

    ions%use_core_electrons = .FALSE.
    IF ( ALLOCATED( ions%core_electrons ) ) CALL errore(sub_name,'Trying to create an already allocated object',1)
    CALL create_environ_density( ions%core )

    RETURN

  END SUBROUTINE create_environ_ions

  SUBROUTINE init_environ_ions_first( nat, ntyp, lsoftcavity, lcoredensity, &
       &  lsmearedions, radius_mode, atom_label, atomicspread, corespread, solvationrad, ions )

    ! First step of ions initialization, cannot initialize everything
    ! because some infos are missing when this routine is called
    ! NEED TO REVISE THE POSITION OF THE CALL INSIDE QE TO MERGE FIRST AND SECOND

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nat, ntyp
    LOGICAL, INTENT(IN) :: lsoftcavity
    LOGICAL, INTENT(IN) :: lcoredensity
    LOGICAL, INTENT(IN) :: lsmearedions
    CHARACTER( LEN=80 ), INTENT(IN) :: radius_mode
    CHARACTER( LEN=3 ), DIMENSION(ntyp), INTENT(IN) :: atom_label
    REAL( DP ), DIMENSION(ntyp), INTENT(IN) :: atomicspread
    REAL( DP ), DIMENSION(ntyp), INTENT(IN) :: corespread
    REAL( DP ), DIMENSION(ntyp), INTENT(IN) :: solvationrad
    TYPE( environ_ions ), INTENT(INOUT) :: ions

    CHARACTER( LEN = 80 ) :: sub_name = 'init_environ_ions_first'

    INTEGER :: i

    ions%center = 0.D0

    ions%number = nat
    ions%ntyp = ntyp

    ! Allocate the basic vectors, cannot initialize them here

    ALLOCATE( ions%tau( 3, nat ) )
    ALLOCATE( ions%ityp( nat ) )

    ! Set ions types, note that also valence charges cannot be initialized here

    ALLOCATE( ions%iontype( ntyp ) )

    DO i = 1, ntyp

       ! Given the label we could set some of the properties with defaults

       CALL set_iontype_defaults( i, atom_label(i), radius_mode, ions%iontype(i) )

       ! Check if values were provided in input and overwrite them

       IF ( atomicspread(i) .GT. 0 ) ions%iontype(i)%atomicspread = atomicspread(i)
       IF ( corespread(i) .GT. 0 ) ions%iontype(i)%corespread = corespread(i)
       IF ( solvationrad(i) .GT. 0 ) ions%iontype(i)%solvationrad = solvationrad(i)

       ! If need cavity defined exclusively on ions, check radius is not zero

       IF ( .NOT. lsoftcavity .AND. ( ions%iontype(i)%solvationrad .EQ. 0.D0 ) ) &
            & CALL errore(sub_name,'Missing solvation radius for one of the atom types',1)

       ! If need smeared ions, check spread is not zero

       IF ( lsmearedions .AND. ( ions%iontype(i)%atomicspread .EQ. 0.D0 ) ) &
            & CALL errore(sub_name,'Missing atomic spread for one of the atom types',1)

    END DO

    ions%use_smeared_ions = lsmearedions

    ions%use_core_electrons = lcoredensity

    RETURN

  END SUBROUTINE init_environ_ions_first

  SUBROUTINE init_environ_ions_second( nat, ntyp, ityp, zv, cell, ions )

    ! Second step of initialization, passing the information on types,
    ! atomic valence charges and whether we need to compute the smeared density or not

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nat, ntyp
    INTEGER, DIMENSION(nat), INTENT(IN) :: ityp
    REAL(DP), DIMENSION(ntyp), INTENT(IN) :: zv
    TYPE( environ_cell), INTENT(IN) :: cell
    TYPE( environ_ions ), INTENT(INOUT) :: ions

    CHARACTER( LEN = 80 ) :: sub_name = 'init_environ_ions_second'

    INTEGER :: i

    ! Check on dimensions, can skip if merged with first step

    IF ( ions%number .NE. nat ) CALL errore(sub_name,'Mismatch in number of atoms',1)
    IF ( ions%ntyp .NE. ntyp ) CALL errore(sub_name,'Mismatch in number of atom types',1)

    ions%ityp = ityp

    DO i = 1, ions%ntyp

       ions%iontype(i)%zv = zv(i)

    ENDDO

    ions%charge = 0.D0

    DO i = 1, ions%number

       ions%charge = ions%charge + ions%iontype(ions%ityp(i))%zv

    ENDDO

    IF ( ions%use_smeared_ions ) THEN
       ! THE FOLLOWING TEST ON ALLOCATION IS ONLY THERE BECAUSE OF WHEN THIS INITIALIZATION
       ! IS CALLED, IF MERGED WITH THE FIRST STEP REMOVE THE TEST
       IF ( .NOT. ALLOCATED( ions%density%of_r ) ) THEN
          CALL init_environ_density( cell, ions%density )
       ELSE
          ions%density%of_r = 0.D0
       ENDIF

       ! Build smeared ions from iontype data

       IF ( .NOT. ALLOCATED( ions%smeared_ions ) ) THEN
          ALLOCATE( ions%smeared_ions( ions%number ) )
          DO i = 1, ions%number
             ions%smeared_ions(i) = environ_functions(1,1,0,0.0_DP,&
                  & ions%iontype(ions%ityp(i))%atomicspread,&
                  & ions%iontype(ions%ityp(i))%zv,ions%tau(:,i))
          ENDDO
       ENDIF

    ENDIF

    IF ( ions%use_core_electrons ) THEN
       ! THE FOLLOWING TEST ON ALLOCATION IS ONLY THERE BECAUSE OF WHEN THIS INITIALIZATION
       ! IS CALLED, IF MERGED WITH THE FIRST STEP REMOVE THE TEST
       IF ( .NOT. ALLOCATED( ions%core%of_r ) ) THEN
          CALL init_environ_density( cell, ions%core )
       ELSE
          ions%core%of_r = 0.D0
       ENDIF

       ! Build core electrons from iontype data

       IF ( .NOT. ALLOCATED( ions%core_electrons ) ) THEN
          ALLOCATE( ions%core_electrons( ions%number ) )
          DO i = 1, ions%number
             ions%core_electrons(i) = environ_functions(1,1,0,0.0_DP,&
                  & ions%iontype(ions%ityp(i))%corespread,&
                  & ions%iontype(ions%ityp(i))%zv,ions%tau(:,i))
          ENDDO
       ENDIF

    ENDIF

    RETURN

  END SUBROUTINE init_environ_ions_second

  SUBROUTINE update_environ_ions( nat, tau, ions )

    ! Update ionic positions and compute derived quantities

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nat
    REAL(DP), DIMENSION(3,nat), INTENT(IN) :: tau
    TYPE( environ_ions ), INTENT(INOUT) :: ions

    INTEGER :: i
    INTEGER :: dim, axis
    REAL(DP) :: charge, spread
    REAL(DP), DIMENSION(3) :: pos
    CHARACTER(LEN=80) :: sub_name = 'update_environ_ions'

    ! Check on dimensions

    IF ( ions%number .NE. nat ) CALL errore(sub_name,'Mismatch in number of atoms',1)

    ! Update positions

    ions%tau = tau

    ! Center of ionic charge used by three sub-modules

    ions%center = 0.D0
    DO i = 1, ions%number

       ions%center(:) = ions%center(:) + ions%tau(:,i)*ions%iontype(ions%ityp(i))%zv

    ENDDO
    ions%center = ions%center / ions%charge

    ! If needed, generate a fictitious ion density using gaussians

    IF ( ions%use_smeared_ions ) CALL density_of_functions(ions%number,ions%smeared_ions,ions%density)

    RETURN

  END SUBROUTINE update_environ_ions

  SUBROUTINE destroy_environ_ions(lflag, ions)

    IMPLICIT NONE

    LOGICAL, INTENT(IN) :: lflag
    TYPE( environ_ions ), INTENT(INOUT) :: ions
    CHARACTER (LEN=80) :: sub_name = 'destroy_environ_ions'

    ! ityp, tau and iontype should have been allocated
    ! raise an error if they are not

    IF ( lflag ) THEN

       ! These components were allocated first, only destroy if lflag = .TRUE.

       IF (.NOT.ALLOCATED(ions%ityp)) &
            & CALL errore(sub_name,'Trying to destroy a non allocated object',1)
       DEALLOCATE( ions%ityp )
       IF (.NOT.ALLOCATED(ions%iontype)) &
            & CALL errore(sub_name,'Trying to destroy a non allocated object',1)
       DEALLOCATE( ions%iontype )

       IF (.NOT.ASSOCIATED(ions%tau)) &
            & CALL errore(sub_name,'Trying to destroy a non associated object',1)
       DEALLOCATE( ions%tau )

    ENDIF

    IF ( ions%use_smeared_ions ) THEN
       CALL destroy_environ_density( ions%density )
       CALL destroy_environ_functions( ions%number, ions%smeared_ions )
    ENDIF

    RETURN

  END SUBROUTINE destroy_environ_ions
!----------------------------------------------------------------------------------------------------------------------------------------
!- ELECTRONS ----------------------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------------------------------
  SUBROUTINE create_environ_electrons(electrons)

    IMPLICIT NONE

    TYPE( environ_electrons ), INTENT(INOUT) :: electrons

    electrons%update = .FALSE.
    electrons%number = 0
    electrons%nspin  = 1
    electrons%charge = 0.D0
    CALL create_environ_density( electrons%density )

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

    INTEGER, INTENT(IN) :: nelec, nspin, nnr
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

    electrons%number = nelec
    electrons%charge = integrate_environ_density( electrons%density )
    IF ( ABS(electrons%charge-DBLE(electrons%number)) .GT. tol ) &
         & CALL errore(sub_name,'Missmatch in integrated electronic charge',1)

    RETURN

  END SUBROUTINE update_environ_electrons

  SUBROUTINE destroy_environ_electrons( electrons )

    IMPLICIT NONE

    TYPE( environ_electrons ), INTENT(INOUT) :: electrons

    CALL destroy_environ_density( electrons%density )

    RETURN

  END SUBROUTINE destroy_environ_electrons
!----------------------------------------------------------------------------------------------------------------------------------------
!- EXTERNALS ----------------------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------------------------------
  SUBROUTINE create_environ_externals(externals)

    IMPLICIT NONE

    TYPE( environ_externals ), INTENT(INOUT) :: externals
    CHARACTER (LEN=80) :: sub_name = 'create_environ_externals'

    externals%update = .FALSE.
    externals%number = 0
    IF ( ALLOCATED( externals%functions ) ) CALL errore(sub_name,'Trying to create an already allocated object',1)
    CALL create_environ_density( externals%density )
    externals%charge = 0.D0

    RETURN

  END SUBROUTINE create_environ_externals

  SUBROUTINE init_environ_externals_first( nexternals, dims, axis, pos, &
       & spreads, charge, externals )

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nexternals
    INTEGER, DIMENSION(nexternals), INTENT(IN) :: dims, axis
    REAL( DP ), DIMENSION(3,nexternals), INTENT(IN) :: pos
    REAL( DP ), DIMENSION(nexternals), INTENT(IN) :: spreads, charge
    TYPE( environ_externals ), INTENT(INOUT) :: externals

    INTEGER :: i

    externals%number = nexternals
    ALLOCATE(externals%functions(externals%number))
    DO i = 1, externals%number
       ALLOCATE(externals%functions(i)%pos(3))
       externals%functions(i)%dim    = dims(i)
       externals%functions(i)%axis   = axis(i)
       externals%functions(i)%pos(:) = pos(:,i)
       externals%functions(i)%spread = spreads(i)
       externals%functions(i)%width  = spreads(i)
       externals%functions(i)%volume = charge(i)
    ENDDO

    RETURN

  END SUBROUTINE init_environ_externals_first

  SUBROUTINE init_environ_externals_second( cell, externals )

    IMPLICIT NONE

    TYPE( environ_cell ), INTENT(IN) :: cell
    TYPE( environ_externals ), INTENT(INOUT) :: externals

    CALL init_environ_density( cell, externals%density )

    RETURN

  END SUBROUTINE init_environ_externals_second

  SUBROUTINE update_environ_externals( externals )

    IMPLICIT NONE

    TYPE( environ_externals ), INTENT(INOUT) :: externals

    CALL density_of_functions( externals%number, externals%functions, externals%density )

    RETURN

  END SUBROUTINE update_environ_externals

  SUBROUTINE destroy_environ_externals( lflag, externals )

    IMPLICIT NONE

    LOGICAL, INTENT(IN) :: lflag
    TYPE( environ_externals ), INTENT(INOUT) :: externals

    IF ( lflag ) CALL destroy_environ_functions( externals%number, externals%functions )

    CALL destroy_environ_density( externals%density )

    RETURN

  END SUBROUTINE destroy_environ_externals
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
!- CHARGES ------------------------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------------------------------
  SUBROUTINE create_environ_charges(charges)

    IMPLICIT NONE

    TYPE( environ_charges ) :: charges

    charges%include_ions = .FALSE.
    charges%nions = 0
    charges%zions = 0.D0
    NULLIFY( charges%ions )

    charges%include_electrons = .FALSE.
    charges%nelectrons = 0
    charges%zelectrons = 0.D0
    NULLIFY( charges%electrons )

    charges%include_externals = .FALSE.
    charges%nexternals = 0
    charges%zexternals = 0.D0
    NULLIFY( charges%externals )

    charges%include_auxiliary = .FALSE.
    charges%nauxiliary = 0
    charges%zauxiliary = 0.D0
    NULLIFY( charges%auxiliary )

    charges%ntot = 0
    charges%ztot = 0.D0
    CALL create_environ_density( charges%density )

    RETURN

  END SUBROUTINE create_environ_charges

  SUBROUTINE init_environ_charges_first( charges, electrons, ions, externals, auxiliary )

    IMPLICIT NONE

    TYPE( environ_charges ), INTENT(INOUT) :: charges
    TYPE( environ_electrons ), OPTIONAL, TARGET, INTENT(IN) :: electrons
    TYPE( environ_ions ),      OPTIONAL, TARGET, INTENT(IN) :: ions
    TYPE( environ_externals ), OPTIONAL, TARGET, INTENT(IN) :: externals
    TYPE( environ_density ),   OPTIONAL, TARGET, INTENT(IN) :: auxiliary

    IF ( PRESENT(ions) ) THEN
       charges%include_ions = .TRUE.
       charges%nions = ions%number
       charges%ions => ions
    END IF

    IF ( PRESENT(electrons) ) THEN
       charges%include_electrons = .TRUE.
       charges%nelectrons = electrons%number
       charges%electrons => electrons
    ENDIF

    IF ( PRESENT(externals) ) THEN
       charges%include_externals = .TRUE.
       charges%nexternals = externals%number
       charges%externals => externals
    ENDIF

    IF ( PRESENT(auxiliary) ) THEN
       charges%include_auxiliary = .TRUE.
       charges%nauxiliary = 1
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

    charges % ntot = 0.D0
    charges % density % of_r = 0.D0

    IF ( charges % include_electrons ) THEN
       charges % ntot = charges % ntot + charges % nelectrons
       IF ( charges % nelectrons .GT. 0 ) charges % density % of_r = &
            & charges % density % of_r + charges % electrons % density % of_r
    ENDIF

    IF ( charges % include_ions ) THEN
       charges % ntot = charges % ntot + charges % nions
       IF ( charges % nions .GT. 0 ) charges % density % of_r = &
            & charges % density % of_r + charges % ions % density % of_r
    ENDIF

    IF ( charges % include_externals ) THEN
       charges % ntot = charges % ntot + charges % nexternals
       IF ( charges % nexternals .GT. 0 ) charges % density % of_r = &
            & charges % density % of_r + charges % externals % density % of_r
    ENDIF

    IF ( charges % include_auxiliary ) THEN
       charges % ntot = charges % ntot + charges % nauxiliary
       IF ( charges % nauxiliary .GT. 0 ) charges % density % of_r = &
            & charges % density % of_r + charges % auxiliary % of_r
    ENDIF

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
!----------------------------------------------------------------------------------------------------------------------------------------
!- BOUNDARY -----------------------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------------------------------
  SUBROUTINE create_environ_boundary(boundary)

    IMPLICIT NONE

    TYPE( environ_boundary ), INTENT(INOUT) :: boundary

    CHARACTER( LEN=80 ) :: sub_name = 'create_environ_boundary'

    boundary%update_status = 0

    CALL create_environ_density( boundary%scaled )
    CALL create_environ_density( boundary%dscaled )
    CALL create_environ_density( boundary%d2scaled )

    boundary%need_electrons = .FALSE.
    NULLIFY( boundary%electrons )
    boundary%need_ions = .FALSE.
    NULLIFY( boundary%ions   )
    boundary%need_theta = .FALSE.
    CALL create_environ_density( boundary%theta )

    IF ( ALLOCATED( boundary%soft_spheres ) ) &
         & CALL errore(sub_name,'Trying to create an already allocated object',1)

    CALL create_environ_density( boundary%density )

    RETURN

  END SUBROUTINE create_environ_boundary

  SUBROUTINE init_environ_boundary_first( mode, constant, type, &
       & rhomax, rhomin, tbeta, need_theta, delta, alpha, &
       & softness, ions, boundary )

    IMPLICIT NONE

    CHARACTER( LEN=80 ), INTENT(IN) :: mode
    REAL( DP ), INTENT(IN) :: constant
    INTEGER, INTENT(IN) :: type
    REAL( DP ), INTENT(IN) :: rhomax, rhomin, tbeta
    LOGICAL, INTENT(IN) :: need_theta
    REAL( DP ), INTENT(IN) :: delta
    REAL( DP ), INTENT(IN) :: alpha
    REAL( DP ), INTENT(IN) :: softness
    TYPE( environ_ions ), TARGET, INTENT(IN) :: ions
    TYPE( environ_boundary ), INTENT(INOUT) :: boundary

    INTEGER :: i
    REAL( DP ) :: radius

    boundary%mode = mode
    boundary%type = type
    boundary%rhomax = rhomax
    boundary%rhomin = rhomin
    boundary%fact = LOG( rhomax / rhomin )
    boundary%rhozero = ( rhomax + rhomin ) * 0.5_DP
    boundary%tbeta = tbeta
    boundary%deltarho = rhomax - rhomin

    IF ( constant .GT. 1.D0 ) THEN
       boundary%constant = constant
       boundary%scaling_factor = -1.D0/(constant-1.D0)
    ELSE
       boundary%constant = 2.D0
       boundary%scaling_factor = -1.D0
    ENDIF

    boundary%need_electrons = ( mode .EQ. 'electronic' ) .OR. ( mode .EQ. 'full' )

    boundary%need_ions = ( mode .EQ. 'ionic' ) .OR. ( mode .EQ. 'full' )
    IF ( boundary%need_ions ) boundary%ions => ions
    boundary%alpha = alpha
    boundary%softness = softness
    IF ( boundary%need_ions .AND. .NOT. boundary%need_electrons ) THEN
       ALLOCATE( boundary%soft_spheres( boundary%ions%number ) )
       DO i = 1, boundary%ions%number
          radius = boundary%ions%iontype(boundary%ions%ityp(i))%solvationrad * boundary%alpha
          boundary%soft_spheres(i) = environ_functions(2,0,1,radius,1.0_DP,boundary%softness,&
                  & boundary%ions%tau(:,i))
       ENDDO
    ENDIF

    boundary%need_theta = need_theta
    boundary%deltatheta = delta

    RETURN

  END SUBROUTINE init_environ_boundary_first

  SUBROUTINE init_environ_boundary_second( cell, boundary )

    IMPLICIT NONE

    TYPE( environ_cell ), INTENT(IN) :: cell
    TYPE( environ_boundary ), INTENT(INOUT) :: boundary

    CALL init_environ_density( cell, boundary%scaled )
    CALL init_environ_density( cell, boundary%dscaled )
    CALL init_environ_density( cell, boundary%d2scaled )

    IF ( boundary%need_electrons ) CALL init_environ_density( cell, boundary%density )
    IF ( boundary%need_theta) CALL init_environ_density( cell, boundary%theta )

    RETURN

  END SUBROUTINE init_environ_boundary_second

  SUBROUTINE destroy_environ_boundary(lflag, boundary)

    IMPLICIT NONE

    LOGICAL, INTENT(IN) :: lflag
    TYPE( environ_boundary ), INTENT(INOUT) :: boundary
    CHARACTER (LEN=80) :: sub_name = 'destroy_environ_boundary'

    IF ( lflag ) THEN

       ! These components were allocated first, destroy only if lflag = .TRUE.

       IF ( boundary%need_ions ) THEN
          IF (.NOT.ASSOCIATED(boundary%ions)) &
               & CALL errore(sub_name,'Trying to destroy a non associated object',1)
          NULLIFY(boundary%ions)
          IF ( .NOT. boundary%need_electrons ) &
               & CALL destroy_environ_functions( boundary%ions%number, boundary%soft_spheres )
       ELSE
          IF (ASSOCIATED(boundary%ions))&
               & CALL errore(sub_name,'Found an unexpected associated object',1)
       ENDIF

       IF ( boundary%need_electrons ) THEN
          IF (ASSOCIATED(boundary%electrons)) NULLIFY(boundary%electrons)
       ENDIF

    ENDIF

    CALL destroy_environ_density( boundary%scaled )
    CALL destroy_environ_density( boundary%dscaled )
    CALL destroy_environ_density( boundary%d2scaled )

    IF ( boundary%need_electrons ) CALL destroy_environ_density( boundary%density )
    IF ( boundary%need_theta ) CALL destroy_environ_density( boundary%theta )

    RETURN

  END SUBROUTINE destroy_environ_boundary
!----------------------------------------------------------------------------------------------------------------------------------------
!- DIELECTRIC ---------------------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------------------------------
  SUBROUTINE create_environ_dielectric(dielectric)

    IMPLICIT NONE

    TYPE( environ_dielectric ), INTENT(INOUT) :: dielectric

    CHARACTER( LEN=80 ) :: sub_name = 'create_environ_dielectric'

    dielectric%constant = 1.0_DP

    IF ( ALLOCATED( dielectric%regions ) ) CALL errore(sub_name,'Trying to create an already allocated object',1)

    CALL create_environ_density( dielectric%background  )
    CALL create_environ_density( dielectric%epsilon     )

    NULLIFY( dielectric%boundary )

    dielectric%need_gradient = .FALSE.
    CALL create_environ_gradient( dielectric%gradient )
    dielectric%need_factsqrt = .FALSE.
    CALL create_environ_density( dielectric%factsqrt )
    dielectric%need_gradlog = .FALSE.
    CALL create_environ_gradient( dielectric%gradlog )
    RETURN

  END SUBROUTINE create_environ_dielectric

  SUBROUTINE create_environ_functions( n, type, dimm, axis, pos, width, spreadd, volume, f )

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: n
    INTEGER, INTENT(IN) :: type
    INTEGER, DIMENSION(:), INTENT(IN) :: dimm, axis
    REAL( DP ), DIMENSION(n), INTENT(IN) :: width, spreadd, volume
    REAL( DP ), DIMENSION(3,n), TARGET, INTENT(IN) :: pos
    TYPE( environ_functions ), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: f

    INTEGER :: i

    ALLOCATE( f(n) )
    DO i = 1, n
       f(i)%type = type
       f(i)%dim  = dimm(i)
       f(i)%axis = axis(i)
       f(i)%spread = spreadd(i)
       f(i)%width = width(i)
       f(i)%volume = volume(i)
       f(i)%pos => pos(:,i)
    ENDDO

    RETURN

  END SUBROUTINE create_environ_functions

  SUBROUTINE destroy_environ_functions( n, f )

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: n
    TYPE( environ_functions ), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: f

    CHARACTER( LEN=80 ) :: sub_name = 'destroy_environ_functions'

    INTEGER :: i

    IF ( .NOT. ALLOCATED( f ) ) &
         & CALL errore(sub_name,'Trying to destroy a non allocated object',1)
    IF ( SIZE(f) .NE. n ) &
         & CALL errore(sub_name,'Inconsistent size of allocated object',1)

    DO i = 1, n
       IF ( .NOT. ASSOCIATED( f(i)%pos ) ) &
            & CALL errore(sub_name,'Trying to destroy a non allocated object',1)
       NULLIFY( f(i)%pos )
    ENDDO

    DEALLOCATE( f )

    RETURN

  END SUBROUTINE destroy_environ_functions

  SUBROUTINE init_environ_dielectric_first( constant, nregions, &
             & epsregion_dim, epsregion_axis, epsregion_pos, epsregion_width, &
             & epsregion_spread, epsregion_eps, boundary, dielectric )

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nregions
    INTEGER, DIMENSION(nregions), INTENT(IN) :: epsregion_dim, epsregion_axis
    REAL( DP ) :: constant
    REAL( DP ), DIMENSION(nregions), INTENT(IN) :: epsregion_width, epsregion_spread, epsregion_eps
    REAL( DP ), DIMENSION(3,nregions), INTENT(IN) :: epsregion_pos
    TYPE( environ_boundary ), TARGET, INTENT(IN) :: boundary
    TYPE( environ_dielectric ), INTENT(INOUT) :: dielectric

    dielectric%constant = constant

    dielectric%nregions = nregions
    IF ( dielectric%nregions .GT. 0 ) &
         & CALL create_environ_functions( dielectric%nregions, 2, epsregion_dim, &
         & epsregion_axis, epsregion_pos, epsregion_spread, epsregion_width, &
         & epsregion_eps, dielectric%regions )

    dielectric%boundary => boundary

    RETURN

  END SUBROUTINE init_environ_dielectric_first

  SUBROUTINE init_environ_dielectric_flags( need_gradient, need_factsqrt, need_gradlog, dielectric )

    IMPLICIT NONE

    LOGICAL, INTENT(IN) :: need_gradient, need_factsqrt, need_gradlog
    TYPE( environ_dielectric ), INTENT(INOUT) :: dielectric

    dielectric%need_gradient = need_gradient
    dielectric%need_factsqrt = need_factsqrt
    dielectric%need_gradlog  = need_gradlog

    RETURN

  END SUBROUTINE init_environ_dielectric_flags

  SUBROUTINE init_environ_dielectric_second( cell, dielectric )

    IMPLICIT NONE

    TYPE( environ_cell ), INTENT(IN) :: cell
    TYPE( environ_dielectric ), INTENT(INOUT) :: dielectric

    CHARACTER( LEN=80 ) :: sub_name = 'init_environ_dielectric_second'

    CALL init_environ_density( cell, dielectric%background )
    CALL init_environ_density( cell, dielectric%epsilon )

    IF ( dielectric%need_gradient ) CALL init_environ_gradient( cell, dielectric%gradient )
    IF ( dielectric%need_factsqrt ) CALL init_environ_density( cell, dielectric%factsqrt )
    IF ( dielectric%need_gradlog ) CALL init_environ_gradient( cell, dielectric%gradlog )

    RETURN

  END SUBROUTINE init_environ_dielectric_second

  SUBROUTINE destroy_environ_dielectric(lflag,dielectric)

    IMPLICIT NONE

    LOGICAL, INTENT(IN) :: lflag
    TYPE( environ_dielectric ), INTENT(INOUT) :: dielectric
    CHARACTER (LEN=80) :: sub_name = 'destroy_environ_dielectric'

    IF ( lflag ) THEN

       ! These components were allocated first, destroy only if lflag = .TRUE.

       IF ( dielectric%nregions .GT. 0 ) THEN
          CALL destroy_environ_functions( dielectric%nregions, dielectric%regions )
       ELSE
          IF ( ALLOCATED(dielectric%regions) ) &
               & CALL errore(sub_name,'Found unexpected allocated object',1)
       END IF

       IF (.NOT.ASSOCIATED(dielectric%boundary)) &
            & CALL errore(sub_name,'Trying to destroy a non associated object',1)
       NULLIFY( dielectric%boundary )

    END IF

    CALL destroy_environ_density( dielectric%background )
    CALL destroy_environ_density( dielectric%epsilon )

    IF ( dielectric%need_gradient ) CALL destroy_environ_gradient( dielectric%gradient )
    IF ( dielectric%need_factsqrt ) CALL destroy_environ_density( dielectric%factsqrt )
    IF ( dielectric%need_gradlog ) CALL destroy_environ_gradient( dielectric%gradlog )

    RETURN

  END SUBROUTINE destroy_environ_dielectric
!----------------------------------------------------------------------------------------------------------------------------------------
!- ELECTROLYTE --------------------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------------------------------
  SUBROUTINE create_environ_electrolyte( electrolyte )

    IMPLICIT NONE

    TYPE( environ_electrolyte ), INTENT(INOUT) :: electrolyte

    CHARACTER( LEN=80 ) :: sub_name = 'create_environ_electrolyte'

    CALL create_environ_boundary( electrolyte%boundary )
    CALL create_environ_density( electrolyte%density )

    IF ( ALLOCATED( electrolyte%ioncctype ) ) CALL errore(sub_name,'Trying to create an already allocated object',1)

    RETURN

  END SUBROUTINE create_environ_electrolyte

  SUBROUTINE init_environ_electrolyte_first( ntyp, mode, stype, rhomax, rhomin, &
             & tbeta, distance, spread, alpha, softness, ions, temperature, cbulk, &
             &  cmax, radius, z, electrolyte )

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: ntyp, stype
    CHARACTER( LEN=80 ), INTENT(IN) :: mode
    REAL( DP ), INTENT(IN) :: rhomax, rhomin, tbeta, distance, spread, alpha, softness, temperature
    REAL( DP ), DIMENSION(ntyp), INTENT(IN) :: cbulk, cmax, radius, z
    TYPE( environ_ions ), INTENT(IN) :: ions
    TYPE( environ_electrolyte ), INTENT(INOUT) :: electrolyte

    INTEGER :: ityp
    REAL( DP ) :: neutral
    CHARACTER( LEN=80 ) :: sub_name = 'init_environ_electrolyte_first'

    electrolyte%ntyp = ntyp

    electrolyte%temperature = temperature

    CALL init_environ_boundary_first( mode, 1.D0, stype, rhomax, rhomin, tbeta, .FALSE., &
         & 0.D0, alpha, softness, ions, electrolyte%boundary )

    ALLOCATE( electrolyte%ioncctype(ntyp) )

    neutral = 0.D0
    DO ityp = 1, ntyp
       ! If the radius is provided in input, compute cmax from it
       electrolyte%ioncctype(ityp)%cmax = cmax(ityp) * bohr_radius_si**3 / amu_si
       IF ( cmax(ityp) .EQ. 0.D0 .AND. radius(ityp) .GT. 0.D0 ) &
            & electrolyte%ioncctype(ityp)%cmax  = 0.64D0 * 3.D0 / fpi / radius(ityp)**3
       ! Double check that the bulk and max concentrations in input are compatible
       IF ( cbulk(ityp) .GT. 0.D0 .AND. cmax(ityp) .LT. cbulk(ityp) ) &
            & call errore (sub_name,'cmax should be at least greater than cbulk',1)
       electrolyte%ioncctype(ityp)%cbulk = cbulk(ityp) * bohr_radius_si**3 / amu_si
       electrolyte%ioncctype(ityp)%radius = radius(ityp)
       electrolyte%ioncctype(ityp)%z = z(ityp)
       neutral = neutral + cbulk(ityp)*z(ityp)
    END DO

    IF ( neutral .GT. 1.D-8 ) CALL errore(sub_name,'Bulk electrolyte is not neutral',1)

    RETURN

  END SUBROUTINE init_environ_electrolyte_first

  SUBROUTINE init_environ_electrolyte_second( cell, electrolyte )

    IMPLICIT NONE

    TYPE( environ_cell ), INTENT(IN) :: cell
    TYPE( environ_electrolyte ), INTENT(INOUT) :: electrolyte

    CALL init_environ_boundary_second( cell, electrolyte%boundary )

    CALL init_environ_density( cell, electrolyte%density )

    RETURN

  END SUBROUTINE init_environ_electrolyte_second

  SUBROUTINE destroy_environ_electrolyte( lflag, electrolyte )

    IMPLICIT NONE

    LOGICAL, INTENT(IN) :: lflag
    TYPE( environ_electrolyte ), INTENT(INOUT) :: electrolyte
    CHARACTER( LEN=80 ) :: sub_name = 'destroy_environ_electrolyte'

    IF ( lflag ) THEN

       ! These components were allocated first, destroy only if lflag = .TRUE.

       IF ( .NOT. ALLOCATED( electrolyte%ioncctype ) ) &
            & CALL errore(sub_name,'Trying to destroy a non allocated object',1)
       DEALLOCATE( electrolyte%ioncctype )

    ENDIF

    CALL destroy_environ_boundary( lflag, electrolyte%boundary )
    CALL destroy_environ_density( electrolyte%density )

    RETURN

  END SUBROUTINE destroy_environ_electrolyte
!----------------------------------------------------------------------------------------------------------------------------------------
!- UTILITIES ----------------------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------------------------------
  SUBROUTINE set_iontype_defaults(index,label,radius_mode,iontype)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: index
    CHARACTER( LEN=3 ), INTENT(IN) :: label
    CHARACTER( LEN=80 ), INTENT(IN) :: radius_mode
    TYPE(environ_iontype), INTENT(INOUT) :: iontype

    CHARACTER( LEN=80 ) :: sub_name = 'set_iontype_defaults'
    REAL( DP ), DIMENSION(92) :: pauling_radii
    DATA pauling_radii/ 1.20_DP, 0.00_DP, & ! H, -
         & 0.00_DP, 0.00_DP, 0.00_DP, 1.50_DP, 1.50_DP, 1.40_DP, 1.35_DP, 0.00_DP, & ! -, -, -, C, N, O, F, -
         & 0.00_DP, 0.00_DP, 0.00_DP, 0.00_DP, 1.90_DP, 1.85_DP, 1.80_DP, 0.00_DP, & ! -, -, -, -, P, S, Cl, -
         & 0.00_DP, 0.00_DP, 10*0.00_DP, 0.00_DP, 0.00_DP, 0.00_DP, 0.00_DP, 1.95_DP, 0.00_DP, & ! -,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,Br,-
         & 0.00_DP, 0.00_DP, 10*0.00_DP, 0.00_DP, 0.00_DP, 0.00_DP, 0.00_DP, 2.15_DP, 0.00_DP, & ! -,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,I,-
         & 38*0.00_DP / ! ...
    REAL( DP ), DIMENSION(92) :: bondi_radii
    DATA bondi_radii/ 1.20_DP, 1.40_DP, & ! H, He
         & 1.82_DP, 1.45_DP, 1.80_DP, 1.70_DP, 1.55_DP, 1.52_DP, 1.47_DP, 1.54_DP, & ! Li, Be, B, C, N, O, F, Ne
         & 2.27_DP, 1.73_DP, 2.30_DP, 2.10_DP, 1.80_DP, 1.80_DP, 1.75_DP, 1.88_DP, & ! Na, Mg, Al, Si, P, S, Cl, Ar
         & 2.75_DP, 0.00_DP, 0.00_DP, 0.00_DP, 0.00_DP, 0.00_DP, 0.00_DP, 1.63_DP, 0.00_DP, & ! K, -, -, -, -, -, -, Ni, -
         & 0.00_DP, 1.40_DP, 1.39_DP, 1.87_DP, 2.19_DP, 1.85_DP, 1.90_DP, 1.85_DP, 2.02_DP, & ! -, Cu, Zn, Ga, Ge, As, Se, Be, Kr
         & 0.00_DP, 0.00_DP, 0.00_DP, 0.00_DP, 0.00_DP, 0.00_DP, 0.00_DP, 0.00_DP, 0.00_DP, & ! -, -, -, -, -, -, -, -, -
         & 1.63_DP, 1.72_DP, 1.58_DP, 1.93_DP, 2.17_DP, 0.00_DP, 2.06_DP, 1.98_DP, 2.16_DP, & ! Pd, Ag, Cd, In, Sn, -, Te, I, Xe
         & 0.00_DP, 0.00_DP, 0.00_DP, 0.00_DP, 0.00_DP, 0.00_DP, 0.00_DP, 0.00_DP, & ! -, -, -, -, -, -, -, -
         & 0.00_DP, 0.00_DP, 0.00_DP, 0.00_DP, 0.00_DP, 0.00_DP, 0.00_DP, 0.00_DP, & ! -, -, -, -, -, -, -, -
         & 0.00_DP, 0.00_DP, 0.00_DP, 0.00_DP, 0.00_DP, 0.00_DP, 0.00_DP, 1.75_DP, & ! -, -, -, -, -, -, -, Pt
         & 1.66_DP, 1.55_DP, 1.96_DP, 2.02_DP, 0.00_DP, 0.00_DP, 0.00_DP, 0.00_DP, & ! Au, Hg, Tl, Pb, -, -, -, -
         & 0.00_DP, 0.00_DP, 0.00_DP, 0.00_DP, 0.00_DP, 1.86_DP / ! -,-,-,-,-,U
    REAL( DP ), DIMENSION(92) :: UFF_diameters
    DATA UFF_diameters/  2.886_DP, 2.362_DP, & ! H, He
         & 2.451_DP, 2.745_DP, 4.083_DP, 3.851_DP, 3.660_DP, 3.500_DP, 3.364_DP, 3.243_DP, & ! Li, Be, B, C, N, O, F, Ne
         & 2.983_DP, 3.021_DP, 4.499_DP, 4.295_DP, 4.147_DP, 4.035_DP, 3.947_DP, 3.868_DP, & ! Na, Mg, Al, Si, P, S, Cl, Ar
         & 3.812_DP, 3.399_DP, 3.295_DP, 3.175_DP, 3.144_DP, 3.023_DP, 2.961_DP, 2.912_DP, 2.872_DP, & ! K, Ca, Sc, Ti, V, Cr, Mn, Ni, Fe
         & 2.834_DP, 3.495_DP, 2.763_DP, 4.383_DP, 4.280_DP, 4.230_DP, 4.205_DP, 4.189_DP, 4.141_DP, & ! Co, Cu, Zn, Ga, Ge, As, Se, Br, Kr
         & 4.114_DP, 3.641_DP, 3.345_DP, 3.124_DP, 3.165_DP, 3.052_DP, 2.998_DP, 2.963_DP, 2.929_DP, & ! Rb, Sr, Y, Zr, Nb, Mo, Tc, Ru, Rh
         & 2.899_DP, 3.148_DP, 2.848_DP, 4.463_DP, 4.392_DP, 4.420_DP, 4.470_DP, 4.500_DP, 4.404_DP, & ! Pd, Ag, Cd, In, Sn, Sb, Te, I, Xe
         & 4.517_DP, 3.703_DP, 3.522_DP, 3.556_DP, 3.606_DP, 3.575_DP, 3.547_DP, 3.520_DP, & ! Cs, Ba, La, Ce, Pr, Nd, Pm, Sm
         & 3.493_DP, 3.368_DP, 3.451_DP, 3.428_DP, 3.409_DP, 3.391_DP, 3.374_DP, 3.355_DP, & ! Eu, Gd, Tb, Dy, Ho, Er, Tm, Yb
         & 3.640_DP, 3.141_DP, 3.170_DP, 3.069_DP, 2.954_DP, 3.120_DP, 2.840_DP, 2.754_DP, & ! Lu, Hf, Ta, W, Re, Os, Ir, Pt
         & 3.293_DP, 2.705_DP, 4.337_DP, 4.297_DP, 4.379_DP, 4.709_DP, 4.750_DP, 4.765_DP, & ! Au, Hg, Tl, Pb, Bi, Po, At, Rn
         & 4.900_DP, 3.677_DP, 3.478_DP, 3.396_DP, 3.424_DP, 3.395_DP / ! Fr, Ra, Ac, Th, Pa, U

    iontype%index = index
    iontype%label = label
    iontype%zv = 0.D0 ! this cannot be initialized here at this time
    iontype%atmnum = get_atmnum(label)

    IF ( iontype%atmnum .EQ. 0 ) &
         & CALL errore(sub_name,'Can not assign the atom type associated with input label',1)

    iontype%atomicspread = 0.5D0
    iontype%corespread = 0.5D0

    SELECT CASE ( radius_mode )

    CASE ( 'pauling' )

       iontype%solvationrad = pauling_radii(iontype%atmnum)

    CASE ( 'bondi' )

       iontype%solvationrad = bondi_radii(iontype%atmnum)

    CASE ( 'uff' )

       iontype%solvationrad = UFF_diameters(iontype%atmnum) * 0.5_DP

    CASE DEFAULT

       CALL errore(sub_name,'Unknown radius_mode',1)

    END SELECT

    RETURN

  END SUBROUTINE set_iontype_defaults

!--------------------------------------------------------------------
  FUNCTION get_atmnum(label)
!wgt  FUNCTION get_atmwgt(label)
!
! original version by O. Andreussi (MIT)
!
!--------------------------------------------------------------------

    INTEGER :: get_atmnum
!wgt    REAL*8 :: get_atmwgt

    CHARACTER*(*), INTENT(IN) :: label
    CHARACTER*3 :: tmplab

    INTEGER :: num
    REAL*8 :: weigth

    tmplab=TRIM(ADJUSTL(label))
    CALL lowcase(tmplab)
    IF (tmplab(1:1).EQ.'a')  THEN
      IF (tmplab(2:2).EQ.'c') THEN
        num=89
        weigth=227.03
      ELSE IF (tmplab(2:2).EQ.'l') THEN
        num=13
        weigth=26.981538
      ELSE IF (tmplab(2:2).EQ.'r') THEN
        num=18
        weigth=39.948
      ELSE IF (tmplab(2:2).EQ.'g') THEN
        num=47
        weigth=107.8682
      ELSE IF (tmplab(2:2).EQ.'s') THEN
        num=33
        weigth=74.9216
      ELSE IF (tmplab(2:2).EQ.'u') THEN
        num=79
        weigth=196.96655
      ELSE IF (tmplab(2:2).EQ.'t') THEN
        num=85
        weigth=210
      ELSE
        num=13
        weigth=26.981538
      ENDIF
    ELSE IF (tmplab(1:1).EQ.'b')  THEN
      IF (tmplab(2:2).EQ.'e') THEN
        num=4
        weigth=9.012182
      ELSE IF (tmplab(2:2).EQ.'r') THEN
        num=35
        weigth=79.904
      ELSE IF (tmplab(2:2).EQ.'a') THEN
        num=56
        weigth=137.327
      ELSE IF (tmplab(2:2).EQ.'i') THEN
        num=83
        weigth=208.98038
      ELSE
        num=5
        weigth=10.811
      ENDIF
    ELSE IF (tmplab(1:1).EQ.'c')  THEN
      IF (tmplab(2:2).EQ.'a') THEN
        num=20
        weigth=40.078
      ELSE IF (tmplab(2:2).EQ.'r') THEN
        num=24
        weigth=51.9961
      ELSE IF (tmplab(2:2).EQ.'o') THEN
        num=27
        weigth=58.9332
      ELSE IF (tmplab(2:2).EQ.'l') THEN
        num=17
        weigth=35.453
      ELSE IF (tmplab(2:2).EQ.'s') THEN
        num=55
        weigth=132.90545
      ELSE IF (tmplab(2:2).EQ.'d') THEN
        num=48
        weigth=112.411
      ELSE IF (tmplab(2:2).EQ.'u') THEN
        num=29
        weigth=63.546
      ELSE IF (tmplab(2:2).EQ.'e') THEN
        num=58
        weigth=140.116
      ELSE
        num=6
        weigth=12.0107
      ENDIF
    ELSE IF (tmplab(1:1).EQ.'f')  THEN
      IF (tmplab(2:2).EQ.'e') THEN
        num=26
        weigth=55.845
      ELSE IF (tmplab(2:2).EQ.'r') THEN
        num=87
        weigth=223
      ELSE
        num=9
        weigth=18.9984032
      ENDIF
    ELSE IF (tmplab(1:1).EQ.'g')  THEN
      IF (tmplab(2:2).EQ.'a') THEN
        num=31
        weigth=69.723
      ELSE IF (tmplab(2:2).EQ.'e') THEN
        num=32
        weigth=72.64
      ELSE
        num=31
        weigth=69.723
      ENDIF
    ELSE IF (tmplab(1:1).EQ.'h')  THEN
      IF (tmplab(2:2).EQ.'e') THEN
        num=2
        weigth=4.002602
      ELSE IF (tmplab(2:2).EQ.'g') THEN
        num=80
        weigth=200.59
      ELSE IF (tmplab(2:2).EQ.'f') THEN
        num=72
        weigth=178.49
      ELSE
        num=1
        weigth=1.00794
      ENDIF
    ELSE IF (tmplab(1:1).EQ.'i')  THEN
      IF (tmplab(2:2).EQ.'n') THEN
        num=49
        weigth=114.818
      ELSE IF (tmplab(2:2).EQ.'r') THEN
        num=77
        weigth=192.217
      ELSE
        num=53
        weigth=126.90447
      ENDIF
    ELSE IF (tmplab(1:1).EQ.'k')  THEN
      IF (tmplab(2:2).EQ.'r') THEN
        num=36
        weigth=83.798
      ELSE
        num=19
        weigth=39.0983
      ENDIF
    ELSE IF (tmplab(1:1).EQ.'l')  THEN
      IF (tmplab(2:2).EQ.'i') THEN
        num=3
        weigth=6.941
      ELSE IF (tmplab(2:2).EQ.'a') THEN
        num=57
        weigth=138.9055
      ELSE
        num=3
        weigth=6.941
      ENDIF
    ELSE IF (tmplab(1:1).EQ.'m')  THEN
      IF (tmplab(2:2).EQ.'g') THEN
        num=12
        weigth=24.3050
      ELSE IF (tmplab(2:2).EQ.'n') THEN
        num=25
        weigth=54.938049
      ELSE IF (tmplab(2:2).EQ.'o') THEN
        num=42
        weigth=95.94
      ELSE
        num=25
        weigth=54.938049
      ENDIF
    ELSE IF (tmplab(1:1).EQ.'n')  THEN
      IF (tmplab(2:2).EQ.'a') THEN
        num=11
        weigth=22.98977
      ELSE IF (tmplab(2:2).EQ.'i') THEN
        num=28
        weigth=58.6934
      ELSE IF (tmplab(2:2).EQ.'b') THEN
        num=41
        weigth=92.90638
      ELSE IF (tmplab(2:2).EQ.'e') THEN
        num=10
        weigth=20.1797
      ELSE IF (tmplab(2:2).EQ.'d') THEN
        num=60
        weigth=144.24
      ELSE
        num=7
        weigth=14.0067
      ENDIF
    ELSE IF (tmplab(1:1).EQ.'o')  THEN
      IF (tmplab(2:2).EQ.'s') THEN
        num=76
        weigth=190.23
      ELSE
        num=8
        weigth=15.9994
      ENDIF
    ELSE IF (tmplab(1:1).EQ.'p')  THEN
      IF (tmplab(2:2).EQ.'a') THEN
        num=91
        weigth=231.04
      ELSE IF (tmplab(2:2).EQ.'d') THEN
        num=46
        weigth=106.42
      ELSE IF (tmplab(2:2).EQ.'t') THEN
        num=78
        weigth=195.078
      ELSE IF (tmplab(2:2).EQ.'b') THEN
        num=82
        weigth=207.2
      ELSE IF (tmplab(2:2).EQ.'o') THEN
        num=84
        weigth=209
      ELSE
        num=15
        weigth=30.973761
      ENDIF
    ELSE IF (tmplab(1:1).EQ.'r')  THEN
      IF (tmplab(2:2).EQ.'b') THEN
        num=37
        weigth=85.4678
      ELSE IF (tmplab(2:2).EQ.'u') THEN
        num=44
        weigth=101.07
      ELSE IF (tmplab(2:2).EQ.'h') THEN
        num=45
        weigth=102.90550
      ELSE IF (tmplab(2:2).EQ.'e') THEN
        num=75
        weigth=186.207
      ELSE IF (tmplab(2:2).EQ.'n') THEN
        num=86
        weigth=222
      ELSE
        num=37
        weigth=85.4678
      ENDIF
    ELSE IF (tmplab(1:1).EQ.'s')  THEN
      IF (tmplab(2:2).EQ.'i') THEN
        num=14
        weigth=28.0855
      ELSE IF (tmplab(2:2).EQ.'e') THEN
        num=34
        weigth=78.96
      ELSE IF (tmplab(2:2).EQ.'c') THEN
        num=21
        weigth=44.955910
      ELSE IF (tmplab(2:2).EQ.'r') THEN
        num=38
        weigth=87.62
      ELSE IF (tmplab(2:2).EQ.'n') THEN
        num=50
        weigth=118.710
      ELSE IF (tmplab(2:2).EQ.'b') THEN
        num=51
        weigth=121.760
      ELSE
        num=16
        weigth=32.065
      ENDIF
    ELSE IF (tmplab(1:1).EQ.'t')  THEN
      IF (tmplab(2:2).EQ.'a') THEN
        num=73
        weigth=180.9479
      ELSE IF (tmplab(2:2).EQ.'l') THEN
        num=81
        weigth=204.3833
      ELSE IF (tmplab(2:2).EQ.'c') THEN
        num=43
        weigth=98
      ELSE IF (tmplab(2:2).EQ.'h') THEN
        num=90
        weigth=232.04
      ELSE IF (tmplab(2:2).EQ.'i') THEN
        num=22
        weigth=47.867
      ELSE IF (tmplab(2:2).EQ.'e') THEN
        num=52
        weigth=127.60
      ELSE
        num=22
        weigth=47.867
      ENDIF
    ELSE IF (tmplab(1:1).EQ.'u')  THEN
      num=92
      weigth=238.02891
    ELSE IF (tmplab(1:1).EQ.'v')  THEN
      num=23
      weigth=50.9415
    ELSE IF (tmplab(1:1).EQ.'w')  THEN
      num=74
      weigth=183.84
    ELSE IF (tmplab(1:1).EQ.'x')  THEN
      num=54
      weigth=131.293
    ELSE IF (tmplab(1:1).EQ.'y')  THEN
      num=39
      weigth=88.90585
    ELSE IF (tmplab(1:1).EQ.'z')  THEN
      IF (tmplab(2:2).EQ.'n') THEN
        num=30
        weigth=65.409
      ELSE IF (tmplab(2:2).EQ.'r') THEN
        num=40
        weigth=91.224
      ELSE
        num=30
        weigth=65.409
      ENDIF
    ELSE
      num=0
      weigth=0
    ENDIF

    get_atmnum=num
!    get_atmwgt=weigth

!--------------------------------------------------------------------
  END FUNCTION get_atmnum
!wgt  END FUNCTION get_atmwgt
!--------------------------------------------------------------------

!--------------------------------------------------------------------
  SUBROUTINE lowcase(string)
!--------------------------------------------------------------------

    CHARACTER*(*), INTENT(inout) :: string

    INTEGER :: i, length
    CHARACTER(1) :: letter

    length=LEN(string)

    DO i = 1,min(255,length)
      letter = string(i:i)
      IF(letter.eq.'A') THEN
        letter = 'a'
      ELSE IF(letter.eq.'B') THEN
        letter = 'b'
      ELSE IF(letter.eq.'C') THEN
        letter = 'c'
      ELSE IF(letter.eq.'D') THEN
        letter = 'd'
      ELSE IF(letter.eq.'E') THEN
        letter = 'e'
      ELSE IF(letter.eq.'F') THEN
        letter = 'f'
      ELSE IF(letter.eq.'G') THEN
        letter = 'g'
      ELSE IF(letter.eq.'H') THEN
        letter = 'h'
      ELSE IF(letter.eq.'I') THEN
        letter = 'i'
      ELSE IF(letter.eq.'J') THEN
        letter = 'j'
      ELSE IF(letter.eq.'K') THEN
        letter = 'k'
      ELSE IF(letter.eq.'L') THEN
        letter = 'l'
      ELSE IF(letter.eq.'M') THEN
        letter = 'm'
      ELSE IF(letter.eq.'N') THEN
        letter = 'n'
      ELSE IF(letter.eq.'O') THEN
        letter = 'o'
      ELSE IF(letter.eq.'P') THEN
        letter = 'p'
      ELSE IF(letter.eq.'Q') THEN
        letter = 'q'
      ELSE IF(letter.eq.'R') THEN
        letter = 'r'
      ELSE IF(letter.eq.'S') THEN
        letter = 's'
      ELSE IF(letter.eq.'T') THEN
        letter = 't'
      ELSE IF(letter.eq.'U') THEN
        letter = 'u'
      ELSE IF(letter.eq.'V') THEN
        letter = 'v'
      ELSE IF(letter.eq.'W') THEN
        letter = 'w'
      ELSE IF(letter.eq.'X') THEN
        letter = 'x'
      ELSE IF(letter.eq.'Y') THEN
        letter = 'y'
      ELSE IF(letter.eq.'Z') THEN
        letter = 'z'
      END IF
      string(i:i) = letter
    END DO

    RETURN
!--------------------------------------------------------------------
  END SUBROUTINE lowcase
!--------------------------------------------------------------------
END MODULE environ_types
