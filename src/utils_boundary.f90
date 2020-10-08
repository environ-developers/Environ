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
!> Module containing the main routines to handle environ_boundary
!! derived data types.
!!
!! Environ_boundary contains all the specifications and the details of
!! the smooth interface between the QM and the continuum regions of the
!! simulation cell. The main interface function is stored in the %scaled
!! component, the type also stores boundary real-space derivatives (gradient,
!! laplacian, dsurface, hessian) and other quantities needed by Environ
!! modules.
!
!----------------------------------------------------------------------------
!  TYPE environ_boundary
!----------------------------------------------------------------------------
!
!     ! Boundary label
!
!     CHARACTER (LEN=80) :: label
!
!     ! Choice of the interface
!
!     CHARACTER (LEN=80) :: mode
!
!     ! Update status
!
!     INTEGER :: update_status = 0
!
!     LOGICAL :: initialized = .FALSE.
!
!     ! Parameters for the electrons-dependent interface
!
!     LOGICAL :: need_electrons
!     TYPE( environ_electrons ), POINTER :: electrons
!
!     ! Parameters for the ions-dependent interface
!
!     LOGICAL :: need_ions
!     TYPE( environ_ions ), POINTER :: ions
!
!     ! Parameters for the system-dependent interface
!
!     LOGICAL :: need_system
!     TYPE( environ_system ), POINTER :: system
!
!     ! scaled switching function of interface
!     ! varying from 1 (QM region) to 0 (environment region)
!
!     TYPE( environ_density ) :: scaled
!
!     INTEGER :: deriv = 0
!     TYPE( environ_gradient ) :: gradient
!     TYPE( environ_density ) :: laplacian
!     TYPE( environ_density ) :: dsurface
!     TYPE( environ_hessian ) :: hessian
!
!     ! global properties of the boundary
!
!     REAL( DP ) :: volume
!     REAL( DP ) :: surface
!
!     ! Components needed for boundary of density
!
!     INTEGER :: type
!     REAL( DP ) :: rhomax, rhomin, fact
!     REAL( DP ) :: rhozero, deltarho, tbeta
!     REAL( DP ) :: const
!     TYPE( environ_density ) :: density
!
!     TYPE( environ_density ) :: dscaled
!     TYPE( environ_density ) :: d2scaled
!
!     ! Components needed for boundary of functions
!
!     REAL( DP ) :: alpha ! solvent-dependent scaling factor
!     REAL( DP ) :: softness ! sharpness of the interface
!     TYPE( environ_functions ), DIMENSION(:), ALLOCATABLE :: soft_spheres
!
!     ! Components needed for boundary of system
!
!     TYPE( environ_functions ) :: simple
!
!     ! Copmonents needed for solvent-aware boundary
!
!     LOGICAL :: solvent_aware
!     TYPE( environ_functions ) :: solvent_probe
!     REAL( DP ) :: filling_threshold, filling_spread
!
!     TYPE( environ_density ) :: local
!     TYPE( environ_density ) :: probe
!     TYPE( environ_density ) :: filling
!     TYPE( environ_density ) :: dfilling
!
!----------------------------------------------------------------------------
!  END TYPE environ_boundary
!----------------------------------------------------------------------------
!
! Authors: Oliviero Andreussi (Department of Physics, UNT)
!          Francesco Nattino  (THEOS and NCCR-MARVEL, EPFL)
!          Nicola Marzari     (THEOS and NCCR-MARVEL, EPFL)
!
!----------------------------------------------------------------------------
MODULE utils_boundary
!----------------------------------------------------------------------------
  !
  USE core_types
  USE environ_types
  USE environ_output
  USE utils_functions
  USE environ_base, ONLY : niter
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC :: create_boundary_core, init_boundary_core, &
       & destroy_boundary_core
  PUBLIC :: create_environ_boundary, init_environ_boundary_first, &
       & init_environ_boundary_second, copy_environ_boundary, &
       & set_soft_spheres, update_environ_boundary, &
       & destroy_environ_boundary
  !
CONTAINS
!--------------------------------------------------------------------
  SUBROUTINE create_boundary_core( core )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    TYPE( boundary_core ), INTENT(INOUT) :: core
    !
    core % type = 'default'
    core % use_fft = .FALSE.
    NULLIFY( core % fft )
    core % use_fd = .FALSE.
    NULLIFY( core % fd )
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE create_boundary_core
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE init_boundary_core( type, core, fft, fd )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    CHARACTER( LEN = 80 ), INTENT(IN) :: type
    TYPE( boundary_core ), INTENT(INOUT) :: core
    TYPE( fft_core ), INTENT(IN), TARGET, OPTIONAL :: fft
    TYPE( fd_core ), INTENT(IN), TARGET, OPTIONAL :: fd
    !
    INTEGER :: number
    CHARACTER( LEN = 80 ) :: sub_name = 'init_boundary_core'
    !
    core % type = type
    !
    ! Assign the selected numerical core
    !
    SELECT CASE ( TRIM( ADJUSTL( type ) ) )
       !
    CASE ( 'analytic', 'fft', 'default' )
       !
       IF ( .NOT. PRESENT( fft ) ) CALL errore(sub_name,'Missing specified core type',1)
       core % use_fft = .TRUE.
       core % fft => fft
       !
    CASE ( 'fd', 'finite differences', 'finite_differences' )
       !
       ! Note: finite differences core only works for gradient, other derivatives still
       ! require fft core
       !
       IF ( .NOT. PRESENT( fd ) ) CALL errore(sub_name,'Missing specified core type',1)
       IF ( .NOT. PRESENT( fft ) ) CALL errore(sub_name,'Missing specified core type',1)
       core % use_fd = .TRUE.
       core % use_fft = .TRUE.
       core % fft => fft
       core % fd => fd
       !
    CASE DEFAULT
       !
       CALL errore(sub_name,'Unexpected keyword for boundary core type',1)
       !
    END SELECT
    !
    ! double check number of active cores
    !
    number = 0
    IF ( core % use_fft ) number = number + 1
    IF ( core % use_fd ) number = number + 1
    IF ( number .LT. 1 ) &
         & CALL errore(sub_name,'Too few cores are active',1)
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE init_boundary_core
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE destroy_boundary_core( lflag, core )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    LOGICAL, INTENT(IN) :: lflag
    TYPE( boundary_core ), INTENT(INOUT) :: core
    !
    IF ( lflag ) THEN
       core % use_fft = .FALSE.
       NULLIFY( core % fft )
       core % use_fd = .FALSE.
       NULLIFY( core % fd )
    END IF
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE destroy_boundary_core
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE create_environ_boundary(boundary, local_label)
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    TYPE( environ_boundary ), INTENT(INOUT) :: boundary
    CHARACTER( LEN=80 ), INTENT(IN) :: local_label
    !
    CHARACTER( LEN=80 ) :: sub_name = 'create_environ_boundary'
    CHARACTER( LEN=80 ) :: label = ' '
    !
    boundary%update_status = 0
    boundary%label = local_label
    !
    label = 'boundary_'//TRIM(ADJUSTL(local_label))
    CALL create_environ_density( boundary%scaled, label )
    boundary%volume = 0.D0
    boundary%surface = 0.D0
    !
    boundary%need_electrons = .FALSE.
    NULLIFY( boundary%electrons )
    boundary%need_ions = .FALSE.
    NULLIFY( boundary%ions )
    boundary%need_system = .FALSE.
    NULLIFY( boundary%system )
    !
    ! Optional components
    !
    boundary%deriv = 0
    label = 'gradboundary_'//TRIM(ADJUSTL(local_label))
    CALL create_environ_gradient( boundary%gradient, label )
    label = 'laplboundary_'//TRIM(ADJUSTL(local_label))
    CALL create_environ_density( boundary%laplacian, label )
    label = 'dsurface_'//TRIM(ADJUSTL(local_label))
    CALL create_environ_density( boundary%dsurface, label )
    label = 'hessboundary_'//TRIM(ADJUSTL(local_label))
    CALL create_environ_hessian( boundary%hessian, label )
    !
    ! Components required for boundary of density
    !
    label = 'boundary_density_'//TRIM(ADJUSTL(local_label))
    CALL create_environ_density( boundary%density, label )
    label = 'dboundary_'//TRIM(ADJUSTL(local_label))
    CALL create_environ_density( boundary%dscaled, label )
    label = 'd2boundary_'//TRIM(ADJUSTL(local_label))
    CALL create_environ_density( boundary%d2scaled, label )
    !
    ! Components required for boundary of functions
    !
    IF ( ALLOCATED( boundary%soft_spheres ) ) &
         & CALL errore(sub_name,'Trying to create an already allocated object',1)
    !
    ! Components required for solvent-aware interface
    !
    boundary%solvent_aware = .FALSE.
    label = 'local_'//TRIM(ADJUSTL(local_label))
    CALL create_environ_density( boundary%local, label )
    label = 'probe_'//TRIM(ADJUSTL(local_label))
    CALL create_environ_density( boundary%probe, label )
    label = 'filling_'//TRIM(ADJUSTL(local_label))
    CALL create_environ_density( boundary%filling, label )
    label = 'dfilling_'//TRIM(ADJUSTL(local_label))
    CALL create_environ_density( boundary%dfilling, label )
    !
    ! Components required for field-aware interface
    !
    boundary%field_aware = .FALSE.
    label = 'normal_field_'//TRIM(ADJUSTL(local_label))
    CALL create_environ_density( boundary%normal_field, label )
    IF ( ALLOCATED( boundary%ion_field ) ) &
         & CALL errore(sub_name,'Trying to create an already allocated object',1)
    IF ( ALLOCATED( boundary%local_spheres ) ) &
         & CALL errore(sub_name,'Trying to create an already allocated object',1)
    IF ( ALLOCATED( boundary%dion_field_drho ) ) &
         & CALL errore(sub_name,'Trying to create an already allocated object',1)
    IF ( ALLOCATED( boundary%partial_of_ion_field ) ) &
         & CALL errore(sub_name,'Trying to create an already allocated object',1)
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE create_environ_boundary
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE init_environ_boundary_first( need_gradient, need_laplacian, &
       & need_hessian, mode, stype, rhomax, rhomin, tbeta, const, alpha, &
       & softness, system_distance, system_spread, solvent_radius, radial_scale, &
       & radial_spread, filling_threshold, filling_spread, field_factor, &
       & charge_asymmetry, field_max, field_min, electrons, ions, system, &
       & core, boundary )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    CHARACTER( LEN=80 ), INTENT(IN) :: mode
    INTEGER, INTENT(IN) :: stype
    REAL( DP ), INTENT(IN) :: rhomax, rhomin, tbeta, const
    LOGICAL, INTENT(IN) :: need_gradient, need_laplacian, need_hessian
    REAL( DP ), INTENT(IN) :: alpha
    REAL( DP ), INTENT(IN) :: softness
    REAL( DP ), INTENT(IN) :: system_distance
    REAL( DP ), INTENT(IN) :: system_spread
    REAL( DP ), INTENT(IN) :: solvent_radius
    REAL( DP ), INTENT(IN) :: radial_scale, radial_spread
    REAL( DP ), INTENT(IN) :: filling_threshold, filling_spread
    REAL( DP ), INTENT(IN) :: field_factor, charge_asymmetry, field_max, field_min
    TYPE( environ_electrons ), TARGET, INTENT(IN) :: electrons
    TYPE( environ_ions ), TARGET, INTENT(IN) :: ions
    TYPE( environ_system ), TARGET, INTENT(IN) :: system
    TYPE( boundary_core ), TARGET, INTENT(IN) :: core
    TYPE( environ_boundary ), INTENT(INOUT) :: boundary
    CHARACTER( LEN=80 ) :: sub_name = 'init_environ_boundary_first'
    !
    IF ( need_hessian ) THEN
       boundary%deriv = 3
    ELSE IF ( need_laplacian ) THEN
       boundary%deriv = 2
    ELSE IF ( need_gradient ) THEN
       boundary%deriv = 1
    ENDIF
    !
    boundary%mode = mode
    !
    boundary%need_electrons = ( mode .EQ. 'electronic' ) .OR. ( mode .EQ. 'full' ) &
      & .OR. ( mode .EQ. 'fa-ionic' ) .OR. ( mode .EQ. 'fa-electronic' )
    IF ( boundary%need_electrons ) boundary%electrons => electrons
    boundary%need_ions = ( mode .EQ. 'ionic' ) .OR. ( mode .EQ. 'full' ) &
      & .OR. ( mode .EQ. 'fa-ionic' ) .OR. ( mode .EQ. 'fa-electronic' )
    IF ( boundary%need_ions ) boundary%ions => ions
    boundary%need_system = ( mode .EQ. 'system' )
    IF ( boundary%need_system ) boundary%system => system
    !
    boundary%type = stype
    boundary%rhomax = rhomax
    boundary%rhomin = rhomin
    boundary%fact = LOG( rhomax / rhomin )
    boundary%rhozero = ( rhomax + rhomin ) * 0.5_DP
    boundary%tbeta = tbeta
    boundary%deltarho = rhomax - rhomin
    !
    IF ( const .EQ. 1.D0 .AND. boundary%need_electrons .AND. stype .EQ. 2 ) &
     & CALL errore(sub_name,'stype=2 boundary requires dielectric constant > 1',1)
    boundary%const = const
    !
    boundary%alpha = alpha
    boundary%softness = softness
    IF ( boundary%mode .EQ. 'ionic' .OR. boundary%mode .EQ. 'fa-ionic' ) ALLOCATE( boundary%soft_spheres( boundary%ions%number ) )
    !
    boundary%simple%type = 4
    boundary%simple%pos => system%pos
    boundary%simple%volume = 1.D0
    boundary%simple%dim = system%dim
    boundary%simple%axis = system%axis
    boundary%simple%width = system_distance
    boundary%simple%spread = system_spread
    !
    boundary%solvent_aware = solvent_radius .GT. 0.D0
    !
    IF( boundary%solvent_aware ) THEN
       boundary%solvent_probe%type = 2
       ALLOCATE(boundary%solvent_probe%pos(3))
       boundary%solvent_probe%pos = 0.D0
       boundary%solvent_probe%volume = 1.D0
       boundary%solvent_probe%dim = 0
       boundary%solvent_probe%axis = 1
       boundary%solvent_probe%spread = radial_spread
       boundary%solvent_probe%width = solvent_radius * radial_scale
    ENDIF
    !
    boundary%filling_threshold = filling_threshold
    boundary%filling_spread = filling_spread
    !
    boundary%core => core
    boundary%field_aware = field_factor .GT. 0.D0
    boundary%field_factor = field_factor
    boundary%charge_asymmetry = charge_asymmetry
    boundary%field_max = field_max
    boundary%field_min = field_min
    IF ( boundary%field_aware .AND. boundary%mode .EQ. 'fa-ionic' ) THEN
       ALLOCATE( boundary%ion_field( boundary%ions%number ) )
       ALLOCATE( boundary%dion_field_drho( boundary%ions%number ) )
       ALLOCATE( boundary%partial_of_ion_field( 3, boundary%ions%number, boundary%ions%number ) )
       ALLOCATE( boundary%local_spheres( boundary%ions%number ) )
     ENDIF
    !
    boundary%initialized = .FALSE.
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE init_environ_boundary_first
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE init_environ_boundary_second( cell, boundary )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    TYPE( environ_cell ), INTENT(IN) :: cell
    TYPE( environ_boundary ), INTENT(INOUT) :: boundary
    !
    CHARACTER( LEN=80 ) :: sub_name = 'init_environ_boundary_second'
    INTEGER :: i
    !
    CALL init_environ_density( cell, boundary%scaled )
    !
    IF ( boundary%mode .EQ. 'electronic' .OR. boundary%mode .EQ. 'full' .OR. &
         & boundary%mode .EQ. 'fa-electronic' .OR. boundary%mode .EQ. 'fa-full' ) THEN
       CALL init_environ_density( cell, boundary%density )
       CALL init_environ_density( cell, boundary%dscaled )
       CALL init_environ_density( cell, boundary%d2scaled )
    END IF
    IF ( boundary%deriv .GE. 1 ) CALL init_environ_gradient( cell, boundary%gradient )
    IF ( boundary%deriv .GE. 2 ) CALL init_environ_density( cell, boundary%laplacian )
    IF ( boundary%deriv .GE. 3 ) CALL init_environ_density( cell, boundary%dsurface )
    !
    IF ( boundary%solvent_aware ) THEN
       CALL init_environ_density( cell, boundary%local )
       CALL init_environ_density( cell, boundary%probe )
       CALL init_environ_density( cell, boundary%filling )
       CALL init_environ_density( cell, boundary%dfilling )
       IF ( boundary%deriv .GE. 3 ) CALL init_environ_hessian( cell, boundary%hessian )
    ENDIF
    !
    IF ( boundary%field_aware ) THEN
       IF ( boundary%mode .EQ. 'fa-electronic' .OR. boundary%mode .EQ. 'fa-full' ) THEN
          CALL init_environ_density( cell, boundary%normal_field )
       ELSE IF ( boundary%mode .EQ. 'fa-ionic' ) THEN
          DO i = 1, boundary%ions%number
             CALL init_environ_density( cell, boundary%dion_field_drho(i) )
          ENDDO
       ELSE
          CALL errore(sub_name,'boundary must be field-aware',1)
       ENDIF
    ENDIF
    !
    boundary%initialized = .TRUE.
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE init_environ_boundary_second
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE copy_environ_boundary( boriginal, bcopy )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    TYPE( environ_boundary ), INTENT(IN) :: boriginal
    TYPE( environ_boundary ), INTENT(OUT) :: bcopy
    !
    INTEGER :: i, n, m
    !
    bcopy % electrons => boriginal % electrons
    bcopy % ions      => boriginal % ions
    bcopy % system    => boriginal % system
    bcopy % core      => boriginal % core
    !
    bcopy % mode              = boriginal % mode
    bcopy % update_status     = boriginal % update_status
    bcopy % need_electrons    = boriginal % need_electrons
    bcopy % need_ions         = boriginal % need_ions
    bcopy % need_system       = boriginal % need_system
    bcopy % deriv             = boriginal % deriv
    bcopy % volume            = boriginal % volume
    bcopy % surface           = boriginal % surface
    bcopy % type              = boriginal % type
    bcopy % rhomax            = boriginal % rhomax
    bcopy % rhomin            = boriginal % rhomin
    bcopy % fact              = boriginal % fact
    bcopy % rhozero           = boriginal % rhozero
    bcopy % deltarho          = boriginal % deltarho
    bcopy % tbeta             = boriginal % tbeta
    bcopy % const             = boriginal % const
    bcopy % alpha             = boriginal % alpha
    bcopy % softness          = boriginal % softness
    bcopy % solvent_aware     = boriginal % solvent_aware
    bcopy % filling_threshold = boriginal % filling_threshold
    bcopy % filling_spread    = boriginal % filling_spread
    bcopy % field_aware       = boriginal % field_aware
    bcopy % field_factor      = boriginal % field_factor
    bcopy % charge_asymmetry  = boriginal % charge_asymmetry
    bcopy % field_max         = boriginal % field_max
    bcopy % field_min         = boriginal % field_min
    bcopy % initialized       = boriginal % initialized
    !
    IF ( ASSOCIATED( boriginal%scaled%cell ) ) CALL copy_environ_density ( boriginal % scaled, bcopy % scaled )
    IF ( ASSOCIATED( boriginal%gradient%cell ) ) CALL copy_environ_gradient ( boriginal % gradient, bcopy % gradient )
    IF ( ASSOCIATED( boriginal%laplacian%cell ) ) CALL copy_environ_density ( boriginal % laplacian, bcopy % laplacian )
    IF ( ASSOCIATED( boriginal%dsurface%cell ) ) CALL copy_environ_density ( boriginal % dsurface, bcopy % dsurface )
    IF ( ASSOCIATED( boriginal%hessian%cell ) ) CALL copy_environ_hessian ( boriginal % hessian, bcopy % hessian )
    IF ( ASSOCIATED( boriginal%density%cell ) ) CALL copy_environ_density ( boriginal % density, bcopy % density )
    IF ( ASSOCIATED( boriginal%dscaled%cell ) ) CALL copy_environ_density ( boriginal % dscaled, bcopy % dscaled )
    IF ( ASSOCIATED( boriginal%d2scaled%cell ) ) CALL copy_environ_density ( boriginal % d2scaled, bcopy % d2scaled )
    CALL copy_environ_functions ( boriginal % simple        , bcopy % simple        )
    CALL copy_environ_functions ( boriginal % solvent_probe , bcopy % solvent_probe )
    IF ( ASSOCIATED( boriginal%local%cell ) ) CALL copy_environ_density ( boriginal % local, bcopy % local )
    IF ( ASSOCIATED( boriginal%probe%cell ) ) CALL copy_environ_density ( boriginal % probe, bcopy % probe )
    IF ( ASSOCIATED( boriginal%filling%cell ) ) CALL copy_environ_density ( boriginal % filling, bcopy % filling )
    IF ( ASSOCIATED( boriginal%dfilling%cell ) ) CALL copy_environ_density ( boriginal % dfilling, bcopy % dfilling )
    IF ( ASSOCIATED( boriginal%normal_field%cell ) ) CALL copy_environ_density ( boriginal % normal_field, bcopy % normal_field )
    !
    IF ( ALLOCATED( boriginal % soft_spheres ) ) THEN
       n = SIZE( boriginal % soft_spheres )
       IF ( ALLOCATED( bcopy % soft_spheres ) ) THEN
          m = SIZE( bcopy % soft_spheres )
          CALL destroy_environ_functions( m, bcopy % soft_spheres )
       ENDIF
       ALLOCATE( bcopy % soft_spheres( n ) )
       DO i = 1, n
          CALL copy_environ_functions ( boriginal % soft_spheres(i), bcopy % soft_spheres(i) )
       ENDDO
    ELSE
       IF ( ALLOCATED( bcopy % soft_spheres ) ) DEALLOCATE( bcopy%soft_spheres )
    ENDIF
    !
    IF ( ALLOCATED( boriginal % ion_field ) ) THEN
       n = SIZE( boriginal % ion_field )
       IF ( ALLOCATED( bcopy % ion_field ) ) DEALLOCATE( bcopy % ion_field )
       IF ( ALLOCATED( bcopy % partial_of_ion_field ) ) DEALLOCATE( bcopy % partial_of_ion_field )
       ALLOCATE( bcopy % ion_field( n ) )
       ALLOCATE( bcopy % partial_of_ion_field( 3, n, n ) )
       IF ( ALLOCATED( bcopy % dion_field_drho ) ) THEN
          m = SIZE( bcopy % dion_field_drho )
          DO i = 1, m
             CALL destroy_environ_density( bcopy % dion_field_drho( i ) )
          ENDDO
          DEALLOCATE( bcopy % dion_field_drho )
       ENDIF
       ALLOCATE( bcopy % dion_field_drho( n ) )
       DO i = 1, n
          CALL copy_environ_density( boriginal % dion_field_drho(i), bcopy % dion_field_drho(i) )
       ENDDO
    ELSE
       IF ( ALLOCATED( bcopy % ion_field ) ) DEALLOCATE( bcopy % ion_field )
       IF ( ALLOCATED( bcopy % partial_of_ion_field ) ) DEALLOCATE( bcopy % partial_of_ion_field )
       IF ( ALLOCATED( bcopy % dion_field_drho ) ) DEALLOCATE( bcopy % dion_field_drho )
    ENDIF
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE copy_environ_boundary
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE set_soft_spheres( boundary, scale )
!--------------------------------------------------------------------
    !
    USE tools_generate_boundary, ONLY : scaling_of_field
    !
    IMPLICIT NONE
    !
    TYPE( environ_boundary ), INTENT(INOUT) :: boundary
    LOGICAL, INTENT(IN), OPTIONAL :: scale
    !
    LOGICAL :: lscale1, lscale2
    INTEGER :: i
    REAL( DP ) :: radius, f
    !
    IF ( .NOT. ( boundary % mode .EQ. 'ionic' .OR. boundary % mode .EQ. 'fa-ionic' ) ) RETURN
    !
    lscale1 = .FALSE.
    lscale2 = .FALSE.
    IF ( PRESENT( scale ) ) THEN
       lscale1 = scale .AND. ( boundary % mode .EQ. 'fa-ionic' )
    ELSE
       lscale2 = ( boundary % mode .EQ. 'fa-ionic' )
    ENDIF
    !
    f = 1.D0
    DO i = 1, boundary%ions%number
       IF ( lscale1 ) f = scaling_of_field(boundary%field_factor,boundary%charge_asymmetry,&
            & boundary%field_max,boundary%field_min,boundary%ion_field(i))
       radius = boundary%ions%iontype(boundary%ions%ityp(i))%solvationrad * boundary%alpha * f
       boundary%soft_spheres(i) = environ_functions(5,1,0,radius,boundary%softness,1.D0,&
            & boundary%ions%tau(:,i))
       IF ( lscale2 ) boundary%local_spheres(i) = boundary%soft_spheres(i)
       IF ( lscale1 .AND. verbose .GE. 1 ) WRITE(environ_unit,6100)i,boundary%ions%iontype(boundary%ions%ityp(i))%label,&
               & boundary%ions%iontype(boundary%ions%ityp(i))%solvationrad,boundary%alpha,&
               & boundary%ion_field(i),f,radius
6100   FORMAT("atom numer = ",i3," atom label = ",a3, &
    & " solvation radius = ",f8.4," scaling = ",f8.4," field flux = ", &
    & f8.4," scaling of field = ",f8.4," final radius = ",f8.4)
    ENDDO
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE set_soft_spheres
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE update_environ_boundary( bound )
!--------------------------------------------------------------------
    !
    USE tools_generate_boundary, ONLY : boundary_of_density, &
         & boundary_of_functions, boundary_of_system, &
         & solvent_aware_boundary, invert_boundary, &
         & compute_ion_field, compute_normal_field, compute_dion_field_drho, &
         & field_aware_density
    !
    IMPLICIT NONE
    !
    TYPE( environ_boundary ), INTENT(INOUT) :: bound
    !
    LOGICAL :: update_anything
    CHARACTER( LEN=80 ) :: sub_name = 'update_environ_boundary'
    !
    INTEGER :: i
    TYPE( environ_cell ), POINTER :: cell
    CHARACTER( LEN=80 ) :: label
    !
    cell => bound%scaled%cell
    !
    update_anything = .FALSE.
    IF ( bound % need_ions ) update_anything = bound % ions % update
    IF ( bound % need_electrons ) update_anything = update_anything .OR. bound % electrons % update
    IF ( bound % need_system ) update_anything = update_anything .OR. bound % system % update
    IF ( .NOT. update_anything ) THEN
       !
       ! ... Nothing is under update, change update_status and exit
       !
       IF ( bound % update_status .EQ. 2 ) bound % update_status = 0
       RETURN
       !
    ENDIF
    !
    SELECT CASE ( bound % mode )
       !
    CASE ( 'full' )
       !
       IF ( bound % ions % update ) THEN
          !
          ! ... Compute the ionic part
          !
          CALL density_of_functions( bound%ions%number, bound%ions%core_electrons, bound%ions%core, .TRUE. )
          !
          bound % update_status = 1 ! waiting to finish update
          !
       ENDIF
       !
       IF ( bound % electrons % update ) THEN
          !
          ! ... Check if the ionic part has been updated
          !
          IF ( bound % update_status .EQ. 0 ) &
               & CALL errore(sub_name,'Wrong update status, possibly missing ionic update',1)
          !
          bound % density % of_r = bound % electrons % density % of_r + bound % ions % core % of_r
          !
          CALL boundary_of_density( bound % density, bound )
          !
          bound % update_status = 2 ! boundary has changed and is ready
          !
       ENDIF
       !
    CASE ( 'electronic' )
       !
       IF ( bound % electrons % update ) THEN
          !
          bound % density % of_r = bound % electrons % density % of_r
          !
          CALL boundary_of_density( bound % density, bound )
          !
          bound % update_status = 2 ! boundary has changes and is ready
          !
          !CALL test_energy_derivatives( 1, bound )
          !
       ELSE
          !
          IF ( bound % update_status .EQ. 2 ) bound % update_status = 0 ! boundary has not changed
          RETURN
          !
       ENDIF
       !
    CASE ( 'ionic' )
       !
       IF ( bound % ions % update ) THEN
          !
          ! ... Only ions are needed, fully update the boundary
          !
          CALL boundary_of_functions( bound%ions%number, bound%soft_spheres, bound )
          !
          bound % update_status = 2 ! boundary has changed and is ready
          !
       ELSE
          !
          IF ( bound % update_status .EQ. 2 ) bound % update_status = 0 ! boundary has not changed
          RETURN
          !
       ENDIF
       !
    CASE ( 'fa-electronic' )
       !
       IF ( bound % ions % update ) THEN
          !
          bound % update_status = 1 ! waiting to finish update
          !
       ENDIF
       !
       IF ( bound % electrons % update ) THEN
          !
          CALL compute_normal_field( bound%ions, bound%electrons, bound%normal_field )
          !
          CALL field_aware_density( bound%electrons, bound )
          !
          CALL boundary_of_density( bound % density, bound )
          !
! ... TO DEBUG FIELD-AWARE: testing energy derivatives
          !
          !CALL extract_boundary_data( bound )
          !CALL test_energy_derivatives( 2, bound )
          !IF ( niter .EQ. 1 ) CALL test_energy_derivatives( 2, bound )
          !CALL test_normal_field_derivatives( bound )
          !
          bound % update_status = 2 ! boundary has changes and is ready
          !
       ENDIF
       !
    CASE ( 'fa-ionic' )
       !
       IF ( bound % ions % update ) THEN
          !
          CALL compute_dion_field_drho( bound%ions%number, bound%local_spheres, bound%dion_field_drho, bound%core%fft )
          !
          bound % update_status = 1 ! waiting to finish update
          !
       ENDIF
       !
       IF ( bound % electrons % update ) THEN
          !
          CALL compute_ion_field( bound%ions%number, bound%local_spheres, bound%ions, bound%electrons, bound%ion_field )
          !
          CALL set_soft_spheres( bound, .TRUE. )
          !
          CALL boundary_of_functions( bound%ions%number, bound%soft_spheres, bound )
          !
! ... TO DEBUG FIELD-AWARE: testing ion_field derivatives
!          !

          !CALL test_ion_field_derivatives( 2, bound )
          !CALL test_energy_derivatives( 2, bound )
!          !
! ... TO DEBUG FIELD-AWARE: testing energy derivatives
!          !
          !IF ( ionode ) WRITE(program_unit,'(1X,a,i14.7)')' niter = ', niter
          !IF ( ionode ) WRITE(environ_unit,'(a,i14.7)')' niter = ', niter
          !IF ( niter .EQ. 32 ) CALL test_energy_derivatives( 2, bound )
          !IF ( niter .EQ. 32 ) CALL test_ion_field_derivatives( 2, bound )
          !
          bound % update_status = 2 ! boundary has changes and is ready
          !
       ENDIF
       !
    CASE ( 'system' )
       !
       IF ( bound % system % update ) THEN
          !
          ! ... Only ions are needed, fully update the boundary
          !
          CALL boundary_of_system( bound % simple, bound )
          !
! ... TO DEBUG SOLVENT-AWARE
!          !
!          CALL invert_boundary( bound )
!          !
!          CALL test_de_dboundary( bound )
          !
          bound % update_status = 2 ! boundary has changed and is ready
          !
       ELSE
          !
          IF ( bound % update_status .EQ. 2 ) bound % update_status = 0 ! boundary has not changed
          RETURN
          !
       ENDIF
       !
    CASE DEFAULT
       !
       CALL errore(sub_name,'Unrecognized boundary mode',1)
       !
    END SELECT
    !
    ! Solvent-aware interface
    !
    IF ( bound % update_status .EQ. 2 .AND. bound % solvent_aware ) CALL solvent_aware_boundary( bound )
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE update_environ_boundary
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE destroy_environ_boundary(lflag, boundary)
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    LOGICAL, INTENT(IN) :: lflag
    TYPE( environ_boundary ), INTENT(INOUT) :: boundary
    CHARACTER (LEN=80) :: sub_name = 'destroy_environ_boundary'
    !
    INTEGER :: i
    !
    IF ( boundary%initialized ) THEN
       !
       CALL destroy_environ_density( boundary%scaled )
       IF ( boundary%mode .EQ. 'electronic' .OR. boundary%mode .EQ. 'full' ) THEN
          CALL destroy_environ_density( boundary%density )
          CALL destroy_environ_density( boundary%dscaled )
          CALL destroy_environ_density( boundary%d2scaled )
       ENDIF
       IF ( boundary%deriv .GE. 1 ) CALL destroy_environ_gradient( boundary%gradient )
       IF ( boundary%deriv .GE. 2 ) CALL destroy_environ_density( boundary%laplacian )
       IF ( boundary%deriv .GE. 3 ) CALL destroy_environ_density( boundary%dsurface )
       !
       IF ( boundary%solvent_aware ) THEN
          CALL destroy_environ_density( boundary%local )
          CALL destroy_environ_density( boundary%probe )
          CALL destroy_environ_density( boundary%filling )
          CALL destroy_environ_density( boundary%dfilling )
          IF ( boundary%deriv .GE. 3 ) CALL destroy_environ_hessian( boundary%hessian )
       ENDIF
       !
       IF ( boundary%field_aware ) THEN
          IF ( boundary%mode .EQ. 'fa-electronic' .OR. boundary%mode .EQ. 'fa-full' ) THEN
             CALL destroy_environ_density( boundary%normal_field )
          ELSE IF ( boundary%mode .EQ. 'fa-ionic' )THEN
             DO i = 1, boundary%ions%number
                CALL destroy_environ_density( boundary%dion_field_drho(i) )
             ENDDO
          ENDIF
       ENDIF
       !
       boundary%initialized = .FALSE.
       !
    END IF
    !
    IF ( lflag ) THEN
       !
       ! These components were allocated first, destroy only if lflag = .TRUE.
       !
       IF ( boundary % need_ions ) THEN
          IF ( boundary%mode .EQ. 'ionic' .OR. boundary%mode .EQ. 'fa-ionic' ) THEN
             CALL destroy_environ_functions( boundary%ions%number, boundary%soft_spheres )
             IF ( boundary%field_aware .AND. boundary%mode .EQ. 'fa-ionic' ) THEN
                DEALLOCATE( boundary%ion_field )
                DEALLOCATE( boundary%partial_of_ion_field )
                CALL destroy_environ_functions( boundary%ions%number, boundary%local_spheres )
                DEALLOCATE( boundary%dion_field_drho )
             ENDIF
          END IF
          IF (.NOT.ASSOCIATED(boundary%ions)) &
               & CALL errore(sub_name,'Trying to destroy a non associated object',1)
          NULLIFY(boundary%ions)
       ELSE
          IF (ASSOCIATED(boundary%ions))&
               & CALL errore(sub_name,'Found an unexpected associated object',1)
       ENDIF
       !
       IF ( boundary%need_electrons ) THEN
          IF (ASSOCIATED(boundary%electrons)) NULLIFY(boundary%electrons)
       ENDIF
       !
       IF ( boundary%solvent_aware ) DEALLOCATE(boundary%solvent_probe%pos)
       !
       IF ( boundary%need_system ) THEN
          IF (ASSOCIATED(boundary%system)) NULLIFY(boundary%system)
       ENDIF
       !
    ENDIF
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE destroy_environ_boundary
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE test_de_dboundary( boundary )
!--------------------------------------------------------------------
   !
    !
    USE tools_generate_boundary, ONLY : solvent_aware_boundary, solvent_aware_de_dboundary
    USE embedding_surface,       ONLY : calc_esurface, calc_desurface_dboundary
    USE embedding_volume,        ONLY : calc_evolume, calc_devolume_dboundary
    !
    IMPLICIT NONE
    !
    ! ... Test functional derivative of energy wrt local boundary
    !
    TYPE( environ_boundary ), INTENT(IN), TARGET :: boundary
    !
    TYPE( environ_boundary ) :: localbound
    TYPE( environ_density ) :: de_dboundary
    TYPE( environ_functions ) :: test_function
    TYPE( environ_density ) :: delta
    TYPE( environ_gradient ) :: graddelta
    !
    TYPE( environ_cell ), POINTER :: cell
    !
    INTEGER :: i
    REAL( DP ) :: localpressure, localsurface_tension
    REAL( DP ) :: eplus, eminus, de_fd, de_analytic, epsilon
    !
    cell => boundary % scaled % cell
    !
    CALL copy_environ_boundary( boundary, localbound )
    !
    CALL init_environ_density( cell, de_dboundary )
    !
    CALL solvent_aware_boundary( localbound )
    !
    localpressure = 100.D0
    CALL calc_devolume_dboundary( localpressure, localbound, de_dboundary )
    !
    localsurface_tension = 100.D0
    CALL calc_desurface_dboundary( localsurface_tension, localbound, de_dboundary )
    !
    CALL solvent_aware_de_dboundary( localbound, de_dboundary )
    !
    test_function % type = 1
    test_function % dim = 0
    test_function % axis = 3
    test_function % spread = 0.3D0
    test_function % width = 0.D0
    test_function % volume = 1.D0
    !
    epsilon = 0.000001
    !
    CALL init_environ_density( cell, delta )
    CALL init_environ_gradient( cell, graddelta )
    !
    ALLOCATE( test_function % pos( 3 ) )
    test_function % pos(1) = 11.79D0 / cell % alat
    test_function % pos(2) = 12.05D0 / cell % alat
    !
    DO i = 1, cell % dfft % nr3
       !
       test_function % pos(3) = DBLE(i-1) * cell % at(3,3) / DBLE( cell % dfft % nr3 )
       CALL density_of_functions( test_function, delta, .TRUE. )
       CALL gradient_of_functions( test_function, graddelta, .TRUE. )
       !
       de_fd = 0.D0
       CALL copy_environ_boundary( boundary, localbound )
       localbound % scaled % of_r = localbound % scaled % of_r + epsilon * delta % of_r
       localbound % volume = integrate_environ_density( localbound % scaled )
       IF ( localbound % deriv .GE. 1 ) THEN
          localbound % gradient % of_r = localbound % gradient % of_r + epsilon * graddelta % of_r
          CALL update_gradient_modulus( localbound % gradient )
          localbound % surface = integrate_environ_density( localbound % gradient % modulus )
       END IF
       !
       CALL solvent_aware_boundary( localbound )
       !
       CALL calc_evolume( localpressure, localbound, eplus )
       de_fd = de_fd + eplus
       !
       CALL calc_esurface( localsurface_tension, localbound, eplus )
       de_fd = de_fd + eplus
       !
       CALL copy_environ_boundary( boundary, localbound )
       localbound % scaled % of_r = localbound % scaled % of_r - epsilon * delta % of_r
       localbound % volume = integrate_environ_density( localbound % scaled )
       IF ( localbound % deriv .GE. 1 ) THEN
          localbound % gradient % of_r = localbound % gradient % of_r - epsilon * graddelta % of_r
          CALL update_gradient_modulus( localbound % gradient )
          localbound % surface = integrate_environ_density( localbound % gradient % modulus )
       END IF
       !
       CALL solvent_aware_boundary( localbound )
       !
       CALL calc_evolume( localpressure, localbound, eminus )
       de_fd = de_fd - eminus
       !
       CALL calc_esurface( localsurface_tension, localbound, eminus )
       de_fd = de_fd - eminus
       !
       de_fd = 0.5D0 * de_fd / epsilon
       !
       de_analytic = scalar_product_environ_density( de_dboundary, delta )
       !
       IF ( ionode ) WRITE(environ_unit,'(1X,a,f20.10,3f20.10)')' z = ',test_function % pos(3) * cell % alat,&
            & de_analytic, de_fd, de_analytic - de_fd
       FLUSH(environ_unit)
       !
    ENDDO
    !
    CALL destroy_environ_density( delta )
    CALL destroy_environ_gradient( graddelta )
    CALL destroy_environ_density( de_dboundary )
    !
    CALL destroy_environ_boundary( .TRUE., localbound )
    !
    STOP
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE test_de_dboundary
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE test_ion_field_derivatives( ideriv, bound )
!--------------------------------------------------------------------
    !
    USE utils_ions, ONLY : update_environ_ions
    USE tools_generate_boundary, ONLY : compute_ion_field, compute_ion_field_partial
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: ideriv ! .EQ. 1/2 for electronic/ionic
    TYPE( environ_boundary ), INTENT(INOUT) :: bound
    !
    TYPE( environ_cell ), POINTER :: cell
    TYPE( environ_density ) :: rho
    TYPE( environ_gradient ) :: field
    !
    INTEGER :: i, ipol
    REAL( DP ) :: dx, x0, epsilon
    REAL( DP ), ALLOCATABLE, DIMENSION(:) :: fd_partial_of_ion_field
    REAL( DP ), ALLOCATABLE, DIMENSION(:,:,:) :: analytic_partial_of_ion_field
    !
    INTEGER :: j
    REAL( DP ) :: tmp
    REAL( DP ), ALLOCATABLE, DIMENSION(:) :: fd_dion_field_drho
    TYPE( environ_density ), ALLOCATABLE, DIMENSION(:) :: analytic_dion_field_drho
    TYPE( environ_functions ) :: test_function
    TYPE( environ_density ) :: delta
    TYPE( environ_electrons ) :: localelectrons
    !
    cell => bound % scaled % cell
    !
    ! Recompute total charge density and field
    !
    CALL init_environ_density( cell, rho )
    rho % of_r = bound % electrons % density % of_r + bound % ions % density % of_r
    !
    CALL init_environ_gradient( cell, field )
    CALL gradv_h_of_rho_r( rho%of_r, field%of_r )
    !
    ! Print out individual and global fluxes
    !
    DO i = 1, bound % ions % number
       WRITE(environ_unit, '(a,i3,a,f14.7)' )&
            & 'flux through soft-sphere number ',i,' equal to ',bound%ion_field(i)
    ENDDO
    !
    rho % of_r = rho % of_r * bound % scaled % of_r
    WRITE( environ_unit, '(a,f14.7,a,f14.7)' )&
         & 'total_charge = ',integrate_environ_density(rho),' total flux throught soft-spheres = ',SUM(bound%ion_field(:))
    !
    CALL scalar_product_environ_gradient( field, bound%gradient, rho )
    WRITE( environ_unit, '(a,f14.7)' )'actual total flux = ',integrate_environ_density(rho)
    !
    IF ( ideriv .EQ. 1 ) THEN
       !
       ! Test functional derivative wrt electronic density of each flux
       !
       ALLOCATE( analytic_dion_field_drho( bound%ions%number ) )
       DO i = 1, bound % ions % number
          CALL copy_environ_density( bound%dion_field_drho(i), analytic_dion_field_drho(i) )
       ENDDO
       !
       ALLOCATE( fd_dion_field_drho( bound%ions%number ) )
       !
       test_function % type = 1
       test_function % dim = 0
       test_function % axis = 3
       test_function % spread = 0.3D0
       test_function % width = 0.D0
       test_function % volume = 1.D0
       !
       epsilon = 0.000001
       !
       CALL init_environ_density( cell, delta )
       !
       CALL init_environ_density( cell, localelectrons%density )
       !
       ! We are only going to check delta functions along the z axis passing throught the O atom
       !
       ALLOCATE( test_function % pos( 3 ) )
       test_function % pos(1) = 11.79D0 / cell % alat
       test_function % pos(2) = 12.05D0 / cell % alat
       !
       DO i = 1, cell % dfft % nr3
          !
          test_function % pos(3) = DBLE(i-1) * cell % at(3,3) / DBLE( cell % dfft % nr3 )
          CALL density_of_functions( test_function, delta, .TRUE. )
          !
          fd_dion_field_drho = 0.D0
          !
          localelectrons % density % of_r = bound % electrons % density % of_r - epsilon * delta%of_r
          !
          CALL compute_ion_field( bound%ions%number, bound%local_spheres, bound%ions, localelectrons, &
               & bound%ion_field )
          !
          fd_dion_field_drho = bound%ion_field
          !
          localelectrons % density % of_r = bound % electrons % density % of_r + epsilon * delta%of_r
          !
          CALL compute_ion_field( bound%ions%number, bound%local_spheres, bound%ions, localelectrons, &
               & bound%ion_field )
          !
          fd_dion_field_drho = bound%ion_field - fd_dion_field_drho
          fd_dion_field_drho = fd_dion_field_drho * 0.5D0 / epsilon
          !
          IF ( ionode ) WRITE( environ_unit, * )' z = ',test_function % pos(3) * cell % alat
          !
          DO j = 1, bound % ions % number
             tmp = scalar_product_environ_density( analytic_dion_field_drho(j), delta )
             IF ( ionode ) WRITE(environ_unit,'(a,i3,3f20.10)')'ion = ',j,&
                  & tmp, fd_dion_field_drho(j), tmp - fd_dion_field_drho(j)
          ENDDO
          !
       ENDDO
       !
       DEALLOCATE( test_function % pos )
       CALL destroy_environ_density( delta )
       !
       DO i = 1, bound % ions % number
          CALL destroy_environ_density( analytic_dion_field_drho(i) )
       ENDDO
       DEALLOCATE( analytic_dion_field_drho )
       DEALLOCATE( fd_dion_field_drho )
       !
    ELSE IF ( ideriv .EQ. 2 ) THEN
       !
       ! Test derivative wrt atomic positions with finite differences
       !
       CALL compute_ion_field_partial( bound%ions%number, bound%local_spheres, bound%ions, bound%electrons, &
            & bound%ion_field, bound%partial_of_ion_field, bound%core%fft )
       !
       ALLOCATE( analytic_partial_of_ion_field( 3, bound%ions%number, bound%ions%number ) )
       analytic_partial_of_ion_field = bound % partial_of_ion_field
       !
       ! dx expected units: BOHR
       dx = 2.0D-3
       ! convert to alat
       dx = dx / cell%alat
       !
       ALLOCATE( fd_partial_of_ion_field( bound%ions%number ) )
       !
       DO i = 1, bound % ions % number
          !
          DO ipol = 1, 3
             !
             x0 = bound % ions % tau( ipol, i )
             bound % ions % tau( ipol, i ) = x0 - dx
             !
             CALL update_environ_ions( bound%ions%number, bound%ions%tau, bound%ions )
             rho % of_r = bound%electrons%density%of_r + bound%ions%density%of_r
             CALL gradv_h_of_rho_r( rho%of_r, field%of_r )
             !
             CALL compute_ion_field( bound%ions%number, bound%local_spheres, bound%ions, bound%electrons, &
                  & bound%ion_field )
             !
             fd_partial_of_ion_field = bound%ion_field
             !
             bound % ions % tau( ipol, i ) = x0 + dx
             !
             CALL update_environ_ions( bound%ions%number, bound%ions%tau, bound%ions )
             rho % of_r = bound%electrons%density%of_r + bound%ions%density%of_r
             CALL gradv_h_of_rho_r( rho%of_r, field%of_r )
             !
             CALL compute_ion_field( bound%ions%number, bound%local_spheres, bound%ions, bound%electrons, &
                  & bound%ion_field )
             !
             fd_partial_of_ion_field = bound%ion_field - fd_partial_of_ion_field
             fd_partial_of_ion_field = fd_partial_of_ion_field / 2.D0 / dx / cell%alat
             !
             WRITE( environ_unit, * )' i  = ',i,' ipol = ',ipol
             WRITE( environ_unit, '(a,10f20.10)' )'analytic     = ',analytic_partial_of_ion_field( ipol, :, i )
             WRITE( environ_unit, '(a,10f20.10)' )'finite-diff  = ',fd_partial_of_ion_field(:)
             WRITE( environ_unit, * )' '
             !
             bound % ions % tau( ipol, i ) = x0
             !
          ENDDO
          !
       ENDDO
       !
       DEALLOCATE( analytic_partial_of_ion_field )
       DEALLOCATE( fd_partial_of_ion_field )
       !
    END IF
    !
    STOP
    !
!--------------------------------------------------------------------
  END SUBROUTINE test_ion_field_derivatives
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE update_test_boundary( bound, electrons )
!--------------------------------------------------------------------
    USE tools_generate_boundary
    !
    IMPLICIT NONE
    !
    TYPE( environ_boundary ), INTENT(INOUT) :: bound
    TYPE( environ_electrons ), INTENT(IN) :: electrons
    INTEGER :: debugcubes = 0
    CHARACTER(len=100) :: label = 'None'
    CHARACTER(len=80) :: sub_name = 'update_test_boundary'
    !
    SELECT CASE ( bound % mode )
      CASE ( 'fa-electronic' )
        !
        CALL compute_normal_field( bound%ions, electrons, bound%normal_field )
        CALL field_aware_density( electrons, bound )
        CALL boundary_of_density( bound % density, bound )
        !
        SELECT CASE ( debugcubes )
        CASE ( 2 )
          label = "standard"
          CALL write_cube( bound % scaled, label=label )
        CASE ( 1 )
          label = "standard"
          CALL write_cube( bound % density, label=label )
        CASE DEFAULT
        END SELECT
        !
      CASE ( 'fa-ionic' )
        !
        CALL compute_ion_field( bound%ions%number, bound%local_spheres, bound%ions, electrons, bound%ion_field )
        CALL set_soft_spheres( bound, .TRUE. )
        CALL boundary_of_functions( bound%ions%number, bound%soft_spheres, bound )
        !
      CASE ( 'electronic' )
        !
        CALL boundary_of_density( electrons % density, bound )
        !
      CASE DEFAULT
        !
        CALL errore(sub_name,'Unrecognized boundary mode',1)
        !
    END SELECT
!--------------------------------------------------------------------
  END SUBROUTINE update_test_boundary
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE print_ionic( ions )
!--------------------------------------------------------------------
    TYPE( environ_ions ), INTENT(IN) :: ions
    INTEGER :: i
    DO i = 1, ions%number
      PRINT *, "ATOM", i, ions%tau(1, i), ions%tau(2, i), ions%tau(3, i)
    ENDDO
!--------------------------------------------------------------------
  END SUBROUTINE print_ionic
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE test_energy_derivatives( ideriv, bound )
!--------------------------------------------------------------------
    !
    USE tools_generate_boundary
    USE utils_ions,              ONLY : update_environ_ions
    USE embedding_surface,       ONLY : calc_esurface, calc_desurface_dboundary
    USE embedding_volume,        ONLY : calc_evolume, calc_devolume_dboundary
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: ideriv
    TYPE( environ_boundary ), INTENT(INOUT) :: bound
    !
    TYPE( environ_cell ), POINTER :: cell
    !
    INTEGER :: i, ipol
    REAL( DP ) :: epsilon
    REAL( DP ) :: localpressure, localsurface_tension
    REAL( DP ) :: de_fd, de_analytic, etmp
    REAL( DP ) :: dx, x0, force(3), ssforce(3)
    REAL( DP ) :: flux
    REAL( DP ), ALLOCATABLE, DIMENSION(:, :) :: tau0
    TYPE( environ_density ) :: de_dboundary
    TYPE( environ_density ) :: vanalytic
    TYPE( environ_functions ) :: test_function
    TYPE( environ_density ) :: delta
    TYPE( environ_electrons ) :: localelectrons
    TYPE( environ_gradient ) :: partial
    !
    CHARACTER( LEN=80 ) :: sub_name = 'test_energy_derivatives'
    !
    cell => bound % scaled % cell
    !
    ! Compute the field and the field-aware interface
    !
    CALL update_test_boundary( bound, bound%electrons )
    !
    ! Compute functional derivative of the energy wrt interface
    !
    CALL init_environ_density( cell, de_dboundary )
    !
    localpressure = -0.35
    CALL calc_devolume_dboundary( localpressure, bound, de_dboundary )
    !
    localsurface_tension = 0.0
    !CALL calc_desurface_dboundary( localsurface_tension, bound, de_dboundary )
    !
    IF ( ideriv .EQ. 1 ) THEN
       !
       ! Compute functional derivative wrt electronic density
       !
       IF ( bound % mode .EQ. 'fa-ionic' ) THEN
          CALL compute_dion_field_drho( bound%ions%number, bound%local_spheres, bound%dion_field_drho, bound%core%fft )
       ENDIF
       !
       CALL init_environ_density( cell, vanalytic )
       !
       IF ( bound % field_aware ) THEN
          CALL field_aware_de_drho( bound, de_dboundary, vanalytic )
       ELSE 
          vanalytic % of_r = bound%dscaled%of_r * de_dboundary%of_r 
       ENDIF
       !
       ! Loop over gridpoints with rhoelec + or - a delta function
       !
       test_function % type = 1
       test_function % dim = 0
       test_function % axis = 3
       test_function % spread = 0.4D0
       test_function % width = 0.D0
       test_function % volume = 1.D0
       !
       epsilon = 0.000008
       !
       CALL init_environ_density( cell, delta )
       !
       CALL init_environ_density( cell, localelectrons%density )
       !
       ! We are only going to check delta functions along the z axis passing throught the O atom
       !
       ALLOCATE( test_function % pos( 3 ) )
       test_function % pos(1) = 0.0 / cell % alat
       test_function % pos(2) = 0.0 / cell % alat
       !
       DO i = 1, cell % dfft % nr3
          !
          de_fd = 0.D0
          !
          test_function % pos(3) = DBLE(i-1) * cell % at(3,3) / DBLE( cell % dfft % nr3 )
          IF ( ionode ) WRITE( program_unit,'(a,f14.7)')' z = ', test_function%pos(3)
          !IF (test_function % pos(3) .LE. 0.365) CYCLE
          !IF (test_function % pos(3) * cell % alat .LE. 4.5 .OR. test_function % pos(3) * cell % alat .GE. 7.0) CYCLE
          CALL density_of_functions( test_function, delta, .TRUE. )
          !
          localelectrons % density % of_r = bound % electrons % density % of_r + epsilon * delta%of_r
          !
          CALL update_test_boundary( bound, localelectrons )
          !
          CALL calc_evolume( localpressure, bound, etmp )
          de_fd = de_fd + etmp
          !
          CALL calc_esurface( localsurface_tension, bound, etmp )
          de_fd = de_fd + etmp
          !IF ( ionode ) WRITE(environ_unit,'(2f20.10)')bound%surface, bound%volume
          !
          localelectrons % density % of_r = bound % electrons % density % of_r - epsilon * delta%of_r
          !
          CALL update_test_boundary( bound, localelectrons )
          !
          CALL calc_evolume( localpressure, bound, etmp )
          de_fd = de_fd - etmp
          !
          CALL calc_esurface( localsurface_tension, bound, etmp )
          de_fd = de_fd - etmp
          !IF ( ionode ) WRITE(environ_unit,'(2f20.10)')bound%surface, bound%volume
          !
          de_fd = de_fd * 0.5D0 / epsilon
          !
          de_analytic = scalar_product_environ_density( vanalytic, delta )
          !
          !
          IF ( ionode ) WRITE(environ_unit,'(a,i14.7,3f20.10)')' z = ', i,&
               & de_analytic, de_fd, de_analytic - de_fd
          FLUSH(environ_unit)
          !
          !STOP
       ENDDO
       !
       DEALLOCATE( test_function % pos )
       CALL destroy_environ_density( delta )
       !
       CALL destroy_environ_density( vanalytic )
       !
    ELSE IF ( ideriv .EQ. 2 ) THEN
      PRINT *, 'ideriv = 2'
       !
       ! Compute partial derivative with respect to ionic positions
       !
       ! CALCULATE ion field partials in advance
       IF ( bound % mode .EQ. 'fa-ionic' ) THEN
          IF ( ionode ) WRITE( program_unit, '(1X,a)' ) 'outside compute_ion_field_partial'
          CALL compute_ion_field_partial( bound%ions%number, bound%local_spheres, bound%ions, bound%electrons, bound%ion_field, &
            & bound%partial_of_ion_field, bound%core%fft )
       ENDIF
       !
       ! dx expected units: BOHR
       dx = 2.0D-3
       ! convert to alat
       dx = dx / cell%alat
       !
       CALL init_environ_gradient( cell, partial )
       ALLOCATE( tau0 ( 3, bound%ions%number ) )
       !
       tau0(:,:) = bound % ions % tau(:,:)
       !
       DO i = 1, bound%ions%number
          !
          ! CALCULATE FORCE FIRST, reset the ionic positions
          bound % ions % tau( :, : ) = tau0( :, : )
          !
          CALL update_environ_ions( bound%ions%number, bound%ions%tau, bound%ions )
          CALL update_test_boundary( bound, bound%electrons )
          !
          IF ( bound % mode .EQ. 'fa-ionic' ) THEN
             CALL calc_dboundary_dions( i, bound, partial )
             ssforce = -scalar_product_environ_gradient_density( partial, de_dboundary )
          ENDIF
          !
          IF ( bound % field_aware ) CALL field_aware_dboundary_dions( i, bound, partial )
          force = - scalar_product_environ_gradient_density( partial, de_dboundary )
          !
          DO ipol = 1, 3
             !
             ! FINITE DIFFERENCE VALUE STORED HERE
             de_fd = 0.D0
             ! RESET IONS
             bound % ions % tau( :, : ) = tau0( :, : )
             !
             ! MINUS dx
             bound % ions % tau( ipol, i ) = tau0( ipol, i ) - dx
             !
             CALL update_environ_ions( bound%ions%number, bound%ions%tau, bound%ions )
             CALL update_test_boundary( bound, bound%electrons )
             !
             CALL calc_evolume( localpressure, bound, etmp )
             de_fd = de_fd - etmp
             CALL calc_esurface( localsurface_tension, bound, etmp )
             de_fd = de_fd - etmp
             !
             ! PLUS dx
             bound % ions % tau( ipol, i ) = tau0( ipol, i ) + dx
             !
             CALL update_environ_ions( bound%ions%number, bound%ions%tau, bound%ions )
             CALL update_test_boundary( bound, bound%electrons)
             !
             CALL calc_evolume( localpressure, bound, etmp )
             de_fd = de_fd + etmp
             CALL calc_esurface( localsurface_tension, bound, etmp )
             de_fd = de_fd + etmp
             !
             ! force is negative of the energy derivative
             de_fd = de_fd * (-0.5D0) / dx / cell%alat
             !
             IF ( bound % mode .EQ. 'fa-ionic') THEN
                IF ( ionode ) WRITE(environ_unit,'(a,i3,a,i3,4f20.10)')' i = ',i,' ipol = ',ipol,&
                     & force(ipol), ssforce(ipol), de_fd, force(ipol) - de_fd
             ELSE
                IF ( ionode ) WRITE(environ_unit,'(a,i3,a,i3,3f20.10)')' i = ',i,' ipol = ',ipol,&
                     & force(ipol), de_fd, force(ipol) - de_fd
             ENDIF
             !STOP
             !
          ENDDO
          !
       END DO
       !
       CALL destroy_environ_gradient( partial )
       !
    ENDIF
    !
    CALL destroy_environ_density( de_dboundary )
    STOP
    !
!--------------------------------------------------------------------
  END SUBROUTINE test_energy_derivatives
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE extract_boundary_data( bound )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    TYPE( environ_boundary ), INTENT(IN) :: bound
    !
    TYPE( environ_cell ), POINTER :: cell
    !
    INTEGER :: i, ipol
    REAL( DP ) :: epsilon
    REAL( DP ) :: localpressure, localsurface_tension
    REAL( DP ) :: de_fd, de_analytic, etmp
    REAL( DP ) :: dx, x0, force(3)
    TYPE( environ_density ) :: de_dboundary
    TYPE( environ_density ) :: vanalytic
    TYPE( environ_functions ) :: test_function
    TYPE( environ_density ) :: delta
    TYPE( environ_electrons ) :: localelectrons
    TYPE( environ_gradient ) :: partial
    !
    CHARACTER(len=100) :: strg = 'e'
    CHARACTER(len=100) :: strl = 'i'
    CHARACTER( LEN=80 ) :: sub_name = 'extract_boundary_data'
    !
    !
    cell => bound % scaled % cell
    !
    IF ( ionode ) WRITE( program_unit, '(1X,a)' ) 'extract_boundary_data'
    !
    ! Compute the field and the field-aware interface
    !
    CALL write_cube( bound % electrons % density, label=strg )
    CALL write_cube( bound % ions % density, label=strl )
    STOP
    !
!--------------------------------------------------------------------
  END SUBROUTINE extract_boundary_data
!--------------------------------------------------------------------
!----------------------------------------------------------------------------
END MODULE utils_boundary
!----------------------------------------------------------------------------
