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
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE create_environ_boundary
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE init_environ_boundary_first( need_gradient, need_laplacian, &
       & need_hessian, mode, stype, rhomax, rhomin, tbeta, const, alpha, &
       & softness, system_distance, system_spread, solvent_radius, radial_scale, &
       & radial_spread, filling_threshold, filling_spread, electrons, ions, system, &
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
    boundary%need_electrons = ( mode .EQ. 'electronic' ) .OR. ( mode .EQ. 'full' )
    IF ( boundary%need_electrons ) boundary%electrons => electrons
    boundary%need_ions = ( mode .EQ. 'ionic' ) .OR. ( mode .EQ. 'full' )
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
    IF ( boundary%need_ions .AND. .NOT. boundary%need_electrons ) &
         & ALLOCATE( boundary%soft_spheres( boundary%ions%number ) )
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
    CALL init_environ_density( cell, boundary%scaled )
    !
    IF ( boundary%mode .NE. 'ionic' ) THEN
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
    INTEGER :: i, n
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
    bcopy % initialized       = boriginal % initialized
    !
    CALL copy_environ_density   ( boriginal % scaled        , bcopy % scaled        )
    CALL copy_environ_gradient  ( boriginal % gradient      , bcopy % gradient      )
    CALL copy_environ_density   ( boriginal % laplacian     , bcopy % laplacian     )
    CALL copy_environ_density   ( boriginal % dsurface      , bcopy % dsurface      )
    CALL copy_environ_hessian   ( boriginal % hessian       , bcopy % hessian       )
    CALL copy_environ_density   ( boriginal % density       , bcopy % density       )
    CALL copy_environ_density   ( boriginal % dscaled       , bcopy % dscaled       )
    CALL copy_environ_density   ( boriginal % d2scaled      , bcopy % d2scaled      )
    CALL copy_environ_functions ( boriginal % simple        , bcopy % simple        )
    CALL copy_environ_functions ( boriginal % solvent_probe , bcopy % solvent_probe )
    CALL copy_environ_density   ( boriginal % local         , bcopy % local         )
    CALL copy_environ_density   ( boriginal % probe         , bcopy % probe         )
    CALL copy_environ_density   ( boriginal % filling       , bcopy % filling       )
    CALL copy_environ_density   ( boriginal % dfilling      , bcopy % dfilling      )
    !
    IF ( ALLOCATED( boriginal % soft_spheres ) ) THEN
       n = SIZE( boriginal % soft_spheres )
       IF ( ALLOCATED( bcopy % soft_spheres ) ) DEALLOCATE( bcopy % soft_spheres )
       ALLOCATE( bcopy % soft_spheres( n ) )
       DO i = 1, n
          CALL copy_environ_functions ( boriginal % soft_spheres(i), bcopy % soft_spheres(i) )
       ENDDO
    ELSE
       IF ( ALLOCATED( bcopy % soft_spheres ) ) DEALLOCATE( bcopy%soft_spheres )
    ENDIF
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE copy_environ_boundary
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE set_soft_spheres( boundary )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    TYPE( environ_boundary ), INTENT(INOUT) :: boundary
    !
    INTEGER :: i
    REAL( DP ) :: radius
    !
    IF ( boundary % mode .NE. 'ionic' ) RETURN
    !
    DO i = 1, boundary%ions%number
       radius = boundary%ions%iontype(boundary%ions%ityp(i))%solvationrad * boundary%alpha
       boundary%soft_spheres(i) = environ_functions(5,1,0,radius,boundary%softness,1.D0,&
            & boundary%ions%tau(:,i))
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
         & solvent_aware_boundary, invert_boundary
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
    TYPE( environ_density ) :: local
    !
    cell => bound%density%cell
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
    CASE ( 'system' )
       !
       IF ( bound % system % update ) THEN
          !
          ! ... Only ions are needed, fully update the boundary
          !
          CALL boundary_of_system( bound % simple, bound )
          !
! ... TO DEBUG SOLVENT-AWARE
!
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
    IF ( lflag ) THEN
       !
       ! These components were allocated first, destroy only if lflag = .TRUE.
       !
       IF ( boundary%need_ions ) THEN
          IF ( .NOT. boundary%need_electrons ) &
               & CALL destroy_environ_functions( boundary%ions%number, boundary%soft_spheres )
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
    IF ( boundary%initialized ) THEN
       !
       CALL destroy_environ_density( boundary%scaled )
       IF ( boundary%mode .NE. 'ionic' ) THEN
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
       boundary%initialized = .FALSE.
       !
    END IF
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
    INTEGER, POINTER :: nnr
    TYPE( environ_cell ), POINTER :: cell
    !
    INTEGER :: i
    REAL( DP ) :: localpressure, localsurface_tension
    REAL( DP ) :: eplus, eminus, de_fd, de_analytic, epsilon
    !
    cell => boundary % scaled % cell
    nnr => boundary % scaled % cell % nnr
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
!----------------------------------------------------------------------------
END MODULE utils_boundary
!----------------------------------------------------------------------------
