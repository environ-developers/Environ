! Copyright (C) 2018 ENVIRON (www.quantum-environment.org)
!
!    This file is part of Environ version 1.0
!
!    Environ 1.0 is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 2 of the License, or
!    (at your option) any later version.
!
!    Environ 1.0 is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more detail, either the file
!    `License' in the root directory of the present distribution, or
!    online at <http://www.gnu.org/licenses/>.
!
! Module containing the main routines to handle
!
!              environ_dielectric
!
! derived data types.
!
! Environ_dielectric is the type to store the details of the dielectric
! embedding. It contains the specifics of externally-defined dielectric
! regions and it links the boundary details. Starting from these quantities,
! It builds the dielectric function in space (stored in %epsilon component)
! and the factors derived from it that are required by the generalized
! Poisson solver (gradient of the logarithm, sqrt factor needed by
! preconditioned conjugate gradient, etc.).
!
!----------------------------------------------------------------------------
!  TYPE environ_dielectric
!----------------------------------------------------------------------------
!
!     ! Update status
!
!     LOGICAL :: update = .FALSE.
!
!     LOGICAL :: initialized = .FALSE.
!
!     ! Basic properties of the dielectric space from input
!
!     INTEGER :: nregions
!     TYPE( environ_functions ), DIMENSION(:), ALLOCATABLE :: regions
!
!     REAL( DP ) :: constant
!     TYPE( environ_density ) :: background
!     TYPE( environ_gradient ) :: gradbackground
!     TYPE( environ_density ) :: laplbackground
!
!     ! Boundary is the pointer to the object controlling
!     ! the interface between the QM and the continuum region
!
!     TYPE( environ_boundary ), POINTER :: boundary
!
!     ! The dielectric function over space is built from the
!     ! boundary of the continuum environment and the basic dielectric
!     ! properties of space
!
!     TYPE( environ_density ) :: epsilon
!     TYPE( environ_density ) :: depsilon ! This is needed in the extra term of kohn-sham/forces
!
!     ! Quantities related to the dielectric permittivity and
!     ! thay may be needed by the different solvers
!
!     LOGICAL :: need_gradient = .FALSE.
!     TYPE( environ_gradient ) :: gradient
!
!     LOGICAL :: need_factsqrt = .FALSE.
!     TYPE( environ_density ) :: factsqrt
!
!     LOGICAL :: need_gradlog = .FALSE.
!     TYPE( environ_gradient ) :: gradlog
!
!     ! Dielectric polarization charges and individual components
!
!     TYPE( environ_density ) :: density
!
!     LOGICAL :: need_auxiliary = .FALSE.
!     TYPE( environ_density ) :: iterative
!
!     REAL( DP ) :: charge = 0.0_DP
!
!----------------------------------------------------------------------------
!  END TYPE environ_dielectric
!----------------------------------------------------------------------------
!
! Authors: Oliviero Andreussi (Department of Physics, UNT)
!          Francesco Nattino  (THEOS and NCCR-MARVEL, EPFL)
!          Ismaila Dabo       (DMSE, Penn State)
!          Nicola Marzari     (THEOS and NCCR-MARVEL, EPFL)
!
!----------------------------------------------------------------------------
MODULE dielectric
  !--------------------------------------------------------------------------
  !
  USE environ_types
  USE environ_output
  USE functions
  USE environ_base, ONLY : e2, add_jellium
  !
  IMPLICIT NONE
  !
  SAVE
  !
  PRIVATE
  !
  PUBLIC :: create_environ_dielectric, init_environ_dielectric_first, &
       & set_dielectric_regions, init_environ_dielectric_second, &
       & update_environ_dielectric, calc_dedielectric_dboundary, &
       & destroy_environ_dielectric, dielectric_of_potential
  !
CONTAINS
  !
!--------------------------------------------------------------------
  SUBROUTINE create_environ_dielectric(dielectric)
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    TYPE( environ_dielectric ), INTENT(INOUT) :: dielectric
    !
    CHARACTER( LEN=80 ) :: sub_name = 'create_environ_dielectric'
    CHARACTER( LEN=80 ) :: label = ' '
    !
    dielectric%constant = 1.0_DP
    !
    IF ( ALLOCATED( dielectric%regions ) ) &
         & CALL errore(sub_name,'Trying to create an already allocated object',1)
    !
    label = 'background'
    CALL create_environ_density( dielectric%background, label  )
    label = 'gradbackground'
    CALL create_environ_gradient( dielectric%gradbackground, label  )
    label = 'laplbackground'
    CALL create_environ_density( dielectric%laplbackground, label  )
    !
    label = 'epsilon'
    CALL create_environ_density( dielectric%epsilon, label )
    label = 'depsilon'
    CALL create_environ_density( dielectric%depsilon, label )
    !
    NULLIFY( dielectric%boundary )
    !
    label = 'epsilon_gradlog'
    CALL create_environ_gradient( dielectric%gradlog, label )
    dielectric%need_gradient = .FALSE.
    label = 'epsilon_gradient'
    CALL create_environ_gradient( dielectric%gradient, label )
    dielectric%need_factsqrt = .FALSE.
    label = 'epsilon_factsqrt'
    CALL create_environ_density( dielectric%factsqrt, label )
    !
    label = 'density'
    CALL create_environ_density( dielectric%density, label )
    dielectric%need_auxiliary = .FALSE.
    label = 'iterative'
    CALL create_environ_density( dielectric%iterative, label )
    !
    dielectric%charge = 0.D0
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE create_environ_dielectric
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE init_environ_dielectric_first( constant, boundary, &
             & need_gradient, need_factsqrt, need_auxiliary, dielectric )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    REAL( DP ) :: constant
    LOGICAL, INTENT(IN) :: need_gradient, need_factsqrt, need_auxiliary
    TYPE( environ_boundary ), TARGET, INTENT(IN) :: boundary
    TYPE( environ_dielectric ), INTENT(INOUT) :: dielectric
    !
    dielectric%constant = constant
    !
    dielectric%boundary => boundary
    !
    dielectric%need_gradient = need_gradient
    dielectric%need_factsqrt = need_factsqrt
    !
    dielectric%need_auxiliary = need_auxiliary
    !
    dielectric%initialized = .FALSE.
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE init_environ_dielectric_first
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE set_dielectric_regions( nregions, epsregion_dim, epsregion_axis, &
       & epsregion_pos, epsregion_width, epsregion_spread, epsregion_eps, dielectric )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: nregions
    INTEGER, DIMENSION(nregions), INTENT(IN) :: epsregion_dim, epsregion_axis
    REAL( DP ), DIMENSION(nregions), INTENT(IN) :: epsregion_width, epsregion_spread, epsregion_eps
    REAL( DP ), DIMENSION(3,nregions), INTENT(IN) :: epsregion_pos
    TYPE( environ_dielectric ), INTENT(INOUT) :: dielectric
    !
    dielectric%nregions = nregions
    IF ( dielectric%nregions .GT. 0 ) &
         & CALL create_environ_functions( dielectric%nregions, 4, epsregion_dim, &
         & epsregion_axis, epsregion_pos, epsregion_width, epsregion_spread, &
         & epsregion_eps, dielectric%regions )
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE set_dielectric_regions
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE init_environ_dielectric_second( cell, dielectric )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    TYPE( environ_cell ), INTENT(INOUT) :: cell
    TYPE( environ_dielectric ), INTENT(INOUT) :: dielectric
    !
    CHARACTER( LEN=80 ) :: sub_name = 'init_environ_dielectric_second'
    !
    INTEGER :: i
    !
    IF ( dielectric % nregions .GT. 0 ) THEN
       DO i = 1, dielectric % nregions
          dielectric % regions(i) % pos = dielectric % regions(i) % pos / cell % alat
       END DO
    END IF
    !
    CALL init_environ_density( cell, dielectric%background )
    dielectric % background % of_r(:) = dielectric % constant
    IF ( dielectric % nregions .GT. 0 ) THEN
       CALL init_environ_gradient( cell, dielectric%gradbackground )
       IF ( dielectric%need_factsqrt ) CALL init_environ_density( cell, dielectric%laplbackground )
    END IF
    !
    CALL init_environ_density( cell, dielectric%epsilon )
    CALL init_environ_density( cell, dielectric%depsilon )
    !
    CALL init_environ_gradient( cell, dielectric%gradlog )
    IF ( dielectric%need_gradient ) CALL init_environ_gradient( cell, dielectric%gradient )
    IF ( dielectric%need_factsqrt ) CALL init_environ_density( cell, dielectric%factsqrt )
    !
    CALL init_environ_density( cell, dielectric%density )
    IF ( dielectric%need_auxiliary ) CALL init_environ_density( cell, dielectric%iterative )
    !
    dielectric%initialized = .TRUE.
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE init_environ_dielectric_second
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE update_environ_dielectric( dielectric )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    TYPE( environ_dielectric ), INTENT(INOUT) :: dielectric
    !
    CALL start_clock( 'dielectric' )
    !
    IF ( dielectric % epsilon % cell % update ) THEN
       !
       ! ... Cells has changed, may need to update the background
       !
       IF ( dielectric % nregions .GT. 0 ) THEN
          !
          ! ... Recompute background dielectric and its derivative
          !
          CALL update_dielectric_background( dielectric )
          !
          ! ... Background has changed, need to update the dielectric when ready
          !
          dielectric % update = .TRUE.
          !
       ENDIF
       !
    ENDIF
    !
    ! ... Check if the boundary is under update (status = 1) or has been fully updated (status = 2)
    !
    IF ( dielectric % boundary % update_status .GT. 0 ) dielectric % update = .TRUE.
    !
    IF ( dielectric % update ) THEN
       !
       ! ... Update the dielectric in space and its derivatives if the boundary is ready
       !
       IF ( dielectric % boundary % update_status .EQ. 2 ) THEN
          !
          CALL dielectric_of_boundary( dielectric )
          !
          dielectric % update = .FALSE.
          !
       ENDIF
       !
    END IF
    !
    CALL stop_clock( 'dielectric' )
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE update_environ_dielectric
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE update_dielectric_background( dielectric )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    TYPE( environ_dielectric ), INTENT(INOUT) :: dielectric
    !
    INTEGER :: i, ipol
    TYPE( environ_density ) :: local
    TYPE( environ_gradient ) :: gradlocal
    TYPE( environ_density ) :: lapllocal
    TYPE( environ_cell ), POINTER :: cell
    !
    dielectric % background % of_r = dielectric % constant
    !
    IF ( dielectric % nregions .LE. 0 ) RETURN
    !
    dielectric % gradbackground % of_r = 0.D0
    IF ( dielectric % need_factsqrt ) dielectric % laplbackground % of_r = 0.D0
    !
    cell => dielectric % background % cell
    !
    CALL init_environ_density( cell , local )
    CALL init_environ_gradient( cell, gradlocal )
    IF ( dielectric % need_factsqrt ) CALL init_environ_density( cell, lapllocal )
    !
    DO i = 1, dielectric % nregions
       !
       CALL density_of_functions( dielectric%regions(i), local, .TRUE. )
       CALL gradient_of_functions( dielectric%regions(i), gradlocal, .TRUE. )
       !
       ! Update background and derivatives in reverse order
       !
       IF ( dielectric % need_factsqrt ) THEN
          CALL scalar_product_environ_gradient( dielectric % gradbackground, gradlocal, lapllocal )
          dielectric % laplbackground % of_r = dielectric % laplbackground % of_r(:) * &
               & ( 1.D0 - local % of_r(:) / dielectric % regions(i) % volume ) - &
               & 2.D0 * lapllocal % of_r / dielectric % regions(i) % volume
          CALL laplacian_of_functions( dielectric%regions(i), lapllocal, .TRUE. )
          dielectric % laplbackground % of_r(:) = dielectric % laplbackground % of_r(:) + &
               & lapllocal % of_r(:) * ( 1.D0 - dielectric % background % of_r(:) / dielectric % regions(i) % volume )
       ENDIF
       DO ipol = 1, 3
          dielectric % gradbackground % of_r(ipol,:) = &
               & dielectric % gradbackground % of_r(ipol,:) * ( 1.D0 - local % of_r(:) / dielectric % regions(i) % volume ) + &
               & gradlocal % of_r (ipol,:) * ( 1.D0 - dielectric % background % of_r(:) / dielectric % regions(i) % volume )
       END DO
       dielectric % background % of_r(:) = dielectric % background % of_r(:) + &
            & local % of_r(:) * ( 1.D0 - dielectric % background % of_r(:) / dielectric % regions(i) % volume )
       !
    ENDDO
    !
    CALL destroy_environ_density( local )
    CALL destroy_environ_gradient( gradlocal )
    IF ( dielectric % need_factsqrt ) CALL destroy_environ_density( lapllocal )
    !
    IF ( verbose .GE. 3 ) CALL print_environ_density( dielectric%background )
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE update_dielectric_background
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE dielectric_of_boundary( dielectric )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    TYPE(environ_dielectric), TARGET, INTENT(INOUT) :: dielectric
    !
    INTEGER :: ipol
    INTEGER, POINTER :: nnr
    REAL( DP ), DIMENSION( : ), POINTER :: factsqrteps, eps, deps, const
    REAL( DP ), DIMENSION( : ), POINTER :: scaled, gradscaledmod, laplscaled
    REAL( DP ), DIMENSION( :, : ), POINTER :: gradeps, gradlogeps, gradscaled
    !
    REAL( DP ), DIMENSION( : ), POINTER :: laplback, gradbackmod
    REAL( DP ), DIMENSION( :, : ), POINTER :: gradback
    !
    REAL( DP ), DIMENSION( : ), ALLOCATABLE :: dlogeps
    REAL( DP ), DIMENSION( : ), ALLOCATABLE :: d2eps
    !
    REAL( DP ), DIMENSION( : ), ALLOCATABLE :: deps_dback
    REAL( DP ), DIMENSION( : ), ALLOCATABLE :: dlogeps_dback
    REAL( DP ), DIMENSION( : ), ALLOCATABLE :: d2eps_dback2
    REAL( DP ), DIMENSION( : ), ALLOCATABLE :: d2eps_dbackdbound
    REAL( DP ), DIMENSION( : ), ALLOCATABLE :: gradepsmod2
    !
    CHARACTER ( LEN=80 ) :: sub_name = 'epsilon_of_boundary'
    !
    ! Aliases and sanity checks
    !
    nnr => dielectric%epsilon%cell%nnr
    eps => dielectric % epsilon % of_r
    deps => dielectric % depsilon % of_r
    const => dielectric % background % of_r
    scaled => dielectric % boundary % scaled % of_r
    !
    IF ( .NOT. ALLOCATED( dielectric % boundary % gradient % of_r ) ) &
         & CALL errore(sub_name,'Missing required gradient of boundary',1)
    gradscaled => dielectric % boundary % gradient % of_r
    gradlogeps => dielectric % gradlog % of_r
    ALLOCATE( dlogeps( nnr ) )
    !
    IF ( dielectric % need_gradient ) THEN
       IF ( .NOT. ALLOCATED( dielectric % boundary % gradient % of_r ) ) &
            & CALL errore(sub_name,'Missing required gradient of boundary',1)
       gradscaled => dielectric % boundary % gradient % of_r
       gradeps => dielectric % gradient % of_r
    ENDIF
    IF ( dielectric % need_factsqrt ) THEN
       IF ( .NOT. ALLOCATED( dielectric % boundary % gradient % of_r ) ) &
            & CALL errore(sub_name,'Missing required gradient of boundary',1)
       gradscaledmod => dielectric % boundary % gradient % modulus % of_r
       IF ( .NOT. ALLOCATED( dielectric % boundary % laplacian % of_r ) ) &
            & CALL errore(sub_name,'Missing required laplacian of boundary',1)
       laplscaled => dielectric % boundary % laplacian % of_r
       factsqrteps => dielectric % factsqrt % of_r
       ALLOCATE( d2eps( nnr ) )
    ENDIF
    !
    IF ( dielectric % nregions .GT. 0 ) THEN
       gradback => dielectric % gradbackground % of_r
       ALLOCATE( deps_dback( nnr ) )
       ALLOCATE( dlogeps_dback( nnr ) )
       IF ( dielectric % need_factsqrt ) THEN
          laplback => dielectric % laplbackground % of_r
          gradbackmod => dielectric % gradbackground % modulus % of_r
          ALLOCATE( gradepsmod2( nnr ) )
          ALLOCATE( d2eps_dbackdbound( nnr ) )
          ALLOCATE( d2eps_dback2( nnr ) )
       END IF
    END IF
    !
    ! Compute epsilon(r) and its derivative wrt boundary
    !
    SELECT CASE ( dielectric%boundary%type )
       !
    CASE ( 0, 2 )
       !
       eps = 1.D0 + ( const - 1.D0 ) * ( 1.D0 - scaled )
       deps = ( 1.D0 - const )
       dlogeps = deps / eps
       IF ( dielectric % need_factsqrt ) d2eps = 0.D0
       IF ( dielectric % nregions .GT. 0 ) THEN
          deps_dback = 1.D0 - scaled
          dlogeps_dback = deps_dback / eps
          IF ( dielectric % need_factsqrt ) THEN
             d2eps_dback2 = 0.D0
             d2eps_dbackdbound = -1.D0
          ENDIF
       ENDIF
       !
    CASE ( 1 )
       !
       eps = EXP( LOG( const ) * ( 1.D0 - scaled ) )
       deps = - eps * LOG( const )
       dlogeps = - LOG( const )
       IF ( dielectric % need_factsqrt ) d2eps = eps * LOG( const )**2
       IF ( dielectric % nregions .GT. 0 ) THEN
          deps_dback = eps * ( 1.D0 - scaled ) / const
          dlogeps_dback = ( 1.D0 - scaled ) / const
          IF ( dielectric % need_factsqrt ) THEN
             d2eps_dback2 = - deps_dback * scaled / const
             d2eps_dbackdbound = eps / const * ( 1.D0 - ( 1.D0 - scaled ) * LOG( const ) )
          ENDIF
       ENDIF
       !
    CASE DEFAULT
       !
       CALL errore(sub_name,'Unkown boundary type',1)
       !
    END SELECT
    !
    ! If needed, compute derived quantites
    !
    DO ipol = 1, 3
       gradlogeps( ipol, : ) = dlogeps( : ) * gradscaled( ipol, : )
       IF ( dielectric % nregions .GT. 0 ) gradlogeps( ipol, : ) = &
            & gradlogeps( ipol, : ) + dlogeps_dback( : ) * gradback( ipol, : )
    ENDDO
    CALL update_gradient_modulus( dielectric%gradlog )
    DEALLOCATE( dlogeps )
    !
    IF ( dielectric % need_gradient ) THEN
       DO ipol = 1, 3
          gradeps( ipol, : ) = deps( : ) * gradscaled( ipol, : )
          IF ( dielectric % nregions .GT. 0 ) gradeps( ipol, : ) = &
               & gradeps( ipol, : ) + deps_dback( : ) * gradback( ipol, : )
       ENDDO
       CALL update_gradient_modulus( dielectric%gradient )
    END IF
    IF ( dielectric % need_factsqrt ) THEN
       IF ( dielectric % nregions .LE. 0 ) THEN
          factsqrteps = ( d2eps - 0.5D0 * deps**2 / eps ) * gradscaledmod**2 + deps * laplscaled
       ELSE
          CALL scalar_product_environ_gradient( dielectric % boundary % gradient, &
               & dielectric % gradbackground, dielectric % factsqrt )
          IF ( dielectric % need_gradient ) THEN
             gradepsmod2 = dielectric % gradient % modulus % of_r**2
          ELSE
             gradepsmod2 = 0.D0
             DO ipol = 1, 3
                gradepsmod2(:) = gradepsmod2(:) + ( deps( : ) * gradscaled( ipol, : ) + &
                     & deps_dback( : ) * gradback( ipol, : ) )**2
             ENDDO
          ENDIF
          factsqrteps = 2.D0 * d2eps_dbackdbound * factsqrteps + &
               & d2eps_dback2 * gradbackmod**2 + deps_dback * laplback + &
               & d2eps * gradscaledmod**2 + deps * laplscaled - &
               & 0.5D0 * gradepsmod2 / eps
       ENDIF
       factsqrteps = factsqrteps * 0.5D0 / e2 / fpi
       DEALLOCATE( d2eps )
    ENDIF
    !
    IF ( dielectric % nregions .GT. 0 ) THEN
       DEALLOCATE( deps_dback )
       DEALLOCATE( dlogeps_dback )
       IF ( dielectric % need_factsqrt ) THEN
          DEALLOCATE( d2eps_dback2 )
          DEALLOCATE( d2eps_dbackdbound )
       ENDIF
    ENDIF
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE dielectric_of_boundary
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE dielectric_of_potential( charges, potential, dielectric )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    TYPE( environ_density ), INTENT(IN) :: charges
    TYPE( environ_density ), INTENT(IN) :: potential
    TYPE( environ_dielectric ), INTENT(INOUT) :: dielectric
    !
    TYPE( environ_cell ), POINTER :: cell
    !
    REAL( DP ) :: jellium
    TYPE( environ_gradient ) :: gradient
    CHARACTER( LEN=80 ) :: sub_name = 'dielectric_of_potential'
    !
    IF ( .NOT. ASSOCIATED( potential % cell, charges % cell ) ) &
         & CALL errore(sub_name,'Missmatch in domains of potential and charges',1)
    IF ( .NOT. ASSOCIATED( potential % cell, dielectric % density % cell ) ) &
         & CALL errore(sub_name,'Missmatch in domains of potential and dielectric',1)
    cell => charges % cell
    !
    CALL init_environ_gradient( cell, gradient )
    CALL external_gradient( potential%of_r, gradient%of_r )
    CALL scalar_product_environ_gradient( dielectric%gradlog, gradient, dielectric % density )
    !
    jellium = 0.D0
    IF ( add_jellium ) jellium = integrate_environ_density( charges ) / cell % omega
    dielectric % density % of_r = dielectric % density % of_r / fpi / e2 + ( charges % of_r - jellium ) * &
         & ( 1.D0 - dielectric % epsilon % of_r ) / dielectric % epsilon % of_r
    !
    CALL destroy_environ_gradient( gradient )
    !
    dielectric % charge = integrate_environ_density( dielectric % density )
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE dielectric_of_potential
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE calc_dedielectric_dboundary( dielectric, velectrostatic, de_dboundary )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    TYPE( environ_dielectric ), INTENT(IN) :: dielectric
    TYPE( environ_density ), INTENT(IN) :: velectrostatic
    TYPE( environ_density ), INTENT(INOUT) :: de_dboundary
    !
    TYPE( environ_cell ), POINTER :: cell
    TYPE( environ_gradient ) :: gradient
    !
    cell => de_dboundary%cell
    !
    CALL init_environ_gradient( cell, gradient )
    CALL external_gradient( velectrostatic%of_r, gradient%of_r )
    CALL update_gradient_modulus( gradient )
    !
    de_dboundary % of_r = de_dboundary % of_r - &
         & gradient % modulus % of_r**2 * dielectric % depsilon % of_r * &
         & 0.5D0 / fpi / e2
    !
    CALL destroy_environ_gradient( gradient )
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE calc_dedielectric_dboundary
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE destroy_environ_dielectric(lflag,dielectric)
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    LOGICAL, INTENT(IN) :: lflag
    TYPE( environ_dielectric ), INTENT(INOUT) :: dielectric
    CHARACTER (LEN=80) :: sub_name = 'destroy_environ_dielectric'
    !
    IF ( lflag ) THEN
       !
       ! These components were allocated first, destroy only if lflag = .TRUE.
       !
       IF ( dielectric%nregions .GT. 0 ) THEN
          CALL destroy_environ_functions( dielectric%nregions, dielectric%regions )
       ELSE
          IF ( ALLOCATED(dielectric%regions) ) &
               & CALL errore(sub_name,'Found unexpected allocated object',1)
       END IF
       !
       IF (.NOT.ASSOCIATED(dielectric%boundary)) &
            & CALL errore(sub_name,'Trying to destroy a non associated object',1)
       NULLIFY( dielectric%boundary )
       !
    END IF
    !
    IF ( dielectric%initialized ) THEN
       !
       CALL destroy_environ_density( dielectric%background )
       IF ( dielectric % nregions .GT. 0 ) THEN
          CALL destroy_environ_gradient( dielectric%gradbackground )
          IF ( dielectric%need_factsqrt ) CALL destroy_environ_density( dielectric%laplbackground )
       END IF
       !
       CALL destroy_environ_density( dielectric%epsilon )
       CALL destroy_environ_density( dielectric%depsilon )
       !
       CALL destroy_environ_gradient( dielectric%gradlog )
       IF ( dielectric%need_gradient ) CALL destroy_environ_gradient( dielectric%gradient )
       IF ( dielectric%need_factsqrt ) CALL destroy_environ_density( dielectric%factsqrt )
       !
       CALL destroy_environ_density( dielectric%density )
       IF ( dielectric%need_auxiliary ) CALL destroy_environ_density( dielectric%iterative )
       !
       dielectric%initialized = .FALSE.
       !
    END IF
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE destroy_environ_dielectric
!--------------------------------------------------------------------
!----------------------------------------------------------------------------
END MODULE dielectric
!----------------------------------------------------------------------------
