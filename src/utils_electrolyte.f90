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
!              environ_electrolyte
!
! derived data types.
!
! Environ_electrolyte contains all the specifications and the details of
! the electrolyte medium and of the ionic distribution in it.
!
!----------------------------------------------------------------------------
!  TYPE environ_electrolyte
!----------------------------------------------------------------------------
!     !
!     ! Update status
!     !
!     LOGICAL :: update = .FALSE.
!     !
!     LOGICAL :: initialized = .FALSE.
!     !
!     CHARACTER( LEN=80 ) :: electrolyte_entropy
!     LOGICAL :: linearized = .FALSE.
!     INTEGER :: ntyp
!     TYPE( environ_ioncctype ), DIMENSION(:), ALLOCATABLE :: ioncctype
!     !
!     REAL( DP ) :: temperature
!     REAL( DP ) :: k2
!     REAL( DP ) :: cionmax
!     REAL( DP ) :: permittivity
!     !
!     TYPE( environ_boundary ) :: boundary
!     !
!     TYPE( environ_density ) :: density
!     !
!     ! The electrolyte switch function and relate quantities
!     !
!     TYPE( environ_density ) :: gamma
!     TYPE( environ_density ) :: dgamma
!     !
!     TYPE( environ_density ) :: de_dboundary_second_order
!     REAL( DP ) :: energy_second_order
!     !
!     REAL( DP ) :: charge = 0.0_DP
!     !
!----------------------------------------------------------------------------
!  END TYPE environ_electrolyte
!----------------------------------------------------------------------------
!
! Authors: Oliviero Andreussi (Department of Physics, UNT)
!          Francesco Nattino  (THEOS and NCCR-MARVEL, EPFL)
!          Nicola Marzari     (THEOS and NCCR-MARVEL, EPFL)
!
!--------------------------------------------------------------------
MODULE utils_electrolyte
!--------------------------------------------------------------------
  !
  USE environ_types
  USE environ_output
  !
  USE utils_boundary
  USE environ_base, ONLY : e2
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC :: create_environ_electrolyte, init_environ_electrolyte_first, &
       & init_environ_electrolyte_second, destroy_environ_electrolyte,  &
       & update_environ_electrolyte, electrolyte_of_potential,          &
       & calc_eelectrolyte, calc_deelectrolyte_dboundary
  !
CONTAINS
!--------------------------------------------------------------------
  SUBROUTINE create_environ_electrolyte( electrolyte )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    TYPE( environ_electrolyte ), INTENT(INOUT) :: electrolyte
    !
    CHARACTER( LEN=80 ) :: sub_name = 'create_environ_electrolyte'
    CHARACTER( LEN=80 ) :: label
    !
    label = 'electrolyte'
    CALL create_environ_boundary( electrolyte%boundary, label )
    CALL create_environ_density( electrolyte%density, label )
    label = 'gamma'
    CALL create_environ_density( electrolyte%gamma  , label )
    label = 'dgamma'
    CALL create_environ_density( electrolyte%dgamma  , label )
    !
    electrolyte % charge = 0.D0
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE create_environ_electrolyte
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE init_environ_electrolyte_first( ntyp, mode, stype, rhomax, rhomin, &
       & tbeta, const, alpha, softness, distance, spread, solvent_radius, radial_scale, &
       & radial_spread, filling_threshold, filling_spread, electrons, ions, system, &
       & temperature, cbulk, cionmax, radius, z, electrolyte_entropy, ion_adsorption, &
       & adsorption_energy, linearized, electrolyte )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    LOGICAL, INTENT(IN) :: linearized
    INTEGER, INTENT(IN) :: ntyp, stype
    CHARACTER( LEN=80 ), INTENT(IN) :: mode, electrolyte_entropy, ion_adsorption
    REAL( DP ), INTENT(IN) :: rhomax, rhomin, tbeta, const, distance, spread, alpha, softness, temperature
    REAL( DP ), INTENT(IN) :: solvent_radius, radial_scale, radial_spread, filling_threshold, filling_spread
    REAL( DP ), INTENT(IN) :: cionmax, radius, adsorption_energy
    REAL( DP ), DIMENSION(ntyp), INTENT(IN) :: cbulk, z
    TYPE( environ_electrons ), INTENT(IN) :: electrons
    TYPE( environ_ions ), INTENT(IN) :: ions
    TYPE( environ_system ), TARGET, INTENT(IN) :: system
    TYPE( environ_electrolyte ), INTENT(INOUT) :: electrolyte
    !
    INTEGER :: ityp
    REAL( DP ) :: neutral, sumcbulk, amin, amax
    CHARACTER( LEN=80 ) :: sub_name = 'init_environ_electrolyte_first'
    CHARACTER( LEN=80 ) :: ityps, label
    !
    ! ... Fix Stern boundary parameter and initiate electrolyte boundary
    !
!    IF (( mode == 'electronic' .OR. mode == 'full') .AND. &
!          distance .GT. 0.D0 .AND. spread .GT. 0.D0 ) THEN
!       !
!       ! ... rhomin, rhomax defined as in Ringe et al. JCTC 12, 4052 (2016)
!       !
!       amin = LOG( rhomin ) + ( LOG( rhomin ) - LOG( rhomax )) * distance
!       amax = LOG( rhomax ) + ( LOG( rhomin ) - LOG( rhomax )) * distance
!       !
!       rhomin_ = EXP( amin + ( amax - amin ) * ( 1.D0 - spread) * 0.5D0 )
!       rhomax_ = EXP( amax - ( amax - amin ) * ( 1.D0 - spread) * 0.5D0 )
!       !
!    ELSE
!       !
!       rhomax_ = rhomin
!       rhomin_ = rhopb
!       !
!    END IF
    !
    CALL init_environ_boundary_first( .TRUE., .TRUE., .FALSE., mode, stype, &
         & rhomax, rhomin, tbeta, const, alpha, softness, distance, spread, &
         & solvent_radius, radial_scale, radial_spread, filling_threshold, &
         & filling_spread, electrons, ions, system, electrolyte%boundary )
    !
    ! ... Setup all electrolyte parameters (with checks)
    !
    electrolyte%initialized = .FALSE.
    electrolyte%linearized = linearized
    electrolyte%ntyp = ntyp
    electrolyte%electrolyte_entropy = TRIM( electrolyte_entropy )
    ! maybe unify these in stern_interaction : exclusion, repulsion, full.
    electrolyte%ion_adsorption = TRIM( ion_adsorption )
    electrolyte%temperature = temperature
    electrolyte%cionmax = 0.D0
    electrolyte%permittivity = const 
    !
    IF ( electrolyte%ion_adsorption .NE. 'none' ) THEN
       !
       ! ... Setup exponential function
       !
       electrolyte%adsorption_energy = adsorption_energy
       electrolyte%function%type = 3
       electrolyte%function%pos => system%pos
!      electrolyte%function%volume = 1.D0 ! 3.D0 * k_boltzmann_ry * electrolyte%temperature
       electrolyte%function%dim = system%dim
       electrolyte%function%axis = system%axis
       electrolyte%function%width = distance
       electrolyte%function%spread = spread
       !
    END IF
    !
    IF ( ALLOCATED( electrolyte%ioncctype ) ) &
       CALL errore(sub_name,'Trying to allocate an already allocated object',1)
    ALLOCATE( electrolyte%ioncctype(ntyp) )
    !
    neutral = 0.D0
    !
    DO ityp = 1, ntyp
       !
       electrolyte%ioncctype(ityp)%index = ityp
       !
       ! ... Convert bulk concentrations in atomic units
       electrolyte%ioncctype(ityp)%cbulk = cbulk(ityp) * bohr_radius_si**3 / amu_si
       !
       electrolyte%ioncctype(ityp)%z = -z(ityp)
       neutral = neutral + cbulk(ityp)*z(ityp)
       !
       ! ... Create density for the local electrolyte concentration and related quantities
       WRITE(ityps,'(I2.2)') ityp
       label = 'c_electrolyte_'//TRIM(ityps)
       CALL create_environ_density( electrolyte%ioncctype(ityp)%c, label)
       label = 'cfactor_electrolyte_'//TRIM(ityps)
       CALL create_environ_density( electrolyte%ioncctype(ityp)%cfactor, label)
       !
       IF ( electrolyte%ion_adsorption .NE. 'none' ) THEN
         label = 'potential_'//TRIM(ityps)
         CALL create_environ_density( electrolyte%ioncctype(ityp)%potential, label)
       END IF
       !
    END DO
    !
    IF ( neutral .GT. 1.D-8 ) CALL errore(sub_name,'Bulk electrolyte is not neutral',1)
    !
    ! ... If cionmax is not provided in input but radius is, calculate cionmax
    electrolyte%cionmax = cionmax * bohr_radius_si**3 / amu_si
    IF ( cionmax .EQ. 0.D0 .AND. radius .GT. 0.D0 ) &
            & electrolyte%cionmax = 0.64D0 * 3.D0 / fpi / radius**3
    !
    ! ... Check suitability of cionmax value
    sumcbulk = SUM( electrolyte%ioncctype(:)%cbulk )
    IF ( electrolyte%cionmax .GT. 0.D0 .AND. electrolyte%cionmax .LE. sumcbulk ) &
          & CALL errore( sub_name,'cionmax should be larger than the sum of cbulks',1)
    !
    electrolyte % energy_second_order = 0.D0
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE init_environ_electrolyte_first
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE init_environ_electrolyte_second( cell, electrolyte )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    TYPE( environ_cell ), INTENT(IN) :: cell
    TYPE( environ_electrolyte ), INTENT(INOUT) :: electrolyte
    !
    INTEGER    :: ityp
    REAL( DP ) :: sum_cz2, arg, kT, e
    !
    sum_cz2 = 0.D0
    kT      = k_boltzmann_ry * electrolyte%temperature
    !
    CALL init_environ_boundary_second( cell, electrolyte%boundary )
    !
    CALL init_environ_density( cell, electrolyte%gamma )
    CALL init_environ_density( cell, electrolyte%dgamma )
    CALL init_environ_density( cell, electrolyte%density )
    !
    DO ityp = 1, electrolyte%ntyp
      !
      CALL init_environ_density( cell, electrolyte%ioncctype(ityp)%c )
      CALL init_environ_density( cell, electrolyte%ioncctype(ityp)%cfactor )
      !
      sum_cz2 = sum_cz2 + electrolyte%ioncctype(ityp)%cbulk * electrolyte%ioncctype(ityp)%z ** 2
      !
      IF ( electrolyte%ion_adsorption .NE. 'none' ) THEN
        !
        CALL init_environ_density( cell, electrolyte%ioncctype(ityp)%potential )
        !
      END IF
      !
    END DO
    !
    ! k ** 2 = eps / lambda_D ** 2
    electrolyte%k2 = sum_cz2 / kT
    electrolyte%k2 = electrolyte%k2 * e2 * fpi
    !
    IF ( electrolyte % linearized ) THEN
      CALL init_environ_density( cell, electrolyte % de_dboundary_second_order )
    END IF
    !
    electrolyte%initialized = .TRUE.
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE init_environ_electrolyte_second
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE update_environ_electrolyte( electrolyte )
!--------------------------------------------------------------------
    !
    USE utils_functions, ONLY: density_of_functions
    !
    IMPLICIT NONE
    !
    TYPE( environ_electrolyte ), INTENT(INOUT) :: electrolyte
    !
    INTEGER :: ityp
    !
    CALL start_clock( 'electrolyte' )
    !
    ! ... Check if the boundary is under update (status = 1) or has been fully updated (status = 2)
    IF ( electrolyte % boundary % update_status .GT. 0 ) electrolyte % update = .TRUE.
    !
    IF ( electrolyte % update ) THEN
       ! ... Update the electrolyte in space if the boundary is ready
       IF ( electrolyte % boundary % update_status .EQ. 2 ) THEN
          !
          CALL electrolyte_of_boundary( electrolyte )
          !
          electrolyte % update = .FALSE.
          !
          SELECT CASE ( electrolyte%ion_adsorption )
          CASE( 'anion' )
             DO ityp = 1, electrolyte%ntyp
                IF ( electrolyte%ioncctype(ityp)%z .GT. 0.D0 ) THEN
                   CALL density_of_functions( electrolyte%function, &
                                              electrolyte%ioncctype(ityp)%potential, .TRUE. )
                   electrolyte%ioncctype(ityp)%potential%of_r(:) = electrolyte%adsorption_energy * &
                                        ( electrolyte%ioncctype(ityp)%potential%of_r(:)**2 - &
                                          electrolyte%ioncctype(ityp)%potential%of_r(:)*2.D0 )
                ELSE
                   electrolyte%ioncctype(ityp)%potential%of_r(:) = 0.D0
                END IF
!                CALL print_environ_density( electrolyte%ioncctype(ityp)%potential )
             END DO
          CASE( 'cation' )
             DO ityp = 1, electrolyte%ntyp
                IF ( electrolyte%ioncctype(ityp)%z .LT. 0.D0 ) THEN
                   CALL density_of_functions( electrolyte%function, &
                                              electrolyte%ioncctype(ityp)%potential, .TRUE. )
                   electrolyte%ioncctype(ityp)%potential%of_r(:) = electrolyte%adsorption_energy * &
                                        ( electrolyte%ioncctype(ityp)%potential%of_r(:)**2 - &
                                          electrolyte%ioncctype(ityp)%potential%of_r(:)*2.D0 )
                ELSE
                   electrolyte%ioncctype(ityp)%potential%of_r(:) = 0.D0
                END IF
             END DO
          CASE( 'repulsion' )
             DO ityp = 1, electrolyte%ntyp
                CALL density_of_functions( electrolyte%function, &
                                           electrolyte%ioncctype(ityp)%potential, .TRUE. )
!                CALL print_environ_density( electrolyte%ioncctype(ityp)%potential )
             END DO
          END SELECT
          !
       END IF
       !
    END IF
    !
    CALL stop_clock( 'electrolyte' )
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE update_environ_electrolyte
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE electrolyte_of_boundary( electrolyte )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    TYPE( environ_electrolyte ), TARGET, INTENT(INOUT) :: electrolyte
    !
    TYPE( environ_density ), POINTER :: gam, dgam, scaled
    !
    CHARACTER ( LEN=80 ) :: sub_name = 'electrolyte_of_boundary'
    !
    ! ... Aliases and init local variables
    gam    => electrolyte % gamma
    dgam   => electrolyte % dgamma
    scaled => electrolyte % boundary % scaled
    !
    ! ... Compute gamma(r)
    gam % of_r = 1.D0 - scaled % of_r
    dgam % of_r = -1.D0
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE electrolyte_of_boundary
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE electrolyte_of_potential( potential, electrolyte )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    TYPE( environ_density ),     TARGET, INTENT(IN)    :: potential
    TYPE( environ_electrolyte ), TARGET, INTENT(INOUT) :: electrolyte
    !
    REAL( DP ), POINTER :: z, cbulk
    REAL( DP ), DIMENSION(:), POINTER :: pot, rho, c, cfactor, gam
    !
    REAL( DP )              :: kT, e, factor, sumcbulk, arg
    INTEGER                 :: ityp, ir
    TYPE( environ_density ) :: denominator
    CHARACTER ( LEN=80 )    :: sub_name = 'calc_electrolyte_density'
    !
    REAL( DP ), PARAMETER   :: exp_arg_limit = 40.D0 !LOG( HUGE(1.0_DP) )
    !
    gam => electrolyte%gamma%of_r
    pot => potential%of_r
    rho => electrolyte%density%of_r
    !
    rho = 0.D0
    kT  = k_boltzmann_ry * electrolyte%temperature
    !
    CALL init_environ_density( potential%cell, denominator )
    denominator%of_r = 1.D0
    !
    DO ityp = 1, electrolyte%ntyp
      !
      cfactor => electrolyte%ioncctype(ityp)%cfactor%of_r
      cbulk => electrolyte%ioncctype(ityp)%cbulk
      z => electrolyte%ioncctype(ityp)%z
      !
      cfactor  = 1.D0
      !
      IF ( electrolyte % linearized ) THEN
         ! Linearized PB and Modified PB
         cfactor = 1.D0 - z*pot / kT
         !
         IF ( electrolyte % cionmax .GT. 0.D0 ) THEN
            !
            factor = cbulk / electrolyte%cionmax
            !
            IF ( electrolyte % electrolyte_entropy == 'ions' ) THEN
               ! Linearized Modified PB with electrolyte entropy 'ions'
               denominator%of_r = denominator%of_r - factor * ( 1.D0 - gam )
               !
            END IF
            !
         END IF
         !
      ELSE
         ! PB
         DO ir = 1, potential % cell % ir_end
            ! numerical problems arise when computing exp( -z*pot/kT )
            ! in regions close to the nuclei (exponent is too large).
            arg = - z * pot(ir) / kT
            IF ( electrolyte%ion_adsorption .NE. 'none' ) &
               & arg = arg - electrolyte%ioncctype(ityp)%potential%of_r(ir) / kT
            IF ( arg .GT. exp_arg_limit ) THEN
                cfactor (ir) = EXP( exp_arg_limit )
            ELSE IF ( arg .LT. -exp_arg_limit ) THEN
                cfactor (ir) = EXP( -exp_arg_limit )
            ELSE
                cfactor (ir) = EXP( arg )
            END IF
         END DO
         !
         IF ( electrolyte % cionmax .NE. 0.D0 ) THEN
            ! Modified PB
            factor = cbulk / electrolyte%cionmax
            !
            SELECT CASE ( electrolyte % electrolyte_entropy )
               !
            CASE ( 'full' )
               !
               denominator%of_r = denominator%of_r - factor * ( 1.D0 - cfactor )
               !
            CASE ( 'ions' )
               !
               denominator%of_r = denominator%of_r - factor * ( 1.D0 - gam * cfactor )
               !
            END SELECT
            !
         END IF
         !
      END IF
      !
      NULLIFY( cfactor )
      NULLIFY( cbulk )
      NULLIFY( z )
      !
   END DO
   !
   DO ityp = 1, electrolyte%ntyp
      !
      c => electrolyte%ioncctype(ityp)%c%of_r
      cfactor => electrolyte%ioncctype(ityp)%cfactor%of_r
      cbulk => electrolyte%ioncctype(ityp)%cbulk
      z => electrolyte%ioncctype(ityp)%z
      !
      c   = gam * cbulk * cfactor / denominator%of_r
      rho = rho + c * z
      !
      NULLIFY( c )
      NULLIFY( cfactor )
      NULLIFY( cbulk )
      NULLIFY( z )
      !
   END DO
   !
   electrolyte % charge = integrate_environ_density( electrolyte % density )
   !
   ! Compute energy and de_dboundary terms that depend on the potential.
   ! These are the second order terms, first order terms are equal to zero.
   !
   IF ( electrolyte % linearized ) THEN
      ! The energy term cancels the electrostatic interaction of the electrolyte, and is
      !  identical for the various implementations of the Stern entropy.
      electrolyte % energy_second_order = &
           &  0.5D0 * scalar_product_environ_density( electrolyte % density, potential )
      !
      IF ( electrolyte % electrolyte_entropy == 'ions' .AND. electrolyte%cionmax .GT. 0.D0 ) THEN
         !
         sumcbulk = SUM( electrolyte%ioncctype(:)%cbulk )
         electrolyte % de_dboundary_second_order % of_r = &
              & -0.5D0 * electrolyte%k2 / e2 / fpi * ( 1.D0 - sumcbulk / electrolyte%cionmax ) * &
              & ( pot / denominator%of_r ) ** 2 * electrolyte % dgamma % of_r
         !
      ELSE
         !
         electrolyte % de_dboundary_second_order % of_r = &
              & -0.5D0 * electrolyte%k2 / e2 / fpi * pot * pot * electrolyte % dgamma % of_r
         !
      END IF
      !
   END IF
   !
   CALL destroy_environ_density( denominator )
   !
   RETURN
   !
!--------------------------------------------------------------------
 END SUBROUTINE electrolyte_of_potential
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE calc_eelectrolyte( electrolyte, energy )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    TYPE( environ_electrolyte ), INTENT(IN)  :: electrolyte
    REAL( DP ),                  INTENT(OUT) :: energy
    !
    REAL( DP )              :: kT, sumcbulk, logterm, integral
    INTEGER                 :: ityp
    TYPE( environ_density ) :: arg, f
    !
    energy = 0.D0
    !
    kT       = k_boltzmann_ry * electrolyte%temperature
    sumcbulk = SUM( electrolyte%ioncctype(:)%cbulk )
    !
    CALL init_environ_density( electrolyte%gamma%cell, arg )
    CALL init_environ_density( electrolyte%gamma%cell, f )
    !
    IF ( electrolyte % linearized ) THEN
       !
       IF ( electrolyte % cionmax .EQ. 0.D0 ) THEN
          ! Linearized PB
          integral = integrate_environ_density( electrolyte%gamma )
          energy   = kT * sumcbulk * ( electrolyte%gamma%cell%omega - integral )
          !
       ELSE
          ! Linearized Modified PB
          SELECT CASE ( electrolyte % electrolyte_entropy )
             !
          CASE ( 'full' )
             !
             integral = integrate_environ_density( electrolyte%gamma )
             logterm  = LOG( 1.D0 - sumcbulk / electrolyte%cionmax )
             !
             energy   = kT * electrolyte%cionmax * logterm * &
                  & ( integral - electrolyte%gamma%cell%omega )
             !
          CASE ( 'ions' )
             !
             arg%of_r = 1.D0 - sumcbulk/electrolyte%cionmax*(1.D0 - electrolyte%gamma%of_r)
             f%of_r   = LOG( arg%of_r )
             integral = integrate_environ_density( f )
             !
             energy   = - kT * electrolyte%cionmax * integral
!             if (ionode) print *, 'energy electrolyte ions'
             !
          END SELECT
          !
       END IF
       !
       energy = energy + electrolyte % energy_second_order
       !
    ELSE
       !
       arg%of_r = 0.D0
       !
       DO ityp = 1, electrolyte%ntyp
          !
          arg%of_r = arg%of_r + electrolyte%ioncctype(ityp)%cfactor%of_r * &
               & electrolyte%ioncctype(ityp)%cbulk
          !
       END DO
       !
       IF ( electrolyte % cionmax .EQ. 0.D0 ) THEN
          ! PB
          f%of_r   = electrolyte%gamma%of_r * arg%of_r - sumcbulk
          integral = integrate_environ_density( f )
          !
          energy   = - kT * integral
          !
       ELSE
          ! Modified PB
          SELECT CASE ( electrolyte % electrolyte_entropy )
             !
          CASE ( 'full' )
             !
             arg%of_r = arg%of_r / ( electrolyte%cionmax - sumcbulk )
             arg%of_r = arg%of_r + 1.D0
             f%of_r   = electrolyte%gamma%of_r * LOG ( arg%of_r )
             integral = integrate_environ_density( f )
             logterm  = LOG( 1.D0 - sumcbulk / electrolyte%cionmax )
             !
             energy   = - kT * electrolyte%cionmax * &
                  & ( integral + logterm * electrolyte%gamma%cell%omega )
             !
          CASE ( 'ions' )
             !
             arg%of_r = arg%of_r / electrolyte%cionmax * electrolyte%gamma%of_r
             arg%of_r = arg%of_r + 1.D0 - sumcbulk/electrolyte%cionmax
             f%of_r   = LOG ( arg%of_r )
             integral = integrate_environ_density( f )
             !
             energy   = - kT * electrolyte%cionmax * integral
             !
          END SELECT
          !
       END IF
       !
    END IF
    !
    CALL destroy_environ_density( arg )
    CALL destroy_environ_density( f )
    !
!--------------------------------------------------------------------
  END SUBROUTINE calc_eelectrolyte
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE calc_deelectrolyte_dboundary( electrolyte, de_dboundary)
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    TYPE( environ_electrolyte ), TARGET, INTENT(IN)    :: electrolyte
    TYPE( environ_density ),     TARGET, INTENT(INOUT) :: de_dboundary
    !
    REAL( DP ), DIMENSION(:), POINTER :: gam
    !
    REAL( DP )              :: kT, sumcbulk
    INTEGER                 :: ityp
    TYPE( environ_density ) :: arg
    !
    gam => electrolyte%gamma%of_r
    !
    kT       = k_boltzmann_ry * electrolyte%temperature
    sumcbulk = SUM( electrolyte % ioncctype(:)%cbulk )
    !
    IF ( electrolyte % linearized ) THEN
       !
       IF ( electrolyte % cionmax .EQ. 0.D0 ) THEN
          ! Linearized PB
          de_dboundary % of_r = de_dboundary % of_r - &
               & electrolyte % dgamma % of_r * kT * sumcbulk
          !
       ELSE
          ! Linearized Modified PB
          SELECT CASE ( electrolyte % electrolyte_entropy )
             !
          CASE ( 'full' )
             !
             de_dboundary % of_r = de_dboundary % of_r + &
                  & electrolyte % dgamma % of_r * kT * electrolyte%cionmax * &
                  & LOG( 1.D0 - sumcbulk / electrolyte%cionmax )
             !
          CASE ( 'ions' )
             !
             de_dboundary % of_r = de_dboundary % of_r - &
                  & electrolyte % dgamma % of_r * kT * sumcbulk / &
                  & ( 1.D0 - sumcbulk / electrolyte%cionmax * (1.D0 - gam ))
             !
          END SELECT
          !
       END IF
       !
       de_dboundary % of_r = de_dboundary % of_r + &
            & electrolyte % de_dboundary_second_order % of_r
       !
    ELSE
       !
       CALL init_environ_density( de_dboundary%cell, arg )
       !
       arg%of_r = 0.D0
       !
       DO ityp = 1, electrolyte%ntyp
          !
          arg%of_r = arg%of_r + electrolyte%ioncctype(ityp)%cfactor%of_r * &
               & electrolyte%ioncctype(ityp)%cbulk
          !
       END DO
       !
       IF ( electrolyte % cionmax .EQ. 0.D0 ) THEN
          ! PB
          de_dboundary % of_r = de_dboundary % of_r - &
               & electrolyte % dgamma % of_r * kT * arg%of_r
          !
       ELSE
          ! Modified PB
          SELECT CASE ( electrolyte % electrolyte_entropy )
             !
          CASE ( 'full' )
             !
             arg%of_r = arg%of_r / ( electrolyte%cionmax - sumcbulk )
             arg%of_r = arg%of_r + 1.D0
             de_dboundary % of_r = de_dboundary % of_r - &
                  & electrolyte % dgamma % of_r * kT * electrolyte%cionmax * LOG ( arg%of_r )
             !
          CASE ( 'ions' )
             !
             de_dboundary % of_r = de_dboundary % of_r - &
                  & electrolyte % dgamma % of_r * kT * arg%of_r / &
                  & ( 1.D0 - ( sumcbulk - arg%of_r * gam ) / electrolyte%cionmax )
             !
          END SELECT
          !
       END IF
       !
       CALL destroy_environ_density( arg )
       !
    END IF
    !
!--------------------------------------------------------------------
  END SUBROUTINE calc_deelectrolyte_dboundary
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE destroy_environ_electrolyte( lflag, electrolyte )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    LOGICAL, INTENT(IN) :: lflag
    TYPE( environ_electrolyte ), INTENT(INOUT) :: electrolyte
    CHARACTER( LEN=80 ) :: sub_name = 'destroy_environ_electrolyte'
    INTEGER :: ityp
    !
    IF ( electrolyte%initialized ) THEN
       !
       DO ityp = 1, electrolyte%ntyp
         CALL destroy_environ_density( electrolyte%ioncctype(ityp)%c )
         CALL destroy_environ_density( electrolyte%ioncctype(ityp)%cfactor )
         IF ( electrolyte%ion_adsorption .NE. 'none' ) THEN
           CALL destroy_environ_density( electrolyte%ioncctype(ityp)%potential )
         END IF
       END DO
       !
       CALL destroy_environ_boundary( lflag, electrolyte%boundary )
       CALL destroy_environ_density( electrolyte%gamma   )
       CALL destroy_environ_density( electrolyte%dgamma  )
       CALL destroy_environ_density( electrolyte%density )
       !
       IF ( electrolyte % linearized ) THEN
          CALL destroy_environ_density( electrolyte%de_dboundary_second_order )
       END IF
       !
       electrolyte%initialized = .FALSE.
       !
    END IF
    !
    IF ( lflag ) THEN
       ! These components were allocated first, destroy only if lflag = .TRUE.
       IF ( .NOT. ALLOCATED( electrolyte%ioncctype ) ) &
            & CALL errore(sub_name,'Trying to destroy a non allocated object',1)
       DEALLOCATE( electrolyte%ioncctype )
    ENDIF
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE destroy_environ_electrolyte
!--------------------------------------------------------------------
!----------------------------------------------------------------------------
END MODULE utils_electrolyte
!----------------------------------------------------------------------------
