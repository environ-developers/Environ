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
! This module contains the main drivers and routines to compute the
! electrostatic potential that is the solution of the full
! Poisson-Boltzmann equation, possibly in a dielectric medium:
!
! \nabla \cdot \epsilon (r) \nabla \phi = -4 \pi (\rho + \gamma \sum_i z_i c_i(\phi) )
!
! Authors: Francesco Nattino  (THEOS and NCCR-MARVEL, EPFL)
!          Oliviero Andreussi (Department of Physics, UNT)
!          Nicola Marzari     (THEOS and NCCR-MARVEL, EPFL)
!
!----------------------------------------------------------------------------
MODULE problem_pb
!----------------------------------------------------------------------------
  !
  USE environ_types
  USE electrostatic_types
  USE environ_output
  USE problem_linearized_pb, ONLY : linearized_pb_gradient!, linearized_pb_gradient_sqrt
  USE problem_generalized, ONLY : generalized_gradient
  USE problem_poisson, ONLY : poisson_direct, poisson_energy
  USE environ_base, ONLY : e2, oldenviron
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC :: pb_nested, pb_energy
  !
  INTERFACE pb_nested
     MODULE PROCEDURE pb_nested_charges, pb_nested_density
  END INTERFACE pb_nested
  !
CONTAINS
!--------------------------------------------------------------------
  SUBROUTINE pb_nested_charges( solver, core, charges, potential, inner_setup )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    TYPE( electrostatic_solver ), INTENT(IN) :: solver
    TYPE( electrostatic_core ), INTENT(IN) :: core
    TYPE( environ_charges ), INTENT(IN) :: charges
    TYPE( environ_density ), INTENT(INOUT) :: potential
    TYPE( electrostatic_setup ), OPTIONAL, INTENT(IN) :: inner_setup
    !
    CHARACTER( LEN=25 ) :: sub_name = 'pb_nested_charges'
    !
    CALL start_clock( 'calc_vpb' )
    !
    SELECT CASE( solver % type ) 
    !
    CASE( 'iterative' )
       !
       IF ( solver % use_iterative ) THEN
          ! 
          IF ( solver % auxiliary .EQ. 'ioncc' ) THEN
             !
             IF ( ASSOCIATED( charges % dielectric ) ) THEN
                !
                IF (.NOT. PRESENT(inner_setup) ) & 
                   & CALL errore( sub_name, 'missing setup for inner generalized problem',1)
                !
                CALL pb_iterative( solver%iterative, core, potential, charges%density, charges%electrolyte, & 
                       charges%dielectric, inner_setup%solver, inner_setup%core )
                !
             ELSE 
                !
                CALL pb_iterative( solver%iterative, core, potential, charges%density, charges%electrolyte )
                !
             END IF
             !
          ELSE
          !
          CALL errore(sub_name,'Option not available',1)
          !
          END IF
          !
       END IF
       !
    CASE( 'newton' )
       !
       IF (.NOT. PRESENT(inner_setup) ) CALL errore( sub_name, 'missing inner setup',1)
       !
       IF ( ASSOCIATED( charges % dielectric ) ) THEN
          !
          CALL pb_newton( solver%newton, core, inner_setup%solver, inner_setup%core, &
                          potential, charges%density, charges%electrolyte, charges%dielectric )
          !
       ELSE 
          !
          CALL pb_newton( solver%newton, core, inner_setup%solver, inner_setup%core, &
                          potential, charges%density, charges%electrolyte ) 
          !
       END IF
       !
    CASE( 'lbfgs' )
       !
       CALL errore( sub_name, 'Option not yet implemented', 1 )
       !
    CASE DEFAULT
       !
       CALL errore( sub_name, 'Option not yet implemented', 1 )
       !
    END SELECT
    !
    CALL stop_clock( 'calc_vpb' )
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE pb_nested_charges
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE pb_nested_density( solver, core, charges, potential, electrolyte, dielectric, inner_setup )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    TYPE( electrostatic_solver ), INTENT(IN) :: solver
    TYPE( electrostatic_core ), INTENT(IN) :: core
    TYPE( environ_density ), INTENT(IN) :: charges
    TYPE( environ_density ), INTENT(INOUT) :: potential
    TYPE( environ_electrolyte ), INTENT(IN) :: electrolyte
    TYPE( environ_dielectric ), OPTIONAL, INTENT(IN) :: dielectric
    TYPE( electrostatic_setup ), OPTIONAL, INTENT(IN) :: inner_setup
    !
    CHARACTER( LEN=25 ) :: sub_name = 'pb_nested_density'
    !
    CALL start_clock( 'calc_vpb' )
    !
    SELECT CASE( solver % type ) 
    !
    CASE( 'iterative' )
       !
       IF ( solver % use_iterative ) THEN
          ! 
          IF ( solver % auxiliary .EQ. 'ioncc' ) THEN
             !
             IF ( PRESENT(dielectric) ) THEN
                !
                IF (.NOT. PRESENT(inner_setup) ) CALL errore( sub_name, 'missing inner setup',1)
                !
                CALL pb_iterative( solver%iterative, core, potential, charges, electrolyte, & 
                       dielectric, inner_setup%solver, inner_setup%core )
                !
             ELSE 
                !
                CALL pb_iterative( solver%iterative, core, potential, charges, electrolyte )
                !
             END IF
             !
          ELSE
          !
          CALL errore(sub_name,'Option not available',1)
          !
          END IF
          !
       END IF
       !
    CASE( 'newton' )
       !
       IF (.NOT. PRESENT(inner_setup) ) CALL errore( sub_name, 'missing inner setup',1)
       !
       IF ( PRESENT(dielectric) ) THEN
          !
          CALL pb_newton( solver%newton, core, inner_setup%solver, inner_setup%core, &
                          potential, charges, electrolyte, dielectric )
          !
       ELSE 
          !
          CALL pb_newton( solver%newton, core, inner_setup%solver, inner_setup%core, &
                          potential, charges, electrolyte ) 
          !
       END IF
       !
    CASE( 'lbfgs' )
       !
       CALL errore( sub_name, 'Option not yet implemented', 1 )
       !
    CASE DEFAULT
       !
       CALL errore( sub_name, 'Option not yet implemented', 1 )
       !
    END SELECT
    !
    CALL stop_clock( 'calc_vpb' )
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE pb_nested_density
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE pb_energy( core, charges, potential, energy )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    TYPE( electrostatic_core ),  INTENT(IN)    :: core
    TYPE( environ_charges ),     INTENT(IN)    :: charges
    TYPE( environ_density ),     INTENT(IN)    :: potential
    REAL( DP ),                  INTENT(OUT)   :: energy
    !
    REAL( DP )                                 :: degauss, eself, eions, aux
    CHARACTER( LEN=25 )                        :: sub_name = 'pb_energy'
    !
    CALL start_clock( 'calc_epb' )
    !
    ! Aliases and sanity checks
    !
    IF ( .NOT. ASSOCIATED( charges % density % cell, potential % cell ) ) &
         & CALL errore(sub_name,'Missmatch in charges and potential domains',1)
    !
    energy     = 0.D0
    eions      = 0.D0
    !
    ! Electrostatic interaction of solute
    !
    CALL poisson_energy( core, charges, potential, energy )
    !
    ! Electrostatic interaction of electrolyte
    !
    eions = -0.5D0 * scalar_product_environ_density( charges % electrolyte % density, potential )
    !
    ! Adding correction for point-like nuclei: only affects simulations of charged
    ! systems, it does not affect forces, but shift the energy depending on the
    ! fictitious Gaussian spread of the nuclei
    !
    eself = 0.D0
    degauss = 0.D0
    !
    IF ( charges % include_ions .AND. charges % ions % use_smeared_ions ) THEN
       !
       ! Compute spurious self-polarization energy
       !
       eself = 0.D0
       !
       IF ( core % use_qe_fft ) THEN
          !
          IF ( core % qe_fft % use_internal_pbc_corr .OR. core % need_correction ) THEN
             !
             degauss = 0.D0
             !
          ELSE
             !
             aux = charges % electrolyte % charge
             IF ( ASSOCIATED( charges % dielectric ) ) aux = aux + charges % dielectric % charge
             degauss = - charges % ions % quadrupole_correction * aux * e2 * pi &
                  & / charges % density % cell % omega
             !
          ENDIF
          !
       ENDIF
       !
    ENDIF
    !
    energy = energy + eions + eself + degauss
    !
    CALL stop_clock( 'calc_epb' )
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE pb_energy
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE pb_iterative( iterative, core, potential, charges, & 
                   electrolyte, dielectric, inner_solver, inner_core ) 
!--------------------------------------------------------------------
!    USE mdiis,     ONLY : mdiis_type, allocate_mdiis, deallocate_mdiis, update_by_mdiis
    !
    IMPLICIT NONE
    !
    TYPE( iterative_solver ), TARGET, INTENT(IN) :: iterative
    TYPE( electrostatic_core ), INTENT(IN) :: core
    TYPE( environ_density ), TARGET, INTENT(IN) :: charges
    TYPE( environ_electrolyte), TARGET, INTENT(IN) :: electrolyte
    TYPE( environ_dielectric ), OPTIONAL, INTENT(IN) :: dielectric
    TYPE( electrostatic_solver ), OPTIONAL, INTENT(IN) :: inner_solver
    TYPE( electrostatic_core ), OPTIONAL, INTENT(IN) :: inner_core
    TYPE( environ_density ), TARGET, INTENT(INOUT) :: potential
    !
    TYPE( environ_cell ), POINTER :: cell
    TYPE( environ_density ), POINTER :: x, gam
    !
    INTEGER :: iter, ityp, ir
    REAL( DP ) :: total, totaux, delta_qm, delta_en, kT, factor, arg
    TYPE( environ_density ) :: residual, rhotot, denominator, rhoaux, cfactor!, rhonew
    !
    CHARACTER( LEN=80 ) :: sub_name = 'pb_nested'
!    TYPE(mdiis_type) :: mdiist
    !
    INTEGER, POINTER :: maxiter, ir_end
    REAL( DP ), POINTER :: tolrhoaux, mix, cbulk, cionmax, z
    !
    REAL( DP ), PARAMETER     :: exp_arg_limit = 25.D0 !LOG( HUGE(1.0_DP) )
    !
    IF ( verbose .GE. 1 .AND. ionode ) WRITE(environ_unit,9000)
9000 FORMAT(/,4('%'),' OUTER LOOP ON POTENTIAL ',50('%'))
    !
    ! ... Various checks
    !
    IF ( PRESENT(dielectric) ) THEN
       IF ( .NOT. ASSOCIATED(charges%cell,dielectric%epsilon%cell) ) &
         & CALL errore(sub_name,'Inconsistent cells of input fields',1)
       IF ( .NOT. (PRESENT(inner_solver) .OR. PRESENT(inner_core)) ) &
         & CALL errore(sub_name,'Missing inner solver',1)
    END IF
    IF ( .NOT. ASSOCIATED(charges%cell,electrolyte%gamma%cell) ) &
         & CALL errore(sub_name,'Inconsistent cells of input fields',1)
    IF ( .NOT. ASSOCIATED(charges%cell,potential%cell) ) &
         & CALL errore(sub_name,'Inconsistent cells for charges and potential',1)
    !
    ! ... Aliases
    !
    cell      => charges%cell
    ir_end    => cell % ir_end
    maxiter   => iterative % maxiter
    mix       => iterative % mix
    tolrhoaux => iterative % tol
    x         => potential
    cionmax => electrolyte%cionmax
    gam => electrolyte % gamma
    !
    ! ... Create local variables
    !
    kT  = k_boltzmann_ry * electrolyte%temperature
    CALL init_environ_density( cell, denominator )
    CALL init_environ_density( cell, rhoaux )
    CALL init_environ_density( cell, rhotot )
    CALL init_environ_density( cell, residual )
    CALL init_environ_density( cell, cfactor )
!    IF ( iterative % mix_type == 'diis' ) & 
!        CALL init_environ_density( cell, rhonew )
    !
    ! ... Starting guess from new input and previous solution(s)
    !
    IF ( x % update ) THEN
       !
       rhoaux % of_r = electrolyte % density % of_r
       !
    ELSE 
       !
       x % of_r = 0.D0
       rhoaux % of_r = 0.D0
       x % update = .TRUE.
       !
    END IF
    !
    ! ... Start iterative algorithm 
    !
!    IF ( iterative % mix_type == 'diis' ) & 
!       CALL allocate_mdiis( mdiist, iterative % ndiis, SIZE(x%of_r), mix, 1)
    !
    DO iter = 1, maxiter
       !
       IF ( verbose .GE. 1 .AND. ionode ) WRITE(environ_unit,9002) iter
9002   FORMAT(' Outer loop iteration # ',i10)
       !
       rhotot % of_r = charges % of_r + rhoaux % of_r
       !
       IF ( PRESENT(dielectric) ) THEN
          !
          CALL generalized_gradient( inner_solver, inner_core, rhotot, dielectric, x ) 
          !
       ELSE 
          !
          CALL poisson_direct( core, rhotot, x )
          !
       END IF
       !
       !  ... Calculate electrolyte charge
       !
       residual%of_r = 0.D0
       denominator%of_r = 1.D0
       !
       DO ityp = 1, electrolyte % ntyp 
          !
          cbulk => electrolyte%ioncctype(ityp)%cbulk 
          z => electrolyte%ioncctype(ityp)%z
          cfactor % of_r = 1.D0
          !
          DO ir = 1, ir_end 
             arg = -z * x%of_r(ir) / kT
             IF ( electrolyte%ion_adsorption .NE. 'none' ) &
                & arg = arg - electrolyte%ioncctype(ityp)%potential%of_r(ir) / kT
             IF ( arg .LT. exp_arg_limit ) THEN
                cfactor % of_r (ir) = EXP( arg )
             END IF
          END DO
          !
          residual % of_r = residual % of_r + z * cbulk * cfactor % of_r
          !
          IF ( cionmax .GT. 0.D0 ) THEN
             factor = cbulk / cionmax
             SELECT CASE ( electrolyte % stern_entropy )
             CASE ( 'full' )
                denominator%of_r = denominator%of_r - factor * ( 1.D0 - cfactor%of_r )
             CASE ( 'ions' )
                denominator%of_r = denominator%of_r - factor * ( 1.D0 - gam%of_r * cfactor%of_r )
             END SELECT
          END IF
          !
          NULLIFY( z )
          NULLIFY( cbulk )
          !
       END DO
       !
       residual % of_r = gam % of_r * residual % of_r / denominator % of_r
       !
       ! ... residual is now the new electrolyte charge
       !
!       SELECT CASE ( iterative % mix_type ) 
!       CASE ( 'diis' )  
!          !
!          rhonew % of_r = residual % of_r
!          residual % of_r = residual % of_r - rhoaux % of_r
!          CALL update_by_mdiis(mdiist, rhonew % of_r, residual % of_r, residual%cell%comm)
!          rhoaux % of_r = rhonew % of_r
!          !
!       CASE( 'linear' )
          !
          residual % of_r = residual % of_r - rhoaux % of_r
          rhoaux % of_r = rhoaux % of_r + mix * residual % of_r
          ! 
!       CASE DEFAULT
!          !
!          CALL errore(sub_name, 'mix_type not implemented' )
!          !
!       END SELECT
       !
       ! ... If residual is small enough exit
       !
       delta_en = euclidean_norm_environ_density( residual )
       IF ( oldenviron ) delta_en = quadratic_mean_environ_density_old( residual )
       delta_qm = quadratic_mean_environ_density( residual )
       totaux = integrate_environ_density( rhoaux )
       IF ( verbose .GE. 1 .AND. ionode ) WRITE(environ_unit,9004)delta_qm,delta_en,tolrhoaux
9004   FORMAT('outer delta_qm = ',E14.6,' delta_en = ',E14.6,' tol = ',E14.6)
       IF ( verbose .GE. 1 .AND. ionode ) WRITE(environ_unit,9003) totaux
9003   FORMAT(' Total iterative electrolyte charge = ',F13.6)
       IF ( delta_en .LT. tolrhoaux .AND. iter .GT. 0 ) THEN
          IF ( verbose .GE. 1 .AND. ionode ) WRITE(environ_unit,9005)
9005      FORMAT(' Outer loop is converged, EXIT')
          EXIT
       ELSE IF ( iter .EQ. maxiter ) THEN
          IF ( ionode ) WRITE(program_unit,9006)
9006      FORMAT(' Warning: Electrolyte charge not converged')
       ENDIF
    ENDDO
    !
    IF (.not.tddfpt.AND.verbose.GE.1.AND.ionode) WRITE(program_unit, 9007) delta_en, iter
9007 FORMAT('     electrolyte accuracy =',1PE8.1,', # of iterations = ',i3)
    !
    ! ... Destroy local variables
    !
    CALL destroy_environ_density(denominator)
    CALL destroy_environ_density(rhoaux)
    CALL destroy_environ_density(rhotot)
    CALL destroy_environ_density(residual)
    CALL destroy_environ_density(cfactor)
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE pb_iterative
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE pb_newton( newton, core, inner_solver, inner_core, & 
                        potential, charges, electrolyte, dielectric)
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    TYPE( newton_solver ), TARGET, INTENT(IN) :: newton 
    TYPE( electrostatic_solver ), INTENT(IN) :: inner_solver
    TYPE( electrostatic_core ), INTENT(IN) :: core, inner_core
    TYPE( environ_density ), TARGET, INTENT(IN) :: charges
    TYPE( environ_electrolyte), TARGET, INTENT(IN) :: electrolyte
    TYPE( environ_dielectric ), OPTIONAL, INTENT(IN) :: dielectric
    TYPE( environ_density ), TARGET, INTENT(INOUT) :: potential
    !
    TYPE( environ_cell ), POINTER :: cell
    TYPE( environ_density ), POINTER :: x, gam
    !
    INTEGER :: iter, itypi, itypj, ir
    REAL( DP ) :: total, totaux, delta_qm, delta_en, kT, z, arg
    TYPE( environ_density ) :: residual, rhotot, numerator, denominator, &
                               cfactor, rhoaux, screening
    !
    CHARACTER( LEN=80 ) :: sub_name = 'pb_newton'
!!!!!
!    CHARACTER( LEN=80 ) :: label
!!!!!
    !
    INTEGER, POINTER :: maxiter, ir_end
    REAL( DP ), POINTER :: tol, mix, cbulk, cbulki, cbulkj, cionmax, zi, zj
    !
    !
    REAL( DP ), PARAMETER     :: exp_arg_limit = 25.D0 !LOG( HUGE(1.0_DP) )
    !
    IF ( verbose .GE. 1 .AND. ionode ) WRITE(environ_unit,9000)
9000 FORMAT(/,4('%'),' COMPUTE ELECTROSTATIC POTENTIAL ',43('%'))
    !
    ! ... Check that fields have the same defintion domain
    !
    IF ( PRESENT(dielectric)  .AND. &
         .NOT. ASSOCIATED(charges%cell,dielectric%epsilon%cell) ) &
         & CALL errore(sub_name,'Inconsistent cells of input fields',1)
    IF ( .NOT. ASSOCIATED(charges%cell,electrolyte%gamma%cell) ) &
         & CALL errore(sub_name,'Inconsistent cells of input fields',1)
    IF ( .NOT. ASSOCIATED(charges%cell,potential%cell) ) &
         & CALL errore(sub_name,'Inconsistent cells for charges and potential',1)
    !
    ! ... Aliases
    !
    cell      => charges%cell
    ir_end    => cell % ir_end
    maxiter   => newton % maxiter
    tol       => newton % tol
    cionmax   => electrolyte % cionmax
    x         => potential
    gam       => electrolyte % gamma
    !
    ! ... Create local variables
    !
    kT  = k_boltzmann_ry * electrolyte%temperature
    !
!!!!!
!    label = 'rhoaux'
!    CALL create_environ_density( rhoaux, label )
!    label = 'screening'
!    CALL create_environ_density( screening, label )
!    label = 'numerator'
!    CALL create_environ_density( numerator, label )
!    label = 'denominator'
!    CALL create_environ_density( denominator, label )
!!!!!
    CALL init_environ_density( cell, cfactor )
    CALL init_environ_density( cell, numerator )
    CALL init_environ_density( cell, denominator )
    CALL init_environ_density( cell, rhoaux )
    CALL init_environ_density( cell, rhotot )
    CALL init_environ_density( cell, residual )
    CALL init_environ_density( cell, screening )
    !
    ! ... Starting guess from new input and previous solution(s)
    !
    IF ( x%update ) THEN
       !
       rhoaux % of_r = electrolyte % density % of_r
       screening % of_r = 0.D0
       denominator%of_r = 1.D0
       !
       DO itypi = 1, electrolyte % ntyp
          !
          cbulki => electrolyte%ioncctype(itypi)%cbulk
          zi => electrolyte%ioncctype(itypi)%z
          !
          cfactor % of_r = 1.D0
          !
          DO ir = 1, ir_end
             arg = - zi*x%of_r(ir) /kT
             IF ( electrolyte%ion_adsorption .NE. 'none' ) &
                & arg = arg - electrolyte%ioncctype(itypi)%potential%of_r(ir) / kT
             IF ( ABS(arg) .LT. exp_arg_limit ) THEN
                cfactor % of_r (ir) = EXP( arg )
             END IF
          END DO
          !
          numerator % of_r = 1.D0
          !
          IF ( cionmax .GT. 0.D0 ) THEN
             !
             SELECT CASE ( electrolyte % stern_entropy )
             !
             CASE ( 'full' )
                !
                denominator%of_r = denominator%of_r - cbulki/cionmax * ( 1.D0 - cfactor%of_r)
                !
                DO itypj = 1, electrolyte % ntyp
                   !
                   zj => electrolyte%ioncctype(itypj)%z
                   cbulkj => electrolyte%ioncctype(itypj)%cbulk
                   !
                   IF ( itypj .EQ. itypi ) THEN
                      !
                      numerator % of_r = numerator % of_r - cbulkj/cionmax
                      !
                   ELSE
                      !
                      numerator % of_r = numerator % of_r - cbulkj/cionmax * &
                                 (1.D0 - (1.D0 - zj/zi) * cfactor % of_r **(zj/zi))
                      !
                   END IF
                   !
                   NULLIFY( zj )
                   NULLIFY( cbulkj )
                   !
                END DO
                !
             CASE ( 'ions' )
                !
                denominator%of_r = denominator%of_r - cbulki/cionmax * &
                           ( 1.D0 - gam%of_r * cfactor%of_r )
                !
                DO itypj = 1, electrolyte % ntyp
                   !
                   zj => electrolyte%ioncctype(itypj)%z
                   cbulkj => electrolyte%ioncctype(itypj)%cbulk
                   !
                   IF ( itypj .EQ. itypi ) THEN
                      !
                      numerator % of_r = numerator % of_r - cbulkj/cionmax
                      !
                   ELSE
                      !
                      numerator % of_r = numerator % of_r - cbulkj/cionmax * &
                                  (1.D0 - (1.D0 - zj/zi) * gam%of_r * cfactor % of_r **(zj/zi))
                      !
                   END IF
                   !
                   NULLIFY( zj )
                   NULLIFY( cbulkj )
                   !
                END DO
                 !
             END SELECT
          END IF
          !
          screening % of_r = screening % of_r + &
              cbulki * zi**2 / kT * cfactor % of_r * numerator % of_r
          !
          NULLIFY( zi )
          NULLIFY( cbulki )
          !
       END DO
       !
       screening % of_r = screening % of_r * gam % of_r / denominator % of_r ** 2
       !
    ELSE
       !
       x % update = .TRUE.
       x % of_r = 0.D0
       rhoaux % of_r = 0.D0
       screening % of_r = electrolyte%k2/e2/fpi * gam%of_r
       !
    ENDIF
    !
    residual % of_r = 0.D0
    !
    ! ... Start Newton's algorithm 
    !
    DO iter = 1, maxiter
       !
       IF ( verbose .GE. 1 .AND. ionode ) WRITE(environ_unit,9002) iter
9002   FORMAT(' Newton step # ',i10)
       !
       rhotot % of_r = charges % of_r + rhoaux % of_r + screening % of_r * x%of_r
       residual % of_r = x%of_r
       !
       IF (PRESENT( dielectric )) THEN
          CALL linearized_pb_gradient( inner_solver, inner_core, rhotot, electrolyte, x, dielectric, screening )
       ELSE
          CALL linearized_pb_gradient( inner_solver, inner_core, rhotot, electrolyte, x, screening=screening ) 
       END IF
       !
       residual % of_r = x%of_r - residual % of_r  
       !
       rhoaux % of_r = 0.D0
       screening % of_r = 0.D0
       denominator%of_r = 1.D0
       !
!       IF ( SIZE( electrolyte%ioncctype ) .EQ. 2 )  THEN
       IF (.FALSE.) THEN
          ! Symmetric electrolyte: use hyperbolic functions
          cbulk => electrolyte%ioncctype(1)%cbulk
          z = ABS(electrolyte%ioncctype(1)%z)
          !
          cfactor % of_r = 1.D0
          !
          DO ir = 1, ir_end
             arg = z * x%of_r(ir) / kT
             IF ( ABS(arg) .LT. exp_arg_limit ) THEN
                rhoaux % of_r (ir) = SINH( arg )
                cfactor % of_r (ir) = COSH( arg ) 
             END IF
          END DO
          !
          rhoaux % of_r  = -2.D0 * z * cbulk * rhoaux % of_r
          !
          IF ( cionmax .GT. 0.D0 ) THEN
             !
             SELECT CASE ( electrolyte % stern_entropy )
             !
             CASE ( 'full' )
                !
                denominator%of_r = 2.D0 * cbulk/cionmax * ( cfactor%of_r - 1.D0 ) + 1.D0
                screening % of_r = electrolyte%k2/e2/fpi * &
                   ( cfactor % of_r * ( 1.D0 - 2.D0*cbulk/cionmax ) + 2.D0*cbulk/cionmax )
                !
             CASE ( 'ions' )
                !
                denominator%of_r = 2.D0 * cbulk/cionmax * ( gam%of_r * cfactor%of_r - 1.D0) + 1.D0
                screening % of_r = electrolyte%k2/e2/fpi * & 
                   ( cfactor % of_r * ( 1.D0 - 2.D0*cbulk/cionmax ) + 2.D0*cbulk/cionmax * gam%of_r )
                !
             END SELECT
             !
          ELSE
             !
             screening % of_r = electrolyte%k2/e2/fpi * cfactor%of_r
             !
          END IF
          !
          NULLIFY( cbulk )
          !
       ELSE 
          ! General solution for symmetric & asymmetric electrolyte
          DO itypi = 1, electrolyte % ntyp
             !
             cbulki => electrolyte%ioncctype(itypi)%cbulk
             zi => electrolyte%ioncctype(itypi)%z
             !
             cfactor % of_r = 1.D0 
             !
             DO ir = 1, ir_end
                arg = - zi*x%of_r(ir) /kT
                IF ( electrolyte%ion_adsorption .NE. 'none' ) &
                   & arg = arg - electrolyte%ioncctype(itypi)%potential%of_r(ir) / kT
                IF ( arg .LT. exp_arg_limit ) THEN
                   cfactor % of_r (ir) = EXP( arg )
                END IF
             END DO
             !
             rhoaux % of_r = rhoaux % of_r + zi * cbulki * cfactor % of_r
             !
             numerator % of_r = 1.D0
             !
             IF ( cionmax .GT. 0.D0 ) THEN
                !
                SELECT CASE ( electrolyte % stern_entropy )
                !
                CASE ( 'full' )
                   !
                   denominator%of_r = denominator%of_r - cbulki/cionmax * ( 1.D0 - cfactor%of_r)
                   !
                   DO itypj = 1, electrolyte % ntyp
                      !
                      zj => electrolyte%ioncctype(itypj)%z
                      cbulkj => electrolyte%ioncctype(itypj)%cbulk
                      !
                      IF ( itypj .EQ. itypi ) THEN
                         !
                         numerator % of_r = numerator % of_r - cbulkj/cionmax
                         !
                      ELSE 
                         !
                         numerator % of_r = numerator % of_r - cbulkj/cionmax * & 
                                     (1.D0 - (1.D0 - zj/zi) * cfactor % of_r **(zj/zi))
                         !
                      END IF
                      !
                      NULLIFY( zj )
                      NULLIFY( cbulkj )
                      !
                   END DO
                   !
                CASE ( 'ions' )
                   !
                   denominator%of_r = denominator%of_r - cbulki/cionmax * & 
                              ( 1.D0 - gam%of_r * cfactor%of_r )
                   !
                   DO itypj = 1, electrolyte % ntyp
                      !
                      zj => electrolyte%ioncctype(itypj)%z
                      cbulkj => electrolyte%ioncctype(itypj)%cbulk
                      !
                      IF ( itypj .EQ. itypi ) THEN
                         !
                         numerator % of_r = numerator % of_r - cbulkj/cionmax
                         !
                      ELSE
                         !
                         numerator % of_r = numerator % of_r - cbulkj/cionmax * &
                                     (1.D0 - (1.D0 - zj/zi) * gam%of_r * cfactor % of_r **(zj/zi))
                         !
                      END IF
                      !
                      NULLIFY( zj )
                      NULLIFY( cbulkj )
                      !
                   END DO

                   !
                END SELECT
                !
             END IF
             !
             screening % of_r = screening % of_r + & 
                 cbulki * zi**2 / kT * cfactor % of_r * numerator % of_r
             !
             NULLIFY( zi )
             NULLIFY( cbulki )
             !
          END DO
          !
       END IF 
       !
       rhoaux % of_r = gam % of_r * rhoaux % of_r / denominator % of_r
       screening % of_r = screening % of_r * gam % of_r / denominator % of_r ** 2
!!!!!!
!       CALL print_environ_density( gam, local_verbose=2 )
!       CALL print_environ_density( numerator, local_verbose=2 )
!       CALL print_environ_density( denominator, local_verbose=2 )
!       CALL print_environ_density( screening, local_verbose=2 )
!       CALL print_environ_density( rhoaux, local_verbose=2 )
!       CALL print_environ_density( x, local_verbose=2 )
!       CALL errore( sub_name, 'stop', 1)
!!!!!!
       !
       ! ... If residual is small enough exit
       !
       delta_en = euclidean_norm_environ_density( residual )
       IF ( oldenviron ) delta_en = quadratic_mean_environ_density_old( residual )
       delta_qm = quadratic_mean_environ_density( residual )
       totaux = integrate_environ_density( rhoaux )
       IF ( verbose .GE. 1 .AND. ionode ) WRITE(environ_unit,9004)delta_qm,delta_en,tol
9004   FORMAT('outer delta_qm = ',E14.6,' delta_en = ',E14.6,' tol = ',E14.6)
       IF ( verbose .GE. 1 .AND. ionode ) WRITE(environ_unit,9003) totaux
9003   FORMAT(' Total iterative electrolyte charge = ',F13.6)
       IF ( delta_en .LT. tol .AND. iter .GT. 0 ) THEN
          IF ( verbose .GE. 1 .AND. ionode ) WRITE(environ_unit,9005)
9005      FORMAT(' Newton steps are converged, EXIT')
          EXIT
       ELSE IF ( iter .EQ. maxiter ) THEN
          IF ( ionode ) WRITE(program_unit,9006)
9006      FORMAT(' Warning: Electrolyte charge not converged')
       ENDIF
    ENDDO
    !
    IF (.not.tddfpt.AND.verbose.GE.1.AND.ionode) WRITE(program_unit, 9007) delta_en, iter
9007 FORMAT('     electrolyte accuracy =',1PE8.1,', # of iterations = ',i3)
    !
    ! ... Destroy local variables
    !
    CALL destroy_environ_density(cfactor)
    CALL destroy_environ_density(numerator)
    CALL destroy_environ_density(denominator)
    CALL destroy_environ_density(rhoaux)
    CALL destroy_environ_density(rhotot)
    CALL destroy_environ_density(residual)
    CALL destroy_environ_density(screening) 
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE pb_newton
!--------------------------------------------------------------------
!----------------------------------------------------------------------------
END MODULE problem_pb
!----------------------------------------------------------------------------
