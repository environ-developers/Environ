!
! Copyright (C) Oliviero Andreussi
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------
MODULE linearized_pb
!--------------------------------------------------------------------

  USE environ_types
  USE electrostatic_types
  USE environ_output
  USE poisson,      ONLY : poisson_direct, poisson_energy
  USE environ_base, ONLY : e2, add_jellium
!  USE periodic, ONLY : calc_v0periodic

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: linearized_pb_gradient, linearized_pb_energy

  INTERFACE linearized_pb_gradient
     MODULE PROCEDURE linearized_pb_gradient_charges, linearized_pb_gradient_density
  END INTERFACE linearized_pb_gradient

CONTAINS

!--------------------------------------------------------------------
SUBROUTINE linearized_pb_gradient_charges( solver, core, charges, potential )
!--------------------------------------------------------------------

  IMPLICIT NONE

  TYPE( electrostatic_solver ), INTENT(IN) :: solver
  TYPE( electrostatic_core ), INTENT(IN) :: core
  TYPE( environ_charges ), INTENT(IN) :: charges
  TYPE( environ_density ), INTENT(INOUT) :: potential

  CHARACTER( LEN=25 ) :: sub_name = 'linearized_pb_gradient'

  CALL start_clock( 'calc_vlinpb' )

  IF ( solver % use_gradient ) THEN

     IF ( solver % auxiliary .EQ. 'none' ) THEN

        IF ( ASSOCIATED( charges % dielectric ) ) THEN

           SELECT CASE ( solver % gradient % preconditioner )

           CASE ( 'sqrt' )

              CALL linearized_pb_gradient_sqrt( solver % gradient, core, charges%density, &
                   & charges%dielectric, charges%electrolyte, potential )

           CASE DEFAULT

              CALL errore( sub_name, 'unexpected preconditioner keyword', 1 )

           END SELECT

        ELSE

           CALL linearized_pb_gradient_vacuum( solver % gradient, core, charges%density, &
                & charges % electrolyte, potential )

        END IF

     ELSE

        CALL errore(sub_name,'Option not yet implemented',1)

     END IF

  ELSE

     CALL errore( sub_name, 'unexpected solver keyword', 1 )

  END IF

  CALL stop_clock( 'calc_vlinpb' )

  RETURN

!--------------------------------------------------------------------
END SUBROUTINE linearized_pb_gradient_charges
!--------------------------------------------------------------------
!--------------------------------------------------------------------
SUBROUTINE linearized_pb_gradient_density( solver, core, charges, electrolyte, potential, dielectric )
!--------------------------------------------------------------------

  IMPLICIT NONE

  TYPE( electrostatic_solver ), INTENT(IN) :: solver
  TYPE( electrostatic_core ), INTENT(IN) :: core
  TYPE( environ_density ), INTENT(IN) :: charges
  TYPE( environ_electrolyte ),INTENT(IN) :: electrolyte
  TYPE( environ_density ), INTENT(INOUT) :: potential
  TYPE( environ_dielectric ), INTENT(IN), OPTIONAL :: dielectric

  CHARACTER( LEN=25 ) :: sub_name = 'linearized_pb_gradient'

  CALL start_clock( 'calc_vlinpb' )

  potential % of_r = 0.D0

  IF ( solver % use_gradient ) THEN

     IF ( solver % auxiliary .EQ. 'none' ) THEN

        IF ( PRESENT( dielectric ) ) THEN

           SELECT CASE ( solver % gradient % preconditioner )

           CASE ( 'sqrt' )

              CALL linearized_pb_gradient_sqrt( solver % gradient, core, charges, dielectric, electrolyte, potential )

           CASE DEFAULT

              CALL errore( sub_name, 'unexpected preconditioner keyword', 1 )

           END SELECT

        ELSE

           CALL linearized_pb_gradient_vacuum( solver % gradient, core, charges, electrolyte, potential )

        END IF

     ELSE

        CALL errore(sub_name,'Option not yet implemented',1)

     END IF

  ELSE

     CALL errore( sub_name, 'unexpected solver keyword', 1 )

  END IF

  CALL stop_clock( 'calc_vlinpb' )

  RETURN

!--------------------------------------------------------------------
END SUBROUTINE linearized_pb_gradient_density
!--------------------------------------------------------------------
!--------------------------------------------------------------------
SUBROUTINE linearized_pb_energy( core, charges, potential, energy )
!--------------------------------------------------------------------

  IMPLICIT NONE

  TYPE( electrostatic_core ),  INTENT(IN)    :: core
  TYPE( environ_charges ),     INTENT(IN)    :: charges
  TYPE( environ_density ),     INTENT(IN)    :: potential
  REAL( DP ),                  INTENT(OUT)   :: energy

  REAL( DP )                                 :: degauss, eself, eions, aux
  CHARACTER( LEN=25 )                        :: sub_name = 'linearized_pb_energy'

  CALL start_clock( 'calc_elinpb' )

  ! Aliases and sanity checks

  IF ( .NOT. ASSOCIATED( charges % density % cell, potential % cell ) ) &
       & CALL errore(sub_name,'Missmatch in charges and potential domains',1)

  energy     = 0.D0
  eions      = 0.D0

  ! Electrostatic interaction of solute

  CALL poisson_energy( core, charges, potential, energy )

  ! Electrostatic interaction of electrolyte

  CALL poisson_energy( core, charges % electrolyte % density, potential, eions )

  ! Adding correction for point-like nuclei: only affects simulations of charged
  ! systems, it does not affect forces, but shift the energy depending on the
  ! fictitious Gaussian spread of the nuclei

  eself = 0.D0
  degauss = 0.D0

  IF ( charges % include_ions .AND. charges % ions % use_smeared_ions ) THEN

     ! Compute spurious self-polarization energy

     eself = 0.D0

     IF ( core % use_qe_fft ) THEN

        IF ( core % qe_fft % use_internal_pbc_corr .OR. core % need_correction ) THEN

           degauss = 0.D0

        ELSE

           aux = charges % electrolyte % charge
           IF ( ASSOCIATED( charges % dielectric ) ) aux = aux + charges % dielectric % charge
           degauss = - charges % ions % quadrupole_correction * aux * e2 * pi &
                & / charges % density % cell % omega

        ENDIF

     ENDIF

  ENDIF

  energy = energy + eions + eself + degauss

  CALL stop_clock( 'calc_elinpb' )

  RETURN

!--------------------------------------------------------------------
END SUBROUTINE linearized_pb_energy
!--------------------------------------------------------------------
!--------------------------------------------------------------------
SUBROUTINE linearized_pb_gradient_sqrt( gradient, core, charges, dielectric, electrolyte, potential )
!--------------------------------------------------------------------

  IMPLICIT NONE

  TYPE( gradient_solver ), TARGET, INTENT(IN) :: gradient
  TYPE( electrostatic_core ), INTENT(IN) :: core
  TYPE( environ_density ), TARGET, INTENT(IN) :: charges
  TYPE( environ_dielectric ), TARGET, INTENT(IN) :: dielectric
  TYPE( environ_electrolyte), TARGET, INTENT(IN) :: electrolyte
  TYPE( environ_density ), TARGET, INTENT(INOUT) :: potential

  TYPE( environ_cell ), POINTER :: cell
  TYPE( environ_density ), POINTER :: x, b, eps, factsqrt, gam
  TYPE( environ_gradient ), POINTER :: gradeps

  INTEGER :: iter
  REAL( DP ) :: rznew, rzold, alpha, beta, pAp, delta_qm, delta_en, jellium, shift
  TYPE( environ_density ) :: r, z, p, Ap, invsqrt

  CHARACTER( LEN=80 ) :: sub_name = 'linearized_pb_gradient_sqrt'

  LOGICAL, POINTER :: lconjugate
  INTEGER, POINTER :: maxstep
  REAL( DP ), POINTER :: tolvelect

  lconjugate => gradient % lconjugate
  maxstep => gradient % maxstep
  tolvelect => gradient % tol

  IF ( verbose .GE. 1 ) WRITE(environ_unit,9000)
9000 FORMAT(/,4('%'),' COMPUTE ELECTROSTATIC POTENTIAL ',43('%'))

  ! ... Check that fields have the same defintion domain

  IF ( .NOT. ASSOCIATED(charges%cell,dielectric%epsilon%cell) ) &
       & CALL errore(sub_name,'Inconsistent cells of input fields',1)
  IF ( .NOT. ASSOCIATED(charges%cell,potential%cell) ) &
       & CALL errore(sub_name,'Inconsistent cells for charges and potential',1)
  cell => charges%cell

  ! ... Aliases

  b => charges
  eps => dielectric % epsilon
  factsqrt => dielectric % factsqrt
  x => potential
  gam => electrolyte % gamma
  CALL init_environ_density( cell, invsqrt )
  invsqrt%of_r = 1.D0 / SQRT(eps%of_r)
  jellium = 0.D0
  IF ( add_jellium ) jellium = integrate_environ_density( charges ) / cell % omega

  ! ... Create and initialize local variables

  CALL init_environ_density( cell, r )
  CALL init_environ_density( cell, z )
  CALL init_environ_density( cell, p )
  CALL init_environ_density( cell, Ap )

  ! ... Starting guess from new input and previous solution(s)

  IF ( x%update ) THEN

     r%of_r = ( b%of_r - jellium ) - factsqrt%of_r * x%of_r

     ! ... Preconditioning step

     z%of_r = r%of_r * invsqrt%of_r
     CALL poisson_direct( core, z, z )
     z%of_r = z%of_r * invsqrt%of_r

     rzold = scalar_product_environ_density( r, z )
     IF ( ABS(rzold) .LT. 1.D-30 ) &
          & CALL errore(sub_name,'Null step in gradient descent iteration',1)

     r%of_r = factsqrt%of_r * ( x%of_r - z%of_r )
     delta_en = euclidean_norm_environ_density( r )
     delta_qm = quadratic_mean_environ_density( r )
     IF ( delta_en .LT. 1.D-02 ) THEN
        IF ( verbose .GE. 1 ) WRITE(environ_unit,9008)delta_en
9008    FORMAT(' Sqrt-preconditioned input guess with residual norm = ',E14.6)
        x%of_r = z%of_r
     ELSE
        IF ( verbose .GE. 1 ) WRITE(environ_unit,9001)delta_en
9001    FORMAT(' Warning: bad guess with residual norm = ',E14.6,', reset to no guess')
        x%update = .FALSE.
     ENDIF

  ENDIF

  IF ( .NOT. x%update ) THEN

     x%update = .TRUE.
     x%of_r = 0.D0
     r%of_r = b%of_r - jellium
     rzold = 0.D0

  ENDIF

  ! ... Start gradient descent

  DO iter = 1, maxstep

       IF ( verbose .GE. 1 ) WRITE(environ_unit,9002) iter
9002   FORMAT(' Iteration # ',i10)

       ! ... Apply preconditioner to new state

       z%of_r = r%of_r * invsqrt%of_r
       CALL poisson_direct( core, z, z )
       z%of_r = z%of_r * invsqrt%of_r

       rznew = scalar_product_environ_density( r, z )
       IF ( ABS(rznew) .LT. 1.D-30 ) &
            & CALL errore(sub_name,'Null step in gradient descent iteration',1)

       ! ... Conjugate gradient or steepest descent input

       IF ( lconjugate .AND. ABS(rzold) .GT. 1.D-30 ) THEN
          beta = rznew / rzold
       ELSE
          beta = 0.D0
       END IF
       IF ( verbose .GE. 2 ) WRITE(environ_unit,*)'rznew = ',rznew,' rzold = ',rzold,' beta = ',beta
       rzold = rznew

       p%of_r = z%of_r + beta * p%of_r

       ! ... Apply operator to conjugate direction

       Ap%of_r = (factsqrt%of_r + electrolyte%k2 * gam%of_r ) * z%of_r + r%of_r + beta * Ap%of_r

       ! ... Step downhill

       pAp = scalar_product_environ_density( p, Ap )
       alpha = rzold / pAp

       IF ( verbose .GE. 1 ) WRITE(environ_unit,*)' pAp = ',pAp,' rzold = ',rzold,' alpha = ',alpha

       x%of_r = x%of_r + alpha * p%of_r
       r%of_r = r%of_r - alpha * Ap%of_r

       ! ... If residual is small enough exit

       delta_qm = quadratic_mean_environ_density( r )
       delta_en = euclidean_norm_environ_density( r )
       IF ( verbose .GE. 1 ) WRITE(environ_unit,9004)delta_qm,delta_en,tolvelect
9004   FORMAT(' delta_qm = ',E14.6,' delta_en = ',E14.6,' tol = ',E14.6)
       IF ( delta_en .LT. tolvelect .AND. iter .GT. 0 ) THEN
          IF ( verbose .GE. 1 ) WRITE(environ_unit,9005)
9005      FORMAT(' Charges are converged, exit!')
          EXIT
       ELSE IF ( iter .EQ. maxstep ) THEN
          WRITE(program_unit,9006)
9006      FORMAT(' Warning: Polarization charge not converged')
       ENDIF

    ENDDO

    ! In PBC the potential need to have zero average

    shift = 0.D0

    IF ( core % use_qe_fft ) THEN

       IF ( .NOT. ( core % qe_fft % use_internal_pbc_corr .OR. core % need_correction ) ) THEN

          shift = - integrate_environ_density( x ) / cell % omega

       END IF

    END IF

    x % of_r = x % of_r + shift

    IF (.not.tddfpt.AND.verbose.GE.1) WRITE(program_unit, 9007) delta_en, iter
9007 FORMAT('     polarization accuracy =',1PE8.1,', # of iterations = ',i3)

    CALL destroy_environ_density( r )
    CALL destroy_environ_density( z )
    CALL destroy_environ_density( p )
    CALL destroy_environ_density( Ap )

    CALL destroy_environ_density(invsqrt)

  RETURN

!--------------------------------------------------------------------
END SUBROUTINE linearized_pb_gradient_sqrt
!--------------------------------------------------------------------
!--------------------------------------------------------------------
SUBROUTINE linearized_pb_gradient_vacuum( gradient, core, charges, electrolyte, potential )
!--------------------------------------------------------------------

  IMPLICIT NONE

  TYPE( gradient_solver ), TARGET, INTENT(IN) :: gradient
  TYPE( electrostatic_core ), INTENT(IN) :: core
  TYPE( environ_density ), TARGET, INTENT(IN) :: charges
  TYPE( environ_electrolyte), TARGET, INTENT(IN) :: electrolyte
  TYPE( environ_density ), TARGET, INTENT(INOUT) :: potential

  TYPE( environ_cell ), POINTER :: cell
  TYPE( environ_density ), POINTER :: x, b, gam

  INTEGER :: iter
  REAL( DP ) :: rznew, rzold, alpha, beta, pAp, delta_qm, delta_en, shift
  TYPE( environ_density ) :: r, z, p, Ap

  CHARACTER( LEN=80 ) :: sub_name = 'linearized_pb_gradient_vacuum'

  LOGICAL, POINTER :: lconjugate
  INTEGER, POINTER :: maxstep
  REAL( DP ), POINTER :: tolvelect

  lconjugate => gradient % lconjugate
  maxstep => gradient % maxstep
  tolvelect => gradient % tol

  IF ( verbose .GE. 1 ) WRITE(environ_unit,9000)
9000 FORMAT(/,4('%'),' COMPUTE ELECTROSTATIC POTENTIAL ',43('%'))

  ! ... Check that fields have the same defintion domain

  IF ( .NOT. ASSOCIATED(charges%cell,electrolyte%density%cell) ) &
       & CALL errore(sub_name,'Inconsistent cells of input fields',1)
  IF ( .NOT. ASSOCIATED(charges%cell,potential%cell) ) &
       & CALL errore(sub_name,'Inconsistent cells for charges and potential',1)
  cell => charges%cell

  ! ... Aliases

  b => charges
  x => potential
  gam => electrolyte % gamma

  ! ... Create and initialize local variables

  CALL init_environ_density( cell, r )
  CALL init_environ_density( cell, z )
  CALL init_environ_density( cell, p )
  CALL init_environ_density( cell, Ap )

  ! ... Starting guess from new input and previous solution(s)

  IF ( x%update ) THEN

     x%update = .FALSE.

  ENDIF

  IF ( .NOT. x%update ) THEN

     x%update = .TRUE.
     x%of_r = 0.D0
     r%of_r = b%of_r
     rzold = 0.D0

  ENDIF

  ! ... Start gradient descent

  DO iter = 1, maxstep

       IF ( verbose .GE. 1 ) WRITE(environ_unit,9002) iter
9002   FORMAT(' Iteration # ',i10)

       ! ... Apply preconditioner to new state

       z%of_r = r%of_r
       CALL poisson_direct( core, z, z )

       rznew = scalar_product_environ_density( r, z )
       IF ( ABS(rznew) .LT. 1.D-30 ) &
            & CALL errore(sub_name,'Null step in gradient descent iteration',1)

       ! ... Conjugate gradient or steepest descent input

       IF ( lconjugate .AND. ABS(rzold) .GT. 1.D-30 ) THEN
          beta = rznew / rzold
       ELSE
          beta = 0.D0
       END IF
       IF ( verbose .GE. 2 ) WRITE(environ_unit,*)'rznew = ',rznew,' rzold = ',rzold,' beta = ',beta
       rzold = rznew

       p%of_r = z%of_r + beta * p%of_r

       ! ... Apply operator to conjugate direction

       Ap%of_r = electrolyte%k2 * gam%of_r * z%of_r + r%of_r + beta * Ap%of_r

       ! ... Step downhill

       pAp = scalar_product_environ_density( p, Ap )
       alpha = rzold / pAp

       IF ( verbose .GE. 1 ) WRITE(environ_unit,*)' pAp = ',pAp,' rzold = ',rzold,' alpha = ',alpha

       x%of_r = x%of_r + alpha * p%of_r
       r%of_r = r%of_r - alpha * Ap%of_r

       ! ... If residual is small enough exit

       delta_qm = quadratic_mean_environ_density( r )
       delta_en = euclidean_norm_environ_density( r )
       IF ( verbose .GE. 1 ) WRITE(environ_unit,9004)delta_qm,delta_en,tolvelect
9004   FORMAT(' delta_qm = ',E14.6,' delta_en = ',E14.6,' tol = ',E14.6)
       IF ( delta_en .LT. tolvelect .AND. iter .GT. 0 ) THEN
          IF ( verbose .GE. 1 ) WRITE(environ_unit,9005)
9005      FORMAT(' Charges are converged, exit!')
          EXIT
       ELSE IF ( iter .EQ. maxstep ) THEN
          WRITE(program_unit,9006)
9006      FORMAT(' Warning: Polarization charge not converged')
       ENDIF

    ENDDO

    ! In PBC the potential need to have zero average

    shift = 0.D0

    IF ( core % use_qe_fft ) THEN

       IF ( .NOT. ( core % qe_fft % use_internal_pbc_corr .OR. core % need_correction ) ) THEN

          shift = - integrate_environ_density( x ) / cell % omega

       END IF

    END IF

    x % of_r = x % of_r + shift

    IF (.not.tddfpt.AND.verbose.GE.1) WRITE(program_unit, 9007) delta_en, iter
9007 FORMAT('     polarization accuracy =',1PE8.1,', # of iterations = ',i3)

    CALL destroy_environ_density( r )
    CALL destroy_environ_density( z )
    CALL destroy_environ_density( p )
    CALL destroy_environ_density( Ap )

  RETURN

!--------------------------------------------------------------------
END SUBROUTINE linearized_pb_gradient_vacuum
!--------------------------------------------------------------------
!--------------------------------------------------------------------
END MODULE linearized_pb
!--------------------------------------------------------------------
