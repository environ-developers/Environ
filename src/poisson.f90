!--------------------------------------------------------------------
MODULE poisson
!--------------------------------------------------------------------

  USE environ_types
  USE electrostatic_types
  USE periodic
  USE environ_base, ONLY : e2, oldenviron

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: poisson_direct, poisson_gradient_direct, poisson_energy

  INTERFACE poisson_direct
     MODULE PROCEDURE poisson_direct_charges, poisson_direct_density
  END INTERFACE poisson_direct

  INTERFACE poisson_gradient_direct
     MODULE PROCEDURE poisson_gradient_direct_charges, poisson_gradient_direct_density
  END INTERFACE poisson_gradient_direct

  INTERFACE poisson_energy
     MODULE PROCEDURE poisson_energy_charges, poisson_energy_density
  END INTERFACE poisson_energy

CONTAINS
!--------------------------------------------------------------------
  SUBROUTINE poisson_direct_charges( core, charges, potential )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    TYPE( electrostatic_core ), INTENT(IN) :: core
    TYPE( environ_charges ), INTENT(IN) :: charges
    TYPE( environ_density ), INTENT(INOUT) :: potential
    !
    REAL( DP ) :: edummy, cdummy
    CHARACTER( LEN = 80 ) :: sub_name = 'poisson_direct_charges'
    !
    ! TO IMPLEMENT THE CASE OF nspin .NE. 1
    !
    potential % of_r = 0.D0
    !
    IF ( core % use_qe_fft ) THEN
       !
       CALL v_h_of_rho_r( charges%density%of_r, edummy, cdummy, potential%of_r )
       !
    ELSE IF ( core % use_oned_analytic ) THEN
       !
       CALL errore(sub_name,'Option not yet implemented',1)
       !
    ELSE
       !
       CALL errore(sub_name,'Unexpected setup of electrostatic core',1)
       !
    ENDIF
    !
    IF ( core % need_correction ) THEN
       !
       SELECT CASE ( TRIM( ADJUSTL( core%correction%type ) ) )
          !
       CASE ( '1da', 'oned_analytic' )
          !
          CALL calc_vperiodic( core%correction%oned_analytic, charges%density, potential )
          !
       CASE DEFAULT
          !
          CALL errore(sub_name,'Unexpected option for pbc correction core',1)
          !
       END SELECT
       !
    ENDIF
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE poisson_direct_charges
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE poisson_direct_density( core, charges, potential )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    TYPE( electrostatic_core ), INTENT(IN) :: core
    TYPE( environ_density ), INTENT(IN) :: charges
    TYPE( environ_density ), INTENT(INOUT) :: potential
    !
    TYPE( environ_cell ), POINTER :: cell
    TYPE( environ_density ) :: local
    !
    REAL( DP ) :: edummy, cdummy
    CHARACTER( LEN=80 ) :: sub_name = 'poisson_direct_density'
    !
    ! TO IMPLEMENT THE CASE OF nspin .NE. 1
    !
    IF ( .NOT. ASSOCIATED(charges%cell,potential%cell) ) &
         & CALL errore(sub_name,'Missmatch in domains of charges and potential',1)
    cell => charges % cell
    !
    ! ... Using a local variable for the potential because the
    !     routine may be called with the same argument for charges and potential
    !
    CALL init_environ_density( cell, local )
    !
    IF ( core % use_qe_fft ) THEN
       !
       CALL v_h_of_rho_r( charges%of_r, edummy, cdummy, local%of_r )
       !
    ELSE IF ( core % use_oned_analytic ) THEN
       !
       CALL errore(sub_name,'Option not yet implemented',1)
       !
    ELSE
       !
       CALL errore(sub_name,'Unexpected setup of electrostatic core',1)
       !
    ENDIF
    !
    ! ... PBC corrections, if needed
    !
    IF ( core % need_correction ) THEN
       !
       SELECT CASE ( TRIM( ADJUSTL( core%correction%type ) ) )
          !
       CASE ( '1da', 'oned_analytic' )
          !
          CALL calc_vperiodic( core%correction%oned_analytic, charges, local )
          !
       CASE DEFAULT
          !
          CALL errore(sub_name,'Unexpected option for pbc correction core',1)
          !
       END SELECT
       !
    ENDIF
    !
    ! ... Only update the potential at the end
    !
    potential % of_r = local % of_r
    !
    CALL destroy_environ_density( local )
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE poisson_direct_density
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE poisson_gradient_direct_charges( core, charges, gradient )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    TYPE( electrostatic_core ), INTENT(IN) :: core
    TYPE( environ_charges ), INTENT(IN) :: charges
    TYPE( environ_gradient ), INTENT(INOUT) :: gradient
    !
    CHARACTER( LEN = 80 ) :: sub_name = 'poisson_gradient_direct_charges'
    !
    ! TO IMPLEMENT THE CASE OF nspin .NE. 1
    !
    IF ( core % use_qe_fft ) THEN
       !
       CALL gradv_h_of_rho_r( charges%density%of_r, gradient%of_r )
       !
    ELSE IF ( core % use_oned_analytic ) THEN
       !
       CALL errore(sub_name,'Option not yet implemented',1)
       !
    ELSE
       !
       CALL errore(sub_name,'Unexpected setup of electrostatic core',1)
       !
    ENDIF
    !
    IF ( core % need_correction ) THEN
       !
       SELECT CASE ( TRIM( ADJUSTL( core%correction%type ) ) )
          !
       CASE ( '1da', 'oned_analytic' )
          !
          CALL calc_gradvperiodic( core%correction%oned_analytic, charges%density, gradient )
          !
       CASE DEFAULT
          !
          CALL errore(sub_name,'Unexpected option for pbc correction core',1)
          !
       END SELECT
       !
    ENDIF
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE poisson_gradient_direct_charges
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE poisson_gradient_direct_density( core, charges, gradient )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    TYPE( electrostatic_core ), INTENT(IN) :: core
    TYPE( environ_density ), INTENT(IN) :: charges
    TYPE( environ_gradient ), INTENT(INOUT) :: gradient
    !
    CHARACTER( LEN = 80 ) :: sub_name = 'poisson_gradient_direct_density'
    !
    ! TO IMPLEMENT THE CASE OF nspin .NE. 1
    !
    IF ( core % use_qe_fft ) THEN
       !
       CALL gradv_h_of_rho_r( charges%of_r, gradient%of_r )
       !
    ELSE IF ( core % use_oned_analytic ) THEN
       !
       CALL errore(sub_name,'Option not yet implemented',1)
       !
    ELSE
       !
       CALL errore(sub_name,'Unexpected setup of electrostatic core',1)
       !
    ENDIF
    !
    IF ( core % need_correction ) THEN
       !
       SELECT CASE ( TRIM( ADJUSTL( core%correction%type ) ) )
          !
       CASE ( '1da', 'oned_analytic' )
          !
          CALL calc_gradvperiodic( core%correction%oned_analytic, charges, gradient )
          !
       CASE DEFAULT
          !
          CALL errore(sub_name,'Unexpected option for pbc correction core',1)
          !
       END SELECT
       !
    ENDIF
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE poisson_gradient_direct_density
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE poisson_energy_charges( core, charges, potential, energy )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    TYPE( electrostatic_core ), INTENT(IN) :: core
    TYPE( environ_charges ), INTENT(INOUT) :: charges
    TYPE( environ_density ), INTENT(IN) :: potential
    REAL( DP ), INTENT(OUT) :: energy
    !
    REAL( DP ) :: degauss, eself
    CHARACTER( LEN = 80 ) :: sub_name = 'poisson_energy_charges'
    !
    energy = 0.d0
    !
    IF ( core % use_qe_fft ) THEN
       !
       energy = 0.5D0 * scalar_product_environ_density( charges%density, potential )
       !
    ELSE IF ( core % use_oned_analytic ) THEN
       !
       CALL errore(sub_name,'Option not yet implemented',1)
       !
    ELSE
       !
       CALL errore(sub_name,'Unexpected setup of electrostatic core',1)
       !
    ENDIF
    !
    ! Remove self-interaction and correct energy of Gaussian ions
    !
    eself = 0.D0
    degauss = 0.D0
    !
    IF ( charges % include_ions .AND. charges % ions % use_smeared_ions ) THEN
       !
       eself = charges % ions % selfenergy_correction * e2
       !
       IF ( core % use_qe_fft ) THEN
          !
          IF ( core % qe_fft % use_internal_pbc_corr .OR. core % need_correction ) THEN
             !
             degauss = 0.D0
             !
          ELSE
             !
             degauss = - charges % ions % quadrupole_correction * charges % charge * e2 * tpi &
                  & / charges % density % cell % omega
             !
          END IF
          !
       ENDIF
       !
       degauss = degauss ! we are missing the difference in interaction of electrons with gaussian vs. pc nuclei
       !
    ENDIF
    !
    energy = energy + eself + degauss
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE poisson_energy_charges
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE poisson_energy_density( core, charges, potential, energy )
!--------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    TYPE( electrostatic_core ), INTENT(IN) :: core
    TYPE( environ_density ), INTENT(IN) :: charges
    TYPE( environ_density ), INTENT(IN) :: potential
    REAL( DP ), INTENT(OUT) :: energy
    !
    REAL( DP ) :: degauss, eself
    CHARACTER( LEN = 80 ) :: sub_name = 'poisson_energy_density'
    !
    energy = 0.D0
    !
    IF ( core % use_qe_fft ) THEN
       !
       energy = 0.5D0 * scalar_product_environ_density( charges, potential )
       !
    ELSE IF ( core % use_oned_analytic ) THEN
       !
       CALL errore(sub_name,'Option not yet implemented',1)
       !
    ELSE
       !
       CALL errore(sub_name,'Unexpected setup of electrostatic core',1)
       !
    ENDIF
    !
    ! Remove self-interaction and correct energy of Gaussian ions
    ! WARNING: THESE CORRECTIONS ARE ONLY POSSIBLE WHEN INFORMATIONS
    ! ABOUT THE POINT-LIKE NUCLEI ARE PASSED
    !
    eself = 0.D0
    !
    degauss = 0.D0
    !
    energy = energy + eself + degauss
    !
    RETURN
    !
!--------------------------------------------------------------------
  END SUBROUTINE poisson_energy_density
!--------------------------------------------------------------------
!--------------------------------------------------------------------
END MODULE poisson
!--------------------------------------------------------------------
