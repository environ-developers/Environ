!
! Copyright (C) 2007-2008 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! original version by O. Andreussi and N. Marzari (MIT)
!
!----------------------------------------------------------------------------
MODULE charges_utils
  !--------------------------------------------------------------------------
  !
  USE environ_types
  USE environ_output
  USE environ_base, ONLY : e2, add_jellium
  !
  SAVE

  PRIVATE

  PUBLIC :: create_environ_charges, init_environ_charges_first, &
       & init_environ_charges_second, update_environ_charges, &
       & charges_of_potential, destroy_environ_charges

CONTAINS

!--------------------------------------------------------------------
  SUBROUTINE create_environ_charges(charges)
!--------------------------------------------------------------------

    IMPLICIT NONE

    TYPE( environ_charges ) :: charges
    CHARACTER( LEN = 80 ) :: label = 'charges'

    charges%include_ions = .FALSE.
    NULLIFY( charges%ions )

    charges%include_electrons = .FALSE.
    NULLIFY( charges%electrons )

    charges%include_externals = .FALSE.
    NULLIFY( charges%externals )

    charges%include_dielectric = .FALSE.
    NULLIFY( charges%dielectric )

    charges%include_electrolyte = .FALSE.
    NULLIFY( charges%electrolyte )

    charges%number = 0
    charges%charge = 0.D0
    CALL create_environ_density( charges%density, label )

    RETURN

!--------------------------------------------------------------------
  END SUBROUTINE create_environ_charges
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE init_environ_charges_first( charges, electrons, ions, externals, dielectric, electrolyte )
!--------------------------------------------------------------------

    IMPLICIT NONE

    TYPE( environ_charges ), INTENT(INOUT) :: charges
    TYPE( environ_electrons ), OPTIONAL, TARGET, INTENT(IN) :: electrons
    TYPE( environ_ions ),      OPTIONAL, TARGET, INTENT(IN) :: ions
    TYPE( environ_externals ), OPTIONAL, TARGET, INTENT(IN) :: externals
    TYPE( environ_dielectric ), OPTIONAL, TARGET, INTENT(IN) :: dielectric
    TYPE( environ_electrolyte ), OPTIONAL, TARGET, INTENT(IN) :: electrolyte

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

    IF ( PRESENT(dielectric) ) THEN
       charges%include_dielectric = .TRUE.
       charges%dielectric => dielectric
    ENDIF

    IF ( PRESENT(electrolyte) ) THEN
       charges%include_electrolyte = .TRUE.
       charges%electrolyte => electrolyte
    ENDIF

!--------------------------------------------------------------------
  END SUBROUTINE init_environ_charges_first
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE init_environ_charges_second( cell, charges )
!--------------------------------------------------------------------

    IMPLICIT NONE

    TYPE( environ_cell ), INTENT( IN ) :: cell
    TYPE( environ_charges ), INTENT( INOUT ) :: charges


    CALL init_environ_density( cell, charges%density )

    RETURN

!--------------------------------------------------------------------
  END SUBROUTINE init_environ_charges_second
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE update_environ_charges( charges, add_externals )
!--------------------------------------------------------------------

    IMPLICIT NONE

    TYPE( environ_charges ), INTENT( INOUT ) :: charges
    LOGICAL, INTENT(IN), OPTIONAL :: add_externals

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
       IF ( PRESENT( add_externals ) .AND. add_externals ) THEN
          charges % number = charges % number + charges % externals % number
          charges % charge = charges % charge + charges % externals % charge
          charges % density % of_r = charges % density % of_r + charges % externals % density % of_r
       END IF
    ENDIF

    local_charge = integrate_environ_density(charges%density)
    IF ( ABS(local_charge-charges%charge) .GT. 1.D-8 ) CALL errore(sub_name,'Inconsistent integral of total charge',1)

    RETURN

  END SUBROUTINE update_environ_charges

  SUBROUTINE charges_of_potential( potential, charges )

    USE dielectric,        ONLY : dielectric_of_potential
    USE electrolyte_utils, ONLY : electrolyte_of_potential

    IMPLICIT NONE

    TYPE( environ_density ), INTENT(IN) :: potential
    TYPE( environ_charges ), INTENT(INOUT) :: charges

    TYPE( environ_density ) :: tot_charge_density
    CHARACTER( LEN = 80 ) :: sub_name = 'charges_of_potential'

    CALL init_environ_density( potential % cell, tot_charge_density )

    tot_charge_density % of_r = charges % density % of_r

    IF ( charges % include_electrolyte ) THEN
       IF ( .NOT. ASSOCIATED( charges % electrolyte ) ) &
            & CALL errore(sub_name,'Missing expected charge component',1)
       ! ELECTROLYTE CHARGES ARE NOT FREE CHARGES, DO NOT ADD THEIR DENSITY TO TOTAL CHARGE
       CALL electrolyte_of_potential( potential, charges%electrolyte )
       tot_charge_density % of_r = tot_charge_density % of_r + charges % electrolyte % density % of_r
    ENDIF

    IF ( charges % include_dielectric ) THEN
       IF ( .NOT. ASSOCIATED( charges % dielectric ) ) &
            & CALL errore(sub_name,'Missing expected charge component',1)
       ! DIELECTRIC CHARGES ARE NOT FREE CHARGES, DO NOT ADD THEIR DENSITY TO TOTAL CHARGE
       CALL dielectric_of_potential( tot_charge_density, potential, charges%dielectric )
    ENDIF

    CALL destroy_environ_density( tot_charge_density )

    RETURN

  END SUBROUTINE charges_of_potential

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

       IF (ASSOCIATED(charges%dielectric)) NULLIFY( charges%dielectric )

       IF (ASSOCIATED(charges%electrolyte)) NULLIFY( charges%electrolyte )

    END IF

    CALL destroy_environ_density( charges%density )

    RETURN

!--------------------------------------------------------------------
  END SUBROUTINE destroy_environ_charges
!--------------------------------------------------------------------
END MODULE charges_utils
